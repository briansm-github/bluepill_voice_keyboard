#include <USBComposite.h> // for usb keyboard emulation
#include "cb.h"

#define WINDOW 240
#define N 80
#define ORDER 10

#define FFT 16
#define FF2 9    
int bin[FF2]={0,8,4,11,2,13,5,14,1};

const int inputPin = 0;  // analog input pin PA0
int8_t s[8000];  // audio sample buffer (roughly 1 second)
int samples=0;
double window[WINDOW],w[WINDOW],ave;
double realtwiddle[16],imtwiddle[16];
double R[ORDER+1];
double k[ORDER+1];
double a[ORDER+1];
char buffer[100];

USBHID HID;
HIDKeyboard Keyboard(HID);

//----------------------------------------------------------------------------------

void autocorrelate(
  double Sn[],  /* frame of Nsam windowed speech samples */
  double Rn[],  /* array of P+1 autocorrelation coefficients */
  int Nsam, /* number of windowed samples to use */
  int order /* order of LPC analysis */
)
{
  int i,j;  /* loop variables */
  double x;

  for(j=0; j<order+1; j++) {
    x = 0.0;
    for(i=0; i<Nsam-j; i++) x+=Sn[i]*Sn[i+j];
    Rn[j]=x;
  }
  if (Rn[0]<1.0) Rn[0]=1.0; // prevent div-by-0 later
}


//----------------------------------------------------------------------------
// Levinson-Durbin function based on Speex 1.0.5 lpc.c....
void wld(
    double       * lpc, /*      [0...p-1] LPC coefficients      */
    const double * ac,  /*  in: [0...p] autocorrelation values  */
    int p
    )
{
   int i, j;  double r, error = ac[0];

   for (i = 0; i < p; i++) {

      /* Sum up this iteration's reflection coefficient.
       */
      r = -ac[i + 1];
      for (j = 0; j < i; j++) r -= lpc[j] * ac[i - j];
      r /= error;

      /*  Update LPC coefficients and total error.
       */
      lpc[i] = r;
      for (j = 0; j < i/2; j++) {
         double tmp  = lpc[j];
         lpc[j]     += r * lpc[i-1-j];
         lpc[i-1-j] += r * tmp;
      }
      if (i % 2) lpc[j] += lpc[j] * r;

      error *= 1.0 - r * r;
   }
}

//----------------------------------------------------------------------
void rootsofunity(double *realtwiddle, double *imtwiddle, unsigned int size)
{
    double twopi = 6.28318530717959;
    unsigned int n;

    for(n=1; n<(size>>1); n++)
    {
       realtwiddle[n] = cos(twopi*n/size);
       imtwiddle[n] = -sin(twopi*n/size);
    }
}

//-----------------------------------------------------------------------
void simple_fft(double *real, double *im, double *realtwiddle, double *imtwiddle, int size)
{
    unsigned int even, odd, span, log=0, rootindex;    // indexes
    double temp;
    log=0;
    int i;

    for(span=size>>1; span; span>>=1, log++)   
    {
       
        for(odd=span; odd<size; odd++)         // iterate over the dual nodes
        {
      
            odd |= span;                    // iterate over odd blocks only
            even = odd ^ span;              // even part of the dual node pair
            //printf("even=%i,odd=%i\n",even,odd);
                       
            temp = real[even] + real[odd];       
            real[odd] = real[even] - real[odd];
            real[even] = temp;
           
            temp = im[even] + im[odd];           
            im[odd] = im[even] - im[odd];
            im[even] = temp;
           
            rootindex = (even<<log) & (size-1); // find root of unity index
            if(rootindex)                    // skip rootindex[0] (has an identity)
            {
                temp=realtwiddle[rootindex]*real[odd]-imtwiddle[rootindex]*im[odd];
                im[odd]=realtwiddle[rootindex]*im[odd]+imtwiddle[rootindex]*real[odd];
                real[odd] = temp;
            }
   
        } // end of loop over n
     
     } // end of loop over FFT stages

} //end of function

//----------------------------------------------------------------------

void process() {
  int i,n,b,c;
 
  double real[FFT],imag[FFT];
  double min,max,snr,prev;
  int segments=1;
  int m;
  int f[FF2];
  int best,d,dist,bestdist;
  char out[2],pho1=0,pho2=0;

    
   printf("\n--------------------\n\n");
   printf("segment %i\n",segments);
   for (n=0; n<samples-WINDOW; n+=N) {
   for (i=0; i<WINDOW; i++) w[i]=(double)s[i+n];
   // run de-emphasis over the window area....
   prev=0; for (i=0; i<WINDOW; i++) {w[i]+=prev/2; prev=w[i];}
   ave=0; for (i=0; i<WINDOW; i++) ave+=w[i]; ave/=WINDOW;
   for (i=0; i<WINDOW; i++) w[i]-=ave;
   for (i=0; i<WINDOW; i++) w[i]*=window[i];
   autocorrelate(w,R,WINDOW,ORDER);
   wld(&a[1],R,ORDER); a[0]=1.0;
   for (i=0; i<FFT; i++) {real[i]=0; imag[i]=0;}
   for (i=0; i<=ORDER; i++) real[i]=a[i];
   simple_fft(real,imag,realtwiddle,imtwiddle,FFT);
   for (i=0; i<FF2; i++) {
           b=bin[i];
           w[i]=log(1.0/sqrt(real[b]*real[b]+imag[b]*imag[b]))*50.0;          
    }
    // find SNR....   
    min=999999; for (i=0; i<FF2; i++) if (w[i]<min) min=w[i];
    max=-1000; for (i=0; i<FF2; i++) if (w[i]>max) max=w[i];
    snr=max-min; 
         
    // ignore low-SNR frames, otherwise normalize...
    if (snr>100) {
      for (i=0; i<FF2; i++) w[i]-=min;  
    } else for (i=0; i<FF2; i++) w[i]=0;    

    // find closest codebook match now.....
    for (i=0; i<FF2; i++) f[i]=(int)w[i];
    bestdist=999999;
    for (c=0; c<CBSIZE; c++) {
      dist=0; for (i=0; i<FF2; i++) {d=f[i]-cb[c][i]; dist+=d*d;}
      if (dist<bestdist) {bestdist=dist; best=c;}
    }
    out[0]=cb[best][FF2]; out[1]=0;
    if (bestdist<2000) {
      if (pho1==0) pho1=cb[best][FF2];
      else if (pho2==0 && cb[best][FF2]!=pho1) pho2=cb[best][FF2];
    }
    
    for (i=0; i<FF2; i++) {Serial.print((int)w[i]); Serial.print(",");}
    Serial.print(out); Serial.print("   "); Serial.print(bestdist);
    Serial.print("\n");
  }
  // lastly, generate keystroke based on pho1/pho2....
  printf("pho1=%c, pho2=%c\n",pho1,pho2);
  for (i=0; i<KEYS; i++) if (pho1==cbp[i][0])
    if (cbp[i][1]==0 || pho2==cbp[i][1]) {
      printf("MATCH = %c (%i)\n",cbp[i][2],i);
      Keyboard.write(cbp[i][2]);
    }
    
  
  Serial.print("--------------------------\n");// end-marker
}
//----------------------------------------------------------------------

void setup() {
  int i,j,a,b,last=2048,ave,e;
  int sound=0,silcount=0;
  float x=1.2345,y=3.456,z;
  
  pinMode(inputPin, INPUT_ANALOG);
  
  HID.begin(HID_KEYBOARD);
  Keyboard.begin();
  
  //Serial.begin(115200); 
  
  for(i=0; i<WINDOW; i++) window[i]=0.5-0.5*cos(2.0*M_PI*i/(WINDOW-1));
  rootsofunity(realtwiddle,imtwiddle,FFT);

  while(1) {
      ave=0;
      for (i=0; i<8; i++) {
        delayMicroseconds(7);          // higher=higher pitch
        a = analogRead(inputPin); // read in the audio signal on PA0
        ave+=a;
      }
      ave/=8;
      
      b=ave-last; last=ave; // pre-emphasis
      //Serial.println(b);
      
      s[samples]=b; samples++;
      
      if (samples==200) { // check energy, reset if too low...
        e=0; for (i=0; i<200; i++) e+=s[i]*s[i];
        //Serial.println(e);

        // wait for high volume to being recording....
        if (e<1000) {for (i=0,j=80; i<120; i++,j++) s[i]=s[j]; samples=120;}
      }
      if (samples>200 && samples%80==40) { // check energy again....
         e=0; for (i=samples-200; i<samples; i++) e+=s[i]*s[i];
         //Serial.println(e); 
         // stop recording when volume drops below tolerance for 30 frames...
         if (e<750) {samples-=80; silcount++;} else silcount=0;
      }
      
      if (samples==10000 || silcount==30) { // dump audio now...
        //for (i=0; i<samples; i++) {Serial.print(s[i]); Serial.print("\n");}
        //Serial.print(65535); Serial.print("\n");// end-marker
        process();
        samples=0; silcount=0;       
      }
  }
}

void loop() {  
 
}
