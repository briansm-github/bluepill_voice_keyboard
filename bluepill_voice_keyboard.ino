#include <USBComposite.h> // for usb keyboard emulation
#include "codebook.h"
#include "phoneme_to_keystroke.h"

#define DEVMODE 1 // dev-mode enabled? 0 or 1.

#define WINDOW 240
#define N 80
#define ORDER 10
#define MAXSAMPLES 8000

#define FFT 32
#define FF2 15    
int bin[FF2]={16, 8,24, 4,20,12,28, 2,18,10,26, 6,22,14,30};

const int inputPin = 0;  // analog input pin PA0
int8_t s[MAXSAMPLES];  // audio sample buffer (roughly 1 second)
int samples=0;
float window[WINDOW],w[WINDOW],ave;
float realtwiddle[16],imtwiddle[16];
float R[ORDER+1];
float k[ORDER+1];
float a[ORDER+1];
char buffer[100];
int enabled=0;

USBHID HID;
HIDKeyboard Keyboard(HID);

//----------------------------------------------------------------------------------

void autocorrelate(
  float Sn[],  /* frame of Nsam windowed speech samples */
  float Rn[],  /* array of P+1 autocorrelation coefficients */
  int Nsam, /* number of windowed samples to use */
  int order /* order of LPC analysis */
)
{
  int i,j;  /* loop variables */
  float x;

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
    float       * lpc, /*      [0...p-1] LPC coefficients      */
    const float * ac,  /*  in: [0...p] autocorrelation values  */
    int p
    )
{
   int i, j;  float r, error = ac[0];

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
         float tmp  = lpc[j];
         lpc[j]     += r * lpc[i-1-j];
         lpc[i-1-j] += r * tmp;
      }
      if (i % 2) lpc[j] += lpc[j] * r;

      error *= 1.0 - r * r;
   }
}

//----------------------------------------------------------------------
void rootsofunity(float *realtwiddle, float *imtwiddle, unsigned int size)
{
    float twopi = 6.28318530717959;
    unsigned int n;

    for(n=1; n<(size>>1); n++)
    {
       realtwiddle[n] = cos(twopi*n/size);
       imtwiddle[n] = -sin(twopi*n/size);
    }
}

//-----------------------------------------------------------------------
void simple_fft(float *real, float *im, float *realtwiddle, float *imtwiddle, int size)
{
    unsigned int even, odd, span, log=0, rootindex;    // indexes
    float temp;
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
 
  float real[FFT],imag[FFT];
  float e;
  int segments=1;
  int min,max;
  int f[FF2];
  int best,d,dist,bestdist,matches;
  char result[10];
  char dump[30];
  char output[100];
  int oplen=0;
  int kscore[KEYS];

   output[0]=0;
   printf("\n--------------------\n\n");
   printf("segment %i\n",segments);
   
   for (n=0; n<samples-WINDOW; n+=N) {
     for (i=0; i<WINDOW; i++) w[i]=(float)s[i+n];
     // run de-emphasis over the window area....
     for (i=WINDOW-1; i>0; i--)  w[i]=w[i]-w[i-1]*0.975;
     // remove average....
     ave=0; for (i=0; i<WINDOW; i++) ave+=w[i]; ave/=WINDOW;
     for (i=0; i<WINDOW; i++) w[i]-=ave;
     
     for (i=0; i<WINDOW; i++) w[i]*=window[i]; // Hann window
     autocorrelate(w,R,WINDOW,ORDER);
     wld(&a[1],R,ORDER); a[0]=1.0;
     e=0;for(i=0; i<=ORDER; i++) e+=a[i]*R[i]; if (e<0) e=0;
     e=powf(e,0.16666666);
     for (i=0; i<FFT; i++) {real[i]=0; imag[i]=0;}
     for (i=0; i<=ORDER; i++) real[i]=a[i];
     simple_fft(real,imag,realtwiddle,imtwiddle,FFT);
     for (i=0; i<FF2; i++) {
        b=bin[i];
        w[i]=log(1.0/sqrt(real[b]*real[b]+imag[b]*imag[b]));
	      w[i]=floor(w[i]*4.0+7.5);
	      f[i]=(int)w[i]; 
	      if (f[i]<0) f[i]=0;
              if (f[i]>15) f[i]=15;	      
     }
     Serial.print(e); Serial.print("  ");
     if (e<2.0) for (i=0; i<FF2; i++) f[i]=0; // ignore quiet frames
 
    // search codebook for closest match(es)...
    bestdist=99999999; matches=0;
    for (c=0; c<CBSIZE; c++) {
      dist=0;
      for (i=0; i<FF2; i++) {
        d=0;
        min=cb[c][i]/16; if (f[i]<min) d=min-f[i];
	      max=cb[c][i]%16; if (f[i]>max) d=f[i]-max;
        dist+=d*d;
      }
      if (dist==0) { // phoneme hit
         result[matches]=cb[c][FF2]; matches++; result[matches]=0;
         // next line prohibits plosives appearing after the first 3 frames...
         if (n>N*5) if (cb[c][FF2]=='b' || cb[c][FF2]=='d' || cb[c][FF2]=='t') {matches--; result[matches]=0;}
      }
      if (dist<bestdist) {bestdist=dist; best=c;}
    }
    //Serial.print(matches); Serial.print("\n");
    if (matches==0) {result[0]=cb[best][FF2]; result[1]=0;}
    if (matches==1 ) {output[oplen]=result[0]; oplen++; output[oplen]=0;}
    if (oplen>1) if (output[oplen-2]==output[oplen-1]) {oplen--; output[oplen]=0;}
     
    for (i=0; i<FF2; i++) if (f[i]<10) dump[i]=f[i]+'0'; else dump[i]=f[i]-10+'a';
    dump[15]=0;
    if (DEVMODE) {
      Serial.print(dump); Serial.print("  ");
      Serial.print(result); Serial.print("   "); Serial.print(bestdist);
      Serial.print("\n");
    }
  }
  // print output hits sequence....
  Serial.print(output); Serial.print("\n");
  
  // generate a score for each keystroke possibility...
  for (i=0; i<KEYS; i++) {
    kscore[i]=0;
    for (n=0; n<oplen; n++) if (output[n]!=cbp[i][0]) kscore[i]--;
    else {kscore[i]+=50; for (n++; n<oplen; n++) if (output[n]==cbp[i][1]) {kscore[i]+=50; n+=10;} else kscore[i]-=5;}
    if (cbp[i][1]!=0 && kscore[i]<80) kscore[i]=0; 
    if (cbp[i][1]==0 && kscore[i]>0) kscore[i]+=30; // balance up single-phoneme items
  }
  d=-100; for (i=0; i<KEYS; i++) if (kscore[i]>d) {d=kscore[i]; best=i;}
  if (cbp[best][2]>2) {output[0]=cbp[best][2]; output[1]=0;}
  else if (cbp[best][2]==0) {output[0]='o'; output[1]='f'; output[2]='f'; output[3]=0;}
  else if (cbp[best][2]==1) {output[0]='o'; output[1]='n'; output[2]=0;}
  if (DEVMODE) Serial.print(output); Serial.print("\n");

  
  // lastly, generate keystroke based on pho1/pho2....
  if (!DEVMODE) {
    if (cbp[best][2]==1) {enabled=1; digitalWrite(PB12, LOW);}
    else if (cbp[best][2]==0) {enabled=0; digitalWrite(PB12, HIGH);}
    else if (enabled) Keyboard.write(cbp[best][2]);
  }    
  
  if (DEVMODE) Serial.print("--------------------------\n");// end-marker
}

//----------------------------------------------------------------------

void setup() {
  int i,j,a,ave,e;
  int sound=0,silcount=0;
 
  pinMode(inputPin, INPUT_ANALOG);
  pinMode(PB12,OUTPUT); digitalWrite(PB12, HIGH);  // LED off
  
  if (DEVMODE) Serial.begin(115200);
  else {HID.begin(HID_KEYBOARD);  Keyboard.begin();}
  
  for(i=0; i<WINDOW; i++) window[i]=0.5-0.5*cos(2.0*M_PI*i/(WINDOW-1));
  rootsofunity(realtwiddle,imtwiddle,FFT);

  while(1) {
      ave=0;
      for (i=0; i<8; i++) {
        delayMicroseconds(7);          // higher=higher pitch
        a = analogRead(inputPin); // read in the audio signal on PA0
        ave+=a;
      }
      ave/=8; ave-=2048;
      if (ave>127) ave=127; if (ave<-127) ave=-127; // avoid clipping
      
      
      //Serial.println(ave);
      
      s[samples]=ave; samples++;
      
      if (samples==200) { // check energy, reset if too low...
       	e=0; for (i=1; i<200; i++) e+=abs(s[i]-s[i-1]);
        //Serial.println(e);

        // wait for high volume to being recording....
        if (e<300) {for (i=0,j=80; i<120; i++,j++) s[i]=s[j]; samples=120;}
      }
      if (samples>200 && samples%80==40) { // check energy again....
         e=0; for (i=samples-199; i<samples; i++) e+=abs(s[i]-s[i-1]);
         //Serial.println(e); 
         // stop recording when volume drops below tolerance for 30 frames...
         if (e<200)  silcount++; else silcount=0;
      }
      
      if (samples==MAXSAMPLES || silcount==30) { // dump audio now...
        samples-=30*80; // remove trailing stuff
        process();
        samples=0; silcount=0;       
      }
  }
}

void loop() {  
 
}
