#include <USBComposite.h> // for usb keyboard emulation
#include "codebook.h"
#include "phoneme_to_keystroke.h"

#define DEVMODE 1 // dev-mode enabled? 0 or 1.

#define WINDOW 200 // 25ms analysis window (at 8000 samples/sec)
#define N 80       // 10ms step size
#define ORDER 10 // LPC order
#define MAXSAMPLES 8000 // 1 second maximum sampling (RAM limited)

#define FFT 32 // size of Fourier transform
#define FF2 15 // number of fft bins to analyse
int bin[FF2]={16, 8,24, 4,20,12,28, 2,18,10,26, 6,22,14,30};

const int inputPin = 0;  // analog input pin PA0
int8_t s[MAXSAMPLES];  // audio sample buffer (roughly 1 second)
int samples=0;
double window[WINDOW],w[WINDOW],ave;
double realtwiddle[16],imtwiddle[16];
double R[ORDER+1]; // autocorrelations  
double a[ORDER+1]; // linear prediction values
char buffer[100];
int enabled=0;

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
// set up FFT twiddles....
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
  double e;
  int segments=1;
  int min,max;
  int f[FF2];
  int best,d,dist,bestdist,matches;
  char result[10];
  char dump[30];
  char output[100];
  int oplen=0;
  int kpos,klen,score,bestscore; // for keystroke scoring

   output[0]=0;
   for (i=0; i<30; i++) dump[i]='_';
   printf("\n--------------------\n\n");
   printf("segment %i\n",segments);
   if (samples>N*30) samples=N*30; // only consider first 30 frames? 
   
   for (n=0; n<samples-WINDOW; n+=N) {
     for (i=0; i<WINDOW; i++) w[i]=(double)s[i+n];
     w[0]=0; for (i=1; i<WINDOW; i++) w[i]+=w[i-1]; // undo 1.0 pre-emph
     for (i=WINDOW-1; i>0; i--) w[i]-=w[i-1]*0.9375; // 15/16 pre-emph
     // remove average....
     ave=0; for (i=0; i<WINDOW; i++) ave+=w[i]; ave/=WINDOW;
     for (i=0; i<WINDOW; i++) w[i]-=ave;
     
     for (i=0; i<WINDOW; i++) w[i]*=window[i]; // Hann window
     autocorrelate(w,R,WINDOW,ORDER);
     wld(&a[1],R,ORDER); a[0]=1.0;
     e=0;for(i=0; i<=ORDER; i++) e+=a[i]*R[i]; if (e<0) e=0;
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
     //Serial.print(e);
     Serial.print("  ");
     
     if (e<100) for (i=0; i<FF2; i++) f[i]=0; // ignore quiet frames
 
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
         // next line gives 'n' sound privilige....
         if (cb[c][FF2]=='n') {matches=1; result[0]='n';}
      }
      if (dist<bestdist) {bestdist=dist; best=c;}
    }
    //Serial.print(matches); Serial.print("\n");
    if (matches==0) {result[0]=cb[best][FF2]; result[1]=0;}
    if (matches==1 ) {output[oplen]=result[0]; oplen++; output[oplen]=0;}
    if (oplen>1) if (output[oplen-2]==output[oplen-1]) {oplen--; output[oplen]=0;}
     
    for (i=0; i<FF2; i++) if (f[i]<10) dump[i+i/4]=f[i]+'0'; else dump[i+i/4]=f[i]-10+'a';
    dump[18]=0;
    if (DEVMODE) {
      Serial.print(dump); Serial.print("  ");
      Serial.print(result); Serial.print("   "); Serial.print(bestdist);
      Serial.print("  r="); Serial.print(cb[best][16]);
      Serial.print("\n");
    }
  }

  // print output phoneme-hits sequence....
  Serial.print(output); Serial.print("\n");

  // search for best keystroke possibility...
  bestscore=0; best='?';  klen=cbk[0]; kpos=1;
  while(klen!=0) {
    score=0; i=0; n=0;
    while(n<oplen) {
       if (output[n]==cbk[kpos]) {
          score+=10; kpos++; i++;
          if (i==klen) {n=oplen-1; if (score>bestscore) {bestscore=score; best=cbk[kpos];}}
       }
       else score--;
       n++; if (n==oplen) kpos-=i; // after over-run of output buffer, rewind kpos
    }
    kpos+=klen+1; klen=cbk[kpos]; kpos++; // move to next entry in cbk
  }
 
  if (best>2 && best!=32) {output[0]=best; output[1]=0;}
  else if (best==0) {output[0]='o'; output[1]='f'; output[2]='f'; output[3]=0;}
  else if (best==1) {output[0]='o'; output[1]='n'; output[2]=0;}
  else if (best==32) {output[0]='_'; output[1]=0;}
  if (DEVMODE) Serial.print(output); Serial.print("\n");

  
  // lastly, generate keystroke (or switch keystroke generation off/on)
  if (!DEVMODE) {
    if (best==1) {enabled=1; digitalWrite(PB12, LOW);}
    else if (best==0) {enabled=0; digitalWrite(PB12, HIGH);}
    else if (enabled && best!='?') Keyboard.write(best);
  }    

  
  if (DEVMODE) Serial.print("--------------------------\n");// end-marker
}

//----------------------------------------------------------------------

void setup() {
  int i,j,a,b,ave,e,last=0;
  int sound=0,silcount=0;
 
  pinMode(inputPin, INPUT_ANALOG);
  pinMode(PB12,OUTPUT); digitalWrite(PB12, HIGH);  // LED off
  
  if (DEVMODE) Serial.begin(115200);
  else {HID.begin(HID_KEYBOARD);  Keyboard.begin();}
  
  for(i=0; i<WINDOW; i++) window[i]=0.5-0.5*cos(2.0*M_PI*i/(WINDOW-1));
  rootsofunity(realtwiddle,imtwiddle,FFT);

  while(1) {
      ave=0;
      for (i=0; i<8; i++) { // read  8 samples and take average value
        delayMicroseconds(7);          // higher=higher pitch
        a = analogRead(inputPin); // read in the audio signal on PA0
        a-=2048;  ave+=a; // 12-bit sample with 11-bit (2048) bias
      }
      ave/=8; 
      
      if (ave-last>127) ave=last+127;  // lowpass + avoid clipping.
      if (ave-last<-127) ave=last-127; 
    
      b=ave-last; last=ave;
    
      //Serial.println(ave);
      
      s[samples]=b; samples++;
      
      if (samples==WINDOW) { // check energy, reset if too low...
       	e=0; for (i=1; i<WINDOW; i++) e+=abs(s[i]);
        //Serial.println(e);

        // wait for high volume to being recording....
        if (e<1000) {for (i=0,j=80; i<120; i++,j++) s[i]=s[j]; samples=120;}
      }
      if (samples>WINDOW && samples%N==40) { // check energy again....
         e=0; for (i=samples-199; i<samples; i++) e+=abs(s[i]);
         //Serial.println(e); 
         // stop recording when volume drops below tolerance for 30 frames...
         if (e<400)  silcount++; else silcount=0;
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
