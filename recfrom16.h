#include <math.h>
#define PI 3.14159
#define pitch 5.0
#define mppcarray 16

/*reconstruct position with weighted mean*/
Double_t recfromst(Double_t *npe){
  Double_t xsum=0;
  Double_t ysum=0;
  for (int ch = 0; ch < mppcarray; ch++) {
    xsum+=npe[ch]*cos(2*(ch+0.5)*PI/mppcarray);
    ysum+=npe[ch]*sin(2*(ch+0.5)*PI/mppcarray);
  }
  Double_t theta= atan(ysum/xsum);
  if (xsum<0) {
    theta=PI+theta;
  }
  if (theta<0) {
    theta=2*PI+theta;
  }
  return theta/(2*PI)*mppcarray*pitch;
}

/*conversion position(0~80) to theta(0~2*PI)*/
Double_t pos2theta(Double_t pos){
  return pos*(2*PI)/(mppcarray*pitch);
}
