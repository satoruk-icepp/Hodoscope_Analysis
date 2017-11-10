#include <math.h>
#define PI 3.14159
#define PITCH 5.0
#define NFIBER 84
#define NCH 16

/*recover dead channels*/
void recvdead(Double_t *EPH){
  Double_t ephdead;
    ephdead=(EPH[(ch+NCH-1)%NCH]+EPH[(ch+1)%NCH])/2.;

  return ephdead;
}

/*extend the array to arbitrary length */
Double_t extarr(Double_t array[NCH], int ch){
  int i=ch%NCH;
  return array[i];
}

Int_t fib2chn(Int_t fiber, Int_t n, Int_t side) {
  //n: AC=0, BD=1, EG=2, FH=3
  //side: ABEF=0, CDGH=1
  if (fiber==72 && n==1 && side==0) fiber=73;
  else if (fiber==73 && n==1 && side==0) fiber=72;
  else if (fiber==75 && n==1 && side==0) fiber=76;
  else if (fiber==76 && n==1 && side==0) fiber=75;
  switch(n) {
    case 0:
      if (side==0)	return fiber%16;
      else		return 31-fiber%16;
    case 1:
      if (side==0)	return 47-fiber%16;
      else		return 48+fiber%16;
    case 2:
      if (side==0)	return fiber%16;
      else		return 16+(fiber+12)%16;
    case 3:
      if (side==0)	return 47-fiber%16;
      else		return 48+fiber%16;
  }
  return 0;
}

/*judge even number*/
Bool_t eventrue(int i){
  return (i%2==0);
}

/*return position of fiber(ch0:2.5mm)*/
Double_t posfib(int ch){
  int i = ch%NCH;
  Double_t position=(i+0.5)*PITCH;
  return position;
}
