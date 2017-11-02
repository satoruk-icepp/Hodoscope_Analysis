#include <math.h>
#define PI 3.14159
#define pitch 5.0
#define fibernum 84
#define mppcarray 16
#define pixelnum 12

/*recover dead channels*/
void recvdead(Double_t *EPH){
  Double_t ephdead;
  if (ch%16==0) {
    ephdead=(EPH[ch+15]+EPH[ch+1])/2.;
  }else if(ch%16==15){
    ephdead=(EPH[ch-1]+EPH[ch-15])/2.;
  }else{
    ephdead=(EPH[ch-1]+EPH[ch+1])/2.;
  }
  return ephdead;
}

/*extend the array to arbitrary length */
Double_t extarr(Double_t array[16], int ch){
  int i=ch%16;
  return array[i];
}

/*to be editted*/
int geometry(int i){
  int cpos;
  switch (i) {
    case 0:
    cpos=4;
    break;
    case 1:
    cpos=6;
    break;
    case 2:
    cpos=8;
    break;
    case 3:
    cpos=2;
    break;
    default:
    cpos=0;
    break;
  }
  return cpos;
}

/*modify array (cyclic permutation, reverse)*/
Double_t modarray(Double_t *array,int shift,int ch,Bool_t flag){
  int revch,truech;
  if (flag==true) {
    revch=15-ch;
  }else{
    revch=ch;
  }
  truech=(revch+shift)%16;
  Double_t value=array[truech];
  return value;
}

/*judge even number*/
Bool_t eventrue(int i){
  Bool_t even=true;
  if (i%2==1) {
    even=false;
  }
  return even;
}

/*return position of fiber(ch0:0mm)*/
Double_t posfib(int ch){
  int i = ch%16;
  Double_t position=(i+0.5)*5.0;
  return position;
}

/*modify the component of array based on geometry*/
void EPHprocess(Double_t *EPH,int m){
  Double_t npe[4][16]={{}};
  for (int ch = 0; ch < 64; ch++) {
    if (ch<16) {
      npe[0][ch]  = EPH[ch];
    }else if(ch<32){
      npe[1][ch%16] = EPH[ch];
    }else if(ch<48){
      npe[2][ch%16] = EPH[ch];
    }else{
      npe[3][ch%16] = EPH[ch];
    }
  }
  Double_t corEPH[64]={};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 16; j++) {
      switch (i) {

        case 0:
        EPH[i*16+j]=modarray(npe[i],0,j,false);
        break;

        case 1:
        EPH[i*16+j]=modarray(npe[i],4,j,true);
        break;
/*bundling error is taken into account*/
        case 2:
        if (m==0) {
          EPH[i*16+j]=modarray(npe[i],0,j,false);
        }else{
          EPH[i*16+j]=modarray(npe[i],4,j,true);
        }
        break;

        case 3:
        EPH[i*16+j]=modarray(npe[i],0,j,true);
        break;
      }
    }
  }
}
