#include <math.h>
#define PI 3.14159
#define pitch 5.0
#define fibernum 84
#define mppcarray 16
#define pixelnum 12

/*recover dead channels*/
Double_t recvdead(Double_t EPH[64],Int_t ch){
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
Double_t modarray(Double_t array[16],int shift,int ch,Bool_t flag){
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

/*reconstruct position with weighted mean*/
Double_t recpos(Double_t npe[16]){
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

/*determine color based on the value(purple ~ red)*/
int colordtm(Double_t value,Double_t vmin, Double_t vmax){
  if (value<vmin) {
    value=vmin;
  }
  if (value>vmax) {
    value=vmax;
  }
  double ratio=(value-vmin)/(vmax-vmin);
  double crange=49;
  double cmin=51;
  double dbcolor=cmin+crange*ratio;
  int color=floor(dbcolor);
  return color;
}

/*make parallelogram*/
void makepo(int color,Double_t xmin,Double_t ymin,Double_t width, Double_t height,Double_t degtheta,Bool_t edge){
  Double_t x[5]={xmin,xmin+height*tan((90-degtheta)*PI/180),xmin+width+height*tan((90-degtheta)*PI/180),xmin+width,xmin};
  Double_t y[5]={ymin,ymin+height,ymin+height,ymin,ymin};
  TPolyLine *pline = new TPolyLine(5,x,y);
  pline->SetFillColor(color);
  if (edge==true) {
    pline->SetLineColor(1);
    pline->SetLineWidth(4);
    pline->Draw();
  }
  pline->Draw("f");

}

/*make square*/
void makesquare(int color,Double_t xmin,Double_t ymin,Double_t size){
  TBox *b1 = new TBox(xmin,ymin,xmin+size,ymin+size);
  b1->SetFillColor(color);
  b1->SetFillStyle(1001);
  b1->Draw();
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
Double_t EPHprocess(Double_t EPH[64],int ch,int m){
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
        corEPH[i*16+j]=modarray(npe[i],0,j,false);
        break;

        case 1:
        corEPH[i*16+j]=modarray(npe[i],0,j,true);
        break;
/*bundling error is taken into account*/
        case 2:
        if (m==0) {
          corEPH[i*16+j]=modarray(npe[i],0,j,false);
        }else{
          corEPH[i*16+j]=modarray(npe[i],4,j,true);
        }
        break;

        case 3:
        corEPH[i*16+j]=modarray(npe[i],0,j,true);
        break;
      }
    }
  }
  return corEPH[ch];
}

/*color event display*/
void colorED(Double_t EPH[64],int runnum,int evenum){
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

  /*define the geometry in the TCanvas*/
  Double_t sizeCB=0.8;
  Double_t boxsize=sizeCB/fibernum;
  Double_t edge[2]={(1.0-0.8)/2.,(1.0+0.8)/2.};

  /*extend npe to 84 fibers*/
  Double_t fibnpe[4][fibernum]={{}};
  for (int i = 0; i < 4; i++) {
    for (int ch = 0; ch < fibernum; ch++) {
      fibnpe[i][ch]=extarr(npe[i],ch);
    }
  }
  Double_t summax=50;

  TString rawtext;
  rawtext.Form("#splitline{Run:%d}{Event:%d}",runnum,evenum);
  TLatex* text=new TLatex(0.15,0.5,rawtext.Data());
  text->SetTextSize(0.15);
  text->DrawClone();

  for (int i = 0; i < 84; i++) {
    for (int j = 0; j < 84; j++) {
      makesquare(colordtm((fibnpe[0][i]+fibnpe[1][i])*(fibnpe[2][j]+fibnpe[3][j]),0,summax*summax),edge[0]+j*boxsize,edge[0]+i*boxsize,boxsize);
    }
  }

  TArrow* arrow[4];
  Double_t arpos[4];
  for (int i = 0; i < 4; i++) {
    arpos[i] = edge[0]+sizeCB*recpos(npe[i])/(pitch*fibernum);
  }
  arrow[0]=new TArrow(0,arpos[0],boxsize,arpos[0],0.01,"|>");
  arrow[1]=new TArrow(1.0,arpos[1],1.-boxsize,arpos[1],0.01,"|>");
  arrow[2]=new TArrow(arpos[2],0,arpos[2],boxsize,0.01,"|>");
  arrow[3]=new TArrow(arpos[3],1,arpos[3],1.-boxsize,0.01,"|>");

  for (int i = 0; i < 4; i++) {
    arrow[i]->SetLineColor(4);
    arrow[i]->SetLineWidth(1);
    arrow[i]->DrawClone();
  }
  TMarker *e = new TMarker((arpos[2]+arpos[3])/2,(arpos[0]+arpos[1])/2,28);
  e->Draw();
}

/*slant event display of one EASIROC module*/
void slantED(Double_t EPH[64],int runnum,int evenum,Double_t edxmin,Double_t edymin, Double_t edwidth,Double_t edheight){
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

  /*define the geometry in the TCanvas*/
  Double_t deg=60;
  Double_t tanvalue=tan((90-deg)*PI/180);
  Double_t widthratio=edwidth/(edwidth+edheight*tanvalue);//ratio<1
  Double_t powidth=edwidth*widthratio/fibernum;
  Double_t poheight=edheight/fibernum;
  Double_t horinc=poheight*tanvalue;

  /*extend npe to 84 fibers*/
  Double_t fibnpe[4][fibernum]={{}};
  for (int i = 0; i < 4; i++) {
    for (int ch = 0; ch < fibernum; ch++) {
      fibnpe[i][ch]=extarr(npe[i],ch);
    }
  }
  Double_t summax=50;

  for (int i = 0; i < 84; i++) {
    for (int j = 0; j < 84; j++) {
      //makesquare(colordtm((fibnpe[0][i]+fibnpe[1][i])*(fibnpe[2][j]+fibnpe[3][j]),0,summax*summax),edge[0]+j*boxsize,edge[0]+i*boxsize,boxsize);
      makepo(colordtm((fibnpe[0][i]+fibnpe[1][i])*(fibnpe[2][j]+fibnpe[3][j]),0,summax*summax),edxmin+j*powidth+i*horinc,edymin+i*poheight,powidth,poheight,deg,false);
    }
  }
}

/*slant event display of one HBU(geometry of each tile should be modified)*/
void slantHBU(int runnum,int evenum,Double_t edxmin,Double_t edymin, Double_t edwidth,Double_t edheight,int HASIC, int HCHAN,int HADC ){
  Double_t deg=60;
  Double_t tanvalue=tan((90-deg)*PI/180);
  Double_t widthratio=edwidth/(edwidth+edheight*tanvalue);//ratio<1
  Double_t powidth=edwidth*widthratio/pixelnum;
  Double_t poheight=edheight/pixelnum;
  Double_t horinc=poheight*tanvalue;
  Int_t ASICcover=pixelnum/2;
  Double_t horincASIC=edheight*tanvalue*ASICcover/pixelnum;
  Double_t inixmin,iniymin;
  Int_t ASICarray[4]={193,194,195,196};
  makepo(0,edxmin,edymin,widthratio*edwidth,edheight,deg,true);
  for (Int_t m = 0; m < 4; m++) {
    switch (ASICarray[m]) {
      case 194:
      inixmin=edxmin;
      iniymin=edymin;
      break;
      case 193:
      inixmin=edxmin+edwidth*widthratio*ASICcover/pixelnum;
      iniymin=edymin;
      break;
      case 195:
      inixmin=edxmin+horincASIC;
      iniymin=edymin+edheight*ASICcover/pixelnum;
      break;
      case 196:
      inixmin=edxmin+horincASIC+edwidth*widthratio*ASICcover/pixelnum;
      iniymin=edymin+edheight*ASICcover/pixelnum;
      break;
    }

    for (Int_t i = 0; i < ASICcover; i++) {
      for (Int_t j = 0; j < ASICcover; j++) {
        if (ASICarray[m]==HASIC&&i==floor((HCHAN)/ASICcover)&&j==(HCHAN)%ASICcover) {
          makepo(colordtm(HADC,350,2000),inixmin+j*powidth+i*horinc,iniymin+i*poheight,powidth,poheight,deg,false);
        }else{
          makepo(0,inixmin+j*powidth+i*horinc,iniymin+i*poheight,powidth,poheight,deg,false);
        }
      }
    }
  }
}

void lateralview(Double_t corEPH[64],int HASIC,int HCHAN,int HADC) {
  TCanvas* canvas= new TCanvas("canvas","canvas",1200,600);
  canvas->Divide(2,1);
  canvas->cd(1);
  Double_t npe[4][16]={{}};
  Double_t barsize=0.8;
  Double_t barind=barsize/fibernum;
  Double_t barxmin=0.1;


  for (int ch = 0; ch < 64; ch++) {
    if (ch<16) {
      npe[0][ch]  = corEPH[ch];
    }else if(ch<32){
      npe[1][ch%16] = corEPH[ch];
    }else if(ch<48){
      npe[2][ch%16] = corEPH[ch];
    }else{
      npe[3][ch%16] = corEPH[ch];
    }
  }
  Double_t fibnpe[4][fibernum]={{}};
  for (int i = 0; i < 4; i++) {
    for (int ch = 0; ch < fibernum; ch++) {
      fibnpe[i][ch]=extarr(npe[i],ch);
    }
  }
  for (int i = 0; i <fibernum ; i++) {
    makesquare(colordtm(fibnpe[0][i],-2,48),barxmin+i*barind,0.1,barind);
  }


}

/*detailed event display of one CR counter*/
void generalED(Double_t corEPH[64],int runnum,int evenum, TCanvas* canvas,int m){

  canvas->Divide(3,3);
  Double_t npe[4][16]={{}};
  for (int ch = 0; ch < 64; ch++) {
    if (ch<16) {
      npe[0][ch]  = corEPH[ch];
    }else if(ch<32){
      npe[1][ch%16] = corEPH[ch];
    }else if(ch<48){
      npe[2][ch%16] = corEPH[ch];
    }else{
      npe[3][ch%16] = corEPH[ch];
    }
  }

  TGraph* siggraph[4];
  canvas->cd(1);
  TString eventtext;
  eventtext.Form("#splitline{Run:%d}{Event:%d}",runnum,evenum);
  TLatex* etext=new TLatex(0.15,0.5,eventtext.Data());
  etext->SetTextSize(0.15);
  etext->DrawClone();
  canvas->cd(3);
  TString HBUtext;
  /*HBUtext.Form("#splitline{ASIC:%d}{Channel:%d}",HASIC,CHAN);
  TLatex* htext=new TLatex(0.15,0.5,HBUtext.Data());
  htext->SetTextSize(0.15);
  htext->DrawClone();*/
  Double_t xpos[16];
  Double_t thpos[16];
  for (int j = 0; j < 16; j++) {
    xpos[j]=posfib(j);
    thpos[j]=(2*j+1)*PI/16;
  }
  for (int i = 0; i < 4; i++) {
    siggraph[i] = new TGraph(16,xpos,npe[i]);
    siggraph[i]->SetTitle("");
    siggraph[i]->GetXaxis()->SetTitle("Fiber position[mm]");
    siggraph[i]->GetYaxis()->SetTitle("[p.e.]");
    siggraph[i]->SetMaximum(50);
    siggraph[i]->SetMinimum(-1);
    siggraph[i]->SetFillColor(2);

    canvas->cd(geometry(i));
    siggraph[i]->Draw("AB");
    TArrow* arrow=new TArrow(recpos(npe[i]),-10,recpos(npe[i]),0,0.02,"|>");
    arrow->SetLineColor(4);
    arrow->SetLineWidth(2);
    arrow->DrawClone();
  }
  canvas->cd(5);
  colorED(corEPH,runnum,evenum);
}

void event_display(int runnum,int evenum){
  /*select file*/
  TString matchfile;
  matchfile.Form("./rootfile/match_run%d.root",runnum);
  TFile *fmatch = new TFile(matchfile.Data(), "read");
  TTree *tmatch = (TTree*)fmatch->Get("correlate");
  Int_t HTDC,HASIC,HADC,HCHAN;
  Int_t ETDC,EADC[2][64];
  Double_t EPH[2][64];

  tmatch->SetBranchAddress("HASIC", &HASIC);
  tmatch->SetBranchAddress("HTDC", &HTDC);
  tmatch->SetBranchAddress("HADC", &HADC);
  tmatch->SetBranchAddress("HCHAN", &HCHAN);
  tmatch->SetBranchAddress("ETDC", &ETDC);
  TString eadcname,ephname;
  for (int m = 0; m < 2; m++) {
    eadcname.Form("EADC%d",m+1);
    ephname.Form("EPH%d",m+1);
    tmatch->SetBranchAddress(eadcname,&EADC[m]);
    tmatch->SetBranchAddress(ephname,&EPH[m]);
  }

  TCanvas *canvas2 = new TCanvas("canvas2","scintillator plate",900,600);

  /*number of photo electron and corrected ...*/

  /*Get event*/
  tmatch->GetEntry(evenum);
  /*recover dead channel*/
  //recvdead(EPH,23);
  if (HADC>400) {
    std::cout<<"yes!!"<<std::endl;
  }
  /*convert calibrated EASIROC data to array*/
  Double_t corEPH[2][64];
  for (int m = 0; m < 2; m++) {
    for (int ch = 0; ch < 64; ch++) {
      corEPH[m][ch]=EPHprocess(EPH[m],ch,m);
    }
  }
  TCanvas* cvgeneral[2];
  for (int i = 0; i < 2; i++) {
    cvgeneral[i]= new TCanvas(Form("event display name run%d event%d easiroc%d",runnum,evenum,i),Form("event display run%d event%d easiroc%d",runnum,evenum,i),600,600);
  }

  canvas2->cd();

  slantED(corEPH[0],runnum,evenum,0.08,0.1,0.84,0.3);
  slantED(corEPH[1],runnum,evenum,0.08,0.6,0.84,0.3);
  slantHBU(runnum,evenum,0.14,0.4,0.72,0.2,HASIC,HCHAN,HADC);
  lateralview(corEPH[0],HASIC,HCHAN,HADC);
  std::cout<<"ASIC:    "<<HASIC<<"    HCHAN:    "<<HCHAN<<"    HADC:    "<<HADC<<std::endl;
  for (int m = 0; m < 2; m++) {
    generalED(corEPH[m],runnum,evenum,cvgeneral[m],m);

  }
}
