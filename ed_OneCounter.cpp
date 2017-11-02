#include <math.h>
#include <recfrom16.h>
#include <colordtm.h>
#include <chassign.h>
#include <makefig.h>
#define PI 3.14159
#define pitch 5.0
#define fibernum 84
#define mppcarray 16
#define pixelnum 12

/*color event display*/
void colorED(Double_t *EPH,int runnum,int evenum){
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
    arpos[i] = edge[0]+sizeCB*recfromst(npe[i])/(pitch*fibernum);
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

/*detailed event display of one CR counter*/
void SquareED(Double_t *corEPH,int runnum,int evenum){
  TCanvas* canvas=new TCanvas("canvas","canvas",600,600);
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
    TArrow* arrow=new TArrow(recfromst(npe[i]),-10,recfromst(npe[i]),0,0.02,"|>");
    arrow->SetLineColor(4);
    arrow->SetLineWidth(2);
    arrow->DrawClone();
  }
  canvas->cd(5);
  colorED(corEPH,runnum,evenum);
}

void ed_OneCounter(int runnum,int module,int evenum){
  /*select file*/
  TString rootfile;
  rootfile.Form("../rootfile/easiroc_%d_%05d.root",module,runnum);
  TFile *fin = new TFile(rootfile.Data(), "read");
  TTree *tin = (TTree*)fin->Get("easiroc");
  Int_t easiroc_ROC,easiroc_TDC,easiroc_ADC[64];
  TString eadcname,ephname;
  for (int m = 0; m < 2; m++) {
    tin->SetBranchAddress("cycle", &easiroc_ROC);
    tin->SetBranchAddress("tdc", &easiroc_TDC);
    tin->SetBranchAddress("easirocadc", easiroc_ADC);
  }


  /*number of photo electron and corrected ...*/

  /*Get event*/
  tin->GetEntry(evenum);
  Double_t temporalEPH[64];
  Double_t MPPCGain[64];
  Double_t MPPCPedestal[64];
  for(int i=0;i<64;i++){
    MPPCGain[i]=30;
    MPPCPedestal[i]=780;
    temporalEPH[i]=((double)easiroc_ADC[i]-MPPCPedestal[i])/MPPCGain[i];
  }
  EPHprocess(temporalEPH,1);
  SquareED(temporalEPH,runnum,evenum);

}
