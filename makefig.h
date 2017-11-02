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
