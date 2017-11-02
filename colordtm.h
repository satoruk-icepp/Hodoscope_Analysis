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
