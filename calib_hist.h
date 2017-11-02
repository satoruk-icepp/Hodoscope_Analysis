#include <math.h>
#include <vector>

double ped[64];
double gain[64];

//calibrate the gain of each MPPC in one counter
//Fit with double gaussian
void calib_hist(TH1I *hist,int ch,int m){
	TString opgausname;
	opgausname.Form("opgaus_%d",0);
	TF1 *fop = new TF1(opgausname, "[2]*exp(-0.5*((x-[1])/[3])^2)+[4]*exp(-0.5*((x-[1]-[0])/[5])^2)");
	//default pedestal value and gain value
	Double_t defpedestal;
	Double_t defopegain;
	Int_t fitmin;
	Int_t fitmax;
	//the second chip of the module ip15 is different from others
	if(m==1&&ch>32){
		defpedestal=815;
		defopegain=35;
		fitmin=800;
		fitmax=870;
	}else{
		defpedestal=800;
		defopegain=35;
		fitmin=780;
		fitmax=850;
	}

	fop->SetParameters(defopegain,defpedestal,500,7,500,7);
	hist->Fit(opgausname,"Q","",fitmin,fitmax);
	Double_t pedestal;
	Double_t opegain;
	Double_t sigma;
	opegain=fop->GetParameter(0);
	pedestal=fop->GetParameter(1);
	sigma=fop->GetParameter(3);
	//if fitting turns to be a failure, default value is returned
	if (pedestal>830||pedestal<770||opegain<15||opegain>50) {
		pedestal=defpedestal;
		opegain=defopegain;
	}
	ped[ch] = pedestal;
	gain[ch] = opegain;
	std::cout<<"pedestal:   "<<ped[ch]<<"   gain:   "<<gain[ch]<<"    sigma:    "<<sigma<<std::endl;

}
