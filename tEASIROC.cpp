#include <math.h>
#include <calib_hist.h>
#include <vector>

void loadeasiroc(TString easirocfinname,TFile* f,int m){
	gROOT->SetStyle("Plain");
	double mean = 0;
	double RMS = 0;
	std::vector<int> cycle_num;
	std::vector<int> tdc_val;
	std::vector<int> cmd_1bit;

	std::ifstream hoge(easirocfinname.Data(),std::ios::binary);
	if(!hoge.is_open()){
		cout << "no file" << endl;
		return;
	}
	f->cd();
	TTree *tree = new TTree("easiroc","easiroc");
	Int_t cycle;
	Int_t tdc;
	Int_t easirocadc[64];
	Int_t accept;
	Double_t pedestal[64];
	tree->Branch("easirocadc",&easirocadc,"easirocadc[64]/I");
	tree->Branch("cycle",&cycle);
	tree->Branch("tdc",&tdc);
	tree->Branch("accept",&accept);
	TH2I* allADC = new TH2I("all_ADC", "all_ADC_title", 64, 0, 64, 300, 700, 1000);
	TH1I *charge[64];
	TString chhistname[64],chhisttitle[64];
	for (int ch = 0; ch < 64; ch++) {
		chhistname[ch].Form("ADC_%d",ch);
		chhisttitle[ch].Form("ADC_histo_%d",ch);
		charge[ch] = new TH1I(chhistname[ch],chhisttitle[ch],300,700,1000);
	}
	while(!hoge.eof()){
		UInt_t val;
		hoge.read((char*)&val, sizeof(int));
		//std::cout <<"header[0]   "<< std::hex << val << std::endl;
		if(val == 0xffff7368){
			hoge.read((char*)&val, sizeof(int));//wordsize
			//std::cout << "wordsize    "<<std::hex << val << std::endl;
			
			hoge.read((char*)&val, sizeof(int));//eventnum
			val = (val >> 16) & 0x0fff;
			//std::cout << "eventnum    "<<std::hex << val << std::endl;
			//tree->GetEntry(val);
			
			hoge.read((char*)&val, sizeof(int));//cycle_num
			cycle = val & 0x00ffffff;
			cycle_num.push_back(cycle);
			//cycle_num.push_back(val);//old firmware
			//std::cout << "cycle_num   "<<std::hex << val << std::endl;
			
			hoge.read((char*)&val, sizeof(int));//tdc_val
			tdc = val & 0x00ffffff;
			tdc_val.push_back(tdc);
			//tdc_val.push_back(val & 0x000fffff);//old firmware
			cmd_1bit.push_back((val >> 31) & 0x1);
			accept = (val >> 31) & 0x1;
			//std::cout << "tdc_val     "<<std::hex << val << std::endl;
			
			for(int i = 0; i<64; ++i){
				hoge.read((char*)&val, sizeof(int));
				//std::cout << std::hex << val << std::endl;
				easirocadc[i] = val & 0xffff;

				charge[i]->Fill(easirocadc[i]);
				allADC->Fill(i,easirocadc[i]);

			}
		}
		tree->Fill();
	}

	for (int ch = 0; ch < 64; ch++) {
		calib_hist(charge[ch],ch,m);
		charge[ch]->Write();
	}
	allADC->Write();
	tree->Write();
	//f->Write();
	f->Close();
}

/*convert ADC to photoelectron*/
void convert(TString easirocfout){
	TFile *feasiroc = new TFile(easirocfout.Data(), "update");
	std::cout<<easirocfout<<std::endl;
	TTree *phtree = new TTree("phtree","phtree");
	TTree *teasiroc = (TTree*)feasiroc->Get("easiroc");
	Int_t easirocadc[64];
	Double_t easirocph[64];
	teasiroc->SetBranchAddress("easirocadc", &easirocadc);
	phtree->Branch("easirocph",&easirocph,"easirocph[64]/D");
	TH1I* photon = new TH1I("photon","photon title",100,0,20);
	int N= teasiroc->GetEntries();
	/*std::cout<<"total events:   "<<N<<"    pedestal:   "<<ped[0]<<"   gain:   "<<gain[0]<<std::endl;*/
	for (int i = 0; i < N; i++) {
		teasiroc->GetEntry(i);
		for (int ch = 0; ch < 64; ch++) {
			easirocph[ch] = ((double)easirocadc[ch]-ped[ch])/gain[ch];
		}
		photon->Fill(easirocph[0]);
		phtree->Fill();
	}
	photon->Write();
	phtree->Write();
	feasiroc->Write();
	feasiroc->Close();
}

//main program
void tEASIROC(int runnum){
	TString easirocfinname[2];
	TString easirocfout[2];
	TFile* f[2];
	for (int i = 0; i < 2; i++) {
		easirocfinname[i].Form("../data/hodoscope%dRaw_%05d.raw",i+1,runnum);
		easirocfout[i].Form("../rootfile/easiroc_%d_%05d.root",i+1,runnum);
		f[i]= new TFile(easirocfout[i].Data(), "RECREATE");
		loadeasiroc(easirocfinname[i],f[i],i);
		convert(easirocfout[i]);
	}
}
