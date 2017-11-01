#include <math.h>
#include <vector>

double ped[64];
double gain[64];

//load the data from HBU
void loadhbufromtxt(TString hbufinname, TString hbufout) {
	FILE *fp;
	TFile* fout;
	fout= new TFile(hbufout,"recreate");
	TTree* tout = new TTree("hbu","hbu");
	fp= fopen(hbufinname.Data(),"r");
	char line[4096];
	int ROC;		/* Readout cycle */
	int BXID;		/* bunch crossing ID */
	int ASIC;		/* ASIC number */
	int MEM;			/* memory cell number*/
	int CHAN;			/* channel number */
	int SPIROCTDC;			/* SPIROC SPIROCTDC */
	int SPIROCADC;			/* SPIROC SPIROCADC */
	int HIT;			/* hit bit */
	int GAIN;			/* gain bit */
	int easiroc_iterator=0;
	int easiroc_previous_iterator=0;
	int global_correlated=0;
	tout->Branch("ROC",&ROC);
	tout->Branch("BXID",&BXID);
	tout->Branch("ASIC",&ASIC);
	tout->Branch("SPIROCTDC",&SPIROCTDC);
	tout->Branch("SPIROCADC",&SPIROCADC);
	tout->Branch("CHAN",&CHAN);

	//printf("#1\t2\t3\t4\t5\t6\t7\t8\t9\n");
	//printf("#ROC\tBXID\tASIC\tMEM\tCHAN\tSPIROCTDC\tSPIROCADC\tHIT\tGAIN\n");
	while (fgets(line,4096,fp)!=NULL) {
		if (line[0]=='#') continue;
		int result;
		result=sscanf(line,"%d %d %d %d %d %d %d %d %d",&ROC,&BXID,&ASIC,&MEM,&CHAN,&SPIROCTDC,&SPIROCADC,&HIT,&GAIN);
		if (result != 9 ) {
			printf("Incorrectly parsed line. only %d values found. Skipping\n",result);
			continue;
		}
		//printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ROC,BXID,ASIC,MEM,CHAN,SPIROCTDC,SPIROCADC,HIT,GAIN);
		tout->Fill();
	}
	tout->Write();
	fout->Close();
}

//calibrate the gain of each MPPC in one counter
//Fit with double gaussian
void calibrate(TH1I *hist,int ch,int m){
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

			val = val >> 16;
			/*std::cout << "eventnum    "<<std::hex << val << std::endl;*/
			tree->GetEntry(val);
			hoge.read((char*)&val, sizeof(int));//cycle_num
			cycle_num.push_back(val & 0x00ffffff);
			//cycle_num.push_back(val & 0xffffffff);//old firmware
			cycle = val & 0xffffffff;
			/*std::cout << "cycle_num   "<<std::hex << val << std::endl;*/
			hoge.read((char*)&val, sizeof(int));//tdc_val
			tdc = val & 0xfffff;
			tdc_val.push_back(val & 0x00ffffff);
			//tdc_val.push_back(val & 0xfffff);//old firmware
			cmd_1bit.push_back(val >> 31);
			//cmd_1bit.push_back(val >> 20);
			accept = val >> 20;
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
		calibrate(charge[ch],ch,m);
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

/*merge rootfiles from two EASIROC modules in order to make it easy to match between EASIROC and HBU.*/
void mergeeasiroc(TString erfout[2],TString mergeeasirocfout){
	TFile* feasiroc[2];
	TTree* teasiroc[2];
	TTree* tph[2];
	Int_t easiroc_ROC[2],easiroc_TDC[2],easiroc_ADC[2][64];
	Double_t easirocph[2][64];
	Int_t N=0;
	for (int i = 0; i < 2; i++) {
		feasiroc[i] = new TFile(erfout[i].Data(),"read");
		teasiroc[i] = (TTree*)feasiroc[i]->Get("easiroc");
		tph[i] = (TTree*)feasiroc[i]->Get("phtree");
		teasiroc[i]->SetBranchAddress("cycle", &easiroc_ROC[i]);
		teasiroc[i]->SetBranchAddress("tdc", &easiroc_TDC[i]);
		teasiroc[i]->SetBranchAddress("easirocadc", &easiroc_ADC[i]);
		tph[i]->SetBranchAddress("easirocph", &easirocph[i]);
		N=teasiroc[i]->GetEntries();

	}

	Int_t EROC,ETDC,EADC[2][64];
	Double_t EPH[2][64];
	TFile* fout= new TFile(mergeeasirocfout.Data(),"recreate");
	TTree* tout = new TTree("measiroc","measiroc");

	tout->Branch("EADC1",&EADC[0],"EADC1[64]/I");
	tout->Branch("EADC2",&EADC[1],"EADC2[64]/I");
	tout->Branch("EPH1",&EPH[0],"EPH1[64]/D");
	tout->Branch("EPH2",&EPH[1],"EPH2[64]/D");
	tout->Branch("ETDC",&ETDC);
	tout->Branch("EROC",&EROC);
	for (int ientry = 0; ientry < N; ientry++) {
		Int_t tmpETDC;
		for (int m = 0; m < 2; m++) {
			teasiroc[m]->GetEntry(ientry);
			tph[m]->GetEntry(ientry);
			for (int ch = 0; ch < 64; ch++) {
				EADC[m][ch]=easiroc_ADC[m][ch];
				EPH[m][ch]=easirocph[m][ch];
			}

			if (m == 2&&tmpETDC==easiroc_TDC[m]) {
				std::cout<<"matching error!!"<<std::endl;
			}
			tmpETDC=easiroc_TDC[m];
			ETDC=tmpETDC;
			EROC=easiroc_ROC[m];
		}
		tout->Fill();
	}
	tout->Write();
	fout->Close();
}

//match between two modules
//If event occurs within the same bunch crossing id, then it is matched.
void match(TString hbufout, TString erfout,TString matchout){
	Int_t cycleoffset=0;

	/*HBU file extraction*/
	int hbu_ROC,hbu_tdc,hbu_BXID,hbu_adc,ASIC, CHAN;
	TFile *fhbu = new TFile(hbufout.Data(), "read");
	TTree *thbu = (TTree*)fhbu->Get("hbu");
	thbu->SetBranchAddress("ROC", &hbu_ROC);
	thbu->SetBranchAddress("BXID", &hbu_BXID);
	thbu->SetBranchAddress("ASIC", &ASIC);
	thbu->SetBranchAddress("SPIROCTDC", &hbu_tdc);
	thbu->SetBranchAddress("SPIROCADC", &hbu_adc);
	thbu->SetBranchAddress("CHAN", &CHAN);

	/*EASIROC file extraction*/
	TFile *feasiroc= new TFile(erfout.Data(),"read");
	TTree *teasiroc = (TTree*)feasiroc->Get("measiroc");
	Int_t EROC,ETDC,EADC[2][64];
	Double_t EPH[2][64];

	teasiroc->SetBranchAddress("EROC",&EROC);
	teasiroc->SetBranchAddress("ETDC", &ETDC);


	TFile* fout= new TFile(matchout.Data(),"recreate");
	TTree* tout = new TTree("correlate","correlate");

	Int_t HTDC,HASIC,HADC,HCHAN;
	Int_t CBXID, CROC;
	//TClonesArray* ermodule= new TClonesArray("")

	tout->Branch("HTDC",&HTDC);
	tout->Branch("HADC",&HADC);
	tout->Branch("HASIC",&HASIC);
	tout->Branch("HCHAN",&HCHAN);
	tout->Branch("CBXID",&CBXID);
	tout->Branch("CROC",&CROC);
	tout->Branch("ETDC",&ETDC);
	TString eadcname,ephname;
	for (int m = 0; m < 2; m++) {
		eadcname.Form("EADC%d",m+1);
		ephname.Form("EPH%d",m+1);
		teasiroc->SetBranchAddress(eadcname.Data(), &EADC[m]);
		teasiroc->SetBranchAddress(ephname.Data(), &EPH[m]);
		tout->Branch(eadcname,&EADC[m],eadcname+"[64]/I");
		tout->Branch(ephname,&EPH[m],ephname+"[64]/D");
	}
	int Nhbu=thbu->GetEntries();
	int Neasiroc=teasiroc->GetEntries();
	int ihbu=0,ieasiroc=0;
	int hbu_tfb=0;
	int startup=2215;

	while (ieasiroc<Neasiroc && ihbu<Nhbu) {
		thbu->GetEntry(ihbu);
		teasiroc->GetEntry(ieasiroc);

		int tmper=0,tmphr=0;
		/*true:increment easiroc counter*/
		Bool_t BXID_match=false;
		Bool_t ROC_match=false;
		Int_t easiroc_BXID=0;
		Bool_t hbuinc=false;
		Int_t eq_er_tdc=0;
		if (EROC-hbu_ROC==cycleoffset) {
			//ROC is matching
			easiroc_BXID=floor(ETDC-startup)/160;
			eq_er_tdc=ETDC-startup-160*easiroc_BXID;
			//restriction of BXID
			if (easiroc_BXID==hbu_BXID&& easiroc_BXID!=0) {
				CBXID=easiroc_BXID;
				CROC=EROC;
				HTDC=hbu_tdc;
				HADC=hbu_adc;
				HASIC=ASIC;
				HCHAN=CHAN;
				tout->Fill();
			}
			if (hbuinc==true) {
				ihbu++;
			}else{
				ieasiroc++;
			}
		}else if(EROC-hbu_ROC>cycleoffset){
			ihbu++;
			hbuinc=false;
		}else if(EROC-hbu_ROC<cycleoffset){
			ieasiroc++;
			hbuinc=true;
		}
	}
	tout->Write();
	fout->Close();
}

//main program
void tEASIROC(int runnum){
	TString matchout;
	TString hbufinname;
	TString hbufout;
	TString mergeeasirocfout;
	TString easirocfinname[2];
	TString easirocfout[2];
	TFile* f[2];
	for (int i = 0; i < 2; i++) {
		easirocfinname[i].Form("./data/hodoscope%dRaw_%05d.raw",i+1,runnum);
		easirocfout[i].Form("./rootfile/easiroc_%d_%05d.root",i+1,runnum);
		f[i]= new TFile(easirocfout[i].Data(), "RECREATE");
		loadeasiroc(easirocfinname[i],f[i],i);
		convert(easirocfout[i]);

	}
	mergeeasirocfout.Form("./rootfile/merge_easiroc_%05d.root",runnum);
	mergeeasiroc(easirocfout,mergeeasirocfout);
	//hbufinname.Form("./data/run%dahcal.txt",runnum);
	//hbufout.Form("./rootfile/hbu%d.root",runnum);
	//loadhbufromtxt(hbufinname,hbufout);

	//matchout.Form("./rootfile/match_run%d.root",runnum);
	//match(hbufout,mergeeasirocfout,matchout);

}
