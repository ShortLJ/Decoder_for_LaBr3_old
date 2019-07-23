double pedsum(Float_t waveform[1024]){
	int pedsum=0;
	for (int iwft=0; iwft<200;iwft++){
		pedsum = pedsum + waveform[iwft];
	}
	return pedsum;
}

double TDCle(Float_t waveform[1024],Float_t waveformtime[1024], double threshold){
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	Double_t pedestal = (double)pedsum(waveform)/200;
	for (int iwft=200; iwft<1024; iwft++){
		if(waveform[iwft]<pedestal-threshold){
			iwfttorecord = (double)iwft+(pedestal-threshold-waveform[iwft])/(waveform[iwft]-waveform[iwft-1]);
			timetorecord = (double)iwfttorecord*0.4;
			break;
		}
	}
	return timetorecord;
}

double pulseheight(Float_t waveform[1024]){
	double ADCmin=4098;
	for (int i=0; i<1024;i++){
		if (ADCmin>waveform[i]){
			ADCmin = waveform[i];
		}
	}
	double pulseheight = (double) pedsum(waveform)/200-ADCmin;
	return pulseheight;
}

double TDChf(Float_t waveform[1024],Float_t waveformtime[1024], double fraction){
	double threshold = fraction * pulseheight(waveform); 
	return TDCle(waveform,waveformtime,threshold);
}	

double TDCcf(Float_t waveform[1024], Float_t waveformtime[1024],Int_t delay, double fraction){
	double pedestal=(double)pedsum(waveform)/200;
	Double_t waveform2[1024-delay];
	for (int i=0; i<1024-delay; i++){
		waveform2[i]=(double) -fraction*(pedestal-waveform[i+delay])+(pedestal-waveform[i]);
	}
	Int_t minwft=0;
	Int_t maxwft=0;
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	for (int iwft=200; iwft<1024-delay; iwft++){
		if(waveform2[iwft]>waveform2[maxwft]){maxwft=iwft;}
		if(waveform2[iwft]<waveform2[minwft]){minwft=iwft;}
	}
	for (int iwft=minwft; iwft<maxwft; iwft++){
		if(waveform2[iwft]>0){
			iwfttorecord = (double)iwft-waveform2[iwft]/(waveform2[iwft]-waveform2[iwft-1]);
			timetorecord = (double)iwfttorecord*0.4;
			break;
		}
	}
	return timetorecord;
}

	
double TDCtrz(Float_t waveform[1024], Float_t waveformtime[1024],int trzup, int trzgap, int trzdown, double fraction){
	int trznloop = 1024-trzup-trzgap-trzdown+1;
	double trzupfactor=(double) -1.0/trzup;
	double trzdownfactor=(double) fraction/trzdown;
	double pedestal=(double)pedsum(waveform)/200;
	Double_t waveformtrz[trznloop];

	for (int itrz=0; itrz<trznloop; itrz++){
		waveformtrz[itrz]=(1.0-fraction)*pedestal;
		for (int iconvolu=0; iconvolu<trzup; iconvolu++){
			waveformtrz[itrz]=(double)waveformtrz[itrz]+trzupfactor*waveform[itrz+iconvolu];
		}
		for (int iconvolu=trzup+trzgap; iconvolu<trzup+trzgap+trzdown; iconvolu++){
			waveformtrz[itrz]=(double)waveformtrz[itrz]+trzdownfactor*waveform[itrz+iconvolu];
		}
	}
	Int_t minwft=0;
	Int_t maxwft=0;
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;

	for (int iwft=200; iwft<trznloop; iwft++){
		if(waveformtrz[iwft]>waveformtrz[maxwft]){maxwft=iwft;}
		if(waveformtrz[iwft]<waveformtrz[minwft]){minwft=iwft;}
	}
	for (int iwft=minwft; iwft<maxwft; iwft++){
		if(waveformtrz[iwft]>0){
			iwfttorecord = (double)iwft-waveformtrz[iwft]/(waveformtrz[iwft]-waveformtrz[iwft-1]);
			timetorecord = (double)iwfttorecord*0.4;
			break;
		}
	}
	return timetorecord;
}


	
float QDC(Float_t waveform[1024]){
	Float_t QDCsum=0;
	for (int iwft=200; iwft<1000; iwft++){
		QDCsum=QDCsum+waveform[iwft];
	}
	Float_t QDC = 4 * pedsum(waveform) - QDCsum;
	return QDC;
}


double energy(float QDC, int Channel){
	double energy;
	if (Channel==0){
		energy = (double) 0.0415113*(double)QDC -1.54326;
	}
	if (Channel==1){
		energy = (double) 0.0377387*(double)QDC -3.99942;
	}
	return energy;
}

void makingCoinTDCtree(){

	char filename[90] = "Co60_test9";
	int numfiles=5;

	TFile *file1[numfiles];
	TTree *tree1[numfiles];
	Int_t Entries1[numfiles];

	Int_t waveform_length = 1024;
	UInt_t EventNumber;
	UInt_t Channel;
	UInt_t TrigTimeTag;
	Float_t waveform[1024];
	Float_t waveformtime[1024];

	for(int ifile=0; ifile<numfiles; ifile++){
		file1[ifile] = new TFile(Form("%s_%d.root", filename, ifile),"read");
		tree1[ifile] = (TTree*) file1[ifile]->Get("wavedata");
		Entries1[ifile] = tree1[ifile]->GetEntries();
		printf("%d entries in %s_%d.root\n", Entries1[ifile], filename, ifile);

		tree1[ifile] -> SetBranchAddress("EventNumber",&EventNumber);
		tree1[ifile] -> SetBranchAddress("Channel", &Channel);
		tree1[ifile] -> SetBranchAddress("TrigTimeTag", &TrigTimeTag);
		tree1[ifile] -> SetBranchAddress("waveform", &waveform);
		tree1[ifile] -> SetBranchAddress("waveformtime", &waveformtime);
	}

	TFile *file2 = new TFile("Co60_test9_pair.root","read");
	TTree *tree2 = (TTree*) file2->Get("Co60_pair");
//	TTree *tree2 = (TTree*) file2->Get("Na22_pair");

	Int_t Entries2 = tree2->GetEntries();
	cout<<"Entries2 " << Entries2 << endl;
	ULong64_t TrigTimeTag2;
	ULong64_t Entry_ch0;
	ULong64_t Entry_ch1;
	Int_t filenum_ch0;
	Int_t filenum_ch1;
	ULong64_t EventNumber_ch0;
	ULong64_t EventNumber_ch1;

	tree2->SetBranchAddress("TrigTimeTag",&TrigTimeTag2);
	tree2->SetBranchAddress("Entry_ch0",&Entry_ch0);
	tree2->SetBranchAddress("Entry_ch1",&Entry_ch1);
	tree2->SetBranchAddress("filenum_ch0",&filenum_ch0);
	tree2->SetBranchAddress("filenum_ch1",&filenum_ch1);
	tree2->SetBranchAddress("EventNumber_ch0",&EventNumber_ch0);
	tree2->SetBranchAddress("EventNumber_ch1",&EventNumber_ch1);


	TFile *file3 = new TFile("Co60_test9_Pair_Coin_TDC.root","recreate");
	TTree *tree3 = new TTree("Co60_Pair_Coin_TDC","Co60_Pair_Coin_TDC");

	ULong64_t TrigTimeTag3;
	Double_t TDChf_ch0;
	Double_t TDChf_ch1;
	Double_t TDCcf_ch0;
	Double_t TDCcf_ch1;
	Double_t QDC_ch0;
	Double_t QDC_ch1;

	tree3->Branch("TrigTimeTag",&TrigTimeTag3);
	tree3->Branch("Entry_ch0",&Entry_ch0);
	tree3->Branch("Entry_ch1",&Entry_ch1);
	tree3->Branch("filenum_ch0",&filenum_ch0);
	tree3->Branch("filenum_ch1",&filenum_ch1);
	tree3->Branch("EventNumber_ch0",&EventNumber_ch0);
	tree3->Branch("EventNumber_ch1",&EventNumber_ch1);
	tree3->Branch("TDChf_ch0",&TDChf_ch0);
	tree3->Branch("TDChf_ch1",&TDChf_ch1);
	tree3->Branch("TDCcf_ch0",&TDCcf_ch0);
	tree3->Branch("TDCcf_ch1",&TDCcf_ch1);
	tree3->Branch("QDC_ch0",&QDC_ch0);
	tree3->Branch("QDC_ch1",&QDC_ch1);

	UInt_t ient=0;
	UInt_t newent=0;


	for (ient=0; ient<Entries2; ient++){
		tree2->GetEntry(ient);

		TrigTimeTag3 = TrigTimeTag2;
		tree1[filenum_ch0]->GetEntry(Entry_ch0);
		QDC_ch0 = QDC(waveform);
		if(QDC_ch0<=10000)	continue;
		TDChf_ch0 = TDChf(waveform, waveformtime, 0.15); 
		TDCcf_ch0 = TDCcf(waveform, waveformtime, 10, 0.15);

		tree1[filenum_ch1]->GetEntry(Entry_ch1);
		QDC_ch1 = QDC(waveform);
		if(QDC_ch1<=10000)	continue;
		TDChf_ch1 = TDChf(waveform, waveformtime, 0.15);
		TDCcf_ch1 = TDCcf(waveform, waveformtime, 10, 0.15);

		if(QDC_ch0>10000 && QDC_ch1>10000){
			tree3->Fill();
			newent++;
		}

		if (ient%2048 == 0){
			printf("\rCalculating TDC, QDC: %u of %d (%.3f)", ient, Entries2, (float) ient/Entries2);
			fflush(stdout);
		}
	}
	printf("\nFinished! New Entries: %u\n", newent);
	tree3->Write();
	file3->Close();		

}
