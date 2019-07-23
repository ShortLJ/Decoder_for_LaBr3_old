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

float QDC(Float_t waveform[1024]){
	Float_t QDCsum=0;
	for (int iwft=200; iwft<1000; iwft++){
		QDCsum=QDCsum+waveform[iwft];
	}
	Float_t QDC = 4 * pedsum(waveform) - QDCsum;
	return QDC;
}

/*
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
*/


//////////////////////main function////////////////////////////////

void datareadingOnline(){
	char filename[90]="Co60_run02";
	int numchannels = 2;

	FILE *datafile = fopen(Form("%s.dat",filename),"rb");
	if(datafile==NULL){
		fputs("File error\n",stderr);
		exit(1);
	}
	fseek(datafile,0,SEEK_END);
	ULong64_t Filesize = ftell(datafile);
	rewind(datafile);
//	cout<<"Filesize = "<<Filesize<<endl;
	uint32_t numevent = Filesize/4120;
	uint32_t remainder = Filesize%4120;
	printf("Filesize: %lld(Bytes), numevent: %u, Error:%u\n", Filesize, numevent, remainder);

	TFile *treefile = new TFile(Form("%s.root",filename),"recreate");
	TTree *tree = new TTree("wavedata","wavedata");

	UInt_t CoinEventNumber=0;
	UInt_t EventSize;
	UInt_t BoardID;
	UInt_t Group;
	UInt_t Channel;
	UInt_t EventNumber[numchannels]; 
	UInt_t TrigTimeTag[numchannels];
	UInt_t waveform_length=1024;
	Float_t waveform[numchannels][waveform_length];
	Float_t waveformtime[waveform_length];
	Float_t QDC_ch0;
	Float_t QDC_ch1;
	Float_t TDC_ch0;
	Float_t TDC_ch1;


	tree->Branch("CoinEventNumber", &CoinEventNumber, "CoinEventNumber/i");
	tree->Branch("EventNumber", &EventNumber ,"EventNumber/i");
	tree->Branch("EventSize", &EventSize, "EventSize/i");
	tree->Branch("BoardID", &BoardID, "BoardID/i");
	tree->Branch("Group", &Group, "Group/i");
	tree->Branch("Channel", &Channel, "Channel/i");
	tree->Branch("TrigTimeTag",&TrigTimeTag, "TrigTimeTag/i");
	tree->Branch("waveform_length", &waveform_length);
	tree->Branch("waveform", waveform,"waveform[Channel][waveform_length]/F");
	tree->Branch("waveformtime", waveformtime,"waveformtime[waveform_length]/F");
	tree->Branch("QDC_ch0", &QDC_ch0, "QDC_ch0/F");
	tree->Branch("QDC_ch1", &QDC_ch1, "QDC_ch1/F");
	tree->Branch("TDC_ch0", &TDC_ch0, "TDC_ch0/F");
	tree->Branch("TDC_ch1", &TDC_ch1, "TDC_ch1/F");

	Float_t samplingrate = 2.5;  // GS/s
	Float_t sampletime = 1.0/samplingrate;  // ns
	for (int iwft=0; iwft<waveform_length; iwft++){
		waveformtime[iwft] = sampletime * iwft;
	}
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);

	uint32_t header[6];
	uint64_t CurrentTime, PrevRateTime, ElapsedTime;

	TCanvas *c1 = new TCanvas("Decoding","Decoding",1500,800);
	c1->Divide(2,2);
	c1->cd(1);
	TH1D *h1_QDC_ch0 = new TH1D("h1_QDC_ch0","h1_QDC_ch0;QDC;count/100QDC", 3000, 0, 300000);
	h1_QDC_ch0->GetXaxis()->SetRangeUser(5000,300000);
	h1_QDC_ch0->Draw();
	c1->cd(2);
	TH1D *h1_QDC_ch1 = new TH1D("h1_QDC_ch1","h1_QDC_ch1;QDC;count/100QDC", 3000, 0, 300000);
	h1_QDC_ch1->GetXaxis()->SetRangeUser(5000,300000);
	h1_QDC_ch1->Draw();
	c1->cd(3);
	TH1D *h1_TDC_ch0 = new TH1D("h1_TDC_ch0","h1_TDC_ch0;TDC (ns);count/0.4ns",1024,0,410);
	h1_TDC_ch0->Draw();
	c1->cd(4);
	TH1D *h1_TDC_ch1 = new TH1D("h1_TDC_ch1","h1_TDC_ch1;TDC (ns);count/0.4ns",1024,0,410);
	h1_TDC_ch1->Draw();



	int ient1=0;
	for (ient1=0; ient1<numevent; ient1++){

		int Nheader=(int)fread(header, 4, 6, datafile);
//		printf("%d %d\n", Nheader, Nwaveform);

		EventSize=header[0];
		BoardID=header[1];
		Group=header[2];
		Channel=header[3];
		EventNumber[Channel]=header[4];
		TrigTimeTag[Channel]=header[5];

		int Nwaveform=(int) fread(waveform[Channel], 4, 1024, datafile);

/*
		printf("header0,EventSize	 = %d\n",header[0]);
		printf("header1,BoardID		 = %d\n",header[1]);
		printf("header2,Group		 = %d\n",header[2]);
		printf("header3,Channel		 = %d\n",header[3]);
		printf("header4,EventCounter	 = %d\n",header[4]);
		printf("header5,TrigTimeTag	 = %d\n",header[5]);
*/
		waveform_length=Nwaveform;

		if(header[0]!=4*(Nheader+Nwaveform)){
			printf("Data Size Error!!  ient1: %d, header[0]: %d, Nheader+Nwaveform: %d\n", ient1, header[0], Nheader+Nwaveform);
			break;
		}
//printf("%f\n", waveform[1]);

		


		if(EventNumber[0]==EventNumber[1] && TrigTimeTag[0]==TrigTimeTag[1]){
			QDC_ch0 = QDC(waveform[0]);
			if(QDC_ch0<=10000)	continue;
			QDC_ch1 = QDC(waveform[1]);
			if(QDC_ch1<=10000)	continue;
			tree->Fill();
			h1_QDC_ch0->Fill(QDC_ch0);
			h1_QDC_ch1->Fill(QDC_ch1);
			h1_TDC_ch0->Fill(TDCcf(waveform[0], waveformtime, 10, 0.15));	//TDCcf(Float_t waveform[1024], Float_t waveformtime[1024],Int_t delay, double fraction)
			h1_TDC_ch1->Fill(TDCcf(waveform[1], waveformtime, 10, 0.15));	//TDCcf(Float_t waveform[1024], Float_t waveformtime[1024],Int_t delay, double fraction)
			CoinEventNumber++;
		}
		


//		CurrentTime = get_time();
//		ElapsedTime = CurrentTime - PrevRateTime;
//		if(ElapsedTime > 1000){
		if(ient1%65536==0){
			c1->cd(1);
			gPad->Modified();			
			gPad->Update();
			c1->cd(2);
			gPad->Modified();			
			gPad->Update();
			c1->cd(3);
			gPad->Modified();			
			gPad->Update();
			c1->cd(4);
			gPad->Modified();			
			gPad->Update();

			printf("\r.Dat->.root decoding... %d of %u (%.3f), CoinEventNumber %u (%.1f%)", ient1, numevent, (float) ient1/numevent, CoinEventNumber-1,(float)((float)(CoinEventNumber-1)/ient1)*100);
//			PrevRateTime = CurrentTime;
			fflush(stdout);
		}


	}
	gPad->Modified();
	gPad->Update();

	printf("\nReading completed.Entries: %d, Event numbers: %u, %u, CoinEventNumber %u\n", ient1, EventNumber[0],EventNumber[1], CoinEventNumber-1);
	tree->Write();
	treefile->Close();

	fclose(datafile);

}
