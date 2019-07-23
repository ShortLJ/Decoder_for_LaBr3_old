//#include <stdio.h>
//#include <stdlib.h>

void datareading(){
	char filename[90]="Co60_run01";


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

	UInt_t EventNumber; 
	UInt_t EventSize;
	UInt_t BoardID;
	UInt_t Group;
	UInt_t Channel;
	UInt_t TrigTimeTag;
	UInt_t waveform_length=1024;
	Float_t waveform[waveform_length];
	Float_t waveformtime[waveform_length];

	tree->Branch("EventNumber", &EventNumber);
	tree->Branch("EventSize", &EventSize);
	tree->Branch("BoardID", &BoardID);
	tree->Branch("Group", &Group);
	tree->Branch("Channel", &Channel);
	tree->Branch("TrigTimeTag",&TrigTimeTag);
	tree->Branch("waveform_length", &waveform_length);
	tree->Branch("waveform", waveform,"waveform[waveform_length]/F");
	tree->Branch("waveformtime", waveformtime,"waveformtime[waveform_length]/F");

	Float_t samplingrate = 2.5;  // GS/s
	Float_t sampletime = 1.0/samplingrate;  // ns
	for (int iwft=0; iwft<waveform_length; iwft++){
		waveformtime[iwft] = sampletime * iwft;
	}
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);

	uint32_t header[6];
	uint64_t CurrentTime, PrevRateTime, ElapsedTime;

	int ient1=0;
	for (ient1=0; ient1<numevent; ient1++){

		int Nheader=(int)fread(header, 4, 6, datafile);
		int Nwaveform=(int) fread(waveform, 4, 1024, datafile);
//		printf("%d %d\n", Nheader, Nwaveform);

		EventSize=header[0];
		BoardID=header[1];
		Group=header[2];
		Channel=header[3];
		EventNumber=header[4];
		TrigTimeTag=header[5];
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
//		if(Channel==0||Channel==1)
			tree->Fill();

//		CurrentTime = get_time();
//		ElapsedTime = CurrentTime - PrevRateTime;
//		if(ElapsedTime > 1000){
		if(ient1%4092==0){
			printf("\r.Dat->.root decoding... %d of %u (%.3f)", ient1, numevent, (float) ient1/numevent);
//			PrevRateTime = CurrentTime;
			fflush(stdout);
		}


	}
	printf("Reading completed.Entries: %d, Event numbers: %d\n", ient1,EventNumber);
	tree->Write();
	treefile->Close();

	fclose(datafile);

}
