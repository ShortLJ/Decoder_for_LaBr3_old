void makingpairplf(){
	char filename[90]="Co60_run01";
	int numfiles=4;


	int pipelength=8;
	UInt_t TrigTimeTag_pipe_ch0[pipelength];
	UInt_t TrigTimeTag_pipe_ch1[pipelength];
	UInt_t Entry_pipe_ch0[pipelength];
	UInt_t Entry_pipe_ch1[pipelength];
	UInt_t EventNumber_pipe_ch0[pipelength];
	UInt_t EventNumber_pipe_ch1[pipelength];
	Int_t pipeout_ch0=0; Int_t pipein_ch0=31;
	Int_t pipeout_ch1=0; Int_t pipein_ch1=31;
	for (int i=0; i<pipelength; i++){
		TrigTimeTag_pipe_ch0[i]=-1;
		TrigTimeTag_pipe_ch1[i]=-1;
		Entry_pipe_ch0[i]=-1;
		Entry_pipe_ch1[i]=-1;
		EventNumber_pipe_ch0[i]=-1;
		EventNumber_pipe_ch1[i]=-1;
	}

	UInt_t EventNumber;
	UInt_t Channel;
	UInt_t TrigTimeTag;

	TH1D *h1_abando = new TH1D("abandoned TrigTimeTag", "abandoned TrigTimeTag;TrigTimeTag;count/640TrigTimeTag", 4000, 0, 2560000);

	TFile *file2 = new TFile(Form("%s_pair.root",filename),"recreate");
	TTree *tree2 = new TTree("Co60_pair","Co60_pair");

	ULong64_t TrigTimeTag1;
	UInt_t filenum_ch0;
	UInt_t filenum_ch1;
	UInt_t Entry_ch0;
	UInt_t Entry_ch1;
	UInt_t EventNumber_ch0;
	UInt_t EventNumber_ch1;

	tree2->Branch("TrigTimeTag",&TrigTimeTag1);
	tree2->Branch("EventNumber_ch0",&EventNumber_ch0);
	tree2->Branch("EventNumber_ch1",&EventNumber_ch1);
	tree2->Branch("filenum_ch0",&filenum_ch0);
	tree2->Branch("filenum_ch1",&filenum_ch1);
	tree2->Branch("Entry_ch0",&Entry_ch0);
	tree2->Branch("Entry_ch1",&Entry_ch1);

	ULong64_t ient=0;
	ULong64_t jent=0;
	UInt_t abando_count=0;
	UInt_t abando=0;

	TFile *file;
	TTree *tree;

	Int_t Entries[numfiles];
	ULong64_t Entriessofar=0;
	ULong64_t totalEntries = 0; 
	for( int ifile=0; ifile<numfiles; ifile++){
		file = new TFile(Form("%s_%d.root",filename,ifile),"read");
		tree = (TTree*) file->Get("wavedata");
		Entries[ifile] = tree->GetEntries();
		totalEntries = totalEntries+Entries[ifile];
		printf("%d Entries in %s_%d.root\n", Entries[ifile], filename, ifile);
	}


	for( int ifile=0; ifile<numfiles; ifile++){
		filenum_ch0=ifile;
		filenum_ch1=ifile;

		file = new TFile(Form("%s_%d.root",filename,ifile),"read");
		tree = (TTree*) file->Get("wavedata");

		tree -> SetBranchAddress("EventNumber",&EventNumber);
		tree -> SetBranchAddress("Channel", &Channel);
		tree -> SetBranchAddress("TrigTimeTag", &TrigTimeTag);


		ient=0;
		jent=0;
		abando=0;

		for (jent=0; ient<Entries[ifile]+40; jent++){
			if(TrigTimeTag_pipe_ch0[pipeout_ch0]==-1){
				pipeout_ch0 = (pipeout_ch0+1)%pipelength;	pipein_ch0 = (pipein_ch0+1)%pipelength;
			}
			if(TrigTimeTag_pipe_ch1[pipeout_ch1]==-1){
				pipeout_ch1 = (pipeout_ch1+1)%pipelength;	pipein_ch1 = (pipein_ch1+1)%pipelength;
			}
		
			if (ient<Entries[ifile]){
				tree->GetEntry(ient);
				if (Channel==0 && TrigTimeTag_pipe_ch0[pipein_ch0]==-1){
					TrigTimeTag_pipe_ch0[pipein_ch0]=TrigTimeTag;
					Entry_pipe_ch0[pipein_ch0]=ient;
					EventNumber_pipe_ch0[pipein_ch0]=EventNumber;
					ient=ient+1;
				}	
				if (Channel==1 && TrigTimeTag_pipe_ch1[pipein_ch1]==-1){
					TrigTimeTag_pipe_ch1[pipein_ch1]=TrigTimeTag;
					Entry_pipe_ch1[pipein_ch1]=ient;
					EventNumber_pipe_ch1[pipein_ch1]=EventNumber;
					ient=ient+1;
				}	
			}
			else if(ient>=Entries[ifile]){ient=ient+1;}

			if(TrigTimeTag_pipe_ch0[pipeout_ch0]!=-1 && TrigTimeTag_pipe_ch1[pipeout_ch1]!=-1){
				if(TrigTimeTag_pipe_ch0[pipeout_ch0]==TrigTimeTag_pipe_ch1[pipeout_ch1]){
					TrigTimeTag1 = TrigTimeTag_pipe_ch0[pipeout_ch0];
					Entry_ch0 = Entry_pipe_ch0[pipeout_ch0];
					Entry_ch1 = Entry_pipe_ch1[pipeout_ch1];
					EventNumber_ch0 = EventNumber_pipe_ch0[pipeout_ch0];
					EventNumber_ch1 = EventNumber_pipe_ch1[pipeout_ch1];

//					cout<<Form("data taken at %llu",TrigTimeTag_pipe_ch0[pipeout_ch0])<<endl;
					tree2->Fill();

					TrigTimeTag_pipe_ch0[pipeout_ch0]=-1;
					TrigTimeTag_pipe_ch1[pipeout_ch1]=-1;
					Entry_pipe_ch0[pipeout_ch0]=-1;
					Entry_pipe_ch1[pipeout_ch1]=-1;
					EventNumber_pipe_ch0[pipeout_ch0]=-1;
					EventNumber_pipe_ch1[pipeout_ch1]=-1;
//					pipeout_ch0 = (pipeout_ch0+1)%pipelength;	pipein_ch0 = (pipein_ch0+1)%pipelength;
//					pipeout_ch1 = (pipeout_ch1+1)%pipelength;	pipein_ch1 = (pipein_ch1+1)%pipelength;
				}
				else if(TrigTimeTag_pipe_ch0[pipeout_ch0] > TrigTimeTag_pipe_ch1[pipeout_ch1]){
					abando_count++;
					abando=TrigTimeTag_pipe_ch1[pipeout_ch1];
					h1_abando->Fill(abando);
//					h1_abando->Fill(TrigTimeTag_pipe_ch1[pipeout_ch1]);
//					cout<<Form("abando at %llu", TrigTimeTag_pipe_ch1[pipeout_ch1])<<endl;
					TrigTimeTag_pipe_ch1[pipeout_ch1]=-1;
					Entry_pipe_ch1[pipeout_ch1]=-1;
					EventNumber_pipe_ch1[pipeout_ch1]=-1;
//					pipeout_ch1 = (pipeout_ch1+1)%pipelength;	pipein_ch1 = (pipein_ch1+1)%pipelength;
				}
				else if(TrigTimeTag_pipe_ch0[pipeout_ch0] < TrigTimeTag_pipe_ch1[pipeout_ch1]){
					abando_count++;
					abando=TrigTimeTag_pipe_ch0[pipeout_ch0];
					h1_abando->Fill(abando);
//					h1_abando->Fill(TrigTimeTag_pipe_ch0[pipeout_ch0]);
//					cout<<Form("abando at %llu", TrigTimeTag_pipe_ch0[pipeout_ch0])<<endl;
					TrigTimeTag_pipe_ch0[pipeout_ch0]=-1;
					Entry_pipe_ch0[pipeout_ch0]=-1;
					EventNumber_pipe_ch0[pipeout_ch0]=-1;
//					pipeout_ch0 = (pipeout_ch0+1)%pipelength;	pipein_ch0 = (pipein_ch0+1)%pipelength;
				}	
			}
			if(ient%16384==0){
				printf("\rfilenum %d: ient %llu of %d (%.3f) jent %llu abando_count %d", ifile, ient, Entries[ifile], (float) ient/Entries[ifile], jent, abando_count);
				printf("\t\ttotal: ient %llu of %llu (%.3f) abando_count %d", Entriessofar+ient, totalEntries, (float)((float)(Entriessofar+ient)/totalEntries), abando_count);
				fflush(stdout);
			}
//			printf("\rient %llu of %d (%.3f) jent %llu abando_count %d", ient, Entries, (float) ient/Entries, jent, abando_count);
		}	//closing for(ient)
	//delete pointer?
	Entriessofar = Entriessofar+Entries[ifile]-40;
	printf("\nfilenum %d:last EventNumber %u, final ient %llu, final jent %llu(%.2f), last TrigTimeTag %u, abando count %d\n", ifile, EventNumber, ient, jent, (double) jent/ient, TrigTimeTag, abando_count);
	}	//closing for(ifile)
	printf("TrigTimeTag_pipe_ch0[pipeout_ch0] %u, TrigTimeTag_pipe_ch1[pipeout_ch1] %u\n", TrigTimeTag_pipe_ch0[pipeout_ch0], TrigTimeTag_pipe_ch0[pipeout_ch1]);
//	gStyle->SetOptStat(0);

	TCanvas *c1234 = new TCanvas ("c1aa", "asdf", 800,600);
cout<<"asdf "<<abando<<endl;
	h1_abando->Draw("");

	tree2->Write();
	file2->Close();



}
