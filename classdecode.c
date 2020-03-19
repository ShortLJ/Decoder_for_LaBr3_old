#include <iostream>
#include <thread>
#include "TTree.h"
#include "TFile.h"
#include "./include/keyb.c"
#include "./include/gettimef.c"
#include "./include/Pulse1024.h"
#include "./include/DecoderCanvas.h"


void KeyboardCommand(int pqrj_flag[3]){		
	if(!kbhit())	return;
	int key = getch();
	switch(key){
		case 'p' :						// 'p'ause
			pqrj_flag[0]^=1; 
			if( pqrj_flag[0]) printf("\nPaused. Press p again to resume\n");
			if(!pqrj_flag[0]) printf("resumed\n");
			break;
		case 'q' : 
			pqrj_flag[1]^=1; break;			// 'q'uit
		case 'Q' : 
			pqrj_flag[1]=2; break;			// 'q'uit
		case 'r' : 						// 'r'enew histograms
			pqrj_flag[2]^=1; 
			if( pqrj_flag[2]) printf("\nRenewing histograms is enabled.\n");
			if(!pqrj_flag[2]) printf("\nRenewing histograms is disabled.\n");
			break;
		case 'j' : 
			pqrj_flag[3]^=1; break;			// saving canvas to filename_canvasname.'j'pg
		default : break;
	}
}
//////////////////////main function////////////////////////////////

void classdecode(){
	//char filename[90]="PMT_9A_primary_Na22_run03_HV_1040_1050_Noatt_2500MHz_1kSa_asym1";
	//char filename[90]="PMT_9A_primary_Na22_run02_HV_1040_1050_Noatt_1GHz_1kSa";
	//char filename[90]="PMT_9A_primary_Na22_run01_HV_1040_1050_Noatt_2.5GHz";
	char filename[90]="PMT_primary_TDC_setvoltage_set2";
	//char filename[90]="PMT_primary_test_setvoltage_set1";
	//char filename[90]="Eu152_att66_newcalib_run02";
	//char filename[90]="NoSource_att66_newcalib_run01";
	//char filename[90]="Ba133_att66_newcalib_run02";
	//char filename[90]="Na22_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run02";
	FILE* datafile = fopen(Form("%s.dat",filename),"rb");	// moved from latter to here
	if(datafile==NULL){
		fputs("File error\n",stderr);
		exit(1);
	}

//tree
	const UInt_t NOpCh = 4;	// number of open channels
	const UInt_t Nthread = 1;	// number of threads
	const UInt_t waveform_length = 1024;
	const Float_t samplingrate = 2.5;  // GS/s
	const Float_t sampletime = 1.0/samplingrate;  // ns
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);

	const int sample_byte = 4*(waveform_length+6);
	fseek(datafile,0,SEEK_END);
	ULong64_t Filesize = ftell(datafile);
	rewind(datafile);
	uint32_t numevent = Filesize/sample_byte;
	uint32_t remainder = Filesize%sample_byte;
	printf("%s.dat Filesize: %lld(Bytes), numevent: %u, Error:%u\n", filename, Filesize, numevent, remainder);

	Float_t waveformBuffer[waveform_length];

	TFile* treefile = new TFile(Form("../../DecodedDataFiles/%s_multithread.root",filename),"recreate");
	TTree* tree = new TTree("wavedata","wavedata");
	UInt_t EventSize;
	UInt_t BoardID=0;
	UInt_t StartIndexCell;
	//UInt_t Group;
	UInt_t NHitCh = NOpCh;	//number of hit channels, maximum NOpCh
	UInt_t Channel[NHitCh];
	UInt_t EventNumber;
	UInt_t TrigTimeTag;
	Float_t waveform[NOpCh][waveform_length];
	//Float_t waveformtime[waveform_length];
	Float_t rPH[NHitCh];
	Float_t QDC[NHitCh];
	Float_t TDC[NHitCh];

	tree->Branch("EventNumber", &EventNumber, "EventNumber/i");
	tree->Branch("TrigTimeTag", &TrigTimeTag, "TrigTimeTag/i");
	tree->Branch("EventSize", &EventSize, "EventSize/i");
	tree->Branch("BoardID", &BoardID, "BoardID/i");
	tree->Branch("StartIndexCell", &StartIndexCell, "StartIndexCell/i");
	//tree->Branch("Group", &Group, "Group/i");
	tree->Branch("NHitCh", &NHitCh, "NHitCh/i");
	tree->Branch("Channel", Channel, "Channel[NHitCh]/i");
//	tree->Branch("waveform_length", &waveform_length, "waveform_length/i");
//	tree->Branch("waveform", waveform,"waveform[NHitCh][waveform_length]/F");
//	tree->Branch("waveformtime", waveformtime,"waveformtime[waveform_length]/F");
	tree->Branch("rPH", rPH, "rPH[NHitCh]/F");
	tree->Branch("QDC", QDC, "QDC[NHitCh]/F");
	tree->Branch("TDC", TDC, "TDC[NHitCh]/F");


//canvas
	const Int_t QDCrange = 300000/2.5*samplingrate;
	Double_t QDCcut[NOpCh][2] = {{8.40333e+04,1.35843e+05},{8.48711e+04,1.22737e+05},{8.48711e+04,1.22737e+05},{8.48711e+04,1.22737e+05}};

	DECODER_CANVAS *canvases = new DECODER_CANVAS(NOpCh);
	canvases->Construct();	

	DECODER_PRIMARY_H1 *h1_primary = new DECODER_PRIMARY_H1(NOpCh, QDCrange, waveform_length, sampletime);
	h1_primary->Draw(canvases);
	DECODER_QDC_H1_COIN *h1_QDC_ch[NOpCh];
	DECODER_QDC_H2 *h2_QDC_xaxis[NOpCh];
	DECODER_TDCDIF_H1 *h1_TDCdif[NOpCh];
	for (int ich=0; ich<NOpCh; ich++){
		h1_QDC_ch[ich] = new DECODER_QDC_H1_COIN(NOpCh, ich, QDCrange);
		h1_QDC_ch[ich]->Draw(canvases);
		h2_QDC_xaxis[ich] = new DECODER_QDC_H2(NOpCh, ich, QDCrange);
		h2_QDC_xaxis[ich]->Draw(canvases);
		h1_TDCdif[ich] = new DECODER_TDCDIF_H1(NOpCh, ich, -3, 3);
		h1_TDCdif[ich]->Draw(canvases);
	}

//////////////////////

	uint32_t header[6];
	int Nheader, Nwaveform;
	uint64_t Curr_Time, Prev_Time, Elap_Time;
	int Curr_ient, Prev_ient, Elap_ient;
	long CurrByte, EndByte;
	int pqrj_flag[4]={0,0,0,0};
	int event_flag=0;
	int fread_flag[Nthread];
	NHitCh=0;
	int nthread=0;
	int ithread;
	EventNumber=100;

	thread thrd[Nthread];
	Pulse1024* pulseBuffer[Nthread];
	Pulse1024* pulse[NOpCh];
	unsigned int threadnumber[NOpCh];

	int leakage=0;

	int ient=0;
	for (pqrj_flag[0]=0; ; ){
		Curr_ient = ient;	Elap_ient = Curr_ient - Prev_ient;
		Curr_Time = get_time();	Elap_Time = Curr_Time - Prev_Time; 

		KeyboardCommand(pqrj_flag);
		if(pqrj_flag[1]) event_flag=1;
		//if(pqrj_flag[1]) break;
		else if(pqrj_flag[0]){sleep(1); continue;}


		nthread=0;
		for (ithread=0; ithread<Nthread; ithread++){
			CurrByte = ftell(datafile);
			fseek(datafile,0,SEEK_END);
			EndByte = ftell(datafile);
			numevent = EndByte/sample_byte;
			fseek(datafile,CurrByte,SEEK_SET);
			if(EndByte-CurrByte>=sample_byte){
				fread_flag[ithread]=1;
				Nheader=(int)fread(header, 4, 6, datafile);
				Nwaveform=(int) fread(waveformBuffer, 4, waveform_length, datafile);
				//printf("EventSize %u, BoardID %u, StartIndexCellch %u, Channel %u, EventNumberch %u, TrigTimeTagch %u \n ", header[0],header[1],header[2],header[3],header[4],header[5]);
				if(header[0]!=4*(Nheader+Nwaveform)){
					printf("\nData Size Error!!  ient: %d, header[0]: %d, Nheader+Nwaveform: %d\n", ient, header[0], Nheader+Nwaveform);
					break;
				}
				pulseBuffer[ithread] = new Pulse1024(header,waveformBuffer);	leakage++;
				//thrd[nthread] = thread(&Pulse1024::Analyze,pulseBuffer[ithread],0,200, 200,950, 10, 0.15, samplingrate, 200,600);
				pulseBuffer[ithread]->Analyze(0,200, 200,950, 10, 0.15, samplingrate, 200,600);
				nthread++;
			}
			else if(EndByte-CurrByte<sample_byte){
				printf("(ONLINE)");
				fflush(stdout);
				fread_flag[ithread]=0;
				continue;
			}
		}

		//printf("thread filled, nthread %d\n", nthread);


		//for (ithread=0; ithread<nthread; ithread++){
		//thrd[ithread].join();
		//}	

		for (ithread=0; ithread<nthread; ithread++){
			if(pulseBuffer[ithread]->EventNumber > EventNumber && NHitCh!=0){
				event_flag=1;
				for (int ihit=0; ihit<NHitCh-1; ihit++){
					if(pulse[ihit]->EventNumber!=pulse[ihit+1]->EventNumber){event_flag=-1; break;}
					if(pulse[ihit]->TrigTimeTag!=pulse[ihit+1]->TrigTimeTag){event_flag=-1; break;}
					if(pulse[ihit]->StartIndexCell!=pulse[ihit+1]->StartIndexCell){event_flag=-1;break;}
				}
			}
			if(event_flag==-1){
				printf("StartIndexCell missmatch!\n", stderr);
				return;
			}
			//printf("event_flag %u @ithread %d\n", event_flag, ithread);

			if(event_flag && NHitCh!=0){
				EventNumber = pulse[0]->EventNumber; TrigTimeTag = pulse[0]->TrigTimeTag;	StartIndexCell=pulse[0]->StartIndexCell;	EventSize=pulse[0]->EventSize;
				for (int i=0; i<NHitCh; i++){
					rPH[i]=pulse[i]->rPH;
					QDC[i]=pulse[i]->QDC;
					TDC[i]=pulse[i]->TDC;
					Channel[i]=pulse[i]->Channel;
					delete pulse[i]; leakage--;
				}
				for (int i=0; i<NHitCh; i++){
					if(rPH[i]>3850) continue;
					h1_primary->rPH[Channel[i]]->Fill(rPH[i]);
					h1_primary->QDC[Channel[i]]->Fill(QDC[i]);
					h1_primary->TDC[Channel[i]]->Fill(TDC[i]);
				}
				for (int i=0; i<NHitCh; i++){
					if(rPH[i]>3850) continue;
					for (int j=0; j<NHitCh; j++){
						if(rPH[j]>3850) continue;
						h1_QDC_ch[Channel[i]]->coin[Channel[j]]->Fill(QDC[i]);
						h2_QDC_xaxis[Channel[i]]->yaxis[Channel[j]]->Fill(QDC[i],QDC[j]);
						if(QDCcut[Channel[i]][0]-2.0*QDCcut[Channel[i]][1]<QDC[i] && QDC[i]<QDCcut[Channel[i]][0]+2.0*QDCcut[Channel[i]][1]){
							if(QDCcut[Channel[j]][0]-2.0*QDCcut[Channel[j]][1]<QDC[j] && QDC[j]<QDCcut[Channel[j]][0]+2.0*QDCcut[Channel[j]][1]){
								if(TDC[i]<0) continue;
								if(TDC[j]<0) continue;
								h1_TDCdif[Channel[i]]->startch[Channel[j]]->Fill(TDC[i]-TDC[j]);
							}
						}
					}
				}	
				tree->Fill();

				NHitCh=0;
				event_flag=0;
			}
			EventNumber=pulseBuffer[ithread]->EventNumber;

			if(fread_flag[ithread] && pulseBuffer[ithread]->Channel<NOpCh){
				//thrd[ithread].join();
				threadnumber[NHitCh]=ithread;
				//printf("thrd[%u].join();\n",ithread);
				pulse[NHitCh]=pulseBuffer[ithread];
				//printf("pulse[NHitCh=%d]=pulseBuffer[ithread=%d];\n", NHitCh, ithread);

				event_flag=0;
				NHitCh++;
			}
			if(pulseBuffer[ithread]->Channel>=NOpCh){
				//thrd[ithread].detach();
				delete pulseBuffer[ithread];	leakage--;
			}
			ient++;

		}

		if(pqrj_flag[1]){
			printf("\nquit decoding\n");
			break;
		}

		ient--;
		if(Elap_Time > 1000){
			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s leakage %d", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time, leakage);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			if(pqrj_flag[2] || pqrj_flag[3]){
				canvases->Modified_Update(0x0F);
			}
			if(pqrj_flag[3]){
				canvases->Modified_Update(0x0F);

				printf("\n");
				//c1->Print(Form("jpg/%s_Decoding_c1.jpg",filename));
				//c2->Print(Form("jpg/%s_Decoding_c2.jpg",filename));
				pqrj_flag[3]=0;
			}
		}
		//if(!fread_flag[ithread]){
		//	sleep(1);
		//	continue;
		//}
		ient++;
	}

	printf("\nReading completed.Entries: %d, Event numbers: %u\n", ient, EventNumber);

	if(pqrj_flag[1]==2)	printf("Force Quit");
	else	canvases->Modified_Update(0x0F);

/*
	TF1 *Na662[NOpCh];
	double fitpar[NOpCh][5] ={	{400, 84000,1000,0,0},
					{300, 83700,1300,0,0},
					{300, 83700,1300,0,0},
					{300, 83700,1300,0,0}};
	for (int ich=0; ich<NOpCh; ich++){
		printf("\nch%d\n",ich);
		c1->cd(ich+1+1*NOpCh);
		Na662[ich] = new TF1(Form("Na662_ch%d",ich),"gaus(0)+pol1(3)", fitpar[ich][1]-3*fitpar[ich][2],fitpar[ich][1]+3*fitpar[ich][2]);
		Na662[ich]->SetParameters(fitpar[ich]);
		//Na662[ich]->SetParLimits(1,100000,130000);
		//Na662[ich]->SetParLimits(2,1000,3000);
		h1_QDC[ich]->Fit(Na662[ich],"M+","",fitpar[ich][1]-3*fitpar[ich][2],fitpar[ich][1]+3*fitpar[ich][2]);
		h1_QDC[ich]->Fit(Na662[ich],"M+","",Na662[ich]->GetParameter(1)-2*Na662[ich]->GetParameter(2),Na662[ich]->GetParameter(1)+2*Na662[ich]->GetParameter(2));
	}
	for (int ich=0; ich<NOpCh; ich++){
		printf("Eres ch%d=%.2f\n", ich, 235*Na662[ich]->GetParameter(2)/Na662[ich]->GetParameter(1));
	}
*/
	tree->Write();
//	c1->Write();
//	c2->Write();
//	treefile->Close();

	fclose(datafile);

	printf("\n");
}
