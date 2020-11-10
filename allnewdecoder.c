#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "./include/keyb.c"
#include "./include/gettimef.c"
#include "./include/Pulse1024.h"
#include "./include/DecoderCanvas.h"


void KeyboardCommand(int pqrjd_flag[5]){		
	if(!kbhit())	return;
	int key = getch();
	switch(key){
		case 'p' :						// 'p'ause
			pqrjd_flag[0]^=1; 
			if( pqrjd_flag[0]) printf("\nPaused. Press p again to resume\n");
			if(!pqrjd_flag[0]) printf("resumed\n");
			break;
		case 'q' : 
			pqrjd_flag[1]^=1; break;			// 'q'uit
		case 'Q' : 
			pqrjd_flag[1]=2; break;			// 'q'uit
		case 'r' : 						// 'r'enew histograms
			pqrjd_flag[2]^=1; 
			if( pqrjd_flag[2]) printf("\nRenewing histograms is enabled.\n");
			if(!pqrjd_flag[2]) printf("\nRenewing histograms is disabled.\n");
			break;
		case 'j' : 
			pqrjd_flag[3]^=1; break;			// saving canvas to filename_canvasname.'j'pg
		case 'd' :
			pqrjd_flag[4]=1;
			printf("\ndisplaying waveform\n"); break;               // pause and 'd'isplay current waveform
		default : break;
	}
}
//////////////////////main function////////////////////////////////

void allnewdecoder(){
	char filename[90]="grouptime_test";
	//char filename[90]="PMT_set1_12ch_Na22_run02_tres";
	//char filename[90]="PMT_set1_12ch_Eu152_run01_efftest";
	//char filename[90]="PMT_set2_12ch_Eu152_run01_efftest";
	FILE* datafile = fopen(Form("%s.dat",filename),"rb");	// moved from latter to here
	if(datafile==NULL){
		fputs("File error\n",stderr);
		exit(1);
	}

//tree
	const UInt_t NOpGr = 2;	// number of open groups
	const UInt_t NOpCh = 6;	// number of open channels
	const UInt_t waveform_length = 1024;
	const Float_t samplingrate = 2.5;  // GS/s
	const Float_t sampletime = 1.0/samplingrate;  // ns
	float waveformtime[waveform_length]; for (int iwft=0; iwft<waveform_length; iwft++) waveformtime[iwft]=sampletime*(float)iwft;
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);

	const int sample_byte = 4*(waveform_length+6);
	fseek(datafile,0,SEEK_END);
	ULong64_t Filesize = ftell(datafile);
	rewind(datafile);
	uint32_t numevent = Filesize/sample_byte;
	uint32_t remainder = Filesize%sample_byte;
	printf("%s.dat Filesize: %lld(Bytes), numevent: %u, Error:%u\n", filename, Filesize, numevent, remainder);

	//Float_t waveformBuffer[waveform_length];

	TFile* treefile = new TFile(Form("../../DecodedDataFiles/%s.root",filename),"recreate");
	TTree* tree = new TTree("wavedata","wavedata");
	UInt_t EventSize;
	//UInt_t BoardID=0;
	UInt_t NHitCh = NOpCh*NOpGr;	//number of hit channels, maximum NOpCh*NOpGr
	UInt_t StartIndexCell[NHitCh];
	//UInt_t Group;
	UInt_t Group[NHitCh];
	UInt_t Channel[NHitCh];
	UInt_t EventNumber;
	UInt_t TrigTimeTag;
	Float_t waveform[NHitCh][waveform_length];
	//Float_t waveformtime[waveform_length];
	Float_t rPH[NHitCh];
	Float_t QDC[NHitCh];
	Float_t TDC[NHitCh];

	tree->Branch("EventNumber", &EventNumber, "EventNumber/i");
	tree->Branch("TrigTimeTag", &TrigTimeTag, "TrigTimeTag/i");
	tree->Branch("EventSize", &EventSize, "EventSize/i");
	//tree->Branch("BoardID", &BoardID, "BoardID/i");
	//tree->Branch("StartIndexCell", &StartIndexCell, "StartIndexCell/i");
	//tree->Branch("Group", &Group, "Group/i");
	tree->Branch("NHitCh", &NHitCh, "NHitCh/i");
	tree->Branch("StartIndexCell", StartIndexCell, "StartIndexCell[NHitCh]/i");
	tree->Branch("Channel", Channel, "Channel[NHitCh]/i");
	tree->Branch("Group", Group, "Group[NHitCh]/i");
//	tree->Branch("waveform_length", &waveform_length, "waveform_length/i");
//	tree->Branch("waveform", waveform,"waveform[NHitCh][waveform_length]/F");
//	tree->Branch("waveformtime", waveformtime,"waveformtime[waveform_length]/F");
	tree->Branch("rPH", rPH, "rPH[NHitCh]/F");
	tree->Branch("QDC", QDC, "QDC[NHitCh]/F");
	tree->Branch("TDC", TDC, "TDC[NHitCh]/F");


//canvas
	const Int_t QDCrange = 300000/2.5*samplingrate;
	Double_t QDCcut[NOpCh][2] = {{8.40333e+04,1.35843e+05},{8.48711e+04,1.22737e+05},{8.48711e+04,1.22737e+05},{8.48711e+04,1.22737e+05},{8.48711e+04,1.22737e+05},{8.48711e+04,1.22737e+05}};

	DECODER_CANVAS *canvases = new DECODER_CANVAS(NOpGr, NOpCh);
	canvases->Construct();	

	DECODER_PRIMARY_H1 *h1_primary[NOpGr];
	for (int igr=0; igr<NOpGr; igr++){
		h1_primary[igr] = new DECODER_PRIMARY_H1(igr, NOpCh, QDCrange, waveform_length, sampletime);
		h1_primary[igr]->Draw(canvases);
	}
	//DECODER_QDC_H1_COIN *h1_QDC_ch[NOpCh];
	//DECODER_QDC_H2 *h2_QDC_xaxis[NOpCh];
	//DECODER_TDCDIF_H1 *h1_TDCdif[NOpCh];
	for (int ich=0; ich<NOpCh; ich++){
		//h1_QDC_ch[ich] = new DECODER_QDC_H1_COIN(NOpCh, ich, QDCrange);
		//h1_QDC_ch[ich]->Draw(canvases);
		//h2_QDC_xaxis[ich] = new DECODER_QDC_H2(NOpCh, ich, QDCrange);
		//h2_QDC_xaxis[ich]->Draw(canvases);
		//h1_TDCdif[ich] = new DECODER_TDCDIF_H1(NOpCh, ich, -3, 3);
		//h1_TDCdif[ich]->Draw(canvases);
	}
	TCanvas *c_waveform = new TCanvas("waveform", "waveform", 800, 600);
//	TGraph *gr_waveform;

//////////////////////

	uint32_t header[6];
	int Nheader, Nwaveform;
	uint64_t Curr_Time, Prev_Time, Elap_Time;
	uint64_t Curr_ient, Prev_ient, Elap_ient;
	uint64_t CurrByte, EndByte;
	//float decoding_speed; uint64_t time_left;
	int pqrjd_flag[5]={0,0,0,0,0};
	int event_flag=0;
	int fread_flag;
	NHitCh=0;

	CurrByte = ftell(datafile);
	Nheader=(int)fread(header, 4, 6, datafile);
	EventNumber=header[4];		//Get first event number.
	fseek(datafile,CurrByte,SEEK_SET);	//return to begining of the file

	Pulse1024* pulseBuffer;
	Pulse1024* pulse[NOpCh*NOpGr];
	int leakage=0;

	uint64_t ient=0;
	for (pqrjd_flag[0]=0; ; ){
		Curr_ient = ient;	Elap_ient = Curr_ient - Prev_ient;
		Curr_Time = get_time();	Elap_Time = Curr_Time - Prev_Time; 

		KeyboardCommand(pqrjd_flag);
		if(pqrjd_flag[1]) event_flag=1;
		//if(pqrjd_flag[1]) break;
		else if(pqrjd_flag[0]){sleep(1); continue;}


		CurrByte = ftell(datafile);
		fseek(datafile,0,SEEK_END);
		EndByte = ftell(datafile);
		numevent = EndByte/sample_byte;
		fseek(datafile,CurrByte,SEEK_SET);
		if(EndByte-CurrByte>=sample_byte){
			fread_flag=1;
			Nheader=(int)fread(header, 4, 6, datafile);
			pulseBuffer = new Pulse1024(header);	leakage++;
			Nwaveform=(int) fread(pulseBuffer->waveform, 4, waveform_length, datafile);
			//Nwaveform=(int) fread(waveformBuffer, 4, waveform_length, datafile);
			//printf("EventSize %u, Group %u, StartIndexCellch %u, Channel %u, EventNumberch %u, TrigTimeTagch %u \n ", header[0],header[1],header[2],header[3],header[4],header[5]);
			if(header[0]!=4*(Nheader+Nwaveform)){
				printf("\nData Size Error!!  ient: %lu, header[0]: %d, Nheader+Nwaveform: %d\n", ient, header[0], Nheader+Nwaveform);
				break;
			}
			if(pqrjd_flag[4]){
				if(c_waveform==NULL) TCanvas *c_waveform = new TCanvas("waveform", "waveform", 800, 600);
				c_waveform->cd();
				//if(gr_waveform!=NULL) gr_waveform->~TGraph();
				//gr_waveform = new TGraph(waveform_length, waveformtime, pulseBuffer->waveform);
				TGraph *gr_waveform = new TGraph(waveform_length, waveformtime, pulseBuffer->waveform);
				gr_waveform->SetName("waveform");
				gr_waveform->SetTitle(Form("waveform %d gr%dch%d", pulseBuffer->EventNumber, pulseBuffer->Group, pulseBuffer->Channel));
				gr_waveform->Draw("AP");
				c_waveform->WaitPrimitive();

				delete gr_waveform;
				//gr_waveform->~TGraph();
				pqrjd_flag[4]=0;
			}
			ient++;
		}
		else if(EndByte-CurrByte<sample_byte){
			printf("(ONLINE)");
			fflush(stdout);
			fread_flag=0;
		}

		if(fread_flag==1 && pulseBuffer->EventNumber!=EventNumber && NHitCh!=0){
			event_flag=1;
			for (int ihit=0; ihit<NHitCh-1; ihit++){
				if(NHitCh==1) break;
				if(pulse[ihit]->EventNumber!=pulse[ihit+1]->EventNumber){event_flag=-1; break;}
				if(pulse[ihit]->TrigTimeTag!=pulse[ihit+1]->TrigTimeTag){event_flag=-2; break;}
				//if(pulse[ihit]->StartIndexCell!=pulse[ihit+1]->StartIndexCell){event_flag=-3;break;}
			}
			if(event_flag==-1){
				printf("EventNumber missmatch @CurrByte %lu!\n", CurrByte, stderr);
				return;
			}
			if(event_flag==-2){
				printf("TrigTimeTag missmatch!\n", stderr);
				return;
			}
			if(event_flag==-3){
				printf("StartIndexCell missmatch!\n", stderr);
				return;
			}
		}

		if(event_flag==1){
			EventNumber = pulse[0]->EventNumber; TrigTimeTag = pulse[0]->TrigTimeTag;	EventSize=pulse[0]->EventSize;
			//EventNumber = pulse[0]->EventNumber; TrigTimeTag = pulse[0]->TrigTimeTag;	StartIndexCell=pulse[0]->StartIndexCell;	EventSize=pulse[0]->EventSize;
			for (int i=0; i<NHitCh; i++){
				rPH[i]=pulse[i]->rPH;
				QDC[i]=pulse[i]->QDC;
				TDC[i]=pulse[i]->TDC;
				Group[i]=pulse[i]->Group;
				Channel[i]=pulse[i]->Channel;
				StartIndexCell[i]=pulse[i]->StartIndexCell;
				delete pulse[i]; leakage--;
			}
			for (int i=0; i<NHitCh; i++){
				if(rPH[i]>3850) continue;
				if(QDC[i]< 1000) continue;
				h1_primary[Group[i]]->rPH[Channel[i]]->Fill(rPH[i]);
				h1_primary[Group[i]]->QDC[Channel[i]]->Fill(QDC[i]);
				h1_primary[Group[i]]->TDC[Channel[i]]->Fill(TDC[i]);
			}
/*			for (int i=0; i<NHitCh; i++){
				if(Group[i]!=0) continue;
				//if(rPH[i]>3850) continue;
				for (int j=0; j<NHitCh; j++){
					if(Group[j]!=0) continue;
					//if(rPH[j]>3850) continue;
					//h1_QDC_ch[Channel[i]]->coin[Channel[j]]->Fill(QDC[i]);
					//h2_QDC_xaxis[Channel[i]]->yaxis[Channel[j]]->Fill(QDC[i],QDC[j]);
					if(QDCcut[Channel[i]][0]-2.0*QDCcut[Channel[i]][1]<QDC[i] && QDC[i]<QDCcut[Channel[i]][0]+2.0*QDCcut[Channel[i]][1]){
						if(QDCcut[Channel[j]][0]-2.0*QDCcut[Channel[j]][1]<QDC[j] && QDC[j]<QDCcut[Channel[j]][0]+2.0*QDCcut[Channel[j]][1]){
							if(TDC[i]<0) continue;
							if(TDC[j]<0) continue;
							//h1_TDCdif[Channel[i]]->startch[Channel[j]]->Fill(TDC[i]-TDC[j]);
						}
					}
				}
			}	
*/			tree->Fill();

			NHitCh=0;
			event_flag=0;
		}

		if(fread_flag && pulseBuffer->Channel<NOpCh){
			EventNumber=pulseBuffer->EventNumber;
			pulseBuffer->Analyze(0,200, 200,950, 10, 0.15, samplingrate, 200,600);
			pulse[NHitCh]=pulseBuffer;

			//event_flag=0;
			NHitCh++;
		}
		if(fread_flag && pulseBuffer->Channel>=NOpCh){
			delete pulseBuffer;	leakage--;
		}


		if(pqrjd_flag[1]){
			printf("\nquit decoding\n");
			break;
		}

		if(Elap_Time > 1000){
			//decoding_speed=(float)Elap_ient*1000.0f/(float)Elap_Time;
			//time_left=(numevent-ient)/(int)decoding_speed;
			//printf("\rdecoding... at %lu of %u (%.3f) %.0fient/s ETF %lusec leakage %d", ient, numevent, (float) ient/numevent, decoding_speed, time_left, leakage);
			printf("\rdecoding... at %lu of %u (%.3f) %.0fient/s leakage %d", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time, leakage);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			if(pqrjd_flag[2] || pqrjd_flag[3]){
				canvases->Modified_Update(0x01);
				//canvases->Modified_Update(0x0F);
			}
			if(pqrjd_flag[3]){
				canvases->Modified_Update(0x0F);
				printf("\n");
				canvases->Print(filename);
				pqrjd_flag[3]=0;
			}
		}
		if(!fread_flag){
			sleep(1);
			continue;
		}
	}

	printf("\nReading completed.Entries: %lu, Event numbers: %u\n", ient, EventNumber);

	if(pqrjd_flag[1]==2)	printf("Force Quit");
	else	canvases->Modified_Update(0x01);
	//else	canvases->Modified_Update(0x0F);


/*
	TF1 *Na662[NOpCh];
	double fitpar[NOpCh][5] ={	{400, 84000,1000,0,0},
					{300, 83700,1300,0,0},
					{300, 83700,1300,0,0},
					{300, 83700,1300,0,0}};
	for (int ich=0; ich<NOpCh; ich++){
		printf("\nch%d\n",ich);
		canvases->Canvas_Decoding->cd(ich+1+1*NOpCh);
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
	treefile->Write();
//	canvases->Canvas_Decoding->Write();
//	for (int i=0; i<3; i++)		canvases->Coincidence[i]->Write();
			///or, maybe I can make directory in TFile, where I save all histograms
//	treefile->Close();

	fclose(datafile);

	printf("\n");
}
