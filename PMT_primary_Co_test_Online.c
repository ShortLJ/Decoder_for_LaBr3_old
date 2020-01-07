#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "./include/keyb.c"
#include "./include/gettimef.c"
#include "./include/pulse_analysis_functions.c"

void waveformcopy(Float_t arrdest[1024],const Float_t arrsour[1024]){
	for (int iwft=0; iwft<1024; iwft++){
		arrdest[iwft]=arrsour[iwft];
	}
}

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
/*
void NewCanvas(TCanvas *c1, TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4);
	c1 = new TCanvas("Decoding","Decoding",1200,800);
	c1->Divide(2,2);
	c1->cd(1);
//	h1_QDC[0]->GetXaxis()->SetRangeUser(5000,QDCrange);
	h1->Draw();
	c1->cd(2);
//	h1_QDC[1]->GetXaxis()->SetRangeUser(5000,QDCrange);
	h2->Draw();
	c1->cd(3);
	h3->Draw();
	c1->cd(4);
	h4->Draw();
}
*/

//////////////////////main function////////////////////////////////

void PMT_primary_Co_test_Online(){
	char filename[90]="PMT_set2_primary_Co60_run04_HV_1160_1040_1050_1110_NoAtt";
	//char filename[90]="PMT_primary_TDC_setvoltage_set2";
	//char filename[90]="PMT_primary_test_setvoltage_set1";
	//char filename[90]="Eu152_att66_newcalib_run02";
	//char filename[90]="NoSource_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run04";
	//char filename[90]="Na22_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run01";
	//char filename[90]="Co60_att66_newcalib_run02";
	FILE *datafile = fopen(Form("%s.dat",filename),"rb");	// moved from latter to here
	if(datafile==NULL){
		fputs("File error\n",stderr);
		exit(1);
	}

//tree

	TFile *treefile = new TFile(Form("../../DecodedDataFiles/%s.root",filename),"recreate");
	//TString str0 = TString("Set_info\n");
	//TString str1 = TString("Co60, attenuator \n");
	//TString str2 = TString(" ch0: Scint-38*38LaBr3 PMT-R13408(Roka numbering1) HV-1500V from HVch0,att6dB\n");
	//TString str3 = TString(" ch1: Scint-38*38LaBr3 PMT-R13408(Roka numvering2) HV-1500V from HVch1,att6dB\n");
	//TString str2 = TString(" ch0: Scint-38*38LaBr3(20170809-1) PMT-R13408(AA0037) HV-1500V from HVch0,att6dB\n");
	//TString str3 = TString(" ch1: Scint-38*38LaBr3(20170809-0) PMT-R13408(AA0038) HV-1500V from HVch1,att6dB\n");
	//TObjString *TObjS = new TObjString(str0+str1+str2+str3);
	//TObjS->Write();


	UInt_t NOpCh = 4;	// number of open channels
	UInt_t OpCh;
	UInt_t waveform_length = 1024;
	Float_t waveformBuffer[waveform_length];
	Float_t waveformch[NOpCh][waveform_length];
	UInt_t StartIndexCellch[NOpCh];
	UInt_t EventNumberch[NOpCh];
	UInt_t TrigTimeTagch[NOpCh];
	Float_t rPHch[NOpCh];
	Float_t QDCch[NOpCh];
	Float_t TDCch[NOpCh];

	TTree *tree = new TTree("wavedata","wavedata");
	UInt_t EventSize;
	UInt_t BoardID;
	UInt_t StartIndexCell;
	//UInt_t Group;
	UInt_t NHitCh = NOpCh;	//number of hit channels, maximum NOpCh
	UInt_t Channel[NHitCh];
	UInt_t EventNumber;
	UInt_t TrigTimeTag;
//	UInt_t waveform_length = 1024;	// moved to upper part with each channel
	Float_t waveform[NOpCh][waveform_length];
	Float_t waveformtime[waveform_length];
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
//	tree->Branch("waveform", waveform,"waveform[NHitCh][1024]/F");
//	tree->Branch("waveformtime", waveformtime,"waveformtime[waveform_length]/F");
	tree->Branch("rPH", rPH, "rPH[NHitCh]/F");
	tree->Branch("QDC", QDC, "QDC[NHitCh]/F");
	tree->Branch("TDC", TDC, "TDC[NHitCh]/F");

	Float_t samplingrate = 2.5;  // GS/s
	Float_t sampletime = 1.0/samplingrate;  // ns
	for (int iwft=0; iwft<waveform_length; iwft++){
		waveformtime[iwft] = sampletime * iwft;
	}
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);

//canvas1

	TCanvas *c1 = new TCanvas("Decoding","Decoding",1800, 1000);
	c1->Divide(NOpCh,3);
	TH1D *h1_rPH[NOpCh];
	TH1D *h1_QDC[NOpCh];
	TH1D *h1_TDC[NOpCh];

	TCanvas *c2 = new TCanvas("Coincidence", "Coincidence",1920, 954);
	c2->Divide(2,1);
	c2->cd(1)->Divide(NOpCh,NOpCh);
	c2->cd(2)->Divide(NOpCh,NOpCh);
	TH1D *h1_QDC_coin[NOpCh][NOpCh];
	TH2D *h2_QDC_coin[NOpCh][NOpCh];
	TH1D *h1_TDCdif_coin[NOpCh][NOpCh];

	Int_t QDCrange = 300000;
	Double_t QDCcut[4][2] = {{8.753975e+04,1.599353e+03},{8.040078e+04,1.289085e+03},{8.346844e+04,1.338430e+03},{8.090783e+04,1.374366e+03}};
	//Double_t QDCcut[NOpCh][2] = {{8.753975e+04,1.599353e+03},{8.040078e+04,1.289085e+03},{8.346844e+04,1.338430e+03},{8.090783e+04,1.374366e+03}};

	for (int ich=0; ich<NOpCh; ich++){
		h1_rPH[ich] = new TH1D(Form("h1_rPH[%d]",ich),Form("h1_rPH[%d];rPH;count/rPH",ich), 4096, 0, 4096);
		h1_QDC[ich] = new TH1D(Form("h1_QDC[%d]",ich),Form("h1_QDC[%d];QDC;count/150QDC",ich), 2000, 0, QDCrange);
		h1_TDC[ich] = new TH1D(Form("h1_TDC[%d]",ich),Form("h1_TDC[%d];TDC;count/0.4ns",ich), 1024, 0, 1024*0.4);
		c1->cd(ich+1+0*NOpCh)->SetGridx();	h1_rPH[ich]->Draw();
		c1->cd(ich+1+1*NOpCh)->SetGridx();	h1_QDC[ich]->Draw();
		c1->cd(ich+1+2*NOpCh)->SetGridx();	h1_TDC[ich]->Draw();
		for (int jch=0; jch<NOpCh; jch++){
			h1_QDC_coin[ich][jch] = new TH1D(Form("QDC_ch%d_coin%d", ich,jch),Form("QDC_ch%d_coin%d;QDC(keV);count per 0.5keV", ich,jch), 4000,0,300000);
			h2_QDC_coin[ich][jch] = new TH2D(Form("QDC_ch%d_QDC_ch%d", ich,jch),Form("QDC_ch%d_QDC_ch%d;QDC_ch%d(keV);QDC%d(keV)", ich,jch,ich,jch), 4000,0,300000, 4000,0,300000);
			h1_TDCdif_coin[ich][jch] = new TH1D(Form("TDCdif_ch%d-ch%d", ich,jch),Form("TDC_ch%d-ch%d;TDC dif(ns);count per 0.02ns", ich,jch), 600,-6,6);
			if(ich==jch) continue;
			c2->cd(1)->cd(1+ich+jch*NOpCh)->SetLogy();
			c2->cd(1)->cd(1+ich+jch*NOpCh)->SetGridx();
			h1_QDC_coin[ich][jch]->Draw();
		}
		for (int jch=0; jch<ich; jch++){
			c2->cd(2)->cd(1+ich+jch*NOpCh)->SetLogz();
			c2->cd(2)->cd(1+ich+jch*NOpCh)->SetGridx();
			c2->cd(2)->cd(1+ich+jch*NOpCh)->SetGridy();
			h2_QDC_coin[ich][jch]->Draw("COLZ");

			c2->cd(2)->cd(1+jch+ich*NOpCh)->SetLogy();
			c2->cd(2)->cd(1+jch+ich*NOpCh)->SetGridx();
			h1_TDCdif_coin[ich][jch]->Draw();
		}
	}


//file

//	FILE *datafile = fopen(Form("%s.dat",filename),"rb");	// moved to upper
//	if(datafile==NULL){
//		fputs("File error\n",stderr);
//		exit(1);
//	}
	fseek(datafile,0,SEEK_END);
	ULong64_t Filesize = ftell(datafile);
	rewind(datafile);
//	cout<<"Filesize = "<<Filesize<<endl;
	uint32_t numevent = Filesize/4120;
	uint32_t remainder = Filesize%4120;
	printf("%s.dat Filesize: %lld(Bytes), numevent: %u, Error:%u\n", filename, Filesize, numevent, remainder);

	uint32_t header[6];
	int Nheader, Nwaveform;
	uint64_t Curr_Time, Prev_Time, Elap_Time;
	int Curr_ient, Prev_ient, Elap_ient;
	long CurrByte, EndByte;
	int pqrj_flag[4]={0,0,1,0};
	int event_flag=0;
	int fread_flag=1;
	NHitCh=0;
	EventNumber=100;

	int ient=0;
//	for (pqrj_flag[0]=0; ient<numevent; ){
	for (pqrj_flag[0]=0; ; ){
		Curr_ient = ient;	Elap_ient = Curr_ient - Prev_ient;
		Curr_Time = get_time();	Elap_Time = Curr_Time - Prev_Time; 

		KeyboardCommand(pqrj_flag);
		if(pqrj_flag[1]) event_flag=1;
		//if(pqrj_flag[1]) break;
		else if(pqrj_flag[0]){sleep(1); continue;}

		CurrByte = ftell(datafile);
		fseek(datafile,0,SEEK_END);
		EndByte = ftell(datafile);
		numevent = EndByte/4120;
		fseek(datafile,CurrByte,SEEK_SET);
		if(EndByte-CurrByte>=4120){
			fread_flag=1;
			Nheader=(int)fread(header, 4, 6, datafile);
			Nwaveform=(int) fread(waveformBuffer, 4, 1024, datafile);
			//printf("EventSize %u, BoardID %u, StartIndexCellch %u, Channel %u, EventNumberch %u, TrigTimeTagch %u \n ", header[0],header[1],header[2],header[3],header[4],header[5]);
			if(header[0]!=4*(Nheader+Nwaveform)){
				printf("\nData Size Error!!  ient: %d, header[0]: %d, Nheader+Nwaveform: %d\n", ient, header[0], Nheader+Nwaveform);
				break;
			}
		}
		else if(EndByte-CurrByte<4120){
			printf("(ONLINE)");
			fflush(stdout);
			fread_flag=0;
		}

		if(header[4]-EventNumber>0 && NHitCh!=0){
			for (int ich=0; ich<NHitCh-1; ich++){
				if(EventNumberch[ich]!=EventNumberch[(ich+1)]){event_flag=-1; break;}
				if(TrigTimeTagch[ich]!=TrigTimeTagch[(ich+1)]){event_flag=-1; break;}
				if(StartIndexCellch[ich]!=StartIndexCellch[(ich+1)]){event_flag=-1;break;}
			}
			event_flag=1;
		}
		if(event_flag==-1){
			printf("StartIndexCell missmatch!\n", stderr);
		}
		if(event_flag && NHitCh!=0){
			EventNumber = EventNumberch[0]; TrigTimeTag = TrigTimeTagch[0];	StartIndexCell=StartIndexCellch[0];
			for (int i=0; i<NHitCh; i++){
				rPH[i] = pulseheight(waveformch[i]);
				QDC[i] = QDCf(waveformch[i]);
				TDC[i] = TDCcf(waveformch[i], waveformtime, 10, 0.15);
				if(rPH[i]>3850) continue;
				h1_rPH[Channel[i]]->Fill(rPH[i]);
				h1_QDC[Channel[i]]->Fill(QDC[i]);
				h1_TDC[Channel[i]]->Fill(TDC[i]);
				for (int j=0; j<i; j++){
					if(rPH[j]>3850) continue;
					h1_QDC_coin[Channel[i]][Channel[j]]->Fill(QDC[i]);
					h1_QDC_coin[Channel[j]][Channel[i]]->Fill(QDC[j]);
					h2_QDC_coin[Channel[i]][Channel[j]]->Fill(QDC[i],QDC[j]);
					if(QDCcut[Channel[i]][0]-1.5*QDCcut[Channel[i]][1]<QDC[i] && QDC[i]<QDCcut[Channel[i]][0]+1.5*QDCcut[Channel[i]][1] && QDCcut[Channel[j]][0]-1.5*QDCcut[Channel[j]][1]<QDC[j] && QDC[j]<QDCcut[Channel[j]][0]+1.5*QDCcut[Channel[j]][1])	h1_TDCdif_coin[Channel[i]][Channel[j]]->Fill(TDC[i]-TDC[j]);
				}
			}	
			tree->Fill();

			NHitCh=0;
			event_flag=0;
		}
		EventNumber=header[4];

		if(fread_flag && header[3]<NOpCh){
			EventSize=header[0];
			BoardID=header[1];
			//OpCh=header[3];
			Channel[NHitCh]=header[3];
			StartIndexCellch[NHitCh]=header[2];
			EventNumberch[NHitCh]=header[4];
			TrigTimeTagch[NHitCh]=header[5];
			waveformcopy(waveformch[NHitCh],waveformBuffer);
			waveform_length=Nwaveform;

			event_flag=0;
			NHitCh++;
		}

		if(pqrj_flag[1]) break;

		if(Elap_Time > 1000){
			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			if(pqrj_flag[2] || pqrj_flag[3]){
				for (int jch=0; jch<NOpCh; jch++){
					c1->cd(3*jch+1);	gPad->Modified();	gPad->Update();
					c1->cd(3*jch+2);	gPad->Modified();	gPad->Update();
					c1->cd(3*jch+3);	gPad->Modified();	gPad->Update();
					//for(int kch=0; kch<NOpCh; kch++){
						//c2->cd(1)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
						//c2->cd(2)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
					//}
				}
			}
			if(pqrj_flag[3]){
				printf("\n");
				c1->Print(Form("jpg/%s_Decoding_c1.jpg",filename));
				c2->Print(Form("jpg/%s_Decoding_c2.jpg",filename));
				pqrj_flag[3]=0;
			}
		}
		if(!fread_flag){
			sleep(1);
			continue;
		}
		ient++;
	}
	for (int jch=0; jch<NOpCh; jch++){
		c1->cd(3*jch+1);	gPad->Modified();	gPad->Update();
		c1->cd(3*jch+2);	gPad->Modified();	gPad->Update();
		c1->cd(3*jch+3);	gPad->Modified();	gPad->Update();
		for(int kch=0; kch<NOpCh; kch++){
			c2->cd(1)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
			c2->cd(2)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
		}
	}
	printf("\nReading completed.Entries: %d, Event numbers: %u\n", ient, EventNumber);

	TF1 *Na662[NOpCh];
	double fitpar[5] = {400,107000,1500,0,0};
	for (int ich=0; ich<NOpCh; ich++){
		printf("\nch%d\n",ich);
		Na662[ich] = new TF1(Form("Na662_ch%d",ich),"gaus(0)+pol1(3)", 100000,120000);
		Na662[ich]->SetParameters(fitpar);
		Na662[ich]->SetParLimits(1,100000,130000);
		Na662[ich]->SetParLimits(2,1000,3000);
		h1_QDC[ich]->Fit(Na662[ich],"M+","",100000,120000);
		h1_QDC[ich]->Fit(Na662[ich],"M+","",Na662[ich]->GetParameter(1)-2*Na662[ich]->GetParameter(2),Na662[ich]->GetParameter(1)+2*Na662[ich]->GetParameter(2));
	}
	for (int ich=0; ich<NOpCh; ich++){
		printf("Eres ch%d=%.2f\n", ich, 235*Na662[ich]->GetParameter(2)/Na662[ich]->GetParameter(1));
	}

	tree->Write();
	c1->Write();
	c2->Write();
//	treefile->Close();

	fclose(datafile);

	printf("\n");
}
