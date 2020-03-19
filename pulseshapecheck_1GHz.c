#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "./include/keyb.c"
#include "./include/gettimef.c"
#include "./include/DRS4_2_pulse_analysis_functions.c"

void Delayf(Float_t waveform[520],Int_t delay, Float_t waveform2[520]){
	for (int i=0; i<520-delay; i++){
		waveform2[i]= waveform[i+delay];
	}
	for (int i=520-delay; i<520; i++){
		waveform2[i]=0;
	}
}

void Compressf(Float_t waveform[520], double fraction, Float_t waveform2[520]){
	double pedestal=(double)pedsum(waveform)/160;
	for (int i=0; i<520; i++){
		waveform2[i]=pedestal-fraction*(pedestal-waveform[i]);
	}
}

void CFDf(Float_t waveform[520],Int_t delay, double fraction, Float_t waveform2[520]){
	double pedestal=(double)pedsum(waveform)/160;
	//Double_t waveform2[520-delay];
	for (int i=0; i<520-delay; i++){
		waveform2[i]=(double) -fraction*(pedestal-waveform[i+delay])+(pedestal-waveform[i])+0.0;
	}
	for (int i=520-delay; i<520; i++){
		waveform2[i]=0.0;
	}	
}

double pulsesum(Float_t waveform[520]){
	Float_t pulsesum=0;
	for (int iwft=0; iwft< 480;iwft++){
		pulsesum = pulsesum + waveform[iwft];
	}
	return pulsesum;
}


void waveformcopy(Float_t arrdest[520], Float_t arrsour[520]){
	for (int iwft=0; iwft<520; iwft++){
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

void pulseshapecheck_1GHz(){
	char filename[90]="PMT_9A_primary_Co60_run01_HV_1040_1050_Noatt_1GHz";
	//char filename[90]="PMT_9A_primary_Na22_run02_HV_1040_1050_Noatt_1GHz";
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

	//TFile *treefile = new TFile(Form("../../DecodedDataFiles/%s.root",filename),"recreate");
	//TString str0 = TString("Set_info\n");
	//TString str1 = TString("Co60, attenuator \n");
	//TString str2 = TString(" ch0: Scint-38*38LaBr3 PMT-R13408(Roka numbering1) HV-1500V from HVch0,att6dB\n");
	//TString str3 = TString(" ch1: Scint-38*38LaBr3 PMT-R13408(Roka numvering2) HV-1500V from HVch1,att6dB\n");
	//TString str2 = TString(" ch0: Scint-38*38LaBr3(20170809-1) PMT-R13408(AA0037) HV-1500V from HVch0,att6dB\n");
	//TString str3 = TString(" ch1: Scint-38*38LaBr3(20170809-0) PMT-R13408(AA0038) HV-1500V from HVch1,att6dB\n");
	//TObjString *TObjS = new TObjString(str0+str1+str2+str3);
	//TObjS->Write();


	UInt_t NOpCh = 2;	// number of open channels
	UInt_t OpCh;
	UInt_t waveform_length = 520;
	Float_t waveformBuffer[waveform_length];
	Float_t waveformch[NOpCh][waveform_length];
	Float_t waveform2ch[NOpCh][waveform_length];
	Float_t waveform3ch[NOpCh][waveform_length];
	Float_t waveform4ch[NOpCh][waveform_length];
	UInt_t StartIndexCellch[NOpCh];
	UInt_t EventNumberch[NOpCh];
	UInt_t TrigTimeTagch[NOpCh];
	Float_t rPHch[NOpCh];
	Float_t QDCch[NOpCh];
	Float_t TDCch[NOpCh];

//	TTree *tree = new TTree("wavedata","wavedata");
	UInt_t EventSize;
	UInt_t BoardID;
	UInt_t StartIndexCell;
	//UInt_t Group;
	UInt_t NHitCh = NOpCh;	//number of hit channels, maximum NOpCh
	UInt_t Channel[NHitCh];
	UInt_t EventNumber;
	UInt_t TrigTimeTag;
	Float_t waveform[NOpCh][waveform_length];
	Float_t waveformtime[waveform_length];
	Float_t rPH[NHitCh];
	Float_t QDC[NHitCh];
	Float_t TDC[NHitCh];

/*	tree->Branch("EventNumber", &EventNumber, "EventNumber/i");
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
*/
	Float_t samplingrate = 1.0;  // GS/s
	Float_t sampletime = 1.0/samplingrate;  // ns
	for (int iwft=0; iwft<waveform_length; iwft++){
		waveformtime[iwft] = sampletime * iwft;
	}
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);

	const int delay=3;
	const float fraction=0.20;

	const int sample_byte = 4*(waveform_length+6);
	fseek(datafile,0,SEEK_END);
	ULong64_t Filesize = ftell(datafile);
	rewind(datafile);
//	cout<<"Filesize = "<<Filesize<<endl;
	uint32_t numevent = Filesize/sample_byte;
	uint32_t remainder = Filesize%sample_byte;
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
	int dummy;

	float pulseavg[NOpCh];
	float pedestal[NOpCh];
	TLine *pedline[NOpCh], *l_tdc1[NOpCh], *l_tdc2[NOpCh];
	TGraph *gr1[NOpCh], *gr2[NOpCh], *gr3[NOpCh];
	TMultiGraph *mg[NOpCh];

	TCanvas *c1 = new TCanvas("pulseshape","pulseshape", 1200,800);
	c1->Divide(NOpCh,2);
	c1->cd(1)->SetGridy();
	c1->cd(2)->SetGridy();
	c1->cd(3)->SetGridy();
	c1->cd(4)->SetGridy();

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
		numevent = EndByte/sample_byte;
		fseek(datafile,CurrByte,SEEK_SET);
		if(EndByte-CurrByte>=sample_byte){
			fread_flag=1;
			Nheader=(int)fread(header, 4, 6, datafile);
			Nwaveform=(int) fread(waveformBuffer, 4, waveform_length, datafile);
			printf("EventSize %u, BoardID %u, StartIndexCellch %u, Channel %u, EventNumberch %u, TrigTimeTagch %u \n ", header[0],header[1],header[2],header[3],header[4],header[5]);
			if(header[0]!=4*(Nheader+Nwaveform)){
				printf("\nData Size Error!!  ient: %d, header[0]: %d, Nheader+Nwaveform: %d\n", ient, header[0], Nheader+Nwaveform);
				break;
			}
		}
		else if(EndByte-CurrByte<sample_byte){
			printf("(ONLINE)");
			fflush(stdout);
			fread_flag=0;
		}

		if(header[4]-EventNumber>0 && NHitCh!=0){
			event_flag=1;
			for (int ich=0; ich<NHitCh-1; ich++){
				if(EventNumberch[ich]!=EventNumberch[(ich+1)]){event_flag=-1; break;}
				if(TrigTimeTagch[ich]!=TrigTimeTagch[(ich+1)]){event_flag=-1; break;}
				if(StartIndexCellch[ich]!=StartIndexCellch[(ich+1)]){event_flag=-1;break;}
			}
		}
		if(event_flag==-1){
			printf("StartIndexCell missmatch!\n", stderr);
		}
		if(event_flag && NHitCh!=0){
			EventNumber = EventNumberch[0]; TrigTimeTag = TrigTimeTagch[0];	StartIndexCell=StartIndexCellch[0];
			for (int i=0; i<NHitCh; i++){
				rPH[i] = pulseheight(waveformch[i]);
				QDC[i] = QDCf(waveformch[i]);
				TDC[i] = TDCcf(waveformch[i], waveformtime, delay, fraction);
				if(NHitCh!=2) continue;
				//if(QDC[i]<3.42934e+04-3*5.27601e+02 || 3.42934e+04+3*5.27601e+02<QDC[i]) continue;
				printf("ient %d, i: %d, rPH[i]=%.1f, QDC[i]=%.1f, TDC[i]=%.1f, Channel[i] %d\n",ient, i,rPH[i],QDC[i], TDC[i], Channel[i]);
				pqrj_flag[0]=1;
				pedestal[i] = pedsum(waveformch[i])/160;
				pulseavg[i] = pulsesum(waveformch[i])/ 480;
				//c1->Clear();
				//c1->cd(1)->Clear();
				gr1[i] = new TGraph( 520,waveformtime,waveformch[i]);
				gr1[i]->SetTitle(Form("waveform ient %d, i: %d;time(ns);ADC", ient, i));
				gr1[i]->GetXaxis()->SetLimits(160,300);
				gr1[i]->GetYaxis()->SetRangeUser(0,4096);
				//gr1->Draw();

				Compressf(waveformch[i], fraction, waveform2ch[i]);
				Delayf(waveform2ch[i],  delay, waveform3ch[i]);
				gr2[i] = new TGraph( 520- delay,waveformtime,waveform3ch[i]);
				gr2[i]->SetLineColor(2);
				gr2[i]->SetTitle(Form("waveform ient %d, i: %d;time(ns);ADC", ient, i));
				//gr2->Draw("SAME");

				c1->cd(1+2*i);
				//c1->cd(1)->Clear();
				mg[i] = new TMultiGraph("concept",Form("waveform ient %d, i: %d;time(ns);ADC", ient, i));
				mg[i]->Add(gr1[i]); mg[i]->Add(gr2[i]);
				mg[i]->GetXaxis()->SetLimits(200,300);
				mg[i]->GetYaxis()->SetRangeUser(0,4096);
				mg[i]->Draw("APL");

				pedline[i] = new TLine(0,pedestal[i],450,pedestal[i]);
				pedline[i]->SetLineColor(3);
				pedline[i]->Draw("SAME");

				l_tdc1[i] = new TLine(TDC[i], 0, TDC[i], 4096);
				l_tdc1[i]->Draw("SAME");

				c1->cd(2+2*i)->Clear();
				CFDf(waveformch[i], delay, fraction, waveform4ch[i]);
				gr3[i] = new TGraph(520, waveformtime, waveform4ch[i]);
				gr3[i]->SetLineColor(4);
				gr3[i]->SetMarkerStyle(8);
				gr3[i]->GetXaxis()->SetLimits(200,300);
				gr3[i]->SetTitle(Form("waveform ient %d, i: %d;time(ns);ADC", ient, i));
				gr3[i]->Draw();
				//gr3->Draw("SAME");

				l_tdc2[i] = new TLine(TDC[i], -500, TDC[i], 500);
				l_tdc2[i]->Draw("SAME");

				c1->cd(1);	gPad->Modified(); gPad->Update();
				c1->cd(2);	gPad->Modified(); gPad->Update();
				c1->cd(3);	gPad->Modified(); gPad->Update();
				c1->cd(4);	gPad->Modified(); gPad->Update();
				cin>>dummy;
				//cin>>pqrj_flag[0];
			}

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

			//printf("EventSize %u, BoardID %u, StartIndexCellch %u, Channel %u, EventNumberch %u, TrigTimeTagch %u \n ", header[0],header[1],header[2],header[3],header[4],header[5]);
			event_flag=0;
			NHitCh++;
		}

		if(pqrj_flag[1]) break;

		if(Elap_Time > 1000){
			//			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time);
			//			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			if(pqrj_flag[2] || pqrj_flag[3]){
				for (int jch=0; jch<NOpCh; jch++){
					//c1->cd(3*jch+1);	gPad->Modified();	gPad->Update();
					//c1->cd(3*jch+2);	gPad->Modified();	gPad->Update();
					//c1->cd(3*jch+3);	gPad->Modified();	gPad->Update();
					/*for(int kch=0; kch<NOpCh; kch++){
					  c2->cd(1)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
					  c2->cd(2)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
					  }*/
				}
			}
			if(pqrj_flag[3]){
				for (int jch=0; jch<NOpCh; jch++){
					for(int kch=0; kch<NOpCh; kch++){
						//c2->cd(1)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
						//c2->cd(2)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
					}
				}
				printf("\n");
				//c1->Print(Form("jpg/%s_Decoding_c1.jpg",filename));
				//c2->Print(Form("jpg/%s_Decoding_c2.jpg",filename));
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
		//c1->cd(3*jch+1);	gPad->Modified();	gPad->Update();
		//c1->cd(3*jch+2);	gPad->Modified();	gPad->Update();
		//c1->cd(3*jch+3);	gPad->Modified();	gPad->Update();
		for(int kch=0; kch<NOpCh; kch++){
			//c2->cd(1)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
			//c2->cd(2)->cd(1+jch+kch*NOpCh);	gPad->Modified();	gPad->Update();
		}
	}
	printf("\nReading completed.Entries: %d, Event numbers: %u\n", ient, EventNumber);



	fclose(datafile);

	printf("\n");
}
