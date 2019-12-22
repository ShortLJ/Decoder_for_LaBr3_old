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

void zerosupdatareadingOnline(){
	char filename[90]="PMT_primary_TDC_setvoltage_set0";
	//char filename[90]="PMT_primary_test_set3_att";
	//char filename[90]="zerosuppressiontest";
	//char filename[90]="Eu152_att66_newcalib_run02";
	//char filename[90]="NoSource_att66_newcalib_run01";
	//char filename[90]="Na22_att66_newcalib_run04";
	//char filename[90]="Na22_att66_newcalib_run01";
	//char filename[90]="Na22_att66_newcalib_run01";
	//char filename[90]="Na22_att66_newcalib_run01";
	//char filename[90]="Na22_att66_newcalib_run02";
	FILE *datafile = fopen(Form("%s.dat",filename),"rb");	// moved from latter to here
	if(datafile==NULL){
		fputs("File error\n",stderr);
		exit(1);
	}

//tree

	TFile *treefile = new TFile(Form("../../DecodedDataFiles/%s.root",filename),"recreate");
	TString str0 = TString("Set_info\n");
	TString str1 = TString("Na22, attenuator \n");
//	TString str2 = TString(" ch0: Scint-38*38LaBr3 PMT-R13408(Roka numbering1) HV-1500V from HVch0,att6dB\n");
//	TString str3 = TString(" ch1: Scint-38*38LaBr3 PMT-R13408(Roka numvering2) HV-1500V from HVch1,att6dB\n");
//	TString str2 = TString(" ch0: Scint-38*38LaBr3(20170809-1) PMT-R13408(AA0037) HV-1500V from HVch0,att6dB\n");
//	TString str3 = TString(" ch1: Scint-38*38LaBr3(20170809-0) PMT-R13408(AA0038) HV-1500V from HVch1,att6dB\n");
	TObjString *TObjS = new TObjString(str0+str1+str2+str3);
	TObjS->Write();


	UInt_t NOpCh = 2;	// number of open channels
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
	UInt_t CoinEventNumber = 0;
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
	tree->Branch("CoinEventNumber", &CoinEventNumber, "CoinEventNumber/i");
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
	TH1D *h1_rPH[NOpCh];
	TH1D *h1_QDC[NOpCh];
	TH1D *h1_TDC[NOpCh];

	Int_t QDCrange = 300000;

	h1_rPH[0] = new TH1D("h1_rPH[0]","h1_rPH[0];rPH;count/rPH", 4096, 0, 4096);
	h1_QDC[0] = new TH1D("h1_QDC[0]","h1_QDC[0];QDC;count/150QDC", 2000, 0, QDCrange);
	h1_TDC[0] = new TH1D("h1_TDC[0]","h1_TDC[0];TDC (ns);count/0.4ns",1024,0,410);
	h1_rPH[1] = new TH1D("h1_rPH[1]","h1_rPH[1];rPH;count/rPH", 4096, 0, 4096);
	h1_QDC[1] = new TH1D("h1_QDC[1]","h1_QDC[1];QDC;count/150QDC", 2000, 0, QDCrange);
	h1_TDC[1] = new TH1D("h1_TDC[1]","h1_TDC[1];TDC (ns);count/0.4ns",1024,0,410);
//	NewCanvas(c1, h1_QDC[0], h1_QDC[1], h1_TDC[0], h1_TDC[1]);
	c1->Divide(3,2);
	c1->cd(1)->SetGridx();	h1_rPH[0]->Draw();
	c1->cd(2)->SetGridx();	h1_QDC[0]->Draw();
	c1->cd(3)->SetGridx();	h1_TDC[0]->Draw();
	c1->cd(4)->SetGridx();	h1_rPH[1]->Draw();
	c1->cd(5)->SetGridx();	h1_QDC[1]->Draw();
	c1->cd(6)->SetGridx();	h1_TDC[1]->Draw();

//canvas2
//	const double ch1Ecut_mean = QDCrange/2;
//	const double ch1Ecut_FWHM = QDCrange;
//	const double ch0Ecut_mean = QDCrange/2;
//	const double ch0Ecut_FWHM = QDCrange;

	const double ch1Ecut_mean = 170000;
	const double ch1Ecut_FWHM = 20000;
	const double ch0Ecut_mean = 245000;
	const double ch0Ecut_FWHM = 15000;


//	const double ch1Ecut_mean = 6.98174e+04;
//	const double ch1Ecut_FWHM = 2.71062e+03*5;
//	const double ch0Ecut_mean = 1.47725e+04;
//	const double ch0Ecut_FWHM = 8.75873e+02*5;

	TCanvas *c2 = new TCanvas("QDCFinal", "QDCFinal", 1800,1000);
//	gStyle->SetOptStat(0);
	c2->Divide(3,2);

	c2->cd(1)->SetLogz();
	c2->cd(1)->SetGridx();
	c2->cd(1)->SetGridy();
	TH2D *h2_QDCch1_QDCch0 = new TH2D("Na22 QDC(ch1) by QDC(ch0)", "Na22 QDC(ch1) by QDC(ch0);QDC,ch0 (qdc);QDC,ch1 (qdc)",2000,0,QDCrange,2000,0,QDCrange);
	h2_QDCch1_QDCch0->Draw("COLZ");
	TLine *lEEvl = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,0,ch0Ecut_mean-ch0Ecut_FWHM/2,QDCrange);
	TLine *lEEvu = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,0,ch0Ecut_mean+ch0Ecut_FWHM/2,QDCrange);
	TLine *lEEhl = new TLine(0,ch1Ecut_mean-ch1Ecut_FWHM/2,QDCrange,ch1Ecut_mean-ch1Ecut_FWHM/2);
	TLine *lEEhu = new TLine(0,ch1Ecut_mean+ch1Ecut_FWHM/2,QDCrange,ch1Ecut_mean+ch1Ecut_FWHM/2);
	lEEvl->Draw("same");
	lEEvu->Draw("same");
	lEEhl->Draw("same");
	lEEhu->Draw("same");

	c2->cd(2)->SetLogy();
	c2->cd(2)->SetGridx();
	TH1D *h1_QDCch0 = new TH1D("QDCch0","QDCch0;QDC (qdc);count/qdc",2000,0,QDCrange);
	h1_QDCch0->Draw();
	TLine *lE0vl = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,0,ch0Ecut_mean-ch0Ecut_FWHM/2,100000);
	TLine *lE0vu = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,0,ch0Ecut_mean+ch0Ecut_FWHM/2,100000);
	lE0vl->Draw("same");
	lE0vu->Draw("same");

	c2->cd(3)->SetLogy();
	c2->cd(3)->SetGridx();
	TH1D *h1_QDCch1	= new TH1D("ch1 QDC range used", "ch1 QDC range used;QDC(qdc);count/qdc", 2000, 0, QDCrange);
	h1_QDCch1->Draw();
	TLine *lE1vl = new TLine(ch1Ecut_mean-ch1Ecut_FWHM/2,0,ch1Ecut_mean-ch1Ecut_FWHM/2,100000);
	TLine *lE1vu = new TLine(ch1Ecut_mean+ch1Ecut_FWHM/2,0,ch1Ecut_mean+ch1Ecut_FWHM/2,100000);
	lE1vl->Draw("same");
	lE1vu->Draw("same");

	c2->cd(4)->SetGridx();
	TH1D *h1_TDCdif	= new TH1D(Form("TDCcf ch0-ch1,ch1 QDCcut %.0fqdc",ch1Ecut_mean), Form("TDCcf dif(ch1 %.0fqdc);TDC ch0-ch1(ns);count/0.05ns",ch1Ecut_mean), 200, -3, 3);	// for Ba 172ps
	//TH1D *h1_TDCdif	= new TH1D(Form("TDCcf ch0-ch1,ch1 QDCcut %.0fqdc",ch1Ecut_mean), Form("TDCcf dif(ch1 %.0fqdc);TDC ch0-ch1(ns);count/0.05ns",ch1Ecut_mean),1000,-15,15);
	h1_TDCdif->Draw();

	c2->cd(5)->SetLogz();
	TH2D *h2_TDCdif_QDCch0 = new TH2D(Form("Na22 TDCcf dif(ch0-ch1) by QDC(ch0), (ch1 %.0fqdc)",ch1Ecut_mean),Form("Na22 TDCcf dif by ch0 QDC (ch1 %.0fqdc);QDC,ch0(qdc);TDC dif, ch0-ch1(ns)",ch1Ecut_mean),2000,0,QDCrange,600,-15,15);
	h2_TDCdif_QDCch0->Draw("COLZ");
	TLine *lE0vl2 = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,-15,ch0Ecut_mean-ch0Ecut_FWHM/2,15);
	TLine *lE0vu2 = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,-15,ch0Ecut_mean+ch0Ecut_FWHM/2,15);
	lE0vl2 -> Draw("same");
	lE0vu2 -> Draw("same");

	c2->cd(6)->SetLogz();
	TH2D *h2_TDCdif_QDCch0_nocut = new TH2D("Na22 TDCcf dif(ch0-ch1) by QDC(ch0),nocut ","Na22 TDCcf dif by ch0 QDC, nocut;QDC,ch0(qdc);TDC dif, ch0-ch1(ns)",2000,0,QDCrange, 600, -15, 15);
	h2_TDCdif_QDCch0_nocut->Draw("COLZ");

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
	UInt_t Curr_CoinN, Prev_CoinN, Elap_CoinN; 
	long CurrByte, EndByte;
	int pqrj_flag[4]={0,0,1,0};
	int event_flag=0;
	int fread_flag=1;
	NHitCh=0;
	EventNumber=100;

	int ient=0;
//	for (pqrj_flag[0]=0; ient<numevent; ){
	for (pqrj_flag[0]=0; CoinEventNumber<9000000; ){
		Curr_ient = ient;	Elap_ient = Curr_ient - Prev_ient;
		Curr_CoinN = CoinEventNumber;	Elap_CoinN = Curr_CoinN - Prev_CoinN;
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

			for (int ich=0; ich<NHitCh; ich++){
				rPH[ich] = pulseheight(waveformch[ich]);
				if(1){
				//if(rPH[ich]>= 80){
				//if(rPHch[ich]>=250){
					QDC[ich] = QDCf(waveformch[ich]);
					TDC[ich] = TDCcf(waveformch[ich], waveformtime, 10, 0.15);
					h1_rPH[Channel[ich]]->Fill(rPH[ich]);
					h1_QDC[Channel[ich]]->Fill(QDC[ich]);
					h1_TDC[Channel[ich]]->Fill(TDC[ich]);
				}
			}	
			tree->Fill();

//			if(TDC[0]-TDC[1]>50){ient++; continue;}
//			if(TDC[1]-TDC[0]>50){ient++; continue;}

			if(NHitCh==2){	// holds only for zero-suppressed data
			//if(rPH[0]>= 80 && rPH[1]>= 80){}
				if((ch1Ecut_mean-ch1Ecut_FWHM/2)<QDC[1] && QDC[1]<(ch1Ecut_mean+ch1Ecut_FWHM/2)){
					h1_QDCch0->Fill(QDC[0]);
					h2_TDCdif_QDCch0->Fill(QDC[0], TDC[0] - TDC[1]);
					if((ch0Ecut_mean-ch0Ecut_FWHM/2)<QDC[0] && QDC[0]<(ch0Ecut_mean+ch0Ecut_FWHM/2)){
						h1_TDCdif->Fill(TDC[0]-TDC[1]);
					}
				}
				h2_TDCdif_QDCch0_nocut->Fill(QDC[0],TDC[0]-TDC[1]);
				h2_QDCch1_QDCch0 -> Fill(QDC[0],QDC[1]);
				h1_QDCch1->Fill(QDC[1]);

				CoinEventNumber++;
			}
			NHitCh=0;
			event_flag=0;
		}
			EventNumber=header[4];

		if(fread_flag && header[3]<NOpCh){
		//if(fread_flag && header[4]-EventNumber<=0){
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
			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s, CoinEventNum:%u (%.1f%%) - %.2fHz ", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time, Curr_CoinN,(float)Curr_CoinN*100.0f/(float)ient, (float)Elap_CoinN*1000.0f/(float)Elap_Time);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			Prev_CoinN = Curr_CoinN;			
			if(pqrj_flag[2] || pqrj_flag[3]){
				c1->cd(1);	gPad->Modified();	gPad->Update();
				c1->cd(2);	gPad->Modified();	gPad->Update();
				c1->cd(3);	gPad->Modified();	gPad->Update();
				c1->cd(4);	gPad->Modified();	gPad->Update();
				c1->cd(5);	gPad->Modified();	gPad->Update();
				c1->cd(6);	gPad->Modified();	gPad->Update();
				c2->cd(1);	gPad->Modified();	gPad->Update();
				c2->cd(2);	gPad->Modified();	gPad->Update();
				c2->cd(3);	gPad->Modified();	gPad->Update();
				c2->cd(4);	gPad->Modified();	gPad->Update();
				c2->cd(5);	gPad->Modified();	gPad->Update();
				c2->cd(6);	gPad->Modified();	gPad->Update();
			}
			if(pqrj_flag[3]){
				printf("\n");
				c1->Print(Form("jpg/%s_Decoding.jpg",filename));
				c2->Print(Form("jpg/%s_QDCFinal.jpg",filename));
				pqrj_flag[3]=0;
			}
		}
		if(!fread_flag){
			sleep(1);
			continue;
		}
		ient++;
	}
	c1->cd(1);	gPad->Modified();	gPad->Update();
	c1->cd(2);	gPad->Modified();	gPad->Update();
	c1->cd(3);	gPad->Modified();	gPad->Update();
	c1->cd(4);	gPad->Modified();	gPad->Update();
	c1->cd(5);	gPad->Modified();	gPad->Update();
	c1->cd(6);	gPad->Modified();	gPad->Update();
	c2->cd(1);	gPad->Modified();	gPad->Update();
	c2->cd(2);	gPad->Modified();	gPad->Update();
	c2->cd(3);	gPad->Modified();	gPad->Update();
	c2->cd(4);	gPad->Modified();	gPad->Update();
	c2->cd(5);	gPad->Modified();	gPad->Update();
	c2->cd(6);	gPad->Modified();	gPad->Update();
//	gPad->Modified();
//	gPad->Update();

	printf("\nReading completed.Entries: %d, Event numbers: %u, CoinEventNumber %u\n", ient, EventNumber, CoinEventNumber);

	c2->cd(4);
	TF1 *fftn = new TF1("fftn", "[3]*exp([0]*([0]*[1]*[1]-2*(x-[2]))/2)*TMath::Erfc(-([0]*[1]*[1]+(x-[2]))/(TMath::Sqrt(2)*[1]))+[4]",-1,1);
	fftn->SetParName(0,"lambda");
	fftn->SetParName(1,"sigma");
	fftn->SetParName(2,"mu");
	fftn->SetParName(3,"amp");
	fftn->SetParName(4,"constant");
	h1_TDCdif->Fit(fftn,"M", "", -1, 1);

	printf("T=log(2)/lambda=%.4fns, sigma/1.414= %.4fns\n", log(2)/(fftn->GetParameter(0)), fftn->GetParameter(1)/1.414);

	tree->Write();
	c1->Write();
	c2->Write();
	treefile->Close();

	fclose(datafile);

	printf("\n");
}
