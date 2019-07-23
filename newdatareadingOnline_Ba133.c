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

void KeyboardCommand(int pqrc_flag[3]){		
	if(!kbhit())	return;
	int key = getch();
	switch(key){
		case 'p' :						// 'p'ause
			pqrc_flag[0]^=1; 
			if( pqrc_flag[0]) printf("\nPaused. Press p again to resume\n");
			if(!pqrc_flag[0]) printf("resumed\n");
			break;
		case 'q' : 
			pqrc_flag[1]^=1; break;			// 'q'uit
		case 'r' : 						// 'r'enew histograms
			pqrc_flag[2]^=1; 
			if( pqrc_flag[2]) printf("\nRenewing histograms is enabled.\n");
			if(!pqrc_flag[2]) printf("\nRenewing histograms is disabled.\n");
			break;
//		case 'c' : 
//			pqrc_flag[3]=1; break;			// new 'c'anvas
		default : break;
	}
}
/*
void NewCanvas(TCanvas *c1, TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4);
	c1 = new TCanvas("Decoding","Decoding",1200,800);
	c1->Divide(2,2);
	c1->cd(1);
//	h1_QDC[0]->GetXaxis()->SetRangeUser(5000,300000);
	h1->Draw();
	c1->cd(2);
//	h1_QDC[1]->GetXaxis()->SetRangeUser(5000,300000);
	h2->Draw();
	c1->cd(3);
	h3->Draw();
	c1->cd(4);
	h4->Draw();
}
*/

//////////////////////main function////////////////////////////////

void newdatareadingOnline_Ba133(){
	char filename[90]="Ba133_att_run01";

//Co60_att_run01
//Na22_att_run04
//Cs137_att_run01
//Ba133_att_run01
//MixSource_att_run01

//tree

	TFile *treefile = new TFile(Form("../../DecodedDataFiles/%s.root",filename),"recreate");
	TString str0 = TString("Set_info\n");
	TString str1 = TString("Ba133, attenuator \n");
	TString str2 = TString(" ch0: Scint-38*38LaBr3(20170809) PMT-R13408(AA0037) HV-1500V,att9dB from HVch1\n");
	TString str3 = TString(" ch1: Scint-38*38LaBr3(UNKNOWN) PMT-R13408(AA0038) HV-1500V,att4dB from HVch0\n");
	TObjString *TObjS = new TObjString(str0+str1+str2+str3);
	TObjS->Write();


	UInt_t NOpCh = 2;	// number of open channels
	UInt_t OpCh;
	UInt_t waveform_length = 1024;
	Float_t waveformch[NOpCh][waveform_length];
	UInt_t EventNumberch[NOpCh];
	UInt_t TrigTimeTagch[NOpCh];
	Float_t rPHch[NOpCh];
	Float_t QDCch[NOpCh];
	Float_t TDCch[NOpCh];

	TTree *tree = new TTree("wavedata","wavedata");
	UInt_t CoinEventNumber = 0;
	UInt_t EventSize;
	UInt_t BoardID;
	UInt_t Group;
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
	tree->Branch("Group", &Group, "Group/i");
	tree->Branch("CoinEventNumber", &CoinEventNumber, "CoinEventNumber/i");
	tree->Branch("NHitCh", &NHitCh, "NHitCh/i");
	tree->Branch("Channel", Channel, "Channel[NHitCh]/i");
	tree->Branch("waveform_length", &waveform_length, "waveform_length/i");
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

	h1_rPH[0] = new TH1D("h1_rPH[0]","h1_rPH[0];rPH;count/rPH", 4096, 0, 4096);
	h1_QDC[0] = new TH1D("h1_QDC[0]","h1_QDC[0];QDC;count/150QDC", 2000, 0, 300000);
	h1_TDC[0] = new TH1D("h1_TDC[0]","h1_TDC[0];TDC (ns);count/0.4ns",1024,0,410);
	h1_rPH[1] = new TH1D("h1_rPH[1]","h1_rPH[1];rPH;count/rPH", 4096, 0, 4096);
	h1_QDC[1] = new TH1D("h1_QDC[1]","h1_QDC[1];QDC;count/150QDC", 2000, 0, 300000);
	h1_TDC[1] = new TH1D("h1_TDC[1]","h1_TDC[1];TDC (ns);count/0.4ns",1024,0,410);
//	NewCanvas(c1, h1_QDC[0], h1_QDC[1], h1_TDC[0], h1_TDC[1]);
	c1->Divide(3,2);
	c1->cd(1);	h1_rPH[0]->Draw();
	c1->cd(2);	h1_QDC[0]->Draw();
	c1->cd(3);	h1_TDC[0]->Draw();
	c1->cd(4);	h1_rPH[1]->Draw();
	c1->cd(5);	h1_QDC[1]->Draw();
	c1->cd(6);	h1_TDC[1]->Draw();

//canvas2
	const double ch1Ecut_mean = 150000;
	const double ch1Ecut_FWHM = 300000;
	const double ch0Ecut_mean = 150000;
	const double ch0Ecut_FWHM = 300000;

//	const double ch1Ecut_mean = 7.41401e+04;
//	const double ch1Ecut_FWHM = 2.23946e+03*5;
//	const double ch0Ecut_mean = 8.09619e+04;
//	const double ch0Ecut_FWHM = 2.21135e+03*5;

	TCanvas *c2 = new TCanvas("QDCFinal", "QDCFinal", 1800,1000);
//	gStyle->SetOptStat(0);
	c2->Divide(3,2);

	c2->cd(1)->SetLogz();
	TH2D *h2_QDCch1_QDCch0 = new TH2D("Ba133 QDC(ch1) by QDC(ch0)", "Ba133 QDC(ch1) by QDC(ch0);QDC,ch0 (qdc);QDC,ch1 (qdc)",2000,0,300000,2000,0,300000);
	h2_QDCch1_QDCch0->Draw("COLZ");
	TLine *lEEvl = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,0,ch0Ecut_mean-ch0Ecut_FWHM/2,300000);
	TLine *lEEvu = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,0,ch0Ecut_mean+ch0Ecut_FWHM/2,300000);
	TLine *lEEhl = new TLine(0,ch1Ecut_mean-ch1Ecut_FWHM/2,300000,ch1Ecut_mean-ch1Ecut_FWHM/2);
	TLine *lEEhu = new TLine(0,ch1Ecut_mean+ch1Ecut_FWHM/2,300000,ch1Ecut_mean+ch1Ecut_FWHM/2);
	lEEvl->Draw("same");
	lEEvu->Draw("same");
	lEEhl->Draw("same");
	lEEhu->Draw("same");

	c2->cd(2)->SetLogy();
	TH1D *h1_QDCch0 = new TH1D("QDCch0","QDCch0;QDC (qdc);count/qdc",2000,0,300000);
	h1_QDCch0->Draw();
	TLine *lE0vl = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,0,ch0Ecut_mean-ch0Ecut_FWHM/2,100000);
	TLine *lE0vu = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,0,ch0Ecut_mean+ch0Ecut_FWHM/2,100000);
	lE0vl->Draw("same");
	lE0vu->Draw("same");

	c2->cd(3)->SetLogy();
	TH1D *h1_QDCch1	= new TH1D("ch1 QDC range used", "ch1 QDC range used;QDC(qdc);count/qdc", 2000, 0, 300000);
	h1_QDCch1->Draw();
	TLine *lE1vl = new TLine(ch1Ecut_mean-ch1Ecut_FWHM/2,0,ch1Ecut_mean-ch1Ecut_FWHM/2,100000);
	TLine *lE1vu = new TLine(ch1Ecut_mean+ch1Ecut_FWHM/2,0,ch1Ecut_mean+ch1Ecut_FWHM/2,100000);
	lE1vl->Draw("same");
	lE1vu->Draw("same");

	c2->cd(4);
	TH1D *h1_TDCdif	= new TH1D(Form("TDCcf ch0-ch1,ch1 QDCcut %.0fqdc",ch1Ecut_mean), Form("TDCcf dif(ch1 %.0fqdc);TDC ch0-ch1(ns);count/0.05ns",ch1Ecut_mean),1000,-5,5);
	h1_TDCdif->Draw();

	c2->cd(5)->SetLogz();
	TH2D *h2_TDCdif_QDCch0 = new TH2D(Form("Ba133 TDCcf dif(ch0-ch1) by QDC(ch0), (ch1 %.0fqdc)",ch1Ecut_mean),Form("Ba133 TDCcf dif by ch0 QDC (ch1 %.0fqdc);QDC,ch0(qdc);TDC dif, ch0-ch1(ns)",ch1Ecut_mean),2000,0,300000,600,-15,15);
	h2_TDCdif_QDCch0->Draw("COLZ");
	TLine *lE0vl2 = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,-15,ch0Ecut_mean-ch0Ecut_FWHM/2,15);
	TLine *lE0vu2 = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,-15,ch0Ecut_mean+ch0Ecut_FWHM/2,15);
	lE0vl2 -> Draw("same");
	lE0vu2 -> Draw("same");

	c2->cd(6)->SetLogz();
	TH2D *h2_TDCdif_QDCch0_nocut = new TH2D("Ba133 TDCcf dif(ch0-ch1) by QDC(ch0),nocut ","Ba133 TDCcf dif by ch0 QDC, nocut;QDC,ch0(qdc);TDC dif, ch0-ch1(ns)",2000,0,300000, 600, -15, 15);
	h2_TDCdif_QDCch0_nocut->Draw("COLZ");

//file

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
	printf("%s.dat Filesize: %lld(Bytes), numevent: %u, Error:%u\n", filename, Filesize, numevent, remainder);

	uint32_t header[6];
	int Nheader, Nwaveform;
	uint64_t Curr_Time, Prev_Time, Elap_Time;
	int Curr_ient, Prev_ient, Elap_ient;
	UInt_t Curr_CoinN, Prev_CoinN, Elap_CoinN; 
	long CurrByte, EndByte;
	int pqrc_flag[3]={0,0,1};
	int event_flag=1;

	int ient=0;
	for (pqrc_flag[0]=0; ient<numevent; ){
//	for (pqrc_flag[0]=0; CoinEventNumber<9000000; ){
		Curr_ient = ient;	Elap_ient = Curr_ient - Prev_ient;
		Curr_CoinN = CoinEventNumber;	Elap_CoinN = Curr_CoinN - Prev_CoinN;
		Curr_Time = get_time();	Elap_Time = Curr_Time - Prev_Time; 

		KeyboardCommand(pqrc_flag);
		if(pqrc_flag[1]) break;
		if(pqrc_flag[0]) continue;

		CurrByte = ftell(datafile);
		fseek(datafile,0,SEEK_END);
		EndByte = ftell(datafile);
		numevent = EndByte/4120;
		fseek(datafile,CurrByte,SEEK_SET);
		if(EndByte-CurrByte<4120){
			printf("(ONLINE)");
			fflush(stdout);
			sleep(1);
			continue;
		}

		Nheader=(int)fread(header, 4, 6, datafile);

		EventSize=header[0];
		BoardID=header[1];
		Group=header[2];
		OpCh=header[3];
		EventNumberch[OpCh]=header[4];
		TrigTimeTagch[OpCh]=header[5];

		int Nwaveform=(int) fread(waveformch[OpCh], 4, 1024, datafile);
//		printf("%d %d\n", Nheader, Nwaveform);
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
			printf("\nData Size Error!!  ient: %d, header[0]: %d, Nheader+Nwaveform: %d\n", ient, header[0], Nheader+Nwaveform);
			break;
		}

		event_flag=1;
		for (int ich=0; ich<NOpCh-1; ich++){
			if(EventNumberch[ich]!=EventNumberch[(ich+1)]){event_flag=0; break;}
			if(TrigTimeTagch[ich]!=TrigTimeTagch[(ich+1)]){event_flag=0; break;}
		}
		if(event_flag){
//		if(EventNumberch[0]==EventNumberch[1] && TrigTimeTagch[0]==TrigTimeTagch[1]){
			EventNumber = EventNumberch[0]; TrigTimeTag = TrigTimeTagch[0];

			NHitCh=0;
			for (int ich=0; ich<NOpCh; ich++){
				rPHch[ich] = pulseheight(waveformch[ich]);
				if(rPHch[ich]>=80){
					QDCch[ich] = QDCf(waveformch[ich]);
					TDCch[ich] = TDCcf(waveformch[ich], waveformtime, 10, 0.15);
					h1_rPH[ich]->Fill(rPHch[ich]);
					h1_QDC[ich]->Fill(QDCch[ich]);
					h1_TDC[ich]->Fill(TDCch[ich]);

					Channel[NHitCh]=ich;
//					waveformcopy(waveform[NHitCh],waveformch[ich]);
					rPH[NHitCh] = rPHch[ich];
					QDC[NHitCh] = QDCch[ich];
					TDC[NHitCh] = TDCch[ich];
					NHitCh=NHitCh+1;
				}
			}	
			if(NHitCh==0){ient++; continue;}
			tree->Fill();

//			if(TDC[0]-TDC[1]>50){ient++; continue;}
//			if(TDC[1]-TDC[0]>50){ient++; continue;}
			if(NHitCh==2){
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
		}

		if(Elap_Time > 1000){
			if(pqrc_flag[2]){
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
			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s, CoinEventNum:%u (%.1f%%) - %.2fHz ", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time, Curr_CoinN,(float)Curr_CoinN*100.0f/(float)ient, (float)Elap_CoinN*1000.0f/(float)Elap_Time);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			Prev_CoinN = Curr_CoinN;			
		}
	ient++;
	}
	gPad->Modified();
	gPad->Update();

	printf("\nReading completed.Entries: %d, Event numbers: %u, CoinEventNumber %u\n", ient, EventNumber, CoinEventNumber);

	c2->cd(4);
	TF1 *fftn = new TF1("fftn", "[3]*exp([0]*([0]*[1]*[1]-2*(x-[2]))/2)*TMath::Erfc(-([0]*[1]*[1]+(x-[2]))/(TMath::Sqrt(2)*[1]))",-1,1);
	fftn->SetParName(0,"lambda");
	fftn->SetParName(1,"sigma");
	fftn->SetParName(2,"mu");
	fftn->SetParName(3,"amp");
	h1_TDCdif->Fit(fftn,"M", "", -1, 1);

	tree->Write();
	c1->Write();
	c2->Write();
	treefile->Close();

	fclose(datafile);

	printf("\n");
}
