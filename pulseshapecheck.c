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

void fastfilter(Float_t waveform[1024], Int_t uftl, Int_t dftl, Float_t waveform2[1024]){
	double pedestal=(double)pedsum(waveform)/200;
	float fraction[uftl+dftl];
	for (int i=0; i<uftl; i++){	fraction[i]=-1.0/(float)uftl;	}
	for (int i=uftl; i<uftl+dftl; i++){	fraction[i]=1.0/(float)dftl;	}

	//printf("fastfilter fraction = %f",fraction[0]);
	for (int i=0; i<1024-uftl-dftl; i++){
		waveform2[i]=0;
		for (int j=0; j<uftl+dftl; j++){
			waveform2[i]=waveform2[i]+fraction[j]*waveform[i+j];
		}
	waveform2[i]=waveform2[i]+pedestal;
	}
	for (int i=1024-uftl-dftl; i<1024;i++){
		waveform2[i]=pedestal;
	}
}
double pulsesum(Float_t waveform[1024]){
	Float_t pulsesum=0;
	for (int iwft=0; iwft<1000;iwft++){
		pulsesum = pulsesum + waveform[iwft];
	}
	return pulsesum;
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

void pulseshapecheck(){
	char filename[90]="terminator_run07_calib1off";
	//char filename[90]="terminator_run04_DRS4";
	//char filename[90]="Ba133_att00_run04_Xsc";
	FILE *datafile = fopen(Form("%s.dat",filename),"rb");	// moved from latter to here
	if(datafile==NULL){
		fputs("File error\n",stderr);
		exit(1);
	}

//tree


	UInt_t NOpCh = 2;	// number of open channels
	UInt_t OpCh;
	UInt_t waveform_length = 1024;
	Float_t waveformch[NOpCh][1024];
	Float_t waveform2ch[NOpCh][1024];
	UInt_t EventNumberch[NOpCh];
	UInt_t TrigTimeTagch[NOpCh];
	Float_t rPHch[NOpCh];
	Float_t QDCch[NOpCh];
	Float_t TDCch[NOpCh];

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
//	Float_t waveform2[NOpCh][waveform_length];
	Float_t waveformtime[1024];
	Float_t rPH[NHitCh];
	Float_t QDC[NHitCh];
	Float_t TDC[NHitCh];

	Float_t samplingrate = 2.5;  // GS/s
	Float_t sampletime = 1.0/samplingrate;  // ns
	for (int iwft=0; iwft<waveform_length; iwft++){
		waveformtime[iwft] = sampletime * iwft;
	}
	printf("samplingrate %f, sampletime %f, waveform_length %d, fullwft %f\n", samplingrate, sampletime, waveform_length, sampletime*waveform_length);


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
	int pqrc_flag[3]={0,0,1};
	int event_flag=1;

	TCanvas *c1 = new TCanvas("pulseshape","pulseshape", 1200,800);
	c1->SetGridy();
	TCanvas *c2 = new TCanvas("pulseshape_overlap","pulseshape_overlap", 1200,800);
	c2->cd()->SetGridy();
	c2->cd()->SetLogz();
	TH2D *h2_wft_ADC = new TH2D("wft_ADC","wft_adc;waveformtime(4sample);ADC(1V/4096)",128,0,1024,100,-100,+100);
	//TH2D *h2_wft_ADC = new TH2D("wft_ADC","wft_adc;waveformtime(4sample);ADC(1V/4096)",128,0,1024,100,3750,3850);
	h2_wft_ADC->Draw("COLZ");

	float pulseavg;
	float pedestal;
	TLine *pedline;
	int ient=0;
//	for (pqrc_flag[0]=0; ient<numevent; ){
	for (pqrc_flag[0]=0; CoinEventNumber<9000000; ){
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
			for (int ich=0; ich<1; ich++){
			//for (int ich=0; ich<NOpCh; ich++){
				rPHch[ich] = pulseheight(waveformch[ich]);
				QDCch[ich] = QDCf(waveformch[ich]);
				TDCch[ich] = TDCcf(waveformch[ich], waveformtime, 10, 0.15);
			}
			//if(rPHch[0]>1000&&rPHch[0]<1200){
			//if(rPHch[0]>150&&rPHch[0]<180){
			//if(QDCch[0]>2000){
			if(1){
			//if(QDCch[0]<-2000){
			//if(QDCch[0]<-700){
			//if(rPHch[0]<300){
				printf("ient %d, rPHch[0]=%.1f, QDCch[0]=%.1f\n",ient,rPHch[0],QDCch[0]);
				pqrc_flag[0]=1;

				pedestal = pedsum(waveformch[0])/200;
				pulseavg = pulsesum(waveformch[0])/1000;
				pedline = new TLine(0,pedestal,450,pedestal);

				//TCanvas *c1 = new TCanvas("pulseshape","pulseshape", 800,600);
				c1->Clear();
				c1->cd();
				TGraph *gr1 = new TGraph(1024,waveformtime,waveformch[0]);
				gr1->GetYaxis()->SetRangeUser(pedestal-50,pedestal+50);
				gr1->Draw();

				fastfilter(waveformch[0],12,12,waveform2ch[0]);
				TGraph *gr2 = new TGraph(1024,waveformtime,waveform2ch[0]);
				gr2->GetYaxis()->SetRangeUser(3900,4096);
				gr2->SetLineColor(2);
				//gr2->Draw("SAME");
				pedline->SetLineColor(3);
				pedline->Draw("SAME");

			//	for(int iwft=0; iwft<1000; iwft++){	
			//		//h2_wft_ADC->Fill(iwft,waveformch[0][iwft]-pulseavg);	
			//		h2_wft_ADC->Fill(iwft,waveformch[0][iwft]-pedestal);	
			//	}
				for(int iwft=0; iwft<1000/4; iwft++){	
					//h2_wft_ADC->Fill(iwft,waveformch[0][iwft]-pulseavg);	
					h2_wft_ADC->Fill(4*iwft,waveformch[0][iwft]+waveformch[0][iwft+1]+waveformch[0][iwft+2]+waveformch[0][iwft+3]-4*pedestal);	
				}


				c1->cd();	gPad->Modified(); gPad->Update();
				c2->cd();	gPad->Modified(); gPad->Update();
			//	cout<<"0 for next pulse"<<endl;
				cin>>pqrc_flag[0];
			}
		}
/*
		if(Elap_Time > 1000){
			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s, CoinEventNum:%u (%.1f%%) - %.2fHz ", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time, Curr_CoinN,(float)Curr_CoinN*100.0f/(float)ient, (float)Elap_CoinN*1000.0f/(float)Elap_Time);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			Prev_CoinN = Curr_CoinN;			
		}
*/	ient++;
	}

	printf("\nReading completed.Entries: %d, Event numbers: %u, CoinEventNumber %u\n", ient, EventNumber, CoinEventNumber);

	printf("\n");
}
