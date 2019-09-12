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
	for (int iwft=0; iwft< 950;iwft++){
		pulsesum = pulsesum + waveform[iwft];
	}
	return pulsesum;
}

void cellindexcorr(Float_t waveform[1024], Int_t StartIndexCell, Float_t waveform2[1024]){
	for (int iwft=0; iwft<1024; iwft++){
		waveform2[(StartIndexCell+iwft)%1024] = waveform[iwft];
	}
//	for (int iwft=900; iwft<1024; iwft++){
//		waveform2[(StartIndexCell+iwft)%1024] = 0;
//	}
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

void cellindexcalib(){
	char filename[90]="testforindex_newcalib_25";
	//char filename[90]="testforindex_nocalib_0";
	//char filename[90]="testforindex_calibon_ch1trig";
	//char filename[90]="testforindex_calibon";
	//char filename[90]="testforindex_caliboff";
	//char filename[90]="testforindex";
	//char filename[90]="terminator_run07_calib1off";
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
	UInt_t StartIndexCellch[NOpCh];
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
	UInt_t StartIndexCell;
	UInt_t NHitCh = NOpCh;	//number of hit channels, maximum NOpCh
	UInt_t Channel[NHitCh];
	UInt_t EventNumber;
	UInt_t TrigTimeTag;
//	UInt_t waveform_length = 1024;	// moved to upper part with each channel
	Float_t waveform[NOpCh][waveform_length];
//	Float_t waveform2[NOpCh][waveform_length];
	Float_t waveformtime[1024];
	Float_t waveformidx[1024];
	Float_t rPH[NHitCh];
	Float_t QDC[NHitCh];
	Float_t TDC[NHitCh];
	Float_t cellcorravg[1024];
	int ahffk=0;

	Float_t samplingrate = 2.5;  // GS/s
	Float_t sampletime = 1.0/samplingrate;  // ns
	for (int iwft=0; iwft<waveform_length; iwft++){
		waveformtime[iwft] = sampletime * iwft;
		waveformidx[iwft]=iwft;
		cellcorravg[iwft]=0;
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
	TCanvas *c2 = new TCanvas("pulseshape_overlap","pulseshape_overlap(onboard calib on)", 1200,800);
	c2->cd()->SetGridy();
	//c2->cd()->SetLogz();
	TH2D *h2_wft_ADC = new TH2D("cellindex_ADC","cellindex_ADCdeviation;cell index;ADC-pulseavg (1V/4096)",1024,0,1024,200,-100, 100);
	h2_wft_ADC->Draw("COLZ");

	TCanvas *c3 = new TCanvas("QDC","QDC", 1200,800);
	c3->SetGridy();
	TH1D *h1_QDC = new TH1D ("QDC","QDC;QDC(ADC*sa);count/100qdc",600,-30000,30000);
	//TH1D *h1_QDC = new TH1D ("QDC","QDC;QDC(ADC*sa);count/100qdc",300,-15000,15000);
	h1_QDC->Draw();

	TCanvas *c4 = new TCanvas("StartIndexCell","StartIndexCell", 1200, 800);
	TH1D *h1_SIC = new TH1D("StartIndexCell","StartIndexCell;StartIndexCell(ch);count", 1024, 0, 1024);
	h1_SIC->Draw();

	TGraph *gr1;
	TGraph *gr2;

	float pulseavg;
	float pedestal;
	TLine *pedline;
	TLine *indexline;
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
		//StartIndexCell=header[2];
		OpCh=header[3];
		StartIndexCellch[OpCh]=header[2];
		EventNumberch[OpCh]=header[4];
		TrigTimeTagch[OpCh]=header[5];

		int Nwaveform=(int) fread(waveformch[OpCh], 4, 1024, datafile);
//		printf("%d %d\n", Nheader, Nwaveform);
/*
		printf("header0,EventSize	 = %d\n",header[0]);
		printf("header1,BoardID		 = %d\n",header[1]);
		printf("header2,StartIndexCell	 = %d\n",header[2]);
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
				//rPHch[ich] = pulseheight(waveformch[ich]);
				QDCch[ich] = QDCf(waveformch[ich]);
				//TDCch[ich] = TDCcf(waveformch[ich], waveformtime, 10, 0.15);
				StartIndexCell=StartIndexCellch[0];
			}
			h1_QDC->Fill(QDCch[0]);
			if(1){
				//printf("ient %d, rPHch[0]=%.1f, QDCch[0]=%.1f\n",ient,rPHch[0],QDCch[0]);
				//pqrc_flag[0]=1;

				pedestal = pedsum(waveformch[0])/200;
				pulseavg = pulsesum(waveformch[0])/ 950;
				cellindexcorr(waveformch[0],StartIndexCell,waveform2ch[0]);

				h1_SIC->Fill(StartIndexCell);

			//	for(int iwft=0; iwft<1000; iwft++){	
			//		//h2_wft_ADC->Fill(iwft,waveformch[0][iwft]-pulseavg);	
			//		h2_wft_ADC->Fill(iwft,waveformch[0][iwft]-pedestal);	
			//	}
				for(int iwft=0; iwft< 950; iwft++){	
					h2_wft_ADC->Fill((StartIndexCell+iwft)%1024,waveformch[0][iwft]-pulseavg);
					cellcorravg[(StartIndexCell+iwft)%1024]=((float)ahffk*cellcorravg[(StartIndexCell+iwft)%1024]+waveformch[0][iwft]-pulseavg)/((float)ahffk+1);
					//cellcorravg[iwft]=((float)EventNumber*cellcorravg[iwft]+waveform2ch[0][iwft]-pulseavg)/((float)EventNumber+1);
					//h2_wft_ADC->Fill((StartIndexCell+iwft)%1024,waveformch[0][iwft]);	
					//h2_wft_ADC->Fill(4*iwft,waveformch[0][iwft]+waveformch[0][iwft+1]+waveformch[0][iwft+2]+waveformch[0][iwft+3]-4*pedestal);	
				}
				ahffk=ahffk+1;

			//	c1->cd();	gPad->Modified(); gPad->Update();
			//	c2->cd();	gPad->Modified(); gPad->Update();
			//	cout<<"0 for next pulse"<<endl;
			//	cin>>pqrc_flag[0];
			}
		}

		if(Elap_Time > 1000){
			c1->cd();
			c1->Clear();
			gr2 = new TGraph(1024,waveformidx,waveform2ch[0]);
			//gr2 = new TGraph( 950,waveformtime,waveform2ch[0]);
			gr2->GetYaxis()->SetRangeUser(pulseavg-100,pulseavg+100);
			gr2->SetLineColor(2);
			gr2->Draw();
			pedline = new TLine(0,pulseavg,1024,pulseavg);
			pedline->SetLineColor(3);
			pedline->Draw("SAME");
			indexline = new TLine(StartIndexCell,pulseavg-100,StartIndexCell,pulseavg+100);
			//indexline = new TLine(StartIndexCell*0.4,pulseavg-100,StartIndexCell*0.4,pulseavg+100);
			indexline->SetLineColor(4);
			indexline->Draw("SAME");


			//printf("\rient %d, rPHch[0]=%.1f, QDCch[0]=%.1f",ient,rPHch[0],QDCch[0]);
			if(pqrc_flag[2]){
				c1->cd();	gPad->Modified(); gPad->Update();
				c2->cd();	gPad->Modified(); gPad->Update();
				c3->cd();	gPad->Modified(); gPad->Update();
				c4->cd();	gPad->Modified(); gPad->Update();
			}
			printf("\rdecoding... at %d of %u (%.3f) %.0fient/s, CoinEventNum:%u (%.1f%%) - %.2fHz ", ient, numevent, (float) ient/numevent, (float)Elap_ient*1000.0f/(float)Elap_Time, Curr_CoinN,(float)Curr_CoinN*100.0f/(float)ient, (float)Elap_CoinN*1000.0f/(float)Elap_Time);
			fflush(stdout);
			Prev_Time = Curr_Time;
			Prev_ient = Curr_ient;
			Prev_CoinN = Curr_CoinN;			
		}
	ient++;
	}
	c2->cd();
	TGraph *gr3 = new TGraph(1024, waveformidx, cellcorravg);
	//gr3->Draw("SAME");

	c1->cd();
	gr2->SetTitle("waveform;cell index;ADC");
	TLatex *tltx1 = new TLatex( 700, pulseavg+50,Form("#splitline{#splitline{%s}{terminator 50#Omega}}{#splitline{DC offset -16%%}{new calib}}",filename));
	tltx1->SetTextSize(0.03);
	tltx1->Draw("SAME");

	c2->cd();
	TLatex *tltx2 = new TLatex( 700, 50, Form("#splitline{#splitline{%s}{terminator 50#Omega}}{#splitline{DC offset -16%%}{new calib}}",filename));
	tltx2->SetTextSize(0.03);
	tltx2->Draw("SAME");

	c3->cd();
	TLatex *tltx3 = new TLatex( 7000,5000, Form("#splitline{#splitline{%s}{terminator 50#Omega}}{#splitline{DC offset -16%%}{new calib}}",filename));
	tltx3->SetTextSize(0.03);
	tltx3->Draw("SAME");

	c4->cd();
	TLatex *tltx4 = new TLatex( 7000,5000, Form("#splitline{#splitline{%s}{terminator 50#Omega}}{#splitline{DC offset -16%%}{new calib}}",filename));
	tltx4->SetTextSize(0.03);
	tltx4->Draw("SAME");




	printf("\nReading completed.Entries: %d, Event numbers: %u, CoinEventNumber %u\n", ient, EventNumber, CoinEventNumber);

	printf("%d\n",ahffk);

	for (int i=0; i<128; i++){
		for (int j=0; j<8; j++){
			//printf("%.0f	%.0f    %.0f    %.0f    %.0f    %.0f    %.0f    %.0f    cell = %d to %d", cellcorravg[i
			printf("%.0f\t", cellcorravg[i*8+j]);
		}
		printf("cell = %d to %d\n", i*8,i*8+7);
	}

}
