float getminimum(Float_t waveform[520]){
	float min=4100.0;
	for (int iwft=0; iwft<520; iwft++){
		if(waveform[iwft]<min) min=waveform[iwft];
	}
	return min;
}

float pedsum(Float_t waveform[520]){
//double pedsum(Float_t waveform[520]){
	Float_t pedsum=0;
	for (int iwft=0; iwft<100;iwft++){
		pedsum = pedsum + waveform[iwft];
	}
	return pedsum;
}

double TDCle(Float_t waveform[520],Float_t waveformtime[520], double threshold){
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	Double_t pedestal = (double)pedsum(waveform)/100;
	for (int iwft=100; iwft<520; iwft++){
		if(waveform[iwft]<pedestal-threshold){
			iwfttorecord = (double)iwft+(pedestal-threshold-waveform[iwft])/(waveform[iwft]-waveform[iwft-1]);
			timetorecord = (double)iwfttorecord*0.4;
			break;
		}
	}
	return timetorecord;
}

double pulseheight(Float_t waveform[520]){
	double ADCmin=4098;
	for (int i=0; i<520;i++){
		if (ADCmin>waveform[i]){
			ADCmin = waveform[i];
		}
	}
	double pulseheight = (double) pedsum(waveform)/200-ADCmin;
	return pulseheight;
}

double TDChf(Float_t waveform[520],Float_t waveformtime[520], double fraction){
	double threshold = fraction * pulseheight(waveform); 
	return TDCle(waveform,waveformtime,threshold);
}	

double TDCcf3(Float_t waveform[520], Float_t waveformtime[520],Int_t delay, double fraction){
	double pedestal=(double)pedsum(waveform)/100;
	Double_t waveform2[520-delay];
	for (int i=0; i<520-delay; i++){
		waveform2[i]=(double) -fraction*(pedestal-waveform[i+delay])+(pedestal-waveform[i]);
	}
	Int_t minwft=0;
	Int_t maxwft=0;
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	for (int iwft=100-delay; iwft<300-delay; iwft++){
//	for (int iwft=200; iwft<520-delay; iwft++){
		if(waveform2[iwft]>waveform2[maxwft]){maxwft=iwft;}
		if(waveform2[iwft]<waveform2[minwft]){minwft=iwft;}
	}
	int iwft=0;
	for (iwft=minwft; iwft<maxwft; iwft++){
		if(waveform2[iwft]<waveform2[iwft+1] && waveform2[iwft+1]<=0 && 0<waveform2[iwft+2] && waveform2[iwft+2]<waveform2[iwft+3]){
			break;
		}
	}
	Double_t y_wf[4] = {waveform2[iwft],waveform2[iwft+1],waveform2[iwft+2],waveform2[iwft+3]};
	Double_t x_wft[4] = {(double)iwft,(double)iwft+1,(double)iwft+2,(double)iwft+3};
	TGraph *gr = new TGraph(4, x_wft,y_wf);
	TF1 *f1 = new TF1("f1", "pol3", iwft, iwft+3);
	gr->Fit(f1, "Q","",iwft,iwft+3);	// for the best? 4 by 4 linear equation can be solved... make some practice!
	iwfttorecord = f1->GetX(0, iwft, iwft+3);

	//Double_t x_wf[4] = {waveform2[iwft],waveform2[iwft+1],waveform2[iwft+2],waveform2[iwft+3]};
	//Double_t y_wft[4] = {(double)iwft,(double)iwft+1,(double)iwft+2,(double)iwft+3};
	//TGraph *gr = new TGraph(4, x_wf,y_wft);
	//TF1 *f1 = new TF1("f1", "pol3", waveform2[iwft], waveform2[iwft+3]+0.0001);
	//gr->Fit(f1, "Q","",waveform2[iwft],waveform2[iwft+3]);	// for the best? 4 by 4 linear equation can be solved... make some practice!
	//iwfttorecord = f1->Eval(0,0,0,0);

	timetorecord = (double)iwfttorecord*0.4;	//0.4ns per sample
	return timetorecord;
}

double TDCcf(Float_t waveform[520], Float_t waveformtime[520],Int_t delay, double fraction){
	double pedestal=(double)pedsum(waveform)/200;
	Double_t waveform2[520-delay];
	for (int i=0; i<520-delay; i++){
		waveform2[i]=(double) -fraction*(pedestal-waveform[i+delay])+(pedestal-waveform[i]);
	}
	Int_t minwft=0;
	Int_t maxwft=0;
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	for (int iwft=100; iwft<250; iwft++){
		if(waveform2[iwft]>waveform2[maxwft]){maxwft=iwft;}
		if(waveform2[iwft]<waveform2[minwft]){minwft=iwft;}
	}
	for (int iwft=minwft; iwft<maxwft; iwft++){
		if(waveform2[iwft-1]<0 && 0<waveform2[iwft] && waveform2[iwft]<waveform2[iwft+1] && waveform2[iwft+1]<waveform2[iwft+2]){
			iwfttorecord = (double)iwft-waveform2[iwft]/(waveform2[iwft]-waveform2[iwft-1]);
			timetorecord = (double)iwfttorecord*0.4;	//0.4ns per sample
			break;
		}
	}
	return timetorecord;
}
/*
float QDCf(Float_t waveform[520]){
	Float_t pedestal=pedsum(waveform)/200.0;
	Float_t QDCsum=0;
	for (int iwft=200; iwft<950; iwft++){
		QDCsum=QDCsum+pedestal-waveform[iwft];
	}
	Float_t QDC = QDCsum/750.0;
	return QDC;
}
*/

float QDCf(Float_t waveform[520]){
	Float_t QDCsum=0;
	for (int iwft=100; iwft<480; iwft++){
		QDCsum=QDCsum+waveform[iwft];
	}
	Float_t QDC = (480.0-100.0)/100.0*pedsum(waveform)-QDCsum;
	return QDC;
}

/*
double energy(float QDC, int Channel){
	double energy;
	if (Channel==0){
		energy = (double) 0.0415113*(double)QDC -1.54326;
	}
	if (Channel==1){
		energy = (double) 0.0377387*(double)QDC -3.99942;
	}
	return energy;
}
*/

