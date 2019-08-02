double getminimum(Float_t waveform[520]){
	double min=4100.0;
	for (int iwft=0; iwft<520; iwft++){
		if(waveform[iwft]<min) min=(float)waveform[iwft];
	}
	return min;
}

double pedsum(Float_t waveform[520]){
	int pedsum=0;
	for (int iwft=0; iwft<80;iwft++){
		pedsum = pedsum + waveform[iwft];
	}
	return pedsum;
}

double TDCle(Float_t waveform[520],Float_t waveformtime[520], double threshold){
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	Double_t pedestal = (double)pedsum(waveform)/80;
	for (int iwft=80; iwft<520; iwft++){
		if(waveform[iwft]<pedestal-threshold){
			iwfttorecord = (double)iwft+(pedestal-threshold-waveform[iwft])/(waveform[iwft]-waveform[iwft-1]);
			timetorecord = (double)iwfttorecord*1.0;
			break;
		}
	}
	return timetorecord;
}

float pulseheight(Float_t waveform[520]){
	float ADCmin=4098;
	for (int i=0; i<520;i++){
		if (ADCmin>waveform[i]){
			ADCmin = waveform[i];
		}
	}
	float pulseheight = (float) pedsum(waveform)/80-ADCmin;
	return pulseheight;
}

double TDChf(Float_t waveform[520],Float_t waveformtime[520], double fraction){
	double threshold = fraction * pulseheight(waveform); 
	return TDCle(waveform,waveformtime,threshold);
}	

double TDCcf(Float_t waveform[520], Float_t waveformtime[520],Int_t delay, double fraction){
	double pedestal=(double)pedsum(waveform)/80;
	Double_t waveform2[520-delay];
	for (int i=0; i<520-delay; i++){
		waveform2[i]=(double) -fraction*(pedestal-waveform[i+delay])+(pedestal-waveform[i]);
	}
	Int_t minwft=0;
	Int_t maxwft=0;
	Double_t iwfttorecord=-1;
	Double_t timetorecord=-2;
	for (int iwft=80; iwft<300; iwft++){
//	for (int iwft=80; iwft<200-delay; iwft++){
		if(waveform2[iwft]>waveform2[maxwft]){maxwft=iwft;}
		if(waveform2[iwft]<waveform2[minwft]){minwft=iwft;}
	}
	for (int iwft=minwft; iwft<maxwft; iwft++){
		if(0<waveform2[iwft] && waveform2[iwft]<waveform2[iwft+1] && waveform2[iwft+1]<waveform2[iwft+2]){
			iwfttorecord = (double)iwft-waveform2[iwft]/(waveform2[iwft]-waveform2[iwft-1]);
			timetorecord = (double)iwfttorecord*1.0;	//1.0ns per sample
			break;
		}
	}
	return timetorecord;
}

float QDCf(Float_t waveform[520]){
	Float_t QDCsum=0;
	for (int iwft=80; iwft<500; iwft++){
		QDCsum=QDCsum+waveform[iwft];
	}
	Float_t QDC = 5.25 * pedsum(waveform) - QDCsum;
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

