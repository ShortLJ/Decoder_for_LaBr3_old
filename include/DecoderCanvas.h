class DECODER_CANVAS{
	private:
		UInt_t NOpCh;
	public:
		TCanvas* Canvas_Decoding;
		TCanvas** Coincidence;

		void Construct();
		void Destruct();
		void Modified_Update(char flag);
		void Print(char _filename[90]);
		//void CanvasWrite(TFile *_treefile);

		DECODER_CANVAS(UInt_t _NOpCh){
			NOpCh=_NOpCh;
			Canvas_Decoding=NULL;
			Coincidence=NULL;
		}
		~DECODER_CANVAS(){
		}
};

void DECODER_CANVAS::Construct(){
	Canvas_Decoding = new TCanvas("Decoding","Decoding",1800, 1000);
	Canvas_Decoding->Divide(NOpCh,3);

	Coincidence = new TCanvas*[3];
	Coincidence[0] = new TCanvas("QDC_H1_COIN", "QDC_H1_COIN", 900, 954);
	Coincidence[1] = new TCanvas("QDC_H2_COIN", "QDC_H1_COIN", 900, 954);
	Coincidence[2] = new TCanvas("TDCdif_H1_COIN", "TDCdif_H1_COIN", 900, 954);
	for(int i=0; i<3; i++)	Coincidence[i]->Divide(NOpCh,NOpCh);

	return;
}

void DECODER_CANVAS::Destruct(){
	if(!Canvas_Decoding) delete Canvas_Decoding;
	if(!Coincidence[0]) delete Coincidence[0];
	if(!Coincidence[1]) delete Coincidence[1];
	if(!Coincidence[2]) delete Coincidence[2];
	if(!Coincidence) delete[] Coincidence;
	return;
}

void DECODER_CANVAS::Modified_Update(char flag){
	for (int ich=0; ich<NOpCh; ich++){
		if(flag & 0x01){
			Canvas_Decoding->cd(1+ich+0*NOpCh); gPad->Modified(); gPad->Update();
			Canvas_Decoding->cd(1+ich+1*NOpCh); gPad->Modified(); gPad->Update();
			Canvas_Decoding->cd(1+ich+2*NOpCh); gPad->Modified(); gPad->Update();
		}
		for (int jch=0; jch<NOpCh; jch++){
			if(flag & 0x02){ 
				Coincidence[0]->cd(1+ich+jch*NOpCh); gPad->Modified(); gPad->Update();
			}
			if(flag & 0x04){ 
				Coincidence[1]->cd(1+ich+jch*NOpCh); gPad->Modified(); gPad->Update();
			}
			if(flag & 0x08){
				Coincidence[2]->cd(1+ich+jch*NOpCh); gPad->Modified(); gPad->Update();
			}
		}
	}
	return;
}

void DECODER_CANVAS::Print(char _filename[90]){
	Canvas_Decoding->Print(Form("jpg/%s_Decoding_c1.jpg",_filename));
	for (int i=0; i<3; i++){
	Coincidence[i]->Print(Form("jpg/%s_Decoding_c%d.jpg",_filename,i+1));
	}
}
/*
void DECODER_CANVAS::CanvasWrite(TFile *_treefile){
	_treefile->cd();
	Canvas_Decoding->Write();
	for (int i=0; i<3; i++){
	Coincidence[i]->Write();
	}
}
*/

////////////////////////////////////////////

class DECODER_PRIMARY_H1{
	private:
		UInt_t NOpCh;

	public:
		TH1D **rPH;
		TH1D **QDC;
		TH1D **TDC;
		void Draw(DECODER_CANVAS *canvas, char flag=0x0F);

		DECODER_PRIMARY_H1(UInt_t _NOpCh, int QDCrange, int waveform_length, float sampletime){
			NOpCh=_NOpCh;
			rPH = new TH1D*[NOpCh];		QDC = new TH1D*[NOpCh];	TDC = new TH1D*[NOpCh];
			for (int ich=0; ich<NOpCh; ich++){
				rPH[ich] = new TH1D(Form("h1_rPH[%d]",ich),Form("h1_rPH[%d];rPH;count/rPH",ich), 4096, 0, 4096);
				QDC[ich] = new TH1D(Form("h1_QDC[%d]",ich),Form("h1_QDC[%d];QDC;count/%dQDC", ich, QDCrange/2000), 2000, 0, QDCrange);
				TDC[ich] = new TH1D(Form("h1_TDC[%d]",ich),Form("h1_TDC[%d];TDC;count/%.1fns",ich, sampletime), waveform_length, 0, waveform_length*sampletime);
			}
		}
		~DECODER_PRIMARY_H1(){
		}

};
inline void DECODER_PRIMARY_H1::Draw(DECODER_CANVAS *canvas, char flag=0x0F){
	if(!(flag & 0x01)) return;
	for (int ich=0; ich<NOpCh; ich++){
		canvas->Canvas_Decoding->cd(1+ich+0*NOpCh)->SetGridx();	rPH[ich]->Draw();
		canvas->Canvas_Decoding->cd(1+ich+1*NOpCh)->SetGridx();	QDC[ich]->Draw();
		canvas->Canvas_Decoding->cd(1+ich+2*NOpCh)->SetGridx();	TDC[ich]->Draw();
	}
	return;
}
//////////////////////////////////////////

class DECODER_QDC_H1_COIN{
	private:
		double binsize=1;
		int nbin=4000;
		int ich; 
		UInt_t NOpCh;
	public:
		TH1D **coin;
		void Draw(DECODER_CANVAS *canvas, char flag=0x0F);


		DECODER_QDC_H1_COIN(UInt_t _NOpCh, int _ich, int QDCrange){
			NOpCh=_NOpCh; ich=_ich;
			binsize=(double)QDCrange/(double)nbin;
			coin = new TH1D*[NOpCh];
			for (int jch=0; jch<NOpCh; jch++){
				coin[jch] = new TH1D(Form("QDC_ch%d_coin%d", ich,jch),Form("QDC_ch%d_coin%d;QDC;count/%.0fQDC", ich,jch, binsize), nbin,0,QDCrange);
			}
			
		}
		~DECODER_QDC_H1_COIN(){
		}
};
inline void DECODER_QDC_H1_COIN::Draw(DECODER_CANVAS *canvas, char flag=0x0F){
	if(!(flag & 0x02)) return;
	for (int jch=0; jch<NOpCh; jch++){
		canvas->Coincidence[0]->cd(1+ich+jch*NOpCh)->SetLogy();
		canvas->Coincidence[0]->cd(1+ich+jch*NOpCh)->SetGridx();
		coin[jch]->Draw();
	}
	return;
}
//////////////////////////////////////////////

class DECODER_QDC_H2{
	private:
		double binsize;
		int nbin=4000;
		int ich; 
		UInt_t NOpCh;
	public:
		TH2D **yaxis;
		void Draw(DECODER_CANVAS *canvas, char flag=0x0F);

		DECODER_QDC_H2(UInt_t _NOpCh, int _ich, int QDCrange){
			NOpCh=_NOpCh; ich=_ich;
			binsize=(double)QDCrange/(double)nbin;
			yaxis = new TH2D*[NOpCh];
			for (int jch=0; jch<NOpCh; jch++){
				yaxis[jch] = new TH2D(Form("QDC_xch%d_QDC_ych%d", ich,jch),Form("QDC_xch%d_QDC_ych%d;QDC_ch%d;QDC%d", ich,jch,ich,jch), 4000,0,QDCrange, 4000,0,QDCrange);
			}
			
		}
		~DECODER_QDC_H2(){
		}
};
inline void DECODER_QDC_H2::Draw(DECODER_CANVAS *canvas, char flag=0x0F){
	if(!(flag & 0x04)) return;
	for (int jch=0; jch<NOpCh; jch++){
		canvas->Coincidence[1]->cd(1+ich+jch*NOpCh)->SetLogz();
		canvas->Coincidence[1]->cd(1+ich+jch*NOpCh)->SetGridx();
		canvas->Coincidence[1]->cd(1+ich+jch*NOpCh)->SetGridy();
		yaxis[jch]->Draw("COLZ");
	}
	return;
}
///////////////////////////////////////////////////////


class DECODER_TDCDIF_H1{
	private:
		double binsize=0.02; // (ns)
		int nbin;
		int ich; 
		UInt_t NOpCh;
	public:
		TH1D **startch;
		void Draw(DECODER_CANVAS *canvas, char flag=0x0F);

		DECODER_TDCDIF_H1(UInt_t _NOpCh, int _ich, double T_lb=-3, double T_ub=3){
			NOpCh=_NOpCh; ich=_ich;
			nbin=(int) ((T_ub-T_lb)/binsize);
			startch = new TH1D*[NOpCh];
			for (int jch=0; jch<NOpCh; jch++){
				startch[jch] = new TH1D(Form("TDCdif_stopch%d-startch%d", ich,jch),Form("TDC_stopch%d-startch%d;TDC dif(ns);count/%.3fns", ich,jch, binsize), nbin, T_lb, T_ub);
			}
		}
		~DECODER_TDCDIF_H1(){
		}
};
inline void DECODER_TDCDIF_H1::Draw(DECODER_CANVAS *canvas, char flag=0x0F){
	if(!(flag & 0x08)) return;
	for (int jch=0; jch<NOpCh; jch++){
		canvas->Coincidence[2]->cd(1+ich+jch*NOpCh)->SetLogy();
		canvas->Coincidence[2]->cd(1+ich+jch*NOpCh)->SetGridx();
		startch[jch]->Draw();
	}
	return;
}




///////////////////////////////////////////////









