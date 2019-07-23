void gainshiftcheck (){
	TFile *file = new TFile("Cs137_att_run01_final.root","read");
	TTree *tree = (TTree*) file->Get("wavedata");
	const Int_t Entries = tree->GetEntries();
	const int Nprd=10;
	const Int_t entdvded = (int) Entries/Nprd;
	printf("%d Entries of whole, %d each period \n", Entries, entdvded);


//Na22
//	const double ch1Ecut_mean = 8.85258e+04;
//	const double ch1Ecut_FWHM = 2.41257e+03*4;
//	const double ch0Ecut_mean = 7.96995e+04;
//	const double ch0Ecut_FWHM = 1.98168e+03*4;

//Co60
//	const double ch1Ecut_mean = 2.26379e+05;
//	const double ch1Ecut_FWHM = 3.46623e+03*2.35;
//	const double ch0Ecut_mean = 1.93290e+05;
//	const double ch0Ecut_FWHM = 3.08465e+03*2.35;


	const double ch1Ecut_mean = 175000;
	const double ch1Ecut_FWHM = 350000;
	const double ch0Ecut_mean = 175000;
	const double ch0Ecut_FWHM = 350000;

//	const double ch1Ecut_mean = 1.16685e+05;
//	const double ch1Ecut_FWHM = 3.96374e+03*2.35;
//	const double ch0Ecut_mean = 175000;
//	const double ch0Ecut_FWHM = 350000;

//	const double ch1Ecut_mean = 75000;
//	const double ch1Ecut_FWHM = 24000;
//	const double ch0Ecut_mean = 81000;
//	const double ch0Ecut_FWHM = 11000;




	UInt_t TrigTimeTag;
	Float_t QDC[2];
	Float_t rPH[2];
	tree -> SetBranchAddress("TrigTimeTag", &TrigTimeTag);
	tree -> SetBranchAddress("QDC", QDC);
	tree -> SetBranchAddress("rPH", rPH);

//	TH2D *h2_QDCch1_QDCch0 = new TH2D("Co60 Energy(ch1) by Energy(ch0)", "Co60 Energy(ch1) by Energy(ch0);Energy,ch0 (keV);Energy,ch1 (keV)",3500,0,350000,3500,0,350000);
//	TH1D *h1_QDCch0 = new TH1D("QDCch0","QDCch0;Energy (keV);count/keV",3500,0,350000);
//	TH1D *h1_QDCch1	= new TH1D("ch1 Energy range used", "ch1 Energy range used;energy(keV);count/keV", 3500, 0, 350000);
//	TH1D *h1_TDCdif	= new TH1D(Form("TDCcf ch0-ch1,ch1 energycut %.0fkeV",ch1Ecut_mean), Form("TDCcf dif(ch1 %.0fkeV);TDC ch0-ch1(ns);count/0.05ns",ch1Ecut_mean),125,-2,2);
//	TH2D *h2_TDCdif_QDCch0 = new TH2D(Form("Co60 TDCcf dif(ch0-ch1) by Energy(ch0), (ch1 %.0fkeV)",ch1Ecut_mean),Form("Co60 TDCcf dif by ch0 Energy (ch1 %.0fkeV);energy,ch0(keV);TDC dif, ch0-ch1(ns)",ch1Ecut_mean),3500,0,350000,600,-15,15);
//	TH2D *h2_TDCdif_QDCch0_nocut = new TH2D("Co60 TDCcf dif(ch0-ch1) by Energy(ch0),nocut ","Co60 TDCcf dif by ch0 Energy, nocut;energy,ch0(keV);TDC dif, ch0-ch1(ns)",3500,0,350000, 600, -15, 15);

	TGraph *gr0 = new TGraph();
	TGraph *gr1 = new TGraph();


	TH1D *h1_QDCch0;
	TH1D *h1_QDCch1;

	TCanvas *c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	
	TF1 *fftn[2];
	Double_t lastfit[2][5] = {{2.68187e+02,7.14518e+04,2.70579e+03,9.54131e+01,-1.05260e-03},{1.90310e+02,7.57921e+04,2.64953e+03,3.15001e+01,6.00249e-06}};
	
	fftn[0] = new TF1("fftn0","gaus(0)+[3]+[4]*x", lastfit[0][1]-3*lastfit[0][2], lastfit[0][1]+3*lastfit[0][2]);
	fftn[1] = new TF1("fftn1","gaus(0)+[3]+[4]*x", lastfit[1][1]-3*lastfit[1][2], lastfit[1][1]+3*lastfit[1][2]);
	fftn[0]->SetParameters(lastfit[0]);
	fftn[1]->SetParameters(lastfit[1]);

	ULong64_t ient=0;
	int iprd = 0;

//	for (iprd=0; iprd<1; iprd++){
	for (iprd=0; iprd<Nprd; iprd++){
	h1_QDCch0 = new TH1D("QDCch0",Form("QDCch0,%d of %d;Energy (keV);count/keV",iprd,Nprd),3500,0,350000);
	h1_QDCch1 = new TH1D("QDCch1",Form("QDCch1,%d of %d;Energy (keV);count/keV",iprd,Nprd),3500,0,350000);
		for (ient=entdvded*iprd; ient<entdvded*(iprd+1); ient++){		//276463
			tree->GetEntry(ient);
			if(rPH[0]>100 && rPH[1]>100){
				c1->cd(1);
				h1_QDCch0->Fill(QDC[0]);
//				h1_QDCch0->Draw();
				c1->cd(2);
				h1_QDCch1->Fill(QDC[1]);
//				h1_QDCch1->Draw();
			}
		}
	h1_QDCch0->Fit("fftn0","MQ","",lastfit[0][1]-3*lastfit[0][2], lastfit[0][1]+3*lastfit[0][2]);
	h1_QDCch1->Fit("fftn1","MQ","",lastfit[1][1]-3*lastfit[1][2], lastfit[1][1]+3*lastfit[1][2]);
	gr0->SetPoint(iprd,iprd,fftn[0]->GetParameter(1));
	gr1->SetPoint(iprd,iprd,fftn[1]->GetParameter(1));
	fftn[0]->GetParameters(lastfit[0]);
	fftn[1]->GetParameters(lastfit[1]);
	printf("Eres ch0 = %.2f%%, ch1 = %.2f%%\n",lastfit[0][2]/lastfit[0][1]*235,lastfit[1][2]/lastfit[1][1]*235);
	
	}
	TCanvas *c2 = new TCanvas("a","d",1200,600);
	c2->Divide(2,1);
	c2->cd(1);
	gr0->Draw();
	gr0->GetYaxis()->SetRangeUser(70000,85000);
	gr0->SetMarkerStyle(8);
	c2->cd(2);
	gr1->Draw();
	gr1->GetYaxis()->SetRangeUser(70000,85000);
	gr1->SetMarkerStyle(8);



}
