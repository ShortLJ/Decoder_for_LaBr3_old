void TDCcf_analysis(){
	TFile *file = new TFile("Ba133_att_run01_final.root","read");
	TTree *tree = (TTree*) file->Get("wavedata");
	const Int_t Entries = tree->GetEntries();
	printf("%d Entries\n", Entries);

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

//Ba133
//	const double ch1Ecut_mean = 1.08998e+05;
//	const double ch1Ecut_FWHM = 3.44549e+03*2.35;
//	const double ch0Ecut_mean = 9.99039e+04;
//	const double ch0Ecut_FWHM = 2.50430e+03*2.35;



	const double ch1Ecut_mean = 175000;
	const double ch1Ecut_FWHM = 350000;
	const double ch0Ecut_mean = 175000;
	const double ch0Ecut_FWHM = 350000;

//	const double ch1Ecut_mean = 1.16685e+05;
//	const double ch1Ecut_FWHM = 3.96374e+03*2.35;
//	const double ch0Ecut_mean = 175000;
//	const double ch0Ecut_FWHM = 350000;






	UInt_t TrigTimeTag;
	Float_t TDC[2];
	Float_t QDC[2];
	Float_t rPH[2];
	tree -> SetBranchAddress("TrigTimeTag", &TrigTimeTag);
	tree -> SetBranchAddress("TDC", TDC);
	tree -> SetBranchAddress("QDC", QDC);
	tree -> SetBranchAddress("rPH", rPH);

	TH2D *h2_QDCch1_QDCch0 = new TH2D("Ba133 Energy(ch1) by Energy(ch0)", "Ba133 Energy(ch1) by Energy(ch0);Energy,ch0 (keV);Energy,ch1 (keV)",3500,0,350000,3500,0,350000);
	TH1D *h1_QDCch0 = new TH1D("QDCch0","QDCch0;Energy (keV);count/keV",3500,0,350000);
	TH1D *h1_QDCch1	= new TH1D("ch1 Energy range used", "ch1 Energy range used;energy(keV);count/keV", 3500, 0, 350000);
	TH1D *h1_TDCdif	= new TH1D(Form("TDCcf ch0-ch1,ch1 energycut %.0fkeV",ch1Ecut_mean), Form("TDCcf dif(ch1 %.0fkeV);TDC ch0-ch1(ns);count/0.05ns",ch1Ecut_mean),500,-5,5);
	TH2D *h2_TDCdif_QDCch0 = new TH2D(Form("Ba133 TDCcf dif(ch0-ch1) by Energy(ch0), (ch1 %.0fkeV)",ch1Ecut_mean),Form("Ba133 TDCcf dif by ch0 Energy (ch1 %.0fkeV);energy,ch0(keV);TDC dif, ch0-ch1(ns)",ch1Ecut_mean),3500,0,350000,600,-15,15);
	TH2D *h2_TDCdif_QDCch0_nocut = new TH2D("Ba133 TDCcf dif(ch0-ch1) by Energy(ch0),nocut ","Ba133 TDCcf dif by ch0 Energy, nocut;energy,ch0(keV);TDC dif, ch0-ch1(ns)",3500,0,350000, 600, -15, 15);


	ULong64_t ient=0;

	for (ient=0; ient<Entries; ient++){
		tree->GetEntry(ient);
//		if(rPH[0]>100 && rPH[1]>100){
//if(TDC[0]-TDC[1]>0.02){continue;}
//if(TDC[1]-TDC[0]>0){continue;}
if(TDC[0]>2){continue;}
if(TDC[1]>2){continue;}
		if((ch1Ecut_mean-ch1Ecut_FWHM/2)<QDC[1] && QDC[1]<(ch1Ecut_mean+ch1Ecut_FWHM/2)){
			h1_QDCch0->Fill(QDC[0]);
			h2_TDCdif_QDCch0->Fill(QDC[0], TDC[0] - TDC[1]);
		if((ch0Ecut_mean-ch0Ecut_FWHM/2)<QDC[0] && QDC[0]<(ch0Ecut_mean+ch0Ecut_FWHM/2)){
			h1_TDCdif->Fill(TDC[0]-TDC[1]);
		}
		}
		h2_TDCdif_QDCch0_nocut->Fill(QDC[0],TDC[0]-TDC[1]);		// test  - double structure
		h2_QDCch1_QDCch0 -> Fill(QDC[0],QDC[1]);
		//h1_QDCch0->Fill(QDC[0]);
		h1_QDCch1->Fill(QDC[1]);

//		h4->Fill(QDC[1],TDC[0]-TDC[1]);
//		h2_TDCdif_Edif->Fill(QDC[0]-QDC[1],TDC[0]-TDC[1]);
	
//		if(ient == 1000000||ient == 2000000||ient==3000000||ient==4000000||ient==5000000){
			//checking
//		}
	}

	TCanvas *c1 = new TCanvas("Ba133 TH1", "Ba133 TH1", 1800,1000);
//	gStyle->SetOptStat(0);
	c1->Divide(3,2);
	
	c1->cd(1)->SetLogz();
	h2_QDCch1_QDCch0->Draw("COLZ");
	TLine *lEEvl = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,0,ch0Ecut_mean-ch0Ecut_FWHM/2,350000);
	TLine *lEEvu = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,0,ch0Ecut_mean+ch0Ecut_FWHM/2,350000);
	TLine *lEEhl = new TLine(0,ch1Ecut_mean-ch1Ecut_FWHM/2,350000,ch1Ecut_mean-ch1Ecut_FWHM/2);
	TLine *lEEhu = new TLine(0,ch1Ecut_mean+ch1Ecut_FWHM/2,350000,ch1Ecut_mean+ch1Ecut_FWHM/2);
	lEEvl->Draw("same");
	lEEvu->Draw("same");
	lEEhl->Draw("same");
	lEEhu->Draw("same");

	c1->cd(2)->SetLogy();
	h1_QDCch0->Draw();
	TLine *lE0vl = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,0,ch0Ecut_mean-ch0Ecut_FWHM/2,100000);
	TLine *lE0vu = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,0,ch0Ecut_mean+ch0Ecut_FWHM/2,100000);
	lE0vl->Draw("same");
	lE0vu->Draw("same");

	c1->cd(3)->SetLogy();
	h1_QDCch1->Draw();
	TLine *lE1vl = new TLine(ch1Ecut_mean-ch1Ecut_FWHM/2,0,ch1Ecut_mean-ch1Ecut_FWHM/2,100000);
	TLine *lE1vu = new TLine(ch1Ecut_mean+ch1Ecut_FWHM/2,0,ch1Ecut_mean+ch1Ecut_FWHM/2,100000);
	lE1vl->Draw("same");
	lE1vu->Draw("same");

	c1->cd(4);
	h1_TDCdif->Draw();
//	TF1 *fftn = new TF1("fftn", "[3]*exp(-[0]*((x-[2])-0.5*[0]*[1]*[1]))*TMath::Erfc(((x-[2])-[0]*[1]*[1])/(TMath::Sqrt(2)*[1]))",-1,1);
//	TF1 *fftn = new TF1("fftn", "[3]*exp([0]*([0]*[1]*[1]-2*(x-[2]))/2)*TMath::Erfc(([0]*[1]*[1]-(x-[2]))/(TMath::Sqrt(2)*[1]))",-1,1);
	TF1 *fftn = new TF1("fftn", "[3]*exp([0]*([0]*[1]*[1]-2*(x-[2]))/2)*TMath::Erfc(((x-[2])-[0]*[1]*[1])/(TMath::Sqrt(2)*[1]))",-1,1);
//	fftn->SetParLimits(1,0.01,0.3);
	fftn->SetParName(0,"lambda");
	fftn->SetParName(1,"sigma");
	fftn->SetParName(2,"mu");
	fftn->SetParName(3,"amp");
	h1_TDCdif->Fit(fftn,"M", "", -1, 1);

	c1->cd(5)->SetLogz();
	h2_TDCdif_QDCch0->Draw("COLZ");
	TLine *lE0vl2 = new TLine(ch0Ecut_mean-ch0Ecut_FWHM/2,-15,ch0Ecut_mean-ch0Ecut_FWHM/2,15);
	TLine *lE0vu2 = new TLine(ch0Ecut_mean+ch0Ecut_FWHM/2,-15,ch0Ecut_mean+ch0Ecut_FWHM/2,15);
	lE0vl2 -> Draw("same");
	lE0vu2 -> Draw("same");

	c1->cd(6)->SetLogz();
	h2_TDCdif_QDCch0_nocut->Draw("COLZ");
//	h2_TDCdif_Edif->Draw("COLZ");





}
