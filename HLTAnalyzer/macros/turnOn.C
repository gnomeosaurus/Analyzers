//calculate rate vs cut for scouting triggers
void setGlobalStyle()
{
  // For the statistics box:                                                                                                                                                                      
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");                                                                                                               
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  tdrStyle->SetPadGridX(true);
  tdrStyle->SetPadGridY(true);

  tdrStyle->SetMarkerSize(0.5);
}

inline
float f_deltaRR (float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = std::abs(phi1-phi2); if (dphi>float(M_PI)) dphi-=float(2*M_PI);  
    return deta*deta + dphi*dphi;
  }
inline
float f_deltaPhi (float phi1, float phi2)
{
  float dphi = std::abs(phi1-phi2); if (dphi>float(M_PI)) dphi-=float(2*M_PI);
  return dphi;
}


int turnOn(std::string fileName)
{
  setGlobalStyle();
  
  //###############################
  unsigned int HT250Calo  = 9; //index of DST_HT250_CaloScouting_v   old:1 ref:9
  float HT250Calo_rate = 1928;

  unsigned int HT410PF = 7; //index of DST_HT410_PFScouting_v   old:3 ref:7
  float HT410PF_rate = 294;

  unsigned int MJJ200Calo = 12; //index of DST_DiCaloWideJetMass200_CaloScouting_v
  
  unsigned int HTT200 = 0; //index if L1_HTT200
  unsigned int HTT240 = 1; //index if L1_HTT240
  unsigned int HTT270 = 2; //index if L1_HTT270
  unsigned int HTT280 = 3; //index if L1_HTT280
  unsigned int DoubleJetC100 = 7; //index if L1_DoubleJetC100
  unsigned int DoubleJetC112 = 8; //index if L1_DoubleJetC112
  unsigned int DoubleIsoTau28er = 11; //index if L1_DoubleJetC112

  unsigned int L1scenario = HTT240;
  //###############################

  TChain* tt = new TChain("MyAnalysis/HLTree");
  tt->Add(fileName.c_str());

  //set branches
  TBranch* b_lumi;
  
  TBranch* b_caloMjj;
  TBranch* b_PFMjj;
  TBranch* b_l1Mjj;

  TBranch* b_hltAccept;
  TBranch* b_l1Accept;
  TBranch* b_l1Names;

  TBranch* b_caloJet1Pt;
  TBranch* b_caloJet2Pt;
  TBranch* b_caloJet1Eta;
  TBranch* b_caloJet2Eta;
  TBranch* b_caloJet1Phi;
  TBranch* b_caloJet2Phi;
  TBranch* b_caloDeltaEta;

  TBranch* b_PFJet1Pt;
  TBranch* b_PFJet2Pt;
  TBranch* b_PFJet1Eta;
  TBranch* b_PFJet2Eta;
  TBranch* b_PFJet1Phi;
  TBranch* b_PFJet2Phi;
  TBranch* b_PFDeltaEta;

  TBranch* b_l1Jet1Pt;
  TBranch* b_l1Jet2Pt;
  TBranch* b_l1Jet1Eta;
  TBranch* b_l1Jet2Eta;
  TBranch* b_l1Jet1Phi;
  TBranch* b_l1Jet2Phi;
  TBranch* b_l1DeltaEta;

  TBranch* b_l1JetPt;
  TBranch* b_l1JetEta;
  TBranch* b_l1JetPhi;

  int lumi = 0;
  float caloMjj = 0;
  float PFMjj = 0;
  float l1Mjj = 0;

  float caloJet1Pt_ = 0;
  float caloJet2Pt_ = 0;
  float caloJet1Eta_ = -999;
  float caloJet2Eta_ = -999;
  float caloJet1Phi_ = -999;
  float caloJet2Phi_ = -999;
  float caloDeltaEta_ = -999;

  float PFJet1Pt_ = 0;
  float PFJet2Pt_ = 0;
  float PFJet1Eta_ = -999;
  float PFJet2Eta_ = -999;
  float PFJet1Phi_ = -999;
  float PFJet2Phi_ = -999;
  float PFDeltaEta_ = -999;

  float l1Jet1Pt_ = 0;
  float l1Jet2Pt_ = 0;
  float l1Jet1Eta_ = -999;
  float l1Jet2Eta_ = -999;
  float l1Jet1Phi_ = -999;
  float l1Jet2Phi_ = -999;
  float l1DeltaEta_ = -999;

  std::vector<float>* l1JetPt_ = 0;
  std::vector<float>* l1JetEta_ = 0;
  std::vector<float>* l1JetPhi_ = 0;
  
  std::vector<int>* hltAccept = 0;
  std::vector<int>* l1Accept = 0;
  std::vector<string>* l1Names = 0;

  tt->SetBranchAddress("lumi", &lumi, &b_lumi);

  tt->SetBranchAddress("caloMjj", &caloMjj, &b_caloMjj);
  tt->SetBranchAddress("PFMjj", &PFMjj, &b_PFMjj);
  tt->SetBranchAddress("l1Mjj", &l1Mjj, &b_l1Mjj);

  tt->SetBranchAddress("caloJet1Pt", &caloJet1Pt_, &b_caloJet1Pt);
  tt->SetBranchAddress("caloJet2Pt", &caloJet2Pt_, &b_caloJet2Pt);
  tt->SetBranchAddress("caloJet1Eta", &caloJet1Eta_, &b_caloJet1Eta);
  tt->SetBranchAddress("caloJet2Eta", &caloJet2Eta_, &b_caloJet2Eta);
  tt->SetBranchAddress("caloJet1Phi", &caloJet1Phi_, &b_caloJet1Phi);
  tt->SetBranchAddress("caloJet2Phi", &caloJet2Phi_, &b_caloJet2Phi);
  tt->SetBranchAddress("caloDeltaEta", &caloDeltaEta_, &b_caloDeltaEta);

  tt->SetBranchAddress("PFJet1Pt", &PFJet1Pt_, &b_PFJet1Pt);
  tt->SetBranchAddress("PFJet2Pt", &PFJet2Pt_, &b_PFJet2Pt);
  tt->SetBranchAddress("PFJet1Eta", &PFJet1Eta_, &b_PFJet1Eta);
  tt->SetBranchAddress("PFJet2Eta", &PFJet2Eta_, &b_PFJet2Eta);
  tt->SetBranchAddress("PFJet1Phi", &PFJet1Phi_, &b_PFJet1Phi);
  tt->SetBranchAddress("PFJet2Phi", &PFJet2Phi_, &b_PFJet2Phi);
  tt->SetBranchAddress("PFDeltaEta", &PFDeltaEta_, &b_PFDeltaEta);

  tt->SetBranchAddress("l1Jet1Pt", &l1Jet1Pt_, &b_l1Jet1Pt);
  tt->SetBranchAddress("l1Jet2Pt", &l1Jet2Pt_, &b_l1Jet2Pt);
  tt->SetBranchAddress("l1Jet1Eta", &l1Jet1Eta_, &b_l1Jet1Eta);
  tt->SetBranchAddress("l1Jet2Eta", &l1Jet2Eta_, &b_l1Jet2Eta);
  tt->SetBranchAddress("l1Jet1Phi", &l1Jet1Phi_, &b_l1Jet1Phi);
  tt->SetBranchAddress("l1Jet2Phi", &l1Jet2Phi_, &b_l1Jet2Phi);
  tt->SetBranchAddress("l1DeltaEta", &l1DeltaEta_, &b_l1DeltaEta);

  tt->SetBranchAddress("l1JetPt", &l1JetPt_, &b_l1JetPt);
  tt->SetBranchAddress("l1JetEta", &l1JetEta_, &b_l1JetEta);
  tt->SetBranchAddress("l1JetPhi", &l1JetPhi_, &b_l1JetPhi);

  tt->SetBranchAddress("hltAccept", &hltAccept, &b_hltAccept);
  tt->SetBranchAddress("l1Accept", &l1Accept, &b_l1Accept);
  tt->SetBranchAddress("l1Names", &l1Names, &b_l1Names);

  int nentries = tt->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;

  //book graphs and plots
  float min = 0.;
  float max = 1000.;
  int nBins = 20;

  TF1* f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])-[0]*TMath::Erf((-x-[1])/[2])",min,max);
  f1->SetParameters(0.5,350,40);  
  f1->FixParameter(0,0.5);
  f1->SetLineWidth(2.);
  f1->SetLineColor(kRed);

  TF1* f2 = (TF1*)f1->Clone("f2");
  f2->SetParameters(0.5,150,10);
  f2->SetLineColor(kBlack);

  TH1F* caloMjjSpectrum = new TH1F("caloMjjSpectrum","caloMjjSpectrum",nBins,min,max);
  
  TEfficiency* caloMjj_eff = new TEfficiency("caloMjj_eff","caloMjj_eff",nBins,min,max);
  caloMjj_eff->SetLineWidth(2);
  caloMjj_eff->SetTitle("turnOn;Mjj [GeV]");
  TEfficiency* caloHT250_eff = new TEfficiency("caloHT250_eff","caloHT250_eff",nBins,min,max);
  caloHT250_eff->SetMarkerColor(kBlue);
  caloHT250_eff->SetLineColor(kBlue);
  TEfficiency* HTT240_eff = new TEfficiency("HTT240_eff","HTT240_eff",nBins,min,max);
  HTT240_eff->SetMarkerColor(kGreen+2);
  HTT240_eff->SetLineColor(kGreen+2);
  TEfficiency* l1Mjj_eff = new TEfficiency("l1Mjj_eff","l1Mjj_eff",nBins,min,max);
  l1Mjj_eff->SetLineWidth(2);
  l1Mjj_eff->SetTitle("turnOn;Mjj [GeV]");
  l1Mjj_eff->SetMarkerColor(kOrange+1);
  l1Mjj_eff->SetLineColor(kOrange+1);

  TH1F* mjj_res = new TH1F("mjj_res","mjj_res",1200,-3.,3.);
  mjj_res->GetXaxis()->SetTitle("2*(l1Mjj-caloMjj)/(l1Mjj+caloMjj)");
  TH1F* mjj_diff = new TH1F("mjj_diff","mjj_diff",1500,-300.,300.);
  mjj_diff->GetXaxis()->SetTitle("l1Mjj-caloMjj");

  TH1F* pt1_res = new TH1F("pt1_res","pt1_res",1200,-3.,3.);
  pt1_res->GetXaxis()->SetTitle("2*(l1Jet1Pt-caloJet1Pt)/(l1Jet1Pt+caloJet1Pt)");
  TH1F* pt2_res = new TH1F("pt2_res","pt2_res",1200,-3.,3.);
  pt2_res->GetXaxis()->SetTitle("2*(l1Jet2Pt-caloJet2Pt)/(l1Jet2Pt+caloJet2Pt)");

  TH1F* deltaR1 = new TH1F("deltaR1","deltaR1",2000,0.,10.);
  deltaR1->GetXaxis()->SetTitle("deltaR1");
  TH1F* deltaR2 = new TH1F("deltaR2","deltaR2",2000,0.,10.);
  deltaR2->GetXaxis()->SetTitle("deltaR2");
  TH1F* deltaPhi = new TH1F("deltaPhi","deltaPhi",350,0.,3.5);
  deltaPhi->GetXaxis()->SetTitle("deltaPhi");

  TH2F* deltaR1_vsMjj = new TH2F("deltaR1_vsMjj","deltaR1_vsMjj",1000,0.,1000.,2000,0.,10.);
  deltaR1_vsMjj->GetXaxis()->SetTitle("Mjj [GeV]");
  deltaR1_vsMjj->GetYaxis()->SetTitle("deltR");
  TH2F* deltaPhi_vsMjj = new TH2F("deltaPhi_vsMjj","deltaPhi_vsMjj",1000,0.,1000.,350,0.,3.5);
  deltaPhi_vsMjj->GetXaxis()->SetTitle("Mjj [GeV]");
  deltaPhi_vsMjj->GetYaxis()->SetTitle("deltPhi");
  TH2F* deltaMjj_vsMjj = new TH2F("deltaMjj_vsMjj","deltaMjj_vsMjj",1000,0.,1000.,600,-300.,300.);
  deltaMjj_vsMjj->GetXaxis()->SetTitle("Mjj [GeV]");
  deltaMjj_vsMjj->GetYaxis()->SetTitle("deltMjj");
  
  TH1F* l1 = new TH1F("l1","l1",14,0.,14.);
  TH1F* l2 = new TH1F("l2","l2",14,0.,14.);
  
  //loop
  for (Long64_t jentry=0; jentry<nentries;++jentry)
    {
      tt->GetEntry(jentry);

      //l1 and hlt rates
      for(unsigned int ii=0; ii<l1Names->size(); ++ii)
	if (l1Accept->at(ii)==1)
	  l1->Fill(ii);


      float l1Jet1MatchedPt = 0.;
      float l1Jet1MatchedEta = -999.;
      float l1Jet1MatchedPhi = -999.;
      float l1Jet2MatchedPt = 0.;
      float l1Jet2MatchedEta = -999.;
      float l1Jet2MatchedPhi = -999.;

      float l1CaloJet1DeltaR = 10000.;
      float l1CaloJet2DeltaR = 10000.;

      for (unsigned int jitr = 0; jitr < l1JetPt_->size(); ++jitr)
	{
	  if(f_deltaRR(l1JetEta_->at(jitr),l1JetPhi_->at(jitr),caloJet1Eta_,caloJet1Phi_) < l1CaloJet1DeltaR)
	    {
	      l1Jet1MatchedPt = l1JetPt_->at(jitr);
	      l1Jet1MatchedEta = l1JetEta_->at(jitr);
	      l1Jet1MatchedPhi = l1JetPhi_->at(jitr);
	    }
	  if(f_deltaRR(l1JetEta_->at(jitr),l1JetPhi_->at(jitr),caloJet2Eta_,caloJet2Phi_) < l1CaloJet2DeltaR)
	    {
	      l1Jet2MatchedPt = l1JetPt_->at(jitr);
	      l1Jet2MatchedEta = l1JetEta_->at(jitr);
	      l1Jet2MatchedPhi = l1JetPhi_->at(jitr);
	    }
	}


      float l1MatchedMjj = sqrt(2*l1Jet1MatchedPt*l1Jet2MatchedPt*(cosh(l1Jet1MatchedEta-l1Jet2MatchedEta) - cos(l1Jet1MatchedPhi-l1Jet2MatchedPhi)));

      
      //resolution and comparison L1 VS Calo (upstream analysis cuts)
      if(l1Mjj > 0 && caloMjj > 0)
	{
	  mjj_res->Fill(2*(l1Mjj-caloMjj)/(l1Mjj+caloMjj));
	  mjj_diff->Fill(l1Mjj-caloMjj);
	  pt1_res->Fill(2*(l1Jet1Pt_-caloJet1Pt_)/(l1Jet1Pt_+caloJet1Pt_));
	  pt2_res->Fill(2*(l1Jet2Pt_-caloJet2Pt_)/(l1Jet2Pt_+caloJet2Pt_));

	  float deltaR1_ = sqrt(f_deltaRR(l1Jet1Eta_,l1Jet1Phi_,caloJet1Eta_,caloJet1Phi_));
	  deltaR1->Fill(deltaR1_);
	  float deltaR2_ = sqrt(f_deltaRR(l1Jet2Eta_,l1Jet2Phi_,caloJet2Eta_,caloJet2Phi_));
	  deltaR2->Fill(deltaR2_);
	  deltaPhi->Fill(f_deltaPhi(l1Jet1Phi_,l1Jet2Phi_));

	  deltaR1_vsMjj->Fill(caloMjj,deltaR1_);
	  deltaPhi_vsMjj->Fill(caloMjj,f_deltaPhi(l1Jet1Phi_,l1Jet2Phi_));
	  deltaMjj_vsMjj->Fill(caloMjj,l1Mjj-caloMjj);
	}



      bool l1Pass = (l1Mjj > 150. &&
		     
      		     l1Jet1Pt_ > 15. &&
      		     l1Jet2Pt_ > 15. &&
      		     fabs(l1Jet1Eta_) < 5.0 &&
      		     fabs(l1Jet2Eta_) < 5.0 &&
      		     l1DeltaEta_ < 5.0
      		     );

      bool l1MatchedPass = (l1MatchedMjj > 150. &&
			    
			    l1Jet1MatchedPt > 15. &&
			    l1Jet2MatchedPt > 15. &&
			    fabs(l1Jet1MatchedEta) < 5.0 &&
			    fabs(l1Jet2MatchedEta) < 5.0 &&
			    fabs(l1Jet1MatchedEta-l1Jet2MatchedEta) < 5.0
			    );

      
      
      
      //analysis cuts needed to compare to the analysis
      //calo analysis
      if (caloJet1Pt_ > 60. &&
	  caloJet2Pt_ > 30. &&
	  fabs(caloJet1Eta_) < 2.5 &&
	  fabs(caloJet2Eta_) < 2.5 &&
	  caloDeltaEta_ < 1.3
	  )//&& l1Mjj > 0) //CAREFUL HERE
	{
	  caloMjjSpectrum->Fill(caloMjj);
	  caloHT250_eff->Fill((hltAccept->at(HT250Calo)==1 && l1Accept->at(L1scenario)==1), caloMjj);

	  //caloMjj_eff->Fill((caloMjj>200 && l1Accept->at(L1scenario)==1) || hltAccept->at(HT250Calo)==1, caloMjj);
	  HTT240_eff->Fill(l1Accept->at(HTT240)==1, caloMjj);

	  
	  // l1Mjj_eff->Fill(l1Pass, caloMjj);
	  // caloMjj_eff->Fill((caloMjj>200 && l1Pass), caloMjj);
	  l1Mjj_eff->Fill(l1MatchedPass, caloMjj);
	  caloMjj_eff->Fill((caloMjj>200 && l1MatchedPass), caloMjj);


	  //l1 and hlt rates
	  for(unsigned int ii=0; ii<l1Names->size(); ++ii)
	    if (l1Accept->at(ii)==1)
	      l2->Fill(ii);
	}

    }

  // caloMjj_eff->Fit(f2,"r");

  caloMjjSpectrum->Scale(1./caloMjjSpectrum->GetBinContent(caloMjjSpectrum->GetMaximumBin()));
			      
  
  TLegend* leg0 = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg0->AddEntry(caloMjj_eff,"caloMjj","P");
  leg0->AddEntry(caloHT250_eff,"HT250_Calo","P");
  leg0->AddEntry(HTT240_eff,"HTT240","P");
  leg0->AddEntry(l1Mjj_eff,"l1Mjj","P");

  TCanvas* c1 = new TCanvas();
  caloMjj_eff->Draw();
  caloHT250_eff->Draw("sames");
  HTT240_eff->Draw("sames");
  l1Mjj_eff->Draw("sames");
  //caloMjjSpectrum->Draw("L,sames");
  leg0->Draw("sames");

  TCanvas* c2 = new TCanvas();
  for(unsigned int ii=0; ii<l1Names->size(); ++ii)
    l1->GetXaxis()->SetBinLabel(ii+1,l1Names->at(ii).c_str());
  //l1->GetYaxis()->SetTitle("L1 Rate @4E33 [Hz]");
  l1->SetMaximum(l1->GetMaximum()+200);
  l2->SetLineColor(kRed);
  
  l1->Draw();
  l2->Draw("same");

  c2->SetLogy();
  c2->Update();

  TCanvas* c3 = new TCanvas();
  mjj_res->Draw();
  TCanvas* c31 = new TCanvas();
  mjj_diff->Draw();
  TCanvas* c4 = new TCanvas();
  deltaR1->Draw();
  TCanvas* c41 = new TCanvas();
  deltaR2->Draw();
  TCanvas* c42 = new TCanvas();
  deltaPhi->Draw();
  TCanvas* c5 = new TCanvas();
  pt1_res->Draw();
  TCanvas* c6 = new TCanvas();
  pt2_res->Draw();
  TCanvas* c7 = new TCanvas();
  deltaR1_vsMjj->Draw("colz");
  TCanvas* c8 = new TCanvas();
  deltaPhi_vsMjj->Draw("colz");
  TCanvas* c9 = new TCanvas();
  deltaMjj_vsMjj->Draw("colz");
  
  c1->Print("turnOn.pdf","pdf");
  c2->Print("L1Rates.pdf","pdf");
  c3->Print("mjj_res.pdf","pdf");
  c31->Print("mjj_diff.pdf","pdf");
  c4->Print("deltaR1.pdf","pdf");
  c41->Print("deltaR2.pdf","pdf");
  c5->Print("pt1_res.pdf","pdf");
  c6->Print("pt2_res.pdf","pdf");
  c7->Print("deltaR1_vsMjj.pdf","pdf");
  c8->Print("deltaPhi_vsMjj.pdf","pdf");
  c9->Print("deltaMjj_vsMjj.pdf","pdf");
  
  return 0;


}
