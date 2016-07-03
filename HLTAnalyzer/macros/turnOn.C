//calculate rate vs cut for scouting triggers
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress (double percentage)
{
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

void setGlobalStyle()
{
  // For the statistics box:                                                                                                                                                                      
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");                                                                                                               
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

  bool matched = true;
  bool wideJets = false;
  //###############################

  TChain* tt = new TChain("MyAnalysis/HLTree");
  tt->Add(fileName.c_str());

  //set branches
  TBranch* b_lumi;
  
  TBranch* b_caloMjj;
  TBranch* b_PFMjj;
  TBranch* b_caloWMjj;
  TBranch* b_PFWMjj;
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
  TBranch* b_caloWJet1Pt;
  TBranch* b_caloWJet2Pt;
  TBranch* b_caloWJet1Eta;
  TBranch* b_caloWJet2Eta;
  TBranch* b_caloWJet1Phi;
  TBranch* b_caloWJet2Phi;
  TBranch* b_caloWDeltaEta;

  TBranch* b_PFJet1Pt;
  TBranch* b_PFJet2Pt;
  TBranch* b_PFJet1Eta;
  TBranch* b_PFJet2Eta;
  TBranch* b_PFJet1Phi;
  TBranch* b_PFJet2Phi;
  TBranch* b_PFDeltaEta;
  TBranch* b_PFWJet1Pt;
  TBranch* b_PFWJet2Pt;
  TBranch* b_PFWJet1Eta;
  TBranch* b_PFWJet2Eta;
  TBranch* b_PFWJet1Phi;
  TBranch* b_PFWJet2Phi;
  TBranch* b_PFWDeltaEta;

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

  int lumi_ = 0;
  float caloMjj_ = 0;
  float PFMjj_ = 0;
  float caloWMjj_ = 0;
  float PFWMjj_ = 0;
  float l1Mjj_ = 0;

  float caloJet1Pt_ = 0;
  float caloJet2Pt_ = 0;
  float caloJet1Eta_ = -999;
  float caloJet2Eta_ = -999;
  float caloJet1Phi_ = -999;
  float caloJet2Phi_ = -999;
  float caloDeltaEta_ = -999;
  float caloWJet1Pt_ = 0;
  float caloWJet2Pt_ = 0;
  float caloWJet1Eta_ = -999;
  float caloWJet2Eta_ = -999;
  float caloWJet1Phi_ = -999;
  float caloWJet2Phi_ = -999;
  float caloWDeltaEta_ = -999;

  float PFJet1Pt_ = 0;
  float PFJet2Pt_ = 0;
  float PFJet1Eta_ = -999;
  float PFJet2Eta_ = -999;
  float PFJet1Phi_ = -999;
  float PFJet2Phi_ = -999;
  float PFDeltaEta_ = -999;
  float PFWJet1Pt_ = 0;
  float PFWJet2Pt_ = 0;
  float PFWJet1Eta_ = -999;
  float PFWJet2Eta_ = -999;
  float PFWJet1Phi_ = -999;
  float PFWJet2Phi_ = -999;
  float PFWDeltaEta_ = -999;

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
  
  std::vector<int>* hltAccept_ = 0;
  std::vector<int>* l1Accept_ = 0;
  std::vector<string>* l1Names_ = 0;

  tt->SetBranchAddress("lumi", &lumi_, &b_lumi);

  tt->SetBranchAddress("caloMjj", &caloMjj_, &b_caloMjj);
  tt->SetBranchAddress("PFMjj", &PFMjj_, &b_PFMjj);
  tt->SetBranchAddress("caloWMjj", &caloWMjj_, &b_caloWMjj);
  tt->SetBranchAddress("PFWMjj", &PFWMjj_, &b_PFWMjj);
  tt->SetBranchAddress("l1Mjj", &l1Mjj_, &b_l1Mjj);

  tt->SetBranchAddress("caloJet1Pt", &caloJet1Pt_, &b_caloJet1Pt);
  tt->SetBranchAddress("caloJet2Pt", &caloJet2Pt_, &b_caloJet2Pt);
  tt->SetBranchAddress("caloJet1Eta", &caloJet1Eta_, &b_caloJet1Eta);
  tt->SetBranchAddress("caloJet2Eta", &caloJet2Eta_, &b_caloJet2Eta);
  tt->SetBranchAddress("caloJet1Phi", &caloJet1Phi_, &b_caloJet1Phi);
  tt->SetBranchAddress("caloJet2Phi", &caloJet2Phi_, &b_caloJet2Phi);
  tt->SetBranchAddress("caloDeltaEta", &caloDeltaEta_, &b_caloDeltaEta);
  tt->SetBranchAddress("caloWJet1Pt", &caloWJet1Pt_, &b_caloWJet1Pt);
  tt->SetBranchAddress("caloWJet2Pt", &caloWJet2Pt_, &b_caloWJet2Pt);
  tt->SetBranchAddress("caloWJet1Eta", &caloWJet1Eta_, &b_caloWJet1Eta);
  tt->SetBranchAddress("caloWJet2Eta", &caloWJet2Eta_, &b_caloWJet2Eta);
  tt->SetBranchAddress("caloWJet1Phi", &caloWJet1Phi_, &b_caloWJet1Phi);
  tt->SetBranchAddress("caloWJet2Phi", &caloWJet2Phi_, &b_caloWJet2Phi);
  tt->SetBranchAddress("caloWDeltaEta", &caloWDeltaEta_, &b_caloWDeltaEta);

  tt->SetBranchAddress("PFJet1Pt", &PFJet1Pt_, &b_PFJet1Pt);
  tt->SetBranchAddress("PFJet2Pt", &PFJet2Pt_, &b_PFJet2Pt);
  tt->SetBranchAddress("PFJet1Eta", &PFJet1Eta_, &b_PFJet1Eta);
  tt->SetBranchAddress("PFJet2Eta", &PFJet2Eta_, &b_PFJet2Eta);
  tt->SetBranchAddress("PFJet1Phi", &PFJet1Phi_, &b_PFJet1Phi);
  tt->SetBranchAddress("PFJet2Phi", &PFJet2Phi_, &b_PFJet2Phi);
  tt->SetBranchAddress("PFDeltaEta", &PFDeltaEta_, &b_PFDeltaEta);
  tt->SetBranchAddress("PFWJet1Pt", &PFWJet1Pt_, &b_PFWJet1Pt);
  tt->SetBranchAddress("PFWJet2Pt", &PFWJet2Pt_, &b_PFWJet2Pt);
  tt->SetBranchAddress("PFWJet1Eta", &PFWJet1Eta_, &b_PFWJet1Eta);
  tt->SetBranchAddress("PFWJet2Eta", &PFWJet2Eta_, &b_PFWJet2Eta);
  tt->SetBranchAddress("PFWJet1Phi", &PFWJet1Phi_, &b_PFWJet1Phi);
  tt->SetBranchAddress("PFWJet2Phi", &PFWJet2Phi_, &b_PFWJet2Phi);
  tt->SetBranchAddress("PFWDeltaEta", &PFWDeltaEta_, &b_PFWDeltaEta);

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

  tt->SetBranchAddress("hltAccept", &hltAccept_, &b_hltAccept);
  tt->SetBranchAddress("l1Accept", &l1Accept_, &b_l1Accept);
  tt->SetBranchAddress("l1Names", &l1Names_, &b_l1Names);

  int nentries = tt->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;

  //book graphs and plots
  float min = 0.;
  float max = 1000.;
  int nBins = 40;

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
  deltaMjj_vsMjj->GetYaxis()->SetTitle("deltaMjj (l1-calo)");
  
  TH1F* l1 = new TH1F("l1","l1",14,0.,14.);
  TH1F* l2 = new TH1F("l2","l2",14,0.,14.);
  
  //loop
  for (Long64_t jentry=0; jentry<nentries;++jentry)
    {
      tt->GetEntry(jentry);
      if(jentry%1000 == 0)
	printProgress((float)jentry/(float)nentries);
    
      //l1 and hlt rates
      for(unsigned int ii=0; ii<l1Names_->size(); ++ii)
	if (l1Accept_->at(ii)==1)
	  l1->Fill(ii);


      float l1Jet1MatchedPt = 0.;
      float l1Jet1MatchedEta = -999.;
      float l1Jet1MatchedPhi = -999.;
      float l1Jet2MatchedPt = 0.;
      float l1Jet2MatchedEta = -999.;
      float l1Jet2MatchedPhi = -999.;

      float l1CaloJet1DeltaR = 100000000.;
      float l1CaloJet2DeltaR = 100000000.;

      int firstItr = -1;
      int secondItr = -1;
      int secondItrTmp = -1;
      for (unsigned int jitr = 0; jitr < l1JetPt_->size(); ++jitr)
	{
	  float l1CaloJet1DeltaR_tmp = f_deltaRR(l1JetEta_->at(jitr),l1JetPhi_->at(jitr),caloJet1Eta_,caloJet1Phi_);
	  if(l1CaloJet1DeltaR_tmp < l1CaloJet1DeltaR)
	    {
	      l1Jet1MatchedPt = l1JetPt_->at(jitr);
	      l1Jet1MatchedEta = l1JetEta_->at(jitr);
	      l1Jet1MatchedPhi = l1JetPhi_->at(jitr);
	      l1CaloJet1DeltaR = l1CaloJet1DeltaR_tmp;
	      firstItr = jitr;
	    }
	  float l1CaloJet2DeltaR_tmp = f_deltaRR(l1JetEta_->at(jitr),l1JetPhi_->at(jitr),caloJet2Eta_,caloJet2Phi_);
	  if(l1CaloJet2DeltaR_tmp < l1CaloJet2DeltaR)
	    {
	      secondItrTmp = secondItr;
	      l1Jet2MatchedPt = l1JetPt_->at(jitr);
	      l1Jet2MatchedEta = l1JetEta_->at(jitr);
	      l1Jet2MatchedPhi = l1JetPhi_->at(jitr);
	      l1CaloJet2DeltaR = l1CaloJet2DeltaR_tmp;
	      secondItr = jitr;
	    }
	}
      //avoid matching both the l1Jets to the same caloJet
      if(firstItr == secondItr && firstItr != -1 && secondItrTmp != -1)
	{
	  l1Jet2MatchedPt = l1JetPt_->at(secondItrTmp);
	  l1Jet2MatchedEta = l1JetEta_->at(secondItrTmp);
	  l1Jet2MatchedPhi = l1JetPhi_->at(secondItrTmp);
	}
	
      float l1MatchedMjj = sqrt(2*l1Jet1MatchedPt*l1Jet2MatchedPt*(cosh(l1Jet1MatchedEta-l1Jet2MatchedEta) - cos(l1Jet1MatchedPhi-l1Jet2MatchedPhi)));


      
      //#######################################
      //do I want to select the matched analysis?
      float l1Mjj = l1Mjj_;
      float l1Jet1Pt = l1Jet1Pt_;
      float l1Jet2Pt = l1Jet2Pt_;
      float l1Jet1Eta = l1Jet1Eta_;
      float l1Jet2Eta = l1Jet2Eta_;
      float l1Jet1Phi = l1Jet1Phi_;
      float l1Jet2Phi = l1Jet2Phi_;
      float l1DeltaEta = l1DeltaEta_;
      if (matched)
	{
	  l1Mjj = l1MatchedMjj;
	  l1Jet1Pt = l1Jet1MatchedPt;
	  l1Jet2Pt = l1Jet2MatchedPt;
	  l1Jet1Eta = l1Jet1MatchedEta;
	  l1Jet2Eta = l1Jet2MatchedEta;
	  l1Jet1Phi = l1Jet1MatchedPhi;
	  l1Jet2Phi = l1Jet2MatchedPhi;
	  l1DeltaEta = fabs(l1Jet1MatchedEta-l1Jet2MatchedEta);
	}
      //#######################################
      //do I want wide jets?
      float caloJet1Pt = caloJet1Pt_;
      float caloJet1Eta = caloJet1Eta_;
      float caloJet1Phi = caloJet1Phi_;
      float caloJet2Pt = caloJet2Pt_;
      float caloJet2Eta = caloJet2Eta_;
      float caloJet2Phi = caloJet2Phi_;
      float caloMjj = caloMjj_;
      float caloDeltaEta = caloDeltaEta_;
      if (wideJets)
	{
	  caloJet1Pt = caloWJet1Pt_;
	  caloJet1Eta = caloWJet1Eta_;
	  caloJet1Phi = caloWJet1Phi_;
	  caloJet2Pt = caloWJet2Pt_;
	  caloJet2Eta = caloWJet2Eta_;
	  caloJet2Phi = caloWJet2Phi_;
	  caloMjj = caloWMjj_;
	  caloDeltaEta = caloWDeltaEta_;
	}
      //#######################################


      
      //resolution and comparison L1 VS Calo (upstream analysis cuts)
      if(l1Mjj > 0 && caloMjj > 0)
	{
	  mjj_res->Fill(2*(l1Mjj-caloMjj)/(l1Mjj+caloMjj));
	  mjj_diff->Fill(l1Mjj-caloMjj);
	  pt1_res->Fill(2*(l1Jet1Pt-caloJet1Pt)/(l1Jet1Pt+caloJet1Pt));
	  pt2_res->Fill(2*(l1Jet2Pt-caloJet2Pt)/(l1Jet2Pt+caloJet2Pt));

	  float deltaR1_ = sqrt(f_deltaRR(l1Jet1Eta,l1Jet1Phi,caloJet1Eta,caloJet1Phi));
	  deltaR1->Fill(deltaR1_);
	  float deltaR2_ = sqrt(f_deltaRR(l1Jet2Eta,l1Jet2Phi,caloJet2Eta,caloJet2Phi));
	  deltaR2->Fill(deltaR2_);
	  deltaPhi->Fill(f_deltaPhi(l1Jet1Phi,l1Jet2Phi));

	  deltaR1_vsMjj->Fill(caloMjj,deltaR1_);
	  deltaPhi_vsMjj->Fill(caloMjj,f_deltaPhi(l1Jet1Phi,l1Jet2Phi));
	  deltaMjj_vsMjj->Fill(caloMjj,l1Mjj-caloMjj);
	}



      bool l1Pass = (l1Mjj > 150. &&
		     
      		     l1Jet1Pt > 15. &&
      		     l1Jet2Pt > 15. &&
      		     fabs(l1Jet1Eta) < 5.0 &&
      		     fabs(l1Jet2Eta) < 5.0 &&
      		     l1DeltaEta < 2.0
      		     );
     
      
      
      //analysis cuts needed to compare to the analysis
      //calo analysis
      if (caloJet1Pt > 60. &&
	  caloJet2Pt > 30. &&
	  fabs(caloJet1Eta) < 2.5 &&
	  fabs(caloJet2Eta) < 2.5 &&
	  caloDeltaEta < 1.3
	  )
	{
	  caloMjjSpectrum->Fill(caloMjj);
	  caloHT250_eff->Fill((hltAccept_->at(HT250Calo)==1 && l1Accept_->at(L1scenario)==1), caloMjj);

	  HTT240_eff->Fill(l1Accept_->at(HTT240)==1, caloMjj);

	  l1Mjj_eff->Fill(l1Pass, caloMjj);
	  caloMjj_eff->Fill((caloMjj>200 && l1Pass), caloMjj);

	  if(!l1Pass && caloMjj > 300)
	    {
	      std::cout << std::endl;
	      std::cout << std::fixed << std::setprecision(2)
			<< "l1Mjj-caloMjj = " << l1Mjj-caloMjj << " l1Mjj        = " << l1Mjj             << " caloMjj    = " << caloMjj << std::endl;
	      std::cout << "  caloJet1Pt= " << caloJet1Pt   << " caloJet1Eta_= " << caloJet1Eta_     << " caloJet1Phi_= " << caloJet1Phi_ << std::endl;
	      std::cout << "  l1Jet1Pt    = " << l1Jet1Pt      << " l1Jet1Eta    = " << l1Jet1Eta        << " l1Jet1Phi    = " << l1Jet1Phi << std::endl;
	      std::cout << std::endl;
	      std::cout << "  caloJet2Pt= " << caloJet2Pt     << " caloJet2Eta= " << caloJet2Eta     << " caloJet2Phi= " << caloJet2Phi << std::endl;
	      std::cout << "  l1Jet2Pt    = " << l1Jet2Pt        << " l1Jet2Eta    = " << l1Jet2Eta        << " l1Jet2Phi    = " << l1Jet2Phi << std::endl;
	      std::cout << "==========================================================" << std::endl;
	    }

	  //l1 and hlt rates
	  for(unsigned int ii=0; ii<l1Names_->size(); ++ii)
	    if (l1Accept_->at(ii)==1)
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
  leg0->Draw("sames");

  TCanvas* c2 = new TCanvas();
  for(unsigned int ii=0; ii<l1Names_->size(); ++ii)
    l1->GetXaxis()->SetBinLabel(ii+1,l1Names_->at(ii).c_str());
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

  std::string folder;
  if(matched)
    folder = "matched";
  else
    folder = "notMatched";
  if(wideJets)
    folder += "WideJets/";
  else
    folder += "CaloJets/";
  c1->Print((folder+"turnOn.pdf").c_str(),"pdf");
  c2->Print((folder+"L1Rates.pdf").c_str(),"pdf");
  c3->Print((folder+"mjj_res.pdf").c_str(),"pdf");
  c31->Print((folder+"mjj_diff.pdf").c_str(),"pdf");
  c4->Print((folder+"deltaR1.pdf").c_str(),"pdf");
  c41->Print((folder+"deltaR2.pdf").c_str(),"pdf");
  c5->Print((folder+"pt1_res.pdf").c_str(),"pdf");
  c6->Print((folder+"pt2_res.pdf").c_str(),"pdf");
  c7->Print((folder+"deltaR1_vsMjj.pdf").c_str(),"pdf");
  c8->Print((folder+"deltaPhi_vsMjj.pdf").c_str(),"pdf");
  c9->Print((folder+"deltaMjj_vsMjj.pdf").c_str(),"pdf");
  
  return 0;


}
