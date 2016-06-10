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


int ntuAnalyzer(std::string fileName)
{
  setGlobalStyle();
  
  //###############################
  //## run274200 ##
  unsigned int HT250Calo  = 1; //index of DST_HT250_CaloScouting_v
  float HT250Calo_rate = 1928;

  unsigned int HT410PF = 3; //index of DST_HT410_PFScouting_v
  float HT410PF_rate = 294;

  unsigned int MJJ200Calo = 10; //index of DST_DiCaloWideJetMass200_CaloScouting_v
  
  unsigned int HTT200 = 0; //index if L1_HTT200
  unsigned int HTT240 = 1; //index if L1_HTT240
  unsigned int HTT270 = 2; //index if L1_HTT270
  unsigned int HTT280 = 3; //index if L1_HTT280
  unsigned int DoubleJetC100 = 7; //index if L1_DoubleJetC100
  unsigned int DoubleJetC112 = 8; //index if L1_DoubleJetC112

  float instLumi = 0.4; //E34
  float targetLumi = 1; //E34
  float lumiScaleFactor = targetLumi/instLumi;

  float PDRate = 59300; //unprescaled rate of HLTPhysics accordingly to: https://cmswbm2.web.cern.ch/cmswbm2/cmsdb/servlet/DatasetSummary?RUN=274200 and prescale of 9000
  

  unsigned int L1scenario = HTT240;
  //###############################

  TChain* tt = new TChain("MyAnalysis/HLTree");
  tt->Add(fileName.c_str());

  //set branches
  TBranch* b_caloMjj;
  TBranch* b_PFMjj;
  TBranch* b_hltAccept;
  TBranch* b_l1Accept;
  TBranch* b_l1Names;

  float caloMjj = 0;
  float PFMjj = 0;
  std::vector<int>* hltAccept = 0;
  std::vector<int>* l1Accept = 0;
  std::vector<string>* l1Names = 0;

  tt->SetBranchAddress("caloMjj", &caloMjj, &b_caloMjj);
  tt->SetBranchAddress("PFMjj", &PFMjj, &b_PFMjj);
  tt->SetBranchAddress("hltAccept", &hltAccept, &b_hltAccept);
  tt->SetBranchAddress("l1Accept", &l1Accept, &b_l1Accept);
  tt->SetBranchAddress("l1Names", &l1Names, &b_l1Names);

  int nentries = tt->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;

  //book graphs and plots
  float min = 150.;
  float max = 1000.;
  int nBins = 17;

  TF1* f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])-[0]*TMath::Erf((-x-[1])/[2])",min,max);
  f1->SetParameters(0.5,350,10);  
  f1->FixParameter(0,0.5);
  f1->SetLineWidth(2.);
  f1->SetLineColor(kRed);

  TF1* f2 = (TF1*)f1->Clone("f2");
  f2->SetParameters(0.5,150,10);
  f2->SetLineColor(kBlack);

  TEfficiency* mjj450_eff = new TEfficiency("mjj450_eff","mjj450_eff",nBins,min,max);
  mjj450_eff->SetMarkerColor(kRed);
  mjj450_eff->SetLineColor(kRed);
  mjj450_eff->SetLineWidth(2);
  mjj450_eff->SetTitle("turnOn;Mjj [GeV]");
  TEfficiency* mjj200_eff = new TEfficiency("mjj200_eff","mjj200_eff",nBins,min,max);
  mjj200_eff->SetLineWidth(2);
  mjj200_eff->SetTitle("turnOn;Mjj [GeV]");
  TEfficiency* pf410_eff = new TEfficiency("pf410_eff","pf410_eff",nBins,min,max);
  pf410_eff->SetMarkerColor(kOrange+1);
  pf410_eff->SetLineColor(kOrange+1);
  TEfficiency* calo250_eff = new TEfficiency("calo250_eff","calo250_eff",nBins,min,max);
  calo250_eff->SetMarkerColor(kBlue);
  calo250_eff->SetLineColor(kBlue);

  TH1F* l1 = new TH1F("l1","l1",9,0.,9.);
  
  //loop
  for (Long64_t jentry=0; jentry<nentries;++jentry)
    {
      tt->GetEntry(jentry);

      mjj450_eff->Fill((caloMjj>450 && l1Accept->at(L1scenario)) || hltAccept->at(HT410PF)==1, PFMjj);
      mjj200_eff->Fill((caloMjj>200 && l1Accept->at(L1scenario)) || hltAccept->at(HT250Calo)==1, caloMjj);
      //for comparison
      pf410_eff->Fill(hltAccept->at(HT410PF)==1, PFMjj);
      calo250_eff->Fill(hltAccept->at(HT250Calo)==1, caloMjj);

      //l1 and hlt rates
      for(unsigned int ii=0; ii<l1Names->size(); ++ii)
	if (l1Accept->at(ii)==1)
	  l1->Fill(ii);
    }

  mjj450_eff->Fit(f1,"r");
  mjj200_eff->Fit(f2,"r");


  TLegend* leg0 = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg0->AddEntry(mjj450_eff,"MJJ450PF || HT410PF","L");
  leg0->AddEntry(mjj200_eff,"MJJ200Calo || HT250Calo","L");
  leg0->AddEntry(pf410_eff,"HT410_PF","P");
  leg0->AddEntry(calo250_eff,"HT250_Calo","P");

  TCanvas* c1 = new TCanvas();
  mjj450_eff->Draw();
  mjj200_eff->Draw("sames");
  pf410_eff->Draw("sames");
  calo250_eff->Draw("sames");
  leg0->Draw("sames");

  TCanvas* c2 = new TCanvas();
  l1->Scale(PDRate/nentries);

  for(unsigned int ii=0; ii<l1Names->size(); ++ii)
    l1->GetXaxis()->SetBinLabel(ii+1,l1Names->at(ii).c_str());
  l1->GetYaxis()->SetTitle("L1 Rate @4E33 [Hz]");
  l1->SetMaximum(l1->GetMaximum()+200);
  
  l1->Draw();
  c2->Update();

  TGaxis *l1axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),
			      l1->GetMinimum()*lumiScaleFactor,
			      l1->GetMaximum()*lumiScaleFactor,510,"+L");

  c2->SetTicky(0);
  l1axis->SetLineColor(kRed);
  l1axis->SetLabelColor(kRed);
  l1axis->SetTextColor(kRed);
  l1axis->SetTitleOffset(1.3);
  l1axis->SetLabelSize(0.03);
  l1axis->SetTitle("L1 Rate @1E34 [Hz]");
  l1axis->Draw();
  
  
  //return 0;

  //##############################################
  //##############################################

  //book graphs and plots
  TGraphErrors* totRateVsCut = new TGraphErrors();
  TGraphErrors* pureRateVsCut450 = new TGraphErrors();
  TGraphErrors* pureRateVsCut280 = new TGraphErrors();

  //loops
  int bin = 0;
  for (int cut = 100; cut < 550; cut=cut+10)
    {
      int mjjPassed = 0;
      int HT250Calo_Passed = 0;
      int excl410_passed = 0;
      int excl250_passed = 0;
      for (Long64_t jentry=0; jentry<nentries;++jentry) 
	{
	  tt->GetEntry(jentry);

	  if (hltAccept->at(HT250Calo) == 1)
	    ++HT250Calo_Passed;

	  //if (caloMjj > cut && !hltAccept->at(HT410PF))
	  if (caloMjj > cut && l1Accept->at(L1scenario) && !hltAccept->at(HT410PF))
	    ++excl410_passed;
	  if (caloMjj > cut && l1Accept->at(L1scenario) && !hltAccept->at(HT250Calo))
	    ++excl250_passed;
	  if (caloMjj > cut && l1Accept->at(L1scenario))
	    ++mjjPassed;

	  // if (hltAccept->at(HT250Calo) == 0 && mjj > cut)
	  //   std::cout << "ref trigger doesn't completely cover cut at " << cut << std::endl;
	}
      // float mjjTotalRate = (float)mjjPassed/(float)HT250Calo_Passed*HT250Calo_rate;
      // float mjjPureRate = (float)exclPassed/(float)HT250Calo_Passed*HT250Calo_rate;

      float sigmaMjjPassed = sqrt((float)mjjPassed);
      float sigmaNentries = sqrt((float)nentries);
      float sigmaExcl410_passed = sqrt((float)excl410_passed);
      float sigmaExcl250_passed = sqrt((float)excl250_passed);

      float mjjTotalRate = (float)mjjPassed/(float)nentries*PDRate;
      float mjjTotalRateE = PDRate*sqrt(pow((sigmaMjjPassed/nentries),2)+pow((sigmaNentries*mjjPassed/nentries/nentries),2));

      float mjj450_PureRate = (float)excl410_passed/(float)nentries*PDRate;
      float mjj450_PureRateE = PDRate*sqrt(pow((sigmaExcl410_passed/nentries),2)+pow((sigmaNentries*excl410_passed/nentries/nentries),2));
      
      float mjj280_PureRate = (float)excl250_passed/(float)nentries*PDRate;
      float mjj280_PureRateE = PDRate*sqrt(pow((sigmaExcl250_passed/nentries),2)+pow((sigmaNentries*excl250_passed/nentries/nentries),2));

      totRateVsCut->SetPoint(bin,cut,mjjTotalRate);
      totRateVsCut->SetPointError(bin,0.,mjjTotalRateE);

      pureRateVsCut450->SetPoint(bin,cut,mjj450_PureRate);
      pureRateVsCut450->SetPointError(bin,0.,mjj450_PureRateE);

      pureRateVsCut280->SetPoint(bin,cut,mjj280_PureRate);
      pureRateVsCut280->SetPointError(bin,0.,mjj280_PureRateE);

      ++bin;
    }

  //plotting and styling
  TLegend* leg = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg->AddEntry(totRateVsCut,"total rate","P");
  leg->AddEntry(pureRateVsCut450,"pure rate wrt HT410PF","P");
  leg->AddEntry(pureRateVsCut280,"pure rate wrt HT250Calo","P");

  totRateVsCut->SetTitle("Rate Ref");

  totRateVsCut->GetXaxis()->SetTitle("Mjj cut threshold [GeV]");
  totRateVsCut->GetYaxis()->SetTitle("Rate @4E33 [Hz]");
  pureRateVsCut450->SetMarkerColor(kRed);
  pureRateVsCut450->SetLineColor(kRed);
  pureRateVsCut280->SetMarkerColor(kOrange+1);
  pureRateVsCut280->SetLineColor(kOrange+1);

  TCanvas* c3 = new TCanvas();
  c3->cd();
  totRateVsCut->Draw("AP");
  pureRateVsCut450->Draw("P,sames");
  pureRateVsCut280->Draw("P,sames");
  leg->Draw("sames");
  c3->Update();

  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),
    			    (totRateVsCut->GetYaxis()->GetBinLowEdge(1))*lumiScaleFactor,
			    (totRateVsCut->GetYaxis()->GetBinLowEdge(totRateVsCut->GetYaxis()->GetNbins())+totRateVsCut->GetYaxis()->GetBinWidth(1))*lumiScaleFactor,510,"+L");

  c3->SetTicky(0);
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetTextColor(kRed);
  axis->SetTitleOffset(1.3);
  axis->SetLabelSize(0.03);
  axis->SetTitle("Rate @1E34 [Hz]");
  axis->Draw();


  return 0;
}
