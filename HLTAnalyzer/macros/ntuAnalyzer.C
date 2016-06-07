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
}


int ntuAnalyzer(std::string fileName)
{
  setGlobalStyle();
  
  //###############################
  unsigned int HT250Calo  = 0; //index of DST_HT250_CaloScouting_v
  float HT250Calo_rate = 1928;

  unsigned int HT410PF = 2; //index of DST_HT410_PFScouting_v
  float HT410PF_rate = 294;

  float instLumi = 0.4; //E34
  float targetLumi = 1; //E34
  float lumiScaleFactor = targetLumi/instLumi;

  float PDRate = 59300; //unprescaled rate of HLTPhysics accordingly to: https://cmswbm2.web.cern.ch/cmswbm2/cmsdb/servlet/DatasetSummary?RUN=274200 and prescale of 9000
  
  //###############################

  TChain* tt = new TChain("MyAnalysis/HLTree");
  tt->Add(fileName.c_str());

  //set branches
  TBranch* b_caloMjj;
  TBranch* b_PFMjj;
  TBranch* b_hltAccept;

  float caloMjj = 0;
  float PFMjj = 0;
  std::vector<int>* hltAccept = 0;

  tt->SetBranchAddress("caloMjj", &caloMjj, &b_caloMjj);
  tt->SetBranchAddress("PFMjj", &PFMjj, &b_PFMjj);
  tt->SetBranchAddress("hltAccept", &hltAccept, &b_hltAccept);

  int nentries = tt->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;

  //book graphs and plots
  float min = 100.;
  float max = 1000.;
  int nBins = 18;

  TF1* f1 = new TF1("f1","[0]*TMath::Erf((x-[1])/[2])-[0]*TMath::Erf((-x-[1])/[2])",min,max);
  f1->SetLineWidth(2.);

  TEfficiency* mjj450_eff = new TEfficiency("mjj450_eff","mjj450_eff",nBins,min,max);
  mjj450_eff->SetMarkerColor(kRed);
  mjj450_eff->SetLineColor(kRed);
  mjj450_eff->SetLineWidth(2.);
  mjj450_eff->SetTitle("turnOn;Mjj [GeV]");
  TEfficiency* mjj280_eff = new TEfficiency("mjj280_eff","mjj280_eff",nBins,min,max);
  mjj280_eff->SetLineWidth(2);
  TEfficiency* pf410_eff = new TEfficiency("pf410_eff","pf410_eff",nBins,min,max);
  pf410_eff->SetMarkerColor(kOrange+1);
  pf410_eff->SetLineColor(kOrange+1);
  TEfficiency* calo250_eff = new TEfficiency("calo250_eff","calo250_eff",nBins,min,max);
  calo250_eff->SetMarkerColor(kBlue);
  calo250_eff->SetLineColor(kBlue);

  //loop
  for (Long64_t jentry=0; jentry<nentries;++jentry)
    {
      tt->GetEntry(jentry);
      mjj450_eff->Fill(caloMjj>450, PFMjj);
      mjj280_eff->Fill(caloMjj>280, caloMjj); //this is illustrative
      pf410_eff->Fill(hltAccept->at(HT410PF)==1, PFMjj);
      calo250_eff->Fill(hltAccept->at(HT250Calo)==1, PFMjj);
    }

  f1->SetParameters(0.5,400,50);  
  mjj450_eff->Fit(f1,"r");

  // f1->SetParameters(0.5,280,10);
  // mjj280_eff->Fit(f1,"r");


  TLegend* leg0 = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg0->AddEntry(mjj450_eff,"MJJ450","L");
  leg0->AddEntry(mjj280_eff,"MJJ280","L");
  leg0->AddEntry(pf410_eff,"HT410_PF","P");
  leg0->AddEntry(calo250_eff,"HT250_Calo","P");

  TCanvas* c0 = new TCanvas();
  mjj450_eff->Draw();
  mjj280_eff->Draw("L,sames");
  pf410_eff->Draw("sames");
  calo250_eff->Draw("sames");
  leg0->Draw("sames");


  //return 0;

  //##############################################
  //##############################################

  //book graphs and plots
  TGraph* totRateVsCut = new TGraph();
  TGraph* pureRateVsCut450 = new TGraph();
  TGraph* pureRateVsCut280 = new TGraph();
  TGraph* totRateVsCut_proj = new TGraph();
  TGraph* pureRateVsCut_proj = new TGraph();

  //loops
  int bin = 0;
  for (int cut = 100; cut < 600; cut=cut+5)
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

	  if (caloMjj > cut && !hltAccept->at(HT410PF))
	    ++excl410_passed;
	  if (caloMjj > cut && !hltAccept->at(HT250Calo))
	    ++excl250_passed;
	  if (caloMjj > cut)
	    ++mjjPassed;

	  // if (hltAccept->at(HT250Calo) == 0 && mjj > cut)
	  //   std::cout << "ref trigger doesn't completely cover cut at " << cut << std::endl;
	}
      // float mjjTotalRate = (float)mjjPassed/(float)HT250Calo_Passed*HT250Calo_rate;
      // float mjjPureRate = (float)exclPassed/(float)HT250Calo_Passed*HT250Calo_rate;
      float mjjTotalRate = (float)mjjPassed/(float)nentries*PDRate;
      float mjj450_PureRate = (float)excl410_passed/(float)nentries*PDRate;
      float mjj280_PureRate = (float)excl250_passed/(float)nentries*PDRate;

      totRateVsCut->SetPoint(bin,cut,mjjTotalRate);
      pureRateVsCut450->SetPoint(bin,cut,mjj450_PureRate);
      pureRateVsCut280->SetPoint(bin,cut,mjj280_PureRate);

      // totRateVsCut_proj->SetPoint(bin,cut,mjjTotalRate/instLumi*targetLumi);
      // pureRateVsCut_proj->SetPoint(bin,cut,mjjPureRate/instLumi*targetLumi);

      ++bin;
    }

  //plotting and styling
  TLegend* leg = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg->AddEntry(totRateVsCut,"total rate","P");
  leg->AddEntry(pureRateVsCut450,"pure rate wrt HT410PF","P");
  leg->AddEntry(pureRateVsCut280,"pure rate wrt HT250Calo","P");

  totRateVsCut->SetTitle("Rate Ref");
  totRateVsCut_proj->SetTitle("Rate target");

  totRateVsCut->GetXaxis()->SetTitle("Mjj cut threshold [GeV]");
  totRateVsCut->GetYaxis()->SetTitle("Rate [Hz]");
  pureRateVsCut450->SetMarkerColor(kRed);
  pureRateVsCut280->SetMarkerColor(kOrange+1);

  totRateVsCut_proj->GetXaxis()->SetTitle("Mjj cut threshold [GeV]");
  totRateVsCut_proj->GetYaxis()->SetTitle("Rate [Hz]");
  pureRateVsCut_proj->SetMarkerColor(kRed);

  TCanvas* c2 = new TCanvas();
  c2->cd();
  totRateVsCut->Draw("AP");
  pureRateVsCut450->Draw("P,sames");
  pureRateVsCut280->Draw("P,sames");
  leg->Draw("sames");


  // TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),
  //  			    (totRateVsCut->GetYaxis()->GetBinLowEdge(1))*lumiScaleFactor,
  // 			    (totRateVsCut->GetYaxis()->GetBinLowEdge(totRateVsCut->GetYaxis()->GetNbins())+totRateVsCut->GetYaxis()->GetBinWidth(1))*lumiScaleFactor,510,"+L");
  //  axis->SetLineColor(kRed);
  //  axis->SetTextColor(kRed);
  //  axis->Draw("same");


  // TCanvas* c3 = new TCanvas();
  // totRateVsCut_proj->Draw("AP");
  // pureRateVsCut_proj->Draw("P,sames");
  // leg->Draw("sames");


  return 0;
}
