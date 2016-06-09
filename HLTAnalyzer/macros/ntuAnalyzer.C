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
      mjj450_eff->Fill(caloMjj>450 || hltAccept->at(HT410PF)==1, PFMjj);
      mjj280_eff->Fill(caloMjj>200 || hltAccept->at(HT250Calo)==1, caloMjj); //this is illustrative
      pf410_eff->Fill(hltAccept->at(HT410PF)==1, PFMjj);
      calo250_eff->Fill(hltAccept->at(HT250Calo)==1, caloMjj);
    }

  mjj450_eff->Fit(f1,"r");
  mjj280_eff->Fit(f2,"r");


  TLegend* leg0 = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg0->AddEntry(mjj450_eff,"MJJ450PF || HT410PF","L");
  leg0->AddEntry(mjj280_eff,"MJJ200Calo || HT250Calo","L");
  leg0->AddEntry(pf410_eff,"HT410_PF","P");
  leg0->AddEntry(calo250_eff,"HT250_Calo","P");

  TCanvas* c0 = new TCanvas();
  mjj450_eff->Draw();
  mjj280_eff->Draw("sames");
  pf410_eff->Draw("sames");
  calo250_eff->Draw("sames");
  leg0->Draw("sames");


  //return 0;

  //##############################################
  //##############################################

  //book graphs and plots
  TGraphErrors* totRateVsCut = new TGraphErrors();
  TGraphErrors* pureRateVsCut450 = new TGraphErrors();
  TGraphErrors* pureRateVsCut280 = new TGraphErrors();
  TGraphErrors* totRateVsCut_proj = new TGraphErrors();
  TGraphErrors* pureRateVsCut_proj = new TGraphErrors();

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
  totRateVsCut->GetYaxis()->SetTitle("Rate @4E33 [Hz]");
  pureRateVsCut450->SetMarkerColor(kRed);
  pureRateVsCut450->SetLineColor(kRed);
  pureRateVsCut280->SetMarkerColor(kOrange+1);
  pureRateVsCut280->SetLineColor(kOrange+1);

  totRateVsCut_proj->GetXaxis()->SetTitle("Mjj cut threshold [GeV]");
  totRateVsCut_proj->GetYaxis()->SetTitle("Rate [Hz]");
  pureRateVsCut_proj->SetMarkerColor(kRed);
  pureRateVsCut_proj->SetLineColor(kRed);

  TCanvas* c2 = new TCanvas();
  c2->cd();
  totRateVsCut->Draw("AP");
  pureRateVsCut450->Draw("P,sames");
  pureRateVsCut280->Draw("P,sames");
  leg->Draw("sames");
  c2->Update();

  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),
    			    (totRateVsCut->GetYaxis()->GetBinLowEdge(1))*lumiScaleFactor,
			    (totRateVsCut->GetYaxis()->GetBinLowEdge(totRateVsCut->GetYaxis()->GetNbins())+totRateVsCut->GetYaxis()->GetBinWidth(1))*lumiScaleFactor,510,"+L");

  c2->SetTicky(0);
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetTextColor(kRed);
  axis->SetTitleOffset(1.3);
  axis->SetLabelSize(0.03);
  axis->SetTitle("Rate @1E34 [Hz]");
  axis->Draw();


  // TCanvas* c3 = new TCanvas();
  // totRateVsCut_proj->Draw("AP");
  // pureRateVsCut_proj->Draw("P,sames");
  // leg->Draw("sames");


  return 0;
}
