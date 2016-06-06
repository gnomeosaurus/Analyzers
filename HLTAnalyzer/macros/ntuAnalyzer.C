//calculate rate vs cut for scouting triggers
void setGlobalStyle()
{
  // For the statistics box:                                                                                                                                                                      
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(000000); // To display the mean and RMS:   SetOptStat("mr");                                                                                                               
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
  unsigned int lowerRef  = 0; //DST_HT250_CaloScouting_v
  float lowerRefRate = 1928;

  unsigned int higherRef = 1; //DST_HT410_PFScouting_v
  float higherRefRate = 294;

  float instLumi = 0.4; //E34
  float targetLumi = 1; //E34
  
  //###############################

  TChain* tt = new TChain("MyAnalysis/HLTree");
  tt->Add(fileName.c_str());

  //set branches
  TBranch* b_mjj;
  TBranch* b_hltAccept;

  float mjj = 0;
  std::vector<int>* hltAccept = 0;

  tt->SetBranchAddress("mjj", &mjj, &b_mjj);
  tt->SetBranchAddress("hltAccept", &hltAccept, &b_hltAccept);

  int nentries = tt->GetEntries();
  std::cout << "Number of entries: " << nentries << std::endl;


  //book of graphs
  TGraph* totRateVsCut = new TGraph();
  TGraph* pureRateVsCut = new TGraph();
  TGraph* totRateVsCut_proj = new TGraph();
  TGraph* pureRateVsCut_proj = new TGraph();

  //loop
  int bin = 0;
  for (int cut = 300; cut < 600; cut=cut+5)
    {
      int mjjPassed = 0;
      int lowerRefPassed = 0;
      int higherRefPassed = 0;
      int exclPassed = 0;
      for (Long64_t jentry=0; jentry<nentries;++jentry) 
	{
	  tt->GetEntry(jentry);

	  if (hltAccept->at(lowerRef) == 1)
	    ++lowerRefPassed;
	  if (hltAccept->at(higherRef) == 1)
	    ++higherRefPassed;

	  if (mjj > cut && hltAccept->at(higherRef) == 0)
	    ++exclPassed;
	  if (mjj > cut)
	    ++mjjPassed;

	  // if (hltAccept->at(lowerRef) == 0 && mjj > cut)
	  //   std::cout << "ref trigger doesn't completely cover cut at " << cut << std::endl;
	}
      float mjjTotalRate = (float)mjjPassed/(float)lowerRefPassed*lowerRefRate;
      float mjjPureRate = (float)exclPassed/(float)higherRefPassed*higherRefRate;

      totRateVsCut->SetPoint(bin,cut,mjjTotalRate);
      pureRateVsCut->SetPoint(bin,cut,mjjPureRate);
      totRateVsCut_proj->SetPoint(bin,cut,mjjTotalRate/instLumi*targetLumi);
      pureRateVsCut_proj->SetPoint(bin,cut,mjjPureRate/instLumi*targetLumi);
      ++bin;
    }

  //plotting and styling
  TLegend* leg = new TLegend(0.62, 0.78, 0.83, 0.89);
  leg->AddEntry(totRateVsCut,"total rate","P");
  leg->AddEntry(pureRateVsCut,"pure rate","P");

  totRateVsCut->SetTitle("Rate Ref");
  totRateVsCut_proj->SetTitle("Rate target");

  totRateVsCut->GetXaxis()->SetTitle("Mjj cut threshold [GeV]");
  totRateVsCut->GetYaxis()->SetTitle("Rate [Hz]");
  pureRateVsCut->SetMarkerColor(kRed);
  pureRateVsCut_proj->SetMarkerColor(kRed);

  TCanvas* c1 = new TCanvas();
  totRateVsCut->Draw("AP");
  pureRateVsCut->Draw("P,sames");
  leg->Draw("sames");

  TCanvas* c2 = new TCanvas();
  totRateVsCut_proj->Draw("AP");
  pureRateVsCut_proj->Draw("P,sames");
  leg->Draw("sames");

  return 0;
}
