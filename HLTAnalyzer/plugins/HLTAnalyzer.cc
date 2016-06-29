#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1Trigger/interface/Jet.h"


#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"



class MyHLTAnalyzer : public edm::EDAnalyzer {
  
public:
  MyHLTAnalyzer(const edm::ParameterSet& cfg);
  virtual ~MyHLTAnalyzer() {};
  
  virtual void analyze (const edm::Event& event, const edm::EventSetup & eventSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  virtual void endRun  (const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  
private:
  
  void fillHltResults(const edm::Handle<edm::TriggerResults> &, 
		      const edm::TriggerNames &);

  void fillL1Jets(const edm::Handle<l1t::JetBxCollection> &,
		  const edm::Event   &);

  template<typename jetCollection>
  void fillMjj(const edm::Handle<jetCollection> &, float &, float &, float &, float &, float &, float &, float &, float &);

  void beginEvent();


  
  /// file service and tree
  edm::Service<TFileService> outfile_;

  TTree* outTree_;
  int run_;
  int evt_;
  int lumi_;
  //  int nVtx_;
  float caloMjj_;
  float PFMjj_;
  float l1Mjj_;

  float caloDeltaEta_;
  float PFDeltaEta_;
  float l1DeltaEta_;

  float caloJet1Pt_;
  float caloJet1Eta_;
  float caloJet1Phi_;
  float caloJet2Pt_;
  float caloJet2Eta_;
  float caloJet2Phi_;

  float PFJet1Pt_;
  float PFJet1Eta_;
  float PFJet1Phi_;
  float PFJet2Pt_;
  float PFJet2Eta_;
  float PFJet2Phi_;

  float l1Jet1Pt_;
  float l1Jet1Eta_;
  float l1Jet1Phi_;
  float l1Jet2Pt_;
  float l1Jet2Eta_;
  float l1Jet2Phi_;

  std::vector<float>* l1JetPt_;
  std::vector<float>* l1JetEta_;
  std::vector<float>* l1JetPhi_;

  std::vector<int>* hltAccept_;
  std::vector<int>* hltWasrun_;
  std::vector<std::string>* hltNames_;

  std::vector<std::string>* l1Names_;
  std::vector<int>* l1Accept_;

  // Trigger process
  edm::InputTag triggerResultTag_1_;
  edm::InputTag triggerResultTag_2_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_1_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_2_;
  edm::InputTag triggerEventTag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> triggerEventToken_;
  edm::EDGetTokenT<l1t::JetBxCollection> l1CandToken_;
  edm::EDGetTokenT<reco::CaloJetCollection>         caloJetToken_;
  edm::EDGetTokenT<reco::PFJetCollection>           PFJetToken_;

  l1t::L1TGlobalUtil *l1GtUtils_;

  double minMass_;
  double fatJetDeltaR_;
  double maxDeltaEta_;
  double maxJetEta_;
  double minJetPt_;
  double maxL1DeltaEta_;
  double maxL1JetEta_;
  double minL1JetPt_;

  std::vector<std::string> hltPaths_;
  edm::EDGetToken algToken_;
  std::vector<std::string> l1Paths_;
};

MyHLTAnalyzer::MyHLTAnalyzer(const edm::ParameterSet& cfg): 
  triggerResultTag_1_       (cfg.getUntrackedParameter<edm::InputTag>("triggerResult_1")),
  triggerResultTag_2_       (cfg.getUntrackedParameter<edm::InputTag>("triggerResult_2")),
  triggerResultToken_1_     (consumes<edm::TriggerResults>(triggerResultTag_1_)),
  triggerResultToken_2_     (consumes<edm::TriggerResults>(triggerResultTag_2_)),
  triggerEventTag_          (cfg.getUntrackedParameter<edm::InputTag>("triggerSummary")),
  triggerEventToken_        (consumes<trigger::TriggerFilterObjectWithRefs>(triggerEventTag_)),
  l1CandToken_              (consumes<l1t::JetBxCollection>(cfg.getUntrackedParameter<edm::InputTag>("l1CandTag"))),
  caloJetToken_             (consumes<reco::CaloJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("caloJetTag"))),
  PFJetToken_               (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJetTag"))),
  //params for wide jet calculation
  minMass_                  (cfg.getUntrackedParameter<double>("minMass")),
  fatJetDeltaR_             (cfg.getUntrackedParameter<double>("fatJetDeltaR")),
  maxDeltaEta_              (cfg.getUntrackedParameter<double>("maxDeltaEta")),
  maxJetEta_                (cfg.getUntrackedParameter<double>("maxJetEta")),
  minJetPt_                 (cfg.getUntrackedParameter<double>("minJetPt")),
  maxL1DeltaEta_            (cfg.getUntrackedParameter<double>("maxL1DeltaEta")),
  maxL1JetEta_              (cfg.getUntrackedParameter<double>("maxL1JetEta")),
  minL1JetPt_               (cfg.getUntrackedParameter<double>("minL1JetPt")),
  
  hltPaths_                 (cfg.getUntrackedParameter<std::vector<std::string> >("hltPaths")),

  algToken_                 (consumes<BXVector<GlobalAlgBlk>>(cfg.getParameter<edm::InputTag>("AlgInputTag"))),
  l1Paths_                  (cfg.getUntrackedParameter<std::vector<std::string> >("l1Paths"))
{
  l1GtUtils_ = new l1t::L1TGlobalUtil(cfg,consumesCollector());
}

void MyHLTAnalyzer::beginEvent()
{
  evt_  = -1;
  lumi_ = -1;
  run_  = -1;

  caloMjj_  = -1;
  PFMjj_  = -1;
  l1Mjj_  = -1;

  caloDeltaEta_ = -999;
  PFDeltaEta_ = -999;
  l1DeltaEta_ = -999;

  caloJet1Pt_ = -1;
  caloJet1Eta_ = -999;
  caloJet1Phi_ = -999;
  caloJet2Pt_ = -1;
  caloJet2Eta_ = -999;
  caloJet2Phi_ = -999;

  PFJet1Pt_ = -1;
  PFJet1Eta_ = -999;
  PFJet1Phi_ = -999;
  PFJet2Pt_ = -1;
  PFJet2Eta_ = -999;
  PFJet2Phi_ = -999;

  l1Jet1Pt_ = -1;
  l1Jet1Eta_ = -999;
  l1Jet1Phi_ = -999;
  l1Jet2Pt_ = -1;
  l1Jet2Eta_ = -999;
  l1Jet2Phi_ = -999;

  l1JetPt_->clear();
  l1JetEta_->clear();
  l1JetPhi_->clear();

  hltNames_->clear();
  hltAccept_->clear();
  hltWasrun_->clear();
  l1Names_->clear();
  l1Accept_->clear();
}
  
void MyHLTAnalyzer::fillHltResults(const edm::Handle<edm::TriggerResults>   & triggerResults, 
				   const edm::TriggerNames                  & triggerNames  )
{

  for (unsigned int itrig=0; itrig < triggerNames.size(); ++itrig) 
    for (unsigned int name=0; name < hltPaths_.size(); ++name)
      if(triggerNames.triggerName(itrig).find(hltPaths_[name]) != std::string::npos)
	{
	  hltNames_->push_back(triggerNames.triggerName(itrig));
	  hltAccept_->push_back(triggerResults->accept(itrig));
	  hltWasrun_->push_back(triggerResults->wasrun(itrig));
	  break;
	}
}

void MyHLTAnalyzer::fillL1Jets(const edm::Handle<l1t::JetBxCollection> & l1cands,
			       const edm::Event                        & event  )
{
  // Selected jets
  std::vector<l1t::Jet> l1jets;
  //select only BX=0
  auto i = l1cands->begin(0);
  for(;i != l1cands->end(0); i++){
    if(std::abs(i->eta()) < maxL1JetEta_ && i->pt() >= minL1JetPt_){
      l1t::Jet jet(i->p4());
      l1jets.push_back(jet);

      l1JetPt_->push_back(i->pt());
      l1JetEta_->push_back(i->eta());
      l1JetPhi_->push_back(i->phi());
      //std::cout << "  et:  "  << jet.et() << "  eta:  "  << jet.eta() << "  phi:  "  << jet.phi() << "\n";
    }
  }

  // events with at least two jets
  if(l1jets.size() < 2) return;

  math::PtEtaPhiMLorentzVector j1(0.1, 0., 0., 0.);
  math::PtEtaPhiMLorentzVector j2(0.1, 0., 0., 0.);
  double jetPt1 = 0.;
  double jetPt2 = 0.;

  // look for the two highest-pT jet
  auto l1jet ( l1jets.begin() );
  for (; l1jet != l1jets.end() ; l1jet++) {
    if(l1jet->pt() > jetPt1) {
      // downgrade the 1st jet to 2nd jet
      j2 = j1;
      jetPt2 = j2.pt();
      // promote this jet to 1st jet
      j1 = l1jet->p4();
      jetPt1 = l1jet->pt();
    } else if(l1jet->pt() > jetPt2) {
      // promote this jet to 2nd jet
      j2 = l1jet->p4();
      jetPt2 = l1jet->pt();
    }
  }

  // apply DeltaEta cut
  double DeltaEta = std::abs(j1.eta() - j2.eta());
  if(DeltaEta > maxL1DeltaEta_) return;

  //FILL HERE THE JET VARIABLES IN THE TUPLE
  l1DeltaEta_ = DeltaEta;
  l1Jet1Pt_ = j1.pt();
  l1Jet1Eta_ = j1.eta();
  l1Jet1Phi_ = j1.phi();
  l1Jet2Pt_ = j2.pt();
  l1Jet2Eta_ = j2.eta();
  l1Jet2Phi_ = j2.phi();

  //as mich as possible similar to L1 math
  l1Mjj_ = sqrt(2*j1.pt()*j2.pt()*(cosh(j1.eta()-j2.eta()) - cos(j1.phi()-j2.phi())));
  return;
}



template<typename jetCollection>
void MyHLTAnalyzer::fillMjj(const edm::Handle<jetCollection> & jets, float & mjj, float & jet1Pt, float & jet2Pt, float & jet1Eta, float & jet2Eta, float & jet1Phi, float & jet2Phi, float & deltaEta)
{
  // Selected jets
  reco::CaloJetCollection recojets;
  typename jetCollection::const_iterator i = jets->begin();
  for(;i != jets->end(); i++){
    if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_){ 
      reco::CaloJet jet(i->p4(), i->vertex(), reco::CaloJet::Specific()); 
      recojets.push_back(jet);
    }
  }
  
  // events with at least two jets
  if(recojets.size() < 2) return;

  math::PtEtaPhiMLorentzVector j1(0.1, 0., 0., 0.);
  math::PtEtaPhiMLorentzVector j2(0.1, 0., 0., 0.);
  double jetPt1 = 0.;
  double jetPt2 = 0.;
  
  // look for the two highest-pT jet
  reco::CaloJetCollection::const_iterator recojet ( recojets.begin() );
  for (; recojet != recojets.end() ; recojet++) {
    if(recojet->pt() > jetPt1) {
      // downgrade the 1st jet to 2nd jet
      j2 = j1;
      jetPt2 = j2.pt();
      // promote this jet to 1st jet
      j1 = recojet->p4();
      jetPt1 = recojet->pt();
    } else if(recojet->pt() > jetPt2) {
      // promote this jet to 2nd jet
      j2 = recojet->p4();
      jetPt2 = recojet->pt();
    }
  }
  
  // apply DeltaEta cut
  double DeltaEta = std::abs(j1.eta() - j2.eta());
  if(DeltaEta > maxDeltaEta_) return;
  
  math::PtEtaPhiMLorentzVector fj1;
  math::PtEtaPhiMLorentzVector fj2;
  
  // apply radiation recovery (WideJets)
  for ( recojet = recojets.begin() ; recojet != recojets.end() ; recojet++) {
    double DeltaR1 = sqrt(pow(recojet->phi()-j1.phi(), 2.)+pow(recojet->eta()-j1.eta(),2.));
    double DeltaR2 = sqrt(pow(recojet->phi()-j2.phi(), 2.)+pow(recojet->eta()-j2.eta(),2.));
    if(DeltaR1 < DeltaR2 && DeltaR1 < fatJetDeltaR_) {
      fj1 += recojet->p4();
    } else if(DeltaR2 < fatJetDeltaR_) {
      fj2 += recojet->p4();
    }
  }

  fj1 += fj2;

  //FILL HERE THE JET VARIABLES IN THE TUPLE
  deltaEta = DeltaEta;
  jet1Pt = j1.pt();
  jet1Eta = j1.eta();
  jet1Phi = j1.phi();
  jet2Pt = j2.pt();
  jet2Eta = j2.eta();
  jet2Phi = j2.phi();


  mjj = fj1.mass();
  return;
}

//===================== beginJob and Analyze =========================
void MyHLTAnalyzer::beginJob() {

  TH1::SetDefaultSumw2() ;
  outTree_ = outfile_-> make<TTree>("HLTree","HLTree");

  outTree_->Branch("run",     &run_,       "run_/I");
  outTree_->Branch("evt",     &evt_,       "evt_/I");
  outTree_->Branch("lumi",    &lumi_,      "lumi_/I");
  //  outTree_->Branch("nvtx",    &nVtx_,      "nVtx_/I");

  outTree_->Branch("caloMjj",   &caloMjj_,     "caloMjj_/F");
  outTree_->Branch("PFMjj",     &PFMjj_,       "PFMjj_/F");
  outTree_->Branch("l1Mjj",     &l1Mjj_,       "l1Mjj_/F");

  outTree_->Branch("caloDeltaEta",  &caloDeltaEta_,    "caloDeltaEta_/F");
  outTree_->Branch("PFDeltaEta",    &PFDeltaEta_,      "PFDeltaEta_/F");
  outTree_->Branch("l1DeltaEta",    &l1DeltaEta_,      "l1DeltaEta_/F");

  outTree_->Branch("caloJet1Pt",   &caloJet1Pt_,     "caloJet1Pt_/F");
  outTree_->Branch("caloJet1Eta",  &caloJet1Eta_,    "caloJet1Eta_/F");
  outTree_->Branch("caloJet1Phi",  &caloJet1Phi_,    "caloJet1Phi_/F");
  outTree_->Branch("caloJet2Pt",   &caloJet2Pt_,     "caloJet2Pt_/F");
  outTree_->Branch("caloJet2Eta",  &caloJet2Eta_,    "caloJet2Eta_/F");
  outTree_->Branch("caloJet2Phi",  &caloJet2Phi_,    "caloJet2Phi_/F");

  outTree_->Branch("PFJet1Pt",   &PFJet1Pt_,     "PFJet1Pt_/F");
  outTree_->Branch("PFJet1Eta",  &PFJet1Eta_,    "PFJet1Eta_/F");
  outTree_->Branch("PFJet1Phi",  &PFJet1Phi_,    "PFJet1Phi_/F");
  outTree_->Branch("PFJet2Pt",   &PFJet2Pt_,     "PFJet2Pt_/F");
  outTree_->Branch("PFJet2Eta",  &PFJet2Eta_,    "PFJet2Eta_/F");
  outTree_->Branch("PFJet2Phi",  &PFJet2Phi_,    "PFJet2Phi_/F");

  outTree_->Branch("l1Jet1Pt",   &l1Jet1Pt_,     "l1Jet1Pt_/F");
  outTree_->Branch("l1Jet1Eta",  &l1Jet1Eta_,    "l1Jet1Eta_/F");
  outTree_->Branch("l1Jet1Phi",  &l1Jet1Phi_,    "l1Jet1Phi_/F");
  outTree_->Branch("l1Jet2Pt",   &l1Jet2Pt_,     "l1Jet2Pt_/F");
  outTree_->Branch("l1Jet2Eta",  &l1Jet2Eta_,    "l1Jet2Eta_/F");
  outTree_->Branch("l1Jet2Phi",  &l1Jet2Phi_,    "l1Jet2Phi_/F");

  l1JetPt_ = new std::vector<float>;
  l1JetEta_ = new std::vector<float>;
  l1JetPhi_ = new std::vector<float>;

  outTree_->Branch("l1JetPt",   "std::vector<float>",      &l1JetPt_);
  outTree_->Branch("l1JetEta",   "std::vector<float>",     &l1JetEta_);
  outTree_->Branch("l1JetPhi",   "std::vector<float>",     &l1JetPhi_);


  hltNames_ = new std::vector<std::string>;
  hltAccept_ = new std::vector<int>;
  hltWasrun_ = new std::vector<int>;

  outTree_->Branch("hltNames",  "std::vector<std::string>", &hltNames_);
  outTree_->Branch("hltAccept", "std::vector<int>", &hltAccept_);
  outTree_->Branch("hltWasrun", "std::vector<int>", &hltWasrun_);

  l1Names_ = new std::vector<std::string>;
  l1Accept_ = new std::vector<int>;

  outTree_->Branch("l1Names",  "std::vector<std::string>", &l1Names_);
  outTree_->Branch("l1Accept", "std::vector<int>", &l1Accept_);

}

void MyHLTAnalyzer::endJob() 
{
  delete hltNames_;
  delete hltAccept_;
  delete hltWasrun_;
  delete l1Names_;
  delete l1Accept_;
}


void MyHLTAnalyzer::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) 
{
  beginEvent();
  
  evt_             = event.id().event();
  lumi_            = event.id().luminosityBlock();
  run_             = event.id().run();
  
  //fill HLT results
  edm::Handle<edm::TriggerResults>   triggerResults_1;
  edm::Handle<edm::TriggerResults>   triggerResults_2;
  edm::Handle<trigger::TriggerFilterObjectWithRefs> triggerEvent;

  if (!event.getByToken(triggerResultToken_1_, triggerResults_1) ||
      !event.getByToken(triggerResultToken_2_, triggerResults_2))
    {
      edm::LogError("") << "Trigger collection not found";
      return;
    }
  
  edm::TriggerNames triggerNames_1 = event.triggerNames(*triggerResults_1);
  edm::TriggerNames triggerNames_2 = event.triggerNames(*triggerResults_2);
  
  fillHltResults(triggerResults_1,triggerNames_1);  
  fillHltResults(triggerResults_2,triggerNames_2);  

  //fill L1 results
  l1GtUtils_->retrieveL1(event,eventSetup,algToken_);
  for( unsigned int iseed = 0; iseed < l1Paths_.size(); iseed++ ) {
    bool l1htbit = 0;
    l1GtUtils_->getFinalDecisionByName(l1Paths_[iseed], l1htbit);
    
    l1Names_ ->push_back(l1Paths_[iseed]);
    l1Accept_->push_back( l1htbit );
  }
  
  //fill L1 jet candidates (not available for all events)
  edm::Handle<l1t::JetBxCollection> l1cands;
  event.getByToken(l1CandToken_, l1cands);
  if(l1cands.isValid())
    fillL1Jets(l1cands, event);

  //fill Jets
  edm::Handle<reco::CaloJetCollection> caloJets;
  event.getByToken(caloJetToken_,caloJets);
  
  edm::Handle<reco::PFJetCollection> PFJets;
  event.getByToken(PFJetToken_,PFJets);

  // not available for every event
  if(caloJets.isValid())
    fillMjj(caloJets,caloMjj_,caloJet1Pt_,caloJet2Pt_,caloJet1Eta_,caloJet2Eta_,caloJet1Phi_,caloJet2Phi_,caloDeltaEta_);
  
  if(PFJets.isValid())
    fillMjj(PFJets,PFMjj_,PFJet1Pt_,PFJet2Pt_,PFJet1Eta_,PFJet2Eta_,PFJet1Phi_,PFJet2Phi_, PFDeltaEta_);
  
  outTree_->Fill();  
}

DEFINE_FWK_MODULE(MyHLTAnalyzer);
