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
  void fillMjj(const edm::Handle<reco::CaloJetCollection> &);
  void beginEvent();


  
  /// file service and tree
  edm::Service<TFileService> outfile_;

  TTree* outTree_;
  int run_;
  int evt_;
  int lumi_;
  //  int nVtx_;
  float mjj_;
  //  float dEta_;
  //  float dPhi_;
  std::vector<int>* hltAccept_;
  std::vector<int>* hltWasrun_;
  std::vector<std::string>* hltNames_;

  // Trigger process
  edm::InputTag triggerResultTag_1_;
  edm::InputTag triggerResultTag_2_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_1_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_2_;
  edm::InputTag triggerEventTag_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> triggerEventToken_;
  edm::EDGetTokenT<reco::CaloJetCollection>         caloJetToken_;

  double minMass_;
  double fatJetDeltaR_;
  double maxDeltaEta_;
  double maxJetEta_;
  double minJetPt_;

  std::vector<std::string> hltPaths_;
};

MyHLTAnalyzer::MyHLTAnalyzer(const edm::ParameterSet& cfg): 
  triggerResultTag_1_       (cfg.getUntrackedParameter<edm::InputTag>("triggerResult_1")),
  triggerResultTag_2_       (cfg.getUntrackedParameter<edm::InputTag>("triggerResult_2")),
  triggerResultToken_1_     (consumes<edm::TriggerResults>(triggerResultTag_1_)),
  triggerResultToken_2_     (consumes<edm::TriggerResults>(triggerResultTag_2_)),
  triggerEventTag_          (cfg.getUntrackedParameter<edm::InputTag>("triggerSummary")),
  triggerEventToken_        (consumes<trigger::TriggerFilterObjectWithRefs>(triggerEventTag_)),
  caloJetToken_             (consumes<reco::CaloJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("caloJetTag"))),
  //params for wide jet calculation
  minMass_                  (cfg.getUntrackedParameter<double>("minMass")),
  fatJetDeltaR_             (cfg.getUntrackedParameter<double>("fatJetDeltaR")),
  maxDeltaEta_              (cfg.getUntrackedParameter<double>("maxDeltaEta")),
  maxJetEta_                (cfg.getUntrackedParameter<double>("maxJetEta")),
  minJetPt_                 (cfg.getUntrackedParameter<double>("minJetPt")),
  
  hltPaths_                 (cfg.getUntrackedParameter<std::vector<std::string> >("hltPaths"))
{
}

void MyHLTAnalyzer::beginEvent()
{
  evt_  = -1;
  lumi_ = -1;
  run_  = -1;
  mjj_  = -1;

  hltNames_->clear();
  hltAccept_->clear();
  hltWasrun_->clear();
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




void MyHLTAnalyzer::fillMjj(const edm::Handle<reco::CaloJetCollection> & jets)
{
  // Selected jets
  reco::CaloJetCollection recojets;
  reco::CaloJetCollection::const_iterator i = jets->begin();
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

  //FILL HERE THE MJJ IN THE TUPLE
  mjj_ = fj1.mass();
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

  outTree_->Branch("mjj",     &mjj_,       "mjj_/F");
  //  outTree_->Branch("dEta",    &dEta_,      "dEta_/F");
  //  outTree_->Branch("dPhi",    &dPhi_,      "dPhi_/F");


  hltNames_ = new std::vector<std::string>;
  hltAccept_ = new std::vector<int>;
  hltWasrun_ = new std::vector<int>;

  outTree_->Branch("hltNames",  "std::vector<std::string>", &hltNames_);
  outTree_->Branch("hltAccept", "std::vector<int>", &hltAccept_);
  outTree_->Branch("hltWasrun", "std::vector<int>", &hltWasrun_);

}

void MyHLTAnalyzer::endJob() 
{
  delete hltNames_;
  delete hltAccept_;
  delete hltWasrun_;
}


void MyHLTAnalyzer::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) 
{
  beginEvent();
  
  evt_             = event.id().event();
  lumi_            = event.id().luminosityBlock();
  run_             = event.id().run();
  

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

  edm::Handle<reco::CaloJetCollection> caloJets;
  event.getByToken(caloJetToken_,caloJets);

  // not available for every event
  if(caloJets.isValid())
    fillMjj(caloJets);

  outTree_->Fill();  
}

DEFINE_FWK_MODULE(MyHLTAnalyzer);
