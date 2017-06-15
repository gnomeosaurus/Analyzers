#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"


#include "DataFormats/EgammaReco/interface/SuperCluster.h" //added because of SuperCluster
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" //added because of SuperCluster (collection)
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" //added because of GsfElectron
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h" //added because of GsfElectron -- QUESTION - do I need to include something similar also for GsfElectronunclean?
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h" 
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"  
#include "DataFormats/METReco/interface/PFMETFwd.h" 
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h" //not entirely sure here
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"

#include <numeric>
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Event.h"

#include <map>
#include <string>
#include <fstream>
#include <iomanip>
#include "TTree.h"
#include "Python.h"



class AODAnalyzer : public edm::EDAnalyzer {
  
public:
  AODAnalyzer(const edm::ParameterSet& cfg);
  virtual ~AODAnalyzer();
  
  virtual void analyze (const edm::Event& event, const edm::EventSetup & eventSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  virtual void endRun  (const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  virtual void beginLuminosityBlock  (const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  virtual void endLuminosityBlock    (const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  
private:
  
  template<typename jetCollection>
  void fillJets(const edm::Handle<jetCollection> &, std::string );

  template<typename PFMETCollection>
  void fillPFChMets(const edm::Handle<PFMETCollection> &);

  template<typename PFMETCollection>
  void fillPFMets(const edm::Handle<PFMETCollection> &);

  template<typename SuperClusterCollection> 
  void fillSC(const edm::Handle<SuperClusterCollection> &); 
  //TODO
  template<typename GsfElectronCollection>
  void fillGsf(const edm::Handle<GsfElectronCollection> &);

  template<typename GsfElectronCollection>
  void fillUNGsf(const edm::Handle<GsfElectronCollection> &);

  template<typename EcalRecHitCollection>
  void fillEBrecHit(const edm::Handle<EcalRecHitCollection> &);

  template<typename EcalRecHitCollection>
  void fillEErecHit(const edm::Handle<EcalRecHitCollection> &);
 
  template<typename EcalRecHitCollection>
  void fillESrecHit(const edm::Handle<EcalRecHitCollection> &);
   //TODO


  void initialize();
  template<typename T>
  void computeQuantiles(std::vector<T>*, std::vector<T>*, std::vector<double>);
  template<typename T>
  void computeMeanAndRms(std::vector<T>*, std::vector<T>*);
  std::map<int, std::vector<std::pair<int, int> > > readJSONFile(const std::string&);
  bool AcceptEventByRunAndLumiSection(const int& , const int& ,std::map<int, std::vector<std::pair<int, int> > >&);




  /// file service and tree
  edm::Service<TFileService> outfile_;

  TTree* outTree_;
  int    runId_;
  int    lumiId_;
  float  lumi_;
  int    isSig_;

  std::vector<float>* PFJetPt_;
  std::vector<float>* PFJetEta_;
  std::vector<float>* PFJetPhi_;
  std::vector<float>* PFChMetPt_;
  std::vector<float>* PFChMetPhi_;
  std::vector<float>* PFMetPt_;
  std::vector<float>* PFMetPhi_;
  std::vector<int>*   nVtx_;

  //std::vector<float>* SuperCluster_; // adding SuperCluster
  std::vector<double>* SCEn_;
  std::vector<double>* SCEta_;
  std::vector<double>* SCPhi_;
  //TODO
  std::vector<float>* SigmaIEta_;
  std::vector<float>* SigmaIPhi_;
  std::vector<float>* r9_;
  std::vector<float>* HadOEm_;
  std::vector<float>* drSumPt_;
  std::vector<float>* drSumEt_;
  std::vector<float>* eSCOP_;
  std::vector<float>* ecEn_;

  std::vector<float>* UNSigmaIEta_;
  std::vector<float>* UNSigmaIPhi_;
  std::vector<float>* UNr9_;
  std::vector<float>* UNHadOEm_;
  std::vector<float>* UNdrSumPt_;
  std::vector<float>* UNdrSumEt_;
  std::vector<float>* UNeSCOP_;
  std::vector<float>* UNecEn_;

  std::vector<float>* EBenergy_;
  std::vector<float>* EBtime_;
  std::vector<float>* EBchi2_;
  std::vector<float>* EEenergy_;
  std::vector<float>* EEtime_;
  std::vector<float>* EEchi2_;
  std::vector<float>* ESenergy_;
  std::vector<float>* EStime_;
  std::vector<float>* ESchi2_;




 
  std::vector<float>* qPFJetPt_;
  std::vector<float>* qPFJetEta_;
  std::vector<float>* qPFJetPhi_;
  std::vector<float>* qPFChMetPt_;
  std::vector<float>* qPFChMetPhi_;
  std::vector<float>* qPFMetPt_;
  std::vector<float>* qPFMetPhi_;
  std::vector<int>*   qNVtx_;

  std::vector<double>* qSCEn_;
  std::vector<double>* qSCEta_;
  std::vector<double>* qSCPhi_;
  //TODO
  std::vector<float>* qSigmaIEta_;
  std::vector<float>* qSigmaIPhi_;
  std::vector<float>* qr9_;
  std::vector<float>* qHadOEm_;
  std::vector<float>* qdrSumPt_;
  std::vector<float>* qdrSumEt_;
  std::vector<float>* qeSCOP_;
  std::vector<float>* qecEn_;
 
  std::vector<float>* qUNSigmaIEta_;
  std::vector<float>* qUNSigmaIPhi_;
  std::vector<float>* qUNr9_;
  std::vector<float>* qUNHadOEm_;
  std::vector<float>* qUNdrSumPt_;
  std::vector<float>* qUNdrSumEt_;
  std::vector<float>* qUNeSCOP_;
  std::vector<float>* qUNecEn_;

  std::vector<float>* qEBenergy_;
  std::vector<float>* qEBtime_;
  std::vector<float>* qEBchi2_;
  std::vector<float>* qEEenergy_;
  std::vector<float>* qEEtime_;
  std::vector<float>* qEEchi2_;
  std::vector<float>* qESenergy_;
  std::vector<float>* qEStime_;
  std::vector<float>* qESchi2_;

  std::vector<float>*   crossSection_;
  std::vector<float>*   pathRates_;
  std::vector<std::string>*   pathNames_;
  std::map<std::string,int> rateMap;


  edm::EDGetTokenT<reco::PFJetCollection>    PFJetToken_;
  edm::EDGetTokenT<reco::PFMETCollection> PFChMETToken_;
  edm::EDGetTokenT<reco::PFMETCollection> PFMETToken_;
  //edm::EDGetTokenT<std::vector<reco::CaloMETCollection> >    CaloMETJetToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerPrescales_;  //pat::PackedTriggerPrescales
  edm::EDGetTokenT<reco::SuperClusterCollection>   SuperClusterToken_;  //adding SuperCluster

  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronUncleanedToken_;
  edm::EDGetTokenT<reco::MuonCollection> MuonToken_;
  edm::EDGetTokenT<reco::PhotonCollection> gedPhotonToken_;
  edm::EDGetTokenT<reco::PhotonCollection> PhotonToken_;   //TWO types of Photons -- one collection?

  
  //TODO -- finish adding collections!  // edm::   //eCALlASERcORRFILTER like

  // edm:: 

  edm::EDGetTokenT<EcalRecHitCollection> ebRHSrcToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRHSrcToken_;
  edm::EDGetTokenT<EcalRecHitCollection> esRHSrcToken_;


  // edm::EDGetTokenT<SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> ebRHSrcToken_;
  // edm::EDGetTokenT<SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> eeRHSrcToken_;





  int eventCounter;

  double maxJetEta_;
  double minJetPt_;
  double maxSCEta_;
  double minSCEn_;

  std::string lumiFile_;
  std::map<int,std::map<int,float> > lumiMap;

  std::vector<double> quantiles_;

  std::vector<std::string> subsystemNames_;
  std::vector<bool>* subsystemQuality_;

  std::vector<std::string> qualityFiles_;

  typedef std::map<int, std::vector<std::pair<int, int> > > jsonMapType;
  jsonMapType myJsonMap;
  std::map<std::string,jsonMapType> qualityMaps;
};



AODAnalyzer::AODAnalyzer(const edm::ParameterSet& cfg): 
  PFJetToken_               (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJetTag"))),  //reco instead of pat
  PFChMETToken_             (consumes<reco::PFMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFChMETTag"))),
  PFMETToken_               (consumes<reco::PFMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFMETTag"))),
  // CaloMETJetToken_          (consumes<reco::CaloMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloMETTag"))),  ///std::vector<reco::CaloMET> 
  vtxToken_                 (consumes<reco::VertexCollection>(cfg.getUntrackedParameter<edm::InputTag>("vtx"))),
  triggerBits_              (consumes<edm::TriggerResults>(cfg.getUntrackedParameter<edm::InputTag>("bits"))),
  triggerPrescales_         (consumes<trigger::TriggerEvent>(cfg.getUntrackedParameter<edm::InputTag>("prescales"))),  // pat::PackedTriggerPrescales
  SuperClusterToken_        (consumes<reco::SuperClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("SuperClusterTag"))), //adding SuperClusterToken_
  GsfElectronToken_         (consumes<reco::GsfElectronCollection>(cfg.getUntrackedParameter<edm::InputTag>("GsfElectronTag"))),
  GsfElectronUncleanedToken_(consumes<reco::GsfElectronCollection>(cfg.getUntrackedParameter<edm::InputTag>("GsfElectronUncleanedTag"))),
  MuonToken_                (consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonTag"))),
  gedPhotonToken_           (consumes<reco::PhotonCollection>(cfg.getUntrackedParameter<edm::InputTag>("gedPhotonTag"))),
  PhotonToken_              (consumes<reco::PhotonCollection>(cfg.getUntrackedParameter<edm::InputTag>("PhotonTag"))),
  
  
  //TODO -- add collections!
  //CollectionEcalRecHitEBTToken_
  ebRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("EBRecHitSourceTag"))),  //ICONFIG -> cfg
  eeRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("EERecHitSourceTag"))),
  esRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("ESRecHitSourceTag"))),





  //params for wide jet calculation
  maxJetEta_                (cfg.getUntrackedParameter<double>("maxJetEta")),
  minJetPt_                 (cfg.getUntrackedParameter<double>("minJetPt")),
  maxSCEta_                 (cfg.getUntrackedParameter<double>("maxSCEta")),
  minSCEn_                  (cfg.getUntrackedParameter<double>("minSCEn")),
  lumiFile_                 (cfg.getUntrackedParameter<std::string>("lumiFile")),
  quantiles_                (cfg.getUntrackedParameter<std::vector<double> >("quantiles")),
  subsystemNames_           (cfg.getUntrackedParameter<std::vector<std::string> >("subsystems")),
  qualityFiles_             (cfg.getUntrackedParameter<std::vector<std::string> >("qualityFiles"))


{}
void AODAnalyzer::initialize()
{
  eventCounter = 0;

  lumiId_ = -1;
  lumi_  = -1;
  runId_ = -1;
  isSig_ = -1;

  PFJetPt_->clear();
  PFJetEta_->clear();
  PFJetPhi_->clear();
  PFChMetPt_->clear();
  PFChMetPhi_->clear();
  PFMetPt_->clear();
  PFMetPhi_->clear();
  nVtx_->clear();
  //SuperCluster_ ->clear(); //adding SuperCluster
  SCEn_ ->clear();
  SCEta_->clear();
  SCPhi_->clear();
  //TODO
  SigmaIEta_->clear();
  SigmaIPhi_->clear();
  r9_->clear();
  HadOEm_->clear();
  drSumPt_->clear();
  drSumEt_->clear();
  eSCOP_->clear();
  ecEn_->clear();

  UNSigmaIEta_->clear();
  UNSigmaIPhi_->clear();
  UNr9_->clear();
  UNHadOEm_->clear();
  UNdrSumPt_->clear();
  UNdrSumEt_->clear();
  UNeSCOP_->clear();
  UNecEn_->clear();
  
  EBenergy_->clear();
  EBtime_->clear();
  EBchi2_->clear();
  EEenergy_->clear();
  EEtime_->clear();
  EEchi2_->clear();
  ESenergy_->clear();
  EStime_->clear();
  ESchi2_->clear();
 
  qPFJetPt_->clear();
  qPFJetEta_->clear();
  qPFJetPhi_->clear();
  qPFChMetPt_->clear();
  qPFChMetPhi_->clear();
  qPFMetPt_->clear();
  qPFMetPhi_->clear();
  qNVtx_->clear();

  //qSuperCluster_ ->clear(); //ASK IF quantiles should be here
  qSCEn_ ->clear();
  qSCEta_->clear();
  qSCPhi_->clear();
  //TODO
  qSigmaIEta_->clear();
  qSigmaIPhi_->clear();  
  qr9_->clear();
  qHadOEm_->clear();
  qdrSumPt_->clear();
  qdrSumEt_->clear();
  qeSCOP_->clear();
  qecEn_->clear();

  qUNSigmaIEta_->clear();
  qUNSigmaIPhi_->clear();  
  qUNr9_->clear();
  qUNHadOEm_->clear();
  qUNdrSumPt_->clear();
  qUNdrSumEt_->clear();
  qUNeSCOP_->clear();
  qUNecEn_->clear();

  qEBenergy_->clear();
  qEBtime_->clear();
  qEBchi2_->clear();
  qEEenergy_->clear();
  qEEtime_->clear();
  qEEchi2_->clear();
  qESenergy_->clear();
  qEStime_->clear();
  qESchi2_->clear();

  crossSection_->clear();
  pathRates_->clear();
  pathNames_->clear();
  
  rateMap.clear();

  subsystemQuality_->clear();
}

template<typename jetCollection>
void AODAnalyzer::fillJets(const edm::Handle<jetCollection> & jets, std::string type)
{
  // Selected jets
  //reco::CaloJetCollection recojets;
  std::cout << "fillJets is called" <<std::endl;
  typename jetCollection::const_iterator i = jets->begin();
  for(;i != jets->end(); i++){
    if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
      {
	
	if(type.compare(std::string("PF")) == 0)
	  {
	    PFJetPt_->push_back(i->pt());
	    PFJetEta_->push_back(i->eta());
	    PFJetPhi_->push_back(i->phi());
      std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
	  }
	
      }
  }
  return;
}

template<typename PFMETCollection>
void AODAnalyzer::fillPFChMets(const edm::Handle<PFMETCollection> & pfchmets)
{
  std::cout << "fillPFChMets is being called!" << std::endl;
  std::cout << pfchmets->size() <<std::endl;
  typename PFMETCollection::const_iterator i = pfchmets->begin();
  for(;i != pfchmets->end(); i++){
    PFChMetPt_->push_back(i->et());
    PFChMetPhi_->push_back(i->phi());

  }
  return;

}

template<typename PFMETCollection>
void AODAnalyzer::fillPFMets(const edm::Handle<PFMETCollection> & pfmets)
{
  std::cout << "fillPFChMets is being called!" << std::endl;
  std::cout << pfmets->size() <<std::endl;
  typename PFMETCollection::const_iterator i = pfmets->begin();
  for(;i != pfmets->end(); i++){
    PFMetPt_->push_back(i->et());
    PFMetPhi_->push_back(i->phi());

  }
  return;

}
template<typename SuperClusterCollection>
void AODAnalyzer::fillSC(const edm::Handle<SuperClusterCollection> & superclusters) //ask for jets analogy //SUPERCLUSTERS
{

  // Selected jets
  //reco::CaloJetCollection recojets;
  typename SuperClusterCollection::const_iterator i = superclusters->begin();
  for(;i != superclusters->end(); i++){
     if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // not sure if needed
      // {
        SCEn_->push_back(i->energy());
        SCEta_->push_back(i->etaWidth());
        SCPhi_->push_back(i->phiWidth());

        std::cout << "ele energy: " << i->energy()   << std::endl; 
        std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;



}

//TODO
template<typename GsfElectronCollection>
void AODAnalyzer::fillGsf(const edm::Handle<GsfElectronCollection> & electrons)
{
  std::cout << "fillGSF is being called!" << std::endl;
  std::cout << electrons->size() <<std::endl;
  typename GsfElectronCollection::const_iterator i = electrons->begin();
  for(;i != electrons->end(); i++){
    SigmaIEta_->push_back(i->sigmaIetaIeta());
    SigmaIPhi_->push_back(i->sigmaIphiIphi());
    r9_       ->push_back(i->r9());
    HadOEm_   ->push_back(i->hadronicOverEm());
    drSumPt_  ->push_back(i->dr03TkSumPt());
    drSumEt_  ->push_back(i->dr03EcalRecHitSumEt());
    eSCOP_    ->push_back(i->eSuperClusterOverP());
    ecEn_     ->push_back(i->ecalEnergy());
    // std::cout << "ele etaeta: "  << i->sigmaIetaIeta()       << std::endl;
    // std::cout << "ele phiphi: "  << i->sigmaIphiIphi()       << std::endl; 
    // std::cout << "ele r9: "      << i->r9()                  << std::endl; 
    // std::cout << "ele hadoem: "  << i->hadronicOverEm()      << std::endl; 
    // std::cout << "ele drsumpt: " << i->dr03TkSumPt()         << std::endl; 
    // std::cout << "ele drsumet: " << i->dr03EcalRecHitSumEt() << std::endl; 
    // std::cout << "ele escop: "   << i->eSuperClusterOverP()  << std::endl; 
    // std::cout << "ele ecen: "    << i->ecalEnergy()          << std::endl; 


  }
  return;

}

template<typename GsfElectronCollection>
void AODAnalyzer::fillUNGsf(const edm::Handle<GsfElectronCollection> & UNelectrons)
{
  
  std::cout << "fillUNGSF is being called!" << std::endl;
  typename GsfElectronCollection::const_iterator i = UNelectrons->begin();
  for(;i != UNelectrons->end(); i++){
    UNSigmaIEta_->push_back(i->sigmaIetaIeta());
    UNSigmaIPhi_->push_back(i->sigmaIphiIphi());
    UNr9_       ->push_back(i->r9());
    UNHadOEm_   ->push_back(i->hadronicOverEm());
    UNdrSumPt_  ->push_back(i->dr03TkSumPt());
    UNdrSumEt_  ->push_back(i->dr03EcalRecHitSumEt());
    UNeSCOP_    ->push_back(i->eSuperClusterOverP());
    UNecEn_     ->push_back(i->ecalEnergy());
    // std::cout << "ele UNetaeta: "  << i->sigmaIetaIeta()       << std::endl;
    // std::cout << "ele UNphiphi: "  << i->sigmaIphiIphi()       << std::endl; 
    // std::cout << "ele UNr9: "      << i->r9()                  << std::endl; 
    // std::cout << "ele UNhadoem: "  << i->hadronicOverEm()      << std::endl; 
    // std::cout << "ele UNdrsumpt: " << i->dr03TkSumPt()         << std::endl; 
    // std::cout << "ele UNdrsumet: " << i->dr03EcalRecHitSumEt() << std::endl; 
    // std::cout << "ele UNescop: "   << i->eSuperClusterOverP()  << std::endl; 
    // std::cout << "ele UNecen: "    << i->ecalEnergy()          << std::endl; 


  }
  return;

}

template<typename EcalRecHitCollection>
void AODAnalyzer::fillEBrecHit(const edm::Handle<EcalRecHitCollection> & EBhits)
{

  std::cout << "fillEBrecHit is being called!" << std::endl;
  typename EcalRecHitCollection::const_iterator i = EBhits->begin();
  for(;i != EBhits->end(); i++){
    EBenergy_ ->push_back(i->energy());
    EBtime_ ->push_back(i->time());
    EBchi2_ ->push_back(i->chi2());

  }
  return;
}

template<typename EcalRecHitCollection>
void AODAnalyzer::fillEErecHit(const edm::Handle<EcalRecHitCollection> & EEhits)
{

  std::cout << "fillEErecHit is being called!" << std::endl;
  typename EcalRecHitCollection::const_iterator i = EEhits->begin();
  for(;i != EEhits->end(); i++){
    EEenergy_ ->push_back(i->energy());
    EEtime_ ->push_back(i->time());
    EEchi2_ ->push_back(i->chi2());

  }
  return;
}

template<typename EcalRecHitCollection>
void AODAnalyzer::fillESrecHit(const edm::Handle<EcalRecHitCollection> & EShits)
{

  std::cout << "fillESrecHit is being called!" << std::endl;
  typename EcalRecHitCollection::const_iterator i = EShits->begin();
  for(;i != EShits->end(); i++){
    ESenergy_ ->push_back(i->energy());
    EStime_ ->push_back(i->time());
    ESchi2_ ->push_back(i->chi2());

  }
  return;
}


template<typename T>
void AODAnalyzer::computeQuantiles(std::vector<T>* myDistr, std::vector<T>* myQuan, std::vector<double> qq)
{
  //need to sort the distr to compute quantiles
  std::vector<T> dummyDistr = *myDistr;
  std::sort(dummyDistr.begin(), dummyDistr.end());

  //scan the vector and find quantiles
  unsigned int qItr = 0;
  for (unsigned int itr = 0; itr < dummyDistr.size(); ++itr)
  {
    if(qItr >= qq.size())
      return;

    float prob = ((float)itr+1.)/(float)dummyDistr.size();
    if( prob >= qq[qItr] )
	  {
  	  //fill root tree
  	  myQuan->push_back(dummyDistr.at(itr));
  	  ++qItr;
  	}
  }
  return;
}


template<typename T>
void AODAnalyzer::computeMeanAndRms(std::vector<T>* myDistr, std::vector<T>* myVect)
{
  double sum = std::accumulate(myDistr->begin(), myDistr->end(), 0.0);
  double mean = sum / myDistr->size();
  myVect->push_back( mean );

  std::vector<double> diff(myDistr->size());
  std::transform(myDistr->begin(), myDistr->end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  myVect->push_back( std::sqrt(sq_sum / myDistr->size()) );
}



//===================== beginJob and Analyze =========================
void AODAnalyzer::beginJob() {

  TH1::SetDefaultSumw2() ;
  outTree_ = outfile_-> make<TTree>("MyTree","MyTree");

  outTree_->Branch("runId",     &runId_,       "runId_/I");
  outTree_->Branch("lumiId",    &lumiId_,      "lumiId_/I");
  outTree_->Branch("lumi",      &lumi_,        "lumi_/F");
  outTree_->Branch("isSig",     &isSig_,       "isSig_/I");

  PFJetPt_   = new std::vector<float>;
  PFJetEta_  = new std::vector<float>;
  PFJetPhi_  = new std::vector<float>;
  PFChMetPt_     = new std::vector<float>;
  PFChMetPhi_    = new std::vector<float>;
  PFMetPt_     = new std::vector<float>;
  PFMetPhi_    = new std::vector<float>;
  nVtx_      = new std::vector<int>;
  
  //SuperCluster_ = new std::vector<float>; //adding SuperCluster
  SCEn_      = new std::vector<double>;
  SCEta_     = new std::vector<double>;
  SCPhi_     = new std::vector<double>;
  //TODO
  SigmaIEta_ = new std::vector<float>;
  SigmaIPhi_ = new std::vector<float>;
  r9_        = new std::vector<float>;
  HadOEm_    = new std::vector<float>;
  drSumPt_   = new std::vector<float>;
  drSumEt_   = new std::vector<float>;
  eSCOP_     = new std::vector<float>;
  ecEn_      = new std::vector<float>;
 
  UNSigmaIEta_ = new std::vector<float>;
  UNSigmaIPhi_ = new std::vector<float>;
  UNr9_        = new std::vector<float>;
  UNHadOEm_    = new std::vector<float>;
  UNdrSumPt_   = new std::vector<float>;
  UNdrSumEt_   = new std::vector<float>;
  UNeSCOP_     = new std::vector<float>;
  UNecEn_      = new std::vector<float>;

  EBenergy_    = new std::vector<float>;
  EBtime_      = new std::vector<float>;
  EBchi2_      = new std::vector<float>;
  EEenergy_    = new std::vector<float>;
  EEtime_      = new std::vector<float>;
  EEchi2_      = new std::vector<float>;
  ESenergy_    = new std::vector<float>;
  EStime_      = new std::vector<float>;
  ESchi2_      = new std::vector<float>;

  // outTree_->Branch("MetPt",     "std::vector<std::float>",     &MetPt_);
  // outTree_->Branch("MetPhi",    "std::vector<std::float>",     &MetPhi_);
  // outTree_->Branch("PFJetPt",     "std::vector<std::float>",     &PFJetPt_);
  // outTree_->Branch("PFJetEta",    "std::vector<std::float>",     &PFJetEta_);
  // outTree_->Branch("PFJetPhi",    "std::vector<std::float>",     &PFJetPhi_);
  // outTree_->Branch("nVtx",           "std::vector<std::int>",       &nVtx_);

  qPFJetPt_  = new std::vector<float>;
  qPFJetEta_ = new std::vector<float>;
  qPFJetPhi_ = new std::vector<float>;
  qPFChMetPt_     = new std::vector<float>;
  qPFChMetPhi_    = new std::vector<float>;
  qPFMetPt_     = new std::vector<float>;
  qPFMetPhi_    = new std::vector<float>;
 
  qSCEn_     = new std::vector<double>;
  qSCEta_    = new std::vector<double>;
  qSCPhi_    = new std::vector<double>;

  qSigmaIEta_= new std::vector<float>;
  qSigmaIPhi_= new std::vector<float>;
  qr9_       = new std::vector<float>;
  qHadOEm_   = new std::vector<float>;
  qdrSumPt_  = new std::vector<float>;
  qdrSumEt_  = new std::vector<float>;
  qeSCOP_    = new std::vector<float>;
  qecEn_     = new std::vector<float>;

  qUNSigmaIEta_= new std::vector<float>;
  qUNSigmaIPhi_= new std::vector<float>;
  qUNr9_       = new std::vector<float>;
  qUNHadOEm_   = new std::vector<float>;
  qUNdrSumPt_  = new std::vector<float>;
  qUNdrSumEt_  = new std::vector<float>;
  qUNeSCOP_    = new std::vector<float>;
  qUNecEn_     = new std::vector<float>;
  //TODO

  qEBenergy_   = new std::vector<float>;
  qEBtime_     = new std::vector<float>;
  qEBchi2_     = new std::vector<float>;
  qEEenergy_   = new std::vector<float>;
  qEEtime_     = new std::vector<float>;
  qEEchi2_     = new std::vector<float>;
  qESenergy_   = new std::vector<float>;
  qEStime_     = new std::vector<float>;
  qESchi2_     = new std::vector<float>;



  qNVtx_       = new std::vector<int>;
  crossSection_= new std::vector<float>;
  pathRates_   = new std::vector<float>;
  pathNames_   = new std::vector<std::string>;
  outTree_->Branch("qPFJetPt",     "std::vector<std::float>",      &qPFJetPt_);
  outTree_->Branch("qPFJetEta",    "std::vector<std::float>",      &qPFJetEta_);
  outTree_->Branch("qPFJetPhi",    "std::vector<std::float>",      &qPFJetPhi_);
  outTree_->Branch("qPFChMetPt",     "std::vector<std::float>",        &qPFChMetPt_);
  outTree_->Branch("qPFChMetPhi",    "std::vector<std::float>",        &qPFChMetPhi_);
  outTree_->Branch("qPFMetPt",     "std::vector<std::float>",        &qPFMetPt_);
  outTree_->Branch("qPFMetPhi",    "std::vector<std::float>",        &qPFMetPhi_);
  outTree_->Branch("qNVtx",        "std::vector<std::int>",        &qNVtx_);

  //outTree_->Branch("SuperCluster","std::vector<std::float>",      &SuperCluster_);
  outTree_->Branch("qSCEn",     "std::vector<std::double>",        &qSCEn_);
  outTree_->Branch("qSCEta",    "std::vector<std::double>",        &qSCEta_);
  outTree_->Branch("qSCPhi",    "std::vector<std::double>",        &qSCPhi_);

  outTree_->Branch("qSigmaIEta",    "std::vector<std::float>",     &qSigmaIEta_);
  outTree_->Branch("qSigmaIPhi",    "std::vector<std::float>",     &qSigmaIPhi_);
  outTree_->Branch("qr9",    "std::vector<std::float>",            &qr9_);
  outTree_->Branch("qHadOEm",    "std::vector<std::float>",        &qHadOEm_);
  outTree_->Branch("qdrSumPt",    "std::vector<std::float>",       &qdrSumPt_);
  outTree_->Branch("qdrSumEt",    "std::vector<std::float>",       &qdrSumEt_);
  outTree_->Branch("qeSCOP",    "std::vector<std::float>",         &qeSCOP_);
  outTree_->Branch("qecEn",    "std::vector<std::float>",          &qecEn_);

  outTree_->Branch("qUNSigmaIEta",    "std::vector<std::float>",     &qUNSigmaIEta_);
  outTree_->Branch("qUNSigmaIPhi",    "std::vector<std::float>",     &qUNSigmaIPhi_);
  outTree_->Branch("qUNr9",    "std::vector<std::float>",            &qUNr9_);
  outTree_->Branch("qUNHadOEm",    "std::vector<std::float>",        &qUNHadOEm_);
  outTree_->Branch("qUNdrSumPt",    "std::vector<std::float>",       &qUNdrSumPt_);
  outTree_->Branch("qUNdrSumEt",    "std::vector<std::float>",       &qUNdrSumEt_);
  outTree_->Branch("qUNeSCOP",    "std::vector<std::float>",         &qUNeSCOP_);
  outTree_->Branch("qUNecEn",    "std::vector<std::float>",          &qUNecEn_);
  //TODO
  outTree_->Branch("qEBenergy",    "std::vector<std::float>",        &qEBenergy_);
  outTree_->Branch("qEBtime",    "std::vector<std::float>",          &qEBtime_);
  outTree_->Branch("qEBchi2",    "std::vector<std::float>",          &qEBchi2_);
  outTree_->Branch("qEEenergy",    "std::vector<std::float>",        &qEEenergy_);
  outTree_->Branch("qEEtime",    "std::vector<std::float>",          &qEEtime_);
  outTree_->Branch("qEEchi2",    "std::vector<std::float>",          &qEEchi2_);
  outTree_->Branch("qESenergy",    "std::vector<std::float>",        &qESenergy_);
  outTree_->Branch("qEStime",    "std::vector<std::float>",          &qEStime_);
  outTree_->Branch("qESchi2",    "std::vector<std::float>",          &qESchi2_);



  outTree_->Branch("crossSection",   "std::vector<std::float>",    &crossSection_);
  outTree_->Branch("pathRates",      "std::vector<std::float>",    &pathRates_);
  outTree_->Branch("pathNames",      "std::vector<std::string>",   &pathNames_);

  subsystemQuality_ = new std::vector<bool>;
  outTree_->Branch("subsystemQuality", "std::vector<bool>",        &subsystemQuality_);
  outTree_->Branch("subsystemNames",   "std::vector<std::string>", &subsystemNames_);


  //open and load luminosity file
  std::ifstream infile;
  infile.open(lumiFile_);
  while(!infile.eof())
    {
      int run, ls;
      float lumi;
      infile >> run >> ls >> lumi;
      lumiMap[run][ls] = lumi;
    }
  infile.close();

  // Load the json files
  for(unsigned int itr = 0; itr < subsystemNames_.size(); ++itr)
    qualityMaps[subsystemNames_.at(itr)] = readJSONFile(qualityFiles_[itr]);

}

void AODAnalyzer::endJob() 
{

  delete PFJetPt_;
  delete PFJetEta_;
  delete PFJetPhi_;
  delete PFChMetPt_;
  delete PFChMetPhi_;
  delete PFMetPt_;
  delete PFMetPhi_;

  delete SCEn_;
  delete SCEta_;
  delete SCPhi_;
  //TODO
  delete SigmaIEta_;
  delete SigmaIPhi_;
  delete r9_;
  delete HadOEm_;
  delete drSumPt_;
  delete drSumEt_;
  delete eSCOP_;
  delete ecEn_;

  delete UNSigmaIEta_;
  delete UNSigmaIPhi_;
  delete UNr9_;
  delete UNHadOEm_;
  delete UNdrSumPt_;
  delete UNdrSumEt_;
  delete UNeSCOP_;
  delete UNecEn_;

  delete EBenergy_;
  delete EBtime_;
  delete EBchi2_;
  delete EEenergy_;
  delete EEtime_;
  delete EEchi2_;
  delete ESenergy_;
  delete EStime_;
  delete ESchi2_;

  //delete SuperCluster_;
  //delete qSuperCluster_;
  delete qPFJetPt_;
  delete qPFJetEta_;
  delete qPFJetPhi_;
  delete qPFChMetPt_;
  delete qPFChMetPhi_;
  delete qPFMetPt_;
  delete qPFMetPhi_;
  delete qSCEn_;
  delete qSCEta_;
  delete qSCPhi_;
  //TODO
  delete qSigmaIEta_;
  delete qSigmaIPhi_;
  delete qr9_;
  delete qHadOEm_;
  delete qdrSumPt_;
  delete qdrSumEt_;
  delete qeSCOP_;
  delete qecEn_;

  delete qUNSigmaIEta_;
  delete qUNSigmaIPhi_;
  delete qUNr9_;
  delete qUNHadOEm_;
  delete qUNdrSumPt_;
  delete qUNdrSumEt_;
  delete qUNeSCOP_;
  delete qUNecEn_;

  delete qEBenergy_;
  delete qEBtime_;
  delete qEBchi2_;
  delete qEEenergy_;
  delete qEEtime_;
  delete qEEchi2_;
  delete qESenergy_;
  delete qEStime_;
  delete qESchi2_;

  delete qNVtx_;
  delete crossSection_;
  delete pathRates_;
  delete pathNames_;

  delete subsystemQuality_;

}

AODAnalyzer::~AODAnalyzer()
{

}

void AODAnalyzer::beginLuminosityBlock (const edm::LuminosityBlock & lumi, const edm::EventSetup &eventSetup) 
{
  initialize();

  lumiId_            = lumi.luminosityBlock();
  runId_             = lumi.run();
  lumi_              = lumiMap[runId_][lumiId_];

  //store quality info
  unsigned int totQuality = 0;
  for(unsigned int itr = 0; itr < subsystemNames_.size(); ++itr)
    {
      bool qq = AcceptEventByRunAndLumiSection(runId_, lumiId_, qualityMaps[subsystemNames_.at(itr)]);
      subsystemQuality_->push_back(qq);
      if (qq == true)
	++totQuality;
    }
  if (totQuality < subsystemNames_.size())
    isSig_ = 0;
  else
    isSig_ = 1;

}

void AODAnalyzer::endLuminosityBlock (const edm::LuminosityBlock & lumi, const edm::EventSetup &eventSetup) 
{

  //compute and store quantiles
  computeMeanAndRms(PFJetPt_, qPFJetPt_);
  computeMeanAndRms(PFJetEta_,qPFJetEta_);
  computeMeanAndRms(PFJetPhi_,qPFJetPhi_);
  computeMeanAndRms(PFChMetPt_, qPFChMetPt_);
  computeMeanAndRms(PFChMetPhi_,  qPFChMetPhi_);
  computeMeanAndRms(PFMetPt_, qPFMetPt_);
  computeMeanAndRms(PFMetPhi_,  qPFMetPhi_);
  computeMeanAndRms(nVtx_,    qNVtx_);
  computeMeanAndRms(SCEn_, qSCEn_);   //adding supercluster pt,eta and phi
  computeMeanAndRms(SCEta_, qSCEta_);  //
  computeMeanAndRms(SCPhi_, qSCPhi_);  //
  //TODO
  computeMeanAndRms(SigmaIEta_, qSigmaIEta_);
  computeMeanAndRms(SigmaIPhi_, qSigmaIPhi_);
  computeMeanAndRms(r9_, qr9_);
  computeMeanAndRms(HadOEm_, qHadOEm_);
  computeMeanAndRms(drSumPt_, qdrSumPt_);
  computeMeanAndRms(drSumEt_, qdrSumEt_);
  computeMeanAndRms(eSCOP_, qeSCOP_);
  computeMeanAndRms(ecEn_, qecEn_);

  computeMeanAndRms(UNSigmaIEta_, qUNSigmaIEta_);
  computeMeanAndRms(UNSigmaIPhi_, qUNSigmaIPhi_);
  computeMeanAndRms(UNr9_, qUNr9_);
  computeMeanAndRms(UNHadOEm_, qUNHadOEm_);
  computeMeanAndRms(UNdrSumPt_, qUNdrSumPt_);
  computeMeanAndRms(UNdrSumEt_, qUNdrSumEt_);
  computeMeanAndRms(UNeSCOP_, qUNeSCOP_);
  computeMeanAndRms(UNecEn_, qUNecEn_);
  computeMeanAndRms(EBenergy_, qEBenergy_);
  computeMeanAndRms(EBtime_, qEBtime_);
  computeMeanAndRms(EBchi2_, qEBchi2_);
  computeMeanAndRms(EEenergy_, qEEenergy_);
  computeMeanAndRms(EEtime_, qEEtime_);
  computeMeanAndRms(EEchi2_, qEEchi2_);
  computeMeanAndRms(ESenergy_, qESenergy_);
  computeMeanAndRms(EStime_, qEStime_);
  computeMeanAndRms(ESchi2_, qESchi2_);



  computeQuantiles(PFJetPt_, qPFJetPt_, quantiles_);
  computeQuantiles(PFJetEta_,qPFJetEta_,quantiles_);
  computeQuantiles(PFJetPhi_,qPFJetPhi_,quantiles_);
  computeQuantiles(PFChMetPt_, qPFChMetPt_,     quantiles_);
  computeQuantiles(PFChMetPhi_,qPFChMetPhi_,    quantiles_);
  computeQuantiles(PFMetPt_, qPFMetPt_,     quantiles_);
  computeQuantiles(PFMetPhi_,qPFMetPhi_,    quantiles_);
  computeQuantiles(nVtx_,    qNVtx_,    quantiles_);
  computeQuantiles(SCEn_, qSCEn_,       quantiles_);
  computeQuantiles(SCEta_, qSCEta_,     quantiles_);
  computeQuantiles(SCPhi_, qSCPhi_,     quantiles_);
  //TODO
  computeQuantiles(SigmaIEta_, qSigmaIEta_, quantiles_);
  computeQuantiles(SigmaIPhi_, qSigmaIPhi_, quantiles_);
  computeQuantiles(r9_, qr9_, quantiles_);
  computeQuantiles(HadOEm_, qHadOEm_, quantiles_);
  computeQuantiles(drSumPt_, qdrSumPt_, quantiles_);
  computeQuantiles(drSumEt_, qdrSumEt_, quantiles_);
  computeQuantiles(eSCOP_, qeSCOP_, quantiles_);
  computeQuantiles(ecEn_, qecEn_, quantiles_);

  computeQuantiles(UNSigmaIEta_, qUNSigmaIEta_, quantiles_);
  computeQuantiles(UNSigmaIPhi_, qUNSigmaIPhi_, quantiles_);
  computeQuantiles(UNr9_, qUNr9_, quantiles_);
  computeQuantiles(UNHadOEm_, qUNHadOEm_, quantiles_);
  computeQuantiles(UNdrSumPt_, qUNdrSumPt_, quantiles_);
  computeQuantiles(UNdrSumEt_, qUNdrSumEt_, quantiles_);
  computeQuantiles(UNeSCOP_, qUNeSCOP_, quantiles_);
  computeQuantiles(UNecEn_, qUNecEn_, quantiles_);

  computeQuantiles(EBenergy_, qEBenergy_, quantiles_);
  computeQuantiles(EBtime_, qEBtime_, quantiles_);
  computeQuantiles(EBchi2_, qEBchi2_, quantiles_);
  computeQuantiles(EEenergy_, qEEenergy_, quantiles_);
  computeQuantiles(EEtime_, qEEtime_, quantiles_);
  computeQuantiles(EEchi2_, qEEchi2_, quantiles_);
  computeQuantiles(ESenergy_, qESenergy_, quantiles_);
  computeQuantiles(EStime_, qEStime_, quantiles_);
  computeQuantiles(ESchi2_, qESchi2_, quantiles_);


  crossSection_->push_back( (float)eventCounter/lumi_ );
  
  for(std::map<std::string,int>::const_iterator itr = rateMap.begin(); itr != rateMap.end(); ++itr)
    {
      pathNames_->push_back(itr->first);
      pathRates_->push_back(itr->second/lumi_);
    }


  //fill tree one event per LS
  outTree_->Fill();


}


void AODAnalyzer::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) 
{
  ++eventCounter;


  //fill Jets
  edm::Handle<reco::PFJetCollection> PFJets;
  event.getByToken(PFJetToken_,PFJets);
  if(PFJets.isValid())
    fillJets(PFJets, std::string("PF")); 

  //   //fill PFChMet

  edm::Handle<reco::PFMETCollection> pfchmetlocalv;
  event.getByToken(PFChMETToken_, pfchmetlocalv);
  // 
  if(pfchmetlocalv.isValid())
    fillPFChMets(pfchmetlocalv);
  // edm::Handle<std::vector<reco::PFMETCollection> > PFChMet;    ////
  // event.getByToken(PFChMETToken_,PFChMet);
  // if(PFChMet.isValid())
  //   {
  //     PFChMetPt_->push_back( (*PFChMet)[0].et() );
  //     PFChMetPhi_->push_back( (*PFChMet)[0].phi() );
  //   }

  // //fill PFMet
  edm::Handle<reco::PFMETCollection> pfmetlocalv;
  event.getByToken(PFMETToken_, pfmetlocalv);
  // 
  if(pfmetlocalv.isValid())
    fillPFMets(pfmetlocalv);
  // edm::Handle<std::vector<reco::PFMETCollection> > PFMet;    ////
  // event.getByToken(PFMETToken_,PFMet);
  // if(PFMet.isValid())
  //   {
  //     PFMetPt_->push_back( (*PFMet)[0].et() );
  //     PFMetPhi_->push_back( (*PFMet)[0].phi() );
  //   } 

  //fill vtx
  edm::Handle<reco::VertexCollection> recVtxs;
  event.getByToken(vtxToken_,recVtxs);
  if(recVtxs.isValid())
    nVtx_->push_back(recVtxs->size());

  //Fill SuperCluster
  edm::Handle<reco::SuperClusterCollection> SuperClusterlocalv;
  event.getByToken(SuperClusterToken_, SuperClusterlocalv);
  // print the size of SuperClusterlocalv
  if(SuperClusterlocalv.isValid())
    fillSC(SuperClusterlocalv);
  

  //TODO --fill Photons, fill Muons
  //fill GsF
  edm::Handle<reco::GsfElectronCollection> GsfElectronlocalv;
  event.getByToken(GsfElectronToken_, GsfElectronlocalv);
  if(GsfElectronlocalv.isValid())
     fillGsf(GsfElectronlocalv);

  //fill UNGsf
  edm::Handle<reco::GsfElectronCollection> UNGsfElectronlocalv;
  event.getByToken(GsfElectronUncleanedToken_, UNGsfElectronlocalv);
  if(UNGsfElectronlocalv.isValid())
     fillUNGsf(UNGsfElectronlocalv);
  

  //fill EcalRec EB    
  edm::Handle<EcalRecHitCollection> ebRHs;    //finish fill ee, eb and so on
  event.getByToken(ebRHSrcToken_, ebRHs);
  if(ebRHs.isValid())
    fillEBrecHit(ebRHs);



  //fill EcalRec EE
  edm::Handle<EcalRecHitCollection> eeRHs;
  event.getByToken(eeRHSrcToken_, eeRHs);
  if(eeRHs.isValid())
    fillEErecHit(eeRHs);

  //fill EcalRec ES
  edm::Handle<EcalRecHitCollection> esRHs;
  event.getByToken(esRHSrcToken_, esRHs);
  if(esRHs.isValid())
    fillESrecHit(esRHs);

  //Lorentz vector
   // Look for MET --LORENTZ VECTOR
  // edm::Handle<reco::GenMETCollection> genMet;
  // iEvent.getByToken(MCMET_, genMet);
  // LorentzVector MET(0.,0.,0.,0.);
  // if(genMet.isValid()){
  //   MET = LorentzVector(
  //       genMet->front().px(),
  //       genMet->front().py(),
  //       0,
  //       genMet->front().pt()
  //   );
  // }     
  // product_MET->push_back(MET);

  //  ANOTHER line I found --   reco::MET::LorentzVector p4(mhx, mhy, 0, sqrt(mhx*mhx + mhy*mhy));
  // ANOTHER LINE i found --    Particle::LorentzVector p,p1,p2;

  //   // Take the SA container       -- LORENTZ
  // LogTrace(metname)<<" Taking the StandAlone muons: "<<theSACollectionLabel;
  // Handle<TrackCollection> tracks; 
  // event.getByToken(tracksToken,tracks);

  // // Create a RecoChargedCandidate collection
  // LogTrace(metname)<<" Creating the RecoChargedCandidate collection";
  // auto candidates = std::make_unique<RecoChargedCandidateCollection>();

  // for (unsigned int i=0; i<tracks->size(); i++) {
  //     TrackRef tkref(tracks,i);
  //     Particle::Charge q = tkref->charge();
  //     Particle::LorentzVector p4(tkref->px(), tkref->py(), tkref->pz(), tkref->p());
  //     Particle::Point vtx(tkref->vx(),tkref->vy(), tkref->vz());
  //     int pid = 13;
  //     if(abs(q)==1) pid = q < 0 ? 13 : -13;
  //     else LogWarning(metname) << "L2MuonCandidate has charge = "<<q;
  //     RecoChargedCandidate cand(q, p4, vtx, pid);
  //     cand.setTrack(tkref);
  //     candidates->push_back(cand);
  // }                                -- LORENTZ


   ///        --lorentz
  // if(selectedMuons.size() == 4){
  //   reco::Candidate::LorentzVector p4CM;
  //   for (pat::MuonCollection::const_iterator muon = selectedMuons.begin();  muon != selectedMuons.end(); ++muon){
  //     p4CM = p4CM + muon->p4();
  //   }
  //   h4MuInvMass->Fill(p4CM.mass());
  // }    --lorentz

  //fill hlt
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<trigger::TriggerEvent> triggerPrescales;

  event.getByToken(triggerBits_, triggerBits);
  event.getByToken(triggerPrescales_, triggerPrescales);
  
 //  const edm::TriggerNames &names = event.triggerNames(*triggerBits);
 //  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) 
 //    {
 //      if(rateMap.find(names.triggerName(i)) != rateMap.end())
	// rateMap[names.triggerName(i)] += triggerPrescales->getPrescaleForIndex(i)*triggerBits->accept(i);
 //     else
	// rateMap[names.triggerName(i)] = triggerPrescales->getPrescaleForIndex(i)*triggerBits->accept(i);

 //      std::cout << names.triggerName(i) << " " << triggerPrescales->getPrescaleForIndex(i) << " " << triggerBits->accept(i) << std::endl;
										     
 //    }
}


//====================================================
//=== functions to parse json files the manual way ===
//====================================================
std::map<int, std::vector<std::pair<int, int> > >
AODAnalyzer::readJSONFile(const std::string& inFileName)
{
  std::ifstream inFile(inFileName.c_str(), std::ios::in);
  
  std::string line;
  while(!inFile.eof())
    {
      std::string buffer;
      inFile >> buffer;
      line += buffer;
    }
  
  
  
  // define map with result
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
    
  
  
  // loop on JSON file
  for(std::string::const_iterator it = line.begin(); it < line.end(); ++it)
    {
      // find run number
      if( (*(it) == '"') && (*(it+7) == '"') )   
	{
	  std::string run(it+1, it+7);
	  //std::cout << "found run " << run << std::endl;
      
      
      
	  // find lumi sections
	  std::vector<std::pair<int, int> > lumisections;
	  for(std::string::const_iterator it2 = it+10; it2 < line.end(); ++it2)
	    {
	      if( (*(it2) == ']') && (*(it2-1) == ']') ) break;
	      if( *(it2) != '[' ) continue;
        
	      std::string::const_iterator it_beg = it2;
	      std::string::const_iterator it_mid;
	      std::string::const_iterator it_end;
        
	      for(std::string::const_iterator it3 = it_beg; it3 < line.end(); ++it3)
		{
		  if( *(it3) == ',' ) it_mid = it3;
		  if( *(it3) == ']' )
		    {
		      it_end = it3;
		      break;
		    }
		}
            
            
            
	      std::string lumi_beg(it_beg+1, it_mid);
	      std::string lumi_end(it_mid+1, it_end);
	      //std::cout << "[" << lumi_beg;
	      //std::cout << ",";
	      //std::cout << lumi_end << "]" << std::endl;
        
	      std::pair<int, int> tempLS(atoi(lumi_beg.c_str()), atoi(lumi_end.c_str()));
	      lumisections.push_back(tempLS);
        
	      it2 = it_end;
	    }
      
      
	  jsonMap[atoi(run.c_str())] = lumisections;
	} // find run number
    
    } // loop on JSON file
  
  
  
  return jsonMap;
}






bool AODAnalyzer::AcceptEventByRunAndLumiSection(const int& runId, const int& lumiId,
						       std::map<int, std::vector<std::pair<int, int> > >& jsonMap)
{
  // select by runId
  if( jsonMap.find(runId) == jsonMap.end() ) return false;
  
  
  
  // select by lumiId
  std::vector<std::pair<int, int> > lumisections = jsonMap[runId];
  
  int skipEvent = true;
  for(unsigned int i = 0; i < lumisections.size(); ++i)
    if( (lumiId >= lumisections.at(i).first) &&
        (lumiId <= lumisections.at(i).second) )
      skipEvent = false;
  
  if( skipEvent == true ) return false;
  
  
  return true;
}









DEFINE_FWK_MODULE(AODAnalyzer);
