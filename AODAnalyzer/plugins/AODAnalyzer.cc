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


#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" //added because of SuperCluster (collection)
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h" 
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"  
#include "DataFormats/METReco/interface/PFMETFwd.h" 

#include "DataFormats/JetReco/interface/CaloJet.h"
// #include "DataFormats/JetReco/interface/CaloJetfwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
// #include "DataFormats/HcalRecHit/interface/HcalRecHitFwd.h"   --gives error
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"

#include "DataFormats/CastorReco/interface/CastorTower.h"
#include "DataFormats/HcalRecHit/interface/CastorRecHit.h"


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

  template<typename PFChMETCollection>
  void fillPFChMets(const edm::Handle<PFChMETCollection> &);

  template<typename PFMETCollection>
  void fillPFMets(const edm::Handle<PFMETCollection> &);

  template<typename CaloJetCollection> //select variables
  void fillCaloJets(const edm::Handle<CaloJetCollection> &);

  template<typename CaloMETCollection>
  void fillCaloMETs(const edm::Handle<CaloMETCollection> &);

  template<typename SuperClusterCollection> 
  void fillSC(const edm::Handle<SuperClusterCollection> &); 

  template<typename SuperClusterhfEMCollection> 
  void fillSChfEM(const edm::Handle<SuperClusterhfEMCollection> &); 

  template<typename CaloClusterhfEMCollection>
  void fillCChfEM(const edm::Handle<CaloClusterhfEMCollection> &); 

  template<typename CaloClusterCollection>
  void fillCC(const edm::Handle<CaloClusterCollection> &); 

  template<typename PhotonCollection>
  void fillPhotons(const edm::Handle<PhotonCollection> &);

  template<typename PhotongedCollection>
  void fillgedPhotons(const edm::Handle<PhotongedCollection> &);

  template<typename MuonCollection>
  void fillMuons(const edm::Handle<MuonCollection> &);

  template<typename GsfElectronCollection>
  void fillGsf(const edm::Handle<GsfElectronCollection> &);

  template<typename UNGsfElectronCollection>
  void fillUNGsf(const edm::Handle<UNGsfElectronCollection> &);

  template<typename EBEcalRecHitCollection>
  void fillEBrecHit(const edm::Handle<EBEcalRecHitCollection> &);

  template<typename EEEcalRecHitCollection>
  void fillEErecHit(const edm::Handle<EEEcalRecHitCollection> &);
 
  template<typename ESEcalRecHitCollection>
  void fillESrecHit(const edm::Handle<ESEcalRecHitCollection> &);

  template<typename HBHERecHitCollection>
  void fillHBHErecHit(const edm::Handle<HBHERecHitCollection> &);

  template<typename HFRecHitCollection>
  void fillHFrecHit(const edm::Handle<HFRecHitCollection> &);

  template<typename HORecHitCollection>
  void fillHOrecHit(const edm::Handle<HORecHitCollection> &);

  template<typename PreshowerClusterCollection>
  void fillPreshowerCluster(const edm::Handle<PreshowerClusterCollection> &);

  template<typename PreshowerClusterCollectionY>
  void fillPreshowerClusterY(const edm::Handle<PreshowerClusterCollectionY> &);

  template<typename CastorTowerCollection>   //TODO   --find variables for LeafCandidate
  void fillCastorTower(const edm::Handle<CastorTowerCollection> &);

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

  //PFJet variables
  std::vector<float>* PFJetPt_;
  std::vector<float>* PFJetEta_;
  std::vector<float>* PFJetPhi_;
  std::vector<float>* PFChMetPt_;
  std::vector<float>* PFChMetPhi_;
  std::vector<float>* PFMetPt_;
  std::vector<float>* PFMetPhi_;
  std::vector<int>*   nVtx_;
  //CaloJet variables
  std::vector<float>* CalJetPt_;
  std::vector<float>* CalJetEta_;
  std::vector<float>* CalJetPhi_;
  std::vector<float>* CalJetEn_;

  //CaloMet variables
  std::vector<float>* CalMETPt_;
  std::vector<float>* CalMETEta_;
  std::vector<float>* CalMETPhi_;
  std::vector<float>* CalMETEn_;

  //std::vector<float>* SuperCluster_; // adding SuperCluster
  std::vector<double>* SCEn_;
  std::vector<double>* SCEta_;
  std::vector<double>* SCPhi_;
  std::vector<double>* SCEnhfEM_;
  std::vector<double>* SCEtahfEM_;
  std::vector<double>* SCPhihfEM_;

  //caloclusters variables
  std::vector<double>* CCEn_;
  std::vector<double>* CCEta_;
  std::vector<double>* CCPhi_;
  //calocluster hfem variables
  std::vector<double>* CCEnhfEM_;
  std::vector<double>* CCEtahfEM_;
  std::vector<double>* CCPhihfEM_;

  //photon variables
  std::vector<float>* PhoPt_;
  std::vector<float>* PhoEta_;
  std::vector<float>* PhoPhi_;
  std::vector<float>* PhoEn_;

  std::vector<float>* Phoe1x5_;
  std::vector<float>* Phoe2x5_;
  std::vector<float>* Phoe3x3_;
  std::vector<float>* Phoe5x5_;
  std::vector<float>* Phomaxenxtal_;
  std::vector<float>* Phosigmaeta_;
  std::vector<float>* PhosigmaIeta_;
  std::vector<float>* Phor1x5_;
  std::vector<float>* Phor2x5_;
  std::vector<float>* Phor9_;

  std::vector<float>* gedPhoPt_;
  std::vector<float>* gedPhoEta_;
  std::vector<float>* gedPhoPhi_;
  // CURRENTLY FILLING ENERGY() AND NOT CORRECTED ENERGY
  std::vector<float>* gedPhoEn_;
  // CURRENTLY FILLING ENERGY() AND NOT CORRECTED ENERGY

  std::vector<float>* gedPhoe1x5_;
  std::vector<float>* gedPhoe2x5_;
  std::vector<float>* gedPhoe3x3_;
  std::vector<float>* gedPhoe5x5_;
  std::vector<float>* gedPhomaxenxtal_;
  std::vector<float>* gedPhosigmaeta_;
  std::vector<float>* gedPhosigmaIeta_;
  std::vector<float>* gedPhor1x5_;
  std::vector<float>* gedPhor2x5_;
  std::vector<float>* gedPhor9_;

  // Muon variables
  std::vector<float>* MuPt_;
  std::vector<float>* MuEta_;
  std::vector<float>* MuPhi_;
  std::vector<float>* MuEn_;
  std::vector<float>* MuCh_;


  // GSF variables
  std::vector<float>* SigmaIEta_;
  std::vector<float>* SigmaIPhi_;
  std::vector<float>* r9_;
  std::vector<float>* HadOEm_;
  std::vector<float>* drSumPt_;
  std::vector<float>* drSumEt_;
  std::vector<float>* eSCOP_;
  std::vector<float>* ecEn_;

  // UNGSF variables
  std::vector<float>* UNSigmaIEta_;
  std::vector<float>* UNSigmaIPhi_;
  std::vector<float>* UNr9_;
  std::vector<float>* UNHadOEm_;
  std::vector<float>* UNdrSumPt_;
  std::vector<float>* UNdrSumEt_;
  std::vector<float>* UNeSCOP_;
  std::vector<float>* UNecEn_;
  
  //EB rechit, EE rechit and ES rechit variables
  std::vector<float>* EBenergy_;
  std::vector<float>* EBtime_;
  std::vector<float>* EBchi2_;
  std::vector<float>* EEenergy_;
  std::vector<float>* EEtime_;
  std::vector<float>* EEchi2_;
  std::vector<float>* ESenergy_;
  std::vector<float>* EStime_;
  std::vector<float>* ESchi2_;

  //HBHE rechit, HF rechit and HO rechit variables
  std::vector<float>* HBHEenergy_;
  std::vector<float>* HBHEtime_;
  // std::vector<float>* HBHEchi2_;
  std::vector<float>* HFenergy_;
  std::vector<float>* HFtime_;
  // std::vector<float>* HFchi2_;
  std::vector<float>* HOenergy_;
  std::vector<float>* HOtime_;
  // std::vector<float>* HOchi2_;
  
  //preshower variables
  std::vector<double>* PreShEn_;
  std::vector<double>* PreShCorrEn_;
  std::vector<double>* PreShEta_;
  std::vector<double>* PreShPhi_;
  std::vector<double>* PreShYEn_;
  std::vector<double>* PreShYCorrEn_;
  std::vector<double>* PreShYEta_;
  std::vector<double>* PreShYPhi_;

  //castor tower variables
  std::vector<double>* CTPt_;
  std::vector<double>* CTEta_;
  std::vector<double>* CTPhi_;

  std::vector<float>* qPFJetPt_;
  std::vector<float>* qPFJetEta_;
  std::vector<float>* qPFJetPhi_;
  std::vector<float>* qPFChMetPt_;
  std::vector<float>* qPFChMetPhi_;
  std::vector<float>* qPFMetPt_;
  std::vector<float>* qPFMetPhi_;
  std::vector<int>*   qNVtx_;

  std::vector<float>* qCalJetPt_;
  std::vector<float>* qCalJetEta_;
  std::vector<float>* qCalJetPhi_;
  std::vector<float>* qCalJetEn_;
  std::vector<float>* qCalMETPt_;
  std::vector<float>* qCalMETEta_;
  std::vector<float>* qCalMETPhi_;
  std::vector<float>* qCalMETEn_;

  std::vector<double>* qSCEn_;
  std::vector<double>* qSCEta_;
  std::vector<double>* qSCPhi_;
  std::vector<double>* qSCEnhfEM_;
  std::vector<double>* qSCEtahfEM_;
  std::vector<double>* qSCPhihfEM_;
  std::vector<double>* qCCEn_;
  std::vector<double>* qCCEta_;
  std::vector<double>* qCCPhi_;
  std::vector<double>* qCCEnhfEM_;
  std::vector<double>* qCCEtahfEM_;
  std::vector<double>* qCCPhihfEM_;

  std::vector<float>* qPhoPt_;
  std::vector<float>* qPhoEta_;
  std::vector<float>* qPhoPhi_;
  std::vector<float>* qPhoEn_;

  std::vector<float>* qPhoe1x5_;
  std::vector<float>* qPhoe2x5_;
  std::vector<float>* qPhoe3x3_;
  std::vector<float>* qPhoe5x5_;
  std::vector<float>* qPhomaxenxtal_;
  std::vector<float>* qPhosigmaeta_;
  std::vector<float>* qPhosigmaIeta_;
  std::vector<float>* qPhor1x5_;
  std::vector<float>* qPhor2x5_;
  std::vector<float>* qPhor9_;

  std::vector<float>* qgedPhoPt_;
  std::vector<float>* qgedPhoEta_;
  std::vector<float>* qgedPhoPhi_;
  std::vector<float>* qgedPhoEn_;

  std::vector<float>* qgedPhoe1x5_;
  std::vector<float>* qgedPhoe2x5_;
  std::vector<float>* qgedPhoe3x3_;
  std::vector<float>* qgedPhoe5x5_;
  std::vector<float>* qgedPhomaxenxtal_;
  std::vector<float>* qgedPhosigmaeta_;
  std::vector<float>* qgedPhosigmaIeta_;
  std::vector<float>* qgedPhor1x5_;
  std::vector<float>* qgedPhor2x5_;
  std::vector<float>* qgedPhor9_;

  std::vector<float>* qMuPt_;
  std::vector<float>* qMuEta_;
  std::vector<float>* qMuPhi_;
  std::vector<float>* qMuEn_;
  std::vector<float>* qMuCh_;

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

  std::vector<float>* qHBHEenergy_;
  std::vector<float>* qHBHEtime_;
  // std::vector<float>* qHBHEchi2_;
  std::vector<float>* qHFenergy_;
  std::vector<float>* qHFtime_;
  // std::vector<float>* qHFchi2_;
  std::vector<float>* qHOenergy_;
  std::vector<float>* qHOtime_;
  // std::vector<float>* qHOchi2_;

  std::vector<double>* qPreShEn_;
  std::vector<double>* qPreShCorrEn_;
  std::vector<double>* qPreShEta_;
  std::vector<double>* qPreShPhi_;
  std::vector<double>* qPreShYEn_;
  std::vector<double>* qPreShYCorrEn_;
  std::vector<double>* qPreShYEta_;
  std::vector<double>* qPreShYPhi_;

  std::vector<double>* qCTPt_;
  std::vector<double>* qCTEta_;
  std::vector<double>* qCTPhi_;

  std::vector<float>*   crossSection_;
  std::vector<float>*   pathRates_;
  std::vector<std::string>*   pathNames_;
  std::map<std::string,int> rateMap;


  edm::EDGetTokenT<reco::PFJetCollection> PFJetToken_;
  edm::EDGetTokenT<reco::PFMETCollection> PFChMETToken_;
  edm::EDGetTokenT<reco::PFMETCollection> PFMETToken_;

  edm::EDGetTokenT<reco::CaloJetCollection> CaloJetToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> CaloMETToken_;      //variables
  edm::EDGetTokenT<reco::VertexCollection>  vtxToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerPrescales_;  //pat::PackedTriggerPrescales
  edm::EDGetTokenT<reco::SuperClusterCollection>   SuperClusterToken_;  //adding SuperCluster
  edm::EDGetTokenT<reco::SuperClusterCollection>    SuperClusterhfEMToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection>   CaloClusterToken_;  //adding SuperCluster
  edm::EDGetTokenT<reco::CaloClusterCollection>   CaloClusterhfEMToken_;

  edm::EDGetTokenT<reco::PhotonCollection> PhotonToken_;   //TWO types of Photons -- one collection?
  edm::EDGetTokenT<reco::PhotonCollection> gedPhotonToken_;
  edm::EDGetTokenT<reco::MuonCollection>   MuonToken_;

  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronUncleanedToken_;


  //TODO -- finish adding collections!  // edm::   //eCALlASERcORRFILTER like

  // edm:: 

  edm::EDGetTokenT<EcalRecHitCollection> ebRHSrcToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRHSrcToken_;
  edm::EDGetTokenT<EcalRecHitCollection> esRHSrcToken_;

  edm::EDGetTokenT<HBHERecHitCollection> hbheRHcToken_;
  edm::EDGetTokenT<HFRecHitCollection>   hfRHcToken_;
  edm::EDGetTokenT<HORecHitCollection>   hoRHcToken_;

  // edm::EDGetTokenT<SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> ebRHSrcToken_;
  // edm::EDGetTokenT<SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> eeRHSrcToken_;
  edm::EDGetTokenT<reco::PreshowerClusterCollection> preshowerXToken_;
  edm::EDGetTokenT<reco::PreshowerClusterCollection> preshowerYToken_;
  edm::EDGetTokenT<reco::CastorTowerCollection> CastorTowerToken_;   //leafcandidate variables  //perhaps reco




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
  CaloJetToken_             (consumes<reco::CaloJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloJetTag"))),
  CaloMETToken_             (consumes<reco::CaloMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloMETTag"))),  ///std::vector<reco::CaloMET> 
  vtxToken_                 (consumes<reco::VertexCollection>(cfg.getUntrackedParameter<edm::InputTag>("vtx"))),
  triggerBits_              (consumes<edm::TriggerResults>(cfg.getUntrackedParameter<edm::InputTag>("bits"))),
  triggerPrescales_         (consumes<trigger::TriggerEvent>(cfg.getUntrackedParameter<edm::InputTag>("prescales"))),  // pat::PackedTriggerPrescales
  SuperClusterToken_        (consumes<reco::SuperClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("SuperClusterTag"))),
  SuperClusterhfEMToken_    (consumes<reco::SuperClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("SuperClusterhfEMTag"))),
  CaloClusterToken_         (consumes<reco::CaloClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloClusterTag"))),
  CaloClusterhfEMToken_     (consumes<reco::CaloClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloClusterhfEMTag"))),
  PhotonToken_              (consumes<reco::PhotonCollection>(cfg.getUntrackedParameter<edm::InputTag>("PhotonTag"))),
  gedPhotonToken_           (consumes<reco::PhotonCollection>(cfg.getUntrackedParameter<edm::InputTag>("gedPhotonTag"))),
  MuonToken_                (consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonTag"))),
  GsfElectronToken_         (consumes<reco::GsfElectronCollection>(cfg.getUntrackedParameter<edm::InputTag>("GsfElectronTag"))),
  GsfElectronUncleanedToken_(consumes<reco::GsfElectronCollection>(cfg.getUntrackedParameter<edm::InputTag>("GsfElectronUncleanedTag"))),

  
  
  //TODO -- add collections!
  //CollectionEcalRecHitEBTToken_
  ebRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("EBRecHitSourceTag"))),  //ICONFIG -> cfg
  eeRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("EERecHitSourceTag"))),
  esRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("ESRecHitSourceTag"))),

  hbheRHcToken_             (consumes<HBHERecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("HBHERecHitTag"))),  //ICONFIG -> cfg
  hfRHcToken_               (consumes<HFRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("HFRecHitTag"))),
  hoRHcToken_               (consumes<HORecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("HORecHitTag"))),
  preshowerXToken_          (consumes<reco::PreshowerClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("PreshowerClusterXTag"))),
  preshowerYToken_          (consumes<reco::PreshowerClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("PreshowerClusterYTag"))),
  CastorTowerToken_         (consumes<reco::CastorTowerCollection>(cfg.getUntrackedParameter<edm::InputTag>("CastorTowerTag"))),  //perhaps reco  https://github.com/cms-sw/cmssw/blob/09c3fce6626f70fd04223e7dacebf0b485f73f54/DataFormats/CastorReco/interface/CastorTower.h



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

  CalJetPt_->clear();
  CalJetEta_->clear();
  CalJetPhi_->clear();
  CalJetEn_->clear();
  CalMETPt_->clear();
  CalMETEta_->clear();
  CalMETPhi_->clear();
  CalMETEn_->clear();
  //SuperCluster_ ->clear(); //adding SuperCluster
  SCEn_ ->clear();
  SCEta_->clear();
  SCPhi_->clear();
  SCEnhfEM_ ->clear();
  SCEtahfEM_->clear();
  SCPhihfEM_->clear();
  CCEn_ ->clear();
  CCEta_->clear();
  CCPhi_->clear();
  CCEnhfEM_ ->clear();
  CCEtahfEM_->clear();
  CCPhihfEM_->clear();

  PhoPt_->clear();
  PhoEta_->clear();
  PhoPhi_->clear();
  PhoEn_->clear();

  Phoe1x5_->clear();
  Phoe2x5_->clear();
  Phoe3x3_->clear();
  Phoe5x5_->clear();
  Phomaxenxtal_->clear();
  Phosigmaeta_->clear();
  PhosigmaIeta_->clear();
  Phor1x5_->clear();
  Phor2x5_->clear();
  Phor9_->clear();

  gedPhoPt_->clear();
  gedPhoEta_->clear();
  gedPhoPhi_->clear();
  gedPhoEn_->clear();

  gedPhoe1x5_->clear();
  gedPhoe2x5_->clear();
  gedPhoe3x3_->clear();
  gedPhoe5x5_->clear();
  gedPhomaxenxtal_->clear();
  gedPhosigmaeta_->clear();
  gedPhosigmaIeta_->clear();
  gedPhor1x5_->clear();
  gedPhor2x5_->clear();
  gedPhor9_->clear();

  MuPt_->clear();
  MuEta_->clear();
  MuPhi_->clear();
  MuEn_->clear();
  MuCh_->clear();

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

  HBHEenergy_->clear();
  HBHEtime_->clear();
  // HBHEchi2_->clear();
  HFenergy_->clear();
  HFtime_->clear();
  // HFchi2_->clear();
  HOenergy_->clear();
  HOtime_->clear();
  // HOchi2_->clear();
 
  PreShEn_->clear();
  PreShCorrEn_->clear();
  PreShEta_->clear();
  PreShPhi_->clear();
  PreShYEn_->clear();
  PreShYCorrEn_->clear();
  PreShYEta_->clear();
  PreShYPhi_->clear();

  CTPt_->clear();
  CTEta_->clear();
  CTPhi_->clear();

  qPFJetPt_->clear();
  qPFJetEta_->clear();
  qPFJetPhi_->clear();
  qPFChMetPt_->clear();
  qPFChMetPhi_->clear();
  qPFMetPt_->clear();
  qPFMetPhi_->clear();
  qNVtx_->clear();

  qCalJetPt_->clear();
  qCalJetEta_->clear();
  qCalJetPhi_->clear();
  qCalJetEn_->clear();
  qCalMETPt_->clear();
  qCalMETEta_->clear();
  qCalMETPhi_->clear();
  qCalMETEn_->clear();
  //qSuperCluster_ ->clear(); //ASK IF quantiles should be here
  qSCEn_ ->clear();
  qSCEta_->clear();
  qSCPhi_->clear();
  qSCEnhfEM_ ->clear();
  qSCEtahfEM_->clear();
  qSCPhihfEM_->clear();
  qCCEn_ ->clear();
  qCCEta_->clear();
  qCCPhi_->clear();
  qCCEnhfEM_ ->clear();
  qCCEtahfEM_->clear();
  qCCPhihfEM_->clear();

  qPhoPt_->clear();
  qPhoEta_->clear();
  qPhoPhi_->clear();
  qPhoEn_->clear();

  qPhoe1x5_->clear();
  qPhoe2x5_->clear();
  qPhoe3x3_->clear();
  qPhoe5x5_->clear();
  qPhomaxenxtal_->clear();
  qPhosigmaeta_->clear();
  qPhosigmaIeta_->clear();
  qPhor1x5_->clear();
  qPhor2x5_->clear();
  qPhor9_->clear();

  qgedPhoPt_->clear();
  qgedPhoEta_->clear();
  qgedPhoPhi_->clear();
  qgedPhoEn_->clear();

  qgedPhoe1x5_->clear();
  qgedPhoe2x5_->clear();
  qgedPhoe3x3_->clear();
  qgedPhoe5x5_->clear();
  qgedPhomaxenxtal_->clear();
  qgedPhosigmaeta_->clear();
  qgedPhosigmaIeta_->clear();
  qgedPhor1x5_->clear();
  qgedPhor2x5_->clear();
  qgedPhor9_->clear();

  qMuPt_->clear();
  qMuEta_->clear();
  qMuPhi_->clear();
  qMuEn_->clear();
  qMuCh_->clear();
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

  qHBHEenergy_->clear();
  qHBHEtime_->clear();
  // qHBHEchi2_->clear();
  qHFenergy_->clear();
  qHFtime_->clear();
  // qHFchi2_->clear();
  qHOenergy_->clear();
  qHOtime_->clear();
  // qHOchi2_->clear();

  qPreShEn_->clear();
  qPreShCorrEn_->clear();
  qPreShEta_->clear();
  qPreShPhi_->clear();
  qPreShYEn_->clear();
  qPreShYCorrEn_->clear();
  qPreShYEta_->clear();
  qPreShYPhi_->clear();

  qCTPt_->clear();
  qCTEta_->clear();
  qCTPhi_->clear();

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
  typename jetCollection::const_iterator i = jets->begin();
  for(;i != jets->end(); i++){
    if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
      {
	
	if(type.compare(std::string("PF")) == 0)
	  {
	    PFJetPt_->push_back(i->pt());
	    PFJetEta_->push_back(i->eta());
	    PFJetPhi_->push_back(i->phi());
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
	  }
	
      }
  }
  return;
}

template<typename PFChMETCollection>
void AODAnalyzer::fillPFChMets(const edm::Handle<PFChMETCollection> & pfchmets)
{
  // std::cout << "fillPFChMets is being called!" << std::endl;
  // std::cout << pfchmets->size() <<std::endl;
  typename PFChMETCollection::const_iterator i = pfchmets->begin();
  for(;i != pfchmets->end(); i++){
    PFChMetPt_->push_back(i->et());
    PFChMetPhi_->push_back(i->phi());

  }
  return;

}

template<typename PFMETCollection>
void AODAnalyzer::fillPFMets(const edm::Handle<PFMETCollection> & pfmets)
{
  // std::cout << "fis is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename PFMETCollection::const_iterator i = pfmets->begin();
  for(;i != pfmets->end(); i++){
    PFMetPt_->push_back(i->et());
    PFMetPhi_->push_back(i->phi());  //try also eta and energy

  }
  return;

}

template<typename CaloJetCollection>
void AODAnalyzer::fillCaloJets(const edm::Handle<CaloJetCollection> & calojets)
{
  // std::cout << "f is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename CaloJetCollection::const_iterator i = calojets->begin();
  for(;i != calojets->end(); i++){
        CalJetPt_->push_back(i->et());
        CalJetEta_->push_back(i->eta());
        CalJetPhi_->push_back(i->phi());
        CalJetEn_->push_back(i->energy());

  }
  return;

}

template<typename CaloMETCollection>
void AODAnalyzer::fillCaloMETs(const edm::Handle<CaloMETCollection> & caloMETs)
{
  // std::cout << "fills is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename CaloMETCollection::const_iterator i = caloMETs->begin();
  for(;i != caloMETs->end(); i++){
        CalMETPt_->push_back(i->et());
        CalMETEta_->push_back(i->eta());
        CalMETPhi_->push_back(i->phi());
        CalMETEn_->push_back(i->energy());

  }
  return;

}
// //TODO
// template<typename CaloJetCollection>
// void AODAnalyzer::fillCaloJets()

// //TODO
// template<typename CaloMETCollection>
// void AODAnalyzer::fillCaloMETs()


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

        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename SuperClusterhfEMCollection>
void AODAnalyzer::fillSChfEM(const edm::Handle<SuperClusterhfEMCollection> & superclustershfEM) //ask for jets analogy //SUPERCLUSTERS
{

  // Selected jets
  //reco::CaloJetCollection recojets;
  typename SuperClusterhfEMCollection::const_iterator i = superclustershfEM->begin();
  for(;i != superclustershfEM->end(); i++){
     if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // not sure if needed
      // {
      SCEnhfEM_->push_back(i->energy());
      SCEtahfEM_->push_back(i->etaWidth());
      SCPhihfEM_->push_back(i->phiWidth());

        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename CaloClusterCollection>
void AODAnalyzer::fillCC(const edm::Handle<CaloClusterCollection> & caloclusters) //ask for jets analogy //SUPERCLUSTERS
{

  // Selected jets
  //reco::CaloJetCollection recojets;
  typename CaloClusterCollection::const_iterator i = caloclusters->begin();
  for(;i != caloclusters->end(); i++){
     //if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_)  //do I need something like maxCCeta and so on? TODO
      // {
      CCEn_->push_back(i->energy());
      CCEta_->push_back(i->eta());
      CCPhi_->push_back(i->phi());

        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename CaloClusterhfEMCollection>
void AODAnalyzer::fillCChfEM(const edm::Handle<CaloClusterhfEMCollection> & caloclustershfEM) //ask for jets analogy //SUPERCLUSTERS
{

  // Selected jets
  //reco::CaloJetCollection recojets;
  typename CaloClusterhfEMCollection::const_iterator i = caloclustershfEM->begin();
  for(;i != caloclustershfEM->end(); i++){
     // if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // do I need something like maxCCeta and so on? TODOd
      // {
      CCEnhfEM_->push_back(i->energy());
      CCEtahfEM_->push_back(i->eta());
      CCPhihfEM_->push_back(i->phi());

        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename PhotonCollection>
void AODAnalyzer::fillPhotons(const edm::Handle<PhotonCollection> & photons)
{
   typename PhotonCollection::const_iterator i = photons->begin();
   for(;i != photons->end(); i++){
     
        PhoPt_->push_back(i->et());
        PhoEta_->push_back(i->eta());
        PhoPhi_->push_back(i->phi());
        PhoEn_->push_back(i->energy());   //GETCORRECTEDENERGY!!
        Phoe1x5_->push_back(i->e1x5());   
        Phoe2x5_->push_back(i->e2x5());
        Phoe3x3_->push_back(i->e3x3());
        Phoe5x5_->push_back(i->e5x5());
        Phomaxenxtal_->push_back(i->maxEnergyXtal());
        Phosigmaeta_->push_back(i->sigmaEtaEta());
        PhosigmaIeta_->push_back(i->sigmaIetaIeta());
        Phor1x5_->push_back(i->r1x5());
        Phor2x5_->push_back(i->r2x5());
        Phor9_->push_back(i->r9());
        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename PhotongedCollection>
void AODAnalyzer::fillgedPhotons(const edm::Handle<PhotongedCollection> & gedphotons)
{
   typename PhotongedCollection::const_iterator i = gedphotons->begin();
   for(;i != gedphotons->end(); i++){
     
        gedPhoPt_->push_back(i->et());
        gedPhoEta_->push_back(i->eta());
        gedPhoPhi_->push_back(i->phi());
        gedPhoEn_->push_back(i->energy());
        gedPhoe1x5_->push_back(i->e1x5());   
        gedPhoe2x5_->push_back(i->e2x5());
        gedPhoe3x3_->push_back(i->e3x3());
        gedPhoe5x5_->push_back(i->e5x5());
        gedPhomaxenxtal_->push_back(i->maxEnergyXtal());
        gedPhosigmaeta_->push_back(i->sigmaEtaEta());
        gedPhosigmaIeta_->push_back(i->sigmaIetaIeta());
        gedPhor1x5_->push_back(i->r1x5());
        gedPhor2x5_->push_back(i->r2x5());
        gedPhor9_->push_back(i->r9());
        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename MuonCollection>
void AODAnalyzer::fillMuons(const edm::Handle<MuonCollection> & muons)
{
   typename MuonCollection::const_iterator i = muons->begin();
   for(;i != muons->end(); i++){
     
        MuPt_->push_back(i->et());
        MuEta_->push_back(i->eta());
        MuPhi_->push_back(i->phi());
        MuEn_->push_back(i->energy());   //GETCORRECTEDENERGY!!
        MuCh_->push_back(i->charge());
        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}


//TODO
template<typename GsfElectronCollection>
void AODAnalyzer::fillGsf(const edm::Handle<GsfElectronCollection> & electrons)
{
  // std::cout << "fillGSF is being called!" << std::endl;
  // std::cout << electrons->size() <<std::endl;
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

template<typename UNGsfElectronCollection>
void AODAnalyzer::fillUNGsf(const edm::Handle<UNGsfElectronCollection> & UNelectrons)
{
  
  // std::cout << "fillUNGSF is being called!" << std::endl;
  typename UNGsfElectronCollection::const_iterator i = UNelectrons->begin();
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

template<typename EBEcalRecHitCollection>
void AODAnalyzer::fillEBrecHit(const edm::Handle<EBEcalRecHitCollection> & EBhits)
{

  // std::cout << "fillEBrecHit is being called!" << std::endl;
  typename EBEcalRecHitCollection::const_iterator i = EBhits->begin();
  for(;i != EBhits->end(); i++){
    EBenergy_ ->push_back(i->energy());
    EBtime_ ->push_back(i->time());
    EBchi2_ ->push_back(i->chi2());

  }
  return;
}

template<typename EEEcalRecHitCollection>
void AODAnalyzer::fillEErecHit(const edm::Handle<EEEcalRecHitCollection> & EEhits)
{

  // std::cout << "fillEErecHit is being called!" << std::endl;
  typename EEEcalRecHitCollection::const_iterator i = EEhits->begin();
  for(;i != EEhits->end(); i++){
    EEenergy_ ->push_back(i->energy());
    EEtime_ ->push_back(i->time());
    EEchi2_ ->push_back(i->chi2());

  }
  return;
}

template<typename ESEcalRecHitCollection>
void AODAnalyzer::fillESrecHit(const edm::Handle<ESEcalRecHitCollection> & EShits)
{

  // std::cout << "fillESrecHit is being called!" << std::endl;
  typename ESEcalRecHitCollection::const_iterator i = EShits->begin();
  for(;i != EShits->end(); i++){
    ESenergy_ ->push_back(i->energy());
    EStime_ ->push_back(i->time());
    ESchi2_ ->push_back(i->chi2()); 

  }
  return;
}
///---------------------------------------- HBHE
template<typename HBHERecHitCollection>
void AODAnalyzer::fillHBHErecHit(const edm::Handle<HBHERecHitCollection> & HBHEhits)
{

  // std::cout << "fillHBHErecHit is being called!" << std::endl;
  typename HBHERecHitCollection::const_iterator i = HBHEhits->begin();
  for(;i != HBHEhits->end(); i++){
    HBHEenergy_ ->push_back(i->energy());
    HBHEtime_ ->push_back(i->time());
    // HBHEchi2_ ->push_back(i->chi2()); const class HBHERecHit' has no member named 'chi2'
    // std::cout << "ele HBHEenergy: " << i->energy()   << std::endl; 
    // std::cout << "ele HBHEtime: "  << i->time() << std::endl;
  }
  return;
}

template<typename HFRecHitCollection>
void AODAnalyzer::fillHFrecHit(const edm::Handle<HFRecHitCollection> & HFhits)
{

  // std::cout << "fillHFrecHit is being called!" << std::endl;
  typename HFRecHitCollection::const_iterator i = HFhits->begin();
  for(;i != HFhits->end(); i++){
    HFenergy_ ->push_back(i->energy());
    HFtime_ ->push_back(i->time());
    // HFchi2_ ->push_back(i->chi2());
    // std::cout << "ele HFenergy: " << i->energy()   << std::endl; 
    // std::cout << "ele HFtime: "  << i->time() << std::endl;
  }
  return;
}

template<typename HORecHitCollection>
void AODAnalyzer::fillHOrecHit(const edm::Handle<HORecHitCollection> & HOhits)
{

  // std::cout << "fillHOrecHit is being called!" << std::endl;
  typename HORecHitCollection::const_iterator i = HOhits->begin();
  for(;i != HOhits->end(); i++){
    HOenergy_ ->push_back(i->energy());
    HOtime_ ->push_back(i->time());
    // HOchi2_ ->push_back(i->chi2());
    std::cout << "ele HOenergy: " << i->energy()   << std::endl; 
    std::cout << "ele HOtime: "  << i->time() << std::endl;
  }
  return;
}

template<typename PreshowerClusterCollection>
void AODAnalyzer::fillPreshowerCluster(const edm::Handle<PreshowerClusterCollection> & preshowerclusterhits)
{

  // std::cout << "fillHOrecHit is being called!" << std::endl;
  typename PreshowerClusterCollection::const_iterator i = preshowerclusterhits->begin();
  for(;i != preshowerclusterhits->end(); i++){

    PreShEn_->push_back(i->energy());
    PreShCorrEn_->push_back(i->correctedEnergy());
    PreShEta_->push_back(i->eta());
    PreShPhi_->push_back(i->phi());

  }
  return;
}

template<typename PreshowerClusterCollectionY>
void AODAnalyzer::fillPreshowerClusterY(const edm::Handle<PreshowerClusterCollectionY> & preshowerclusterYhits)
{

  // std::cout << "fillHOrecHit is being called!" << std::endl;
  typename PreshowerClusterCollectionY::const_iterator i = preshowerclusterYhits->begin();
  for(;i != preshowerclusterYhits->end(); i++){

    PreShYEn_->push_back(i->energy());
    PreShYCorrEn_->push_back(i->correctedEnergy());
    PreShYEta_->push_back(i->eta());
    PreShYPhi_->push_back(i->phi());

  }
  return;
}


template<typename CastorTowerCollection>
void AODAnalyzer::fillCastorTower(const edm::Handle<CastorTowerCollection> & castors)
{

  // std::cout << "fillHOrecHit is being called!" << std::endl;
  typename CastorTowerCollection::const_iterator i = castors->begin();
  for(;i != castors->end(); i++){

    CTPt_->push_back(i->pt());
    CTEta_->push_back(i->eta());
    CTPhi_->push_back(i->phi());

  }
  return;
}

// ----------------------------------- HO

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
  
  CalJetPt_   = new std::vector<float>;
  CalJetEta_   = new std::vector<float>;
  CalJetPhi_   = new std::vector<float>;
  CalJetEn_   = new std::vector<float>;
  CalMETPt_   = new std::vector<float>;
  CalMETEta_   = new std::vector<float>;
  CalMETPhi_   = new std::vector<float>;
  CalMETEn_   = new std::vector<float>;
  //SuperCluster_ = new std::vector<float>; //adding SuperCluster
  SCEn_      = new std::vector<double>;
  SCEta_     = new std::vector<double>;
  SCPhi_     = new std::vector<double>;
  SCEnhfEM_      = new std::vector<double>;
  SCEtahfEM_     = new std::vector<double>;
  SCPhihfEM_     = new std::vector<double>;
  CCEn_      = new std::vector<double>;
  CCEta_     = new std::vector<double>;
  CCPhi_     = new std::vector<double>;
  CCEnhfEM_      = new std::vector<double>;
  CCEtahfEM_     = new std::vector<double>;
  CCPhihfEM_     = new std::vector<double>;

  PhoPt_     = new std::vector<float>;
  PhoEta_    = new std::vector<float>;
  PhoPhi_    = new std::vector<float>;
  PhoEn_ = new std::vector<float>;

  Phoe1x5_   = new std::vector<float>;
  Phoe2x5_   = new std::vector<float>;
  Phoe3x3_   = new std::vector<float>;
  Phoe5x5_   = new std::vector<float>;
  Phomaxenxtal_ = new std::vector<float>;
  Phosigmaeta_  = new std::vector<float>;
  PhosigmaIeta_ = new std::vector<float>;
  Phor1x5_   = new std::vector<float>;
  Phor2x5_   = new std::vector<float>;
  Phor9_     = new std::vector<float>;

  gedPhoPt_     = new std::vector<float>;
  gedPhoEta_    = new std::vector<float>;
  gedPhoPhi_    = new std::vector<float>;
  gedPhoEn_ = new std::vector<float>;

  gedPhoe1x5_   = new std::vector<float>;
  gedPhoe2x5_   = new std::vector<float>;
  gedPhoe3x3_   = new std::vector<float>;
  gedPhoe5x5_   = new std::vector<float>;
  gedPhomaxenxtal_ = new std::vector<float>;
  gedPhosigmaeta_  = new std::vector<float>;
  gedPhosigmaIeta_ = new std::vector<float>;
  gedPhor1x5_   = new std::vector<float>;
  gedPhor2x5_   = new std::vector<float>;
  gedPhor9_     = new std::vector<float>;

  MuPt_         = new std::vector<float>;
  MuEta_        = new std::vector<float>;
  MuPhi_        = new std::vector<float>;
  MuEn_         = new std::vector<float>;
  MuCh_         = new std::vector<float>;
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

  HBHEenergy_  = new std::vector<float>;
  HBHEtime_    = new std::vector<float>;
  // HBHEchi2_    = new std::vector<float>;
  HFenergy_    = new std::vector<float>;
  HFtime_      = new std::vector<float>;
  // HFchi2_      = new std::vector<float>;
  HOenergy_    = new std::vector<float>;
  HOtime_      = new std::vector<float>;
  // HOchi2_      = new std::vector<float>;


  PreShEn_    = new std::vector<double>;
  PreShCorrEn_= new std::vector<double>;
  PreShEta_   = new std::vector<double>;
  PreShPhi_   = new std::vector<double>;

  PreShYEn_    = new std::vector<double>;
  PreShYCorrEn_= new std::vector<double>;
  PreShYEta_   = new std::vector<double>;
  PreShYPhi_   = new std::vector<double>;

  CTPt_        = new std::vector<double>;
  CTEta_       = new std::vector<double>;
  CTPhi_       = new std::vector<double>;
  
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
 
  qCalJetPt_   = new std::vector<float>;
  qCalJetEta_   = new std::vector<float>;
  qCalJetPhi_   = new std::vector<float>;
  qCalJetEn_   = new std::vector<float>;
  qCalMETPt_   = new std::vector<float>;
  qCalMETEta_   = new std::vector<float>;
  qCalMETPhi_   = new std::vector<float>;
  qCalMETEn_   = new std::vector<float>;

  qSCEn_     = new std::vector<double>;
  qSCEta_    = new std::vector<double>;
  qSCPhi_    = new std::vector<double>;
  qSCEnhfEM_     = new std::vector<double>;
  qSCEtahfEM_    = new std::vector<double>;
  qSCPhihfEM_    = new std::vector<double>;
  qCCEn_     = new std::vector<double>;
  qCCEta_    = new std::vector<double>;
  qCCPhi_    = new std::vector<double>;
  qCCEnhfEM_     = new std::vector<double>;
  qCCEtahfEM_    = new std::vector<double>;
  qCCPhihfEM_    = new std::vector<double>;

  qPhoPt_     = new std::vector<float>;
  qPhoEta_    = new std::vector<float>;
  qPhoPhi_    = new std::vector<float>;
  qPhoEn_ = new std::vector<float>;

  qPhoe1x5_   = new std::vector<float>;
  qPhoe2x5_   = new std::vector<float>;
  qPhoe3x3_   = new std::vector<float>;
  qPhoe5x5_   = new std::vector<float>;
  qPhomaxenxtal_ = new std::vector<float>;
  qPhosigmaeta_  = new std::vector<float>;
  qPhosigmaIeta_ = new std::vector<float>;
  qPhor1x5_   = new std::vector<float>;
  qPhor2x5_   = new std::vector<float>;
  qPhor9_     = new std::vector<float>;

  qgedPhoPt_     = new std::vector<float>;
  qgedPhoEta_    = new std::vector<float>;
  qgedPhoPhi_    = new std::vector<float>;
  qgedPhoEn_ = new std::vector<float>;

  qgedPhoe1x5_   = new std::vector<float>;
  qgedPhoe2x5_   = new std::vector<float>;
  qgedPhoe3x3_   = new std::vector<float>;
  qgedPhoe5x5_   = new std::vector<float>;
  qgedPhomaxenxtal_ = new std::vector<float>;
  qgedPhosigmaeta_  = new std::vector<float>;
  qgedPhosigmaIeta_ = new std::vector<float>;
  qgedPhor1x5_   = new std::vector<float>;
  qgedPhor2x5_   = new std::vector<float>;
  qgedPhor9_     = new std::vector<float>;

  qMuPt_   = new std::vector<float>;
  qMuEta_   = new std::vector<float>;
  qMuPhi_   = new std::vector<float>;
  qMuEn_   = new std::vector<float>;
  qMuCh_   = new std::vector<float>;

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

  qHBHEenergy_  = new std::vector<float>;
  qHBHEtime_    = new std::vector<float>;
  // qHBHEchi2_    = new std::vector<float>;
  qHFenergy_    = new std::vector<float>;
  qHFtime_      = new std::vector<float>;
  // qHFchi2_      = new std::vector<float>;
  qHOenergy_    = new std::vector<float>;
  qHOtime_      = new std::vector<float>;
  // qHOchi2_      = new std::vector<float>;
//FINISH doing HBHE<HF<HO

  qPreShEn_     = new std::vector<double>;
  qPreShCorrEn_ = new std::vector<double>;
  qPreShEta_    = new std::vector<double>;
  qPreShPhi_    = new std::vector<double>;
  qPreShYEn_    = new std::vector<double>;
  qPreShYCorrEn_= new std::vector<double>;
  qPreShYEta_   = new std::vector<double>;
  qPreShYPhi_   = new std::vector<double>;

  qCTPt_        = new std::vector<double>;
  qCTEta_       = new std::vector<double>;
  qCTPhi_       = new std::vector<double>;

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

  outTree_->Branch("qCalJetPt",     "std::vector<std::float>",        &qCalJetPt_);
  outTree_->Branch("qCalJetEta",    "std::vector<std::float>",        &qCalJetEta_);
  outTree_->Branch("qCalJetPhi",    "std::vector<std::float>",        &qCalJetPhi_);
  outTree_->Branch("qCalJetEn",    "std::vector<std::float>",        &qCalJetEn_);

  outTree_->Branch("qCalMETPt",     "std::vector<std::float>",        &qCalMETPt_);
  outTree_->Branch("qCalMETEta",    "std::vector<std::float>",        &qCalMETEta_);
  outTree_->Branch("qCalMETPhi",    "std::vector<std::float>",        &qCalMETPhi_);
  outTree_->Branch("qCalMETEn",    "std::vector<std::float>",        &qCalMETEn_);
  //outTree_->Branch("SuperCluster","std::vector<std::float>",      &SuperCluster_);
  outTree_->Branch("qSCEn",     "std::vector<std::double>",        &qSCEn_);
  outTree_->Branch("qSCEta",    "std::vector<std::double>",        &qSCEta_);
  outTree_->Branch("qSCPhi",    "std::vector<std::double>",        &qSCPhi_);
  outTree_->Branch("qSCEnhfEM",     "std::vector<std::double>",        &qSCEnhfEM_);
  outTree_->Branch("qSCEtahfEM",    "std::vector<std::double>",        &qSCEtahfEM_);
  outTree_->Branch("qSCPhihfEM",    "std::vector<std::double>",        &qSCPhihfEM_);
  outTree_->Branch("qCCEn",     "std::vector<std::double>",        &qCCEn_);
  outTree_->Branch("qCCEta",    "std::vector<std::double>",        &qCCEta_);
  outTree_->Branch("qCCPhi",    "std::vector<std::double>",        &qCCPhi_);
  outTree_->Branch("qCCEnhfEM",     "std::vector<std::double>",        &qCCEnhfEM_);
  outTree_->Branch("qCCEtahfEM",    "std::vector<std::double>",        &qCCEtahfEM_);
  outTree_->Branch("qCCPhihfEM",    "std::vector<std::double>",        &qCCPhihfEM_);

  outTree_->Branch("qPhoPt",     "std::vector<std::float>",        &qPhoPt_);
  outTree_->Branch("qPhoEta",    "std::vector<std::float>",        &qPhoEta_);
  outTree_->Branch("qPhoPhi",    "std::vector<std::float>",        &qPhoPhi_);
  outTree_->Branch("qPhoEn_",    "std::vector<std::float>",    &qPhoEn_);

  outTree_->Branch("qPhoe1x5_",     "std::vector<std::float>",        &qPhoe1x5_);
  outTree_->Branch("qPhoe2x5_",    "std::vector<std::float>",        &qPhoe2x5_);
  outTree_->Branch("qPhoe3x3_",    "std::vector<std::float>",        &qPhoe3x3_);
  outTree_->Branch("qPhoe5x5_",    "std::vector<std::float>",    &qPhoe5x5_);
  outTree_->Branch("qPhomaxenxtal_",     "std::vector<std::float>",        &qPhomaxenxtal_);
  outTree_->Branch("qPhosigmaeta_",    "std::vector<std::float>",        &qPhosigmaeta_);
  outTree_->Branch("qPhosigmaIeta_",    "std::vector<std::float>",        &qPhosigmaIeta_);
  outTree_->Branch("qPhor1x5_",    "std::vector<std::float>",    &qPhor1x5_);
  outTree_->Branch("qPhor2x5_",    "std::vector<std::float>",        &qPhor2x5_);
  outTree_->Branch("qPhor9_",    "std::vector<std::float>",    &qPhor9_);

  outTree_->Branch("qgedPhoPt",     "std::vector<std::float>",     &qgedPhoPt_);
  outTree_->Branch("qgedPhoEta",    "std::vector<std::float>",     &qgedPhoEta_);
  outTree_->Branch("qgedPhoPhi",    "std::vector<std::float>",     &qgedPhoPhi_);
  outTree_->Branch("qgedPhoEn_",    "std::vector<std::float>", &qgedPhoEn_);

  outTree_->Branch("qgedPhoe1x5_",     "std::vector<std::float>",        &qgedPhoe1x5_);
  outTree_->Branch("qgedPhoe2x5_",    "std::vector<std::float>",        &qgedPhoe2x5_);
  outTree_->Branch("qgedPhoe3x3_",    "std::vector<std::float>",        &qgedPhoe3x3_);
  outTree_->Branch("qgedPhoe5x5_",    "std::vector<std::float>",    &qgedPhoe5x5_);
  outTree_->Branch("qgedPhomaxenxtal_",     "std::vector<std::float>",        &qgedPhomaxenxtal_);
  outTree_->Branch("qgedPhosigmaeta_",    "std::vector<std::float>",        &qgedPhosigmaeta_);
  outTree_->Branch("qgedPhosigmaIeta_",    "std::vector<std::float>",        &qgedPhosigmaIeta_);
  outTree_->Branch("qgedPhor1x5_",    "std::vector<std::float>",    &qgedPhor1x5_);
  outTree_->Branch("qgedPhor2x5_",    "std::vector<std::float>",        &qgedPhor2x5_);
  outTree_->Branch("qgedPhor9_",    "std::vector<std::float>",    &qgedPhor9_);

  outTree_->Branch("qMuPt",     "std::vector<std::float>",        &qMuPt_);
  outTree_->Branch("qMuEta",    "std::vector<std::float>",        &qMuEta_);
  outTree_->Branch("qMuPhi",    "std::vector<std::float>",        &qMuPhi_);
  outTree_->Branch("qMuEn_",    "std::vector<std::float>",        &qMuEn_);
  outTree_->Branch("qMuCh_",    "std::vector<std::float>",        &qMuCh_);

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

  outTree_->Branch("qHBHEenergy",    "std::vector<std::float>",        &qHBHEenergy_);
  outTree_->Branch("qHBHEtime",    "std::vector<std::float>",          &qHBHEtime_);
  // outTree_->Branch("qHBHEchi2",    "std::vector<std::float>",          &qHBHEchi2_);
  outTree_->Branch("qHFenergy",    "std::vector<std::float>",        &qHFenergy_);
  outTree_->Branch("qHFtime",    "std::vector<std::float>",          &qHFtime_);
  // outTree_->Branch("qHFchi2",    "std::vector<std::float>",          &qHFchi2_);
  outTree_->Branch("qHOenergy",    "std::vector<std::float>",        &qHOenergy_);
  outTree_->Branch("qHOtime",    "std::vector<std::float>",          &qHOtime_);
  // outTree_->Branch("qHOchi2",    "std::vector<std::float>",          &qHOchi2_);

  outTree_->Branch("qPreShEn",     "std::vector<std::double>",        &qPreShEn_);
  outTree_->Branch("qPreShCorrEn", "std::vector<std::double>",        &qPreShCorrEn_);
  outTree_->Branch("qPreShEta",    "std::vector<std::double>",        &qPreShEta_);
  outTree_->Branch("qPreShPhi",    "std::vector<std::double>",        &qPreShPhi_);
  outTree_->Branch("qPreShYEn",     "std::vector<std::double>",        &qPreShYEn_);
  outTree_->Branch("qPreShYCorrEn", "std::vector<std::double>",        &qPreShYCorrEn_);
  outTree_->Branch("qPreShYEta",    "std::vector<std::double>",        &qPreShYEta_);
  outTree_->Branch("qPreShYPhi",    "std::vector<std::double>",        &qPreShYPhi_);

  outTree_->Branch("qCTPt",     "std::vector<std::double>",        &qCTPt_);
  outTree_->Branch("qCTEta",    "std::vector<std::double>",        &qCTEta_);
  outTree_->Branch("qCTPhi",    "std::vector<std::double>",        &qCTPhi_);

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

  delete CalJetPt_;
  delete CalJetEta_;
  delete CalJetPhi_;
  delete CalJetEn_;
  delete CalMETPt_;
  delete CalMETEta_;
  delete CalMETPhi_;
  delete CalMETEn_;

  delete SCEn_;
  delete SCEta_;
  delete SCPhi_;
  delete SCEnhfEM_;
  delete SCEtahfEM_;
  delete SCPhihfEM_;
  delete CCEn_;
  delete CCEta_;
  delete CCPhi_;
  delete CCEnhfEM_;
  delete CCEtahfEM_;
  delete CCPhihfEM_;

  delete PhoPt_;
  delete PhoEta_;
  delete PhoPhi_;
  delete PhoEn_;

  delete Phoe1x5_;
  delete Phoe2x5_;
  delete Phoe3x3_;
  delete Phoe5x5_;
  delete Phomaxenxtal_;
  delete Phosigmaeta_;
  delete PhosigmaIeta_;
  delete Phor1x5_;
  delete Phor2x5_;
  delete Phor9_;

  delete gedPhoPt_;
  delete gedPhoEta_;
  delete gedPhoPhi_;
  delete gedPhoEn_;

  delete gedPhoe1x5_;
  delete gedPhoe2x5_;
  delete gedPhoe3x3_;
  delete gedPhoe5x5_;
  delete gedPhomaxenxtal_;
  delete gedPhosigmaeta_;
  delete gedPhosigmaIeta_;
  delete gedPhor1x5_;
  delete gedPhor2x5_;
  delete gedPhor9_;

  delete MuPt_;
  delete MuEta_;
  delete MuPhi_;
  delete MuEn_;
  delete MuCh_;
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

  delete HBHEenergy_ ;
  delete HBHEtime_ ;
  //delete HBHEchi2_ ;
  delete HFenergy_ ;
  delete HFtime_ ;
  //delete HFchi2_  ;
  delete HOenergy_ ;
  delete HOtime_ ;
  // HOchi2_    
  delete PreShEn_;
  delete PreShCorrEn_;
  delete PreShEta_;
  delete PreShPhi_;
  delete PreShYEn_;
  delete PreShYCorrEn_;
  delete PreShYEta_;
  delete PreShYPhi_;
  delete CTPt_;
  delete CTEta_;
  delete CTPhi_;
  //delete SuperCluster_;
  //delete qSuperCluster_;
  delete qPFJetPt_;
  delete qPFJetEta_;
  delete qPFJetPhi_;
  delete qPFChMetPt_;
  delete qPFChMetPhi_;
  delete qPFMetPt_;
  delete qPFMetPhi_;

  delete qCalJetPt_;
  delete qCalJetEta_;
  delete qCalJetPhi_;
  delete qCalJetEn_;
  delete qCalMETPt_;
  delete qCalMETEta_;
  delete qCalMETPhi_;
  delete qCalMETEn_;

  delete qSCEn_;
  delete qSCEta_;
  delete qSCPhi_;
  delete qSCEnhfEM_;
  delete qSCEtahfEM_;
  delete qSCPhihfEM_;
  delete qCCEn_;
  delete qCCEta_;
  delete qCCPhi_;
  delete qCCEnhfEM_;
  delete qCCEtahfEM_;
  delete qCCPhihfEM_;

  delete qPhoPt_;
  delete qPhoEta_;
  delete qPhoPhi_;
  delete qPhoEn_;

  delete qPhoe1x5_;
  delete qPhoe2x5_;
  delete qPhoe3x3_;
  delete qPhoe5x5_;
  delete qPhomaxenxtal_;
  delete qPhosigmaeta_;
  delete qPhosigmaIeta_;
  delete qPhor1x5_;
  delete qPhor2x5_;
  delete qPhor9_;

  delete qgedPhoPt_;
  delete qgedPhoEta_;
  delete qgedPhoPhi_;
  delete qgedPhoEn_;

  delete qgedPhoe1x5_;
  delete qgedPhoe2x5_;
  delete qgedPhoe3x3_;
  delete qgedPhoe5x5_;
  delete qgedPhomaxenxtal_;
  delete qgedPhosigmaeta_;
  delete qgedPhosigmaIeta_;
  delete qgedPhor1x5_;
  delete qgedPhor2x5_;
  delete qgedPhor9_;

  delete qMuPt_;
  delete qMuEta_;
  delete qMuPhi_;
  delete qMuEn_;
  delete qMuCh_;
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

  delete qHBHEenergy_;
  delete qHBHEtime_;
  // delete qHBHEchi2_;
  delete qHFenergy_;
  delete qHFtime_;
  // delete qHFchi2_;
  delete qHOenergy_;
  delete qHOtime_;
  // delete qHOchi2_;

  delete qPreShEn_;
  delete qPreShCorrEn_;
  delete qPreShEta_;
  delete qPreShPhi_;
  delete qPreShYEn_;
  delete qPreShYCorrEn_;
  delete qPreShYEta_;
  delete qPreShYPhi_;

  delete qCTPt_;
  delete qCTEta_;
  delete qCTPhi_;

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

  computeMeanAndRms(CalJetPt_, qCalJetPt_);
  computeMeanAndRms(CalJetEta_,qCalJetEta_);
  computeMeanAndRms(CalJetPhi_,qCalJetPhi_);
  computeMeanAndRms(CalJetEn_,qCalJetEn_);

  computeMeanAndRms(CalMETPt_, qCalMETPt_);
  computeMeanAndRms(CalMETEta_,qCalMETEta_);
  computeMeanAndRms(CalMETPhi_,qCalMETPhi_);
  computeMeanAndRms(CalMETEn_,qCalMETEn_);

  computeMeanAndRms(SCEn_, qSCEn_);   
  computeMeanAndRms(SCEta_, qSCEta_);  
  computeMeanAndRms(SCPhi_, qSCPhi_);
  computeMeanAndRms(SCEnhfEM_, qSCEnhfEM_);   
  computeMeanAndRms(SCEtahfEM_, qSCEtahfEM_);  
  computeMeanAndRms(SCPhihfEM_, qSCPhihfEM_); 
  computeMeanAndRms(CCEn_, qCCEn_);   
  computeMeanAndRms(CCEta_, qCCEta_);  
  computeMeanAndRms(CCPhi_, qCCPhi_);
  computeMeanAndRms(CCEnhfEM_, qCCEnhfEM_);   
  computeMeanAndRms(CCEtahfEM_, qCCEtahfEM_);  
  computeMeanAndRms(CCPhihfEM_, qCCPhihfEM_); 


  computeMeanAndRms(PhoPt_, qPhoPt_);
  computeMeanAndRms(PhoEta_,qPhoEta_);
  computeMeanAndRms(PhoPhi_,qPhoPhi_);
  computeMeanAndRms(PhoEn_,qPhoEn_);

  computeMeanAndRms(Phoe1x5_, qPhoe1x5_);
  computeMeanAndRms(Phoe2x5_,qPhoe2x5_);
  computeMeanAndRms(Phoe3x3_,qPhoe3x3_);
  computeMeanAndRms(Phoe5x5_,qPhoe5x5_);
  computeMeanAndRms(Phomaxenxtal_, qPhomaxenxtal_);
  computeMeanAndRms(Phosigmaeta_,  qPhosigmaeta_);
  computeMeanAndRms(PhosigmaIeta_, qPhosigmaIeta_);
  computeMeanAndRms(Phor1x5_, qPhor1x5_);
  computeMeanAndRms(Phor2x5_, qPhor2x5_);
  computeMeanAndRms(Phor9_,   qPhor9_);


  computeMeanAndRms(gedPhoPt_, qgedPhoPt_);
  computeMeanAndRms(gedPhoEta_,qgedPhoEta_);
  computeMeanAndRms(gedPhoPhi_,qgedPhoPhi_);
  computeMeanAndRms(gedPhoEn_,qgedPhoEn_);

  computeMeanAndRms(gedPhoe1x5_, qgedPhoe1x5_);
  computeMeanAndRms(gedPhoe2x5_,qgedPhoe2x5_);
  computeMeanAndRms(gedPhoe3x3_,qgedPhoe3x3_);
  computeMeanAndRms(gedPhoe5x5_,qgedPhoe5x5_);
  computeMeanAndRms(gedPhomaxenxtal_, qgedPhomaxenxtal_);
  computeMeanAndRms(gedPhosigmaeta_,  qgedPhosigmaeta_);
  computeMeanAndRms(gedPhosigmaIeta_, qgedPhosigmaIeta_);
  computeMeanAndRms(gedPhor1x5_, qgedPhor1x5_);
  computeMeanAndRms(gedPhor2x5_, qgedPhor2x5_);
  computeMeanAndRms(gedPhor9_,   qgedPhor9_);

  computeMeanAndRms(MuPt_, qMuPt_);
  computeMeanAndRms(MuEta_,qMuEta_);
  computeMeanAndRms(MuPhi_,qMuPhi_);
  computeMeanAndRms(MuEn_,qMuEn_);
  computeMeanAndRms(MuCh_,qMuCh_);
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
  computeMeanAndRms(HBHEenergy_, qHBHEenergy_);
  computeMeanAndRms(HBHEtime_, qHBHEtime_);
  // computeMeanAndRms(HBHEchi2_, qHBHEchi2_);
  computeMeanAndRms(HFenergy_, qHFenergy_);
  computeMeanAndRms(HFtime_, qHFtime_);
  // computeMeanAndRms(HFchi2_, qHFchi2_);
  computeMeanAndRms(HOenergy_, qHOenergy_);
  computeMeanAndRms(HOtime_, qHOtime_);
  // computeMeanAndRms(HOchi2_, qHOchi2_);

  computeMeanAndRms(PreShEn_, qPreShEn_);   
  computeMeanAndRms(PreShCorrEn_, qPreShCorrEn_);   
  computeMeanAndRms(PreShEta_, qPreShEta_);  //
  computeMeanAndRms(PreShPhi_, qPreShPhi_);  //
  computeMeanAndRms(PreShYEn_, qPreShYEn_);   
  computeMeanAndRms(PreShYCorrEn_, qPreShYCorrEn_);   
  computeMeanAndRms(PreShYEta_, qPreShYEta_);  //
  computeMeanAndRms(PreShYPhi_, qPreShYPhi_);  //

  computeMeanAndRms(CTPt_, qCTPt_);   
  computeMeanAndRms(CTEta_, qCTEta_);  
  computeMeanAndRms(CTPhi_, qCTPhi_);

  computeQuantiles(PFJetPt_, qPFJetPt_, quantiles_);
  computeQuantiles(PFJetEta_,qPFJetEta_,quantiles_);
  computeQuantiles(PFJetPhi_,qPFJetPhi_,quantiles_);
  computeQuantiles(PFChMetPt_, qPFChMetPt_,     quantiles_);
  computeQuantiles(PFChMetPhi_,qPFChMetPhi_,    quantiles_);
  computeQuantiles(PFMetPt_, qPFMetPt_,     quantiles_);
  computeQuantiles(PFMetPhi_,qPFMetPhi_,    quantiles_);
  computeQuantiles(nVtx_,    qNVtx_,    quantiles_);

  computeQuantiles(CalJetPt_, qCalJetPt_, quantiles_);
  computeQuantiles(CalJetEta_,qCalJetEta_,quantiles_);
  computeQuantiles(CalJetPhi_,qCalJetPhi_,quantiles_);
  computeQuantiles(CalJetEn_,qCalJetEn_,  quantiles_);
  computeQuantiles(CalMETPt_, qCalMETPt_, quantiles_);
  computeQuantiles(CalMETEta_,qCalMETEta_,quantiles_);
  computeQuantiles(CalMETPhi_,qCalMETPhi_,quantiles_);
  computeQuantiles(CalMETEn_,qCalMETEn_,  quantiles_);

  computeQuantiles(SCEn_, qSCEn_,       quantiles_);
  computeQuantiles(SCEta_, qSCEta_,     quantiles_);
  computeQuantiles(SCPhi_, qSCPhi_,     quantiles_);
  computeQuantiles(SCEnhfEM_, qSCEnhfEM_,       quantiles_);
  computeQuantiles(SCEtahfEM_, qSCEtahfEM_,     quantiles_);
  computeQuantiles(SCPhihfEM_, qSCPhihfEM_,     quantiles_);
  computeQuantiles(CCEn_, qCCEn_,       quantiles_);
  computeQuantiles(CCEta_, qCCEta_,     quantiles_);
  computeQuantiles(CCPhi_, qCCPhi_,     quantiles_);
  computeQuantiles(CCEnhfEM_, qCCEnhfEM_,       quantiles_);
  computeQuantiles(CCEtahfEM_, qCCEtahfEM_,     quantiles_);
  computeQuantiles(CCPhihfEM_, qCCPhihfEM_,     quantiles_);

  computeQuantiles(PhoPt_, qPhoPt_, quantiles_);
  computeQuantiles(PhoEta_,qPhoEta_,quantiles_);
  computeQuantiles(PhoPhi_,qPhoPhi_,quantiles_);
  computeQuantiles(PhoEn_,qPhoEn_,quantiles_);

  computeQuantiles(Phoe1x5_, qPhoe1x5_, quantiles_);
  computeQuantiles(Phoe2x5_,qPhoe2x5_,quantiles_);
  computeQuantiles(Phoe3x3_,qPhoe3x3_,quantiles_);
  computeQuantiles(Phoe5x5_,qPhoe5x5_,quantiles_);
  computeQuantiles(Phomaxenxtal_, qPhomaxenxtal_, quantiles_);
  computeQuantiles(Phosigmaeta_,qPhosigmaeta_,quantiles_);
  computeQuantiles(PhosigmaIeta_,qPhosigmaIeta_,quantiles_);
  computeQuantiles(Phor1x5_,qPhor1x5_,quantiles_);
  computeQuantiles(Phor2x5_, qPhor2x5_, quantiles_);
  computeQuantiles(Phor9_,qPhor9_,quantiles_);

  computeQuantiles(gedPhoPt_, qgedPhoPt_, quantiles_);
  computeQuantiles(gedPhoEta_,qgedPhoEta_,quantiles_);
  computeQuantiles(gedPhoPhi_,qgedPhoPhi_,quantiles_);
  computeQuantiles(gedPhoEn_,qgedPhoEn_,quantiles_);

  computeQuantiles(gedPhoe1x5_, qgedPhoe1x5_, quantiles_);
  computeQuantiles(gedPhoe2x5_,qgedPhoe2x5_,quantiles_);
  computeQuantiles(gedPhoe3x3_,qgedPhoe3x3_,quantiles_);
  computeQuantiles(gedPhoe5x5_,qgedPhoe5x5_,quantiles_);
  computeQuantiles(gedPhomaxenxtal_, qgedPhomaxenxtal_, quantiles_);
  computeQuantiles(gedPhosigmaeta_,qgedPhosigmaeta_,quantiles_);
  computeQuantiles(gedPhosigmaIeta_,qgedPhosigmaIeta_,quantiles_);
  computeQuantiles(gedPhor1x5_,qgedPhor1x5_,quantiles_);
  computeQuantiles(gedPhor2x5_, qgedPhor2x5_, quantiles_);
  computeQuantiles(gedPhor9_,qgedPhor9_,quantiles_);

  computeQuantiles(MuPt_, qMuPt_, quantiles_);
  computeQuantiles(MuEta_,qMuEta_,quantiles_);
  computeQuantiles(MuPhi_,qMuPhi_,quantiles_);
  computeQuantiles(MuEn_,qMuEn_,quantiles_);
  computeQuantiles(MuCh_,qMuCh_,quantiles_);
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

  computeQuantiles(HBHEenergy_, qHBHEenergy_, quantiles_);
  computeQuantiles(HBHEtime_, qHBHEtime_, quantiles_);
  // computeQuantiles(HBHEchi2_, qHBHEchi2_, quantiles_);
  computeQuantiles(HFenergy_, qHFenergy_, quantiles_);
  computeQuantiles(HFtime_, qHFtime_, quantiles_);
  // computeQuantiles(HFchi2_, qHFchi2_, quantiles_);
  computeQuantiles(HOenergy_, qHOenergy_, quantiles_);
  computeQuantiles(HOtime_, qHOtime_, quantiles_);
  // computeQuantiles(HOchi2_, qHOchi2_, quantiles_);


  computeQuantiles(PreShEn_, qPreShEn_,       quantiles_);
  computeQuantiles(PreShCorrEn_, qPreShCorrEn_,   quantiles_);
  computeQuantiles(PreShEta_, qPreShEta_,     quantiles_);
  computeQuantiles(PreShPhi_, qPreShPhi_,     quantiles_);
  computeQuantiles(PreShYEn_, qPreShYEn_,       quantiles_);
  computeQuantiles(PreShYCorrEn_, qPreShYCorrEn_,   quantiles_);
  computeQuantiles(PreShYEta_, qPreShYEta_,     quantiles_);
  computeQuantiles(PreShYPhi_, qPreShYPhi_,     quantiles_);

  computeQuantiles(CTPt_, qCTPt_, quantiles_);
  computeQuantiles(CTEta_,qCTEta_,quantiles_);
  computeQuantiles(CTPhi_,qCTPhi_,quantiles_);

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
  if(pfmetlocalv.isValid())
    fillPFMets(pfmetlocalv);
  // edm::Handle<std::vector<reco::PFMETCollection> > PFMet;    ////
  // event.getByToken(PFMETToken_,PFMet);
  // if(PFMet.isValid())
  //   {
  //     PFMetPt_->push_back( (*PFMet)[0].et() );
  //     PFMetPhi_->push_back( (*PFMet)[0].phi() );
  //   } 

  edm::Handle<reco::CaloJetCollection> calojetlocalv;
  event.getByToken(CaloJetToken_, calojetlocalv);
  if(calojetlocalv.isValid())
    fillCaloJets(calojetlocalv);

  edm::Handle<reco::CaloMETCollection> caloMETlocalv;
  event.getByToken(CaloMETToken_, caloMETlocalv);
  if(caloMETlocalv.isValid())
    fillCaloMETs(caloMETlocalv);

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

  edm::Handle<reco::SuperClusterCollection> SuperClusterhfEMlocalv;
  event.getByToken(SuperClusterhfEMToken_, SuperClusterhfEMlocalv);
  // print the size of SuperClusterlocalv
  if(SuperClusterhfEMlocalv.isValid())
    fillSChfEM(SuperClusterhfEMlocalv);

    //Fill CaloCluster
  edm::Handle<reco::CaloClusterCollection> CaloClusterlocalv;
  event.getByToken(CaloClusterToken_, CaloClusterlocalv);
  // print the size of SuperClusterlocalv
  if(CaloClusterlocalv.isValid())
    fillCC(CaloClusterlocalv);

  edm::Handle<reco::CaloClusterCollection> CaloClusterhfEMlocalv;
  event.getByToken(CaloClusterhfEMToken_, CaloClusterhfEMlocalv);
  // print the size of SuperClusterlocalv
  if(CaloClusterhfEMlocalv.isValid())
    fillCChfEM(CaloClusterhfEMlocalv);


  edm::Handle<reco::PhotonCollection> photonlocalv;
  event.getByToken(PhotonToken_, photonlocalv);
  if(photonlocalv.isValid())
    fillPhotons(photonlocalv);

  edm::Handle<reco::PhotonCollection> gedphotonlocalv;
  event.getByToken(gedPhotonToken_, gedphotonlocalv);
  if(gedphotonlocalv.isValid())
    fillgedPhotons(gedphotonlocalv);

  edm::Handle<reco::MuonCollection> muonlocalv;
  event.getByToken(MuonToken_, muonlocalv);
  if(muonlocalv.isValid())
    fillMuons(muonlocalv);

  //TODO --fill fill Muons
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

  //fill hbhereco  
  edm::Handle<HBHERecHitCollection> hbheRHs;    //finish fill ee, eb and so on
  event.getByToken(hbheRHcToken_, hbheRHs);
  if(hbheRHs.isValid())
    fillHBHErecHit(hbheRHs);

  //fill hfreco
  edm::Handle<HFRecHitCollection> hfRHs;
  event.getByToken(hfRHcToken_, hfRHs);
  if(hfRHs.isValid())
    fillHFrecHit(hfRHs);

  //fill horeco
  edm::Handle<HORecHitCollection> hoRHs;
  event.getByToken(hoRHcToken_, hoRHs);
  if(hoRHs.isValid())
    fillHOrecHit(hoRHs);

  edm::Handle<reco::PreshowerClusterCollection> prShs;
  event.getByToken(preshowerXToken_, prShs);
  if(prShs.isValid())
    fillPreshowerCluster(prShs);

  edm::Handle<reco::PreshowerClusterCollection> prShYs;
  event.getByToken(preshowerYToken_, prShYs);
  if(prShYs.isValid())
    fillPreshowerClusterY(prShYs);

  edm::Handle<reco::CastorTowerCollection> castorsz;
  event.getByToken(CastorTowerToken_, castorsz);
  if(castorsz.isValid())
    fillCastorTower(castorsz);
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
