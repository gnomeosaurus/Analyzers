#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

//new includes for ECAL and HCAL mapping


#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

//#include "DQMServices/Core/interface/MonitorElement.h" gives error NoSocket.h

//end

//PTComparator for sorting by pt
#include "CommonTools/Utils/interface/PtComparator.h"
#include "Validation/EventGenerator/interface/HepMCValidationHelper.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"


#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/HLTGlobalStatus.h"

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
// #include "DataFormats/HcalRecHit/interface/HORecHit.h"
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
  virtual void beginRun(const edm::Run & run,    const edm::EventSetup & eventSetup);
  virtual void endRun  (const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  virtual void beginLuminosityBlock  (const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  virtual void endLuminosityBlock    (const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  
private:
  template<typename jetCollection>
  void fillJets(const edm::Handle<jetCollection> &, std::string );

  template<typename PFJet4CHSCollection>
  void fill4CHSJets(const edm::Handle<PFJet4CHSCollection> &, std::string );

  // //---------------------------Sorted variables start here.

  // template<typename SortedjetCollection>
  // void fillSortedJets(const edm::Handle<SortedjetCollection> &, std::string );

  // template<typename PFJet4CHSSortedCollection>
  // void fill4CHSSortedJets(const edm::Handle<PFJet4CHSSortedCollection> &, std::string );  

  // template<typename PFJet8CHSSortedCollection>
  // void fill8CHSSortedJets(const edm::Handle<PFJet8CHSSortedCollection> &, std::string );

  // template<typename PFJetEISortedCollection>
  // void fillEISortedJets(const edm::Handle<PFJetEISortedCollection> &, std::string );

  // template<typename PFJet8CHSSoftDropSortedCollection>
  // void fill8CHSoftDropSortedJets(const edm::Handle<PFJet8CHSSoftDropSortedCollection> &, std::string );

  // template<typename PFJetTopCHSSortedCollection>
  // void fillTopCHSSortedJets(const edm::Handle<PFJetTopCHSSortedCollection> &, std::string );

  // template<typename CaloJetSortedCollection> 
  // void fillCaloSortedJets(const edm::Handle<CaloJetSortedCollection> &);

  // template<typename PhotonSortedCollection>
  // void fillSortedPhotons(const edm::Handle<PhotonSortedCollection> &);

  // template<typename PhotongedSortedCollection>
  // void fillSortedgedPhotons(const edm::Handle<PhotongedSortedCollection> &);

  // template<typename MuonSortedCollection>
  // void fillSortedMuons(const edm::Handle<MuonSortedCollection> &);

  // template<typename MuonCosmSortedCollection>
  // void fillCosmSortedMuons(const edm::Handle<MuonCosmSortedCollection> &);

  // template<typename MuonCosmLegSortedCollection>
  // void fillCosmLegSortedMuons(const edm::Handle<MuonCosmLegSortedCollection> &);


  //___________________________Sorted variables end here.

   template<typename PFJet8CHSCollection>
  void fill8CHSJets(const edm::Handle<PFJet8CHSCollection> &, std::string );

  template<typename PFJetEICollection>
  void fillEIJets(const edm::Handle<PFJetEICollection> &, std::string );

  template<typename PFJet8CHSSoftDropCollection>
  void fill8CHSoftDropJets(const edm::Handle<PFJet8CHSSoftDropCollection> &, std::string );

  template<typename PFJetTopCHSCollection>
  void fillTopCHSJets(const edm::Handle<PFJetTopCHSCollection> &, std::string );

  template<typename PFChMETCollection>
  void fillPFChMets(const edm::Handle<PFChMETCollection> &);

  template<typename PFMETCollection>
  void fillPFMets(const edm::Handle<PFMETCollection> &);

  template<typename CaloJetCollection> 
  void fillCaloJets(const edm::Handle<CaloJetCollection> &);

  template<typename CaloMETCollection>
  void fillCaloMETs(const edm::Handle<CaloMETCollection> &);

  template<typename CaloMETBECollection>
  void fillCaloMETBEs(const edm::Handle<CaloMETBECollection> &);

  template<typename CaloMETBEFOCollection>
  void fillCaloMETBEFOs(const edm::Handle<CaloMETBEFOCollection> &);

  template<typename CaloMETMCollection>
  void fillCaloMETMs(const edm::Handle<CaloMETMCollection> &);

  template<typename SuperClusterCollection> 
  void fillSC(const edm::Handle<SuperClusterCollection> &); 

  template<typename SuperClusterhfEMCollection> 
  void fillSChfEM(const edm::Handle<SuperClusterhfEMCollection> &); 

  template<typename SuperCluster5x5Collection> 
  void fillSC5x5(const edm::Handle<SuperCluster5x5Collection> &); 

  template<typename CaloCluster5x5Collection>
  void fillCC5x5(const edm::Handle<CaloCluster5x5Collection> &); 

  template<typename CaloClusterCollection>
  void fillCC(const edm::Handle<CaloClusterCollection> &); 

  template<typename PhotonCollection>
  void fillPhotons(const edm::Handle<PhotonCollection> &);

  template<typename PhotongedCollection>
  void fillgedPhotons(const edm::Handle<PhotongedCollection> &);

  template<typename MuonCollection>
  void fillMuons(const edm::Handle<MuonCollection> &);

  template<typename MuonCosmCollection>
  void fillCosmMuons(const edm::Handle<MuonCosmCollection> &);

  template<typename MuonCosmLegCollection>
  void fillCosmLegMuons(const edm::Handle<MuonCosmLegCollection> &);

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

  // template<typename HORecHitCollection>
  // void fillHOrecHit(const edm::Handle<HORecHitCollection> &);

  template<typename PreshowerClusterCollection>
  void fillPreshowerCluster(const edm::Handle<PreshowerClusterCollection> &);

  template<typename PreshowerClusterCollectionY>
  void fillPreshowerClusterY(const edm::Handle<PreshowerClusterCollectionY> &);

  // template<typename CastorTowerCollection>   //TODO   doesnt output any values for eta,en,phi
  // void fillCastorTower(const edm::Handle<CastorTowerCollection> &);

  void initialize();
  template<typename T>
  void computeQuantiles(std::vector<T>*, std::vector<T>*, std::vector<double>); 
  template<typename T>
  void computeMeanAndRms(std::vector<T>*, std::vector<T>*);
  std::map<int, std::vector<std::pair<int, int> > > readJSONFile(const std::string&);
  bool AcceptEventByRunAndLumiSection(const int& , const int& ,std::map<int, std::vector<std::pair<int, int> > >&);



  /// file service and tree
  edm::Service<TFileService> outfile_;
  HLTConfigProvider hltConfigProvider_;

  TTree* outTree_;
  int    runId_;
  int    lumiId_;
  float  lumi_;
  int    isSig_;


  //PFJet variables
  std::vector<float>* PFJetPt_;
  std::vector<float>* PFJetEta_;
  std::vector<float>* PFJetPhi_;
 
  //PFJet sorted variables
  std::vector<float>* PFJet0Pt_; 
  std::vector<float>* PFJet1Pt_; 
  std::vector<float>* PFJet2Pt_; 
  std::vector<float>* PFJet3Pt_; 
  std::vector<float>* PFJet4Pt_; 
  std::vector<float>* PFJet5Pt_; 

  std::vector<float>* PFJet0Eta_;
  std::vector<float>* PFJet1Eta_;
  std::vector<float>* PFJet2Eta_;
  std::vector<float>* PFJet3Eta_;
  std::vector<float>* PFJet4Eta_;
  std::vector<float>* PFJet5Eta_;

  std::vector<float>* PFJet0Phi_;
  std::vector<float>* PFJet1Phi_;
  std::vector<float>* PFJet2Phi_;
  std::vector<float>* PFJet3Phi_;
  std::vector<float>* PFJet4Phi_;
  std::vector<float>* PFJet5Phi_;


  //PF4CHS variables
  std::vector<float>* PFJet4CHSPt_;
  std::vector<float>* PFJet4CHSEta_;
  std::vector<float>* PFJet4CHSPhi_;


  //PFJet4CHS sorted variables
  std::vector<float>* PFJet4CHS0Pt_;
  std::vector<float>* PFJet4CHS1Pt_;
  std::vector<float>* PFJet4CHS2Pt_;
  std::vector<float>* PFJet4CHS3Pt_;
  std::vector<float>* PFJet4CHS4Pt_;
  std::vector<float>* PFJet4CHS5Pt_;

  std::vector<float>* PFJet4CHS0Eta_;
  std::vector<float>* PFJet4CHS1Eta_;
  std::vector<float>* PFJet4CHS2Eta_;
  std::vector<float>* PFJet4CHS3Eta_;
  std::vector<float>* PFJet4CHS4Eta_;
  std::vector<float>* PFJet4CHS5Eta_;

  std::vector<float>* PFJet4CHS0Phi_;
  std::vector<float>* PFJet4CHS1Phi_;
  std::vector<float>* PFJet4CHS2Phi_;
  std::vector<float>* PFJet4CHS3Phi_;
  std::vector<float>* PFJet4CHS4Phi_;
  std::vector<float>* PFJet4CHS5Phi_;


  //PF8CHS variables
  std::vector<float>* PFJet8CHSPt_;
  std::vector<float>* PFJet8CHSEta_;
  std::vector<float>* PFJet8CHSPhi_;


  //PF8CHS sorted variables
  std::vector<float>* PFJet8CHS0Pt_;
  std::vector<float>* PFJet8CHS1Pt_;
  std::vector<float>* PFJet8CHS2Pt_;
  std::vector<float>* PFJet8CHS3Pt_;
  std::vector<float>* PFJet8CHS4Pt_;
  std::vector<float>* PFJet8CHS5Pt_;

  std::vector<float>* PFJet8CHS0Eta_;
  std::vector<float>* PFJet8CHS1Eta_;
  std::vector<float>* PFJet8CHS2Eta_;
  std::vector<float>* PFJet8CHS3Eta_;
  std::vector<float>* PFJet8CHS4Eta_;
  std::vector<float>* PFJet8CHS5Eta_;

  std::vector<float>* PFJet8CHS0Phi_;
  std::vector<float>* PFJet8CHS1Phi_;
  std::vector<float>* PFJet8CHS2Phi_;
  std::vector<float>* PFJet8CHS3Phi_;
  std::vector<float>* PFJet8CHS4Phi_;
  std::vector<float>* PFJet8CHS5Phi_;


  //PFJetEI variables
  std::vector<float>* PFJetEIPt_;
  std::vector<float>* PFJetEIEta_;
  std::vector<float>* PFJetEIPhi_;


  //PFJetEI sorted variables
  std::vector<float>* PFJetEI0Pt_;
  std::vector<float>* PFJetEI1Pt_;
  std::vector<float>* PFJetEI2Pt_;
  std::vector<float>* PFJetEI3Pt_;
  std::vector<float>* PFJetEI4Pt_;
  std::vector<float>* PFJetEI5Pt_;

  std::vector<float>* PFJetEI0Eta_;
  std::vector<float>* PFJetEI1Eta_;
  std::vector<float>* PFJetEI2Eta_;
  std::vector<float>* PFJetEI3Eta_;
  std::vector<float>* PFJetEI4Eta_;
  std::vector<float>* PFJetEI5Eta_;

  std::vector<float>* PFJetEI0Phi_;
  std::vector<float>* PFJetEI1Phi_;
  std::vector<float>* PFJetEI2Phi_;
  std::vector<float>* PFJetEI3Phi_;
  std::vector<float>* PFJetEI4Phi_;
  std::vector<float>* PFJetEI5Phi_;


  //8CHSSoftDrop variables
  std::vector<float>* PFJet8CHSSDPt_;
  std::vector<float>* PFJet8CHSSDEta_;
  std::vector<float>* PFJet8CHSSDPhi_;


  // //8CHSSoftDrop sorted variables
  // std::vector<float>* PFJet8CHSSD0Pt_;
  // std::vector<float>* PFJet8CHSSD1Pt_;
  // std::vector<float>* PFJet8CHSSD2Pt_;
  // std::vector<float>* PFJet8CHSSD3Pt_;
  // std::vector<float>* PFJet8CHSSD4Pt_;
  // std::vector<float>* PFJet8CHSSD5Pt_;

  // std::vector<float>* PFJet8CHSSD0Eta_;
  // std::vector<float>* PFJet8CHSSD1Eta_;
  // std::vector<float>* PFJet8CHSSD2Eta_;
  // std::vector<float>* PFJet8CHSSD3Eta_;
  // std::vector<float>* PFJet8CHSSD4Eta_;
  // std::vector<float>* PFJet8CHSSD5Eta_;

  // std::vector<float>* PFJet8CHSSD0Phi_;
  // std::vector<float>* PFJet8CHSSD1Phi_;
  // std::vector<float>* PFJet8CHSSD2Phi_;
  // std::vector<float>* PFJet8CHSSD3Phi_;
  // std::vector<float>* PFJet8CHSSD4Phi_;
  // std::vector<float>* PFJet8CHSSD5Phi_;


  //TopCHS variables
  std::vector<float>* PFJetTopCHSPt_;
  std::vector<float>* PFJetTopCHSEta_;
  std::vector<float>* PFJetTopCHSPhi_;


  // //TopCHS sorted variables
  // std::vector<float>* PFJetTopCHS0Pt_;
  // std::vector<float>* PFJetTopCHS1Pt_;
  // std::vector<float>* PFJetTopCHS2Pt_;
  // std::vector<float>* PFJetTopCHS3Pt_;
  // std::vector<float>* PFJetTopCHS4Pt_;
  // std::vector<float>* PFJetTopCHS5Pt_;

  // std::vector<float>* PFJetTopCHS0Eta_;
  // std::vector<float>* PFJetTopCHS1Eta_;
  // std::vector<float>* PFJetTopCHS2Eta_;
  // std::vector<float>* PFJetTopCHS3Eta_;
  // std::vector<float>* PFJetTopCHS4Eta_;
  // std::vector<float>* PFJetTopCHS5Eta_;

  // std::vector<float>* PFJetTopCHS0Phi_;
  // std::vector<float>* PFJetTopCHS1Phi_;
  // std::vector<float>* PFJetTopCHS2Phi_;
  // std::vector<float>* PFJetTopCHS3Phi_;
  // std::vector<float>* PFJetTopCHS4Phi_;
  // std::vector<float>* PFJetTopCHS5Phi_;


  //PFChMet variables
  std::vector<float>* PFChMetPt_;
  // std::vector<float>* PFChMetEta_;  //MET doesn't have ETA
  std::vector<float>* PFChMetPhi_;
  std::vector<float>* PFMetPt_;
  // std::vector<float>* PFMetEta_;  //MET doesn't have ETA
  std::vector<float>* PFMetPhi_;
  std::vector<int>*   nVtx_;

  //CaloJet variables
  std::vector<float>* CalJetPt_;
  std::vector<float>* CalJetEta_;
  std::vector<float>* CalJetPhi_;
  std::vector<float>* CalJetEn_;


  //CaloJet sorted variables
  std::vector<float>* CalJet0Pt_;
  std::vector<float>* CalJet1Pt_;
  std::vector<float>* CalJet2Pt_;
  std::vector<float>* CalJet3Pt_;
  std::vector<float>* CalJet4Pt_;
  std::vector<float>* CalJet5Pt_;

  std::vector<float>* CalJet0Eta_;
  std::vector<float>* CalJet1Eta_;
  std::vector<float>* CalJet2Eta_;
  std::vector<float>* CalJet3Eta_;
  std::vector<float>* CalJet4Eta_;
  std::vector<float>* CalJet5Eta_;

  std::vector<float>* CalJet0Phi_;
  std::vector<float>* CalJet1Phi_;
  std::vector<float>* CalJet2Phi_;
  std::vector<float>* CalJet3Phi_;
  std::vector<float>* CalJet4Phi_;
  std::vector<float>* CalJet5Phi_;

  std::vector<float>* CalJet0En_;
  std::vector<float>* CalJet1En_;
  std::vector<float>* CalJet2En_;
  std::vector<float>* CalJet3En_;
  std::vector<float>* CalJet4En_;
  std::vector<float>* CalJet5En_;

  //CaloMet variables
  std::vector<float>* CalMETPt_;
  // std::vector<float>* CalMETEta_;  //MET doesn't have ETA
  std::vector<float>* CalMETPhi_;
  std::vector<float>* CalMETEn_;

  //CaloMet variables BE
  std::vector<float>* CalMETBEPt_;
  // std::vector<float>* CalMETBEEta_;   //MET doesn't have ETA
  std::vector<float>* CalMETBEPhi_;
  std::vector<float>* CalMETBEEn_;

    //CaloMet variables  BEFO
  std::vector<float>* CalMETBEFOPt_;
  // std::vector<float>* CalMETBEFOEta_;  //MET doesn't have ETA
  std::vector<float>* CalMETBEFOPhi_;
  std::vector<float>* CalMETBEFOEn_;

    //CaloMet variables  M
  std::vector<float>* CalMETMPt_;
  // std::vector<float>* CalMETMEta_;   //MET doesn't have ETA
  std::vector<float>* CalMETMPhi_;
  std::vector<float>* CalMETMEn_;

  // SuperCluster variables
  std::vector<float>* SCEn_;
  std::vector<float>* SCEta_;
  std::vector<float>* SCPhi_;
  std::vector<float>* SCEtaWidth_;
  std::vector<float>* SCPhiWidth_;
  std::vector<float>* SCEnhfEM_;
  std::vector<float>* SCEtahfEM_;
  std::vector<float>* SCPhihfEM_;
  // std::vector<float>* SCEtaWidthhfEM_;  //ouputted 0
  // std::vector<float>* SCPhiWidthhfEM_;  //outputted 0
  std::vector<float>* SCEn5x5_;
  std::vector<float>* SCEta5x5_;
  std::vector<float>* SCPhi5x5_;
  std::vector<float>* SCEtaWidth5x5_;
  std::vector<float>* SCPhiWidth5x5_; 

  //caloclusters variables
  std::vector<float>* CCEn_;
  std::vector<float>* CCEta_;
  std::vector<float>* CCPhi_;

  //calocluster hfem variables
  std::vector<float>* CCEn5x5_;
  std::vector<float>* CCEta5x5_;
  std::vector<float>* CCPhi5x5_;

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


  //photon sorted variables (not all of them)
  std::vector<float>* Pho0Pt_;
  std::vector<float>* Pho1Pt_;
  std::vector<float>* Pho2Pt_;
  std::vector<float>* Pho3Pt_;
  std::vector<float>* Pho4Pt_;
  std::vector<float>* Pho5Pt_;

  std::vector<float>* Pho0Eta_;
  std::vector<float>* Pho1Eta_;
  std::vector<float>* Pho2Eta_;
  std::vector<float>* Pho3Eta_;
  std::vector<float>* Pho4Eta_;
  std::vector<float>* Pho5Eta_;

  std::vector<float>* Pho0Phi_;
  std::vector<float>* Pho1Phi_;
  std::vector<float>* Pho2Phi_;
  std::vector<float>* Pho3Phi_;
  std::vector<float>* Pho4Phi_;
  std::vector<float>* Pho5Phi_;

  std::vector<float>* Pho0En_;
  std::vector<float>* Pho1En_;
  std::vector<float>* Pho2En_;
  std::vector<float>* Pho3En_;
  std::vector<float>* Pho4En_;
  std::vector<float>* Pho5En_;


  //ged photon variables
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


  //ged photons sorted variables (not all of them)
  std::vector<float>* gedPho0Pt_;
  std::vector<float>* gedPho1Pt_;
  std::vector<float>* gedPho2Pt_;
  std::vector<float>* gedPho3Pt_;
  std::vector<float>* gedPho4Pt_;
  std::vector<float>* gedPho5Pt_;
 
  std::vector<float>* gedPho0Eta_;
  std::vector<float>* gedPho1Eta_;
  std::vector<float>* gedPho2Eta_;
  std::vector<float>* gedPho3Eta_;
  std::vector<float>* gedPho4Eta_;
  std::vector<float>* gedPho5Eta_;

  std::vector<float>* gedPho0Phi_;
  std::vector<float>* gedPho1Phi_;
  std::vector<float>* gedPho2Phi_;
  std::vector<float>* gedPho3Phi_;
  std::vector<float>* gedPho4Phi_;
  std::vector<float>* gedPho5Phi_;

  std::vector<float>* gedPho0En_;
  std::vector<float>* gedPho1En_;
  std::vector<float>* gedPho2En_;
  std::vector<float>* gedPho3En_;
  std::vector<float>* gedPho4En_;
  std::vector<float>* gedPho5En_;


  // Muon variables
  std::vector<float>* MuPt_;
  std::vector<float>* MuEta_;
  std::vector<float>* MuPhi_;
  std::vector<float>* MuEn_;
  std::vector<float>* MuCh_;
  std::vector<float>* MuChi2_;


  //Muon sorted variables (not all of them)
  std::vector<float>* Mu0Pt_;
  std::vector<float>* Mu1Pt_;
  std::vector<float>* Mu2Pt_;
  std::vector<float>* Mu3Pt_;
  std::vector<float>* Mu4Pt_;
  std::vector<float>* Mu5Pt_;

  std::vector<float>* Mu0Eta_;
  std::vector<float>* Mu1Eta_;
  std::vector<float>* Mu2Eta_;
  std::vector<float>* Mu3Eta_;
  std::vector<float>* Mu4Eta_;
  std::vector<float>* Mu5Eta_;

  std::vector<float>* Mu0Phi_;
  std::vector<float>* Mu1Phi_;
  std::vector<float>* Mu2Phi_;
  std::vector<float>* Mu3Phi_;
  std::vector<float>* Mu4Phi_;
  std::vector<float>* Mu5Phi_;

  std::vector<float>* Mu0En_;
  std::vector<float>* Mu1En_;
  std::vector<float>* Mu2En_;
  std::vector<float>* Mu3En_;
  std::vector<float>* Mu4En_;
  std::vector<float>* Mu5En_;

  //Muon Cosmic variables
  std::vector<float>* MuCosmPt_;
  std::vector<float>* MuCosmEta_;
  std::vector<float>* MuCosmPhi_;
  std::vector<float>* MuCosmEn_;
  std::vector<float>* MuCosmCh_;
  std::vector<float>* MuCosmChi2_;


  //Muon sorted cosmic variables (not all of them)
  std::vector<float>* MuCosm0Pt_;
  std::vector<float>* MuCosm1Pt_;
  std::vector<float>* MuCosm2Pt_;
  std::vector<float>* MuCosm3Pt_;
  std::vector<float>* MuCosm4Pt_;
  std::vector<float>* MuCosm5Pt_;

  std::vector<float>* MuCosm0Eta_;
  std::vector<float>* MuCosm1Eta_;
  std::vector<float>* MuCosm2Eta_;
  std::vector<float>* MuCosm3Eta_;
  std::vector<float>* MuCosm4Eta_;
  std::vector<float>* MuCosm5Eta_;

  std::vector<float>* MuCosm0Phi_;
  std::vector<float>* MuCosm1Phi_;
  std::vector<float>* MuCosm2Phi_;
  std::vector<float>* MuCosm3Phi_;
  std::vector<float>* MuCosm4Phi_;
  std::vector<float>* MuCosm5Phi_;

  std::vector<float>* MuCosm0En_;
  std::vector<float>* MuCosm1En_;
  std::vector<float>* MuCosm2En_;
  std::vector<float>* MuCosm3En_;
  std::vector<float>* MuCosm4En_;
  std::vector<float>* MuCosm5En_;


  //Muon Cosmic1Leg variables
  std::vector<float>* MuCosmLegPt_;
  std::vector<float>* MuCosmLegEta_;
  std::vector<float>* MuCosmLegPhi_;
  std::vector<float>* MuCosmLegEn_;
  std::vector<float>* MuCosmLegCh_;
  std::vector<float>* MuCosmLegChi2_;


  //Muon Cosmic1Leg sorted variables (not all of them)
  std::vector<float>* MuCosmLeg0Pt_;
  std::vector<float>* MuCosmLeg1Pt_;
  std::vector<float>* MuCosmLeg2Pt_;
  std::vector<float>* MuCosmLeg3Pt_;
  std::vector<float>* MuCosmLeg4Pt_;
  std::vector<float>* MuCosmLeg5Pt_;

  std::vector<float>* MuCosmLeg0Eta_;
  std::vector<float>* MuCosmLeg1Eta_;
  std::vector<float>* MuCosmLeg2Eta_;
  std::vector<float>* MuCosmLeg3Eta_;
  std::vector<float>* MuCosmLeg4Eta_;
  std::vector<float>* MuCosmLeg5Eta_;

  std::vector<float>* MuCosmLeg0Phi_;
  std::vector<float>* MuCosmLeg1Phi_;
  std::vector<float>* MuCosmLeg2Phi_;
  std::vector<float>* MuCosmLeg3Phi_;
  std::vector<float>* MuCosmLeg4Phi_;
  std::vector<float>* MuCosmLeg5Phi_;

  std::vector<float>* MuCosmLeg0En_;
  std::vector<float>* MuCosmLeg1En_;
  std::vector<float>* MuCosmLeg2En_;
  std::vector<float>* MuCosmLeg3En_;
  std::vector<float>* MuCosmLeg4En_;
  std::vector<float>* MuCosmLeg5En_;

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
  std::vector<float>* EBiEta_;
  std::vector<float>* EBiPhi_;   

  std::vector<float>* EEenergy_;
  std::vector<float>* EEtime_;
  std::vector<float>* EEchi2_;
  std::vector<float>* EEix_;
  std::vector<float>* EEiy_; 

  std::vector<float>* ESenergy_;
  std::vector<float>* EStime_;
  // std::vector<float>* ESchi2_;
  std::vector<float>* ESix_;
  std::vector<float>* ESiy_; 
  //HBHE rechit, HF rechit and HO rechit variables
  std::vector<float>* HBHEenergy_;
  std::vector<float>* HBHEtime_;
  std::vector<float>* HBHEauxe_;
  std::vector<float>* HBHEieta_;
  std::vector<float>* HBHEiphi_;

  std::vector<float>* HFenergy_;
  std::vector<float>* HFtime_;
  // std::vector<float>* HFchi2_;
  std::vector<float>* HFieta_;
  std::vector<float>* HFiphi_;

  // std::vector<float>* HOenergy_;
  // std::vector<float>* HOtime_;
  // // std::vector<float>* HOchi2_;
  // std::vector<float>* HOieta_;
  // std::vector<float>* HOiphi_;

  
  //preshower variables
  std::vector<float>* PreShEn_;
  // std::vector<float>* PreShCorrEn_;
  std::vector<float>* PreShEta_;
  std::vector<float>* PreShPhi_;
  std::vector<float>* PreShYEn_;
  // std::vector<float>* PreShYCorrEn_;
  std::vector<float>* PreShYEta_;
  std::vector<float>* PreShYPhi_;

  //castor tower variables
  // std::vector<float>* CTPt_;
  // std::vector<float>* CTEta_;
  // std::vector<float>* CTPhi_;



  std::vector<float>* qPFJetPt_;
  std::vector<float>* qPFJetEta_;
  std::vector<float>* qPFJetPhi_;


  //PFJet sorted variables
  std::vector<float>* qPFJet0Pt_; 
  std::vector<float>* qPFJet1Pt_; 
  std::vector<float>* qPFJet2Pt_; 
  std::vector<float>* qPFJet3Pt_; 
  std::vector<float>* qPFJet4Pt_; 
  std::vector<float>* qPFJet5Pt_; 

  std::vector<float>* qPFJet0Eta_;
  std::vector<float>* qPFJet1Eta_;
  std::vector<float>* qPFJet2Eta_;
  std::vector<float>* qPFJet3Eta_;
  std::vector<float>* qPFJet4Eta_;
  std::vector<float>* qPFJet5Eta_;

  std::vector<float>* qPFJet0Phi_;
  std::vector<float>* qPFJet1Phi_;
  std::vector<float>* qPFJet2Phi_;
  std::vector<float>* qPFJet3Phi_;
  std::vector<float>* qPFJet4Phi_;
  std::vector<float>* qPFJet5Phi_;

  //PFJet4CHS sorted variables
  std::vector<float>* qPFJet4CHS0Pt_;
  std::vector<float>* qPFJet4CHS1Pt_;
  std::vector<float>* qPFJet4CHS2Pt_;
  std::vector<float>* qPFJet4CHS3Pt_;
  std::vector<float>* qPFJet4CHS4Pt_;
  std::vector<float>* qPFJet4CHS5Pt_;

  std::vector<float>* qPFJet4CHS0Eta_;
  std::vector<float>* qPFJet4CHS1Eta_;
  std::vector<float>* qPFJet4CHS2Eta_;
  std::vector<float>* qPFJet4CHS3Eta_;
  std::vector<float>* qPFJet4CHS4Eta_;
  std::vector<float>* qPFJet4CHS5Eta_;

  std::vector<float>* qPFJet4CHS0Phi_;
  std::vector<float>* qPFJet4CHS1Phi_;
  std::vector<float>* qPFJet4CHS2Phi_;
  std::vector<float>* qPFJet4CHS3Phi_;
  std::vector<float>* qPFJet4CHS4Phi_;
  std::vector<float>* qPFJet4CHS5Phi_;

  //PF8CHS sorted variables
  std::vector<float>* qPFJet8CHS0Pt_;
  std::vector<float>* qPFJet8CHS1Pt_;
  std::vector<float>* qPFJet8CHS2Pt_;
  std::vector<float>* qPFJet8CHS3Pt_;
  std::vector<float>* qPFJet8CHS4Pt_;
  std::vector<float>* qPFJet8CHS5Pt_;

  std::vector<float>* qPFJet8CHS0Eta_;
  std::vector<float>* qPFJet8CHS1Eta_;
  std::vector<float>* qPFJet8CHS2Eta_;
  std::vector<float>* qPFJet8CHS3Eta_;
  std::vector<float>* qPFJet8CHS4Eta_;
  std::vector<float>* qPFJet8CHS5Eta_;

  std::vector<float>* qPFJet8CHS0Phi_;
  std::vector<float>* qPFJet8CHS1Phi_;
  std::vector<float>* qPFJet8CHS2Phi_;
  std::vector<float>* qPFJet8CHS3Phi_;
  std::vector<float>* qPFJet8CHS4Phi_;
  std::vector<float>* qPFJet8CHS5Phi_;

  //PFJetEI sorted variables
  std::vector<float>* qPFJetEI0Pt_;
  std::vector<float>* qPFJetEI1Pt_;
  std::vector<float>* qPFJetEI2Pt_;
  std::vector<float>* qPFJetEI3Pt_;
  std::vector<float>* qPFJetEI4Pt_;
  std::vector<float>* qPFJetEI5Pt_;

  std::vector<float>* qPFJetEI0Eta_;
  std::vector<float>* qPFJetEI1Eta_;
  std::vector<float>* qPFJetEI2Eta_;
  std::vector<float>* qPFJetEI3Eta_;
  std::vector<float>* qPFJetEI4Eta_;
  std::vector<float>* qPFJetEI5Eta_;

  std::vector<float>* qPFJetEI0Phi_;
  std::vector<float>* qPFJetEI1Phi_;
  std::vector<float>* qPFJetEI2Phi_;
  std::vector<float>* qPFJetEI3Phi_;
  std::vector<float>* qPFJetEI4Phi_;
  std::vector<float>* qPFJetEI5Phi_;

  // //8CHSSoftDrop sorted variables
  // std::vector<float>* qPFJet8CHSSD0Pt_;
  // std::vector<float>* qPFJet8CHSSD1Pt_;
  // std::vector<float>* qPFJet8CHSSD2Pt_;
  // std::vector<float>* qPFJet8CHSSD3Pt_;
  // std::vector<float>* qPFJet8CHSSD4Pt_;
  // std::vector<float>* qPFJet8CHSSD5Pt_;

  // std::vector<float>* qPFJet8CHSSD0Eta_;
  // std::vector<float>* qPFJet8CHSSD1Eta_;
  // std::vector<float>* qPFJet8CHSSD2Eta_;
  // std::vector<float>* qPFJet8CHSSD3Eta_;
  // std::vector<float>* qPFJet8CHSSD4Eta_;
  // std::vector<float>* qPFJet8CHSSD5Eta_;

  // std::vector<float>* qPFJet8CHSSD0Phi_;
  // std::vector<float>* qPFJet8CHSSD1Phi_;
  // std::vector<float>* qPFJet8CHSSD2Phi_;
  // std::vector<float>* qPFJet8CHSSD3Phi_;
  // std::vector<float>* qPFJet8CHSSD4Phi_;
  // std::vector<float>* qPFJet8CHSSD5Phi_;

  // //TopCHS sorted variables
  // std::vector<float>* qPFJetTopCHS0Pt_;
  // std::vector<float>* qPFJetTopCHS1Pt_;
  // std::vector<float>* qPFJetTopCHS2Pt_;
  // std::vector<float>* qPFJetTopCHS3Pt_;
  // std::vector<float>* qPFJetTopCHS4Pt_;
  // std::vector<float>* qPFJetTopCHS5Pt_;

  // std::vector<float>* qPFJetTopCHS0Eta_;
  // std::vector<float>* qPFJetTopCHS1Eta_;
  // std::vector<float>* qPFJetTopCHS2Eta_;
  // std::vector<float>* qPFJetTopCHS3Eta_;
  // std::vector<float>* qPFJetTopCHS4Eta_;
  // std::vector<float>* qPFJetTopCHS5Eta_;

  // std::vector<float>* qPFJetTopCHS0Phi_;
  // std::vector<float>* qPFJetTopCHS1Phi_;
  // std::vector<float>* qPFJetTopCHS2Phi_;
  // std::vector<float>* qPFJetTopCHS3Phi_;
  // std::vector<float>* qPFJetTopCHS4Phi_;
  // std::vector<float>* qPFJetTopCHS5Phi_;

//______________________________________
  std::vector<float>* qPFJet4CHSPt_;
  std::vector<float>* qPFJet4CHSEta_;
  std::vector<float>* qPFJet4CHSPhi_;

  std::vector<float>* qPFJet8CHSPt_;
  std::vector<float>* qPFJet8CHSEta_;
  std::vector<float>* qPFJet8CHSPhi_;

  std::vector<float>* qPFJetEIPt_;
  std::vector<float>* qPFJetEIEta_;
  std::vector<float>* qPFJetEIPhi_;

  std::vector<float>* qPFJet8CHSSDPt_;
  std::vector<float>* qPFJet8CHSSDEta_;
  std::vector<float>* qPFJet8CHSSDPhi_;

  std::vector<float>* qPFJetTopCHSPt_;
  std::vector<float>* qPFJetTopCHSEta_;
  std::vector<float>* qPFJetTopCHSPhi_;

  std::vector<float>* qPFChMetPt_;
  // std::vector<float>* qPFChMetEta_;   
  std::vector<float>* qPFChMetPhi_;
  std::vector<float>* qPFMetPt_;
  // std::vector<float>* qPFMetEta_;  
  std::vector<float>* qPFMetPhi_;
  std::vector<int>*   qNVtx_;

  std::vector<float>* qCalJetPt_;
  std::vector<float>* qCalJetEta_;
  std::vector<float>* qCalJetPhi_;
  std::vector<float>* qCalJetEn_;


  //CaloJet sorted variables
  std::vector<float>* qCalJet0Pt_;
  std::vector<float>* qCalJet1Pt_;
  std::vector<float>* qCalJet2Pt_;
  std::vector<float>* qCalJet3Pt_;
  std::vector<float>* qCalJet4Pt_;
  std::vector<float>* qCalJet5Pt_;

  std::vector<float>* qCalJet0Eta_;
  std::vector<float>* qCalJet1Eta_;
  std::vector<float>* qCalJet2Eta_;
  std::vector<float>* qCalJet3Eta_;
  std::vector<float>* qCalJet4Eta_;
  std::vector<float>* qCalJet5Eta_;

  std::vector<float>* qCalJet0Phi_;
  std::vector<float>* qCalJet1Phi_;
  std::vector<float>* qCalJet2Phi_;
  std::vector<float>* qCalJet3Phi_;
  std::vector<float>* qCalJet4Phi_;
  std::vector<float>* qCalJet5Phi_;

  std::vector<float>* qCalJet0En_;
  std::vector<float>* qCalJet1En_;
  std::vector<float>* qCalJet2En_;
  std::vector<float>* qCalJet3En_;
  std::vector<float>* qCalJet4En_;
  std::vector<float>* qCalJet5En_;


  std::vector<float>* qCalMETPt_;
  // std::vector<float>* qCalMETEta_;
  std::vector<float>* qCalMETPhi_;
  std::vector<float>* qCalMETEn_;


  std::vector<float>* qCalMETBEPt_;
  // std::vector<float>* qCalMETBEEta_;
  std::vector<float>* qCalMETBEPhi_;
  std::vector<float>* qCalMETBEEn_;


  std::vector<float>* qCalMETBEFOPt_;
  // std::vector<float>* qCalMETBEFOEta_;
  std::vector<float>* qCalMETBEFOPhi_;
  std::vector<float>* qCalMETBEFOEn_;


  std::vector<float>* qCalMETMPt_;
  // std::vector<float>* qCalMETMEta_;
  std::vector<float>* qCalMETMPhi_;
  std::vector<float>* qCalMETMEn_;


  std::vector<float>* qSCEn_;
  std::vector<float>* qSCEta_;
  std::vector<float>* qSCPhi_;
  std::vector<float>* qSCEtaWidth_;
  std::vector<float>* qSCPhiWidth_;  
  std::vector<float>* qSCEnhfEM_;
  std::vector<float>* qSCEtahfEM_;
  std::vector<float>* qSCPhihfEM_;
  // std::vector<float>* qSCEtaWidthhfEM_;
  // std::vector<float>* qSCPhiWidthhfEM_; 
  std::vector<float>* qSCEn5x5_;
  std::vector<float>* qSCEta5x5_;
  std::vector<float>* qSCPhi5x5_;
  std::vector<float>* qSCEtaWidth5x5_;
  std::vector<float>* qSCPhiWidth5x5_;    
  std::vector<float>* qCCEn_;
  std::vector<float>* qCCEta_;
  std::vector<float>* qCCPhi_;
  std::vector<float>* qCCEn5x5_;
  std::vector<float>* qCCEta5x5_;
  std::vector<float>* qCCPhi5x5_;

  std::vector<float>* qPhoPt_;
  std::vector<float>* qPhoEta_;
  std::vector<float>* qPhoPhi_;
  std::vector<float>* qPhoEn_;


  //qPhoton sorted variables (not all of them)
  std::vector<float>* qPho0Pt_;
  std::vector<float>* qPho1Pt_;
  std::vector<float>* qPho2Pt_;
  std::vector<float>* qPho3Pt_;
  std::vector<float>* qPho4Pt_;
  std::vector<float>* qPho5Pt_;

  std::vector<float>* qPho0Eta_;
  std::vector<float>* qPho1Eta_;
  std::vector<float>* qPho2Eta_;
  std::vector<float>* qPho3Eta_;
  std::vector<float>* qPho4Eta_;
  std::vector<float>* qPho5Eta_;

  std::vector<float>* qPho0Phi_;
  std::vector<float>* qPho1Phi_;
  std::vector<float>* qPho2Phi_;
  std::vector<float>* qPho3Phi_;
  std::vector<float>* qPho4Phi_;
  std::vector<float>* qPho5Phi_;

  std::vector<float>* qPho0En_;
  std::vector<float>* qPho1En_;
  std::vector<float>* qPho2En_;
  std::vector<float>* qPho3En_;
  std::vector<float>* qPho4En_;
  std::vector<float>* qPho5En_;




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

  //ged qPhotons sorted variables (not all of them)
  std::vector<float>* qgedPho0Pt_;
  std::vector<float>* qgedPho1Pt_;
  std::vector<float>* qgedPho2Pt_;
  std::vector<float>* qgedPho3Pt_;
  std::vector<float>* qgedPho4Pt_;
  std::vector<float>* qgedPho5Pt_;
 
  std::vector<float>* qgedPho0Eta_;
  std::vector<float>* qgedPho1Eta_;
  std::vector<float>* qgedPho2Eta_;
  std::vector<float>* qgedPho3Eta_;
  std::vector<float>* qgedPho4Eta_;
  std::vector<float>* qgedPho5Eta_;

  std::vector<float>* qgedPho0Phi_;
  std::vector<float>* qgedPho1Phi_;
  std::vector<float>* qgedPho2Phi_;
  std::vector<float>* qgedPho3Phi_;
  std::vector<float>* qgedPho4Phi_;
  std::vector<float>* qgedPho5Phi_;

  std::vector<float>* qgedPho0En_;
  std::vector<float>* qgedPho1En_;
  std::vector<float>* qgedPho2En_;
  std::vector<float>* qgedPho3En_;
  std::vector<float>* qgedPho4En_;
  std::vector<float>* qgedPho5En_;


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
  std::vector<float>* qMuChi2_;
  std::vector<float>* qMuCosmPt_;
  std::vector<float>* qMuCosmEta_;
  std::vector<float>* qMuCosmPhi_;
  std::vector<float>* qMuCosmEn_;
  std::vector<float>* qMuCosmCh_;
  std::vector<float>* qMuCosmChi2_;
  std::vector<float>* qMuCosmLegPt_;
  std::vector<float>* qMuCosmLegEta_;
  std::vector<float>* qMuCosmLegPhi_;
  std::vector<float>* qMuCosmLegEn_;
  std::vector<float>* qMuCosmLegCh_;
  std::vector<float>* qMuCosmLegChi2_; 
  

  //Muon sorted variables (not all of them)
  std::vector<float>* qMu0Pt_;
  std::vector<float>* qMu1Pt_;
  std::vector<float>* qMu2Pt_;
  std::vector<float>* qMu3Pt_;
  std::vector<float>* qMu4Pt_;
  std::vector<float>* qMu5Pt_;

  std::vector<float>* qMu0Eta_;
  std::vector<float>* qMu1Eta_;
  std::vector<float>* qMu2Eta_;
  std::vector<float>* qMu3Eta_;
  std::vector<float>* qMu4Eta_;
  std::vector<float>* qMu5Eta_;

  std::vector<float>* qMu0Phi_;
  std::vector<float>* qMu1Phi_;
  std::vector<float>* qMu2Phi_;
  std::vector<float>* qMu3Phi_;
  std::vector<float>* qMu4Phi_;
  std::vector<float>* qMu5Phi_;

  std::vector<float>* qMu0En_;
  std::vector<float>* qMu1En_;
  std::vector<float>* qMu2En_;
  std::vector<float>* qMu3En_;
  std::vector<float>* qMu4En_;
  std::vector<float>* qMu5En_;
  //Muon sorted cosmic variables (not all of them)
  std::vector<float>* qMuCosm0Pt_;
  std::vector<float>* qMuCosm1Pt_;
  std::vector<float>* qMuCosm2Pt_;
  std::vector<float>* qMuCosm3Pt_;
  std::vector<float>* qMuCosm4Pt_;
  std::vector<float>* qMuCosm5Pt_;

  std::vector<float>* qMuCosm0Eta_;
  std::vector<float>* qMuCosm1Eta_;
  std::vector<float>* qMuCosm2Eta_;
  std::vector<float>* qMuCosm3Eta_;
  std::vector<float>* qMuCosm4Eta_;
  std::vector<float>* qMuCosm5Eta_;

  std::vector<float>* qMuCosm0Phi_;
  std::vector<float>* qMuCosm1Phi_;
  std::vector<float>* qMuCosm2Phi_;
  std::vector<float>* qMuCosm3Phi_;
  std::vector<float>* qMuCosm4Phi_;
  std::vector<float>* qMuCosm5Phi_;

  std::vector<float>* qMuCosm0En_;
  std::vector<float>* qMuCosm1En_;
  std::vector<float>* qMuCosm2En_;
  std::vector<float>* qMuCosm3En_;
  std::vector<float>* qMuCosm4En_;
  std::vector<float>* qMuCosm5En_;

  //Muon Cosmic1Leg sorted variables (not all of them)
  std::vector<float>* qMuCosmLeg0Pt_;
  std::vector<float>* qMuCosmLeg1Pt_;
  std::vector<float>* qMuCosmLeg2Pt_;
  std::vector<float>* qMuCosmLeg3Pt_;
  std::vector<float>* qMuCosmLeg4Pt_;
  std::vector<float>* qMuCosmLeg5Pt_;

  std::vector<float>* qMuCosmLeg0Eta_;
  std::vector<float>* qMuCosmLeg1Eta_;
  std::vector<float>* qMuCosmLeg2Eta_;
  std::vector<float>* qMuCosmLeg3Eta_;
  std::vector<float>* qMuCosmLeg4Eta_;
  std::vector<float>* qMuCosmLeg5Eta_;

  std::vector<float>* qMuCosmLeg0Phi_;
  std::vector<float>* qMuCosmLeg1Phi_;
  std::vector<float>* qMuCosmLeg2Phi_;
  std::vector<float>* qMuCosmLeg3Phi_;
  std::vector<float>* qMuCosmLeg4Phi_;
  std::vector<float>* qMuCosmLeg5Phi_;

  std::vector<float>* qMuCosmLeg0En_;
  std::vector<float>* qMuCosmLeg1En_;
  std::vector<float>* qMuCosmLeg2En_;
  std::vector<float>* qMuCosmLeg3En_;
  std::vector<float>* qMuCosmLeg4En_;
  std::vector<float>* qMuCosmLeg5En_;



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
  std::vector<float>* qEBiEta_;
  std::vector<float>* qEBiPhi_;   

  std::vector<float>* qEEenergy_;
  std::vector<float>* qEEtime_;
  std::vector<float>* qEEchi2_;
  std::vector<float>* qEEix_;
  std::vector<float>* qEEiy_; 
   
  std::vector<float>* qESenergy_;
  std::vector<float>* qEStime_;
  // std::vector<float>* qESchi2_;
  std::vector<float>* qESix_;
  std::vector<float>* qESiy_; 

  std::vector<float>* qHBHEenergy_;
  std::vector<float>* qHBHEtime_;
  std::vector<float>* qHBHEauxe_;
  std::vector<float>* qHBHEieta_;
  std::vector<float>* qHBHEiphi_;

  std::vector<float>* qHFenergy_;
  std::vector<float>* qHFtime_;
  // std::vector<float>* qHFchi2_;
  std::vector<float>* qHFieta_;
  std::vector<float>* qHFiphi_;

  // std::vector<float>* qHOenergy_;
  // std::vector<float>* qHOtime_;
  // // std::vector<float>* qHOchi2_;
  // std::vector<float>* qHOieta_;
  // std::vector<float>* qHOiphi_;

  std::vector<float>* qPreShEn_;
  // std::vector<float>* qPreShCorrEn_;
  std::vector<float>* qPreShEta_;
  std::vector<float>* qPreShPhi_;
  std::vector<float>* qPreShYEn_;
  // std::vector<float>* qPreShYCorrEn_;
  std::vector<float>* qPreShYEta_;
  std::vector<float>* qPreShYPhi_;

  // std::vector<float>* qCTPt_;
  // std::vector<float>* qCTEta_;
  // std::vector<float>* qCTPhi_;

  std::vector<float>*   crossSection_;
  std::vector<float>*   pathRates_;
  std::vector<std::string>*   pathNames_;
  std::map<std::string,int> rateMap;




  edm::EDGetTokenT<reco::PFJetCollection> PFJetToken_;

  //----------------------------------------Sorted variables start here.




  //________________________________________Sorted variables end here.

  edm::EDGetTokenT<reco::PFJetCollection> PFJet4CHSToken_;
  edm::EDGetTokenT<reco::PFJetCollection> PFJet8CHSToken_;
  edm::EDGetTokenT<reco::PFJetCollection> PFJetEIToken_;
  edm::EDGetTokenT<reco::PFJetCollection> PFJet8CHSSoftDropToken_;
  edm::EDGetTokenT<reco::PFJetCollection> PFJetTopCHSToken_;
  edm::EDGetTokenT<reco::PFMETCollection> PFChMETToken_;
  edm::EDGetTokenT<reco::PFMETCollection> PFMETToken_;

  edm::EDGetTokenT<reco::CaloJetCollection> CaloJetToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> CaloMETToken_;      
  edm::EDGetTokenT<reco::CaloMETCollection> CaloMETBEToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> CaloMETBEFOToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> CaloMETMToken_;
  edm::EDGetTokenT<reco::VertexCollection>  vtxToken_;

  edm::EDGetTokenT<edm::TriggerResults>     triggerBits_;
  // edm::EDGetTokenT<edm::TriggerResults>  triggerResultsToken_;
  //edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;  //pat::PackedTriggerPrescales is only present in miniAOD. Have to find alternative using TriggerResults
  

  edm::EDGetTokenT<reco::SuperClusterCollection>   SuperClusterToken_;  
  edm::EDGetTokenT<reco::SuperClusterCollection>    SuperClusterhfEMToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection>   SuperCluster5x5Token_;  

  edm::EDGetTokenT<reco::CaloClusterCollection>   CaloClusterToken_;  
  edm::EDGetTokenT<reco::CaloClusterCollection>   CaloCluster5x5Token_;

  edm::EDGetTokenT<reco::PhotonCollection> PhotonToken_;  
  edm::EDGetTokenT<reco::PhotonCollection> gedPhotonToken_;
  edm::EDGetTokenT<reco::MuonCollection>   MuonToken_;
  edm::EDGetTokenT<reco::MuonCollection>   MuonCosmToken_;
  edm::EDGetTokenT<reco::MuonCollection>   MuonCosmLegToken_;

  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> GsfElectronUncleanedToken_;

  edm::EDGetTokenT<EcalRecHitCollection> ebRHSrcToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRHSrcToken_;
  edm::EDGetTokenT<EcalRecHitCollection> esRHSrcToken_;

  edm::EDGetTokenT<HBHERecHitCollection> hbheRHcToken_;
  edm::EDGetTokenT<HFRecHitCollection>   hfRHcToken_;
  // edm::EDGetTokenT<HORecHitCollection>   hoRHcToken_;

  // edm::EDGetTokenT<SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> ebRHSrcToken_;
  // edm::EDGetTokenT<SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> eeRHSrcToken_;
  edm::EDGetTokenT<reco::PreshowerClusterCollection> preshowerXToken_;
  edm::EDGetTokenT<reco::PreshowerClusterCollection> preshowerYToken_;
  // edm::EDGetTokenT<reco::CastorTowerCollection> CastorTowerToken_;   //leafcandidate variables  //perhaps reco




  int eventCounter;

  // double maxJetEta_;
  // double minJetPt_;
  // double maxSCEta_;
  // double minSCEn_;

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
  PFJetToken_               (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJetTag"))),  
  PFJet4CHSToken_           (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJet4CHSTag"))),
  PFJet8CHSToken_           (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJet8CHSTag"))),
  PFJetEIToken_             (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJetEITag"))),
  PFJet8CHSSoftDropToken_   (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJet8CHSSoftDropTag"))),
  PFJetTopCHSToken_         (consumes<reco::PFJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJetTopCHSTag"))),
  PFChMETToken_             (consumes<reco::PFMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFChMETTag"))),
  PFMETToken_               (consumes<reco::PFMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFMETTag"))),
  CaloJetToken_             (consumes<reco::CaloJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloJetTag"))),
  CaloMETToken_             (consumes<reco::CaloMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloMETTag"))),  
  CaloMETBEToken_           (consumes<reco::CaloMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloMETBETag"))),  
  CaloMETBEFOToken_         (consumes<reco::CaloMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloMETBEFOTag"))),  
  CaloMETMToken_            (consumes<reco::CaloMETCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloMETMTag"))), 
 
  vtxToken_                 (consumes<reco::VertexCollection>(cfg.getUntrackedParameter<edm::InputTag>("vtx"))),
  triggerBits_              (consumes<edm::TriggerResults>(cfg.getUntrackedParameter<edm::InputTag>("bits"))),
  // triggerResultsToken_      (consumes<edm::TriggerResults>(cfg.getUntrackedParameter<edm::InputTag>("triggerResultsTag"))),  // alternative for pat::PackedTriggerPrescales
  //triggerPrescales_         (consumes<pat::PackedTriggerPrescales>(cfg.getUntrackedParameter<edm::InputTag>("prescales"))),  // pat::PackedTriggerPrescales is not in AOD
  SuperClusterToken_        (consumes<reco::SuperClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("SuperClusterTag"))),
  SuperClusterhfEMToken_    (consumes<reco::SuperClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("SuperClusterhfEMTag"))),
  SuperCluster5x5Token_     (consumes<reco::SuperClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("SuperCluster5x5Tag"))),
  CaloClusterToken_         (consumes<reco::CaloClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloClusterTag"))),
  CaloCluster5x5Token_      (consumes<reco::CaloClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("CaloCluster5x5Tag"))),
  PhotonToken_              (consumes<reco::PhotonCollection>(cfg.getUntrackedParameter<edm::InputTag>("PhotonTag"))),
  gedPhotonToken_           (consumes<reco::PhotonCollection>(cfg.getUntrackedParameter<edm::InputTag>("gedPhotonTag"))),
  MuonToken_                (consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonTag"))),
  MuonCosmToken_            (consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonCosmTag"))),
  MuonCosmLegToken_         (consumes<reco::MuonCollection>(cfg.getUntrackedParameter<edm::InputTag>("MuonCosmLegTag"))),

  GsfElectronToken_         (consumes<reco::GsfElectronCollection>(cfg.getUntrackedParameter<edm::InputTag>("GsfElectronTag"))),
  GsfElectronUncleanedToken_(consumes<reco::GsfElectronCollection>(cfg.getUntrackedParameter<edm::InputTag>("GsfElectronUncleanedTag"))),

  ebRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("EBRecHitSourceTag"))),  
  eeRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("EERecHitSourceTag"))),
  esRHSrcToken_             (consumes<EcalRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("ESRecHitSourceTag"))),

  hbheRHcToken_             (consumes<HBHERecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("HBHERecHitTag"))),
  hfRHcToken_               (consumes<HFRecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("HFRecHitTag"))),
  // hoRHcToken_               (consumes<HORecHitCollection>(cfg.getUntrackedParameter<edm::InputTag>("HORecHitTag"))),
  preshowerXToken_          (consumes<reco::PreshowerClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("PreshowerClusterXTag"))),
  preshowerYToken_          (consumes<reco::PreshowerClusterCollection>(cfg.getUntrackedParameter<edm::InputTag>("PreshowerClusterYTag"))),
  // CastorTowerToken_         (consumes<reco::CastorTowerCollection>(cfg.getUntrackedParameter<edm::InputTag>("CastorTowerTag"))),  //perhaps reco  https://github.com/cms-sw/cmssw/blob/09c3fce6626f70fd04223e7dacebf0b485f73f54/DataFormats/CastorReco/interface/CastorTower.h

  // hltPrescaleProvider_(cfg, consumesCollector(), *this),

  //params for wide jet calculation
  // maxJetEta_                (cfg.getUntrackedParameter<double>("maxJetEta")),
  // minJetPt_                 (cfg.getUntrackedParameter<double>("minJetPt")),
  // maxSCEta_                 (cfg.getUntrackedParameter<double>("maxSCEta")),
  // minSCEn_                  (cfg.getUntrackedParameter<double>("minSCEn")),
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
//harambe
  // hits_->clear();


  PFJetPt_->clear();
  PFJetEta_->clear();
  PFJetPhi_->clear();

//---------------------- All sorted variables
//PFJet sorted variables
  PFJet0Pt_->clear(); 
  PFJet1Pt_->clear(); 
  PFJet2Pt_->clear(); 
  PFJet3Pt_->clear(); 
  PFJet4Pt_->clear(); 
  PFJet5Pt_->clear(); 

  PFJet0Eta_->clear();
  PFJet1Eta_->clear();
  PFJet2Eta_->clear();
  PFJet3Eta_->clear();
  PFJet4Eta_->clear();
  PFJet5Eta_->clear();

  PFJet0Phi_->clear();
  PFJet1Phi_->clear();
  PFJet2Phi_->clear();
  PFJet3Phi_->clear();
  PFJet4Phi_->clear();
  PFJet5Phi_->clear();
//PFJet4CHS sorted variables
  PFJet4CHS0Pt_->clear();
  PFJet4CHS1Pt_->clear();
  PFJet4CHS2Pt_->clear();
  PFJet4CHS3Pt_->clear();
  PFJet4CHS4Pt_->clear();
  PFJet4CHS5Pt_->clear();

  PFJet4CHS0Eta_->clear();
  PFJet4CHS1Eta_->clear();
  PFJet4CHS2Eta_->clear();
  PFJet4CHS3Eta_->clear();
  PFJet4CHS4Eta_->clear();
  PFJet4CHS5Eta_->clear();

  PFJet4CHS0Phi_->clear();
  PFJet4CHS1Phi_->clear();
  PFJet4CHS2Phi_->clear();
  PFJet4CHS3Phi_->clear();
  PFJet4CHS4Phi_->clear();
  PFJet4CHS5Phi_->clear();
 //PF8CHS sorted variables
  PFJet8CHS0Pt_->clear();
  PFJet8CHS1Pt_->clear();
  PFJet8CHS2Pt_->clear();
  PFJet8CHS3Pt_->clear();
  PFJet8CHS4Pt_->clear();
  PFJet8CHS5Pt_->clear();

  PFJet8CHS0Eta_->clear();
  PFJet8CHS1Eta_->clear();
  PFJet8CHS2Eta_->clear();
  PFJet8CHS3Eta_->clear();
  PFJet8CHS4Eta_->clear();
  PFJet8CHS5Eta_->clear();

  PFJet8CHS0Phi_->clear();
  PFJet8CHS1Phi_->clear();
  PFJet8CHS2Phi_->clear();
  PFJet8CHS3Phi_->clear();
  PFJet8CHS4Phi_->clear();
  PFJet8CHS5Phi_->clear();
  //PFJetEI sorted variables
  PFJetEI0Pt_->clear();
  PFJetEI1Pt_->clear();
  PFJetEI2Pt_->clear();
  PFJetEI3Pt_->clear();
  PFJetEI4Pt_->clear();
  PFJetEI5Pt_->clear();

  PFJetEI0Eta_->clear();
  PFJetEI1Eta_->clear();
  PFJetEI2Eta_->clear();
  PFJetEI3Eta_->clear();
  PFJetEI4Eta_->clear();
  PFJetEI5Eta_->clear();

  PFJetEI0Phi_->clear();
  PFJetEI1Phi_->clear();
  PFJetEI2Phi_->clear();
  PFJetEI3Phi_->clear();
  PFJetEI4Phi_->clear();
  PFJetEI5Phi_->clear();

  // //8CHSSoftDrop sorted variables
  // PFJet8CHSSD0Pt_->clear();
  // PFJet8CHSSD1Pt_->clear();
  // PFJet8CHSSD2Pt_->clear();
  // PFJet8CHSSD3Pt_->clear();
  // PFJet8CHSSD4Pt_->clear();
  // PFJet8CHSSD5Pt_->clear();

  // PFJet8CHSSD0Eta_->clear();
  // PFJet8CHSSD1Eta_->clear();
  // PFJet8CHSSD2Eta_->clear();
  // PFJet8CHSSD3Eta_->clear();
  // PFJet8CHSSD4Eta_->clear();
  // PFJet8CHSSD5Eta_->clear();

  // PFJet8CHSSD0Phi_->clear();
  // PFJet8CHSSD1Phi_->clear();
  // PFJet8CHSSD2Phi_->clear();
  // PFJet8CHSSD3Phi_->clear();
  // PFJet8CHSSD4Phi_->clear();
  // PFJet8CHSSD5Phi_->clear();
  // //TopCHS sorted variables
  // PFJetTopCHS0Pt_->clear();
  // PFJetTopCHS1Pt_->clear();
  // PFJetTopCHS2Pt_->clear();
  // PFJetTopCHS3Pt_->clear();
  // PFJetTopCHS4Pt_->clear();
  // PFJetTopCHS5Pt_->clear();

  // PFJetTopCHS0Eta_->clear();
  // PFJetTopCHS1Eta_->clear();
  // PFJetTopCHS2Eta_->clear();
  // PFJetTopCHS3Eta_->clear();
  // PFJetTopCHS4Eta_->clear();
  // PFJetTopCHS5Eta_->clear();

  // PFJetTopCHS0Phi_->clear();
  // PFJetTopCHS1Phi_->clear();
  // PFJetTopCHS2Phi_->clear();
  // PFJetTopCHS3Phi_->clear();
  // PFJetTopCHS4Phi_->clear();
  // PFJetTopCHS5Phi_->clear();


  //CaloJet sorted variables
  CalJet0Pt_->clear();
  CalJet1Pt_->clear();
  CalJet2Pt_->clear();
  CalJet3Pt_->clear();
  CalJet4Pt_->clear();
  CalJet5Pt_->clear();

  CalJet0Eta_->clear();
  CalJet1Eta_->clear();
  CalJet2Eta_->clear();
  CalJet3Eta_->clear();
  CalJet4Eta_->clear();
  CalJet5Eta_->clear();

  CalJet0Phi_->clear();
  CalJet1Phi_->clear();
  CalJet2Phi_->clear();
  CalJet3Phi_->clear();
  CalJet4Phi_->clear();
  CalJet5Phi_->clear();

  CalJet0En_->clear();
  CalJet1En_->clear();
  CalJet2En_->clear();
  CalJet3En_->clear();
  CalJet4En_->clear();
  CalJet5En_->clear();


  //photon sorted variables (not all of them)
  Pho0Pt_->clear();
  Pho1Pt_->clear();
  Pho2Pt_->clear();
  Pho3Pt_->clear();
  Pho4Pt_->clear();
  Pho5Pt_->clear();

  Pho0Eta_->clear();
  Pho1Eta_->clear();
  Pho2Eta_->clear();
  Pho3Eta_->clear();
  Pho4Eta_->clear();
  Pho5Eta_->clear();

  Pho0Phi_->clear();
  Pho1Phi_->clear();
  Pho2Phi_->clear();
  Pho3Phi_->clear();
  Pho4Phi_->clear();
  Pho5Phi_->clear();

  Pho0En_->clear();
  Pho1En_->clear();
  Pho2En_->clear();
  Pho3En_->clear();
  Pho4En_->clear();
  Pho5En_->clear();

  //ged Photons sorted variables (not all of them)
  gedPho0Pt_->clear();
  gedPho1Pt_->clear();
  gedPho2Pt_->clear();
  gedPho3Pt_->clear();
  gedPho4Pt_->clear();
  gedPho5Pt_->clear();
 
  gedPho0Eta_->clear();
  gedPho1Eta_->clear();
  gedPho2Eta_->clear();
  gedPho3Eta_->clear();
  gedPho4Eta_->clear();
  gedPho5Eta_->clear();

  gedPho0Phi_->clear();
  gedPho1Phi_->clear();
  gedPho2Phi_->clear();
  gedPho3Phi_->clear();
  gedPho4Phi_->clear();
  gedPho5Phi_->clear();

  gedPho0En_->clear();
  gedPho1En_->clear();
  gedPho2En_->clear();
  gedPho3En_->clear();
  gedPho4En_->clear();
  gedPho5En_->clear();


    //Muon sorted variables (not all of them)
  Mu0Pt_->clear();
  Mu1Pt_->clear();
  Mu2Pt_->clear();
  Mu3Pt_->clear();
  Mu4Pt_->clear();
  Mu5Pt_->clear();

  Mu0Eta_->clear();
  Mu1Eta_->clear();
  Mu2Eta_->clear();
  Mu3Eta_->clear();
  Mu4Eta_->clear();
  Mu5Eta_->clear();

  Mu0Phi_->clear();
  Mu1Phi_->clear();
  Mu2Phi_->clear();
  Mu3Phi_->clear();
  Mu4Phi_->clear();
  Mu5Phi_->clear();

  Mu0En_->clear();
  Mu1En_->clear();
  Mu2En_->clear();
  Mu3En_->clear();
  Mu4En_->clear();
  Mu5En_->clear();
  //Muon sorted cosmic variables (not all of them)
  MuCosm0Pt_->clear();
  MuCosm1Pt_->clear();
  MuCosm2Pt_->clear();
  MuCosm3Pt_->clear();
  MuCosm4Pt_->clear();
  MuCosm5Pt_->clear();

  MuCosm0Eta_->clear();
  MuCosm1Eta_->clear();
  MuCosm2Eta_->clear();
  MuCosm3Eta_->clear();
  MuCosm4Eta_->clear();
  MuCosm5Eta_->clear();

  MuCosm0Phi_->clear();
  MuCosm1Phi_->clear();
  MuCosm2Phi_->clear();
  MuCosm3Phi_->clear();
  MuCosm4Phi_->clear();
  MuCosm5Phi_->clear();

  MuCosm0En_->clear();
  MuCosm1En_->clear();
  MuCosm2En_->clear();
  MuCosm3En_->clear();
  MuCosm4En_->clear();
  MuCosm5En_->clear();


  //Muon Cosmic1Leg sorted variables (not all of them)
  MuCosmLeg0Pt_->clear();
  MuCosmLeg1Pt_->clear();
  MuCosmLeg2Pt_->clear();
  MuCosmLeg3Pt_->clear();
  MuCosmLeg4Pt_->clear();
  MuCosmLeg5Pt_->clear();

  MuCosmLeg0Eta_->clear();
  MuCosmLeg1Eta_->clear();
  MuCosmLeg2Eta_->clear();
  MuCosmLeg3Eta_->clear();
  MuCosmLeg4Eta_->clear();
  MuCosmLeg5Eta_->clear();

  MuCosmLeg0Phi_->clear();
  MuCosmLeg1Phi_->clear();
  MuCosmLeg2Phi_->clear();
  MuCosmLeg3Phi_->clear();
  MuCosmLeg4Phi_->clear();
  MuCosmLeg5Phi_->clear();

  MuCosmLeg0En_->clear();
  MuCosmLeg1En_->clear();
  MuCosmLeg2En_->clear();
  MuCosmLeg3En_->clear();
  MuCosmLeg4En_->clear();
  MuCosmLeg5En_->clear();

//___________________________All sorted variables end.



  PFJet4CHSPt_->clear();
  PFJet4CHSEta_->clear();
  PFJet4CHSPhi_->clear();

  PFJet8CHSPt_->clear();
  PFJet8CHSEta_->clear();
  PFJet8CHSPhi_->clear();
  
  PFJetEIPt_->clear();
  PFJetEIEta_->clear();
  PFJetEIPhi_->clear();

  PFJet8CHSSDPt_->clear();
  PFJet8CHSSDEta_->clear();
  PFJet8CHSSDPhi_->clear();

  PFJetTopCHSPt_->clear();
  PFJetTopCHSEta_->clear();
  PFJetTopCHSPhi_->clear();       

  PFChMetPt_->clear();
  // PFChMetEta_->clear();  //MET doesn't have ETA
  PFChMetPhi_->clear();
  PFMetPt_->clear();
  // PFMetEta_->clear();  //MET doesn't have ETA
  PFMetPhi_->clear();
  nVtx_->clear();

  CalJetPt_->clear();
  CalJetEta_->clear();
  CalJetPhi_->clear();
  CalJetEn_->clear();

  CalMETPt_->clear();
  // CalMETEta_->clear();  //MET doesn't have ETA
  CalMETPhi_->clear();
  CalMETEn_->clear();

  CalMETBEPt_->clear();
  // CalMETBEEta_->clear();  //MET doesn't have ETA
  CalMETBEPhi_->clear();
  CalMETBEEn_->clear();

  CalMETBEFOPt_->clear();
  // CalMETBEFOEta_->clear();  //MET doesn't have ETA
  CalMETBEFOPhi_->clear();
  CalMETBEFOEn_->clear();

  CalMETMPt_->clear();
  // CalMETMEta_->clear();  //MET doesn't have ETA
  CalMETMPhi_->clear();
  CalMETMEn_->clear();

  SCEn_ ->clear();
  SCEta_->clear();
  SCPhi_->clear();
  SCEtaWidth_->clear();
  SCPhiWidth_->clear();  
  SCEnhfEM_ ->clear();
  SCEtahfEM_->clear();
  SCPhihfEM_->clear();
  // SCEtaWidthhfEM_->clear();  //outputted 0
  // SCPhiWidthhfEM_->clear();  //outputted 0
  SCEn5x5_ ->clear();
  SCEta5x5_->clear();
  SCPhi5x5_->clear();
  SCEtaWidth5x5_->clear();
  SCPhiWidth5x5_->clear();    
  CCEn_ ->clear();
  CCEta_->clear();
  CCPhi_->clear();
  CCEn5x5_ ->clear();
  CCEta5x5_->clear();
  CCPhi5x5_->clear();

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
  MuChi2_->clear();
  MuCosmPt_->clear();
  MuCosmEta_->clear();
  MuCosmPhi_->clear();
  MuCosmEn_->clear();
  MuCosmCh_->clear();
  MuCosmChi2_->clear();  
  MuCosmLegPt_->clear();
  MuCosmLegEta_->clear();
  MuCosmLegPhi_->clear();
  MuCosmLegEn_->clear();
  MuCosmLegCh_->clear();
  MuCosmLegChi2_->clear();
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
  EBiEta_->clear();
  EBiPhi_->clear();

  EEenergy_->clear();
  EEtime_->clear();
  EEchi2_->clear();
  EEix_->clear();
  EEiy_->clear();

  ESenergy_->clear();
  EStime_->clear();
  // ESchi2_->clear();
  ESix_->clear();
  ESiy_->clear();

  HBHEenergy_->clear();
  HBHEtime_->clear();
  HBHEauxe_->clear();
  HBHEieta_->clear();
  HBHEiphi_->clear();

  HFenergy_->clear();
  HFtime_->clear();
  // HFchi2_->clear();
  HFieta_->clear();
  HFiphi_->clear();

  // HOenergy_->clear();
  // HOtime_->clear();
  // // HOchi2_->clear();
  // HOieta_->clear();
  // HOiphi_->clear();

 
  PreShEn_->clear();
  // PreShCorrEn_->clear();
  PreShEta_->clear();
  PreShPhi_->clear();
  PreShYEn_->clear();
  // PreShYCorrEn_->clear();
  PreShYEta_->clear();
  PreShYPhi_->clear();

  // CTPt_->clear();
  // CTEta_->clear();
  // CTPhi_->clear();

  qPFJetPt_->clear();
  qPFJetEta_->clear();
  qPFJetPhi_->clear();

//---------------------- All sorted variables
//PFJet sorted variables
  qPFJet0Pt_->clear(); 
  qPFJet1Pt_->clear(); 
  qPFJet2Pt_->clear(); 
  qPFJet3Pt_->clear(); 
  qPFJet4Pt_->clear(); 
  qPFJet5Pt_->clear(); 

  qPFJet0Eta_->clear();
  qPFJet1Eta_->clear();
  qPFJet2Eta_->clear();
  qPFJet3Eta_->clear();
  qPFJet4Eta_->clear();
  qPFJet5Eta_->clear();

  qPFJet0Phi_->clear();
  qPFJet1Phi_->clear();
  qPFJet2Phi_->clear();
  qPFJet3Phi_->clear();
  qPFJet4Phi_->clear();
  qPFJet5Phi_->clear();
//PFJet4CHS sorted variables
  qPFJet4CHS0Pt_->clear();
  qPFJet4CHS1Pt_->clear();
  qPFJet4CHS2Pt_->clear();
  qPFJet4CHS3Pt_->clear();
  qPFJet4CHS4Pt_->clear();
  qPFJet4CHS5Pt_->clear();
  
  qPFJet4CHS0Eta_->clear();
  qPFJet4CHS1Eta_->clear();
  qPFJet4CHS2Eta_->clear();
  qPFJet4CHS3Eta_->clear();
  qPFJet4CHS4Eta_->clear();
  qPFJet4CHS5Eta_->clear();

  qPFJet4CHS0Phi_->clear();
  qPFJet4CHS1Phi_->clear();
  qPFJet4CHS2Phi_->clear();
  qPFJet4CHS3Phi_->clear();
  qPFJet4CHS4Phi_->clear();
  qPFJet4CHS5Phi_->clear();
 //PF8CHS sorted variables
  qPFJet8CHS0Pt_->clear();
  qPFJet8CHS1Pt_->clear();
  qPFJet8CHS2Pt_->clear();
  qPFJet8CHS3Pt_->clear();
  qPFJet8CHS4Pt_->clear();
  qPFJet8CHS5Pt_->clear();

  qPFJet8CHS0Eta_->clear();
  qPFJet8CHS1Eta_->clear();
  qPFJet8CHS2Eta_->clear();
  qPFJet8CHS3Eta_->clear();
  qPFJet8CHS4Eta_->clear();
  qPFJet8CHS5Eta_->clear();

  qPFJet8CHS0Phi_->clear();
  qPFJet8CHS1Phi_->clear();
  qPFJet8CHS2Phi_->clear();
  qPFJet8CHS3Phi_->clear();
  qPFJet8CHS4Phi_->clear();
  qPFJet8CHS5Phi_->clear();
  //PFJetEI sorted variables
  qPFJetEI0Pt_->clear();
  qPFJetEI1Pt_->clear();
  qPFJetEI2Pt_->clear();
  qPFJetEI3Pt_->clear();
  qPFJetEI4Pt_->clear();
  qPFJetEI5Pt_->clear();

  qPFJetEI0Eta_->clear();
  qPFJetEI1Eta_->clear();
  qPFJetEI2Eta_->clear();
  qPFJetEI3Eta_->clear();
  qPFJetEI4Eta_->clear();
  qPFJetEI5Eta_->clear();

  qPFJetEI0Phi_->clear();
  qPFJetEI1Phi_->clear();
  qPFJetEI2Phi_->clear();
  qPFJetEI3Phi_->clear();
  qPFJetEI4Phi_->clear();
  qPFJetEI5Phi_->clear();

  // //8CHSSoftDrop sorted variables
  // qPFJet8CHSSD0Pt_->clear();
  // qPFJet8CHSSD1Pt_->clear();
  // qPFJet8CHSSD2Pt_->clear();
  // qPFJet8CHSSD3Pt_->clear();
  // qPFJet8CHSSD4Pt_->clear();
  // qPFJet8CHSSD5Pt_->clear();

  // qPFJet8CHSSD0Eta_->clear();
  // qPFJet8CHSSD1Eta_->clear();
  // qPFJet8CHSSD2Eta_->clear();
  // qPFJet8CHSSD3Eta_->clear();
  // qPFJet8CHSSD4Eta_->clear();
  // qPFJet8CHSSD5Eta_->clear();

  // qPFJet8CHSSD0Phi_->clear();
  // qPFJet8CHSSD1Phi_->clear();
  // qPFJet8CHSSD2Phi_->clear();
  // qPFJet8CHSSD3Phi_->clear();
  // qPFJet8CHSSD4Phi_->clear();
  // qPFJet8CHSSD5Phi_->clear();
  // //TopCHS sorted variables
  // qPFJetTopCHS0Pt_->clear();
  // qPFJetTopCHS1Pt_->clear();
  // qPFJetTopCHS2Pt_->clear();
  // qPFJetTopCHS3Pt_->clear();
  // qPFJetTopCHS4Pt_->clear();
  // qPFJetTopCHS5Pt_->clear();

  // qPFJetTopCHS0Eta_->clear();
  // qPFJetTopCHS1Eta_->clear();
  // qPFJetTopCHS2Eta_->clear();
  // qPFJetTopCHS3Eta_->clear();
  // qPFJetTopCHS4Eta_->clear();
  // qPFJetTopCHS5Eta_->clear();

  // qPFJetTopCHS0Phi_->clear();
  // qPFJetTopCHS1Phi_->clear();
  // qPFJetTopCHS2Phi_->clear();
  // qPFJetTopCHS3Phi_->clear();
  // qPFJetTopCHS4Phi_->clear();
  // qPFJetTopCHS5Phi_->clear();


  //CaloJet sorted variables
  qCalJet0Pt_->clear();
  qCalJet1Pt_->clear();
  qCalJet2Pt_->clear();
  qCalJet3Pt_->clear();
  qCalJet4Pt_->clear();
  qCalJet5Pt_->clear();

  qCalJet0Eta_->clear();
  qCalJet1Eta_->clear();
  qCalJet2Eta_->clear();
  qCalJet3Eta_->clear();
  qCalJet4Eta_->clear();
  qCalJet5Eta_->clear();

  qCalJet0Phi_->clear();
  qCalJet1Phi_->clear();
  qCalJet2Phi_->clear();
  qCalJet3Phi_->clear();
  qCalJet4Phi_->clear();
  qCalJet5Phi_->clear();

  qCalJet0En_->clear();
  qCalJet1En_->clear();
  qCalJet2En_->clear();
  qCalJet3En_->clear();
  qCalJet4En_->clear();
  qCalJet5En_->clear();


  //photon sorted variables (not all of them)
  qPho0Pt_->clear();
  qPho1Pt_->clear();
  qPho2Pt_->clear();
  qPho3Pt_->clear();
  qPho4Pt_->clear();
  qPho5Pt_->clear();

  qPho0Eta_->clear();
  qPho1Eta_->clear();
  qPho2Eta_->clear();
  qPho3Eta_->clear();
  qPho4Eta_->clear();
  qPho5Eta_->clear();

  qPho0Phi_->clear();
  qPho1Phi_->clear();
  qPho2Phi_->clear();
  qPho3Phi_->clear();
  qPho4Phi_->clear();
  qPho5Phi_->clear();

  qPho0En_->clear();
  qPho1En_->clear();
  qPho2En_->clear();
  qPho3En_->clear();
  qPho4En_->clear();
  qPho5En_->clear();

  //ged qPhotons sorted variables (not all of them)
  qgedPho0Pt_->clear();
  qgedPho1Pt_->clear();
  qgedPho2Pt_->clear();
  qgedPho3Pt_->clear();
  qgedPho4Pt_->clear();
  qgedPho5Pt_->clear();
 
  qgedPho0Eta_->clear();
  qgedPho1Eta_->clear();
  qgedPho2Eta_->clear();
  qgedPho3Eta_->clear();
  qgedPho4Eta_->clear();
  qgedPho5Eta_->clear();

  qgedPho0Phi_->clear();
  qgedPho1Phi_->clear();
  qgedPho2Phi_->clear();
  qgedPho3Phi_->clear();
  qgedPho4Phi_->clear();
  qgedPho5Phi_->clear();

  qgedPho0En_->clear();
  qgedPho1En_->clear();
  qgedPho2En_->clear();
  qgedPho3En_->clear();
  qgedPho4En_->clear();
  qgedPho5En_->clear();


    //Muon sorted variables (not all of them)
  qMu0Pt_->clear();
  qMu1Pt_->clear();
  qMu2Pt_->clear();
  qMu3Pt_->clear();
  qMu4Pt_->clear();
  qMu5Pt_->clear();

  qMu0Eta_->clear();
  qMu1Eta_->clear();
  qMu2Eta_->clear();
  qMu3Eta_->clear();
  qMu4Eta_->clear();
  qMu5Eta_->clear();

  qMu0Phi_->clear();
  qMu1Phi_->clear();
  qMu2Phi_->clear();
  qMu3Phi_->clear();
  qMu4Phi_->clear();
  qMu5Phi_->clear();

  qMu0En_->clear();
  qMu1En_->clear();
  qMu2En_->clear();
  qMu3En_->clear();
  qMu4En_->clear();
  qMu5En_->clear();
  //Muon sorted cosmic variables (not all of them)
  qMuCosm0Pt_->clear();
  qMuCosm1Pt_->clear();
  qMuCosm2Pt_->clear();
  qMuCosm3Pt_->clear();
  qMuCosm4Pt_->clear();
  qMuCosm5Pt_->clear();

  qMuCosm0Eta_->clear();
  qMuCosm1Eta_->clear();
  qMuCosm2Eta_->clear();
  qMuCosm3Eta_->clear();
  qMuCosm4Eta_->clear();
  qMuCosm5Eta_->clear();

  qMuCosm0Phi_->clear();
  qMuCosm1Phi_->clear();
  qMuCosm2Phi_->clear();
  qMuCosm3Phi_->clear();
  qMuCosm4Phi_->clear();
  qMuCosm5Phi_->clear();

  qMuCosm0En_->clear();
  qMuCosm1En_->clear();
  qMuCosm2En_->clear();
  qMuCosm3En_->clear();
  qMuCosm4En_->clear();
  qMuCosm5En_->clear();


  //Muon Cosmic1Leg sorted variables (not all of them)
  qMuCosmLeg0Pt_->clear();
  qMuCosmLeg1Pt_->clear();
  qMuCosmLeg2Pt_->clear();
  qMuCosmLeg3Pt_->clear();
  qMuCosmLeg4Pt_->clear();
  qMuCosmLeg5Pt_->clear();

  qMuCosmLeg0Eta_->clear();
  qMuCosmLeg1Eta_->clear();
  qMuCosmLeg2Eta_->clear();
  qMuCosmLeg3Eta_->clear();
  qMuCosmLeg4Eta_->clear();
  qMuCosmLeg5Eta_->clear();

  qMuCosmLeg0Phi_->clear();
  qMuCosmLeg1Phi_->clear();
  qMuCosmLeg2Phi_->clear();
  qMuCosmLeg3Phi_->clear();
  qMuCosmLeg4Phi_->clear();
  qMuCosmLeg5Phi_->clear();

  qMuCosmLeg0En_->clear();
  qMuCosmLeg1En_->clear();
  qMuCosmLeg2En_->clear();
  qMuCosmLeg3En_->clear();
  qMuCosmLeg4En_->clear();
  qMuCosmLeg5En_->clear();

//___________________________All sorted qvariables end.



  qPFJet4CHSPt_->clear();
  qPFJet4CHSEta_->clear();
  qPFJet4CHSPhi_->clear();

  qPFJet8CHSPt_->clear();
  qPFJet8CHSEta_->clear();
  qPFJet8CHSPhi_->clear();
  
  qPFJetEIPt_->clear();
  qPFJetEIEta_->clear();
  qPFJetEIPhi_->clear();
  
  qPFJet8CHSSDPt_->clear();
  qPFJet8CHSSDEta_->clear();
  qPFJet8CHSSDPhi_->clear();

  qPFJetTopCHSPt_->clear();
  qPFJetTopCHSEta_->clear();
  qPFJetTopCHSPhi_->clear(); 

  qPFChMetPt_->clear();
  // qPFChMetEta_->clear();  
  qPFChMetPhi_->clear();
  qPFMetPt_->clear();
  // qPFMetEta_->clear();  
  qPFMetPhi_->clear();
  qNVtx_->clear();

  qCalJetPt_->clear();
  qCalJetEta_->clear();
  qCalJetPhi_->clear();
  qCalJetEn_->clear();
  qCalMETPt_->clear();
  // qCalMETEta_->clear();
  qCalMETPhi_->clear();
  qCalMETEn_->clear();
  qCalMETBEPt_->clear();
  // qCalMETBEEta_->clear();
  qCalMETBEPhi_->clear();
  qCalMETBEEn_->clear();
  qCalMETBEFOPt_->clear();
  // qCalMETBEFOEta_->clear();
  qCalMETBEFOPhi_->clear();
  qCalMETBEFOEn_->clear();
  qCalMETMPt_->clear();
  // qCalMETMEta_->clear();
  qCalMETMPhi_->clear();
  qCalMETMEn_->clear();
  qSCEn_ ->clear();
  qSCEta_->clear();
  qSCPhi_->clear();
  qSCEtaWidth_->clear();
  qSCPhiWidth_->clear();  
  qSCEnhfEM_ ->clear();
  qSCEtahfEM_->clear();
  qSCPhihfEM_->clear();
  // qSCEtaWidthhfEM_->clear();
  // qSCPhiWidthhfEM_->clear();  
  qSCEn5x5_ ->clear();
  qSCEta5x5_->clear();
  qSCPhi5x5_->clear();
  qSCEtaWidth5x5_->clear();
  qSCPhiWidth5x5_->clear();  
  qCCEn_ ->clear();
  qCCEta_->clear();
  qCCPhi_->clear();
  qCCEn5x5_ ->clear();
  qCCEta5x5_->clear();
  qCCPhi5x5_->clear();

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
  qMuChi2_->clear();
  qMuCosmPt_->clear();
  qMuCosmEta_->clear();
  qMuCosmPhi_->clear();
  qMuCosmEn_->clear();
  qMuCosmCh_->clear();
  qMuCosmChi2_->clear();  
  qMuCosmLegPt_->clear();
  qMuCosmLegEta_->clear();
  qMuCosmLegPhi_->clear();
  qMuCosmLegEn_->clear();
  qMuCosmLegCh_->clear();
  qMuCosmLegChi2_->clear();    

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
  qEBiEta_->clear();
  qEBiPhi_->clear();

  qEEenergy_->clear();
  qEEtime_->clear();
  qEEchi2_->clear();
  qEEix_->clear();
  qEEiy_->clear();

  qESenergy_->clear();
  qEStime_->clear();
  // qESchi2_->clear();
  qESix_->clear();
  qESiy_->clear();

  qHBHEenergy_->clear();
  qHBHEtime_->clear();
  qHBHEauxe_->clear();
  qHBHEieta_->clear();
  qHBHEiphi_->clear();

  qHFenergy_->clear();
  qHFtime_->clear();
  // qHFchi2_->clear();
  qHFieta_->clear();
  qHFiphi_->clear();

  // qHOenergy_->clear();
  // qHOtime_->clear();
  // // qHOchi2_->clear();
  // qHOieta_->clear();
  // qHOiphi_->clear();

  qPreShEn_->clear();
  // qPreShCorrEn_->clear();
  qPreShEta_->clear();
  qPreShPhi_->clear();
  qPreShYEn_->clear();
  // qPreShYCorrEn_->clear();
  qPreShYEta_->clear();
  qPreShYPhi_->clear();

  // qCTPt_->clear();
  // qCTEta_->clear();
  // qCTPhi_->clear();

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
    // if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
    //   {
	
	if(type.compare(std::string("PF")) == 0)
	  {
	    PFJetPt_->push_back(i->pt());
	    PFJetEta_->push_back(i->eta());
	    PFJetPhi_->push_back(i->phi()); 
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
      // std::cout << "ele Jet pt: " << i->pt()   << std::endl; //TODO
      //       std::cout << "ele Jet eta: " << i->eta()   << std::endl; //TODO
      //       std::cout << "ele Jet phi: " << i->phi()   << std::endl; //TODO

	  }
	 


  }
        // std::cout << "end of Jets collision "  << std::endl; //TODO
  PFJet0Pt_->push_back(jets->begin()->pt()); 
          // std::cout << "ele Jet 0pt: " << (jets->begin()->pt())   << std::endl;

  PFJet1Pt_->push_back((jets->begin()+jets->size()/5)->pt()); 
        // std::cout << "ele Jet 1pt: " << ((jets->begin()+jets->size()/5)->pt())    << std::endl;

  PFJet2Pt_->push_back((jets->begin()+jets->size()*2/5)->pt());
        // std::cout << "ele Jet 2pt: " << ((jets->begin()+jets->size()*2/5)->pt())  << std::endl;
  
  PFJet3Pt_->push_back((jets->begin()+jets->size()*3/5)->pt());
        // std::cout << "ele Jet 3pt: " << ((jets->begin()+jets->size()*3/5)->pt())  << std::endl;

  PFJet4Pt_->push_back((jets->begin()+jets->size()*4/5)->pt());
        // std::cout << "ele Jet 4pt: " << ((jets->begin()+jets->size()*4/5)->pt())  << std::endl; 

  PFJet5Pt_->push_back((jets->begin()+jets->size()-1)->pt());
       // std::cout << "ele Jet 5pt: " << ((jets->begin()+jets->size()-1)->pt())   << std::endl;

  PFJet0Eta_->push_back(jets->begin()->eta()); 
          // std::cout << "ele Jet 0eta: " << (jets->begin()->eta())   << std::endl;

  PFJet1Eta_->push_back((jets->begin()+jets->size()/5)->eta()); 
        // std::cout << "ele Jet 1eta: " << ((jets->begin()+jets->size()/5)->eta())    << std::endl;

  PFJet2Eta_->push_back((jets->begin()+jets->size()*2/5)->eta());
        // std::cout << "ele Jet 2eta: " << ((jets->begin()+jets->size()*2/5)->eta())  << std::endl;
  
  PFJet3Eta_->push_back((jets->begin()+jets->size()*3/5)->eta());
        // std::cout << "ele Jet 3eta: " << ((jets->begin()+jets->size()*3/5)->eta())  << std::endl;

  PFJet4Eta_->push_back((jets->begin()+jets->size()*4/5)->eta());
        // std::cout << "ele Jet 4eta: " << ((jets->begin()+jets->size()*4/5)->eta())  << std::endl; 

  PFJet5Eta_->push_back((jets->begin()+jets->size()-1)->eta());
       // std::cout << "ele Jet 5eta: " << ((jets->begin()+jets->size()-1)->eta())   << std::endl; 

  PFJet0Phi_->push_back(jets->begin()->phi()); 
          // std::cout << "ele Jet 0Phi: " << (jets->begin()->phi())   << std::endl;

  PFJet1Phi_->push_back((jets->begin()+jets->size()/5)->phi()); 
        // std::cout << "ele Jet 1Phi: " << ((jets->begin()+jets->size()/5)->phi())    << std::endl;

  PFJet2Phi_->push_back((jets->begin()+jets->size()*2/5)->phi());
        // std::cout << "ele Jet 2Phi: " << ((jets->begin()+jets->size()*2/5)->phi())  << std::endl;
  
  PFJet3Phi_->push_back((jets->begin()+jets->size()*3/5)->phi());
        // std::cout << "ele Jet 3Phi: " << ((jets->begin()+jets->size()*3/5)->phi())  << std::endl;

  PFJet4Phi_->push_back((jets->begin()+jets->size()*4/5)->phi());
        // std::cout << "ele Jet 4Phi: " << ((jets->begin()+jets->size()*4/5)->phi())  << std::endl; 

  PFJet5Phi_->push_back((jets->begin()+jets->size()-1)->phi());
       // std::cout << "ele Jet 5Phi: " << ((jets->begin()+jets->size()-1)->phi())   << std::endl;             
     
  return;
}




template<typename PFJet4CHSCollection>
void AODAnalyzer::fill4CHSJets(const edm::Handle<PFJet4CHSCollection> & jets4CHS, std::string type)
{
  // Selected jets
  //reco::CaloJetCollection recojets;
  typename PFJet4CHSCollection::const_iterator i = jets4CHS->begin();
  for(;i != jets4CHS->end(); i++){
    // if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
    //   {
  
  if(type.compare(std::string("PF")) == 0)
    {
      PFJet4CHSPt_->push_back(i->pt());
      PFJet4CHSEta_->push_back(i->eta());
      PFJet4CHSPhi_->push_back(i->phi());
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
            // std::cout << "ele Jet pt: " << i->pt()   << std::endl; //TODO
            // std::cout << "ele Jet eta: " << i->eta()   << std::endl; //TODO
            // std::cout << "ele Jet phi: " << i->phi()   << std::endl; //TODO
    }
  
      // }
  }

        // std::cout << "end of Jets collision "  << std::endl; //TODO
  PFJet4CHS0Pt_->push_back(jets4CHS->begin()->pt()); 
          // std::cout << "ele Jet 0pt: " << (jets4CHS->begin()->pt())   << std::endl;

  PFJet4CHS1Pt_->push_back((jets4CHS->begin()+jets4CHS->size()/5)->pt()); 
        // std::cout << "ele Jet 1pt: " << ((jets4CHS->begin()+jets4CHS->size()/5)->pt())    << std::endl;

  PFJet4CHS2Pt_->push_back((jets4CHS->begin()+jets4CHS->size()*2/5)->pt());
        // std::cout << "ele Jet 2pt: " << ((jets4CHS->begin()+jets4CHS->size()*2/5)->pt())  << std::endl;
  
  PFJet4CHS3Pt_->push_back((jets4CHS->begin()+jets4CHS->size()*3/5)->pt());
        // std::cout << "ele Jet 3pt: " << ((jets4CHS->begin()+jets4CHS->size()*3/5)->pt())  << std::endl;

  PFJet4CHS4Pt_->push_back((jets4CHS->begin()+jets4CHS->size()*4/5)->pt());
        // std::cout << "ele Jet 4pt: " << ((jets4CHS->begin()+jets4CHS->size()*4/5)->pt())  << std::endl; 

  PFJet4CHS5Pt_->push_back((jets4CHS->begin()+jets4CHS->size()-1)->pt());
       // std::cout << "ele Jet 5pt: " << ((jets4CHS->begin()+jets4CHS->size()-1)->pt())   << std::endl;

  PFJet4CHS0Eta_->push_back(jets4CHS->begin()->eta()); 
          // std::cout << "ele Jet 0eta: " << (jets4CHS->begin()->eta())   << std::endl;

  PFJet4CHS1Eta_->push_back((jets4CHS->begin()+jets4CHS->size()/5)->eta()); 
        // std::cout << "ele Jet 1eta: " << ((jets4CHS->begin()+jets4CHS->size()/5)->eta())    << std::endl;

  PFJet4CHS2Eta_->push_back((jets4CHS->begin()+jets4CHS->size()*2/5)->eta());
        // std::cout << "ele Jet 2eta: " << ((jets4CHS->begin()+jets4CHS->size()*2/5)->eta())  << std::endl;
  
  PFJet4CHS3Eta_->push_back((jets4CHS->begin()+jets4CHS->size()*3/5)->eta());
        // std::cout << "ele Jet 3eta: " << ((jets4CHS->begin()+jets4CHS->size()*3/5)->eta())  << std::endl;

  PFJet4CHS4Eta_->push_back((jets4CHS->begin()+jets4CHS->size()*4/5)->eta());
        // std::cout << "ele Jet 4eta: " << ((jets4CHS->begin()+jets4CHS->size()*4/5)->eta())  << std::endl; 

  PFJet4CHS5Eta_->push_back((jets4CHS->begin()+jets4CHS->size()-1)->eta());
       // std::cout << "ele Jet 5eta: " << ((jets4CHS->begin()+jets4CHS->size()-1)->eta())   << std::endl; 

  PFJet4CHS0Phi_->push_back(jets4CHS->begin()->phi()); 
          // std::cout << "ele Jet 0Phi: " << (jets4CHS->begin()->phi())   << std::endl;

  PFJet4CHS1Phi_->push_back((jets4CHS->begin()+jets4CHS->size()/5)->phi()); 
        // std::cout << "ele Jet 1Phi: " << ((jets4CHS->begin()+jets4CHS->size()/5)->phi())    << std::endl;

  PFJet4CHS2Phi_->push_back((jets4CHS->begin()+jets4CHS->size()*2/5)->phi());
        // std::cout << "ele Jet 2Phi: " << ((jets4CHS->begin()+jets4CHS->size()*2/5)->phi())  << std::endl;
  
  PFJet4CHS3Phi_->push_back((jets4CHS->begin()+jets4CHS->size()*3/5)->phi());
        // std::cout << "ele Jet 3Phi: " << ((jets4CHS->begin()+jets4CHS->size()*3/5)->phi())  << std::endl;

  PFJet4CHS4Phi_->push_back((jets4CHS->begin()+jets4CHS->size()*4/5)->phi());
        // std::cout << "ele Jet 4Phi: " << ((jets4CHS->begin()+jets4CHS->size()*4/5)->phi())  << std::endl; 

  PFJet4CHS5Phi_->push_back((jets4CHS->begin()+jets4CHS->size()-1)->phi());
       // std::cout << "ele Jet 5Phi: " << ((jets4CHS->begin()+jets4CHS->size()-1)->phi())   << std::endl; 



  return;
}

template<typename PFJet8CHSCollection>
void AODAnalyzer::fill8CHSJets(const edm::Handle<PFJet8CHSCollection> & jets8CHS, std::string type)
{
  // Selected jets
  //reco::CaloJetCollection recojets;
  typename PFJet8CHSCollection::const_iterator i = jets8CHS->begin();
  for(;i != jets8CHS->end(); i++){
    // if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
    //   {
  
  if(type.compare(std::string("PF")) == 0)
    {
      PFJet8CHSPt_->push_back(i->pt());
      PFJet8CHSEta_->push_back(i->eta());
      PFJet8CHSPhi_->push_back(i->phi());
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
    }
  
      // }
  }

       // std::cout << "end of Jets collision "  << std::endl; //TODO
  PFJet8CHS0Pt_->push_back(jets8CHS->begin()->pt()); 
          // std::cout << "ele Jet 0pt: " << (jets8CHS->begin()->pt())   << std::endl;

  PFJet8CHS1Pt_->push_back((jets8CHS->begin()+jets8CHS->size()/5)->pt()); 
        // std::cout << "ele Jet 1pt: " << ((jets8CHS->begin()+jets8CHS->size()/5)->pt())    << std::endl;

  PFJet8CHS2Pt_->push_back((jets8CHS->begin()+jets8CHS->size()*2/5)->pt());
        // std::cout << "ele Jet 2pt: " << ((jets8CHS->begin()+jets8CHS->size()*2/5)->pt())  << std::endl;
  
  PFJet8CHS3Pt_->push_back((jets8CHS->begin()+jets8CHS->size()*3/5)->pt());
        // std::cout << "ele Jet 3pt: " << ((jets8CHS->begin()+jets8CHS->size()*3/5)->pt())  << std::endl;

  PFJet8CHS4Pt_->push_back((jets8CHS->begin()+jets8CHS->size()*4/5)->pt());
        // std::cout << "ele Jet 4pt: " << ((jets8CHS->begin()+jets8CHS->size()*4/5)->pt())  << std::endl; 

  PFJet8CHS5Pt_->push_back((jets8CHS->begin()+jets8CHS->size()-1)->pt());
       // std::cout << "ele Jet 5pt: " << ((jets8CHS->begin()+jets8CHS->size()-1)->pt())   << std::endl;

  PFJet8CHS0Eta_->push_back(jets8CHS->begin()->eta()); 
          // std::cout << "ele Jet 0eta: " << (jets8CHS->begin()->eta())   << std::endl;

  PFJet8CHS1Eta_->push_back((jets8CHS->begin()+jets8CHS->size()/5)->eta()); 
        // std::cout << "ele Jet 1eta: " << ((jets8CHS->begin()+jets8CHS->size()/5)->eta())    << std::endl;

  PFJet8CHS2Eta_->push_back((jets8CHS->begin()+jets8CHS->size()*2/5)->eta());
        // std::cout << "ele Jet 2eta: " << ((jets8CHS->begin()+jets8CHS->size()*2/5)->eta())  << std::endl;
  
  PFJet8CHS3Eta_->push_back((jets8CHS->begin()+jets8CHS->size()*3/5)->eta());
        // std::cout << "ele Jet 3eta: " << ((jets8CHS->begin()+jets8CHS->size()*3/5)->eta())  << std::endl;

  PFJet8CHS4Eta_->push_back((jets8CHS->begin()+jets8CHS->size()*4/5)->eta());
        // std::cout << "ele Jet 4eta: " << ((jets8CHS->begin()+jets8CHS->size()*4/5)->eta())  << std::endl; 

  PFJet8CHS5Eta_->push_back((jets8CHS->begin()+jets8CHS->size()-1)->eta());
       // std::cout << "ele Jet 5eta: " << ((jets8CHS->begin()+jets8CHS->size()-1)->eta())   << std::endl; 

  PFJet8CHS0Phi_->push_back(jets8CHS->begin()->phi()); 
          // std::cout << "ele Jet 0Phi: " << (jets8CHS->begin()->phi())   << std::endl;

  PFJet8CHS1Phi_->push_back((jets8CHS->begin()+jets8CHS->size()/5)->phi()); 
        // std::cout << "ele Jet 1Phi: " << ((jets8CHS->begin()+jets8CHS->size()/5)->phi())    << std::endl;

  PFJet8CHS2Phi_->push_back((jets8CHS->begin()+jets8CHS->size()*2/5)->phi());
        // std::cout << "ele Jet 2Phi: " << ((jets8CHS->begin()+jets8CHS->size()*2/5)->phi())  << std::endl;
  
  PFJet8CHS3Phi_->push_back((jets8CHS->begin()+jets8CHS->size()*3/5)->phi());
        // std::cout << "ele Jet 3Phi: " << ((jets8CHS->begin()+jets8CHS->size()*3/5)->phi())  << std::endl;

  PFJet8CHS4Phi_->push_back((jets8CHS->begin()+jets8CHS->size()*4/5)->phi());
        // std::cout << "ele Jet 4Phi: " << ((jets8CHS->begin()+jets8CHS->size()*4/5)->phi())  << std::endl; 

  PFJet8CHS5Phi_->push_back((jets8CHS->begin()+jets8CHS->size()-1)->phi());
       // std::cout << "ele Jet 5Phi: " << ((jets8CHS->begin()+jets8CHS->size()-1)->phi())   << std::endl; 


  return;
}

template<typename PFJetEICollection>
void AODAnalyzer::fillEIJets(const edm::Handle<PFJetEICollection> & jetsEI, std::string type)
{
  // Selected jets
  //reco::CaloJetCollection recojets;
  typename PFJetEICollection::const_iterator i = jetsEI->begin();
  for(;i != jetsEI->end(); i++){
    // if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
    //   {
  
  if(type.compare(std::string("PF")) == 0)
    {
      PFJetEIPt_->push_back(i->pt());
      PFJetEIEta_->push_back(i->eta());
      PFJetEIPhi_->push_back(i->phi());
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
    }
  
      // }
  }
     // std::cout << "end of Jets collision "  << std::endl; //TODO
  PFJetEI0Pt_->push_back(jetsEI->begin()->pt()); 
          // std::cout << "ele Jet 0pt: " << (jetsEI->begin()->pt())   << std::endl;

  PFJetEI1Pt_->push_back((jetsEI->begin()+jetsEI->size()/5)->pt()); 
        // std::cout << "ele Jet 1pt: " << ((jetsEI->begin()+jetsEI->size()/5)->pt())    << std::endl;

  PFJetEI2Pt_->push_back((jetsEI->begin()+jetsEI->size()*2/5)->pt());
        // std::cout << "ele Jet 2pt: " << ((jetsEI->begin()+jetsEI->size()*2/5)->pt())  << std::endl;
  
  PFJetEI3Pt_->push_back((jetsEI->begin()+jetsEI->size()*3/5)->pt());
        // std::cout << "ele Jet 3pt: " << ((jetsEI->begin()+jetsEI->size()*3/5)->pt())  << std::endl;

  PFJetEI4Pt_->push_back((jetsEI->begin()+jetsEI->size()*4/5)->pt());
        // std::cout << "ele Jet 4pt: " << ((jetsEI->begin()+jetsEI->size()*4/5)->pt())  << std::endl; 

  PFJetEI5Pt_->push_back((jetsEI->begin()+jetsEI->size()-1)->pt());
       // std::cout << "ele Jet 5pt: " << ((jetsEI->begin()+jetsEI->size()-1)->pt())   << std::endl;

  PFJetEI0Eta_->push_back(jetsEI->begin()->eta()); 
          // std::cout << "ele Jet 0eta: " << (jetsEI->begin()->eta())   << std::endl;

  PFJetEI1Eta_->push_back((jetsEI->begin()+jetsEI->size()/5)->eta()); 
        // std::cout << "ele Jet 1eta: " << ((jetsEI->begin()+jetsEI->size()/5)->eta())    << std::endl;

  PFJetEI2Eta_->push_back((jetsEI->begin()+jetsEI->size()*2/5)->eta());
        // std::cout << "ele Jet 2eta: " << ((jetsEI->begin()+jetsEI->size()*2/5)->eta())  << std::endl;
  
  PFJetEI3Eta_->push_back((jetsEI->begin()+jetsEI->size()*3/5)->eta());
        // std::cout << "ele Jet 3eta: " << ((jetsEI->begin()+jetsEI->size()*3/5)->eta())  << std::endl;

  PFJetEI4Eta_->push_back((jetsEI->begin()+jetsEI->size()*4/5)->eta());
        // std::cout << "ele Jet 4eta: " << ((jetsEI->begin()+jetsEI->size()*4/5)->eta())  << std::endl; 

  PFJetEI5Eta_->push_back((jetsEI->begin()+jetsEI->size()-1)->eta());
       // std::cout << "ele Jet 5eta: " << ((jetsEI->begin()+jetsEI->size()-1)->eta())   << std::endl; 

  PFJetEI0Phi_->push_back(jetsEI->begin()->phi()); 
          // std::cout << "ele Jet 0Phi: " << (jetsEI->begin()->phi())   << std::endl;

  PFJetEI1Phi_->push_back((jetsEI->begin()+jetsEI->size()/5)->phi()); 
        // std::cout << "ele Jet 1Phi: " << ((jetsEI->begin()+jetsEI->size()/5)->phi())    << std::endl;

  PFJetEI2Phi_->push_back((jetsEI->begin()+jetsEI->size()*2/5)->phi());
        // std::cout << "ele Jet 2Phi: " << ((jetsEI->begin()+jetsEI->size()*2/5)->phi())  << std::endl;
  
  PFJetEI3Phi_->push_back((jetsEI->begin()+jetsEI->size()*3/5)->phi());
        // std::cout << "ele Jet 3Phi: " << ((jetsEI->begin()+jetsEI->size()*3/5)->phi())  << std::endl;

  PFJetEI4Phi_->push_back((jetsEI->begin()+jetsEI->size()*4/5)->phi());
        // std::cout << "ele Jet 4Phi: " << ((jetsEI->begin()+jetsEI->size()*4/5)->phi())  << std::endl; 

  PFJetEI5Phi_->push_back((jetsEI->begin()+jetsEI->size()-1)->phi());
       // std::cout << "ele Jet 5Phi: " << ((jetsEI->begin()+jetsEI->size()-1)->phi())   << std::endl; 


  return;
}

template<typename PFJet8CHSSoftDropCollection>
void AODAnalyzer::fill8CHSoftDropJets(const edm::Handle<PFJet8CHSSoftDropCollection> & jets8CHSSoftDrop, std::string type)
{
  // Selected jets
  //reco::CaloJetCollection recojets;
  typename PFJet8CHSSoftDropCollection::const_iterator i = jets8CHSSoftDrop->begin();
  for(;i != jets8CHSSoftDrop->end(); i++){
    // if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
    //   {
  
  if(type.compare(std::string("PF")) == 0)
    {
      PFJet8CHSSDPt_->push_back(i->pt());
      PFJet8CHSSDEta_->push_back(i->eta());
      PFJet8CHSSDPhi_->push_back(i->phi());
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
    }
  
      // }
  }
     // std::cout << "end of Jets collision "  << std::endl; //TODO
  // PFJet8CHSSD0Pt_->push_back(jets8CHSSoftDrop->begin()->pt()); 
  // std::cout << "ele Jet 0pt: " << (jets8CHSSoftDrop->begin()->pt())   << std::endl;

  // PFJet8CHSSD1Pt_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()/5)->pt()); 
         // std::cout << "ele Jet 1pt: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()/5)->pt())    << std::endl;

  // PFJet8CHSSD2Pt_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*2/5)->pt());
  //       // std::cout << "ele Jet 2pt: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*2/5)->pt())  << std::endl;
  
  // PFJet8CHSSD3Pt_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*3/5)->pt());
  //       // std::cout << "ele Jet 3pt: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*3/5)->pt())  << std::endl;

  // PFJet8CHSSD4Pt_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*4/5)->pt());
  //       // std::cout << "ele Jet 4pt: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*4/5)->pt())  << std::endl; 

  // PFJet8CHSSD5Pt_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()-1)->pt());
  //      // std::cout << "ele Jet 5pt: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()-1)->pt())   << std::endl;

  // PFJet8CHSSD0Eta_->push_back(jets8CHSSoftDrop->begin()->eta()); 
  //         // std::cout << "ele Jet 0eta: " << (jets8CHSSoftDrop->begin()->eta())   << std::endl;

  // PFJet8CHSSD1Eta_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()/5)->eta()); 
  //       // std::cout << "ele Jet 1eta: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()/5)->eta())    << std::endl;

  // PFJet8CHSSD2Eta_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*2/5)->eta());
  //       // std::cout << "ele Jet 2eta: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*2/5)->eta())  << std::endl;
  
  // PFJet8CHSSD3Eta_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*3/5)->eta());
  //       // std::cout << "ele Jet 3eta: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*3/5)->eta())  << std::endl;

  // PFJet8CHSSD4Eta_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*4/5)->eta());
  //       // std::cout << "ele Jet 4eta: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*4/5)->eta())  << std::endl; 

  // PFJet8CHSSD5Eta_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()-1)->eta());
  //      // std::cout << "ele Jet 5eta: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()-1)->eta())   << std::endl; 

  // PFJet8CHSSD0Phi_->push_back(jets8CHSSoftDrop->begin()->phi()); 
  //         // std::cout << "ele Jet 0Phi: " << (jets8CHSSoftDrop->begin()->phi())   << std::endl;

  // PFJet8CHSSD1Phi_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()/5)->phi()); 
  //       // std::cout << "ele Jet 1Phi: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()/5)->phi())    << std::endl;

  // PFJet8CHSSD2Phi_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*2/5)->phi());
  //       // std::cout << "ele Jet 2Phi: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*2/5)->phi())  << std::endl;
  
  // PFJet8CHSSD3Phi_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*3/5)->phi());
  //       // std::cout << "ele Jet 3Phi: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*3/5)->phi())  << std::endl;

  // PFJet8CHSSD4Phi_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*4/5)->phi());
  //       // std::cout << "ele Jet 4Phi: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()*4/5)->phi())  << std::endl; 

  // PFJet8CHSSD5Phi_->push_back((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()-1)->phi());
  //      // std::cout << "ele Jet 5Phi: " << ((jets8CHSSoftDrop->begin()+jets8CHSSoftDrop->size()-1)->phi())   << std::endl; 

  return;
}

template<typename PFJetTopCHSCollection>
void AODAnalyzer::fillTopCHSJets(const edm::Handle<PFJetTopCHSCollection> & jetsTopCHS, std::string type)
{
  // Selected jets
  //reco::CaloJetCollection recojets;
  typename PFJetTopCHSCollection::const_iterator i = jetsTopCHS->begin();
  for(;i != jetsTopCHS->end(); i++){
    // if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
    //   {
  
  if(type.compare(std::string("PF")) == 0)
    {
      PFJetTopCHSPt_->push_back(i->pt());
      PFJetTopCHSEta_->push_back(i->eta());
      PFJetTopCHSPhi_->push_back(i->phi());
      // std::cout << "ele pt: " << i->pt() << std::endl; //TEST -- works
      // std::cout << "ele eta: " << i->eta() << std::endl;   //TEST --works
      // std::cout << "ele phi: " << i->phi() << std::endl;   //TEST --works
    }
  
      // }
  }

  //    // std::cout << "end of Jets collision "  << std::endl; //TODO
  // PFJetTopCHS0Pt_->push_back(jetsTopCHS->begin()->pt()); 
  //         // std::cout << "ele Jet 0pt: " << (jetsTopCHS->begin()->pt())   << std::endl;

  // PFJetTopCHS1Pt_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()/5)->pt()); 
  //       // std::cout << "ele Jet 1pt: " << ((jetsTopCHS->begin()+jetsTopCHS->size()/5)->pt())    << std::endl;

  // PFJetTopCHS2Pt_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*2/5)->pt());
  //       // std::cout << "ele Jet 2pt: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*2/5)->pt())  << std::endl;
  
  // PFJetTopCHS3Pt_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*3/5)->pt());
  //       // std::cout << "ele Jet 3pt: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*3/5)->pt())  << std::endl;

  // PFJetTopCHS4Pt_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*4/5)->pt());
  //       // std::cout << "ele Jet 4pt: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*4/5)->pt())  << std::endl; 

  // PFJetTopCHS5Pt_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()-1)->pt());
  //      // std::cout << "ele Jet 5pt: " << ((jetsTopCHS->begin()+jetsTopCHS->size()-1)->pt())   << std::endl;

  // PFJetTopCHS0Eta_->push_back(jetsTopCHS->begin()->eta()); 
  //         // std::cout << "ele Jet 0eta: " << (jetsTopCHS->begin()->eta())   << std::endl;

  // PFJetTopCHS1Eta_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()/5)->eta()); 
  //       // std::cout << "ele Jet 1eta: " << ((jetsTopCHS->begin()+jetsTopCHS->size()/5)->eta())    << std::endl;

  // PFJetTopCHS2Eta_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*2/5)->eta());
  //       // std::cout << "ele Jet 2eta: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*2/5)->eta())  << std::endl;
  
  // PFJetTopCHS3Eta_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*3/5)->eta());
  //       // std::cout << "ele Jet 3eta: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*3/5)->eta())  << std::endl;

  // PFJetTopCHS4Eta_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*4/5)->eta());
  //       // std::cout << "ele Jet 4eta: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*4/5)->eta())  << std::endl; 

  // PFJetTopCHS5Eta_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()-1)->eta());
  //      // std::cout << "ele Jet 5eta: " << ((jetsTopCHS->begin()+jetsTopCHS->size()-1)->eta())   << std::endl; 

  // PFJetTopCHS0Phi_->push_back(jetsTopCHS->begin()->phi()); 
  //         // std::cout << "ele Jet 0Phi: " << (jetsTopCHS->begin()->phi())   << std::endl;

  // PFJetTopCHS1Phi_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()/5)->phi()); 
  //       // std::cout << "ele Jet 1Phi: " << ((jetsTopCHS->begin()+jetsTopCHS->size()/5)->phi())    << std::endl;

  // PFJetTopCHS2Phi_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*2/5)->phi());
  //       // std::cout << "ele Jet 2Phi: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*2/5)->phi())  << std::endl;
  
  // PFJetTopCHS3Phi_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*3/5)->phi());
  //       // std::cout << "ele Jet 3Phi: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*3/5)->phi())  << std::endl;

  // PFJetTopCHS4Phi_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()*4/5)->phi());
  //       // std::cout << "ele Jet 4Phi: " << ((jetsTopCHS->begin()+jetsTopCHS->size()*4/5)->phi())  << std::endl; 

  // PFJetTopCHS5Phi_->push_back((jetsTopCHS->begin()+jetsTopCHS->size()-1)->phi());
  //      // std::cout << "ele Jet 5Phi: " << ((jetsTopCHS->begin()+jetsTopCHS->size()-1)->phi())   << std::endl; 
  return;
}

template<typename PFChMETCollection>
void AODAnalyzer::fillPFChMets(const edm::Handle<PFChMETCollection> & pfchmets)
{
  // std::cout << "fillPFChMets is being called!" << std::endl;
  // std::cout << pfchmets->size() <<std::endl;
  typename PFChMETCollection::const_iterator i = pfchmets->begin();
  for(;i != pfchmets->end(); i++){
    PFChMetPt_->push_back(i->pt());
    // PFChMetEta_->push_back(i->eta());
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
    PFMetPt_->push_back(i->pt());
    // PFMetEta_->push_back(i->eta());
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
        CalJetPt_->push_back(i->pt());
        CalJetEta_->push_back(i->eta());
        CalJetPhi_->push_back(i->phi());
        CalJetEn_->push_back(i->energy());
      // std::cout << "ele CaloJets pt: " << i->pt()   << std::endl; //TODO
      // std::cout << "ele CaloJets eta: " << i->eta()   << std::endl; //TODO
      // std::cout << "ele CaloJets phi: " << i->phi()   << std::endl; //TODO
      // std::cout << "ele CaloJets energy: " << i->energy()   << std::endl; //TODO


  }
     // sort muons by pt
       // std::cout << "end of calojets collision "  << std::endl; //TODO
  CalJet0Pt_->push_back(calojets->begin()->pt()); 
          // std::cout << "ele Jet 0pt: " << (calojets->begin()->pt())   << std::endl;

  CalJet1Pt_->push_back((calojets->begin()+calojets->size()/5)->pt()); 
        // std::cout << "ele Jet 1pt: " << ((calojets->begin()+calojets->size()/5)->pt())    << std::endl;

  CalJet2Pt_->push_back((calojets->begin()+calojets->size()*2/5)->pt());
        // std::cout << "ele Jet 2pt: " << ((calojets->begin()+calojets->size()*2/5)->pt())  << std::endl;
  
  CalJet3Pt_->push_back((calojets->begin()+calojets->size()*3/5)->pt());
        // std::cout << "ele Jet 3pt: " << ((calojets->begin()+calojets->size()*3/5)->pt())  << std::endl;

  CalJet4Pt_->push_back((calojets->begin()+calojets->size()*4/5)->pt());
        // std::cout << "ele Jet 4pt: " << ((calojets->begin()+calojets->size()*4/5)->pt())  << std::endl; 

  CalJet5Pt_->push_back((calojets->begin()+calojets->size()-1)->pt());
       // std::cout << "ele Jet 5pt: " << ((calojets->begin()+calojets->size()-1)->pt())   << std::endl;

  CalJet0Eta_->push_back(calojets->begin()->eta()); 
          // std::cout << "ele Jet 0eta: " << (calojets->begin()->eta())   << std::endl;

  CalJet1Eta_->push_back((calojets->begin()+calojets->size()/5)->eta()); 
        // std::cout << "ele Jet 1eta: " << ((calojets->begin()+calojets->size()/5)->eta())    << std::endl;

  CalJet2Eta_->push_back((calojets->begin()+calojets->size()*2/5)->eta());
        // std::cout << "ele Jet 2eta: " << ((calojets->begin()+calojets->size()*2/5)->eta())  << std::endl;
  
  CalJet3Eta_->push_back((calojets->begin()+calojets->size()*3/5)->eta());
        // std::cout << "ele Jet 3eta: " << ((calojets->begin()+calojets->size()*3/5)->eta())  << std::endl;

  CalJet4Eta_->push_back((calojets->begin()+calojets->size()*4/5)->eta());
        // std::cout << "ele Jet 4eta: " << ((calojets->begin()+calojets->size()*4/5)->eta())  << std::endl; 

  CalJet5Eta_->push_back((calojets->begin()+calojets->size()-1)->eta());
       // std::cout << "ele Jet 5eta: " << ((calojets->begin()+calojets->size()-1)->eta())   << std::endl; 

  CalJet0Phi_->push_back(calojets->begin()->phi()); 
        // std::cout << "ele Jet 0Phi: " << (calojets->begin()->phi())   << std::endl;

  CalJet1Phi_->push_back((calojets->begin()+calojets->size()/5)->phi()); 
        // std::cout << "ele Jet 1Phi: " << ((calojets->begin()+calojets->size()/5)->phi())    << std::endl;

  CalJet2Phi_->push_back((calojets->begin()+calojets->size()*2/5)->phi());
        // std::cout << "ele Jet 2Phi: " << ((calojets->begin()+calojets->size()*2/5)->phi())  << std::endl;
  
  CalJet3Phi_->push_back((calojets->begin()+calojets->size()*3/5)->phi());
        // std::cout << "ele Jet 3Phi: " << ((calojets->begin()+calojets->size()*3/5)->phi())  << std::endl;

  CalJet4Phi_->push_back((calojets->begin()+calojets->size()*4/5)->phi());
        // std::cout << "ele Jet 4Phi: " << ((calojets->begin()+calojets->size()*4/5)->phi())  << std::endl; 

  CalJet5Phi_->push_back((calojets->begin()+calojets->size()-1)->phi());
       // std::cout << "ele Jet 5Phi: " << ((calojets->begin()+calojets->size()-1)->phi())   << std::endl; 

  CalJet0En_->push_back(calojets->begin()->energy()); 
          // std::cout << "ele Jet 0En: " << (calojets->begin()->energy())   << std::endl;

  CalJet1En_->push_back((calojets->begin()+calojets->size()/5)->energy()); 
        // std::cout << "ele Jet 1En: " << ((calojets->begin()+calojets->size()/5)->energy())    << std::endl;

  CalJet2En_->push_back((calojets->begin()+calojets->size()*2/5)->energy());
        // std::cout << "ele Jet 2En: " << ((calojets->begin()+calojets->size()*2/5)->energy())  << std::endl;
  
  CalJet3En_->push_back((calojets->begin()+calojets->size()*3/5)->energy());
        // std::cout << "ele Jet 3En: " << ((calojets->begin()+calojets->size()*3/5)->energy())  << std::endl;

  CalJet4En_->push_back((calojets->begin()+calojets->size()*4/5)->energy());
        // std::cout << "ele Jet 4En: " << ((calojets->begin()+calojets->size()*4/5)->energy())  << std::endl; 

  CalJet5En_->push_back((calojets->begin()+calojets->size()-1)->energy());
       // std::cout << "ele Jet 5En: " << ((calojets->begin()+calojets->size()-1)->energy())   << std::endl;             

  return;

}

template<typename CaloMETCollection>
void AODAnalyzer::fillCaloMETs(const edm::Handle<CaloMETCollection> & caloMETs)
{
  // std::cout << "fills is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename CaloMETCollection::const_iterator i = caloMETs->begin();
  for(;i != caloMETs->end(); i++){
        CalMETPt_->push_back(i->pt());
        // CalMETEta_->push_back(i->eta());
        CalMETPhi_->push_back(i->phi());
        CalMETEn_->push_back(i->energy());

  }
  return;

}


template<typename CaloMETBECollection>
void AODAnalyzer::fillCaloMETBEs(const edm::Handle<CaloMETBECollection> & caloMETBEs)
{
  // std::cout << "fills is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename CaloMETBECollection::const_iterator i = caloMETBEs->begin();
  for(;i != caloMETBEs->end(); i++){
        CalMETBEPt_->push_back(i->pt());
        // CalMETBEEta_->push_back(i->eta());
        CalMETBEPhi_->push_back(i->phi());
        CalMETBEEn_->push_back(i->energy());

  }
  return;

}

template<typename CaloMETBEFOCollection>
void AODAnalyzer::fillCaloMETBEFOs(const edm::Handle<CaloMETBEFOCollection> & caloMETBEFOs)
{
  // std::cout << "fills is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename CaloMETBEFOCollection::const_iterator i = caloMETBEFOs->begin();
  for(;i != caloMETBEFOs->end(); i++){
        CalMETBEFOPt_->push_back(i->pt());
        // CalMETBEFOEta_->push_back(i->eta());
        CalMETBEFOPhi_->push_back(i->phi());
        CalMETBEFOEn_->push_back(i->energy());

  }
  return;

}

template<typename CaloMETMCollection>
void AODAnalyzer::fillCaloMETMs(const edm::Handle<CaloMETMCollection> & caloMETMs)
{
  // std::cout << "fills is being called!" << std::endl;
  // std::cout << pfmets->size() <<std::endl;
  typename CaloMETMCollection::const_iterator i = caloMETMs->begin();
  for(;i != caloMETMs->end(); i++){
        CalMETMPt_->push_back(i->pt());
        // CalMETMEta_->push_back(i->eta());
        CalMETMPhi_->push_back(i->phi());
        CalMETMEn_->push_back(i->energy());

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
     // if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // not sure if needed
      // {
      SCEn_->push_back(i->energy());
      SCEta_->push_back(i->eta());
      SCPhi_->push_back(i->phi());
      SCEtaWidth_->push_back(i->etaWidth());
      SCPhiWidth_->push_back(i->phiWidth());
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
     // if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // not sure if needed
      // {
      SCEnhfEM_->push_back(i->energy());
      SCEtahfEM_->push_back(i->eta());
      SCPhihfEM_->push_back(i->phi());
      // SCEtaWidthhfEM_->push_back(i->etaWidth());
      // SCPhiWidthhfEM_->push_back(i->phiWidth());
        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
      // }
  }
  return;


}

template<typename SuperCluster5x5Collection>
void AODAnalyzer::fillSC5x5(const edm::Handle<SuperCluster5x5Collection> & superclusters5x5) //ask for jets analogy //SUPERCLUSTERS
{

  // Selected jets
  //reco::CaloJetCollection recojets;
  typename SuperCluster5x5Collection::const_iterator i = superclusters5x5->begin();
  for(;i != superclusters5x5->end(); i++){
     // if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // not sure if needed
      // {
      SCEn5x5_->push_back(i->energy());
      SCEta5x5_->push_back(i->eta());
      SCPhi5x5_->push_back(i->phi());
      SCEtaWidth5x5_->push_back(i->etaWidth());
      SCPhiWidth5x5_->push_back(i->phiWidth());
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

template<typename CaloCluster5x5Collection>
void AODAnalyzer::fillCC5x5(const edm::Handle<CaloCluster5x5Collection> & caloclusters5x5) //ask for jets analogy //SUPERCLUSTERS
{

  // Selected jets
  //reco::CaloJetCollection recojets;
  typename CaloCluster5x5Collection::const_iterator i = caloclusters5x5->begin();
  for(;i != caloclusters5x5->end(); i++){
     // if(std::abs(i->eta()) < maxSCEta_ && i->energy() >= minSCEn_) // do I need something like maxCCeta and so on? TODOd
      // {
      CCEn5x5_->push_back(i->energy());
      CCEta5x5_->push_back(i->eta());
      CCPhi5x5_->push_back(i->phi());

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
     
        PhoPt_->push_back(i->pt());
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
                // std::cout << "ele Photons pt: " << i->pt()   << std::endl; //TODO

  }

        // std::cout << "end of Photons collision "  << std::endl; //TODO

  return;


}

template<typename PhotongedCollection>
void AODAnalyzer::fillgedPhotons(const edm::Handle<PhotongedCollection> & gedphotons)
{
   typename PhotongedCollection::const_iterator i = gedphotons->begin();
   for(;i != gedphotons->end(); i++){
     
        gedPhoPt_->push_back(i->pt());
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
   // https://github.com/cms-sw/cmssw/blob/1f04ccd3079678a4be345749f9ad13ea2edce89e/DataFormats/CTPPSDigi/src/CTPPSPixelDigiCollection.cc  -- Change muons to temporary
   // https://github.com/cms-sw/cmssw/blob/09c3fce6626f70fd04223e7dacebf0b485f73f54/SimCalorimetry/CaloSimAlgos/interface/CaloDigiCollectionSorter.h
   // https://github.com/cms-sw/cmssw/blob/09c3fce6626f70fd04223e7dacebf0b485f73f54/RecoTauTag/RecoTau/src/RecoTauCommonUtilities.cc 

   typename MuonCollection::const_iterator i = muons->begin();
   
   for(;i != muons->end(); i++){

        // Mu0Pt_->push_back(i->pt());
        MuPt_->push_back(i->pt());
        MuEta_->push_back(i->eta());
        MuPhi_->push_back(i->phi());
        MuEn_->push_back(i->energy());   //GETCORRECTEDENERGY!!
        MuCh_->push_back(i->charge());
        MuChi2_->push_back(i->vertexNormalizedChi2());  
        std::cout << "ele Muon pt: " << i->pt()   << std::endl; //TODO
        std::cout << "ele Muon eta: "  << i->eta() << std::endl;
        std::cout << "ele Muon phi: "  << i->phi() << std::endl;
        std::cout << "ele Muon energy: "  << i->energy() << std::endl;
        std::cout << "ele Muon charge: "  << i->charge() << std::endl;
        std::cout << "ele Muon vertexNormalizedChi2: "  << i->vertexNormalizedChi2() << std::endl;

        // std::sort(MuPt_->begin(), MuPt_->end()); 
        // std::cout << "ele Muon in the loop pt: " << i->pt()   << std::endl; //TODO
      // }
  }
        std::cout << "end of Muon collision "  << std::endl; //TODO
    

  // typename MuonCollection::const_iterator muonsBegin = muons->begin();
  // typename MuonCollection::const_iterator muonsEnd = muons->end().
    // std::vector<MuonCollection> muonsvector;
    // muonsvector.clear();
    // muonsvector.push_back(muons);
    // std::sort(muonsvector.begin(), muonsvector.end(), GreaterByPt<reco::Muon>()); 
             
  // muons->clear();
  std::sort( muons->begin(), muons->end(), GreaterByPt<reco::Muon>());
  //   [] ( const reco::Muon & muon1, const reco::Muon & muon2 ) { 
  //   return (muon1.pt() > muon2.pt());
  // } );


  Mu0Pt_->push_back(muons->begin()->pt()); 
          std::cout << "ele Jet 0pt: " << (muons->begin()->pt())   << std::endl;

  Mu1Pt_->push_back((muons->begin()+muons->size()/5)->pt()); 
        std::cout << "ele Jet 1pt: " << ((muons->begin()+muons->size()/5)->pt())    << std::endl;

  Mu2Pt_->push_back((muons->begin()+muons->size()*2/5)->pt());
        std::cout << "ele Jet 2pt: " << ((muons->begin()+muons->size()*2/5)->pt())  << std::endl;
  
  Mu3Pt_->push_back((muons->begin()+muons->size()*3/5)->pt());
        std::cout << "ele Jet 3pt: " << ((muons->begin()+muons->size()*3/5)->pt())  << std::endl;

  Mu4Pt_->push_back((muons->begin()+muons->size()*4/5)->pt());
        std::cout << "ele Jet 4pt: " << ((muons->begin()+muons->size()*4/5)->pt())  << std::endl; 

  Mu5Pt_->push_back((muons->begin()+muons->size()-1)->pt());
       std::cout << "ele Jet 5pt: " << ((muons->begin()+muons->size()-1)->pt())   << std::endl;
  //        // Mu0Eta_
         // Mu1Eta_
         // Mu2Eta_
         // Mu3Eta_
         // Mu4Eta_
         // Mu5Eta_

         // Mu0Phi_
         // Mu1Phi_
         // Mu2Phi_
         // Mu3Phi_
         // Mu4Phi_
         // Mu5Phi_

         // Mu0En_
         // Mu1En_
         // Mu2En_
         // Mu3En_
         // Mu4En_
         // Mu5En_






        // temporary.insert(std::end(temporary), sort_begin, sort_end);

        // std::sort(MuPt_->begin(), MuPt_->end()); 
        // Mu0Pt_->push_back(muons->begin()->pt());
        //     std::cout << "ele Mu0pt: " << (muons->begin()->pt())   << std::endl;
        // //         std::cout << "ele Muon pt begin: " << i->pt()   << std::endl; //TODO



  //-------------------------------------------------------------
  // typename CaloJetCollection::const_iterator i = calojets->begin();
  // for(;i != calojets->end(); i++){
  //       CalJetPt_->push_back(i->pt());
  //       CalJetEta_->push_back(i->eta());
  //       CalJetPhi_->push_back(i->phi());
  //       CalJetEn_->push_back(i->energy());
  //     // std::cout << "ele CaloJets pt: " << i->pt()   << std::endl; //TODO
  //     // std::cout << "ele CaloJets eta: " << i->eta()   << std::endl; //TODO
  //     // std::cout << "ele CaloJets phi: " << i->phi()   << std::endl; //TODO
  //     // std::cout << "ele CaloJets energy: " << i->energy()   << std::endl; //TODO


  // }

  //      // std::cout << "end of calojets collision "  << std::endl; //TODO
  // CalJet0Pt_->push_back(calojets->begin()->pt()); 
  //         // std::cout << "ele Jet 0pt: " << (calojets->begin()->pt())   << std::endl;

  //____________________________________________________________________






  return;


}


template<typename MuonCosmCollection>
void AODAnalyzer::fillCosmMuons(const edm::Handle<MuonCosmCollection> & muonsCosm)
{
   typename MuonCosmCollection::const_iterator i = muonsCosm->begin();
   for(;i != muonsCosm->end(); i++){
     
        MuCosmPt_->push_back(i->pt());
        MuCosmEta_->push_back(i->eta());
        MuCosmPhi_->push_back(i->phi());
        MuCosmEn_->push_back(i->energy());   //GETCORRECTEDENERGY!!
        MuCosmCh_->push_back(i->charge());
        MuCosmChi2_->push_back(i->vertexNormalizedChi2());
        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
         //  std::cout << "ele cosm charge: "  << i->charge() << std::endl;

      // }
  }
  return;


}



template<typename MuonCosmLegCollection>
void AODAnalyzer::fillCosmLegMuons(const edm::Handle<MuonCosmLegCollection> & muonsCosmLeg)
{
   typename MuonCosmLegCollection::const_iterator i = muonsCosmLeg->begin();
   for(;i != muonsCosmLeg->end(); i++){
     
        MuCosmLegPt_->push_back(i->pt());
        MuCosmLegEta_->push_back(i->eta());
        MuCosmLegPhi_->push_back(i->phi());
        MuCosmLegEn_->push_back(i->energy());   //GETCORRECTEDENERGY!!
        MuCosmLegCh_->push_back(i->charge());
        MuCosmLegChi2_->push_back(i->vertexNormalizedChi2());
        // std::cout << "ele energy: " << i->energy()   << std::endl; 
        // std::cout << "ele SCeta: "  << i->etaWidth() << std::endl;
        // std::cout << "ele SCphi: "  << i->phiWidth() << std::endl;
         //std::cout << "ele cosmlegcharge: "  << i->charge() << std::endl;

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
    // EBenergy_ ->push_back( 3.14 ); // E.g
    EBtime_ ->push_back(i->time());
    EBchi2_ ->push_back(i->chi2());
    // EBdetid_ = (i->detid());
    // DetId id();
    // float eta = geomH->getGeometry(i->detid())->getPosition().eta();
    EBiEta_->push_back( EBDetId( i->detid() ).ieta() );
    EBiPhi_->push_back( EBDetId( i->detid() ).iphi() );
     // std::cout << "ele EB energy: "  << i->energy()       << std::endl;


  }
       // std::cout << "ele EB end of collision"  << std::endl;

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
    EEix_->push_back( EEDetId( i->detid() ).ix() );
    EEiy_->push_back( EEDetId( i->detid() ).iy() );
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
    // ESchi2_ ->push_back(i->chi2()); 
    ESix_->push_back( ESDetId( i->detid() ).six() );
    ESiy_->push_back( ESDetId( i->detid() ).siy() );
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
    //

    HBHEenergy_ ->push_back(i->energy());
    HBHEtime_ ->push_back(i->time());
    HBHEauxe_ ->push_back(i->eaux()); //const class HBHERecHit' has no member named 'chi2'

    HBHEieta_ ->push_back( HcalDetId( i->detid() ).ieta() );
    HBHEiphi_ ->push_back( HcalDetId( i->detid() ).iphi() );
    // std::cout << "ele HBHEenergy: " << i->energy()   << std::endl; 
    //std::cout << "ele HBHEchi: "  << i->eaux() << std::endl;
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
    HFieta_ ->push_back(HcalDetId( i->detid() ).ieta());
    HFiphi_ ->push_back(HcalDetId( i->detid() ).iphi());
    // std::cout << "ele HFenergy: " << i->energy()   << std::endl; 
    // std::cout << "ele HFtime: "  << i->time() << std::endl;
    // HFchi2_ ->push_back(i->chi2());
    // std::cout << "ele HFenergy: " << i->energy()   << std::endl; 
    // std::cout << "ele HFtime: "  << i->time() << std::endl;
  }
  return;
}

// template<typename HORecHitCollection>
// void AODAnalyzer::fillHOrecHit(const edm::Handle<HORecHitCollection> & HOhits)
// {

//   // std::cout << "fillHOrecHit is being called!" << std::endl;
//   typename HORecHitCollection::const_iterator i = HOhits->begin();
//   for(;i != HOhits->end(); i++){
//     HOenergy_ ->push_back(i->energy());
//     HOtime_ ->push_back(i->time());
//     HOieta_ ->push_back(HcalDetId( i->detid() ).ieta());
//     HOiphi_ ->push_back(HcalDetId( i->detid() ).iphi());
//     // HOchi2_ ->push_back(i->chi2());
//     std::cout << "ele HOenergy: " << i->energy()   << std::endl; 
//     std::cout << "ele HOtime: "  << i->time() << std::endl;
//   }
//   return;
// }

template<typename PreshowerClusterCollection>
void AODAnalyzer::fillPreshowerCluster(const edm::Handle<PreshowerClusterCollection> & preshowerclusterhits)
{

  // std::cout << "fillHOrecHit is being called!" << std::endl;
  typename PreshowerClusterCollection::const_iterator i = preshowerclusterhits->begin();
  for(;i != preshowerclusterhits->end(); i++){

    PreShEn_->push_back(i->energy());
    // PreShCorrEn_->push_back(i->correctedEnergy());
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
    // PreShYCorrEn_->push_back(i->correctedEnergy());
    PreShYEta_->push_back(i->eta());
    PreShYPhi_->push_back(i->phi());

  }
  return;
}


// template<typename CastorTowerCollection>
// void AODAnalyzer::fillCastorTower(const edm::Handle<CastorTowerCollection> & castors)
// {

//   // std::cout << "fillHOrecHit is being called!" << std::endl;
//   typename CastorTowerCollection::const_iterator i = castors->begin();
//   for(;i != castors->end(); i++){

//     CTPt_->push_back(i->pt());
//     CTEta_->push_back(i->eta());
//     CTPhi_->push_back(i->phi());

//   }
//   return;
// }

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

//---------------------All sorted variables start.
//PFJet sorted variables  
  PFJet0Pt_   = new std::vector<float>;
  PFJet1Pt_   = new std::vector<float>;
  PFJet2Pt_   = new std::vector<float>;
  PFJet3Pt_   = new std::vector<float>;
  PFJet4Pt_   = new std::vector<float>;
  PFJet5Pt_   = new std::vector<float>;

  PFJet0Eta_  = new std::vector<float>;
  PFJet1Eta_  = new std::vector<float>;
  PFJet2Eta_  = new std::vector<float>;
  PFJet3Eta_  = new std::vector<float>;
  PFJet4Eta_  = new std::vector<float>;
  PFJet5Eta_  = new std::vector<float>;

  PFJet0Phi_  = new std::vector<float>;
  PFJet1Phi_  = new std::vector<float>;
  PFJet2Phi_  = new std::vector<float>;
  PFJet3Phi_  = new std::vector<float>;
  PFJet4Phi_  = new std::vector<float>;
  PFJet5Phi_  = new std::vector<float>;

//PFJet4CHS sorted variables
  PFJet4CHS0Pt_  = new std::vector<float>;
  PFJet4CHS1Pt_  = new std::vector<float>;
  PFJet4CHS2Pt_  = new std::vector<float>;
  PFJet4CHS3Pt_  = new std::vector<float>;
  PFJet4CHS4Pt_  = new std::vector<float>;
  PFJet4CHS5Pt_  = new std::vector<float>;

  PFJet4CHS0Eta_  = new std::vector<float>;
  PFJet4CHS1Eta_  = new std::vector<float>;
  PFJet4CHS2Eta_  = new std::vector<float>;
  PFJet4CHS3Eta_  = new std::vector<float>;
  PFJet4CHS4Eta_  = new std::vector<float>;
  PFJet4CHS5Eta_  = new std::vector<float>;

  PFJet4CHS0Phi_  = new std::vector<float>;
  PFJet4CHS1Phi_  = new std::vector<float>;
  PFJet4CHS2Phi_  = new std::vector<float>;
  PFJet4CHS3Phi_  = new std::vector<float>;
  PFJet4CHS4Phi_  = new std::vector<float>;
  PFJet4CHS5Phi_  = new std::vector<float>;
 //PF8CHS sorted variables 
  PFJet8CHS0Pt_  = new std::vector<float>;
  PFJet8CHS1Pt_  = new std::vector<float>;
  PFJet8CHS2Pt_  = new std::vector<float>;
  PFJet8CHS3Pt_  = new std::vector<float>;
  PFJet8CHS4Pt_  = new std::vector<float>;
  PFJet8CHS5Pt_  = new std::vector<float>;

  PFJet8CHS0Eta_  = new std::vector<float>;
  PFJet8CHS1Eta_  = new std::vector<float>;
  PFJet8CHS2Eta_  = new std::vector<float>;
  PFJet8CHS3Eta_  = new std::vector<float>;
  PFJet8CHS4Eta_  = new std::vector<float>;
  PFJet8CHS5Eta_  = new std::vector<float>;

  PFJet8CHS0Phi_  = new std::vector<float>;
  PFJet8CHS1Phi_  = new std::vector<float>;
  PFJet8CHS2Phi_  = new std::vector<float>;
  PFJet8CHS3Phi_  = new std::vector<float>;
  PFJet8CHS4Phi_  = new std::vector<float>;
  PFJet8CHS5Phi_  = new std::vector<float>;
  //PFJetEI sorted variables 
  PFJetEI0Pt_  = new std::vector<float>;
  PFJetEI1Pt_  = new std::vector<float>;
  PFJetEI2Pt_  = new std::vector<float>;
  PFJetEI3Pt_  = new std::vector<float>;
  PFJetEI4Pt_  = new std::vector<float>;
  PFJetEI5Pt_  = new std::vector<float>;

  PFJetEI0Eta_  = new std::vector<float>;
  PFJetEI1Eta_  = new std::vector<float>;
  PFJetEI2Eta_  = new std::vector<float>;
  PFJetEI3Eta_  = new std::vector<float>;
  PFJetEI4Eta_  = new std::vector<float>;
  PFJetEI5Eta_  = new std::vector<float>;

  PFJetEI0Phi_  = new std::vector<float>;
  PFJetEI1Phi_  = new std::vector<float>;
  PFJetEI2Phi_  = new std::vector<float>;
  PFJetEI3Phi_  = new std::vector<float>;
  PFJetEI4Phi_  = new std::vector<float>;
  PFJetEI5Phi_  = new std::vector<float>;

  // //8CHSSoftDrop sorted variables
  // PFJet8CHSSD0Pt_  = new std::vector<float>;
  // PFJet8CHSSD1Pt_  = new std::vector<float>;
  // PFJet8CHSSD2Pt_  = new std::vector<float>;
  // PFJet8CHSSD3Pt_  = new std::vector<float>;
  // PFJet8CHSSD4Pt_  = new std::vector<float>;
  // PFJet8CHSSD5Pt_  = new std::vector<float>;

  // PFJet8CHSSD0Eta_  = new std::vector<float>;
  // PFJet8CHSSD1Eta_  = new std::vector<float>;
  // PFJet8CHSSD2Eta_  = new std::vector<float>;
  // PFJet8CHSSD3Eta_  = new std::vector<float>;
  // PFJet8CHSSD4Eta_  = new std::vector<float>;
  // PFJet8CHSSD5Eta_  = new std::vector<float>;

  // PFJet8CHSSD0Phi_  = new std::vector<float>;
  // PFJet8CHSSD1Phi_  = new std::vector<float>;
  // PFJet8CHSSD2Phi_  = new std::vector<float>;
  // PFJet8CHSSD3Phi_  = new std::vector<float>;
  // PFJet8CHSSD4Phi_  = new std::vector<float>;
  // PFJet8CHSSD5Phi_  = new std::vector<float>;
  // //TopCHS sorted variables
  // PFJetTopCHS0Pt_  = new std::vector<float>;
  // PFJetTopCHS1Pt_  = new std::vector<float>;
  // PFJetTopCHS2Pt_  = new std::vector<float>;
  // PFJetTopCHS3Pt_  = new std::vector<float>;
  // PFJetTopCHS4Pt_  = new std::vector<float>;
  // PFJetTopCHS5Pt_  = new std::vector<float>;

  // PFJetTopCHS0Eta_  = new std::vector<float>;
  // PFJetTopCHS1Eta_  = new std::vector<float>;
  // PFJetTopCHS2Eta_  = new std::vector<float>;
  // PFJetTopCHS3Eta_  = new std::vector<float>;
  // PFJetTopCHS4Eta_  = new std::vector<float>;
  // PFJetTopCHS5Eta_  = new std::vector<float>;

  // PFJetTopCHS0Phi_  = new std::vector<float>;
  // PFJetTopCHS1Phi_  = new std::vector<float>;
  // PFJetTopCHS2Phi_  = new std::vector<float>;
  // PFJetTopCHS3Phi_  = new std::vector<float>;
  // PFJetTopCHS4Phi_  = new std::vector<float>;
  // PFJetTopCHS5Phi_  = new std::vector<float>;


  //CaloJet sorted variables
  CalJet0Pt_  = new std::vector<float>;
  CalJet1Pt_  = new std::vector<float>;
  CalJet2Pt_  = new std::vector<float>;
  CalJet3Pt_  = new std::vector<float>;
  CalJet4Pt_  = new std::vector<float>;
  CalJet5Pt_  = new std::vector<float>;

  CalJet0Eta_  = new std::vector<float>;
  CalJet1Eta_  = new std::vector<float>;
  CalJet2Eta_  = new std::vector<float>;
  CalJet3Eta_  = new std::vector<float>;
  CalJet4Eta_  = new std::vector<float>;
  CalJet5Eta_  = new std::vector<float>;

  CalJet0Phi_  = new std::vector<float>;
  CalJet1Phi_  = new std::vector<float>;
  CalJet2Phi_  = new std::vector<float>;
  CalJet3Phi_  = new std::vector<float>;
  CalJet4Phi_  = new std::vector<float>;
  CalJet5Phi_  = new std::vector<float>;

  CalJet0En_  = new std::vector<float>;
  CalJet1En_  = new std::vector<float>;
  CalJet2En_  = new std::vector<float>;
  CalJet3En_  = new std::vector<float>;
  CalJet4En_  = new std::vector<float>;
  CalJet5En_  = new std::vector<float>;


  //photon sorted variables (not all of them)
  Pho0Pt_  = new std::vector<float>;
  Pho1Pt_  = new std::vector<float>;
  Pho2Pt_  = new std::vector<float>;
  Pho3Pt_  = new std::vector<float>;
  Pho4Pt_  = new std::vector<float>;
  Pho5Pt_  = new std::vector<float>;

  Pho0Eta_  = new std::vector<float>;
  Pho1Eta_  = new std::vector<float>;
  Pho2Eta_  = new std::vector<float>;
  Pho3Eta_  = new std::vector<float>;
  Pho4Eta_  = new std::vector<float>;
  Pho5Eta_  = new std::vector<float>;

  Pho0Phi_  = new std::vector<float>;
  Pho1Phi_  = new std::vector<float>;
  Pho2Phi_  = new std::vector<float>;
  Pho3Phi_  = new std::vector<float>;
  Pho4Phi_  = new std::vector<float>;
  Pho5Phi_  = new std::vector<float>;

  Pho0En_  = new std::vector<float>;
  Pho1En_  = new std::vector<float>;
  Pho2En_  = new std::vector<float>;
  Pho3En_  = new std::vector<float>;
  Pho4En_  = new std::vector<float>;
  Pho5En_  = new std::vector<float>;

  //ged qPhotons sorted variables (not all of them)
  gedPho0Pt_  = new std::vector<float>;
  gedPho1Pt_  = new std::vector<float>;
  gedPho2Pt_  = new std::vector<float>;
  gedPho3Pt_  = new std::vector<float>;
  gedPho4Pt_  = new std::vector<float>;
  gedPho5Pt_  = new std::vector<float>;

  gedPho0Eta_  = new std::vector<float>;
  gedPho1Eta_  = new std::vector<float>;
  gedPho2Eta_  = new std::vector<float>;
  gedPho3Eta_  = new std::vector<float>;
  gedPho4Eta_  = new std::vector<float>;
  gedPho5Eta_  = new std::vector<float>;

  gedPho0Phi_  = new std::vector<float>;
  gedPho1Phi_  = new std::vector<float>;
  gedPho2Phi_  = new std::vector<float>;
  gedPho3Phi_  = new std::vector<float>;
  gedPho4Phi_  = new std::vector<float>;
  gedPho5Phi_  = new std::vector<float>;

  gedPho0En_  = new std::vector<float>;
  gedPho1En_  = new std::vector<float>;
  gedPho2En_  = new std::vector<float>;
  gedPho3En_  = new std::vector<float>;
  gedPho4En_  = new std::vector<float>;
  gedPho5En_  = new std::vector<float>;


    //Muon sorted variables (not all of them)
  Mu0Pt_  = new std::vector<float>;
  Mu1Pt_  = new std::vector<float>;
  Mu2Pt_  = new std::vector<float>;
  Mu3Pt_  = new std::vector<float>;
  Mu4Pt_  = new std::vector<float>;
  Mu5Pt_  = new std::vector<float>;

  Mu0Eta_  = new std::vector<float>;
  Mu1Eta_  = new std::vector<float>;
  Mu2Eta_  = new std::vector<float>;
  Mu3Eta_  = new std::vector<float>;
  Mu4Eta_  = new std::vector<float>;
  Mu5Eta_  = new std::vector<float>;

  Mu0Phi_  = new std::vector<float>;
  Mu1Phi_  = new std::vector<float>;
  Mu2Phi_  = new std::vector<float>;
  Mu3Phi_  = new std::vector<float>;
  Mu4Phi_  = new std::vector<float>;
  Mu5Phi_  = new std::vector<float>;

  Mu0En_  = new std::vector<float>;
  Mu1En_  = new std::vector<float>;
  Mu2En_  = new std::vector<float>;
  Mu3En_  = new std::vector<float>;
  Mu4En_  = new std::vector<float>;
  Mu5En_  = new std::vector<float>;

  //Muon sorted cosmic variables (not all of them)
  MuCosm0Pt_  = new std::vector<float>;
  MuCosm1Pt_  = new std::vector<float>;
  MuCosm2Pt_  = new std::vector<float>;
  MuCosm3Pt_  = new std::vector<float>;
  MuCosm4Pt_  = new std::vector<float>;
  MuCosm5Pt_  = new std::vector<float>;

  MuCosm0Eta_  = new std::vector<float>;
  MuCosm1Eta_  = new std::vector<float>;
  MuCosm2Eta_  = new std::vector<float>;
  MuCosm3Eta_  = new std::vector<float>;
  MuCosm4Eta_  = new std::vector<float>;
  MuCosm5Eta_  = new std::vector<float>;

  MuCosm0Phi_  = new std::vector<float>;
  MuCosm1Phi_  = new std::vector<float>;
  MuCosm2Phi_  = new std::vector<float>;
  MuCosm3Phi_  = new std::vector<float>;
  MuCosm4Phi_  = new std::vector<float>;
  MuCosm5Phi_  = new std::vector<float>;

  MuCosm0En_  = new std::vector<float>;
  MuCosm1En_  = new std::vector<float>;
  MuCosm2En_  = new std::vector<float>;
  MuCosm3En_  = new std::vector<float>;
  MuCosm4En_  = new std::vector<float>;
  MuCosm5En_  = new std::vector<float>;


  //Muon Cosmic1Leg sorted variables (not all of them)
  MuCosmLeg0Pt_  = new std::vector<float>;
  MuCosmLeg1Pt_  = new std::vector<float>;
  MuCosmLeg2Pt_  = new std::vector<float>;
  MuCosmLeg3Pt_  = new std::vector<float>;
  MuCosmLeg4Pt_  = new std::vector<float>;
  MuCosmLeg5Pt_  = new std::vector<float>;

  MuCosmLeg0Eta_  = new std::vector<float>;
  MuCosmLeg1Eta_  = new std::vector<float>;
  MuCosmLeg2Eta_  = new std::vector<float>;
  MuCosmLeg3Eta_  = new std::vector<float>;
  MuCosmLeg4Eta_  = new std::vector<float>;
  MuCosmLeg5Eta_  = new std::vector<float>;

  MuCosmLeg0Phi_  = new std::vector<float>;
  MuCosmLeg1Phi_  = new std::vector<float>;
  MuCosmLeg2Phi_  = new std::vector<float>;
  MuCosmLeg3Phi_  = new std::vector<float>;
  MuCosmLeg4Phi_  = new std::vector<float>;
  MuCosmLeg5Phi_  = new std::vector<float>;

  MuCosmLeg0En_  = new std::vector<float>;
  MuCosmLeg1En_  = new std::vector<float>;
  MuCosmLeg2En_  = new std::vector<float>;
  MuCosmLeg3En_  = new std::vector<float>;
  MuCosmLeg4En_  = new std::vector<float>;
  MuCosmLeg5En_  = new std::vector<float>;




//_____________________All sorted variables end.  

  PFJet4CHSPt_   = new std::vector<float>;
  PFJet4CHSEta_  = new std::vector<float>;
  PFJet4CHSPhi_  = new std::vector<float>;

  PFJet8CHSPt_   = new std::vector<float>;
  PFJet8CHSEta_  = new std::vector<float>;
  PFJet8CHSPhi_  = new std::vector<float>;

  PFJetEIPt_   = new std::vector<float>;
  PFJetEIEta_  = new std::vector<float>;
  PFJetEIPhi_  = new std::vector<float>;

  PFJet8CHSSDPt_   = new std::vector<float>;
  PFJet8CHSSDEta_  = new std::vector<float>;
  PFJet8CHSSDPhi_  = new std::vector<float>;

  PFJetTopCHSPt_   = new std::vector<float>;
  PFJetTopCHSEta_  = new std::vector<float>;
  PFJetTopCHSPhi_  = new std::vector<float>;        

  PFChMetPt_     = new std::vector<float>;
  // PFChMetEta_    = new std::vector<float>; //MET doesn't have ETA  
  PFChMetPhi_    = new std::vector<float>;
  PFMetPt_     = new std::vector<float>;
  // PFMetEta_    = new std::vector<float>; //MET doesn't have ETA 
  PFMetPhi_    = new std::vector<float>;
  nVtx_      = new std::vector<int>;
  
  CalJetPt_   = new std::vector<float>;
  CalJetEta_   = new std::vector<float>;
  CalJetPhi_   = new std::vector<float>;
  CalJetEn_   = new std::vector<float>;
  CalMETPt_   = new std::vector<float>;
  // CalMETEta_   = new std::vector<float>; //MET doesn't have ETA 
  CalMETPhi_   = new std::vector<float>;
  CalMETEn_   = new std::vector<float>;

  CalMETBEPt_   = new std::vector<float>;
  // CalMETBEEta_   = new std::vector<float>; //MET doesn't have ETA 
  CalMETBEPhi_   = new std::vector<float>;
  CalMETBEEn_   = new std::vector<float>;
  CalMETBEFOPt_   = new std::vector<float>;
  // CalMETBEFOEta_   = new std::vector<float>; //MET doesn't have ETA 
  CalMETBEFOPhi_   = new std::vector<float>;
  CalMETBEFOEn_   = new std::vector<float>;
  CalMETMPt_   = new std::vector<float>;
  // CalMETMEta_   = new std::vector<float>;  //MET doesn't have ETA 
  CalMETMPhi_   = new std::vector<float>;
  CalMETMEn_   = new std::vector<float>;


  SCEn_      = new std::vector<float>;
  SCEta_     = new std::vector<float>;
  SCPhi_     = new std::vector<float>;
  SCEtaWidth_     = new std::vector<float>;
  SCPhiWidth_     = new std::vector<float>;  
  SCEnhfEM_      = new std::vector<float>;
  SCEtahfEM_     = new std::vector<float>;
  SCPhihfEM_     = new std::vector<float>;
  // SCEtaWidthhfEM_     = new std::vector<float>;
  // SCPhiWidthhfEM_     = new std::vector<float>;  
  SCEn5x5_      = new std::vector<float>;
  SCEta5x5_     = new std::vector<float>;
  SCPhi5x5_     = new std::vector<float>;
  SCEtaWidth5x5_     = new std::vector<float>;
  SCPhiWidth5x5_     = new std::vector<float>;  
  CCEn_      = new std::vector<float>;
  CCEta_     = new std::vector<float>;
  CCPhi_     = new std::vector<float>;
  CCEn5x5_      = new std::vector<float>;
  CCEta5x5_     = new std::vector<float>;
  CCPhi5x5_     = new std::vector<float>;

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
  MuChi2_         = new std::vector<float>;  
  MuCosmPt_         = new std::vector<float>;
  MuCosmEta_        = new std::vector<float>;
  MuCosmPhi_        = new std::vector<float>;
  MuCosmEn_         = new std::vector<float>;
  MuCosmCh_         = new std::vector<float>;
  MuCosmChi2_         = new std::vector<float>;  
  MuCosmLegPt_         = new std::vector<float>;
  MuCosmLegEta_        = new std::vector<float>;
  MuCosmLegPhi_        = new std::vector<float>;
  MuCosmLegEn_         = new std::vector<float>;
  MuCosmLegCh_         = new std::vector<float>;
  MuCosmLegChi2_         = new std::vector<float>;  
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
  EBiEta_       = new std::vector<float>;
  EBiPhi_       = new std::vector<float>;

  EEenergy_    = new std::vector<float>;
  EEtime_      = new std::vector<float>;
  EEchi2_      = new std::vector<float>;
  EEix_       = new std::vector<float>;
  EEiy_       = new std::vector<float>;

  ESenergy_    = new std::vector<float>;
  EStime_      = new std::vector<float>;
  // ESchi2_      = new std::vector<float>;
  ESix_       = new std::vector<float>;
  ESiy_       = new std::vector<float>;

  HBHEenergy_  = new std::vector<float>;
  HBHEtime_    = new std::vector<float>;
  HBHEauxe_    = new std::vector<float>;
  HBHEieta_    = new std::vector<float>;
  HBHEiphi_    = new std::vector<float>;

  HFenergy_    = new std::vector<float>;
  HFtime_      = new std::vector<float>;
  HFieta_      = new std::vector<float>;
  HFiphi_      = new std::vector<float>;
  // HFchi2_      = new std::vector<float>;
  // HOenergy_    = new std::vector<float>;
  // HOtime_      = new std::vector<float>;
  // HOieta_      = new std::vector<float>;
  // HOiphi_      = new std::vector<float>;  
  // // HOchi2_      = new std::vector<float>;


  PreShEn_    = new std::vector<float>;
  // PreShCorrEn_= new std::vector<float>;
  PreShEta_   = new std::vector<float>;
  PreShPhi_   = new std::vector<float>;

  PreShYEn_    = new std::vector<float>;
  // PreShYCorrEn_= new std::vector<float>;
  PreShYEta_   = new std::vector<float>;
  PreShYPhi_   = new std::vector<float>;

  // CTPt_        = new std::vector<float>;
  // CTEta_       = new std::vector<float>;
  // CTPhi_       = new std::vector<float>;
  
  // outTree_->Branch("MetPt",     "std::vector<std::float>",     &MetPt_);
  // outTree_->Branch("MetPhi",    "std::vector<std::float>",     &MetPhi_);
  // outTree_->Branch("PFJetPt",     "std::vector<std::float>",     &PFJetPt_);
  // outTree_->Branch("PFJetEta",    "std::vector<std::float>",     &PFJetEta_);
  // outTree_->Branch("PFJetPhi",    "std::vector<std::float>",     &PFJetPhi_);
  // outTree_->Branch("nVtx",           "std::vector<std::int>",       &nVtx_);

  qPFJetPt_  = new std::vector<float>;
  qPFJetEta_ = new std::vector<float>;
  qPFJetPhi_ = new std::vector<float>;


//---------------------All sorted variables start.
//PFJet sorted variables  
  qPFJet0Pt_   = new std::vector<float>;
  qPFJet1Pt_   = new std::vector<float>;
  qPFJet2Pt_   = new std::vector<float>;
  qPFJet3Pt_   = new std::vector<float>;
  qPFJet4Pt_   = new std::vector<float>;
  qPFJet5Pt_   = new std::vector<float>;

  qPFJet0Eta_  = new std::vector<float>;
  qPFJet1Eta_  = new std::vector<float>;
  qPFJet2Eta_  = new std::vector<float>;
  qPFJet3Eta_  = new std::vector<float>;
  qPFJet4Eta_  = new std::vector<float>;
  qPFJet5Eta_  = new std::vector<float>;

  qPFJet0Phi_  = new std::vector<float>;
  qPFJet1Phi_  = new std::vector<float>;
  qPFJet2Phi_  = new std::vector<float>;
  qPFJet3Phi_  = new std::vector<float>;
  qPFJet4Phi_  = new std::vector<float>;
  qPFJet5Phi_  = new std::vector<float>;

//PFJet4CHS sorted variables
  qPFJet4CHS0Pt_  = new std::vector<float>;
  qPFJet4CHS1Pt_  = new std::vector<float>;
  qPFJet4CHS2Pt_  = new std::vector<float>;
  qPFJet4CHS3Pt_  = new std::vector<float>;
  qPFJet4CHS4Pt_  = new std::vector<float>;
  qPFJet4CHS5Pt_  = new std::vector<float>;

  qPFJet4CHS0Eta_  = new std::vector<float>;
  qPFJet4CHS1Eta_  = new std::vector<float>;
  qPFJet4CHS2Eta_  = new std::vector<float>;
  qPFJet4CHS3Eta_  = new std::vector<float>;
  qPFJet4CHS4Eta_  = new std::vector<float>;
  qPFJet4CHS5Eta_  = new std::vector<float>;

  qPFJet4CHS0Phi_  = new std::vector<float>;
  qPFJet4CHS1Phi_  = new std::vector<float>;
  qPFJet4CHS2Phi_  = new std::vector<float>;
  qPFJet4CHS3Phi_  = new std::vector<float>;
  qPFJet4CHS4Phi_  = new std::vector<float>;
  qPFJet4CHS5Phi_  = new std::vector<float>;
 //PF8CHS sorted variables 
  qPFJet8CHS0Pt_  = new std::vector<float>;
  qPFJet8CHS1Pt_  = new std::vector<float>;
  qPFJet8CHS2Pt_  = new std::vector<float>;
  qPFJet8CHS3Pt_  = new std::vector<float>;
  qPFJet8CHS4Pt_  = new std::vector<float>;
  qPFJet8CHS5Pt_  = new std::vector<float>;

  qPFJet8CHS0Eta_  = new std::vector<float>;
  qPFJet8CHS1Eta_  = new std::vector<float>;
  qPFJet8CHS2Eta_  = new std::vector<float>;
  qPFJet8CHS3Eta_  = new std::vector<float>;
  qPFJet8CHS4Eta_  = new std::vector<float>;
  qPFJet8CHS5Eta_  = new std::vector<float>;

  qPFJet8CHS0Phi_  = new std::vector<float>;
  qPFJet8CHS1Phi_  = new std::vector<float>;
  qPFJet8CHS2Phi_  = new std::vector<float>;
  qPFJet8CHS3Phi_  = new std::vector<float>;
  qPFJet8CHS4Phi_  = new std::vector<float>;
  qPFJet8CHS5Phi_  = new std::vector<float>;

  //PFJetEI sorted variables 
  qPFJetEI0Pt_  = new std::vector<float>;
  qPFJetEI1Pt_  = new std::vector<float>;
  qPFJetEI2Pt_  = new std::vector<float>;
  qPFJetEI3Pt_  = new std::vector<float>;
  qPFJetEI4Pt_  = new std::vector<float>;
  qPFJetEI5Pt_  = new std::vector<float>;

  qPFJetEI0Eta_  = new std::vector<float>;
  qPFJetEI1Eta_  = new std::vector<float>;
  qPFJetEI2Eta_  = new std::vector<float>;
  qPFJetEI3Eta_  = new std::vector<float>;
  qPFJetEI4Eta_  = new std::vector<float>;
  qPFJetEI5Eta_  = new std::vector<float>;

  qPFJetEI0Phi_  = new std::vector<float>;
  qPFJetEI1Phi_  = new std::vector<float>;
  qPFJetEI2Phi_  = new std::vector<float>;
  qPFJetEI3Phi_  = new std::vector<float>;
  qPFJetEI4Phi_  = new std::vector<float>;
  qPFJetEI5Phi_  = new std::vector<float>;

  // //8CHSSoftDrop sorted variables
  // qPFJet8CHSSD0Pt_  = new std::vector<float>;
  // qPFJet8CHSSD1Pt_  = new std::vector<float>;
  // qPFJet8CHSSD2Pt_  = new std::vector<float>;
  // qPFJet8CHSSD3Pt_  = new std::vector<float>;
  // qPFJet8CHSSD4Pt_  = new std::vector<float>;
  // qPFJet8CHSSD5Pt_  = new std::vector<float>;

  // qPFJet8CHSSD0Eta_  = new std::vector<float>;
  // qPFJet8CHSSD1Eta_  = new std::vector<float>;
  // qPFJet8CHSSD2Eta_  = new std::vector<float>;
  // qPFJet8CHSSD3Eta_  = new std::vector<float>;
  // qPFJet8CHSSD4Eta_  = new std::vector<float>;
  // qPFJet8CHSSD5Eta_  = new std::vector<float>;

  // qPFJet8CHSSD0Phi_  = new std::vector<float>;
  // qPFJet8CHSSD1Phi_  = new std::vector<float>;
  // qPFJet8CHSSD2Phi_  = new std::vector<float>;
  // qPFJet8CHSSD3Phi_  = new std::vector<float>;
  // qPFJet8CHSSD4Phi_  = new std::vector<float>;
  // qPFJet8CHSSD5Phi_  = new std::vector<float>;
  // //TopCHS sorted variables
  // qPFJetTopCHS0Pt_  = new std::vector<float>;
  // qPFJetTopCHS1Pt_  = new std::vector<float>;
  // qPFJetTopCHS2Pt_  = new std::vector<float>;
  // qPFJetTopCHS3Pt_  = new std::vector<float>;
  // qPFJetTopCHS4Pt_  = new std::vector<float>;
  // qPFJetTopCHS5Pt_  = new std::vector<float>;

  // qPFJetTopCHS0Eta_  = new std::vector<float>;
  // qPFJetTopCHS1Eta_  = new std::vector<float>;
  // qPFJetTopCHS2Eta_  = new std::vector<float>;
  // qPFJetTopCHS3Eta_  = new std::vector<float>;
  // qPFJetTopCHS4Eta_  = new std::vector<float>;
  // qPFJetTopCHS5Eta_  = new std::vector<float>;

  // qPFJetTopCHS0Phi_  = new std::vector<float>;
  // qPFJetTopCHS1Phi_  = new std::vector<float>;
  // qPFJetTopCHS2Phi_  = new std::vector<float>;
  // qPFJetTopCHS3Phi_  = new std::vector<float>;
  // qPFJetTopCHS4Phi_  = new std::vector<float>;
  // qPFJetTopCHS5Phi_  = new std::vector<float>;


  //CaloJet sorted variables
  qCalJet0Pt_  = new std::vector<float>;
  qCalJet1Pt_  = new std::vector<float>;
  qCalJet2Pt_  = new std::vector<float>;
  qCalJet3Pt_  = new std::vector<float>;
  qCalJet4Pt_  = new std::vector<float>;
  qCalJet5Pt_  = new std::vector<float>;

  qCalJet0Eta_  = new std::vector<float>;
  qCalJet1Eta_  = new std::vector<float>;
  qCalJet2Eta_  = new std::vector<float>;
  qCalJet3Eta_  = new std::vector<float>;
  qCalJet4Eta_  = new std::vector<float>;
  qCalJet5Eta_  = new std::vector<float>;

  qCalJet0Phi_  = new std::vector<float>;
  qCalJet1Phi_  = new std::vector<float>;
  qCalJet2Phi_  = new std::vector<float>;
  qCalJet3Phi_  = new std::vector<float>;
  qCalJet4Phi_  = new std::vector<float>;
  qCalJet5Phi_  = new std::vector<float>;

  qCalJet0En_  = new std::vector<float>;
  qCalJet1En_  = new std::vector<float>;
  qCalJet2En_  = new std::vector<float>;
  qCalJet3En_  = new std::vector<float>;
  qCalJet4En_  = new std::vector<float>;
  qCalJet5En_  = new std::vector<float>;


  //photon sorted variables (not all of them)
  qPho0Pt_  = new std::vector<float>;
  qPho1Pt_  = new std::vector<float>;
  qPho2Pt_  = new std::vector<float>;
  qPho3Pt_  = new std::vector<float>;
  qPho4Pt_  = new std::vector<float>;
  qPho5Pt_  = new std::vector<float>;

  qPho0Eta_  = new std::vector<float>;
  qPho1Eta_  = new std::vector<float>;
  qPho2Eta_  = new std::vector<float>;
  qPho3Eta_  = new std::vector<float>;
  qPho4Eta_  = new std::vector<float>;
  qPho5Eta_  = new std::vector<float>;

  qPho0Phi_  = new std::vector<float>;
  qPho1Phi_  = new std::vector<float>;
  qPho2Phi_  = new std::vector<float>;
  qPho3Phi_  = new std::vector<float>;
  qPho4Phi_  = new std::vector<float>;
  qPho5Phi_  = new std::vector<float>;

  qPho0En_  = new std::vector<float>;
  qPho1En_  = new std::vector<float>;
  qPho2En_  = new std::vector<float>;
  qPho3En_  = new std::vector<float>;
  qPho4En_  = new std::vector<float>;
  qPho5En_  = new std::vector<float>;

  //ged qPhotons sorted variables (not all of them)
  qgedPho0Pt_  = new std::vector<float>;
  qgedPho1Pt_  = new std::vector<float>;
  qgedPho2Pt_  = new std::vector<float>;
  qgedPho3Pt_  = new std::vector<float>;
  qgedPho4Pt_  = new std::vector<float>;
  qgedPho5Pt_  = new std::vector<float>;

  qgedPho0Eta_  = new std::vector<float>;
  qgedPho1Eta_  = new std::vector<float>;
  qgedPho2Eta_  = new std::vector<float>;
  qgedPho3Eta_  = new std::vector<float>;
  qgedPho4Eta_  = new std::vector<float>;
  qgedPho5Eta_  = new std::vector<float>;

  qgedPho0Phi_  = new std::vector<float>;
  qgedPho1Phi_  = new std::vector<float>;
  qgedPho2Phi_  = new std::vector<float>;
  qgedPho3Phi_  = new std::vector<float>;
  qgedPho4Phi_  = new std::vector<float>;
  qgedPho5Phi_  = new std::vector<float>;

  qgedPho0En_  = new std::vector<float>;
  qgedPho1En_  = new std::vector<float>;
  qgedPho2En_  = new std::vector<float>;
  qgedPho3En_  = new std::vector<float>;
  qgedPho4En_  = new std::vector<float>;
  qgedPho5En_  = new std::vector<float>;


    //Muon sorted variables (not all of them)
  qMu0Pt_  = new std::vector<float>;
  qMu1Pt_  = new std::vector<float>;
  qMu2Pt_  = new std::vector<float>;
  qMu3Pt_  = new std::vector<float>;
  qMu4Pt_  = new std::vector<float>;
  qMu5Pt_  = new std::vector<float>;

  qMu0Eta_  = new std::vector<float>;
  qMu1Eta_  = new std::vector<float>;
  qMu2Eta_  = new std::vector<float>;
  qMu3Eta_  = new std::vector<float>;
  qMu4Eta_  = new std::vector<float>;
  qMu5Eta_  = new std::vector<float>;

  qMu0Phi_  = new std::vector<float>;
  qMu1Phi_  = new std::vector<float>;
  qMu2Phi_  = new std::vector<float>;
  qMu3Phi_  = new std::vector<float>;
  qMu4Phi_  = new std::vector<float>;
  qMu5Phi_  = new std::vector<float>;

  qMu0En_  = new std::vector<float>;
  qMu1En_  = new std::vector<float>;
  qMu2En_  = new std::vector<float>;
  qMu3En_  = new std::vector<float>;
  qMu4En_  = new std::vector<float>;
  qMu5En_  = new std::vector<float>;

  //Muon sorted cosmic variables (not all of them)
  qMuCosm0Pt_  = new std::vector<float>;
  qMuCosm1Pt_  = new std::vector<float>;
  qMuCosm2Pt_  = new std::vector<float>;
  qMuCosm3Pt_  = new std::vector<float>;
  qMuCosm4Pt_  = new std::vector<float>;
  qMuCosm5Pt_  = new std::vector<float>;

  qMuCosm0Eta_  = new std::vector<float>;
  qMuCosm1Eta_  = new std::vector<float>;
  qMuCosm2Eta_  = new std::vector<float>;
  qMuCosm3Eta_  = new std::vector<float>;
  qMuCosm4Eta_  = new std::vector<float>;
  qMuCosm5Eta_  = new std::vector<float>;

  qMuCosm0Phi_  = new std::vector<float>;
  qMuCosm1Phi_  = new std::vector<float>;
  qMuCosm2Phi_  = new std::vector<float>;
  qMuCosm3Phi_  = new std::vector<float>;
  qMuCosm4Phi_  = new std::vector<float>;
  qMuCosm5Phi_  = new std::vector<float>;

  qMuCosm0En_  = new std::vector<float>;
  qMuCosm1En_  = new std::vector<float>;
  qMuCosm2En_  = new std::vector<float>;
  qMuCosm3En_  = new std::vector<float>;
  qMuCosm4En_  = new std::vector<float>;
  qMuCosm5En_  = new std::vector<float>;


  //Muon Cosmic1Leg sorted variables (not all of them)
  qMuCosmLeg0Pt_  = new std::vector<float>;
  qMuCosmLeg1Pt_  = new std::vector<float>;
  qMuCosmLeg2Pt_  = new std::vector<float>;
  qMuCosmLeg3Pt_  = new std::vector<float>;
  qMuCosmLeg4Pt_  = new std::vector<float>;
  qMuCosmLeg5Pt_  = new std::vector<float>;

  qMuCosmLeg0Eta_  = new std::vector<float>;
  qMuCosmLeg1Eta_  = new std::vector<float>;
  qMuCosmLeg2Eta_  = new std::vector<float>;
  qMuCosmLeg3Eta_  = new std::vector<float>;
  qMuCosmLeg4Eta_  = new std::vector<float>;
  qMuCosmLeg5Eta_  = new std::vector<float>;

  qMuCosmLeg0Phi_  = new std::vector<float>;
  qMuCosmLeg1Phi_  = new std::vector<float>;
  qMuCosmLeg2Phi_  = new std::vector<float>;
  qMuCosmLeg3Phi_  = new std::vector<float>;
  qMuCosmLeg4Phi_  = new std::vector<float>;
  qMuCosmLeg5Phi_  = new std::vector<float>;

  qMuCosmLeg0En_  = new std::vector<float>;
  qMuCosmLeg1En_  = new std::vector<float>;
  qMuCosmLeg2En_  = new std::vector<float>;
  qMuCosmLeg3En_  = new std::vector<float>;
  qMuCosmLeg4En_  = new std::vector<float>;
  qMuCosmLeg5En_  = new std::vector<float>;



//_____________________All sorted variables end.  


  qPFJet4CHSPt_   = new std::vector<float>;
  qPFJet4CHSEta_  = new std::vector<float>;
  qPFJet4CHSPhi_  = new std::vector<float>;

  qPFJet8CHSPt_   = new std::vector<float>;
  qPFJet8CHSEta_  = new std::vector<float>;
  qPFJet8CHSPhi_  = new std::vector<float>;

  qPFJetEIPt_   = new std::vector<float>;
  qPFJetEIEta_  = new std::vector<float>;
  qPFJetEIPhi_  = new std::vector<float>;

  qPFJet8CHSSDPt_   = new std::vector<float>;
  qPFJet8CHSSDEta_  = new std::vector<float>;
  qPFJet8CHSSDPhi_  = new std::vector<float>;

  qPFJetTopCHSPt_   = new std::vector<float>;
  qPFJetTopCHSEta_  = new std::vector<float>;
  qPFJetTopCHSPhi_  = new std::vector<float>;

  qPFChMetPt_     = new std::vector<float>;
  // qPFChMetEta_    = new std::vector<float>;  
  qPFChMetPhi_    = new std::vector<float>;
  qPFMetPt_     = new std::vector<float>;
  // qPFMetEta_    = new std::vector<float>;
  qPFMetPhi_    = new std::vector<float>;
 
  qCalJetPt_   = new std::vector<float>;
  qCalJetEta_   = new std::vector<float>;
  qCalJetPhi_   = new std::vector<float>;
  qCalJetEn_   = new std::vector<float>;
  qCalMETPt_   = new std::vector<float>;
  // qCalMETEta_   = new std::vector<float>;
  qCalMETPhi_   = new std::vector<float>;
  qCalMETEn_   = new std::vector<float>;

  qCalMETBEPt_   = new std::vector<float>;
  // qCalMETBEEta_   = new std::vector<float>;
  qCalMETBEPhi_   = new std::vector<float>;
  qCalMETBEEn_   = new std::vector<float>;
  qCalMETBEFOPt_   = new std::vector<float>;
  // qCalMETBEFOEta_   = new std::vector<float>;
  qCalMETBEFOPhi_   = new std::vector<float>;
  qCalMETBEFOEn_   = new std::vector<float>;
  qCalMETMPt_   = new std::vector<float>;
  // qCalMETMEta_   = new std::vector<float>;
  qCalMETMPhi_   = new std::vector<float>;
  qCalMETMEn_   = new std::vector<float>;

  qSCEn_      = new std::vector<float>;
  qSCEta_     = new std::vector<float>;
  qSCPhi_     = new std::vector<float>;
  qSCEtaWidth_     = new std::vector<float>;
  qSCPhiWidth_     = new std::vector<float>;  
  qSCEnhfEM_      = new std::vector<float>;
  qSCEtahfEM_     = new std::vector<float>;
  qSCPhihfEM_     = new std::vector<float>;
  // qSCEtaWidthhfEM_     = new std::vector<float>;
  // qSCPhiWidthhfEM_     = new std::vector<float>;  
  qSCEn5x5_      = new std::vector<float>;
  qSCEta5x5_     = new std::vector<float>;
  qSCPhi5x5_     = new std::vector<float>;
  qSCEtaWidth5x5_     = new std::vector<float>;
  qSCPhiWidth5x5_     = new std::vector<float>; 
  qCCEn_     = new std::vector<float>;
  qCCEta_    = new std::vector<float>;
  qCCPhi_    = new std::vector<float>;
  qCCEn5x5_     = new std::vector<float>;
  qCCEta5x5_    = new std::vector<float>;
  qCCPhi5x5_    = new std::vector<float>;

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

  qMuPt_         = new std::vector<float>;
  qMuEta_        = new std::vector<float>;
  qMuPhi_        = new std::vector<float>;
  qMuEn_         = new std::vector<float>;
  qMuCh_         = new std::vector<float>;
  qMuChi2_         = new std::vector<float>;  
  qMuCosmPt_         = new std::vector<float>;
  qMuCosmEta_        = new std::vector<float>;
  qMuCosmPhi_        = new std::vector<float>;
  qMuCosmEn_         = new std::vector<float>;
  qMuCosmCh_         = new std::vector<float>;
  qMuCosmChi2_         = new std::vector<float>;  
  qMuCosmLegPt_         = new std::vector<float>;
  qMuCosmLegEta_        = new std::vector<float>;
  qMuCosmLegPhi_        = new std::vector<float>;
  qMuCosmLegEn_         = new std::vector<float>;
  qMuCosmLegCh_         = new std::vector<float>;
  qMuCosmLegChi2_         = new std::vector<float>;
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


  qEBenergy_    = new std::vector<float>;
  qEBtime_      = new std::vector<float>;
  qEBchi2_      = new std::vector<float>;
  qEBiEta_       = new std::vector<float>;
  qEBiPhi_       = new std::vector<float>;

  qEEenergy_    = new std::vector<float>;
  qEEtime_      = new std::vector<float>;
  qEEchi2_      = new std::vector<float>;
  qEEix_       = new std::vector<float>;
  qEEiy_       = new std::vector<float>;

  qESenergy_    = new std::vector<float>;
  qEStime_      = new std::vector<float>;
  // qESchi2_      = new std::vector<float>;
  qESix_       = new std::vector<float>;
  qESiy_       = new std::vector<float>;

  qHBHEenergy_  = new std::vector<float>;
  qHBHEtime_    = new std::vector<float>;
  qHBHEauxe_    = new std::vector<float>;
  qHBHEieta_    = new std::vector<float>;
  qHBHEiphi_    = new std::vector<float>;

  qHFenergy_    = new std::vector<float>;
  qHFtime_      = new std::vector<float>;
  qHFieta_      = new std::vector<float>;
  qHFiphi_      = new std::vector<float>;
  // qHFchi2_      = new std::vector<float>;
//   qHOenergy_    = new std::vector<float>;
//   qHOtime_      = new std::vector<float>;
//   qHOieta_      = new std::vector<float>;
//   qHOiphi_      = new std::vector<float>;  
//   // qHOchi2_      = new std::vector<float>;

  qPreShEn_     = new std::vector<float>;
  // qPreShCorrEn_ = new std::vector<float>;
  qPreShEta_    = new std::vector<float>;
  qPreShPhi_    = new std::vector<float>;
  qPreShYEn_    = new std::vector<float>;
  // qPreShYCorrEn_= new std::vector<float>;
  qPreShYEta_   = new std::vector<float>;
  qPreShYPhi_   = new std::vector<float>;

  // qCTPt_        = new std::vector<float>;
  // qCTEta_       = new std::vector<float>;
  // qCTPhi_       = new std::vector<float>;

  qNVtx_       = new std::vector<int>;
  crossSection_= new std::vector<float>;
  pathRates_   = new std::vector<float>;
  pathNames_   = new std::vector<std::string>;
  outTree_->Branch("qPFJetPt",     "std::vector<std::float>",      &qPFJetPt_);
  outTree_->Branch("qPFJetEta",    "std::vector<std::float>",      &qPFJetEta_);
  outTree_->Branch("qPFJetPhi",    "std::vector<std::float>",      &qPFJetPhi_);

//-------------------------------------------Sorted variables start here.

//PFJet sorted variables
  outTree_->Branch("qPFJet0Pt",    "std::vector<std::float>",       qPFJet0Pt_);
  outTree_->Branch("qPFJet1Pt",    "std::vector<std::float>",       qPFJet1Pt_); 
  outTree_->Branch("qPFJet2Pt",    "std::vector<std::float>",       qPFJet2Pt_); 
  outTree_->Branch("qPFJet3Pt",    "std::vector<std::float>",       qPFJet3Pt_); 
  outTree_->Branch("qPFJet4Pt",    "std::vector<std::float>",       qPFJet4Pt_); 
  outTree_->Branch("qPFJet5Pt",    "std::vector<std::float>",       qPFJet5Pt_); 

  outTree_->Branch("qPFJet0Eta",     "std::vector<std::float>",       qPFJet0Eta_);
  outTree_->Branch("qPFJet1Eta",     "std::vector<std::float>",       qPFJet1Eta_);
  outTree_->Branch("qPFJet2Eta",     "std::vector<std::float>",       qPFJet2Eta_);
  outTree_->Branch("qPFJet3Eta",     "std::vector<std::float>",       qPFJet3Eta_);
  outTree_->Branch("qPFJet4Eta",     "std::vector<std::float>",       qPFJet4Eta_);
  outTree_->Branch("qPFJet5Eta",     "std::vector<std::float>",       qPFJet5Eta_);

  outTree_->Branch("qPFJet0Phi",     "std::vector<std::float>",       qPFJet0Phi_);
  outTree_->Branch("qPFJet1Phi",     "std::vector<std::float>",       qPFJet1Phi_);
  outTree_->Branch("qPFJet2Phi",     "std::vector<std::float>",       qPFJet2Phi_);
  outTree_->Branch("qPFJet3Phi",     "std::vector<std::float>",       qPFJet3Phi_);
  outTree_->Branch("qPFJet4Phi",     "std::vector<std::float>",       qPFJet4Phi_);
  outTree_->Branch("qPFJet5Phi",     "std::vector<std::float>",       qPFJet5Phi_);

  //PFJet4CHS sorted variables
  outTree_->Branch("qPFJet4CHS0Pt",   "std::vector<std::float>",       qPFJet4CHS0Pt_);
  outTree_->Branch("qPFJet4CHS1Pt",   "std::vector<std::float>",       qPFJet4CHS1Pt_);
  outTree_->Branch("qPFJet4CHS2Pt",   "std::vector<std::float>",       qPFJet4CHS2Pt_);
  outTree_->Branch("qPFJet4CHS3Pt",   "std::vector<std::float>",       qPFJet4CHS3Pt_);
  outTree_->Branch("qPFJet4CHS4Pt",   "std::vector<std::float>",       qPFJet4CHS4Pt_);
  outTree_->Branch("qPFJet4CHS5Pt",   "std::vector<std::float>",       qPFJet4CHS5Pt_);

  outTree_->Branch("qPFJet4CHS0Eta",     "std::vector<std::float>",       qPFJet4CHS0Eta_);
  outTree_->Branch("qPFJet4CHS1Eta",     "std::vector<std::float>",       qPFJet4CHS1Eta_);
  outTree_->Branch("qPFJet4CHS2Eta",     "std::vector<std::float>",       qPFJet4CHS2Eta_);
  outTree_->Branch("qPFJet4CHS3Eta",     "std::vector<std::float>",       qPFJet4CHS3Eta_);
  outTree_->Branch("qPFJet4CHS4Eta",     "std::vector<std::float>",       qPFJet4CHS4Eta_);
  outTree_->Branch("qPFJet4CHS5Eta",     "std::vector<std::float>",       qPFJet4CHS5Eta_);

  outTree_->Branch("qPFJet4CHS0Phi",     "std::vector<std::float>",       qPFJet4CHS0Phi_);
  outTree_->Branch("qPFJet4CHS1Phi",     "std::vector<std::float>",       qPFJet4CHS1Phi_);
  outTree_->Branch("qPFJet4CHS2Phi",     "std::vector<std::float>",       qPFJet4CHS2Phi_);
  outTree_->Branch("qPFJet4CHS3Phi",     "std::vector<std::float>",       qPFJet4CHS3Phi_);
  outTree_->Branch("qPFJet4CHS4Phi",     "std::vector<std::float>",       qPFJet4CHS4Phi_);
  outTree_->Branch("qPFJet4CHS5Phi",     "std::vector<std::float>",       qPFJet4CHS5Phi_);

 //PF8CHS sorted variables
  outTree_->Branch("qPFJet8CHS0Pt",     "std::vector<std::float>",       qPFJet8CHS0Pt_);
  outTree_->Branch("qPFJet8CHS1Pt",     "std::vector<std::float>",       qPFJet8CHS1Pt_);
  outTree_->Branch("qPFJet8CHS2Pt",     "std::vector<std::float>",       qPFJet8CHS2Pt_);
  outTree_->Branch("qPFJet8CHS3Pt",     "std::vector<std::float>",       qPFJet8CHS3Pt_);
  outTree_->Branch("qPFJet8CHS4Pt",     "std::vector<std::float>",       qPFJet8CHS4Pt_);
  outTree_->Branch("qPFJet8CHS5Pt",     "std::vector<std::float>",       qPFJet8CHS5Pt_);

  outTree_->Branch("qPFJet8CHS0Eta",     "std::vector<std::float>",       qPFJet8CHS0Eta_);
  outTree_->Branch("qPFJet8CHS1Eta",     "std::vector<std::float>",       qPFJet8CHS1Eta_);
  outTree_->Branch("qPFJet8CHS2Eta",     "std::vector<std::float>",       qPFJet8CHS2Eta_);
  outTree_->Branch("qPFJet8CHS3Eta",     "std::vector<std::float>",       qPFJet8CHS3Eta_);
  outTree_->Branch("qPFJet8CHS4Eta",     "std::vector<std::float>",       qPFJet8CHS4Eta_);
  outTree_->Branch("qPFJet8CHS5Eta",     "std::vector<std::float>",       qPFJet8CHS5Eta_);

  outTree_->Branch("qPFJet8CHS0Phi",     "std::vector<std::float>",       qPFJet8CHS0Phi_);
  outTree_->Branch("qPFJet8CHS1Phi",     "std::vector<std::float>",       qPFJet8CHS1Phi_);
  outTree_->Branch("qPFJet8CHS2Phi",     "std::vector<std::float>",       qPFJet8CHS2Phi_);
  outTree_->Branch("qPFJet8CHS3Phi",     "std::vector<std::float>",       qPFJet8CHS3Phi_);
  outTree_->Branch("qPFJet8CHS4Phi",     "std::vector<std::float>",       qPFJet8CHS4Phi_);
  outTree_->Branch("qPFJet8CHS5Phi",     "std::vector<std::float>",       qPFJet8CHS5Phi_);

//PFJetEI sorted variables

  outTree_->Branch("qPFJetEI0Pt",     "std::vector<std::float>",       qPFJetEI0Pt_);
  outTree_->Branch("qPFJetEI1Pt",     "std::vector<std::float>",       qPFJetEI1Pt_);
  outTree_->Branch("qPFJetEI2Pt",     "std::vector<std::float>",       qPFJetEI2Pt_);
  outTree_->Branch("qPFJetEI3Pt",     "std::vector<std::float>",       qPFJetEI3Pt_);
  outTree_->Branch("qPFJetEI4Pt",     "std::vector<std::float>",       qPFJetEI4Pt_);
  outTree_->Branch("qPFJetEI5Pt",     "std::vector<std::float>",       qPFJetEI5Pt_);

  outTree_->Branch("qPFJetEI0Eta",     "std::vector<std::float>",       qPFJetEI0Eta_);
  outTree_->Branch("qPFJetEI1Eta",     "std::vector<std::float>",       qPFJetEI1Eta_);
  outTree_->Branch("qPFJetEI2Eta",     "std::vector<std::float>",       qPFJetEI2Eta_);
  outTree_->Branch("qPFJetEI3Eta",     "std::vector<std::float>",       qPFJetEI3Eta_);
  outTree_->Branch("qPFJetEI4Eta",     "std::vector<std::float>",       qPFJetEI4Eta_);
  outTree_->Branch("qPFJetEI5Eta",     "std::vector<std::float>",       qPFJetEI5Eta_);

  outTree_->Branch("qPFJetEI0Phi",     "std::vector<std::float>",       qPFJetEI0Phi_);
  outTree_->Branch("qPFJetEI1Phi",     "std::vector<std::float>",       qPFJetEI1Phi_);
  outTree_->Branch("qPFJetEI2Phi",     "std::vector<std::float>",       qPFJetEI2Phi_);
  outTree_->Branch("qPFJetEI3Phi",     "std::vector<std::float>",       qPFJetEI3Phi_);
  outTree_->Branch("qPFJetEI4Phi",     "std::vector<std::float>",       qPFJetEI4Phi_);
  outTree_->Branch("qPFJetEI5Phi",     "std::vector<std::float>",       qPFJetEI5Phi_);

// //8CHSSoftDrop sorted variables

//   outTree_->Branch("qPFJet8CHSSD0Pt",     "std::vector<std::float>",       qPFJet8CHSSD0Pt_);
//   outTree_->Branch("qPFJet8CHSSD1Pt",     "std::vector<std::float>",       qPFJet8CHSSD1Pt_);
//   outTree_->Branch("qPFJet8CHSSD2Pt",     "std::vector<std::float>",       qPFJet8CHSSD2Pt_);
//   outTree_->Branch("qPFJet8CHSSD3Pt",     "std::vector<std::float>",       qPFJet8CHSSD3Pt_);
//   outTree_->Branch("qPFJet8CHSSD4Pt",     "std::vector<std::float>",       qPFJet8CHSSD4Pt_);
//   outTree_->Branch("qPFJet8CHSSD5Pt",     "std::vector<std::float>",       qPFJet8CHSSD5Pt_);

//   outTree_->Branch("qPFJet8CHSSD0Eta",     "std::vector<std::float>",       qPFJet8CHSSD0Eta_);
//   outTree_->Branch("qPFJet8CHSSD1Eta",     "std::vector<std::float>",       qPFJet8CHSSD1Eta_);
//   outTree_->Branch("qPFJet8CHSSD2Eta",     "std::vector<std::float>",       qPFJet8CHSSD2Eta_);
//   outTree_->Branch("qPFJet8CHSSD3Eta",     "std::vector<std::float>",       qPFJet8CHSSD3Eta_);
//   outTree_->Branch("qPFJet8CHSSD4Eta",     "std::vector<std::float>",       qPFJet8CHSSD4Eta_);
//   outTree_->Branch("qPFJet8CHSSD5Eta",     "std::vector<std::float>",       qPFJet8CHSSD5Eta_);

//   outTree_->Branch("qPFJet8CHSSD0Phi",     "std::vector<std::float>",       qPFJet8CHSSD0Phi_);
//   outTree_->Branch("qPFJet8CHSSD1Phi",     "std::vector<std::float>",       qPFJet8CHSSD1Phi_);
//   outTree_->Branch("qPFJet8CHSSD2Phi",     "std::vector<std::float>",       qPFJet8CHSSD2Phi_);
//   outTree_->Branch("qPFJet8CHSSD3Phi",     "std::vector<std::float>",       qPFJet8CHSSD3Phi_);
//   outTree_->Branch("qPFJet8CHSSD4Phi",     "std::vector<std::float>",       qPFJet8CHSSD4Phi_);
//   outTree_->Branch("qPFJet8CHSSD5Phi",     "std::vector<std::float>",       qPFJet8CHSSD5Phi_);

// //TopCHS sorted variables

//   outTree_->Branch("qPFJetTopCHS0Pt",    "std::vector<std::float>",       qPFJetTopCHS0Pt_);
//   outTree_->Branch("qPFJetTopCHS1Pt",    "std::vector<std::float>",       qPFJetTopCHS1Pt_);
//   outTree_->Branch("qPFJetTopCHS2Pt",    "std::vector<std::float>",       qPFJetTopCHS2Pt_);
//   outTree_->Branch("qPFJetTopCHS3Pt",    "std::vector<std::float>",       qPFJetTopCHS3Pt_);
//   outTree_->Branch("qPFJetTopCHS4Pt",    "std::vector<std::float>",       qPFJetTopCHS4Pt_);
//   outTree_->Branch("qPFJetTopCHS5Pt",    "std::vector<std::float>",       qPFJetTopCHS5Pt_);

//   outTree_->Branch("qPFJetTopCHS0Eta",     "std::vector<std::float>",       qPFJetTopCHS0Eta_);
//   outTree_->Branch("qPFJetTopCHS1Eta",     "std::vector<std::float>",       qPFJetTopCHS1Eta_);
//   outTree_->Branch("qPFJetTopCHS2Eta",     "std::vector<std::float>",       qPFJetTopCHS2Eta_);
//   outTree_->Branch("qPFJetTopCHS3Eta",     "std::vector<std::float>",       qPFJetTopCHS3Eta_);
//   outTree_->Branch("qPFJetTopCHS4Eta",     "std::vector<std::float>",       qPFJetTopCHS4Eta_);
//   outTree_->Branch("qPFJetTopCHS5Eta",     "std::vector<std::float>",       qPFJetTopCHS5Eta_);

//   outTree_->Branch("qPFJetTopCHS0Phi",     "std::vector<std::float>",       qPFJetTopCHS0Phi_);
//   outTree_->Branch("qPFJetTopCHS1Phi",     "std::vector<std::float>",       qPFJetTopCHS1Phi_);
//   outTree_->Branch("qPFJetTopCHS2Phi",     "std::vector<std::float>",       qPFJetTopCHS2Phi_);
//   outTree_->Branch("qPFJetTopCHS3Phi",     "std::vector<std::float>",       qPFJetTopCHS3Phi_);
//   outTree_->Branch("qPFJetTopCHS4Phi",     "std::vector<std::float>",       qPFJetTopCHS4Phi_);
//   outTree_->Branch("qPFJetTopCHS5Phi",     "std::vector<std::float>",       qPFJetTopCHS5Phi_);


//CaloJet sorted variables

  outTree_->Branch("qCalJet0Pt",     "std::vector<std::float>",       qCalJet0Pt_);
  outTree_->Branch("qCalJet1Pt",     "std::vector<std::float>",       qCalJet1Pt_);
  outTree_->Branch("qCalJet2Pt",     "std::vector<std::float>",       qCalJet2Pt_);
  outTree_->Branch("qCalJet3Pt",     "std::vector<std::float>",       qCalJet3Pt_);
  outTree_->Branch("qCalJet4Pt",     "std::vector<std::float>",       qCalJet4Pt_);
  outTree_->Branch("qCalJet5Pt",     "std::vector<std::float>",       qCalJet5Pt_);

  outTree_->Branch("qCalJet0Eta",     "std::vector<std::float>",       qCalJet0Eta_);
  outTree_->Branch("qCalJet1Eta",     "std::vector<std::float>",       qCalJet1Eta_);
  outTree_->Branch("qCalJet2Eta",     "std::vector<std::float>",       qCalJet2Eta_);
  outTree_->Branch("qCalJet3Eta",     "std::vector<std::float>",       qCalJet3Eta_);
  outTree_->Branch("qCalJet4Eta",     "std::vector<std::float>",       qCalJet4Eta_);
  outTree_->Branch("qCalJet5Eta",     "std::vector<std::float>",       qCalJet5Eta_);

  outTree_->Branch("qCalJet0Phi",     "std::vector<std::float>",       qCalJet0Phi_);
  outTree_->Branch("qCalJet1Phi",     "std::vector<std::float>",       qCalJet1Phi_);
  outTree_->Branch("qCalJet2Phi",     "std::vector<std::float>",       qCalJet2Phi_);
  outTree_->Branch("qCalJet3Phi",     "std::vector<std::float>",       qCalJet3Phi_);
  outTree_->Branch("qCalJet4Phi",     "std::vector<std::float>",       qCalJet4Phi_);
  outTree_->Branch("qCalJet5Phi",     "std::vector<std::float>",       qCalJet5Phi_);

  outTree_->Branch("qCalJet0En",     "std::vector<std::float>",       qCalJet0En_);
  outTree_->Branch("qCalJet1En",     "std::vector<std::float>",       qCalJet1En_);
  outTree_->Branch("qCalJet2En",     "std::vector<std::float>",       qCalJet2En_);
  outTree_->Branch("qCalJet3En",     "std::vector<std::float>",       qCalJet3En_);
  outTree_->Branch("qCalJet4En",     "std::vector<std::float>",       qCalJet4En_);
  outTree_->Branch("qCalJet5En",     "std::vector<std::float>",       qCalJet5En_);


//photon sorted variables (not all of them)

  outTree_->Branch("qPho0Pt",     "std::vector<std::float>",       qPho0Pt_);
  outTree_->Branch("qPho1Pt",     "std::vector<std::float>",       qPho1Pt_);
  outTree_->Branch("qPho2Pt",     "std::vector<std::float>",       qPho2Pt_);
  outTree_->Branch("qPho3Pt",     "std::vector<std::float>",       qPho3Pt_);
  outTree_->Branch("qPho4Pt",     "std::vector<std::float>",       qPho4Pt_);
  outTree_->Branch("qPho5Pt",     "std::vector<std::float>",       qPho5Pt_);

  outTree_->Branch("qPho0Eta",     "std::vector<std::float>",       qPho0Eta_);
  outTree_->Branch("qPho1Eta",     "std::vector<std::float>",       qPho1Eta_);
  outTree_->Branch("qPho2Eta",     "std::vector<std::float>",       qPho2Eta_);
  outTree_->Branch("qPho3Eta",     "std::vector<std::float>",       qPho3Eta_);
  outTree_->Branch("qPho4Eta",     "std::vector<std::float>",       qPho4Eta_);
  outTree_->Branch("qPho5Eta",     "std::vector<std::float>",       qPho5Eta_);

  outTree_->Branch("qPho0Phi",     "std::vector<std::float>",       qPho0Phi_);
  outTree_->Branch("qPho1Phi",     "std::vector<std::float>",       qPho1Phi_);
  outTree_->Branch("qPho2Phi",     "std::vector<std::float>",       qPho2Phi_);
  outTree_->Branch("qPho3Phi",     "std::vector<std::float>",       qPho3Phi_);
  outTree_->Branch("qPho4Phi",     "std::vector<std::float>",       qPho4Phi_);
  outTree_->Branch("qPho5Phi",     "std::vector<std::float>",       qPho5Phi_);

  outTree_->Branch("qPho0En",     "std::vector<std::float>",       qPho0En_);
  outTree_->Branch("qPho1En",     "std::vector<std::float>",       qPho1En_);
  outTree_->Branch("qPho2En",     "std::vector<std::float>",       qPho2En_);
  outTree_->Branch("qPho3En",     "std::vector<std::float>",       qPho3En_);
  outTree_->Branch("qPho4En",     "std::vector<std::float>",       qPho4En_);
  outTree_->Branch("qPho5En",     "std::vector<std::float>",       qPho5En_);

//ged qPhotons sorted variables (not all of them)

  outTree_->Branch("qgedPho0Pt",     "std::vector<std::float>",       qgedPho0Pt_);
  outTree_->Branch("qgedPho1Pt",     "std::vector<std::float>",       qgedPho1Pt_);
  outTree_->Branch("qgedPho2Pt",     "std::vector<std::float>",       qgedPho2Pt_);
  outTree_->Branch("qgedPho3Pt",     "std::vector<std::float>",       qgedPho3Pt_);
  outTree_->Branch("qgedPho4Pt",     "std::vector<std::float>",       qgedPho4Pt_);
  outTree_->Branch("qgedPho5Pt",     "std::vector<std::float>",       qgedPho5Pt_);
        
  outTree_->Branch("qgedPho0Eta",     "std::vector<std::float>",       qgedPho0Eta_);
  outTree_->Branch("qgedPho1Eta",     "std::vector<std::float>",       qgedPho1Eta_);
  outTree_->Branch("qgedPho2Eta",     "std::vector<std::float>",       qgedPho2Eta_);
  outTree_->Branch("qgedPho3Eta",     "std::vector<std::float>",       qgedPho3Eta_);
  outTree_->Branch("qgedPho4Eta",     "std::vector<std::float>",       qgedPho4Eta_);
  outTree_->Branch("qgedPho5Eta",     "std::vector<std::float>",       qgedPho5Eta_);

  outTree_->Branch("qgedPho0Phi",     "std::vector<std::float>",       qgedPho0Phi_);
  outTree_->Branch("qgedPho1Phi",     "std::vector<std::float>",       qgedPho1Phi_);
  outTree_->Branch("qgedPho2Phi",     "std::vector<std::float>",       qgedPho2Phi_);
  outTree_->Branch("qgedPho3Phi",     "std::vector<std::float>",       qgedPho3Phi_);
  outTree_->Branch("qgedPho4Phi",     "std::vector<std::float>",       qgedPho4Phi_);
  outTree_->Branch("qgedPho5Phi",     "std::vector<std::float>",       qgedPho5Phi_);

  outTree_->Branch("qgedPho0En",     "std::vector<std::float>",       qgedPho0En_);
  outTree_->Branch("qgedPho1En",     "std::vector<std::float>",       qgedPho1En_);
  outTree_->Branch("qgedPho2En",     "std::vector<std::float>",       qgedPho2En_);
  outTree_->Branch("qgedPho3En",     "std::vector<std::float>",       qgedPho3En_);
  outTree_->Branch("qgedPho4En",     "std::vector<std::float>",       qgedPho4En_);
  outTree_->Branch("qgedPho5En",     "std::vector<std::float>",       qgedPho5En_);


//Muon sorted variables (not all of them)

  outTree_->Branch("qMu0Pt",     "std::vector<std::float>",       qMu0Pt_);
  outTree_->Branch("qMu1Pt",     "std::vector<std::float>",       qMu1Pt_);
  outTree_->Branch("qMu2Pt",     "std::vector<std::float>",       qMu2Pt_);
  outTree_->Branch("qMu3Pt",     "std::vector<std::float>",       qMu3Pt_);
  outTree_->Branch("qMu4Pt",     "std::vector<std::float>",       qMu4Pt_);
  outTree_->Branch("qMu5Pt",     "std::vector<std::float>",       qMu5Pt_);

  outTree_->Branch("qMu0Eta",     "std::vector<std::float>",       qMu0Eta_);
  outTree_->Branch("qMu1Eta",     "std::vector<std::float>",       qMu1Eta_);
  outTree_->Branch("qMu2Eta",     "std::vector<std::float>",       qMu2Eta_);
  outTree_->Branch("qMu3Eta",     "std::vector<std::float>",       qMu3Eta_);
  outTree_->Branch("qMu4Eta",     "std::vector<std::float>",       qMu4Eta_);
  outTree_->Branch("qMu5Eta",     "std::vector<std::float>",       qMu5Eta_);

  outTree_->Branch("qMu0Phi",     "std::vector<std::float>",       qMu0Phi_);
  outTree_->Branch("qMu1Phi",     "std::vector<std::float>",       qMu1Phi_);
  outTree_->Branch("qMu2Phi",     "std::vector<std::float>",       qMu2Phi_);
  outTree_->Branch("qMu3Phi",     "std::vector<std::float>",       qMu3Phi_);
  outTree_->Branch("qMu4Phi",     "std::vector<std::float>",       qMu4Phi_);
  outTree_->Branch("qMu5Phi",     "std::vector<std::float>",       qMu5Phi_);

  outTree_->Branch("qMu0En",     "std::vector<std::float>",       qMu0En_);
  outTree_->Branch("qMu1En",     "std::vector<std::float>",       qMu1En_);
  outTree_->Branch("qMu2En",     "std::vector<std::float>",       qMu2En_);
  outTree_->Branch("qMu3En",     "std::vector<std::float>",       qMu3En_);
  outTree_->Branch("qMu4En",     "std::vector<std::float>",       qMu4En_);
  outTree_->Branch("qMu5En",     "std::vector<std::float>",       qMu5En_);

//Muon sorted cosmic variables (not all of them)

  outTree_->Branch("qMuCosm0Pt",     "std::vector<std::float>",       qMuCosm0Pt_);
  outTree_->Branch("qMuCosm1Pt",     "std::vector<std::float>",       qMuCosm1Pt_);
  outTree_->Branch("qMuCosm2Pt",     "std::vector<std::float>",       qMuCosm2Pt_);
  outTree_->Branch("qMuCosm3Pt",     "std::vector<std::float>",       qMuCosm3Pt_);
  outTree_->Branch("qMuCosm4Pt",     "std::vector<std::float>",       qMuCosm4Pt_);
  outTree_->Branch("qMuCosm5Pt",     "std::vector<std::float>",       qMuCosm5Pt_);

  outTree_->Branch("qMuCosm0Eta",     "std::vector<std::float>",       qMuCosm0Eta_);
  outTree_->Branch("qMuCosm1Eta",     "std::vector<std::float>",       qMuCosm1Eta_);
  outTree_->Branch("qMuCosm2Eta",     "std::vector<std::float>",       qMuCosm2Eta_);
  outTree_->Branch("qMuCosm3Eta",     "std::vector<std::float>",       qMuCosm3Eta_);
  outTree_->Branch("qMuCosm4Eta",     "std::vector<std::float>",       qMuCosm4Eta_);
  outTree_->Branch("qMuCosm5Eta",     "std::vector<std::float>",       qMuCosm5Eta_);

  outTree_->Branch("qMuCosm0Phi",     "std::vector<std::float>",       qMuCosm0Phi_);
  outTree_->Branch("qMuCosm1Phi",     "std::vector<std::float>",       qMuCosm1Phi_);
  outTree_->Branch("qMuCosm2Phi",     "std::vector<std::float>",       qMuCosm2Phi_);
  outTree_->Branch("qMuCosm3Phi",     "std::vector<std::float>",       qMuCosm3Phi_);
  outTree_->Branch("qMuCosm4Phi",     "std::vector<std::float>",       qMuCosm4Phi_);
  outTree_->Branch("qMuCosm5Phi",     "std::vector<std::float>",       qMuCosm5Phi_);

  outTree_->Branch("qMuCosm0En",     "std::vector<std::float>",       qMuCosm0En_);
  outTree_->Branch("qMuCosm1En",     "std::vector<std::float>",       qMuCosm1En_);
  outTree_->Branch("qMuCosm2En",     "std::vector<std::float>",       qMuCosm2En_);
  outTree_->Branch("qMuCosm3En",     "std::vector<std::float>",       qMuCosm3En_);
  outTree_->Branch("qMuCosm4En",     "std::vector<std::float>",       qMuCosm4En_);
  outTree_->Branch("qMuCosm5En",     "std::vector<std::float>",       qMuCosm5En_);


 //Muon Cosmic1Leg sorted variables (not all of them)

  outTree_->Branch("qMuCosmLeg0Pt",     "std::vector<std::float>",       qMuCosmLeg0Pt_);
  outTree_->Branch("qMuCosmLeg1Pt",     "std::vector<std::float>",       qMuCosmLeg1Pt_);
  outTree_->Branch("qMuCosmLeg2Pt",     "std::vector<std::float>",       qMuCosmLeg2Pt_);
  outTree_->Branch("qMuCosmLeg3Pt",     "std::vector<std::float>",       qMuCosmLeg3Pt_);
  outTree_->Branch("qMuCosmLeg4Pt",     "std::vector<std::float>",       qMuCosmLeg4Pt_);
  outTree_->Branch("qMuCosmLeg5Pt",     "std::vector<std::float>",       qMuCosmLeg5Pt_);

  outTree_->Branch("qMuCosmLeg0Eta",     "std::vector<std::float>",       qMuCosmLeg0Eta_);
  outTree_->Branch("qMuCosmLeg1Eta",     "std::vector<std::float>",       qMuCosmLeg1Eta_);
  outTree_->Branch("qMuCosmLeg2Eta",     "std::vector<std::float>",       qMuCosmLeg2Eta_);
  outTree_->Branch("qMuCosmLeg3Eta",     "std::vector<std::float>",       qMuCosmLeg3Eta_);
  outTree_->Branch("qMuCosmLeg4Eta",     "std::vector<std::float>",       qMuCosmLeg4Eta_);
  outTree_->Branch("qMuCosmLeg5Eta",     "std::vector<std::float>",       qMuCosmLeg5Eta_);

  outTree_->Branch("qMuCosmLeg0Phi",     "std::vector<std::float>",       qMuCosmLeg0Phi_);
  outTree_->Branch("qMuCosmLeg1Phi",     "std::vector<std::float>",       qMuCosmLeg1Phi_);
  outTree_->Branch("qMuCosmLeg2Phi",     "std::vector<std::float>",       qMuCosmLeg2Phi_);
  outTree_->Branch("qMuCosmLeg3Phi",     "std::vector<std::float>",       qMuCosmLeg3Phi_);
  outTree_->Branch("qMuCosmLeg4Phi",     "std::vector<std::float>",       qMuCosmLeg4Phi_);
  outTree_->Branch("qMuCosmLeg5Phi",     "std::vector<std::float>",       qMuCosmLeg5Phi_);

  outTree_->Branch("qMuCosmLeg0En",     "std::vector<std::float>",       qMuCosmLeg0En_);
  outTree_->Branch("qMuCosmLeg1En",     "std::vector<std::float>",       qMuCosmLeg1En_);
  outTree_->Branch("qMuCosmLeg2En",     "std::vector<std::float>",       qMuCosmLeg2En_);
  outTree_->Branch("qMuCosmLeg3En",     "std::vector<std::float>",       qMuCosmLeg3En_);
  outTree_->Branch("qMuCosmLeg4En",     "std::vector<std::float>",       qMuCosmLeg4En_);
  outTree_->Branch("qMuCosmLeg5En",     "std::vector<std::float>",       qMuCosmLeg5En_);





//___________________________________________Sorted variables end here.

  outTree_->Branch("qPFJet4CHSPt",     "std::vector<std::float>",      &qPFJet4CHSPt_);
  outTree_->Branch("qPFJet4CHSEta",    "std::vector<std::float>",      &qPFJet4CHSEta_);
  outTree_->Branch("qPFJet4CHSPhi",    "std::vector<std::float>",      &qPFJet4CHSPhi_);

  outTree_->Branch("qPFJet8CHSPt",     "std::vector<std::float>",      &qPFJet8CHSPt_);
  outTree_->Branch("qPFJet8CHSEta",    "std::vector<std::float>",      &qPFJet8CHSEta_);
  outTree_->Branch("qPFJet8CHSPhi",    "std::vector<std::float>",      &qPFJet8CHSPhi_);

  outTree_->Branch("qPFJetEIPt",     "std::vector<std::float>",      &qPFJetEIPt_);
  outTree_->Branch("qPFJetEIEta",    "std::vector<std::float>",      &qPFJetEIEta_);
  outTree_->Branch("qPFJetEIPhi",    "std::vector<std::float>",      &qPFJetEIPhi_);

  outTree_->Branch("qPFJet8CHSSDPt",     "std::vector<std::float>",      &qPFJet8CHSSDPt_);
  outTree_->Branch("qPFJet8CHSSDEta",    "std::vector<std::float>",      &qPFJet8CHSSDEta_);
  outTree_->Branch("qPFJet8CHSSDPhi",    "std::vector<std::float>",      &qPFJet8CHSSDPhi_);

  outTree_->Branch("qPFJetTopCHSPt",     "std::vector<std::float>",      &qPFJetTopCHSPt_);
  outTree_->Branch("qPFJetTopCHSEta",    "std::vector<std::float>",      &qPFJetTopCHSEta_);
  outTree_->Branch("qPFJetTopCHSPhi",    "std::vector<std::float>",      &qPFJetTopCHSPhi_);

  outTree_->Branch("qPFChMetPt",     "std::vector<std::float>",        &qPFChMetPt_);
  // outTree_->Branch("qPFChMetEta",    "std::vector<std::float>",        &qPFChMetEta_);  
  outTree_->Branch("qPFChMetPhi",    "std::vector<std::float>",        &qPFChMetPhi_);
  outTree_->Branch("qPFMetPt",     "std::vector<std::float>",        &qPFMetPt_);
  // outTree_->Branch("qPFMetEta",    "std::vector<std::float>",        &qPFMetEta_);
  outTree_->Branch("qPFMetPhi",    "std::vector<std::float>",        &qPFMetPhi_);
  outTree_->Branch("qNVtx",        "std::vector<std::int>",        &qNVtx_);

  outTree_->Branch("qCalJetPt",     "std::vector<std::float>",        &qCalJetPt_);
  outTree_->Branch("qCalJetEta",    "std::vector<std::float>",        &qCalJetEta_);
  outTree_->Branch("qCalJetPhi",    "std::vector<std::float>",        &qCalJetPhi_);
  outTree_->Branch("qCalJetEn",    "std::vector<std::float>",        &qCalJetEn_);

  outTree_->Branch("qCalMETPt",     "std::vector<std::float>",        &qCalMETPt_);
  // outTree_->Branch("qCalMETEta",    "std::vector<std::float>",        &qCalMETEta_);
  outTree_->Branch("qCalMETPhi",    "std::vector<std::float>",        &qCalMETPhi_);
  outTree_->Branch("qCalMETEn",    "std::vector<std::float>",        &qCalMETEn_);

  outTree_->Branch("qCalMETBEPt",     "std::vector<std::float>",        &qCalMETBEPt_);
  // outTree_->Branch("qCalMETBEEta",    "std::vector<std::float>",        &qCalMETBEEta_);
  outTree_->Branch("qCalMETBEPhi",    "std::vector<std::float>",        &qCalMETBEPhi_);
  outTree_->Branch("qCalMETBEEn",    "std::vector<std::float>",         &qCalMETBEEn_);

  outTree_->Branch("qCalMETBEFOPt",     "std::vector<std::float>",        &qCalMETBEFOPt_);
  // outTree_->Branch("qCalMETBEFOEta",    "std::vector<std::float>",        &qCalMETBEFOEta_);
  outTree_->Branch("qCalMETBEFOPhi",    "std::vector<std::float>",        &qCalMETBEFOPhi_);
  outTree_->Branch("qCalMETBEFOEn",    "std::vector<std::float>",         &qCalMETBEFOEn_);

  outTree_->Branch("qCalMETMPt",     "std::vector<std::float>",        &qCalMETMPt_);
  // outTree_->Branch("qCalMETMEta",    "std::vector<std::float>",        &qCalMETMEta_);
  outTree_->Branch("qCalMETMPhi",    "std::vector<std::float>",        &qCalMETMPhi_);
  outTree_->Branch("qCalMETMEn",    "std::vector<std::float>",         &qCalMETMEn_);

  outTree_->Branch("qSCEn",     "std::vector<std::float>",        &qSCEn_);
  outTree_->Branch("qSCEta",    "std::vector<std::float>",        &qSCEta_);
  outTree_->Branch("qSCPhi",    "std::vector<std::float>",        &qSCPhi_);
  outTree_->Branch("qSCEtaWidth",    "std::vector<std::float>",        &qSCEtaWidth_);
  outTree_->Branch("qSCPhiWidth",    "std::vector<std::float>",        &qSCPhiWidth_);  
  outTree_->Branch("qSCEnhfEM",     "std::vector<std::float>",        &qSCEnhfEM_);
  outTree_->Branch("qSCEtahfEM",    "std::vector<std::float>",        &qSCEtahfEM_);
  outTree_->Branch("qSCPhihfEM",    "std::vector<std::float>",        &qSCPhihfEM_);
  // outTree_->Branch("qSCEtaWidthhfEM",    "std::vector<std::float>",        &qSCEtaWidthhfEM_);
  // outTree_->Branch("qSCPhiWidthhfEM",    "std::vector<std::float>",        &qSCPhiWidthhfEM_);  
  outTree_->Branch("qSCEn5x5",     "std::vector<std::float>",        &qSCEn5x5_);
  outTree_->Branch("qSCEta5x5",    "std::vector<std::float>",        &qSCEta5x5_);
  outTree_->Branch("qSCPhi5x5",    "std::vector<std::float>",        &qSCPhi5x5_);
  outTree_->Branch("qSCEtaWidth5x5",    "std::vector<std::float>",        &qSCEtaWidth5x5_);
  outTree_->Branch("qSCPhiWidth5x5",    "std::vector<std::float>",        &qSCPhiWidth5x5_);  
  outTree_->Branch("qCCEn",     "std::vector<std::float>",        &qCCEn_);
  outTree_->Branch("qCCEta",    "std::vector<std::float>",        &qCCEta_);
  outTree_->Branch("qCCPhi",    "std::vector<std::float>",        &qCCPhi_);
  outTree_->Branch("qCCEn5x5",     "std::vector<std::float>",        &qCCEn5x5_);
  outTree_->Branch("qCCEta5x5",    "std::vector<std::float>",        &qCCEta5x5_);
  outTree_->Branch("qCCPhi5x5",    "std::vector<std::float>",        &qCCPhi5x5_);

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
  outTree_->Branch("qMuChi2_",    "std::vector<std::float>",        &qMuChi2_);  
  outTree_->Branch("qMuCosmPt",     "std::vector<std::float>",        &qMuCosmPt_);
  outTree_->Branch("qMuCosmEta",    "std::vector<std::float>",        &qMuCosmEta_);
  outTree_->Branch("qMuCosmPhi",    "std::vector<std::float>",        &qMuCosmPhi_);
  outTree_->Branch("qMuCosmEn_",    "std::vector<std::float>",        &qMuCosmEn_);
  outTree_->Branch("qMuCosmCh_",    "std::vector<std::float>",        &qMuCosmCh_);
  outTree_->Branch("qMuCosmChi2_",    "std::vector<std::float>",        &qMuCosmChi2_);    
  outTree_->Branch("qMuCosmLegPt",     "std::vector<std::float>",        &qMuCosmLegPt_);
  outTree_->Branch("qMuCosmLegEta",    "std::vector<std::float>",        &qMuCosmLegEta_);
  outTree_->Branch("qMuCosmLegPhi",    "std::vector<std::float>",        &qMuCosmLegPhi_);
  outTree_->Branch("qMuCosmLegEn_",    "std::vector<std::float>",        &qMuCosmLegEn_);
  outTree_->Branch("qMuCosmLegCh_",    "std::vector<std::float>",        &qMuCosmLegCh_);
  outTree_->Branch("qMuCosmLegChi2_",    "std::vector<std::float>",        &qMuCosmLegChi2_);   
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
  outTree_->Branch("qEBiEta",    "std::vector<std::float>",          &qEBiEta_);
  outTree_->Branch("qEBiPhi",    "std::vector<std::float>",          &qEBiPhi_);

  outTree_->Branch("qEEenergy",    "std::vector<std::float>",        &qEEenergy_);
  outTree_->Branch("qEEtime",    "std::vector<std::float>",          &qEEtime_);
  outTree_->Branch("qEEchi2",    "std::vector<std::float>",          &qEEchi2_);
  outTree_->Branch("qEEix",    "std::vector<std::float>",          &qEEix_);
  outTree_->Branch("qEEiy",    "std::vector<std::float>",          &qEEiy_);

  outTree_->Branch("qESenergy",    "std::vector<std::float>",        &qESenergy_);
  outTree_->Branch("qEStime",    "std::vector<std::float>",          &qEStime_);
  // outTree_->Branch("qESchi2",    "std::vector<std::float>",          &qESchi2_);
  outTree_->Branch("qESix",    "std::vector<std::float>",          &qESix_);
  outTree_->Branch("qESiy",    "std::vector<std::float>",          &qESiy_);

  outTree_->Branch("qHBHEenergy",    "std::vector<std::float>",        &qHBHEenergy_);
  outTree_->Branch("qHBHEtime",    "std::vector<std::float>",          &qHBHEtime_);
  outTree_->Branch("qHBHEauxe",    "std::vector<std::float>",          &qHBHEauxe_);
  outTree_->Branch("qHBHEieta",    "std::vector<std::float>",        &qHBHEieta_);
  outTree_->Branch("qHBHEiphi",    "std::vector<std::float>",          &qHBHEiphi_);

  outTree_->Branch("qHFenergy",    "std::vector<std::float>",        &qHFenergy_);
  outTree_->Branch("qHFtime",    "std::vector<std::float>",          &qHFtime_);
  // outTree_->Branch("qHFchi2",    "std::vector<std::float>",          &qHFchi2_);
  outTree_->Branch("qHFieta",    "std::vector<std::float>",          &qHFieta_);
  outTree_->Branch("qHFiphi",    "std::vector<std::float>",          &qHFiphi_);

  // outTree_->Branch("qHOenergy",    "std::vector<std::float>",        &qHOenergy_);
  // outTree_->Branch("qHOtime",    "std::vector<std::float>",          &qHOtime_);
  // // outTree_->Branch("qHOchi2",    "std::vector<std::float>",          &qHOchi2_);
  // outTree_->Branch("qHOieta",    "std::vector<std::float>",          &qHOieta_);
  // outTree_->Branch("qHOiphi",    "std::vector<std::float>",          &qHOiphi_);

  outTree_->Branch("qPreShEn",     "std::vector<std::float>",        &qPreShEn_);
  // outTree_->Branch("qPreShCorrEn", "std::vector<std::float>",        &qPreShCorrEn_);
  outTree_->Branch("qPreShEta",    "std::vector<std::float>",        &qPreShEta_);
  outTree_->Branch("qPreShPhi",    "std::vector<std::float>",        &qPreShPhi_);
  outTree_->Branch("qPreShYEn",     "std::vector<std::float>",        &qPreShYEn_);
  // outTree_->Branch("qPreShYCorrEn", "std::vector<std::float>",        &qPreShYCorrEn_);
  outTree_->Branch("qPreShYEta",    "std::vector<std::float>",        &qPreShYEta_);
  outTree_->Branch("qPreShYPhi",    "std::vector<std::float>",        &qPreShYPhi_);

  // outTree_->Branch("qCTPt",     "std::vector<std::float>",        &qCTPt_);
  // outTree_->Branch("qCTEta",    "std::vector<std::float>",        &qCTEta_);
  // outTree_->Branch("qCTPhi",    "std::vector<std::float>",        &qCTPhi_);

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

 //-------------------------Sorted variables start here.

//PFJet sorted variables
delete PFJet0Pt_;
delete PFJet1Pt_; 
delete PFJet2Pt_; 
delete PFJet3Pt_; 
delete PFJet4Pt_; 
delete PFJet5Pt_; 

delete PFJet0Eta_;
delete PFJet1Eta_;
delete PFJet2Eta_;
delete PFJet3Eta_;
delete PFJet4Eta_;
delete PFJet5Eta_;

delete PFJet0Phi_;
delete PFJet1Phi_;
delete PFJet2Phi_;
delete PFJet3Phi_;
delete PFJet4Phi_;
delete PFJet5Phi_;
//PFJet4CHS sorted variables
delete PFJet4CHS0Pt_;
delete PFJet4CHS1Pt_;
delete PFJet4CHS2Pt_;
delete PFJet4CHS3Pt_;
delete PFJet4CHS4Pt_;
delete PFJet4CHS5Pt_;

delete PFJet4CHS0Eta_;
delete PFJet4CHS1Eta_;
delete PFJet4CHS2Eta_;
delete PFJet4CHS3Eta_;
delete PFJet4CHS4Eta_;
delete PFJet4CHS5Eta_;

delete PFJet4CHS0Phi_;
delete PFJet4CHS1Phi_;
delete PFJet4CHS2Phi_;
delete PFJet4CHS3Phi_;
delete PFJet4CHS4Phi_;
delete PFJet4CHS5Phi_;
 //PF8CHS sorted variables
delete PFJet8CHS0Pt_;
delete PFJet8CHS1Pt_;
delete PFJet8CHS2Pt_;
delete PFJet8CHS3Pt_;
delete PFJet8CHS4Pt_;
delete PFJet8CHS5Pt_;

delete PFJet8CHS0Eta_;
delete PFJet8CHS1Eta_;
delete PFJet8CHS2Eta_;
delete PFJet8CHS3Eta_;
delete PFJet8CHS4Eta_;
delete PFJet8CHS5Eta_;

delete PFJet8CHS0Phi_;
delete PFJet8CHS1Phi_;
delete PFJet8CHS2Phi_;
delete PFJet8CHS3Phi_;
delete PFJet8CHS4Phi_;
delete PFJet8CHS5Phi_;
  //PFJetEI sorted variables
delete PFJetEI0Pt_;
delete PFJetEI1Pt_;
delete PFJetEI2Pt_;
delete PFJetEI3Pt_;
delete PFJetEI4Pt_;
delete PFJetEI5Pt_;

delete PFJetEI0Eta_;
delete PFJetEI1Eta_;
delete PFJetEI2Eta_;
delete PFJetEI3Eta_;
delete PFJetEI4Eta_;
delete PFJetEI5Eta_;

delete PFJetEI0Phi_;
delete PFJetEI1Phi_;
delete PFJetEI2Phi_;
delete PFJetEI3Phi_;
delete PFJetEI4Phi_;
delete PFJetEI5Phi_;

//   //8CHSSoftDrop sorted variables
// delete PFJet8CHSSD0Pt_;
// delete PFJet8CHSSD1Pt_;
// delete PFJet8CHSSD2Pt_;
// delete PFJet8CHSSD3Pt_;
// delete PFJet8CHSSD4Pt_;
// delete PFJet8CHSSD5Pt_;

// delete PFJet8CHSSD0Eta_;
// delete PFJet8CHSSD1Eta_;
// delete PFJet8CHSSD2Eta_;
// delete PFJet8CHSSD3Eta_;
// delete PFJet8CHSSD4Eta_;
// delete PFJet8CHSSD5Eta_;

// delete PFJet8CHSSD0Phi_;
// delete PFJet8CHSSD1Phi_;
// delete PFJet8CHSSD2Phi_;
// delete PFJet8CHSSD3Phi_;
// delete PFJet8CHSSD4Phi_;
// delete PFJet8CHSSD5Phi_;
//   //TopCHS sorted variables
// delete PFJetTopCHS0Pt_;
// delete PFJetTopCHS1Pt_;
// delete PFJetTopCHS2Pt_;
// delete PFJetTopCHS3Pt_;
// delete PFJetTopCHS4Pt_;
// delete PFJetTopCHS5Pt_;

// delete PFJetTopCHS0Eta_;
// delete PFJetTopCHS1Eta_;
// delete PFJetTopCHS2Eta_;
// delete PFJetTopCHS3Eta_;
// delete PFJetTopCHS4Eta_;
// delete PFJetTopCHS5Eta_;

// delete PFJetTopCHS0Phi_;
// delete PFJetTopCHS1Phi_;
// delete PFJetTopCHS2Phi_;
// delete PFJetTopCHS3Phi_;
// delete PFJetTopCHS4Phi_;
// delete PFJetTopCHS5Phi_;


  //CaloJet sorted variables
delete CalJet0Pt_;
delete CalJet1Pt_;
delete CalJet2Pt_;
delete CalJet3Pt_;
delete CalJet4Pt_;
delete CalJet5Pt_;

delete CalJet0Eta_;
delete CalJet1Eta_;
delete CalJet2Eta_;
delete CalJet3Eta_;
delete CalJet4Eta_;
delete CalJet5Eta_;

delete CalJet0Phi_;
delete CalJet1Phi_;
delete CalJet2Phi_;
delete CalJet3Phi_;
delete CalJet4Phi_;
delete CalJet5Phi_;

delete CalJet0En_;
delete CalJet1En_;
delete CalJet2En_;
delete CalJet3En_;
delete CalJet4En_;
delete CalJet5En_;


  //photon sorted variables (not all of them)
delete Pho0Pt_;
delete Pho1Pt_;
delete Pho2Pt_;
delete Pho3Pt_;
delete Pho4Pt_;
delete Pho5Pt_;

delete Pho0Eta_;
delete Pho1Eta_;
delete Pho2Eta_;
delete Pho3Eta_;
delete Pho4Eta_;
delete Pho5Eta_;

delete Pho0Phi_;
delete Pho1Phi_;
delete Pho2Phi_;
delete Pho3Phi_;
delete Pho4Phi_;
delete Pho5Phi_;

delete Pho0En_;
delete Pho1En_;
delete Pho2En_;
delete Pho3En_;
delete Pho4En_;
delete Pho5En_;

  //ged qPhotons sorted variables (not all of them)
delete gedPho0Pt_;
delete gedPho1Pt_;
delete gedPho2Pt_;
delete gedPho3Pt_;
delete gedPho4Pt_;
delete gedPho5Pt_;
 
delete gedPho0Eta_;
delete gedPho1Eta_;
delete gedPho2Eta_;
delete gedPho3Eta_;
delete gedPho4Eta_;
delete gedPho5Eta_;

delete gedPho0Phi_;
delete gedPho1Phi_;
delete gedPho2Phi_;
delete gedPho3Phi_;
delete gedPho4Phi_;
delete gedPho5Phi_;

delete gedPho0En_;
delete gedPho1En_;
delete gedPho2En_;
delete gedPho3En_;
delete gedPho4En_;
delete gedPho5En_;


    //Muon sorted variables (not all of them)
delete Mu0Pt_;
delete Mu1Pt_;
delete Mu2Pt_;
delete Mu3Pt_;
delete Mu4Pt_;
delete Mu5Pt_;

delete Mu0Eta_;
delete Mu1Eta_;
delete Mu2Eta_;
delete Mu3Eta_;
delete Mu4Eta_;
delete Mu5Eta_;

delete Mu0Phi_;
delete Mu1Phi_;
delete Mu2Phi_;
delete Mu3Phi_;
delete Mu4Phi_;
delete Mu5Phi_;

delete Mu0En_;
delete Mu1En_;
delete Mu2En_;
delete Mu3En_;
delete Mu4En_;
delete Mu5En_;
  //Muon sorted cosmic variables (not all of them)
delete MuCosm0Pt_;
delete MuCosm1Pt_;
delete MuCosm2Pt_;
delete MuCosm3Pt_;
delete MuCosm4Pt_;
delete MuCosm5Pt_;

delete MuCosm0Eta_;
delete MuCosm1Eta_;
delete MuCosm2Eta_;
delete MuCosm3Eta_;
delete MuCosm4Eta_;
delete MuCosm5Eta_;

delete MuCosm0Phi_;
delete MuCosm1Phi_;
delete MuCosm2Phi_;
delete MuCosm3Phi_;
delete MuCosm4Phi_;
delete MuCosm5Phi_;

delete MuCosm0En_;
delete MuCosm1En_;
delete MuCosm2En_;
delete MuCosm3En_;
delete MuCosm4En_;
delete MuCosm5En_;


  //Muon Cosmic1Leg sorted variables (not all of them)
delete MuCosmLeg0Pt_;
delete MuCosmLeg1Pt_;
delete MuCosmLeg2Pt_;
delete MuCosmLeg3Pt_;
delete MuCosmLeg4Pt_;
delete MuCosmLeg5Pt_;

delete MuCosmLeg0Eta_;
delete MuCosmLeg1Eta_;
delete MuCosmLeg2Eta_;
delete MuCosmLeg3Eta_;
delete MuCosmLeg4Eta_;
delete MuCosmLeg5Eta_;

delete MuCosmLeg0Phi_;
delete MuCosmLeg1Phi_;
delete MuCosmLeg2Phi_;
delete MuCosmLeg3Phi_;
delete MuCosmLeg4Phi_;
delete MuCosmLeg5Phi_;

delete MuCosmLeg0En_;
delete MuCosmLeg1En_;
delete MuCosmLeg2En_;
delete MuCosmLeg3En_;
delete MuCosmLeg4En_;
delete MuCosmLeg5En_;


 //_________________________Sorted variables end here. 

  delete PFJet4CHSPt_;
  delete PFJet4CHSEta_;
  delete PFJet4CHSPhi_;

  delete PFJet8CHSPt_;
  delete PFJet8CHSEta_;
  delete PFJet8CHSPhi_;

  delete PFJetEIPt_;
  delete PFJetEIEta_;
  delete PFJetEIPhi_;

  delete PFJet8CHSSDPt_;
  delete PFJet8CHSSDEta_;
  delete PFJet8CHSSDPhi_;

  delete PFJetTopCHSPt_;
  delete PFJetTopCHSEta_;
  delete PFJetTopCHSPhi_;

  delete PFChMetPt_;
  // delete PFChMetEta_;  
  delete PFChMetPhi_;
  delete PFMetPt_;
  // delete PFMetEta_;
  delete PFMetPhi_;

  delete CalJetPt_;
  delete CalJetEta_;
  delete CalJetPhi_;
  delete CalJetEn_;
  delete CalMETPt_;
  // delete CalMETEta_;
  delete CalMETPhi_;
  delete CalMETEn_;

  delete CalMETBEPt_;
  // delete CalMETBEEta_;
  delete CalMETBEPhi_;
  delete CalMETBEEn_;
  delete CalMETBEFOPt_;
  // delete CalMETBEFOEta_;
  delete CalMETBEFOPhi_;
  delete CalMETBEFOEn_;
  delete CalMETMPt_;
  // delete CalMETMEta_;
  delete CalMETMPhi_;
  delete CalMETMEn_;


  delete SCEn_;
  delete SCEta_;
  delete SCPhi_;
  delete SCEtaWidth_;
  delete SCPhiWidth_;  
  delete SCEnhfEM_;
  delete SCEtahfEM_;
  delete SCPhihfEM_;
  // delete SCEtaWidthhfEM_;
  // delete SCPhiWidthhfEM_; 
  delete SCEn5x5_;
  delete SCEta5x5_;
  delete SCPhi5x5_;
  delete SCEtaWidth5x5_;
  delete SCPhiWidth5x5_;  
  delete CCEn_;
  delete CCEta_;
  delete CCPhi_;
  delete CCEn5x5_;
  delete CCEta5x5_;
  delete CCPhi5x5_;

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
  delete MuChi2_;  
  delete MuCosmPt_;
  delete MuCosmEta_;
  delete MuCosmPhi_;
  delete MuCosmEn_;
  delete MuCosmCh_;
  delete MuCosmChi2_;    
  delete MuCosmLegPt_;
  delete MuCosmLegEta_;
  delete MuCosmLegPhi_;
  delete MuCosmLegEn_;
  delete MuCosmLegCh_;    
  delete MuCosmLegChi2_;  
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
  delete EBiEta_;
  delete EBiPhi_;
  delete EEenergy_;
  delete EEtime_;
  delete EEchi2_;
  delete EEix_;
  delete EEiy_;
  delete ESenergy_;
  delete EStime_;
  // delete ESchi2_;
  delete ESix_;
  delete ESiy_;

  delete HBHEenergy_ ;
  delete HBHEtime_ ;
  delete HBHEauxe_ ;
  delete HBHEieta_;
  delete HBHEiphi_;  
  delete HFenergy_ ;
  delete HFtime_ ;
  //delete HFchi2_  ;
  delete HFieta_;
  delete HFiphi_;
  // delete HOenergy_ ;
  // delete HOtime_;
  // delete HOieta_;
  // delete HOiphi_;
  // HOchi2_    
  delete PreShEn_;
  // delete PreShCorrEn_;
  delete PreShEta_;
  delete PreShPhi_;
  delete PreShYEn_;
  // delete PreShYCorrEn_;
  delete PreShYEta_;
  delete PreShYPhi_;
  // delete CTPt_;
  // delete CTEta_;
  // delete CTPhi_;

  delete qPFJetPt_;
  delete qPFJetEta_;
  delete qPFJetPhi_;

//-------------------------Sorted qvariables start here.

//PFJet sorted variables
delete qPFJet0Pt_;
delete qPFJet1Pt_; 
delete qPFJet2Pt_; 
delete qPFJet3Pt_; 
delete qPFJet4Pt_; 
delete qPFJet5Pt_; 

delete qPFJet0Eta_;
delete qPFJet1Eta_;
delete qPFJet2Eta_;
delete qPFJet3Eta_;
delete qPFJet4Eta_;
delete qPFJet5Eta_;

delete qPFJet0Phi_;
delete qPFJet1Phi_;
delete qPFJet2Phi_;
delete qPFJet3Phi_;
delete qPFJet4Phi_;
delete qPFJet5Phi_;
//PFJet4CHS sorted variables
delete qPFJet4CHS0Pt_;
delete qPFJet4CHS1Pt_;
delete qPFJet4CHS2Pt_;
delete qPFJet4CHS3Pt_;
delete qPFJet4CHS4Pt_;
delete qPFJet4CHS5Pt_;

delete qPFJet4CHS0Eta_;
delete qPFJet4CHS1Eta_;
delete qPFJet4CHS2Eta_;
delete qPFJet4CHS3Eta_;
delete qPFJet4CHS4Eta_;
delete qPFJet4CHS5Eta_;

delete qPFJet4CHS0Phi_;
delete qPFJet4CHS1Phi_;
delete qPFJet4CHS2Phi_;
delete qPFJet4CHS3Phi_;
delete qPFJet4CHS4Phi_;
delete qPFJet4CHS5Phi_;
 //PF8CHS sorted variables
delete qPFJet8CHS0Pt_;
delete qPFJet8CHS1Pt_;
delete qPFJet8CHS2Pt_;
delete qPFJet8CHS3Pt_;
delete qPFJet8CHS4Pt_;
delete qPFJet8CHS5Pt_;

delete qPFJet8CHS0Eta_;
delete qPFJet8CHS1Eta_;
delete qPFJet8CHS2Eta_;
delete qPFJet8CHS3Eta_;
delete qPFJet8CHS4Eta_;
delete qPFJet8CHS5Eta_;

delete qPFJet8CHS0Phi_;
delete qPFJet8CHS1Phi_;
delete qPFJet8CHS2Phi_;
delete qPFJet8CHS3Phi_;
delete qPFJet8CHS4Phi_;
delete qPFJet8CHS5Phi_;
  //PFJetEI sorted variables
delete qPFJetEI0Pt_;
delete qPFJetEI1Pt_;
delete qPFJetEI2Pt_;
delete qPFJetEI3Pt_;
delete qPFJetEI4Pt_;
delete qPFJetEI5Pt_;

delete qPFJetEI0Eta_;
delete qPFJetEI1Eta_;
delete qPFJetEI2Eta_;
delete qPFJetEI3Eta_;
delete qPFJetEI4Eta_;
delete qPFJetEI5Eta_;

delete qPFJetEI0Phi_;
delete qPFJetEI1Phi_;
delete qPFJetEI2Phi_;
delete qPFJetEI3Phi_;
delete qPFJetEI4Phi_;
delete qPFJetEI5Phi_;

//   //8CHSSoftDrop sorted variables
// delete qPFJet8CHSSD0Pt_;
// delete qPFJet8CHSSD1Pt_;
// delete qPFJet8CHSSD2Pt_;
// delete qPFJet8CHSSD3Pt_;
// delete qPFJet8CHSSD4Pt_;
// delete qPFJet8CHSSD5Pt_;

// delete qPFJet8CHSSD0Eta_;
// delete qPFJet8CHSSD1Eta_;
// delete qPFJet8CHSSD2Eta_;
// delete qPFJet8CHSSD3Eta_;
// delete qPFJet8CHSSD4Eta_;
// delete qPFJet8CHSSD5Eta_;

// delete qPFJet8CHSSD0Phi_;
// delete qPFJet8CHSSD1Phi_;
// delete qPFJet8CHSSD2Phi_;
// delete qPFJet8CHSSD3Phi_;
// delete qPFJet8CHSSD4Phi_;
// delete qPFJet8CHSSD5Phi_;
//   //TopCHS sorted variables
// delete qPFJetTopCHS0Pt_;
// delete qPFJetTopCHS1Pt_;
// delete qPFJetTopCHS2Pt_;
// delete qPFJetTopCHS3Pt_;
// delete qPFJetTopCHS4Pt_;
// delete qPFJetTopCHS5Pt_;

// delete qPFJetTopCHS0Eta_;
// delete qPFJetTopCHS1Eta_;
// delete qPFJetTopCHS2Eta_;
// delete qPFJetTopCHS3Eta_;
// delete qPFJetTopCHS4Eta_;
// delete qPFJetTopCHS5Eta_;

// delete qPFJetTopCHS0Phi_;
// delete qPFJetTopCHS1Phi_;
// delete qPFJetTopCHS2Phi_;
// delete qPFJetTopCHS3Phi_;
// delete qPFJetTopCHS4Phi_;
// delete qPFJetTopCHS5Phi_;


  //CaloJet sorted variables
delete qCalJet0Pt_;
delete qCalJet1Pt_;
delete qCalJet2Pt_;
delete qCalJet3Pt_;
delete qCalJet4Pt_;
delete qCalJet5Pt_;

delete qCalJet0Eta_;
delete qCalJet1Eta_;
delete qCalJet2Eta_;
delete qCalJet3Eta_;
delete qCalJet4Eta_;
delete qCalJet5Eta_;

delete qCalJet0Phi_;
delete qCalJet1Phi_;
delete qCalJet2Phi_;
delete qCalJet3Phi_;
delete qCalJet4Phi_;
delete qCalJet5Phi_;

delete qCalJet0En_;
delete qCalJet1En_;
delete qCalJet2En_;
delete qCalJet3En_;
delete qCalJet4En_;
delete qCalJet5En_;


  //photon sorted variables (not all of them)
delete qPho0Pt_;
delete qPho1Pt_;
delete qPho2Pt_;
delete qPho3Pt_;
delete qPho4Pt_;
delete qPho5Pt_;

delete qPho0Eta_;
delete qPho1Eta_;
delete qPho2Eta_;
delete qPho3Eta_;
delete qPho4Eta_;
delete qPho5Eta_;

delete qPho0Phi_;
delete qPho1Phi_;
delete qPho2Phi_;
delete qPho3Phi_;
delete qPho4Phi_;
delete qPho5Phi_;

delete qPho0En_;
delete qPho1En_;
delete qPho2En_;
delete qPho3En_;
delete qPho4En_;
delete qPho5En_;

  //ged qPhotons sorted variables (not all of them)
delete qgedPho0Pt_;
delete qgedPho1Pt_;
delete qgedPho2Pt_;
delete qgedPho3Pt_;
delete qgedPho4Pt_;
delete qgedPho5Pt_;
 
delete qgedPho0Eta_;
delete qgedPho1Eta_;
delete qgedPho2Eta_;
delete qgedPho3Eta_;
delete qgedPho4Eta_;
delete qgedPho5Eta_;

delete qgedPho0Phi_;
delete qgedPho1Phi_;
delete qgedPho2Phi_;
delete qgedPho3Phi_;
delete qgedPho4Phi_;
delete qgedPho5Phi_;

delete qgedPho0En_;
delete qgedPho1En_;
delete qgedPho2En_;
delete qgedPho3En_;
delete qgedPho4En_;
delete qgedPho5En_;


    //Muon sorted variables (not all of them)
delete qMu0Pt_;
delete qMu1Pt_;
delete qMu2Pt_;
delete qMu3Pt_;
delete qMu4Pt_;
delete qMu5Pt_;

delete qMu0Eta_;
delete qMu1Eta_;
delete qMu2Eta_;
delete qMu3Eta_;
delete qMu4Eta_;
delete qMu5Eta_;

delete qMu0Phi_;
delete qMu1Phi_;
delete qMu2Phi_;
delete qMu3Phi_;
delete qMu4Phi_;
delete qMu5Phi_;

delete qMu0En_;
delete qMu1En_;
delete qMu2En_;
delete qMu3En_;
delete qMu4En_;
delete qMu5En_;
  //Muon sorted cosmic variables (not all of them)
delete qMuCosm0Pt_;
delete qMuCosm1Pt_;
delete qMuCosm2Pt_;
delete qMuCosm3Pt_;
delete qMuCosm4Pt_;
delete qMuCosm5Pt_;

delete qMuCosm0Eta_;
delete qMuCosm1Eta_;
delete qMuCosm2Eta_;
delete qMuCosm3Eta_;
delete qMuCosm4Eta_;
delete qMuCosm5Eta_;

delete qMuCosm0Phi_;
delete qMuCosm1Phi_;
delete qMuCosm2Phi_;
delete qMuCosm3Phi_;
delete qMuCosm4Phi_;
delete qMuCosm5Phi_;

delete qMuCosm0En_;
delete qMuCosm1En_;
delete qMuCosm2En_;
delete qMuCosm3En_;
delete qMuCosm4En_;
delete qMuCosm5En_;


  //Muon Cosmic1Leg sorted variables (not all of them)
delete qMuCosmLeg0Pt_;
delete qMuCosmLeg1Pt_;
delete qMuCosmLeg2Pt_;
delete qMuCosmLeg3Pt_;
delete qMuCosmLeg4Pt_;
delete qMuCosmLeg5Pt_;

delete qMuCosmLeg0Eta_;
delete qMuCosmLeg1Eta_;
delete qMuCosmLeg2Eta_;
delete qMuCosmLeg3Eta_;
delete qMuCosmLeg4Eta_;
delete qMuCosmLeg5Eta_;

delete qMuCosmLeg0Phi_;
delete qMuCosmLeg1Phi_;
delete qMuCosmLeg2Phi_;
delete qMuCosmLeg3Phi_;
delete qMuCosmLeg4Phi_;
delete qMuCosmLeg5Phi_;

delete qMuCosmLeg0En_;
delete qMuCosmLeg1En_;
delete qMuCosmLeg2En_;
delete qMuCosmLeg3En_;
delete qMuCosmLeg4En_;
delete qMuCosmLeg5En_;


 
//_________________________Sorted qvariables end here.  

  delete qPFJet4CHSPt_;
  delete qPFJet4CHSEta_;
  delete qPFJet4CHSPhi_;

  delete qPFJet8CHSPt_;
  delete qPFJet8CHSEta_;
  delete qPFJet8CHSPhi_;

  delete qPFJetEIPt_;
  delete qPFJetEIEta_;
  delete qPFJetEIPhi_;

  delete qPFJet8CHSSDPt_;
  delete qPFJet8CHSSDEta_;
  delete qPFJet8CHSSDPhi_;

  delete qPFJetTopCHSPt_;
  delete qPFJetTopCHSEta_;
  delete qPFJetTopCHSPhi_;

  delete qPFChMetPt_;
  // delete qPFChMetEta_;  
  delete qPFChMetPhi_;
  delete qPFMetPt_;
  // delete qPFMetEta_;
  delete qPFMetPhi_;

  delete qCalJetPt_;
  delete qCalJetEta_;
  delete qCalJetPhi_;
  delete qCalJetEn_;
  delete qCalMETPt_;
  // delete qCalMETEta_;
  delete qCalMETPhi_;
  delete qCalMETEn_;

  delete qCalMETBEPt_;
  // delete qCalMETBEEta_;
  delete qCalMETBEPhi_;
  delete qCalMETBEEn_;
  delete qCalMETBEFOPt_;
  // delete qCalMETBEFOEta_;
  delete qCalMETBEFOPhi_;
  delete qCalMETBEFOEn_;
  delete qCalMETMPt_;
  // delete qCalMETMEta_;
  delete qCalMETMPhi_;
  delete qCalMETMEn_;

  delete qSCEn_;
  delete qSCEta_;
  delete qSCPhi_;
  delete qSCEtaWidth_;
  delete qSCPhiWidth_;  
  delete qSCEnhfEM_;
  delete qSCEtahfEM_;
  delete qSCPhihfEM_;
  // delete qSCEtaWidthhfEM_;
  // delete qSCPhiWidthhfEM_; 
  delete qSCEn5x5_;
  delete qSCEta5x5_;
  delete qSCPhi5x5_;
  delete qSCEtaWidth5x5_;
  delete qSCPhiWidth5x5_; 
  delete qCCEn_;
  delete qCCEta_;
  delete qCCPhi_;
  delete qCCEn5x5_;
  delete qCCEta5x5_;
  delete qCCPhi5x5_;

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
  delete qMuChi2_;  
  delete qMuCosmPt_;
  delete qMuCosmEta_;
  delete qMuCosmPhi_;
  delete qMuCosmEn_;
  delete qMuCosmCh_;
  delete qMuCosmChi2_;    
  delete qMuCosmLegPt_;
  delete qMuCosmLegEta_;
  delete qMuCosmLegPhi_;
  delete qMuCosmLegEn_;
  delete qMuCosmLegCh_;    
  delete qMuCosmLegChi2_;    
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
  delete qEBiEta_;
  delete qEBiPhi_;
  delete qEEenergy_;
  delete qEEtime_;
  delete qEEchi2_;
  delete qEEix_;
  delete qEEiy_;
  delete qESenergy_;
  delete qEStime_;
  // delete qESchi2_;
  delete qESix_;
  delete qESiy_;

  delete qHBHEenergy_ ;
  delete qHBHEtime_ ;
  delete qHBHEauxe_ ;
  delete qHBHEieta_;
  delete qHBHEiphi_;  
  delete qHFenergy_ ;
  delete qHFtime_;
  //delete qHFchi2_  ;
  delete qHFieta_;
  delete qHFiphi_;
  // delete qHOenergy_ ;
  // delete qHOtime_;
  // delete qHOieta_;
  // delete qHOiphi_;
  // //delete qHOchi2_;    

  delete qPreShEn_;
  // delete qPreShCorrEn_;
  delete qPreShEta_;
  delete qPreShPhi_;
  delete qPreShYEn_;
  // delete qPreShYCorrEn_;
  delete qPreShYEta_;
  delete qPreShYPhi_;

  // delete qCTPt_;
  // delete qCTEta_;
  // delete qCTPhi_;

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

//---------------------Sorted variables start here.

//PFJet sorted variables
  computeMeanAndRms(  PFJet0Pt_,  qPFJet0Pt_);
  computeMeanAndRms(  PFJet1Pt_,  qPFJet1Pt_); 
  computeMeanAndRms(  PFJet2Pt_,  qPFJet2Pt_); 
  computeMeanAndRms(  PFJet3Pt_,  qPFJet3Pt_); 
  computeMeanAndRms(  PFJet4Pt_,  qPFJet4Pt_); 
  computeMeanAndRms(  PFJet5Pt_,  qPFJet5Pt_); 

  computeMeanAndRms(  PFJet0Eta_, qPFJet0Eta_);
  computeMeanAndRms(  PFJet1Eta_, qPFJet1Eta_);
  computeMeanAndRms(  PFJet2Eta_, qPFJet2Eta_);
  computeMeanAndRms(  PFJet3Eta_, qPFJet3Eta_);
  computeMeanAndRms(  PFJet4Eta_, qPFJet4Eta_);
  computeMeanAndRms(  PFJet5Eta_ ,qPFJet5Eta_);

  computeMeanAndRms(  PFJet0Phi_, qPFJet0Phi_);
  computeMeanAndRms(  PFJet1Phi_, qPFJet1Phi_);
  computeMeanAndRms(  PFJet2Phi_, qPFJet2Phi_);
  computeMeanAndRms(  PFJet3Phi_, qPFJet3Phi_);
  computeMeanAndRms(  PFJet4Phi_, qPFJet4Phi_);
  computeMeanAndRms(  PFJet5Phi_, qPFJet5Phi_);
//PFJet4CHS sorted variables
  computeMeanAndRms(  PFJet4CHS0Pt_, qPFJet4CHS0Pt_);
  computeMeanAndRms(  PFJet4CHS1Pt_, qPFJet4CHS1Pt_);
  computeMeanAndRms(  PFJet4CHS2Pt_, qPFJet4CHS2Pt_);
  computeMeanAndRms(  PFJet4CHS3Pt_, qPFJet4CHS3Pt_);
  computeMeanAndRms(  PFJet4CHS4Pt_, qPFJet4CHS4Pt_);
  computeMeanAndRms(  PFJet4CHS5Pt_, qPFJet4CHS5Pt_);

  computeMeanAndRms(  PFJet4CHS0Eta_, qPFJet4CHS0Eta_);
  computeMeanAndRms(  PFJet4CHS1Eta_, qPFJet4CHS1Eta_);
  computeMeanAndRms(  PFJet4CHS2Eta_, qPFJet4CHS2Eta_);
  computeMeanAndRms(  PFJet4CHS3Eta_, qPFJet4CHS3Eta_);
  computeMeanAndRms(  PFJet4CHS4Eta_, qPFJet4CHS4Eta_);
  computeMeanAndRms(  PFJet4CHS5Eta_, qPFJet4CHS5Eta_);

  computeMeanAndRms(  PFJet4CHS0Phi_, qPFJet4CHS0Phi_);
  computeMeanAndRms(  PFJet4CHS1Phi_, qPFJet4CHS1Phi_);
  computeMeanAndRms(  PFJet4CHS2Phi_, qPFJet4CHS2Phi_);
  computeMeanAndRms(  PFJet4CHS3Phi_, qPFJet4CHS3Phi_);
  computeMeanAndRms(  PFJet4CHS4Phi_, qPFJet4CHS4Phi_);
  computeMeanAndRms(  PFJet4CHS5Phi_, qPFJet4CHS5Phi_);
 //PF8CHS sorted variables
  computeMeanAndRms(  PFJet8CHS0Pt_, qPFJet8CHS0Pt_);
  computeMeanAndRms(  PFJet8CHS1Pt_, qPFJet8CHS1Pt_);
  computeMeanAndRms(  PFJet8CHS2Pt_, qPFJet8CHS2Pt_);
  computeMeanAndRms(  PFJet8CHS3Pt_, qPFJet8CHS3Pt_);
  computeMeanAndRms(  PFJet8CHS4Pt_, qPFJet8CHS4Pt_);
  computeMeanAndRms(  PFJet8CHS5Pt_, qPFJet8CHS5Pt_);

  computeMeanAndRms(  PFJet8CHS0Eta_, qPFJet8CHS0Eta_);
  computeMeanAndRms(  PFJet8CHS1Eta_, qPFJet8CHS1Eta_);
  computeMeanAndRms(  PFJet8CHS2Eta_, qPFJet8CHS2Eta_);
  computeMeanAndRms(  PFJet8CHS3Eta_, qPFJet8CHS3Eta_);
  computeMeanAndRms(  PFJet8CHS4Eta_, qPFJet8CHS4Eta_);
  computeMeanAndRms(  PFJet8CHS5Eta_, qPFJet8CHS5Eta_);

  computeMeanAndRms(  PFJet8CHS0Phi_, qPFJet8CHS0Phi_);
  computeMeanAndRms(  PFJet8CHS1Phi_, qPFJet8CHS1Phi_);
  computeMeanAndRms(  PFJet8CHS2Phi_, qPFJet8CHS2Phi_);
  computeMeanAndRms(  PFJet8CHS3Phi_, qPFJet8CHS3Phi_);
  computeMeanAndRms(  PFJet8CHS4Phi_, qPFJet8CHS4Phi_);
  computeMeanAndRms(  PFJet8CHS5Phi_, qPFJet8CHS5Phi_);
  //PFJetEI sorted variables
  computeMeanAndRms(  PFJetEI0Pt_, qPFJetEI0Pt_);
  computeMeanAndRms(  PFJetEI1Pt_, qPFJetEI1Pt_);
  computeMeanAndRms(  PFJetEI2Pt_, qPFJetEI2Pt_);
  computeMeanAndRms(  PFJetEI3Pt_, qPFJetEI3Pt_);
  computeMeanAndRms(  PFJetEI4Pt_, qPFJetEI4Pt_);
  computeMeanAndRms(  PFJetEI5Pt_, qPFJetEI5Pt_);

  computeMeanAndRms(  PFJetEI0Eta_, qPFJetEI0Eta_);
  computeMeanAndRms(  PFJetEI1Eta_, qPFJetEI1Eta_);
  computeMeanAndRms(  PFJetEI2Eta_, qPFJetEI2Eta_);
  computeMeanAndRms(  PFJetEI3Eta_, qPFJetEI3Eta_);
  computeMeanAndRms(  PFJetEI4Eta_, qPFJetEI4Eta_);
  computeMeanAndRms(  PFJetEI5Eta_, qPFJetEI5Eta_);

  computeMeanAndRms(  PFJetEI0Phi_, qPFJetEI0Phi_);
  computeMeanAndRms(  PFJetEI1Phi_, qPFJetEI1Phi_);
  computeMeanAndRms(  PFJetEI2Phi_, qPFJetEI2Phi_);
  computeMeanAndRms(  PFJetEI3Phi_, qPFJetEI3Phi_);
  computeMeanAndRms(  PFJetEI4Phi_, qPFJetEI4Phi_);
  computeMeanAndRms(  PFJetEI5Phi_, qPFJetEI5Phi_);

  // //8CHSSoftDrop sorted variables
  // computeMeanAndRms(  PFJet8CHSSD0Pt_, qPFJet8CHSSD0Pt_);
  // computeMeanAndRms(  PFJet8CHSSD1Pt_, qPFJet8CHSSD1Pt_);
  // computeMeanAndRms(  PFJet8CHSSD2Pt_, qPFJet8CHSSD2Pt_);
  // computeMeanAndRms(  PFJet8CHSSD3Pt_, qPFJet8CHSSD3Pt_);
  // computeMeanAndRms(  PFJet8CHSSD4Pt_, qPFJet8CHSSD4Pt_);
  // computeMeanAndRms(  PFJet8CHSSD5Pt_, qPFJet8CHSSD5Pt_);

  // computeMeanAndRms(  PFJet8CHSSD0Eta_, qPFJet8CHSSD0Eta_);
  // computeMeanAndRms(  PFJet8CHSSD1Eta_, qPFJet8CHSSD1Eta_);
  // computeMeanAndRms(  PFJet8CHSSD2Eta_, qPFJet8CHSSD2Eta_);
  // computeMeanAndRms(  PFJet8CHSSD3Eta_, qPFJet8CHSSD3Eta_);
  // computeMeanAndRms(  PFJet8CHSSD4Eta_, qPFJet8CHSSD4Eta_);
  // computeMeanAndRms(  PFJet8CHSSD5Eta_, qPFJet8CHSSD5Eta_);

  // computeMeanAndRms(  PFJet8CHSSD0Phi_, qPFJet8CHSSD0Phi_);
  // computeMeanAndRms(  PFJet8CHSSD1Phi_, qPFJet8CHSSD1Phi_);
  // computeMeanAndRms(  PFJet8CHSSD2Phi_, qPFJet8CHSSD2Phi_);
  // computeMeanAndRms(  PFJet8CHSSD3Phi_, qPFJet8CHSSD3Phi_);
  // computeMeanAndRms(  PFJet8CHSSD4Phi_, qPFJet8CHSSD4Phi_);
  // computeMeanAndRms(  PFJet8CHSSD5Phi_, qPFJet8CHSSD5Phi_);
  // //TopCHS sorted variables
  // computeMeanAndRms(  PFJetTopCHS0Pt_, qPFJetTopCHS0Pt_);
  // computeMeanAndRms(  PFJetTopCHS1Pt_, qPFJetTopCHS1Pt_);
  // computeMeanAndRms(  PFJetTopCHS2Pt_, qPFJetTopCHS2Pt_);
  // computeMeanAndRms(  PFJetTopCHS3Pt_, qPFJetTopCHS3Pt_);
  // computeMeanAndRms(  PFJetTopCHS4Pt_, qPFJetTopCHS4Pt_);
  // computeMeanAndRms(  PFJetTopCHS5Pt_, qPFJetTopCHS5Pt_);

  // computeMeanAndRms(  PFJetTopCHS0Eta_,  qPFJetTopCHS0Eta_);
  // computeMeanAndRms(  PFJetTopCHS1Eta_,  qPFJetTopCHS1Eta_);
  // computeMeanAndRms(  PFJetTopCHS2Eta_,  qPFJetTopCHS2Eta_);
  // computeMeanAndRms(  PFJetTopCHS3Eta_,  qPFJetTopCHS3Eta_);
  // computeMeanAndRms(  PFJetTopCHS4Eta_,  qPFJetTopCHS4Eta_);
  // computeMeanAndRms(  PFJetTopCHS5Eta_,  qPFJetTopCHS5Eta_);

  // computeMeanAndRms(  PFJetTopCHS0Phi_, qPFJetTopCHS0Phi_);
  // computeMeanAndRms(  PFJetTopCHS1Phi_, qPFJetTopCHS1Phi_);
  // computeMeanAndRms(  PFJetTopCHS2Phi_, qPFJetTopCHS2Phi_);
  // computeMeanAndRms(  PFJetTopCHS3Phi_, qPFJetTopCHS3Phi_);
  // computeMeanAndRms(  PFJetTopCHS4Phi_, qPFJetTopCHS4Phi_);
  // computeMeanAndRms(  PFJetTopCHS5Phi_, qPFJetTopCHS5Phi_);


  //CaloJet sorted variables
  computeMeanAndRms(  CalJet0Pt_, qCalJet0Pt_);
  computeMeanAndRms(  CalJet1Pt_, qCalJet1Pt_);
  computeMeanAndRms(  CalJet2Pt_, qCalJet2Pt_);
  computeMeanAndRms(  CalJet3Pt_, qCalJet3Pt_);
  computeMeanAndRms(  CalJet4Pt_, qCalJet4Pt_);
  computeMeanAndRms(  CalJet5Pt_, qCalJet5Pt_);

  computeMeanAndRms(  CalJet0Eta_, qCalJet0Eta_);
  computeMeanAndRms(  CalJet1Eta_, qCalJet1Eta_);
  computeMeanAndRms(  CalJet2Eta_, qCalJet2Eta_);
  computeMeanAndRms(  CalJet3Eta_, qCalJet3Eta_);
  computeMeanAndRms(  CalJet4Eta_, qCalJet4Eta_);
  computeMeanAndRms(  CalJet5Eta_, qCalJet5Eta_);

  computeMeanAndRms(  CalJet0Phi_, qCalJet0Phi_);
  computeMeanAndRms(  CalJet1Phi_, qCalJet1Phi_);
  computeMeanAndRms(  CalJet2Phi_, qCalJet2Phi_);
  computeMeanAndRms(  CalJet3Phi_, qCalJet3Phi_);
  computeMeanAndRms(  CalJet4Phi_, qCalJet4Phi_);
  computeMeanAndRms(  CalJet5Phi_, qCalJet5Phi_);

  computeMeanAndRms(  CalJet0En_, qCalJet0En_);
  computeMeanAndRms(  CalJet1En_, qCalJet1En_);
  computeMeanAndRms(  CalJet2En_, qCalJet2En_);
  computeMeanAndRms(  CalJet3En_, qCalJet3En_);
  computeMeanAndRms(  CalJet4En_, qCalJet4En_);
  computeMeanAndRms(  CalJet5En_, qCalJet5En_);


  //photon sorted variables (not all of them)
  computeMeanAndRms(  Pho0Pt_, qPho0Pt_);
  computeMeanAndRms(  Pho1Pt_, qPho1Pt_);
  computeMeanAndRms(  Pho2Pt_, qPho2Pt_);
  computeMeanAndRms(  Pho3Pt_, qPho3Pt_);
  computeMeanAndRms(  Pho4Pt_, qPho4Pt_);
  computeMeanAndRms(  Pho5Pt_, qPho5Pt_);

  computeMeanAndRms(  Pho0Eta_, qPho0Eta_);
  computeMeanAndRms(  Pho1Eta_, qPho1Eta_);
  computeMeanAndRms(  Pho2Eta_, qPho2Eta_);
  computeMeanAndRms(  Pho3Eta_, qPho3Eta_);
  computeMeanAndRms(  Pho4Eta_, qPho4Eta_);
  computeMeanAndRms(  Pho5Eta_, qPho5Eta_);

  computeMeanAndRms(  Pho0Phi_, qPho0Phi_);
  computeMeanAndRms(  Pho1Phi_, qPho1Phi_);
  computeMeanAndRms(  Pho2Phi_, qPho2Phi_);
  computeMeanAndRms(  Pho3Phi_, qPho3Phi_);
  computeMeanAndRms(  Pho4Phi_, qPho4Phi_);
  computeMeanAndRms(  Pho5Phi_, qPho5Phi_);

  computeMeanAndRms(  Pho0En_, qPho0En_);
  computeMeanAndRms(  Pho1En_, qPho1En_);
  computeMeanAndRms(  Pho2En_, qPho2En_);
  computeMeanAndRms(  Pho3En_, qPho3En_);
  computeMeanAndRms(  Pho4En_, qPho4En_);
  computeMeanAndRms(  Pho5En_, qPho5En_);

  //ged computeMeanAndRms(qPhotons sorted variables (not all of them)
  computeMeanAndRms(  gedPho0Pt_, qgedPho0Pt_);
  computeMeanAndRms(  gedPho1Pt_, qgedPho1Pt_);
  computeMeanAndRms(  gedPho2Pt_, qgedPho2Pt_);
  computeMeanAndRms(  gedPho3Pt_, qgedPho3Pt_);
  computeMeanAndRms(  gedPho4Pt_, qgedPho4Pt_);
  computeMeanAndRms(  gedPho5Pt_, qgedPho5Pt_);
 
  computeMeanAndRms(  gedPho0Eta_, qgedPho0Eta_);
  computeMeanAndRms(  gedPho1Eta_, qgedPho1Eta_);
  computeMeanAndRms(  gedPho2Eta_, qgedPho2Eta_);
  computeMeanAndRms(  gedPho3Eta_, qgedPho3Eta_);
  computeMeanAndRms(  gedPho4Eta_, qgedPho4Eta_);
  computeMeanAndRms(  gedPho5Eta_, qgedPho5Eta_);

  computeMeanAndRms(  gedPho0Phi_, qgedPho0Phi_);
  computeMeanAndRms(  gedPho1Phi_, qgedPho1Phi_);
  computeMeanAndRms(  gedPho2Phi_, qgedPho2Phi_);
  computeMeanAndRms(  gedPho3Phi_, qgedPho3Phi_);
  computeMeanAndRms(  gedPho4Phi_, qgedPho4Phi_);
  computeMeanAndRms(  gedPho5Phi_, qgedPho5Phi_);

  computeMeanAndRms(  gedPho0En_, qgedPho0En_);
  computeMeanAndRms(  gedPho1En_, qgedPho1En_);
  computeMeanAndRms(  gedPho2En_, qgedPho2En_);
  computeMeanAndRms(  gedPho3En_, qgedPho3En_);
  computeMeanAndRms(  gedPho4En_, qgedPho4En_);
  computeMeanAndRms(  gedPho5En_, qgedPho5En_);


    //Muon sorted variables (not all of them)
  computeMeanAndRms(  Mu0Pt_, qMu0Pt_);
  computeMeanAndRms(  Mu1Pt_, qMu1Pt_);
  computeMeanAndRms(  Mu2Pt_, qMu2Pt_);
  computeMeanAndRms(  Mu3Pt_, qMu3Pt_);
  computeMeanAndRms(  Mu4Pt_, qMu4Pt_);
  computeMeanAndRms(  Mu5Pt_, qMu5Pt_);

  computeMeanAndRms(  Mu0Eta_, qMu0Eta_);
  computeMeanAndRms(  Mu1Eta_, qMu1Eta_);
  computeMeanAndRms(  Mu2Eta_, qMu2Eta_);
  computeMeanAndRms(  Mu3Eta_, qMu3Eta_);
  computeMeanAndRms(  Mu4Eta_, qMu4Eta_);
  computeMeanAndRms(  Mu5Eta_, qMu5Eta_);

  computeMeanAndRms(  Mu0Phi_, qMu0Phi_);
  computeMeanAndRms(  Mu1Phi_, qMu1Phi_);
  computeMeanAndRms(  Mu2Phi_, qMu2Phi_);
  computeMeanAndRms(  Mu3Phi_, qMu3Phi_);
  computeMeanAndRms(  Mu4Phi_, qMu4Phi_);
  computeMeanAndRms(  Mu5Phi_, qMu5Phi_);

  computeMeanAndRms(  Mu0En_, qMu0En_);
  computeMeanAndRms(  Mu1En_, qMu1En_);
  computeMeanAndRms(  Mu2En_, qMu2En_);
  computeMeanAndRms(  Mu3En_, qMu3En_);
  computeMeanAndRms(  Mu4En_, qMu4En_);
  computeMeanAndRms(  Mu5En_, qMu5En_);
  //Muon sorted cosmic variables (not all of them)
  computeMeanAndRms(  MuCosm0Pt_, qMuCosm0Pt_);
  computeMeanAndRms(  MuCosm1Pt_, qMuCosm1Pt_);
  computeMeanAndRms(  MuCosm2Pt_, qMuCosm2Pt_);
  computeMeanAndRms(  MuCosm3Pt_, qMuCosm3Pt_);
  computeMeanAndRms(  MuCosm4Pt_, qMuCosm4Pt_);
  computeMeanAndRms(  MuCosm5Pt_, qMuCosm5Pt_);

  computeMeanAndRms(  MuCosm0Eta_, qMuCosm0Eta_);
  computeMeanAndRms(  MuCosm1Eta_, qMuCosm1Eta_);
  computeMeanAndRms(  MuCosm2Eta_, qMuCosm2Eta_);
  computeMeanAndRms(  MuCosm3Eta_, qMuCosm3Eta_);
  computeMeanAndRms(  MuCosm4Eta_, qMuCosm4Eta_);
  computeMeanAndRms(  MuCosm5Eta_, qMuCosm5Eta_);

  computeMeanAndRms(  MuCosm0Phi_, qMuCosm0Phi_);
  computeMeanAndRms(  MuCosm1Phi_, qMuCosm1Phi_);
  computeMeanAndRms(  MuCosm2Phi_, qMuCosm2Phi_);
  computeMeanAndRms(  MuCosm3Phi_, qMuCosm3Phi_);
  computeMeanAndRms(  MuCosm4Phi_, qMuCosm4Phi_);
  computeMeanAndRms(  MuCosm5Phi_, qMuCosm5Phi_);

  computeMeanAndRms(  MuCosm0En_, qMuCosm0En_);
  computeMeanAndRms(  MuCosm1En_, qMuCosm1En_);
  computeMeanAndRms(  MuCosm2En_, qMuCosm2En_);
  computeMeanAndRms(  MuCosm3En_, qMuCosm3En_);
  computeMeanAndRms(  MuCosm4En_, qMuCosm4En_);
  computeMeanAndRms(  MuCosm5En_, qMuCosm5En_);


  //Muon Cosmic1Leg sorted variables (not all of them)
  computeMeanAndRms(  MuCosmLeg0Pt_, qMuCosmLeg0Pt_);
  computeMeanAndRms(  MuCosmLeg1Pt_, qMuCosmLeg1Pt_);
  computeMeanAndRms(  MuCosmLeg2Pt_, qMuCosmLeg2Pt_);
  computeMeanAndRms(  MuCosmLeg3Pt_, qMuCosmLeg3Pt_);
  computeMeanAndRms(  MuCosmLeg4Pt_, qMuCosmLeg4Pt_);
  computeMeanAndRms(  MuCosmLeg5Pt_, qMuCosmLeg5Pt_);

  computeMeanAndRms(  MuCosmLeg0Eta_, qMuCosmLeg0Eta_);
  computeMeanAndRms(  MuCosmLeg1Eta_, qMuCosmLeg1Eta_);
  computeMeanAndRms(  MuCosmLeg2Eta_, qMuCosmLeg2Eta_);
  computeMeanAndRms(  MuCosmLeg3Eta_, qMuCosmLeg3Eta_);
  computeMeanAndRms(  MuCosmLeg4Eta_, qMuCosmLeg4Eta_);
  computeMeanAndRms(  MuCosmLeg5Eta_, qMuCosmLeg5Eta_);

  computeMeanAndRms(  MuCosmLeg0Phi_, qMuCosmLeg0Phi_);
  computeMeanAndRms(  MuCosmLeg1Phi_, qMuCosmLeg1Phi_);
  computeMeanAndRms(  MuCosmLeg2Phi_, qMuCosmLeg2Phi_);
  computeMeanAndRms(  MuCosmLeg3Phi_, qMuCosmLeg3Phi_);
  computeMeanAndRms(  MuCosmLeg4Phi_, qMuCosmLeg4Phi_);
  computeMeanAndRms(  MuCosmLeg5Phi_, qMuCosmLeg5Phi_);

  computeMeanAndRms(  MuCosmLeg0En_, qMuCosmLeg0En_);
  computeMeanAndRms(  MuCosmLeg1En_, qMuCosmLeg1En_);
  computeMeanAndRms(  MuCosmLeg2En_, qMuCosmLeg2En_);
  computeMeanAndRms(  MuCosmLeg3En_, qMuCosmLeg3En_);
  computeMeanAndRms(  MuCosmLeg4En_, qMuCosmLeg4En_);
  computeMeanAndRms(  MuCosmLeg5En_, qMuCosmLeg5En_);


//_____________________Sorted variables end here.  

  computeMeanAndRms(PFJet4CHSPt_, qPFJet4CHSPt_);
  computeMeanAndRms(PFJet4CHSEta_,qPFJet4CHSEta_);
  computeMeanAndRms(PFJet4CHSPhi_,qPFJet4CHSPhi_);

  computeMeanAndRms(PFJet8CHSPt_, qPFJet8CHSPt_);
  computeMeanAndRms(PFJet8CHSEta_,qPFJet8CHSEta_);
  computeMeanAndRms(PFJet8CHSPhi_,qPFJet8CHSPhi_);

  computeMeanAndRms(PFJetEIPt_, qPFJetEIPt_);
  computeMeanAndRms(PFJetEIEta_,qPFJetEIEta_);
  computeMeanAndRms(PFJetEIPhi_,qPFJetEIPhi_);

  computeMeanAndRms(PFJet8CHSSDPt_, qPFJet8CHSSDPt_);
  computeMeanAndRms(PFJet8CHSSDEta_,qPFJet8CHSSDEta_);
  computeMeanAndRms(PFJet8CHSSDPhi_,qPFJet8CHSSDPhi_);

  computeMeanAndRms(PFJetTopCHSPt_, qPFJetTopCHSPt_);
  computeMeanAndRms(PFJetTopCHSEta_,qPFJetTopCHSEta_);
  computeMeanAndRms(PFJetTopCHSPhi_,qPFJetTopCHSPhi_);

  computeMeanAndRms(PFChMetPt_, qPFChMetPt_);
  // computeMeanAndRms(PFChMetEta_,  qPFChMetEta_);  
  computeMeanAndRms(PFChMetPhi_,  qPFChMetPhi_);
  computeMeanAndRms(PFMetPt_, qPFMetPt_);
  // computeMeanAndRms(PFMetEta_,  qPFMetEta_);
  computeMeanAndRms(PFMetPhi_,  qPFMetPhi_);
  computeMeanAndRms(nVtx_,    qNVtx_);

  computeMeanAndRms(CalJetPt_, qCalJetPt_);
  computeMeanAndRms(CalJetEta_,qCalJetEta_);
  computeMeanAndRms(CalJetPhi_,qCalJetPhi_);
  computeMeanAndRms(CalJetEn_,qCalJetEn_);

  computeMeanAndRms(CalMETPt_, qCalMETPt_);
  // computeMeanAndRms(CalMETEta_,qCalMETEta_);
  computeMeanAndRms(CalMETPhi_,qCalMETPhi_);
  computeMeanAndRms(CalMETEn_,qCalMETEn_);

  computeMeanAndRms(CalMETBEPt_, qCalMETBEPt_);
  // computeMeanAndRms(CalMETBEEta_,qCalMETBEEta_);
  computeMeanAndRms(CalMETBEPhi_,qCalMETBEPhi_);
  computeMeanAndRms(CalMETBEEn_, qCalMETBEEn_);

  computeMeanAndRms(CalMETBEFOPt_, qCalMETBEFOPt_);
  // computeMeanAndRms(CalMETBEFOEta_,qCalMETBEFOEta_);
  computeMeanAndRms(CalMETBEFOPhi_,qCalMETBEFOPhi_);
  computeMeanAndRms(CalMETBEFOEn_, qCalMETBEFOEn_);

  computeMeanAndRms(CalMETMPt_, qCalMETMPt_);
  // computeMeanAndRms(CalMETMEta_,qCalMETMEta_);
  computeMeanAndRms(CalMETMPhi_,qCalMETMPhi_);
  computeMeanAndRms(CalMETMEn_, qCalMETMEn_);


  computeMeanAndRms(SCEn_, qSCEn_);   
  computeMeanAndRms(SCEta_, qSCEta_);  
  computeMeanAndRms(SCPhi_, qSCPhi_);
  computeMeanAndRms(SCEtaWidth_, qSCEtaWidth_); 
  computeMeanAndRms(SCPhiWidth_, qSCPhiWidth_);  

  computeMeanAndRms(SCEnhfEM_, qSCEnhfEM_);   
  computeMeanAndRms(SCEtahfEM_, qSCEtahfEM_);  
  computeMeanAndRms(SCPhihfEM_, qSCPhihfEM_); 
  // computeMeanAndRms(SCEtaWidthhfEM_, qSCEtaWidthhfEM_);  
  // computeMeanAndRms(SCPhiWidthhfEM_, qSCPhiWidthhfEM_); 

  computeMeanAndRms(SCEn5x5_, qSCEn5x5_);   
  computeMeanAndRms(SCEta5x5_, qSCEta5x5_);  
  computeMeanAndRms(SCPhi5x5_, qSCPhi5x5_);
  computeMeanAndRms(SCEtaWidth5x5_, qSCEtaWidth5x5_);  
  computeMeanAndRms(SCPhiWidth5x5_, qSCPhiWidth5x5_);  
  computeMeanAndRms(CCEn_, qCCEn_);   
  computeMeanAndRms(CCEta_, qCCEta_);  
  computeMeanAndRms(CCPhi_, qCCPhi_);
  computeMeanAndRms(CCEn5x5_, qCCEn5x5_);   
  computeMeanAndRms(CCEta5x5_, qCCEta5x5_);  
  computeMeanAndRms(CCPhi5x5_, qCCPhi5x5_); 


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
  computeMeanAndRms(MuChi2_,qMuChi2_);  
  computeMeanAndRms(MuCosmPt_, qMuCosmPt_);
  computeMeanAndRms(MuCosmEta_,qMuCosmEta_);
  computeMeanAndRms(MuCosmPhi_,qMuCosmPhi_);
  computeMeanAndRms(MuCosmEn_, qMuCosmEn_);
  computeMeanAndRms(MuCosmCh_, qMuCosmCh_);
  computeMeanAndRms(MuCosmChi2_,qMuCosmChi2_);   
  computeMeanAndRms(MuCosmLegPt_, qMuCosmLegPt_);
  computeMeanAndRms(MuCosmLegEta_,qMuCosmLegEta_);
  computeMeanAndRms(MuCosmLegPhi_,qMuCosmLegPhi_);
  computeMeanAndRms(MuCosmLegEn_, qMuCosmLegEn_);
  computeMeanAndRms(MuCosmLegCh_, qMuCosmLegCh_); 
  computeMeanAndRms(MuCosmLegChi2_,qMuCosmLegChi2_);      
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
  computeMeanAndRms(EBiEta_, qEBiEta_);
  computeMeanAndRms(EBiPhi_, qEBiPhi_);

  computeMeanAndRms(EEenergy_, qEEenergy_);
  computeMeanAndRms(EEtime_, qEEtime_);
  computeMeanAndRms(EEchi2_, qEEchi2_);
  computeMeanAndRms(EEix_, qEEix_);
  computeMeanAndRms(EEiy_, qEEiy_);

  computeMeanAndRms(ESenergy_, qESenergy_);
  computeMeanAndRms(EStime_, qEStime_);
  // computeMeanAndRms(ESchi2_, qESchi2_);
  computeMeanAndRms(ESix_, qESix_);
  computeMeanAndRms(ESiy_, qESiy_);

  computeMeanAndRms(HBHEenergy_, qHBHEenergy_);
  computeMeanAndRms(HBHEtime_, qHBHEtime_);
  computeMeanAndRms(HBHEauxe_, qHBHEauxe_);
  computeMeanAndRms(HBHEieta_, qHBHEieta_);
  computeMeanAndRms(HBHEiphi_, qHBHEiphi_);

  computeMeanAndRms(HFenergy_, qHFenergy_);
  computeMeanAndRms(HFtime_, qHFtime_);
  // computeMeanAndRms(HFchi2_, qHFchi2_);
  computeMeanAndRms(HFieta_, qHFieta_);
  computeMeanAndRms(HFiphi_, qHFiphi_);

  // computeMeanAndRms(HOenergy_, qHOenergy_);
  // computeMeanAndRms(HOtime_, qHOtime_);
  // // computeMeanAndRms(HOchi2_, qHOchi2_);
  // computeMeanAndRms(HOieta_, qHOieta_);
  // computeMeanAndRms(HOiphi_, qHOiphi_);


  computeMeanAndRms(PreShEn_, qPreShEn_);   
  // computeMeanAndRms(PreShCorrEn_, qPreShCorrEn_);   
  computeMeanAndRms(PreShEta_, qPreShEta_);  //
  computeMeanAndRms(PreShPhi_, qPreShPhi_);  //
  computeMeanAndRms(PreShYEn_, qPreShYEn_);   
  // computeMeanAndRms(PreShYCorrEn_, qPreShYCorrEn_);   
  computeMeanAndRms(PreShYEta_, qPreShYEta_);  //
  computeMeanAndRms(PreShYPhi_, qPreShYPhi_);  //

  // computeMeanAndRms(CTPt_, qCTPt_);   
  // computeMeanAndRms(CTEta_, qCTEta_);  
  // computeMeanAndRms(CTPhi_, qCTPhi_);


  computeQuantiles(PFJetPt_, qPFJetPt_, quantiles_);
  computeQuantiles(PFJetEta_,qPFJetEta_,quantiles_);
  computeQuantiles(PFJetPhi_,qPFJetPhi_,quantiles_);

//------------------------------------------Sorted variables start here.


//PFJet sorted variables
  computeQuantiles(  PFJet0Pt_ , qPFJet0Pt_, quantiles_);
  computeQuantiles(  PFJet1Pt_ , qPFJet1Pt_, quantiles_); 
  computeQuantiles(  PFJet2Pt_ , qPFJet2Pt_, quantiles_); 
  computeQuantiles(  PFJet3Pt_ , qPFJet3Pt_, quantiles_); 
  computeQuantiles(  PFJet4Pt_ , qPFJet4Pt_, quantiles_); 
  computeQuantiles(  PFJet5Pt_ , qPFJet5Pt_, quantiles_); 

  computeQuantiles(  PFJet0Eta_, qPFJet0Eta_, quantiles_);
  computeQuantiles(  PFJet1Eta_, qPFJet1Eta_, quantiles_);
  computeQuantiles(  PFJet2Eta_, qPFJet2Eta_, quantiles_);
  computeQuantiles(  PFJet3Eta_, qPFJet3Eta_, quantiles_);
  computeQuantiles(  PFJet4Eta_, qPFJet4Eta_, quantiles_);
  computeQuantiles(  PFJet5Eta_, qPFJet5Eta_, quantiles_);

  computeQuantiles(  PFJet0Phi_, qPFJet0Phi_, quantiles_);
  computeQuantiles(  PFJet1Phi_, qPFJet1Phi_, quantiles_);
  computeQuantiles(  PFJet2Phi_, qPFJet2Phi_, quantiles_);
  computeQuantiles(  PFJet3Phi_, qPFJet3Phi_, quantiles_);
  computeQuantiles(  PFJet4Phi_, qPFJet4Phi_, quantiles_);
  computeQuantiles(  PFJet5Phi_, qPFJet5Phi_, quantiles_);
//PFJet4CHS sorted variables
  computeQuantiles(  PFJet4CHS0Pt_, qPFJet4CHS0Pt_, quantiles_);
  computeQuantiles(  PFJet4CHS1Pt_, qPFJet4CHS1Pt_, quantiles_);
  computeQuantiles(  PFJet4CHS2Pt_, qPFJet4CHS2Pt_, quantiles_);
  computeQuantiles(  PFJet4CHS3Pt_, qPFJet4CHS3Pt_, quantiles_);
  computeQuantiles(  PFJet4CHS4Pt_, qPFJet4CHS4Pt_, quantiles_);
  computeQuantiles(  PFJet4CHS5Pt_, qPFJet4CHS5Pt_, quantiles_);

  computeQuantiles(  PFJet4CHS0Eta_, qPFJet4CHS0Eta_, quantiles_);
  computeQuantiles(  PFJet4CHS1Eta_, qPFJet4CHS1Eta_, quantiles_);
  computeQuantiles(  PFJet4CHS2Eta_, qPFJet4CHS2Eta_, quantiles_);
  computeQuantiles(  PFJet4CHS3Eta_, qPFJet4CHS3Eta_, quantiles_);
  computeQuantiles(  PFJet4CHS4Eta_, qPFJet4CHS4Eta_, quantiles_);
  computeQuantiles(  PFJet4CHS5Eta_, qPFJet4CHS5Eta_, quantiles_);

  computeQuantiles(  PFJet4CHS0Phi_, qPFJet4CHS0Phi_, quantiles_);
  computeQuantiles(  PFJet4CHS1Phi_, qPFJet4CHS1Phi_, quantiles_);
  computeQuantiles(  PFJet4CHS2Phi_, qPFJet4CHS2Phi_, quantiles_);
  computeQuantiles(  PFJet4CHS3Phi_, qPFJet4CHS3Phi_, quantiles_);
  computeQuantiles(  PFJet4CHS4Phi_, qPFJet4CHS4Phi_, quantiles_);
  computeQuantiles(  PFJet4CHS5Phi_, qPFJet4CHS5Phi_, quantiles_);

 //PF8CHS sorted variables
 
  computeQuantiles(  PFJet8CHS0Pt_, qPFJet8CHS0Pt_, quantiles_);
  computeQuantiles(  PFJet8CHS1Pt_, qPFJet8CHS1Pt_, quantiles_);
  computeQuantiles(  PFJet8CHS2Pt_, qPFJet8CHS2Pt_, quantiles_);
  computeQuantiles(  PFJet8CHS3Pt_, qPFJet8CHS3Pt_, quantiles_);
  computeQuantiles(  PFJet8CHS4Pt_, qPFJet8CHS4Pt_, quantiles_);
  computeQuantiles(  PFJet8CHS5Pt_, qPFJet8CHS5Pt_, quantiles_);

  computeQuantiles(  PFJet8CHS0Eta_, qPFJet8CHS0Eta_, quantiles_);
  computeQuantiles(  PFJet8CHS1Eta_, qPFJet8CHS1Eta_, quantiles_);
  computeQuantiles(  PFJet8CHS2Eta_, qPFJet8CHS2Eta_, quantiles_);
  computeQuantiles(  PFJet8CHS3Eta_, qPFJet8CHS3Eta_, quantiles_);
  computeQuantiles(  PFJet8CHS4Eta_, qPFJet8CHS4Eta_, quantiles_);
  computeQuantiles(  PFJet8CHS5Eta_, qPFJet8CHS5Eta_, quantiles_);

  computeQuantiles(  PFJet8CHS0Phi_, qPFJet8CHS0Phi_, quantiles_);
  computeQuantiles(  PFJet8CHS1Phi_, qPFJet8CHS1Phi_, quantiles_);
  computeQuantiles(  PFJet8CHS2Phi_, qPFJet8CHS2Phi_, quantiles_);
  computeQuantiles(  PFJet8CHS3Phi_, qPFJet8CHS3Phi_, quantiles_);
  computeQuantiles(  PFJet8CHS4Phi_, qPFJet8CHS4Phi_, quantiles_);
  computeQuantiles(  PFJet8CHS5Phi_, qPFJet8CHS5Phi_, quantiles_);
  //PFJetEI sorted variables
  computeQuantiles(  PFJetEI0Pt_, qPFJetEI0Pt_, quantiles_);
  computeQuantiles(  PFJetEI1Pt_, qPFJetEI1Pt_, quantiles_);
  computeQuantiles(  PFJetEI2Pt_, qPFJetEI2Pt_, quantiles_);
  computeQuantiles(  PFJetEI3Pt_, qPFJetEI3Pt_, quantiles_);
  computeQuantiles(  PFJetEI4Pt_, qPFJetEI4Pt_, quantiles_);
  computeQuantiles(  PFJetEI5Pt_, qPFJetEI5Pt_, quantiles_);

  computeQuantiles(  PFJetEI0Eta_, qPFJetEI0Eta_, quantiles_);
  computeQuantiles(  PFJetEI1Eta_, qPFJetEI1Eta_, quantiles_);
  computeQuantiles(  PFJetEI2Eta_, qPFJetEI2Eta_, quantiles_);
  computeQuantiles(  PFJetEI3Eta_, qPFJetEI3Eta_, quantiles_);
  computeQuantiles(  PFJetEI4Eta_, qPFJetEI4Eta_, quantiles_);
  computeQuantiles(  PFJetEI5Eta_, qPFJetEI5Eta_, quantiles_);

  computeQuantiles(  PFJetEI0Phi_, qPFJetEI0Phi_, quantiles_);
  computeQuantiles(  PFJetEI1Phi_, qPFJetEI1Phi_, quantiles_);
  computeQuantiles(  PFJetEI2Phi_, qPFJetEI2Phi_, quantiles_);
  computeQuantiles(  PFJetEI3Phi_, qPFJetEI3Phi_, quantiles_);
  computeQuantiles(  PFJetEI4Phi_, qPFJetEI4Phi_, quantiles_);
  computeQuantiles(  PFJetEI5Phi_, qPFJetEI5Phi_, quantiles_);

  // //8CHSSoftDrop sorted variables
  // computeQuantiles(  PFJet8CHSSD0Pt_, qPFJet8CHSSD0Pt_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD1Pt_, qPFJet8CHSSD1Pt_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD2Pt_, qPFJet8CHSSD2Pt_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD3Pt_, qPFJet8CHSSD3Pt_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD4Pt_, qPFJet8CHSSD4Pt_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD5Pt_, qPFJet8CHSSD5Pt_, quantiles_);

  // computeQuantiles(  PFJet8CHSSD0Eta_, qPFJet8CHSSD0Eta_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD1Eta_, qPFJet8CHSSD1Eta_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD2Eta_, qPFJet8CHSSD2Eta_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD3Eta_, qPFJet8CHSSD3Eta_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD4Eta_, qPFJet8CHSSD4Eta_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD5Eta_, qPFJet8CHSSD5Eta_, quantiles_);

  // computeQuantiles(  PFJet8CHSSD0Phi_, qPFJet8CHSSD0Phi_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD1Phi_, qPFJet8CHSSD1Phi_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD2Phi_, qPFJet8CHSSD2Phi_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD3Phi_, qPFJet8CHSSD3Phi_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD4Phi_, qPFJet8CHSSD4Phi_, quantiles_);
  // computeQuantiles(  PFJet8CHSSD5Phi_, qPFJet8CHSSD5Phi_, quantiles_);
  // //TopCHS sorted variables
  // computeQuantiles(  PFJetTopCHS0Pt_, qPFJetTopCHS0Pt_, quantiles_);
  // computeQuantiles(  PFJetTopCHS1Pt_, qPFJetTopCHS1Pt_, quantiles_);
  // computeQuantiles(  PFJetTopCHS2Pt_, qPFJetTopCHS2Pt_, quantiles_);
  // computeQuantiles(  PFJetTopCHS3Pt_, qPFJetTopCHS3Pt_, quantiles_);
  // computeQuantiles(  PFJetTopCHS4Pt_, qPFJetTopCHS4Pt_, quantiles_);
  // computeQuantiles(  PFJetTopCHS5Pt_, qPFJetTopCHS5Pt_, quantiles_);

  // computeQuantiles(  PFJetTopCHS0Eta_, qPFJetTopCHS0Eta_, quantiles_);
  // computeQuantiles(  PFJetTopCHS1Eta_, qPFJetTopCHS1Eta_, quantiles_);
  // computeQuantiles(  PFJetTopCHS2Eta_, qPFJetTopCHS2Eta_, quantiles_);
  // computeQuantiles(  PFJetTopCHS3Eta_, qPFJetTopCHS3Eta_, quantiles_);
  // computeQuantiles(  PFJetTopCHS4Eta_, qPFJetTopCHS4Eta_, quantiles_);
  // computeQuantiles(  PFJetTopCHS5Eta_, qPFJetTopCHS5Eta_, quantiles_);

  // computeQuantiles(  PFJetTopCHS0Phi_, qPFJetTopCHS0Phi_, quantiles_);
  // computeQuantiles(  PFJetTopCHS1Phi_, qPFJetTopCHS1Phi_, quantiles_);
  // computeQuantiles(  PFJetTopCHS2Phi_, qPFJetTopCHS2Phi_, quantiles_);
  // computeQuantiles(  PFJetTopCHS3Phi_, qPFJetTopCHS3Phi_, quantiles_);
  // computeQuantiles(  PFJetTopCHS4Phi_, qPFJetTopCHS4Phi_, quantiles_);
  // computeQuantiles(  PFJetTopCHS5Phi_, qPFJetTopCHS5Phi_, quantiles_);


  //CaloJet sorted variables
  computeQuantiles(  CalJet0Pt_, qCalJet0Pt_, quantiles_);
  computeQuantiles(  CalJet1Pt_, qCalJet1Pt_, quantiles_);
  computeQuantiles(  CalJet2Pt_, qCalJet2Pt_, quantiles_);
  computeQuantiles(  CalJet3Pt_, qCalJet3Pt_, quantiles_);
  computeQuantiles(  CalJet4Pt_, qCalJet4Pt_, quantiles_);
  computeQuantiles(  CalJet5Pt_, qCalJet5Pt_, quantiles_);

  computeQuantiles(  CalJet0Eta_, qCalJet0Eta_, quantiles_);
  computeQuantiles(  CalJet1Eta_, qCalJet1Eta_, quantiles_);
  computeQuantiles(  CalJet2Eta_, qCalJet2Eta_, quantiles_);
  computeQuantiles(  CalJet3Eta_, qCalJet3Eta_, quantiles_);
  computeQuantiles(  CalJet4Eta_, qCalJet4Eta_, quantiles_);
  computeQuantiles(  CalJet5Eta_, qCalJet5Eta_, quantiles_);

  computeQuantiles(  CalJet0Phi_, qCalJet0Phi_, quantiles_);
  computeQuantiles(  CalJet1Phi_, qCalJet1Phi_, quantiles_);
  computeQuantiles(  CalJet2Phi_, qCalJet2Phi_, quantiles_);
  computeQuantiles(  CalJet3Phi_, qCalJet3Phi_, quantiles_);
  computeQuantiles(  CalJet4Phi_, qCalJet4Phi_, quantiles_);
  computeQuantiles(  CalJet5Phi_, qCalJet5Phi_, quantiles_);

  computeQuantiles(  CalJet0En_, qCalJet0En_, quantiles_);
  computeQuantiles(  CalJet1En_, qCalJet1En_, quantiles_);
  computeQuantiles(  CalJet2En_, qCalJet2En_, quantiles_);
  computeQuantiles(  CalJet3En_, qCalJet3En_, quantiles_);
  computeQuantiles(  CalJet4En_, qCalJet4En_, quantiles_);
  computeQuantiles(  CalJet5En_, qCalJet5En_, quantiles_);


  //photon sorted variables (not all of them)
  computeQuantiles(  Pho0Pt_, qPho0Pt_, quantiles_);
  computeQuantiles(  Pho1Pt_, qPho1Pt_, quantiles_);
  computeQuantiles(  Pho2Pt_, qPho2Pt_, quantiles_);
  computeQuantiles(  Pho3Pt_, qPho3Pt_, quantiles_);
  computeQuantiles(  Pho4Pt_, qPho4Pt_, quantiles_);
  computeQuantiles(  Pho5Pt_, qPho5Pt_, quantiles_);

  computeQuantiles(  Pho0Eta_, qPho0Eta_, quantiles_);
  computeQuantiles(  Pho1Eta_, qPho1Eta_, quantiles_);
  computeQuantiles(  Pho2Eta_, qPho2Eta_, quantiles_);
  computeQuantiles(  Pho3Eta_, qPho3Eta_, quantiles_);
  computeQuantiles(  Pho4Eta_, qPho4Eta_, quantiles_);
  computeQuantiles(  Pho5Eta_, qPho5Eta_, quantiles_);

  computeQuantiles(  Pho0Phi_, qPho0Phi_, quantiles_);
  computeQuantiles(  Pho1Phi_, qPho1Phi_, quantiles_);
  computeQuantiles(  Pho2Phi_, qPho2Phi_, quantiles_);
  computeQuantiles(  Pho3Phi_, qPho3Phi_, quantiles_);
  computeQuantiles(  Pho4Phi_, qPho4Phi_, quantiles_);
  computeQuantiles(  Pho5Phi_, qPho5Phi_, quantiles_);

  computeQuantiles(  Pho0En_, qPho0En_, quantiles_);
  computeQuantiles(  Pho1En_, qPho1En_, quantiles_);
  computeQuantiles(  Pho2En_, qPho2En_, quantiles_);
  computeQuantiles(  Pho3En_, qPho3En_, quantiles_);
  computeQuantiles(  Pho4En_, qPho4En_, quantiles_);
  computeQuantiles(  Pho5En_, qPho5En_, quantiles_);

  //ged computeQuantiles(, qPhotons sorted variables (not all of them)
  computeQuantiles(  gedPho0Pt_, qgedPho0Pt_, quantiles_);
  computeQuantiles(  gedPho1Pt_, qgedPho1Pt_, quantiles_);
  computeQuantiles(  gedPho2Pt_, qgedPho2Pt_, quantiles_);
  computeQuantiles(  gedPho3Pt_, qgedPho3Pt_, quantiles_);
  computeQuantiles(  gedPho4Pt_, qgedPho4Pt_, quantiles_);
  computeQuantiles(  gedPho5Pt_, qgedPho5Pt_, quantiles_);
 
  computeQuantiles(  gedPho0Eta_, qgedPho0Eta_, quantiles_);
  computeQuantiles(  gedPho1Eta_, qgedPho1Eta_, quantiles_);
  computeQuantiles(  gedPho2Eta_, qgedPho2Eta_, quantiles_);
  computeQuantiles(  gedPho3Eta_, qgedPho3Eta_, quantiles_);
  computeQuantiles(  gedPho4Eta_, qgedPho4Eta_, quantiles_);
  computeQuantiles(  gedPho5Eta_, qgedPho5Eta_, quantiles_);

  computeQuantiles(  gedPho0Phi_, qgedPho0Phi_, quantiles_);
  computeQuantiles(  gedPho1Phi_, qgedPho1Phi_, quantiles_);
  computeQuantiles(  gedPho2Phi_, qgedPho2Phi_, quantiles_);
  computeQuantiles(  gedPho3Phi_, qgedPho3Phi_, quantiles_);
  computeQuantiles(  gedPho4Phi_, qgedPho4Phi_, quantiles_);
  computeQuantiles(  gedPho5Phi_, qgedPho5Phi_, quantiles_);

  computeQuantiles(  gedPho0En_, qgedPho0En_, quantiles_);
  computeQuantiles(  gedPho1En_, qgedPho1En_, quantiles_);
  computeQuantiles(  gedPho2En_, qgedPho2En_, quantiles_);
  computeQuantiles(  gedPho3En_, qgedPho3En_, quantiles_);
  computeQuantiles(  gedPho4En_, qgedPho4En_, quantiles_);
  computeQuantiles(  gedPho5En_, qgedPho5En_, quantiles_);


    //Muon sorted variables (not all of them)
  computeQuantiles(  Mu0Pt_, qMu0Pt_, quantiles_);
  computeQuantiles(  Mu1Pt_, qMu1Pt_, quantiles_);
  computeQuantiles(  Mu2Pt_, qMu2Pt_, quantiles_);
  computeQuantiles(  Mu3Pt_, qMu3Pt_, quantiles_);
  computeQuantiles(  Mu4Pt_, qMu4Pt_, quantiles_);
  computeQuantiles(  Mu5Pt_, qMu5Pt_, quantiles_);

  computeQuantiles(  Mu0Eta_, qMu0Eta_, quantiles_);
  computeQuantiles(  Mu1Eta_, qMu1Eta_, quantiles_);
  computeQuantiles(  Mu2Eta_, qMu2Eta_, quantiles_);
  computeQuantiles(  Mu3Eta_, qMu3Eta_, quantiles_);
  computeQuantiles(  Mu4Eta_, qMu4Eta_, quantiles_);
  computeQuantiles(  Mu5Eta_, qMu5Eta_, quantiles_);

  computeQuantiles(  Mu0Phi_, qMu0Phi_, quantiles_);
  computeQuantiles(  Mu1Phi_, qMu1Phi_, quantiles_);
  computeQuantiles(  Mu2Phi_, qMu2Phi_, quantiles_);
  computeQuantiles(  Mu3Phi_, qMu3Phi_, quantiles_);
  computeQuantiles(  Mu4Phi_, qMu4Phi_, quantiles_);
  computeQuantiles(  Mu5Phi_, qMu5Phi_, quantiles_);

  computeQuantiles(  Mu0En_, qMu0En_, quantiles_);
  computeQuantiles(  Mu1En_, qMu1En_, quantiles_);
  computeQuantiles(  Mu2En_, qMu2En_, quantiles_);
  computeQuantiles(  Mu3En_, qMu3En_, quantiles_);
  computeQuantiles(  Mu4En_, qMu4En_, quantiles_);
  computeQuantiles(  Mu5En_, qMu5En_, quantiles_);
  //Muon sorted cosmic variables (not all of them)
  computeQuantiles(  MuCosm0Pt_, qMuCosm0Pt_, quantiles_);
  computeQuantiles(  MuCosm1Pt_, qMuCosm1Pt_, quantiles_);
  computeQuantiles(  MuCosm2Pt_, qMuCosm2Pt_, quantiles_);
  computeQuantiles(  MuCosm3Pt_, qMuCosm3Pt_, quantiles_);
  computeQuantiles(  MuCosm4Pt_, qMuCosm4Pt_, quantiles_);
  computeQuantiles(  MuCosm5Pt_, qMuCosm5Pt_, quantiles_);

  computeQuantiles(  MuCosm0Eta_, qMuCosm0Eta_, quantiles_);
  computeQuantiles(  MuCosm1Eta_, qMuCosm1Eta_, quantiles_);
  computeQuantiles(  MuCosm2Eta_, qMuCosm2Eta_, quantiles_);
  computeQuantiles(  MuCosm3Eta_, qMuCosm3Eta_, quantiles_);
  computeQuantiles(  MuCosm4Eta_, qMuCosm4Eta_, quantiles_);
  computeQuantiles(  MuCosm5Eta_, qMuCosm5Eta_, quantiles_);

  computeQuantiles(  MuCosm0Phi_, qMuCosm0Phi_, quantiles_);
  computeQuantiles(  MuCosm1Phi_, qMuCosm1Phi_, quantiles_);
  computeQuantiles(  MuCosm2Phi_, qMuCosm2Phi_, quantiles_);
  computeQuantiles(  MuCosm3Phi_, qMuCosm3Phi_, quantiles_);
  computeQuantiles(  MuCosm4Phi_, qMuCosm4Phi_, quantiles_);
  computeQuantiles(  MuCosm5Phi_, qMuCosm5Phi_, quantiles_);

  computeQuantiles(  MuCosm0En_, qMuCosm0En_, quantiles_);
  computeQuantiles(  MuCosm1En_, qMuCosm1En_, quantiles_);
  computeQuantiles(  MuCosm2En_, qMuCosm2En_, quantiles_);
  computeQuantiles(  MuCosm3En_, qMuCosm3En_, quantiles_);
  computeQuantiles(  MuCosm4En_, qMuCosm4En_, quantiles_);
  computeQuantiles(  MuCosm5En_, qMuCosm5En_, quantiles_);


  //Muon Cosmic1Leg sorted variables (not all of them)
  computeQuantiles(  MuCosmLeg0Pt_, qMuCosmLeg0Pt_, quantiles_);
  computeQuantiles(  MuCosmLeg1Pt_, qMuCosmLeg1Pt_, quantiles_);
  computeQuantiles(  MuCosmLeg2Pt_, qMuCosmLeg2Pt_, quantiles_);
  computeQuantiles(  MuCosmLeg3Pt_, qMuCosmLeg3Pt_, quantiles_);
  computeQuantiles(  MuCosmLeg4Pt_, qMuCosmLeg4Pt_, quantiles_);
  computeQuantiles(  MuCosmLeg5Pt_, qMuCosmLeg5Pt_, quantiles_);

  computeQuantiles(  MuCosmLeg0Eta_, qMuCosmLeg0Eta_, quantiles_);
  computeQuantiles(  MuCosmLeg1Eta_, qMuCosmLeg1Eta_, quantiles_);
  computeQuantiles(  MuCosmLeg2Eta_, qMuCosmLeg2Eta_, quantiles_);
  computeQuantiles(  MuCosmLeg3Eta_, qMuCosmLeg3Eta_, quantiles_);
  computeQuantiles(  MuCosmLeg4Eta_, qMuCosmLeg4Eta_, quantiles_);
  computeQuantiles(  MuCosmLeg5Eta_, qMuCosmLeg5Eta_, quantiles_);

  computeQuantiles(  MuCosmLeg0Phi_, qMuCosmLeg0Phi_, quantiles_);
  computeQuantiles(  MuCosmLeg1Phi_, qMuCosmLeg1Phi_, quantiles_);
  computeQuantiles(  MuCosmLeg2Phi_, qMuCosmLeg2Phi_, quantiles_);
  computeQuantiles(  MuCosmLeg3Phi_, qMuCosmLeg3Phi_, quantiles_);
  computeQuantiles(  MuCosmLeg4Phi_, qMuCosmLeg4Phi_, quantiles_);
  computeQuantiles(  MuCosmLeg5Phi_, qMuCosmLeg5Phi_, quantiles_);

  computeQuantiles(  MuCosmLeg0En_, qMuCosmLeg0En_, quantiles_);
  computeQuantiles(  MuCosmLeg1En_, qMuCosmLeg1En_, quantiles_);
  computeQuantiles(  MuCosmLeg2En_, qMuCosmLeg2En_, quantiles_);
  computeQuantiles(  MuCosmLeg3En_, qMuCosmLeg3En_, quantiles_);
  computeQuantiles(  MuCosmLeg4En_, qMuCosmLeg4En_, quantiles_);
  computeQuantiles(  MuCosmLeg5En_, qMuCosmLeg5En_, quantiles_);




//__________________________________________Sorted variables end here.  

  computeQuantiles(PFJet4CHSPt_, qPFJet4CHSPt_, quantiles_);
  computeQuantiles(PFJet4CHSEta_,qPFJet4CHSEta_,quantiles_);
  computeQuantiles(PFJet4CHSPhi_,qPFJet4CHSPhi_,quantiles_);

  computeQuantiles(PFJet8CHSPt_, qPFJet8CHSPt_, quantiles_);
  computeQuantiles(PFJet8CHSEta_,qPFJet8CHSEta_,quantiles_);
  computeQuantiles(PFJet8CHSPhi_,qPFJet8CHSPhi_,quantiles_);

  computeQuantiles(PFJetEIPt_, qPFJetEIPt_, quantiles_);
  computeQuantiles(PFJetEIEta_,qPFJetEIEta_,quantiles_);
  computeQuantiles(PFJetEIPhi_,qPFJetEIPhi_,quantiles_);

  computeQuantiles(PFJet8CHSSDPt_, qPFJet8CHSSDPt_, quantiles_);
  computeQuantiles(PFJet8CHSSDEta_,qPFJet8CHSSDEta_,quantiles_);
  computeQuantiles(PFJet8CHSSDPhi_,qPFJet8CHSSDPhi_,quantiles_);

  computeQuantiles(PFJetTopCHSPt_, qPFJetTopCHSPt_, quantiles_);
  computeQuantiles(PFJetTopCHSEta_,qPFJetTopCHSEta_,quantiles_);
  computeQuantiles(PFJetTopCHSPhi_,qPFJetTopCHSPhi_,quantiles_);      

  computeQuantiles(PFChMetPt_, qPFChMetPt_,     quantiles_);
  // computeQuantiles(PFChMetEta_,qPFChMetEta_,    quantiles_);  
  computeQuantiles(PFChMetPhi_,qPFChMetPhi_,    quantiles_);
  computeQuantiles(PFMetPt_, qPFMetPt_,     quantiles_);
  // computeQuantiles(PFMetEta_,qPFMetEta_,    quantiles_);
  computeQuantiles(PFMetPhi_,qPFMetPhi_,    quantiles_);
  computeQuantiles(nVtx_,    qNVtx_,    quantiles_);

  computeQuantiles(CalJetPt_, qCalJetPt_, quantiles_);
  computeQuantiles(CalJetEta_,qCalJetEta_,quantiles_);
  computeQuantiles(CalJetPhi_,qCalJetPhi_,quantiles_);
  computeQuantiles(CalJetEn_,qCalJetEn_,  quantiles_);
  computeQuantiles(CalMETPt_, qCalMETPt_, quantiles_);
  // computeQuantiles(CalMETEta_,qCalMETEta_,quantiles_);
  computeQuantiles(CalMETPhi_,qCalMETPhi_,quantiles_);
  computeQuantiles(CalMETEn_,qCalMETEn_,  quantiles_);

  computeQuantiles(CalMETBEPt_, qCalMETBEPt_, quantiles_);
  // computeQuantiles(CalMETBEEta_,qCalMETBEEta_,quantiles_);
  computeQuantiles(CalMETBEPhi_,qCalMETBEPhi_,quantiles_);
  computeQuantiles(CalMETBEEn_, qCalMETBEEn_,  quantiles_);

  computeQuantiles(CalMETBEFOPt_, qCalMETBEFOPt_, quantiles_);
  // computeQuantiles(CalMETBEFOEta_,qCalMETBEFOEta_,quantiles_);
  computeQuantiles(CalMETBEFOPhi_,qCalMETBEFOPhi_,quantiles_);
  computeQuantiles(CalMETBEFOEn_, qCalMETBEFOEn_,  quantiles_);

  computeQuantiles(CalMETMPt_, qCalMETMPt_, quantiles_);
  // computeQuantiles(CalMETMEta_,qCalMETMEta_,quantiles_);
  computeQuantiles(CalMETMPhi_,qCalMETMPhi_,quantiles_);
  computeQuantiles(CalMETMEn_, qCalMETMEn_,  quantiles_);

  computeQuantiles(SCEn_, qSCEn_,       quantiles_);
  computeQuantiles(SCEta_, qSCEta_,     quantiles_);
  computeQuantiles(SCPhi_, qSCPhi_,     quantiles_);
  computeQuantiles(SCEtaWidth_, qSCEtaWidth_,     quantiles_);
  computeQuantiles(SCPhiWidth_, qSCPhiWidth_,     quantiles_);

  computeQuantiles(SCEnhfEM_, qSCEnhfEM_,       quantiles_);
  computeQuantiles(SCEtahfEM_, qSCEtahfEM_,     quantiles_);
  computeQuantiles(SCPhihfEM_, qSCPhihfEM_,     quantiles_);
  // computeQuantiles(SCEtaWidthhfEM_, qSCEtaWidthhfEM_,     quantiles_);
  // computeQuantiles(SCPhiWidthhfEM_, qSCPhiWidthhfEM_,     quantiles_);

  computeQuantiles(SCEn5x5_, qSCEn5x5_,       quantiles_);
  computeQuantiles(SCEta5x5_, qSCEta5x5_,     quantiles_);
  computeQuantiles(SCPhi5x5_, qSCPhi5x5_,     quantiles_);
  computeQuantiles(SCEtaWidth5x5_, qSCEtaWidth5x5_,     quantiles_);
  computeQuantiles(SCPhiWidth5x5_, qSCPhiWidth5x5_,     quantiles_);

  computeQuantiles(CCEn_, qCCEn_,       quantiles_);
  computeQuantiles(CCEta_, qCCEta_,     quantiles_);
  computeQuantiles(CCPhi_, qCCPhi_,     quantiles_);
  computeQuantiles(CCEn5x5_, qCCEn5x5_,       quantiles_);
  computeQuantiles(CCEta5x5_, qCCEta5x5_,     quantiles_);
  computeQuantiles(CCPhi5x5_, qCCPhi5x5_,     quantiles_);

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
  computeQuantiles(MuChi2_,qMuChi2_,quantiles_);  
  computeQuantiles(MuCosmPt_, qMuCosmPt_, quantiles_);
  computeQuantiles(MuCosmEta_,qMuCosmEta_,quantiles_);
  computeQuantiles(MuCosmPhi_,qMuCosmPhi_,quantiles_);
  computeQuantiles(MuCosmEn_, qMuCosmEn_,quantiles_);
  computeQuantiles(MuCosmCh_, qMuCosmCh_,quantiles_);
  computeQuantiles(MuCosmChi2_,qMuCosmChi2_,quantiles_);   
  computeQuantiles(MuCosmLegPt_, qMuCosmLegPt_, quantiles_);
  computeQuantiles(MuCosmLegEta_,qMuCosmLegEta_,quantiles_);
  computeQuantiles(MuCosmLegPhi_,qMuCosmLegPhi_,quantiles_);
  computeQuantiles(MuCosmLegEn_, qMuCosmLegEn_,quantiles_);
  computeQuantiles(MuCosmLegCh_, qMuCosmLegCh_,quantiles_);
  computeQuantiles(MuCosmLegChi2_,qMuCosmLegChi2_,quantiles_);       
  //TODOCosm
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
  computeQuantiles(EBiEta_, qEBiEta_, quantiles_);
  computeQuantiles(EBiPhi_, qEBiPhi_, quantiles_);


  computeQuantiles(EEenergy_, qEEenergy_, quantiles_);
  computeQuantiles(EEtime_, qEEtime_, quantiles_);
  computeQuantiles(EEchi2_, qEEchi2_, quantiles_);
  computeQuantiles(EEix_, qEEix_, quantiles_);
  computeQuantiles(EEiy_, qEEiy_, quantiles_);

  computeQuantiles(ESenergy_, qESenergy_, quantiles_);
  computeQuantiles(EStime_, qEStime_, quantiles_);
  // computeQuantiles(ESchi2_, qESchi2_, quantiles_);
  computeQuantiles(ESix_, qESix_, quantiles_);
  computeQuantiles(ESiy_, qESiy_, quantiles_);

  computeQuantiles(HBHEenergy_, qHBHEenergy_, quantiles_);
  computeQuantiles(HBHEtime_, qHBHEtime_, quantiles_);
  computeQuantiles(HBHEauxe_, qHBHEauxe_, quantiles_);
  computeQuantiles(HBHEieta_, qHBHEieta_, quantiles_);
  computeQuantiles(HBHEiphi_, qHBHEiphi_, quantiles_);

  computeQuantiles(HFenergy_, qHFenergy_, quantiles_);
  computeQuantiles(HFtime_, qHFtime_, quantiles_);
  // computeQuantiles(HFchi2_, qHFchi2_, quantiles_);
  computeQuantiles(HFieta_, qHFieta_, quantiles_);
  computeQuantiles(HFiphi_, qHFiphi_, quantiles_);

  // computeQuantiles(HOenergy_, qHOenergy_, quantiles_);
  // computeQuantiles(HOtime_, qHOtime_, quantiles_);
  // // computeQuantiles(HOchi2_, qHOchi2_, quantiles_);
  // computeQuantiles(HOieta_, qHOieta_, quantiles_);
  // computeQuantiles(HOiphi_, qHOiphi_, quantiles_);  


  computeQuantiles(PreShEn_, qPreShEn_,         quantiles_);
  // computeQuantiles(PreShCorrEn_, qPreShCorrEn_,   quantiles_);
  computeQuantiles(PreShEta_, qPreShEta_,       quantiles_);
  computeQuantiles(PreShPhi_, qPreShPhi_,       quantiles_);
  computeQuantiles(PreShYEn_, qPreShYEn_,       quantiles_);
  // computeQuantiles(PreShYCorrEn_, qPreShYCorrEn_,   quantiles_);
  computeQuantiles(PreShYEta_, qPreShYEta_,     quantiles_);
  computeQuantiles(PreShYPhi_, qPreShYPhi_,     quantiles_);

  // computeQuantiles(CTPt_, qCTPt_, quantiles_);
  // computeQuantiles(CTEta_,qCTEta_,quantiles_);
  // computeQuantiles(CTPhi_,qCTPhi_,quantiles_);

  crossSection_->push_back( (float)eventCounter/lumi_ );
  
  for(std::map<std::string,int>::const_iterator itr = rateMap.begin(); itr != rateMap.end(); ++itr)
    {
      pathNames_->push_back(itr->first);
      pathRates_->push_back(itr->second/lumi_);
       // std::cout << "pathRates_: " << itr->second/lumi_ << std::endl; //TEST -- works
       // std::cout << "pathNames_: " << itr->first << std::endl; //TEST -- works


    }


  //fill tree one event per LS
  outTree_->Fill();


}



void AODAnalyzer::beginRun (const edm::Run &run, const edm::EventSetup &eventSetup) 
{
 bool isConfigChanged = false;
 edm::InputTag myTag("TriggerResults::HLT");
 hltConfigProvider_.init(run, eventSetup, myTag.process(), isConfigChanged);
}


void AODAnalyzer::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) 
{
  ++eventCounter;



  //fill Jets
  edm::Handle<reco::PFJetCollection> PFJets;
  event.getByToken(PFJetToken_,PFJets);
  if(PFJets.isValid())
    fillJets(PFJets, std::string("PF")); 

//---------------------------------Sorted variables start here.

  // //fill SortedJets
  // edm::Handle<reco::PFJetCollection> PFSortedJets;
  // event.getByToken(PFJetToken_,PFSortedJets);
  // if(PFSortedJets.isValid())  
  //   fillSortedJets(PFSortedJets, std::string("PF"));





//_________________________________Sorted variables end here.




  //fill PFJet4CHS Jets
  edm::Handle<reco::PFJetCollection> PF4CHSJets;
  event.getByToken(PFJet4CHSToken_,PF4CHSJets);
  if(PF4CHSJets.isValid())
    fill4CHSJets(PF4CHSJets, std::string("PF"));

    //fill PFJet8CHS Jets
  edm::Handle<reco::PFJetCollection> PF8CHSJets;
  event.getByToken(PFJet8CHSToken_,PF8CHSJets);
  if(PF8CHSJets.isValid())
    fill8CHSJets(PF8CHSJets, std::string("PF"));

    //fill PFJetEI Jets
  edm::Handle<reco::PFJetCollection> PFEIJets;
  event.getByToken(PFJetEIToken_,PFEIJets);
  if(PFEIJets.isValid())
    fillEIJets(PFEIJets, std::string("PF"));

    //fill PFJet8CHSSoftDrop Jets
  edm::Handle<reco::PFJetCollection> PF8CHSSoftDropJets;
  event.getByToken(PFJet8CHSSoftDropToken_,PF8CHSSoftDropJets);
  if(PF8CHSSoftDropJets.isValid())
    fill8CHSoftDropJets(PF8CHSSoftDropJets, std::string("PF"));

    //fill PFJetTopCHS Jets
  edm::Handle<reco::PFJetCollection> PFTopCHSJets;
  event.getByToken(PFJetTopCHSToken_,PFTopCHSJets);
  if(PFTopCHSJets.isValid())
    fillTopCHSJets(PFTopCHSJets, std::string("PF"));
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

  //fill calomet
  edm::Handle<reco::CaloMETCollection> caloMETlocalv;
  event.getByToken(CaloMETToken_, caloMETlocalv);
  if(caloMETlocalv.isValid())
    fillCaloMETs(caloMETlocalv);

  //fill calomet BE
  edm::Handle<reco::CaloMETCollection> caloMETBElocalv;
  event.getByToken(CaloMETBEToken_, caloMETBElocalv);
  if(caloMETBElocalv.isValid())
    fillCaloMETBEs(caloMETBElocalv);

    //fill calomet BEFO
  edm::Handle<reco::CaloMETCollection> caloMETBEFOlocalv;
  event.getByToken(CaloMETBEFOToken_, caloMETBEFOlocalv);
  if(caloMETBEFOlocalv.isValid())
    fillCaloMETBEFOs(caloMETBEFOlocalv);

    //fill calomet M
  edm::Handle<reco::CaloMETCollection> caloMETMlocalv;
  event.getByToken(CaloMETMToken_, caloMETMlocalv);
  if(caloMETMlocalv.isValid())
    fillCaloMETMs(caloMETMlocalv);

  //fill vtx
  edm::Handle<reco::VertexCollection> recVtxs;
  event.getByToken(vtxToken_,recVtxs);
  if(recVtxs.isValid())
    nVtx_->push_back(recVtxs->size());


  //Bits
  edm::Handle<edm::TriggerResults> triggerBits;
  event.getByToken(triggerBits_, triggerBits);


  // Prescales (alternative for pat::PackedTriggerPrescales in order to get pathRates)
  // edm::Handle<edm::TriggerResults> triggerResults;
  // event.getByToken(triggerResultsToken_, triggerResults);

  // edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  // event.getByToken(triggerPrescales_, triggerPrescales);

  // edm::Handle<trigger::TriggerEvent> triggerEvent;
  // event.getByToken( "patTriggerEvent", triggerEvent );
  

  //Fill SuperCluster
  edm::Handle<reco::SuperClusterCollection> SuperClusterlocalv;
  event.getByToken(SuperClusterToken_, SuperClusterlocalv);
  if(SuperClusterlocalv.isValid())
    fillSC(SuperClusterlocalv);

  edm::Handle<reco::SuperClusterCollection> SuperClusterhfEMlocalv;
  event.getByToken(SuperClusterhfEMToken_, SuperClusterhfEMlocalv);
  if(SuperClusterhfEMlocalv.isValid())
    fillSChfEM(SuperClusterhfEMlocalv);

  edm::Handle<reco::SuperClusterCollection> SuperCluster5x5localv;
  event.getByToken(SuperCluster5x5Token_, SuperCluster5x5localv);
  if(SuperCluster5x5localv.isValid())
    fillSC5x5(SuperCluster5x5localv);

    //Fill CaloCluster
  edm::Handle<reco::CaloClusterCollection> CaloClusterlocalv;
  event.getByToken(CaloClusterToken_, CaloClusterlocalv);
  if(CaloClusterlocalv.isValid())
    fillCC(CaloClusterlocalv);

  edm::Handle<reco::CaloClusterCollection> CaloCluster5x5localv;
  event.getByToken(CaloCluster5x5Token_, CaloCluster5x5localv);
  if(CaloCluster5x5localv.isValid())
    fillCC5x5(CaloCluster5x5localv);

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

  edm::Handle<reco::MuonCollection> muonCosmlocalv;
  event.getByToken(MuonCosmToken_, muonCosmlocalv);
  if(muonCosmlocalv.isValid())
    fillCosmMuons(muonCosmlocalv);

  edm::Handle<reco::MuonCollection> muonCosmLeglocalv;
  event.getByToken(MuonCosmLegToken_, muonCosmLeglocalv);
  if(muonCosmLeglocalv.isValid())
    fillCosmLegMuons(muonCosmLeglocalv);

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
  edm::Handle<EcalRecHitCollection> ebRHs;   
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
  edm::Handle<HBHERecHitCollection> hbheRHs;    
  event.getByToken(hbheRHcToken_, hbheRHs);
  if(hbheRHs.isValid())
    fillHBHErecHit(hbheRHs);

  //fill hfreco
  edm::Handle<HFRecHitCollection> hfRHs;
  event.getByToken(hfRHcToken_, hfRHs);
  if(hfRHs.isValid())
    fillHFrecHit(hfRHs);

  //fill horeco
  // edm::Handle<HORecHitCollection> hoRHs;
  // event.getByToken(hoRHcToken_, hoRHs);
  // if(hoRHs.isValid())
  //   fillHOrecHit(hoRHs);

  // fill PreshowerCluster
  edm::Handle<reco::PreshowerClusterCollection> prShs;
  event.getByToken(preshowerXToken_, prShs);
  if(prShs.isValid())
    fillPreshowerCluster(prShs);

  edm::Handle<reco::PreshowerClusterCollection> prShYs;
  event.getByToken(preshowerYToken_, prShYs);
  if(prShYs.isValid())
    fillPreshowerClusterY(prShYs);

  // edm::Handle<reco::CastorTowerCollection> castorsz;
  // event.getByToken(CastorTowerToken_, castorsz);
  // if(castorsz.isValid())
  //   fillCastorTower(castorsz);


  // -----------------------------------------------------------------
  // --------Hunt for pathRates---------------------------------------
    // This was in https://github.com/deguio/Analyzers/blob/master/miniAODAnalyzer/plugins/miniAODAnalyzer.cc#L518
   //=====================================================================================================================================
  // const edm::TriggerNames &names = event.triggerNames(*triggerBits);
  // for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) 
  //   {
  //     if(rateMap.find(names.triggerName(i)) != rateMap.end())
  // rateMap[names.triggerName(i)] += triggerPrescales->getPrescaleForIndex(i)*triggerBits->accept(i);
  //     else
  // rateMap[names.triggerName(i)] = triggerPrescales->getPrescaleForIndex(i)*triggerBits->accept(i);

  //     std::cout << names.triggerName(i) << " " << triggerPrescales->getPrescaleForIndex(i) << " " << triggerBits->accept(i) << std::endl;
  // ======================================================================================================================================

 //PLAYGROUND-------------------------

   const edm::TriggerNames &names = event.triggerNames(*triggerBits);
   const std::vector<std::string>& triggerNames = hltConfigProvider_.triggerNames();
   for (size_t ts = 0; ts < triggerNames.size(); ts++) {
       std::string trig = triggerNames[ts]; //--adding because of test

        // std::cout << "HLT name " << trig;  //adding because of test

    const unsigned int prescaleSize = hltConfigProvider_.prescaleSize();
        for (unsigned int ps = 0; ps < prescaleSize; ps++) 
          {
            const unsigned int prescaleValue = hltConfigProvider_.prescaleValue(ps, trig);
               if(rateMap.find(names.triggerName(ts)) != rateMap.end())
                  rateMap[names.triggerName(ts)] += prescaleValue;
               else
                  rateMap[names.triggerName(ts)] = prescaleValue;
                // std::cout << prescaleValue << " ";  //adding because of test
          }
          // std::cout << std::endl; //adding because of test
      }





 //PLAYGROUND ENDS...-------------------------------------
  // FEDES CODE
  //  const std::vector<std::string>& triggerNames =
  //       hltConfigProvider_.triggerNames();
  //   for (size_t ts = 0; ts < triggerNames.size(); ts++) {
  //     std::string trig = triggerNames[ts];

  //       std::cout << "HLT name " << trig;
  //       // See if the trigger is prescaled;
  //       /// number of prescale sets available
  //       const unsigned int prescaleSize = hltConfigProvider_.prescaleSize();
  //       for (unsigned int ps = 0; ps < prescaleSize; ps++) 
  //         {
  //           const unsigned int prescaleValue = hltConfigProvider_.prescaleValue(ps, trig);
  //               std::cout << prescaleValue << " ";
  //         //  std::cout<< " prescaleValue[" << ps << "] =" << prescaleValue
  //           //<<std::endl;
    
  //          }
  //         std::cout << std::endl;
  // }
  // FEDES CODE
  // --------------------------------------------------------------
  // ---------------------------------------------------------------
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
