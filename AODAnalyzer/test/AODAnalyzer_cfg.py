import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                                   #"/store/data/Run2015D/SingleElectron/MINIAOD/16Dec2015-v1/20000/00050EF1-F9A6-E511-86B2-0025905A48D0.root", #going for SuperCluster in MINIAOD
				                           # '/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v3/000/284/040/00000/76F4EC3E-879F-E611-8792-02163E012498.root', #randomly choosing file regarding  electrons
                                  #'/store/data/Run2016H/DoubleEG/AOD/PromptReco-v1/000/281/130/00000/70B9D8E1-7680-E611-ABA1-02163E014716.root', #exercise2 - going for other variables in AOD
                                   # '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/08080C9F-C078-E611-AEF1-FA163E5647FC.root',
                                   # '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/084B4FBE-C078-E611-9740-02163E013993.root',
                                   # '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/0A6B1BD0-C078-E611-A332-02163E014113.root',
                                   "/store/data/Run2016H/DoubleEG/AOD/PromptReco-v3/000/284/044/00000/42E8B572-A29F-E611-A76B-FA163E4C7CB3.root",
                        ),
                    secondaryFileNames = cms.untracked.vstring(),
#                   lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v11'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


basePath = '/afs/cern.ch/user/d/deguio/public/qualityPerSub_2016/'
subsystemList = ['L1tcalo','L1tmu','Hlt','Pix','Strip','Ecal','Hcal','Dt','Rpc','Es','Csc','Track','Egamma','Muon','Jetmet','Lumi']
fileList = []
for sub in subsystemList:
    fileList.append(basePath+'Cert_13TeV_2016_'+sub+'.txt')

#vector<reco::SuperCluster>            "reducedEgamma"             "reducedSuperClusters"   "RECO"     --output of edmDumpEventContent for path #going for Supercluster
process.MyAnalysis =cms.EDAnalyzer("AODAnalyzer",

                       
                       PFJetTag                = cms.untracked.InputTag("ak8PFJetsCHS"), #ak4PFJets
                       PFChMETTag              = cms.untracked.InputTag("pfChMet"),
                       PFMETTag                = cms.untracked.InputTag("pfMet"),
                       CaloJetTag              = cms.untracked.InputTag("ak4CaloJets"),
                       CaloMETTag              = cms.untracked.InputTag("caloMet"),  #calometTag

                       vtx                     = cms.untracked.InputTag("offlinePrimaryVertices"),
                       bits                    = cms.untracked.InputTag("TriggerResults","","HLT"),  #ASK IF HLT OR RECO!!!!!!!
                       prescales               = cms.untracked.InputTag("hltTriggerSummaryAOD"), # //PROBABLY get from https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/TriggerEvent.h  ... ask about how

                       SuperClusterTag         = cms.untracked.InputTag("particleFlowEGamma"),
                       SuperClusterhfEMTag     = cms.untracked.InputTag("hfEMClusters"),
                       SuperCluster5x5Tag      = cms.untracked.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                       CaloClusterTag          = cms.untracked.InputTag("particleFlowEGamma","EBEEClusters"),  #ESClusters also possible instead of EBEEClusters
                       CaloCluster5x5Tag       = cms.untracked.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters"),    #CHANGING HEFM TO MULTI5X5  HFEM was empty!

                       PhotonTag               = cms.untracked.InputTag("photons"),
                       gedPhotonTag            = cms.untracked.InputTag("gedPhotons"),
                       MuonTag                 = cms.untracked.InputTag("muons"),
                       MuonCosmTag             = cms.untracked.InputTag("muonsFromCosmics"),
                       MuonCosmLegTag          = cms.untracked.InputTag("muonsFromCosmics1Leg"),

                       GsfElectronTag          = cms.untracked.InputTag("gedGsfElectrons"),
                       GsfElectronUncleanedTag = cms.untracked.InputTag("uncleanedOnlyGsfElectrons"),


                       EBRecHitSourceTag       = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                       EERecHitSourceTag       = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                       ESRecHitSourceTag       = cms.untracked.InputTag("reducedEcalRecHitsES"),

                       HBHERecHitTag           = cms.untracked.InputTag("reducedHcalRecHits","hbhereco"),
                       HFRecHitTag             = cms.untracked.InputTag("reducedHcalRecHits","hfreco"),
                       HORecHitTag             = cms.untracked.InputTag("reducedHcalRecHits","horeco"),
                       PreshowerClusterXTag    = cms.untracked.InputTag("multi5x5SuperClustersWithPreshower","preshowerXClusters"),  
                       PreshowerClusterYTag    = cms.untracked.InputTag("multi5x5SuperClustersWithPreshower","preshowerYClusters"),
                       CastorTowerTag          = cms.untracked.InputTag("CastorTowerReco"),  #add LeafCandidate variables
        
                       maxJetEta               = cms.untracked.double(5.0),   #is this ok?
                       minJetPt                = cms.untracked.double(10.0),  #is this ok?

                       maxSCEta                = cms.untracked.double(3.0),   # IS THIS OK?
                       minSCEn                 = cms.untracked.double(8.0),  # IS THIS OK?

                       lumiFile                = cms.untracked.string(basePath+'run_ls_lumi_2016.txt'),
                       subsystems              = cms.untracked.vstring(subsystemList),
                       qualityFiles            = cms.untracked.vstring(fileList),
                       quantiles               = cms.untracked.vdouble(0.,0.25,0.50,0.75,1.)

                        )


process.mypath  = cms.Path(process.MyAnalysis)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("AODTree.root"),
                                   closeFileFast = cms.untracked.bool(False),
                                   )
 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options =  cms.untracked.PSet(
                   #allowUnscheduled = cms.untracked.bool(True),
                   wantSummary = cms.untracked.bool(True),
                   #SkipEvent = cms.untracked.vstring('ProductNotFound'), #!! only for testing
                   )


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
    debugModules  = cms.untracked.vstring('*')
)
