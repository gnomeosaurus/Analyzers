import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                                   #"/store/data/Run2015D/SingleElectron/MINIAOD/16Dec2015-v1/20000/00050EF1-F9A6-E511-86B2-0025905A48D0.root", #going for SuperCluster in MINIAOD
				                           # '/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v3/000/284/040/00000/76F4EC3E-879F-E611-8792-02163E012498.root', #randomly choosing file regarding  electrons
                                  '/store/data/Run2016H/DoubleEG/AOD/PromptReco-v1/000/281/130/00000/70B9D8E1-7680-E611-ABA1-02163E014716.root', #exercise2 - going for other variables in AOD
                                   # '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/08080C9F-C078-E611-AEF1-FA163E5647FC.root',
                                   # '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/084B4FBE-C078-E611-9740-02163E013993.root',
                                   # '/store/data/Run2016G/ZeroBias/MINIAOD/PromptReco-v1/000/280/385/00000/0A6B1BD0-C078-E611-A332-02163E014113.root'
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
process.MyAnalysis =cms.EDAnalyzer("MyMiniAODAnalyzer",

                       #caloJetTag              = cms.untracked.InputTag("slimmedJets"),
                       PFJetTag                = cms.untracked.InputTag("slimmedJets"),
                       metTag                  = cms.untracked.InputTag("slimmedMETs"),
                       vtx                     = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                       bits                    = cms.untracked.InputTag("TriggerResults","","HLT"),
                       prescales               = cms.untracked.InputTag("patTrigger"), #this is giving us a crash...

                       SuperClusterTag         = cms.untracked.InputTag("reducedEgamma","reducedSuperClusters"), #adding SuperCluster --needs to be changed because miniAOD!=AOD

                       #44-57 are new tags
                       GsfElectronTag          = cms.untracked.InputTag("gedGsfElectrons"),
                       GsfElectronUncleanedTag = cms.untracked.InputTag("uncleanedOnlyGsfElectrons"),
                       MuonTag                 = cms.untracked.InputTag("muons"),
                       gedPhotonTag            = cms.untracked.InputTag("gedPhotons"),
                       PhotonTag               = cms.untracked.InputTag("photons"),
                       ChPFMETTag              = cms.untracked.InputTag("pfChMet"),
                       PFMETTag                = cms.untracked.InputTag("pfMet"),

                       CollectionEcalRecHitEBTag  = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                       CollectionEcalRecHitEETag  = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                       CollectionEcalRecHitESTag  = cms.untracked.InputTag("reducedEcalRecHitsES"),
                       CollectionHBHERecHitTag    = cms.untracked.InputTag("reducedHcalRecHits","hbhereco"),
                       CollectionHFRecHitTag      = cms.untracked.InputTag("reducedHcalRecHits","hfreco"),
                       CollectionHORecHitTag      = cms.untracked.InputTag("reducedHcalRecHits","horeco"),

                                   
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
                                   fileName = cms.string("miniAODTree.root"),
                                   closeFileFast = cms.untracked.bool(False),
                                   )
 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options =  cms.untracked.PSet(
                   #allowUnscheduled = cms.untracked.bool(True),
                   wantSummary = cms.untracked.bool(True),
                   SkipEvent = cms.untracked.vstring('ProductNotFound'), #!!! only for testing
                   )


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
    debugModules  = cms.untracked.vstring('*')
)
