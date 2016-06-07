import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

import glob
fileList = glob.glob("/tmp/deguio/HLTPhysics_Run274200/myResult*.root")
finalList = []
for fil in fileList:
    finalList.append("file:"+fil)

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                        finalList

                        #'file:/afs/cern.ch/work/d/deguio/Analysis/DiJetScouting/triggerStudies/CMSSW_8_0_9_Mjj_trigger/src/OpenPaths_0606/myResults.root',
                        #'file:/tmp/deguio/HLTPhysics_Run274200/HLTPhysics_Run274200.root',
                        #'file:/tmp/deguio/HLTPhysics_Run274200/myResults_600.root'
                    ),
                    secondaryFileNames = cms.untracked.vstring(),
#                     lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_dataRun2_HLT_v12'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


process.MyAnalysis =cms.EDAnalyzer("MyHLTAnalyzer",
                       triggerResult_1         = cms.untracked.InputTag("TriggerResults::HLT"),
                       triggerResult_2         = cms.untracked.InputTag("TriggerResults::TEST"),
                       triggerSummary          = cms.untracked.InputTag("hltDiCaloWideJetMass200::TEST"),
                       caloJetTag              = cms.untracked.InputTag("hltAK4CaloJetsCorrectedIDPassed::TEST"),
                       PFJetTag                = cms.untracked.InputTag("hltAK4PFJetsCorrected::TEST"),

                       #params for wide jet mass calculation
                       minMass                 = cms.untracked.double(200),  
                       fatJetDeltaR            = cms.untracked.double(1.1),
                       maxDeltaEta             = cms.untracked.double(1.5),
                       maxJetEta               = cms.untracked.double(3.0),
                       minJetPt                = cms.untracked.double(30.0),

                       #paths to monitor
                       hltPaths                = cms.untracked.vstring(
                                                                       #from TriggerResults::HLT
                                                                       'DST_HT250_CaloScouting_v',
                                                                       'DST_HT250_CaloBTagScouting_v',
                                                                       'DST_HT410_PFScouting_v',
                                                                       'DST_HT410_BTagScouting_v',
                                                                       #from TriggerResults::TEST
                                                                       'DST_HT250_CaloScouting_v2_ref',
                                                                       'DST_HT250_CaloBTagScouting_v1_ref',
                                                                       'DST_HT410_PFScouting_v1_ref',
                                                                       'DST_HT410_BTagScouting_v1_ref',
                                                                       # from TriggerResults::TEST
                                                                       'DST_DiPFWideJetMass200_PFScouting_v',
                                                                       'DST_DiPFWideJetMass200_BTagScouting_v',
                                                                       'DST_DiCaloWideJetMass200_CaloScouting_v',
                                                                       'DST_DiCaloWideJetMass200_CaloBTagScouting_v'
                                                                      )
                       )
                       
process.mypath  = cms.Path(process.MyAnalysis)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hlTree.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))   


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
    debugModules  = cms.untracked.vstring('*')
)
