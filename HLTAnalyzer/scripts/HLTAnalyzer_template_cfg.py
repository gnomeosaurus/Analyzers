import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(INPUTLIST),
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
                       l1CandTag               = cms.untracked.InputTag("hltCaloStage2Digis","Jet","TEST"),
                       caloJetTag              = cms.untracked.InputTag("hltAK4CaloJetsCorrectedIDPassed::TEST"),
                       PFJetTag                = cms.untracked.InputTag("hltAK4PFJetsCorrected::TEST"),
                                   

                       #params for wide jet mass calculation
                       minMass                 = cms.untracked.double(200),  
                       fatJetDeltaR            = cms.untracked.double(1.1),

                       maxDeltaEta             = cms.untracked.double(1.5),
                       maxJetEta               = cms.untracked.double(3.0),
                       minJetPt                = cms.untracked.double(30.0),

                       maxL1DeltaEta           = cms.untracked.double(100.0),
                       maxL1JetEta             = cms.untracked.double(100.0),
                       minL1JetPt              = cms.untracked.double(1.0),

                       #paths to monitor
                       hltPaths                = cms.untracked.vstring(
                                                                       #from TriggerResults::HLT
                                                                       'HLT_ZeroBias_v',
                                                                       'DST_HT250_CaloScouting_v',
                                                                       'DST_HT250_CaloBTagScouting_v',
                                                                       'DST_HT410_PFScouting_v',
                                                                       'DST_HT410_BTagScouting_v',
                                                                       #from TriggerResults::TEST
                                                                       'HLT_ZeroBias_v2_ref',
                                                                       'DST_HT250_CaloScouting_v2_ref',
                                                                       'DST_HT250_CaloBTagScouting_v1_ref',
                                                                       'DST_HT410_PFScouting_v1_ref',
                                                                       'DST_HT410_BTagScouting_v1_ref',
                                                                       # from TriggerResults::TEST
                                                                       'DST_DiPFWideJetMass200_PFScouting_v',
                                                                       'DST_DiPFWideJetMass200_BTagScouting_v',
                                                                       'DST_DiCaloWideJetMass200_CaloScouting_v',
                                                                       'DST_DiCaloWideJetMass200_CaloBTagScouting_v',
                                                                       'DST_L1HTT_CaloScouting_PFScouting_v'
                                                                      ),
                       l1Paths                 = cms.untracked.vstring('L1_HTT200',
                                                                       'L1_HTT240',
                                                                       'L1_HTT270',
                                                                       'L1_HTT280',
                                                                       'L1_HTT300',
                                                                       'L1_HTT320',
                                                                       'L1_ZeroBias',
                                                                       'L1_DoubleJetC80',
                                                                       'L1_DoubleJetC100',
                                                                       'L1_DoubleJetC112',
                                                                       'L1_DoubleIsoTau26er',
                                                                       'L1_DoubleIsoTau28er',
                                                                       'L1_DoubleIsoTau30er',
                                                                       'L1_DoubleIsoTau32er'),
                        AlgInputTag = cms.InputTag("hltGtStage2Digis"),
                        l1tAlgBlkInputTag = cms.InputTag("hltGtStage2Digis"),
                        l1tExtBlkInputTag = cms.InputTag("hltGtStage2Digis")
                       )


#needed if I want to run on raw directly
process.hltCaloStage2Digis = cms.EDProducer("L1TRawToDigi",
    CTP7 = cms.untracked.bool(False),
    FWId = cms.uint32(0),
    FWOverride = cms.bool(False),
    FedIds = cms.vint32(1360, 1366),
    InputLabel = cms.InputTag("rawDataCollector"),
    MTF7 = cms.untracked.bool(False),
    Setup = cms.string('stage2::CaloSetup'),
    debug = cms.untracked.bool(False),
    lenAMC13Header = cms.untracked.int32(8),
    lenAMC13Trailer = cms.untracked.int32(8),
    lenAMCHeader = cms.untracked.int32(8),
    lenAMCTrailer = cms.untracked.int32(0),
    lenSlinkHeader = cms.untracked.int32(8),
    lenSlinkTrailer = cms.untracked.int32(8)
)
                       
process.mypath  = cms.Path(process.MyAnalysis)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(OUTPUTFILE),
                                   closeFileFast = cms.untracked.bool(False)
                                   )
 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options =  cms.untracked.PSet(
                   #allowUnscheduled = cms.untracked.bool(True),
                   wantSummary = cms.untracked.bool(True),
                   )


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
    debugModules  = cms.untracked.vstring('*')
)
