from WMCore.Configuration import Configuration
# from IPython import embed

config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.workArea = 'crab3final5thsubmission'
config.General.requestName = ''
#config.General.requestName = 'Run2016allruns'

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/user/f/fsiroky/CMSSW_8_0_28/src/Analyzers/AODAnalyzer/test/AODAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['AODTree.root']
config.JobType.maxJobRuntimeMin = 2750 


config.section_('Data')
# config.Data.inputDataset = '/DoubleEG/Run2016H-PromptReco-v1/AOD'
config.Data.inputDataset = ''
#config.Data.inputDataset = '/DoubleEG/Run2016H-18Apr2017-v1/AOD'   #RUN2016A-H  DoubleEG/*/AOD
config.Data.unitsPerJob = 40 #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
#config.Data.runRange = '275656-284044'  
# config.Data.runRange = '281613-284044' //2016H range   //275656 C range starts   // 284044 H range ends
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
#config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
#config.Data.outLFNDirBase = '/store/group/alca_global/filipTemporary/'
config.Data.outLFNDirBase = '/store/user/fsiroky/data_fifthrun/'
config.Data.ignoreLocality = True

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

# WHITELIST DAS LINK: https://cmsweb.cern.ch/das/request?view=list&limit=100&instance=prod%2Fglobal&input=site+dataset%3D%2F*%2FRun2016*-18Apr2017-v1%2FAOD
# config.Site.whitelist = ['T1_DE_KIT',
#                          'T1_ES_PIC',
#                          'T1_IT_CNAF',
#                          'T1_RU_JINR',
#                          'T1_UK_RAL',
#                          'T1_US_FNAL',                       
#                          'T2_AT_Vienna',
#                          'T2_BE_IIHE',
#                          'T2_BE_UCL',
#                          'T2_BR_SPRACE',
#                          'T2_CH_CERN',
#                          'T2_CH_CSCS',
#                          'T2_CN_Beijing',
#                          'T2_DE_DESY',
#                          'T2_DE_RWTH',
#                          'T2_EE_Estonia',
#                          'T2_ES_CIEMAT', 
#                          'T2_FI_HIP',
#                          'T2_FR_IPHC',
#                          'T2_IT_Legnaro',
#                          'T2_IT_Pisa',
#                          'T2_PL_Swierk',
#                          'T2_RU_JINR',
#                          'T2_US_Caltech',
#                          'T2_US_MIT',
#                          'T2_US_Nebraska',
#                          'T2_US_Purdue',
#                          'T2_US_Vanderbilt',
#                          'T2_US_Wisconsin',
#                          'T3_IT_Trieste',
#                          'T2_UK_SGrid_RALPP']
#config.Site.whitelist = ['T2_RU_JINR','T1_RU_JINR', 'T1_US_FNAL', 'T2_CH_CERN', 'T2_US_MIT']
#config.Site.whitelist = ['T2_UK_SGrid_RALPP', 'T2_CH_CERN ']  #FINISHED ON CRAB3FINAL
#config.Site.blacklist = ['T2_US_Wisconsin','T2_US_Purdue']
#config.Site.ignoreGlobalBlacklist = True


if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    from multiprocessing import Process
    
    def submit(config):
        res = crabCommand('submit', config = config)
        
    #########From now on that's what users should modify: this is the a-la-CRAB2 configuration part.

    requestNameList = ['BTagCSVRun_2016C-18Apr2017-v1',
                        'BTagCSVRun_2016D-18Apr2017-v1',
                        'BTagCSVRun_2016E-18Apr2017-v1',
                        'BTagCSVRun_2016F-18Apr2017-v1',
                        'BTagCSVRun_2016G-18Apr2017-v1',
                        'BTagCSVRun_2016H-18Apr2017-v1',
                        'BTagMuRun_2016C-18Apr2017-v1',
                        'BTagMuRun_2016D-18Apr2017-v1',
                        'BTagMuRun_2016E-18Apr2017-v1',
                        'BTagMuRun_2016F-18Apr2017-v1',
                        'BTagMuRun_2016G-18Apr2017-v1',
                        'BTagMuRun_2016H-18Apr2017-v1',
                        'CharmoniumRun_2016C-18Apr2017-v1',
                        'CharmoniumRun_2016D-18Apr2017-v1',
                        'CharmoniumRun_2016E-18Apr2017-v1',
                        'CharmoniumRun_2016F-18Apr2017-v1',
                        'CharmoniumRun_2016G-18Apr2017-v1',
                        'CharmoniumRun_2016H-18Apr2017-v1',
                        'DisplacedJetRun_2016C-18Apr2017-v1',
                        'DisplacedJetRun_2016D-18Apr2017-v1',
                        'DisplacedJetRun_2016E-18Apr2017-v1',
                        'DisplacedJetRun_2016F-18Apr2017-v1',
                        'DisplacedJetRun_2016G-18Apr2017-v1',
                        'DisplacedJetRun_2016H-18Apr2017-v1',
                        'DoubleEGRun_2016C-18Apr2017-v1',
                        'DoubleEGRun_2016D-18Apr2017-v1',
                        'DoubleEGRun_2016E-18Apr2017-v1',
                        'DoubleEGRun_2016F-18Apr2017-v1',
                        'DoubleEGRun_2016G-18Apr2017-v1',
                        'DoubleEGRun_2016H-18Apr2017-v1',
                        'DoubleMuonRun_2016C-18Apr2017-v1',
                        'DoubleMuonRun_2016D-18Apr2017-v1',
                        'DoubleMuonRun_2016E-18Apr2017-v1',
                        'DoubleMuonRun_2016F-18Apr2017-v1',
                        'DoubleMuonRun_2016G-18Apr2017-v1',
                        'DoubleMuonRun_2016H-18Apr2017-v1',
                        'DoubleMuonLowMassRun_2016C-18Apr2017-v1',
                        'DoubleMuonLowMassRun_2016D-18Apr2017-v1',
                        'DoubleMuonLowMassRun_2016E-18Apr2017-v1',
                        'DoubleMuonLowMassRun_2016F-18Apr2017-v1',
                        'DoubleMuonLowMassRun_2016G-18Apr2017-v1',
                        'DoubleMuonLowMassRun_2016H-18Apr2017-v1',
                        'FSQJetsRun_2016C-18Apr2017-v1',
                        'FSQJetsRun_2016G-18Apr2017-v1',
                        'FSQJetsRun_2016H-18Apr2017-v1',
                        'HTMHTRun_2016C-18Apr2017-v1',
                        'HTMHTRun_2016D-18Apr2017-v1',
                        'HTMHTRun_2016E-18Apr2017-v1',
                        'HTMHTRun_2016F-18Apr2017-v1',
                        'HTMHTRun_2016G-18Apr2017-v1',
                        'HTMHTRun_2016H-18Apr2017-v1',
                        'HighMultiplicityEOFRun_2016G-18Apr2017-v1',
                        'HighMultiplicityEOFRun_2016H-18Apr2017-v1',
                        'JetHTRun_2016C-18Apr2017-v1',
                        'JetHTRun_2016D-18Apr2017-v1',
                        'JetHTRun_2016E-18Apr2017-v1',
                        'JetHTRun_2016F-18Apr2017-v1',
                        'JetHTRun_2016G-18Apr2017-v1',
                        'JetHTRun_2016H-18Apr2017-v1',
                        'METRun_2016C-18Apr2017-v1',
                        'METRun_2016D-18Apr2017-v1',
                        'METRun_2016E-18Apr2017-v1',
                        'METRun_2016F-18Apr2017-v1',
                        'METRun_2016G-18Apr2017-v1',
                        'METRun_2016H-18Apr2017-v1',
                        'MinimumBiasRun_2016F-18Apr2017-v1',
                        'MuOniaRun_2016C-18Apr2017-v1',
                        'MuOniaRun_2016D-18Apr2017-v1',
                        'MuOniaRun_2016E-18Apr2017-v1',
                        'MuOniaRun_2016F-18Apr2017-v1',
                        'MuOniaRun_2016G-18Apr2017-v1',
                        'MuOniaRun_2016H-18Apr2017-v1',
                        'MuonEGRun_2016C-18Apr2017-v1',
                        'MuonEGRun_2016D-18Apr2017-v1',
                        'MuonEGRun_2016E-18Apr2017-v1',
                        'MuonEGRun_2016F-18Apr2017-v1',
                        'MuonEGRun_2016G-18Apr2017-v1',
                        'MuonEGRun_2016H-18Apr2017-v1',
                        'NoBPTXRun_2016C-18Apr2017-v1',
                        'NoBPTXRun_2016D-18Apr2017-v1',
                        'NoBPTXRun_2016E-18Apr2017-v1',
                        'NoBPTXRun_2016F-18Apr2017-v1',
                        'NoBPTXRun_2016G-18Apr2017-v1',
                        'NoBPTXRun_2016H-18Apr2017-v1',
                        'SingleElectronRun_2016C-18Apr2017-v1',
                        'SingleElectronRun_2016D-18Apr2017-v1',
                        'SingleElectronRun_2016E-18Apr2017-v1',
                        'SingleElectronRun_2016F-18Apr2017-v1',
                        'SingleElectronRun_2016G-18Apr2017-v1',
                        'SingleElectronRun_2016H-18Apr2017-v1',
                        'SingleMuonRun_2016C-18Apr2017-v1',
                        'SingleMuonRun_2016D-18Apr2017-v1',
                        'SingleMuonRun_2016E-18Apr2017-v1',
                        'SingleMuonRun_2016F-18Apr2017-v1',
                        'SingleMuonRun_2016G-18Apr2017-v1',
                        'SingleMuonRun_2016H-18Apr2017-v1',
                        'SinglePhotonRun_2016C-18Apr2017-v1',
                        'SinglePhotonRun_2016D-18Apr2017-v1',
                        'SinglePhotonRun_2016E-18Apr2017-v1',
                        'SinglePhotonRun_2016F-18Apr2017-v1',
                        'SinglePhotonRun_2016G-18Apr2017-v1',
                        'SinglePhotonRun_2016H-18Apr2017-v1',
                        'TauRun_2016C-18Apr2017-v1',
                        'TauRun_2016D-18Apr2017-v1',
                        'TauRun_2016E-18Apr2017-v1',
                        'TauRun_2016F-18Apr2017-v1',
                        'TauRun_2016G-18Apr2017-v1',
                        'TauRun_2016H-18Apr2017-v1',
                        'ZeroBiasRun_2016C-18Apr2017-v1',
                        'ZeroBiasRun_2016D-18Apr2017-v1',
                        'ZeroBiasRun_2016E-18Apr2017-v1',
                        'ZeroBiasRun_2016F-18Apr2017-v1',
                        'ZeroBiasRun_2016G-18Apr2017-v1',
                        'ZeroBiasRun_2016H-18Apr2017-v1']


    inputDatasetList = ['/BTagCSV/Run2016C-18Apr2017-v1/AOD',
                        '/BTagCSV/Run2016D-18Apr2017-v1/AOD',
                        '/BTagCSV/Run2016E-18Apr2017-v1/AOD',
                        '/BTagCSV/Run2016F-18Apr2017-v1/AOD',
                        '/BTagCSV/Run2016G-18Apr2017-v1/AOD',
                        '/BTagCSV/Run2016H-18Apr2017-v1/AOD',
                        '/BTagMu/Run2016C-18Apr2017-v1/AOD',
                        '/BTagMu/Run2016D-18Apr2017-v1/AOD',
                        '/BTagMu/Run2016E-18Apr2017-v1/AOD',
                        '/BTagMu/Run2016F-18Apr2017-v1/AOD',
                        '/BTagMu/Run2016G-18Apr2017-v1/AOD',
                        '/BTagMu/Run2016H-18Apr2017-v1/AOD',
                        '/Charmonium/Run2016C-18Apr2017-v1/AOD',
                        '/Charmonium/Run2016D-18Apr2017-v1/AOD',
                        '/Charmonium/Run2016E-18Apr2017-v1/AOD',
                        '/Charmonium/Run2016F-18Apr2017-v1/AOD',
                        '/Charmonium/Run2016G-18Apr2017-v1/AOD',
                        '/Charmonium/Run2016H-18Apr2017-v1/AOD',
                        '/DisplacedJet/Run2016C-18Apr2017-v1/AOD',
                        '/DisplacedJet/Run2016D-18Apr2017-v1/AOD',
                        '/DisplacedJet/Run2016E-18Apr2017-v1/AOD',
                        '/DisplacedJet/Run2016F-18Apr2017-v1/AOD',
                        '/DisplacedJet/Run2016G-18Apr2017-v1/AOD',
                        '/DisplacedJet/Run2016H-18Apr2017-v1/AOD',
                        '/DoubleEG/Run2016C-18Apr2017-v1/AOD',
                        '/DoubleEG/Run2016D-18Apr2017-v1/AOD',
                        '/DoubleEG/Run2016E-18Apr2017-v1/AOD',
                        '/DoubleEG/Run2016F-18Apr2017-v1/AOD',
                        '/DoubleEG/Run2016G-18Apr2017-v1/AOD',
                        '/DoubleEG/Run2016H-18Apr2017-v1/AOD',
                        '/DoubleMuon/Run2016C-18Apr2017-v1/AOD',
                        '/DoubleMuon/Run2016D-18Apr2017-v1/AOD',
                        '/DoubleMuon/Run2016E-18Apr2017-v1/AOD',
                        '/DoubleMuon/Run2016F-18Apr2017-v1/AOD',
                        '/DoubleMuon/Run2016G-18Apr2017-v1/AOD',
                        '/DoubleMuon/Run2016H-18Apr2017-v1/AOD',
                        '/DoubleMuonLowMass/Run2016C-18Apr2017-v1/AOD',
                        '/DoubleMuonLowMass/Run2016D-18Apr2017-v1/AOD',
                        '/DoubleMuonLowMass/Run2016E-18Apr2017-v1/AOD',
                        '/DoubleMuonLowMass/Run2016F-18Apr2017-v1/AOD',
                        '/DoubleMuonLowMass/Run2016G-18Apr2017-v1/AOD',
                        '/DoubleMuonLowMass/Run2016H-18Apr2017-v1/AOD',
                        '/FSQJets/Run2016C-18Apr2017-v1/AOD',
                        '/FSQJets/Run2016G-18Apr2017-v1/AOD',
                        '/FSQJets/Run2016H-18Apr2017-v1/AOD',
                        '/HTMHT/Run2016C-18Apr2017-v1/AOD',
                        '/HTMHT/Run2016D-18Apr2017-v1/AOD',
                        '/HTMHT/Run2016E-18Apr2017-v1/AOD',
                        '/HTMHT/Run2016F-18Apr2017-v1/AOD',
                        '/HTMHT/Run2016G-18Apr2017-v1/AOD',
                        '/HTMHT/Run2016H-18Apr2017-v1/AOD',
                        '/HighMultiplicityEOF/Run2016G-18Apr2017-v1/AOD',
                        '/HighMultiplicityEOF/Run2016H-18Apr2017-v1/AOD',
                        '/JetHT/Run2016C-18Apr2017-v1/AOD',
                        '/JetHT/Run2016D-18Apr2017-v1/AOD',
                        '/JetHT/Run2016E-18Apr2017-v1/AOD',
                        '/JetHT/Run2016F-18Apr2017-v1/AOD',
                        '/JetHT/Run2016G-18Apr2017-v1/AOD',
                        '/JetHT/Run2016H-18Apr2017-v1/AOD',
                        '/MET/Run2016C-18Apr2017-v1/AOD',
                        '/MET/Run2016D-18Apr2017-v1/AOD',
                        '/MET/Run2016E-18Apr2017-v1/AOD',
                        '/MET/Run2016F-18Apr2017-v1/AOD',
                        '/MET/Run2016G-18Apr2017-v1/AOD',
                        '/MET/Run2016H-18Apr2017-v1/AOD',
                        '/MinimumBias/Run2016F-18Apr2017-v1/AOD',
                        '/MuOnia/Run2016C-18Apr2017-v1/AOD',
                        '/MuOnia/Run2016D-18Apr2017-v1/AOD',
                        '/MuOnia/Run2016E-18Apr2017-v1/AOD',
                        '/MuOnia/Run2016F-18Apr2017-v1/AOD',
                        '/MuOnia/Run2016G-18Apr2017-v1/AOD',
                        '/MuOnia/Run2016H-18Apr2017-v1/AOD',
                        '/MuonEG/Run2016C-18Apr2017-v1/AOD',
                        '/MuonEG/Run2016D-18Apr2017-v1/AOD',
                        '/MuonEG/Run2016E-18Apr2017-v1/AOD',
                        '/MuonEG/Run2016F-18Apr2017-v1/AOD',
                        '/MuonEG/Run2016G-18Apr2017-v1/AOD',
                        '/MuonEG/Run2016H-18Apr2017-v1/AOD',
                        '/NoBPTX/Run2016C-18Apr2017-v1/AOD',
                        '/NoBPTX/Run2016D-18Apr2017-v1/AOD',
                        '/NoBPTX/Run2016E-18Apr2017-v1/AOD',
                        '/NoBPTX/Run2016F-18Apr2017-v1/AOD',
                        '/NoBPTX/Run2016G-18Apr2017-v1/AOD',
                        '/NoBPTX/Run2016H-18Apr2017-v1/AOD',
                        '/SingleElectron/Run2016C-18Apr2017-v1/AOD',
                        '/SingleElectron/Run2016D-18Apr2017-v1/AOD',
                        '/SingleElectron/Run2016E-18Apr2017-v1/AOD',
                        '/SingleElectron/Run2016F-18Apr2017-v1/AOD',
                        '/SingleElectron/Run2016G-18Apr2017-v1/AOD',
                        '/SingleElectron/Run2016H-18Apr2017-v1/AOD',
                        '/SingleMuon/Run2016C-18Apr2017-v1/AOD',
                        '/SingleMuon/Run2016D-18Apr2017-v1/AOD',
                        '/SingleMuon/Run2016E-18Apr2017-v1/AOD',
                        '/SingleMuon/Run2016F-18Apr2017-v1/AOD',
                        '/SingleMuon/Run2016G-18Apr2017-v1/AOD',
                        '/SingleMuon/Run2016H-18Apr2017-v1/AOD',
                        '/SinglePhoton/Run2016C-18Apr2017-v1/AOD',
                        '/SinglePhoton/Run2016D-18Apr2017-v1/AOD',
                        '/SinglePhoton/Run2016E-18Apr2017-v1/AOD',
                        '/SinglePhoton/Run2016F-18Apr2017-v1/AOD',
                        '/SinglePhoton/Run2016G-18Apr2017-v1/AOD',
                        '/SinglePhoton/Run2016H-18Apr2017-v1/AOD',
                        '/Tau/Run2016C-18Apr2017-v1/AOD',
                        '/Tau/Run2016D-18Apr2017-v1/AOD',
                        '/Tau/Run2016E-18Apr2017-v1/AOD',
                        '/Tau/Run2016F-18Apr2017-v1/AOD',
                        '/Tau/Run2016G-18Apr2017-v1/AOD',
                        '/Tau/Run2016H-18Apr2017-v1/AOD',
                        '/ZeroBias/Run2016C-18Apr2017-v1/AOD',
                        '/ZeroBias/Run2016D-18Apr2017-v1/AOD',
                        '/ZeroBias/Run2016E-18Apr2017-v1/AOD',
                        '/ZeroBias/Run2016F-18Apr2017-v1/AOD',
                        '/ZeroBias/Run2016G-18Apr2017-v1/AOD',
                        '/ZeroBias/Run2016H-18Apr2017-v1/AOD']



    for req,dataset in zip(requestNameList,inputDatasetList):
        config.General.requestName = req
        config.Data.inputDataset = dataset
        if "Run2016C" in dataset:
          config.Data.runRange = '275656-276283'  
        elif "Run2016D" in dataset:
          config.Data.runRange = '276315-276811'  
        elif "Run2016E" in dataset:
          config.Data.runRange = '276831-277420'
        elif "Run2016F" in dataset:
          config.Data.runRange = '277932-278808' 
        elif "Run2016G" in dataset:
          config.Data.runRange = '278820-280385' 
        elif "Run2016H" in dataset:
          config.Data.runRange = '281613-284044'           

        # print('REQUEST:', req, 'DATASET:',dataset)
        #config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
