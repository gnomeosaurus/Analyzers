from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab3Submission'
config.General.requestName = 'JetHT_2016G'

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/d/deguio/DQM/ML-DQM/CMSSW_8_0_20/src/Analyzers/miniAODAnalyzer/test/miniAODAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['miniAODTree.root']

config.section_('Data')
config.Data.inputDataset = '/JetHT/Run2016G-PromptReco-v1/MINIAOD'
config.Data.unitsPerJob = 100 #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
#config.Data.runRange = ''
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
config.Data.ignoreLocality = True

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_CH_CERN']
