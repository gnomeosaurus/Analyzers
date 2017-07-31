from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'crab3thirdsubmission'
config.General.requestName = 'Run2016H387runs'

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/user/f/fsiroky/CMSSW_8_0_28/src/Analyzers/AODAnalyzer/test/AODAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['AODTree.root']
config.JobType.maxJobRuntimeMin = 2750 #45 h


config.section_('Data')
# config.Data.inputDataset = '/DoubleEG/Run2016H-PromptReco-v1/AOD'
config.Data.inputDataset = '/DoubleEG/Run2016H-18Apr2017-v1/AOD'
config.Data.unitsPerJob = 40 #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.runRange = '281613-282000'
# config.Data.runRange = '281613-284044' //2016H range
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
#config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
config.Data.outLFNDirBase = '/store/group/alca_global/filipTemporary/'
config.Data.ignoreLocality = False

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.Site.whitelist = ['T3_IT_Trieste','T2_UK_SGrid_RALPP','T2_RU_JINR','T1_RU_JINR']
config.Site.blacklist = ['T2_US_Wisconsin','T2_US_Purdue']



# if __name__ == '__main__':
#     from CRABAPI.RawCommand import crabCommand
#     from multiprocessing import Process
    
#     def submit(config):
#         res = crabCommand('submit', config = config)
        
#     #########From now on that's what users should modify: this is the a-la-CRAB2 configuration part.

#     requestNameList = ['DoubleEG_2016H_v1',]
#                        #'DoubleEG_2016H_v2',
#                        #'DoubleEG_2016H_v3']

#     inputDatasetList = ['/DoubleEG/Run2016H-PromptReco-v1/AOD',]
#                         #'/DoubleEG/Run2016H-PromptReco-v2/AOD',
#                         #'/DoubleEG/Run2016H-PromptReco-v3/AOD']

#     for req,dataset in zip(requestNameList,inputDatasetList):
#         config.General.requestName = req
#         config.Data.inputDataset = dataset
#         print 'REQUEST:', req, 'DATASET:',dataset
#         #config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
#         p = Process(target=submit, args=(config,))
#         p.start()
#         p.join()
        
