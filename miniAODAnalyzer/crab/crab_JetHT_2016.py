from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = False
config.General.workArea = 'crab3Submission_0310'
config.General.requestName = ''

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/d/deguio/DQM/ML-DQM/CMSSW_8_0_20/src/Analyzers/miniAODAnalyzer/test/miniAODAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['miniAODTree.root']
config.JobType.maxJobRuntimeMin = 2750 #45 h


config.section_('Data')
config.Data.inputDataset = ''
config.Data.unitsPerJob = 100 #without '' since it must be an int
config.Data.splitting = 'LumiBased'
config.Data.publication = False
#config.Data.runRange = ''
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
#config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/deguio/ML-DQM/'
config.Data.ignoreLocality = False

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.whitelist = ['T2_CH_CERN']



if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    from multiprocessing import Process
    
    def submit(config):
        res = crabCommand('submit', config = config)
        
    #########From now on that's what users should modify: this is the a-la-CRAB2 configuration part.

    requestNameList = ['JetHT_2016B_v1',
                       'JetHT_2016B_v2',
                       'JetHT_2016C_v2',
                       'JetHT_2016D_v2',
                       'JetHT_2016E_v2',
                       'JetHT_2016F_v1',
                       'JetHT_2016G_v1',
                       'JetHT_2016H_v1',
                       'JetHT_2016H_v2']

    inputDatasetList = ['/JetHT/Run2016B-PromptReco-v1/MINIAOD',
                        '/JetHT/Run2016B-PromptReco-v2/MINIAOD',
                        '/JetHT/Run2016C-PromptReco-v2/MINIAOD',
                        '/JetHT/Run2016D-PromptReco-v2/MINIAOD',
                        '/JetHT/Run2016E-PromptReco-v2/MINIAOD',
                        '/JetHT/Run2016F-PromptReco-v1/MINIAOD',
                        '/JetHT/Run2016G-PromptReco-v1/MINIAOD',
                        '/JetHT/Run2016H-PromptReco-v1/MINIAOD',
                        '/JetHT/Run2016H-PromptReco-v2/MINIAOD']

    for req,dataset in zip(requestNameList,inputDatasetList):
        config.General.requestName = req
        config.Data.inputDataset = dataset
        print 'REQUEST:', req, 'DATASET:',dataset
        #config.Data.outLFNDirBase = '/store/user/deguio/ML-DQM/'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
