import CRABClient
from CRABClient.UserUtilities import config

config = config()
# config.General.requestName = f'{run_era}_{run_number}_07Mar2023_BKG_analysis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step3_RAW2DIGI_RECO.py'
# config.JobType.psetName = 'step3_flower.py'

# config.Data.inputDataset = f'/ZeroBias/{run_era}-PromptReco-v1/AOD'
# config.Data.inputDataset = f'/ZeroBias/{run_era}-v1/RAW'
config.Data.inputDBS = 'global'
config.Site.storageSite = 'T3_KR_UOS'
config.Site.outLFNDirBase = '/store/user/jheo/'
# config.Data.splitting = 'LumiBased'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
# config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt'
# config.Data.runRange = '359691-359694'
# config.Data.runRange = f'{run_number}'
config.Data.publication = False
# config.Data.outputDatasetTag = f'{fill}_{run_number}'
