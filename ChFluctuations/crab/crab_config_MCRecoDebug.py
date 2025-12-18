###For a description of the crabConfig.py parameters. See:
###https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.workArea        = 'Test_Run3_NCF_MC'
config.General.requestName     = 'Test_Run3_NCF_MC'
config.General.transferLogs    = False 
config.General.transferOutputs = True
################################
config.section_('JobType')
config.JobType.pluginName      = 'Analysis'
config.JobType.psetName        = '../cfg/subash_config_cfg.py'
config.JobType.inputFiles      = ['../cfg/CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_Nominal.db']
################################
config.section_('Data')
#config.Data.inputDataset       = '/MinBias_PbPb_5p36TeV_Hydjet_v1/sarteaga-MinBias_PbPb_5p36TeV_Hydjet_RECODEBUG_v5-0e6c11377ba727d4466887a72ad361ed/USER'
#config.Data.inputDataset       = '/MinBias_PbPb_5p36TeV_Hydjet_v1/sarteaga-MinBias_PbPb_5p36TeV_Hydjet_MINIAOD_v5-45d226c4498e665b330b684cce3e6ab4/USER'
config.Data.inputDataset       = '/MinBias_Drum5F_5p36TeV_hydjet/HINPbPbSpring23MiniAOD-NoPU_132X_mcRun3_2023_realistic_HI_v9-v2/MINIAODSIM'
config.Data.splitting          = 'FileBased'
config.Data.unitsPerJob        = 10
config.Data.totalUnits         = -1
config.Data.publication        = False
config.Data.inputDBS           = 'global'
config.Data.outLFNDirBase      = '/store/user/anarzary/pbpb'
config.Data.outputDatasetTag   = 'Test_Run3_NCF_MC'
################################
config.section_('Site')
config.Site.storageSite        = 'T3_CH_CERNBOX'


