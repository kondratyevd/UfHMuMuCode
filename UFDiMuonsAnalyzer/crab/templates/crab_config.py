
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.requestName = 'STR'
config.General.workArea = 'logs'
config.General.transferOutputs = True  ## Do output root files
config.General.transferLogs = True

config.section_('JobType')
config.JobType.psetName = 'STR'
config.JobType.outputFiles = ['tuple.root'] ## Must be the same as the output file in process.TFileService in config.JobType.psetName python file
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDBS = 'global'

config.Data.inputDataset = 'STR'
# config.Data.lumiMask = 'STR'

config.Data.useParent = False
config.Data.splitting = 'STR'
config.Data.unitsPerJob = NUM  ## Should take ~10 minutes for 100k events

config.Data.outLFNDirBase = '/store/user/dkondrat'
config.Data.publication = False
config.Data.outputDatasetTag = 'STR'
config.Data.ignoreLocality = True 

config.section_('User')

config.section_('Site')
config.Site.blacklist = ['T2_US_*']
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T2_US_Purdue'
