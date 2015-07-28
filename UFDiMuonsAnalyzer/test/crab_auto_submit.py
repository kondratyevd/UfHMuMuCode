from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'dy_jetsToLL_PU40bx50'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dimuanalyzer_for_crab.py'
#config.JobType.outputFiles = ['outputfile.root']

config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased' # Must use this with JSON file
#config.Data.lumiMask = 's.jsonfiles[1]'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/acarnes/'
config.Data.publication = False
config.Data.publishDataName = 'dy_jetsToLL_PU40bx50'

config.Site.storageSite = 'T2_US_Florida'
