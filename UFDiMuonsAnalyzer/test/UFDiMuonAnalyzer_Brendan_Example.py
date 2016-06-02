# =============================================================#
# UFDiMuonAnalyzer                                             #
# =============================================================#
# Makes stage1 trees.                                          #
#                                                              #
################################################################ 

# /////////////////////////////////////////////////////////////
# Load some things
# /////////////////////////////////////////////////////////////

import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##process.MessageLogger.destinations.append("detailedInfo")
##process.MessageLogger.detailedInfo = cms.untracked.PSet(
##    threshold = cms.untracked.string("INFO"),
##    categories = cms.untracked.vstring("UFHLTTests")
##)


## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond

# /////////////////////////////////////////////////////////////
# INPUT INFO ABOUT YOUR SAMPLE
# /////////////////////////////////////////////////////////////

isData = False
sname = "Name of mc sample. Used to name the output files."
sdir =  "NO DAS DIR. Locally created samples."               ### You don't need to worry about this.
sjsonfile = "/no/json/for/mc"
sglobaltag = "YOUR_GLOBAL_TAG"                               ### Put the correct global tag here
sfilename = "file:/LOCATION/OF/SAMPLE/SAMPLENAME.root"       ### The MINIAOD root file you want to run over

# /////////////////////////////////////////////////////////////
# Set up and run the analyzer
# /////////////////////////////////////////////////////////////

thisIsData = isData

print ""
print ""
if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'

print "Sample Name:    " +  sname
print "Sample DAS DIR: " +  sdir
if thisIsData:
    print "Sample JSON:    " +  sjsonfile
print ""

# /////////////////////////////////////////////////////////////
# global tag
# /////////////////////////////////////////////////////////////

globalTag = sglobaltag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

print 'Loading Global Tag: ' + globalTag
process.GlobalTag.globaltag = globalTag

# /////////////////////////////////////////////////////////////
# ------------ PoolSource -------------
# /////////////////////////////////////////////////////////////
readFiles = cms.untracked.vstring();
readFiles.append(sfilename);

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = readFiles)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()

# use a JSON file when locally executing cmsRun
if thisIsData:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = sjsonfile).getVLuminosityBlockRange()



# /////////////////////////////////////////////////////////////
# Save output with TFileService
# /////////////////////////////////////////////////////////////

process.TFileService = cms.Service("TFileService", fileName = cms.string("stage_1_"+sname+".root") )

# /////////////////////////////////////////////////////////////
# Load UFDiMuonAnalyzer
# /////////////////////////////////////////////////////////////

if thisIsData:
  process.load("UfHMuMuCode.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_cff")
else:
  process.load("UfHMuMuCode.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_MC_cff")

process.dimuons = process.DiMuons.clone()

# /////////////////////////////////////////////////////////////
# Set the order of operations
# /////////////////////////////////////////////////////////////

process.p = cms.Path(process.dimuons)
