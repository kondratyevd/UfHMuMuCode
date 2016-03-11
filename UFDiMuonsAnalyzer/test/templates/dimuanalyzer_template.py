# =============================================================#
# UFDiMuonAnalyzer                                             #
# =============================================================#
# Makes stage1 trees.                                          #
# Adds a cleaner vector of jets to each event.                 #
# Originally Made by Justin Hugon. Edited by Andrew Carnes.    #
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
# Get a sample from our collection of samples
# /////////////////////////////////////////////////////////////

thisIsData = s.isData

print ""
print ""
if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'

print "Sample Name:    " +  s.name
print "Sample DAS DIR: " +  s.dir
if thisIsData:
    print "Sample JSON:    " +  s.jsonfiles[1]
print ""

# /////////////////////////////////////////////////////////////
# global tag, automatically retrieved from the imported sample
# /////////////////////////////////////////////////////////////

globalTag = s.globaltag

# The updated FrontierConditions_GlobalTag load needed for 2015 13TeV data does not like the ::All at the end of the tag
#if not thisIsData:
#    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#    globalTag+="::All"
#else:
#    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") 

# Newer MC global tags are consistent with data so we don't need the separate conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

print 'Loading Global Tag: ' + globalTag
process.GlobalTag.globaltag = globalTag

# /////////////////////////////////////////////////////////////
# ------------ PoolSource -------------
# /////////////////////////////////////////////////////////////
readFiles = cms.untracked.vstring();

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = readFiles)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()

# use a JSON file when locally executing cmsRun
if thisIsData:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = s.jsonfiles[1]).getVLuminosityBlockRange()


# /////////////////////////////////////////////////////////////
# -------- Add a cleaner vector jets to each event -----------
# /////////////////////////////////////////////////////////////

# Clean the Jets from good muons, apply loose jet Id
ccMuPreSel = "pt > 10. && isGlobalMuon "
ccMuPreSel += " && globalTrack().normalizedChi2 < 10 "
ccMuPreSel += " && isPFMuon "
ccMuPreSel += " && innerTrack().hitPattern().trackerLayersWithMeasurement > 5 "
ccMuPreSel += " && innerTrack().hitPattern().numberOfValidPixelHits > 0 "
ccMuPreSel += " && globalTrack().hitPattern().numberOfValidMuonHits > 0 "
ccMuPreSel += " && numberOfMatchedStations > 1 && dB < 0.2 && abs(eta) < 2.4 "
ccMuPreSel += " && ( chargedHadronIso + max(0.,neutralHadronIso + photonIso - 0.5*puChargedHadronIso) ) < 0.25 * pt"

jetSelection = 'neutralEmEnergy/energy < 0.99 '
jetSelection += ' && neutralHadronEnergy/energy < 0.99 '
jetSelection += ' && (chargedMultiplicity + neutralMultiplicity) > 1 '
jetSelection += ' && ((abs(eta)>2.4) || (chargedMultiplicity > 0 '
jetSelection += ' && chargedHadronEnergy/energy > 0.0'
jetSelection += ' && chargedEmEnergy/energy < 0.99))'

process.cleanJets = cms.EDProducer("PATJetCleaner",
          src = cms.InputTag("slimmedJets"),
          preselection = cms.string(jetSelection),
          checkOverlaps = cms.PSet(
             muons = cms.PSet(
               src       = cms.InputTag("slimmedMuons"),
               algorithm = cms.string("byDeltaR"),
               preselection        = cms.string(ccMuPreSel),
               deltaR              = cms.double(0.3),
               checkRecoComponents = cms.bool(False),
               pairCut             = cms.string(""),
               requireNoOverlaps   = cms.bool(True),
             ),
             #electrons = cms.PSet(
             #  src       = cms.InputTag("slimmedElectrons"),
             #  algorithm = cms.string("byDeltaR"),
             #  preselection        = cms.string(ccElePreSel),
             #  deltaR              = cms.double(0.5),
             #  checkRecoComponents = cms.bool(False),
             #  pairCut             = cms.string(""),
             #  requireNoOverlaps   = cms.bool(True),
             #),
         ),
         finalCut = cms.string('')
)

# /////////////////////////////////////////////////////////////
# Save output with TFileService
# /////////////////////////////////////////////////////////////

process.TFileService = cms.Service("TFileService", fileName = cms.string("stage_1_"+s.name+".root") )

# /////////////////////////////////////////////////////////////
# Load UFDiMuonAnalyzer
# /////////////////////////////////////////////////////////////

if thisIsData:
  process.load("UfHMuMuCode.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_cff")
else:
  process.load("UfHMuMuCode.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_MC_cff")

process.dimuons = process.DiMuons.clone()
process.dimuons.pfJetsTag = cms.InputTag("cleanJets")


# /////////////////////////////////////////////////////////////
# Set the order of operations
# /////////////////////////////////////////////////////////////

process.p = cms.Path(
                     process.cleanJets*
                     process.dimuons
                     )
