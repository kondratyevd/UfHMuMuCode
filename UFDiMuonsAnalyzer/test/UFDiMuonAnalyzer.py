import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

thisIsData = False

if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##process.MessageLogger.destinations.append("detailedInfo")
##process.MessageLogger.detailedInfo = cms.untracked.PSet(
##    threshold = cms.untracked.string("INFO"),
##    categories = cms.untracked.vstring("UFHLTTests")
##)

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

## Geometry and Detector Conditions (needed for a few patTuple production steps)

process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.EventContent.EventContent_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond


# global tag
globalTag = "POSTLS170_V5"
print 'Loading Global Tag: '+globalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag+"::All"


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2-Spring14miniaod-PU20bx25_POSTLS170_V5-v1.root'))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
# -------- PoolSource END -------------

#===============================================================================

# Clean the Jets from good muons, apply loose jet Id
ccMuPreSel = "pt > 15. && isGlobalMuon "
ccMuPreSel += " && globalTrack().normalizedChi2 < 10 "
ccMuPreSel += " && isPFMuon "
ccMuPreSel += " && innerTrack().hitPattern().trackerLayersWithMeasurement > 5 "
ccMuPreSel += " && innerTrack().hitPattern().numberOfValidPixelHits > 0 "
ccMuPreSel += " && globalTrack().hitPattern().numberOfValidMuonHits > 0 "
ccMuPreSel += " && numberOfMatchedStations > 1 && dB < 0.2 && abs(eta) < 2.4 "
ccMuPreSel += " && ( chargedHadronIso + max(0.,neutralHadronIso + photonIso - 0.5*puChargedHadronIso) ) < 0.12 * pt"

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
               deltaR              = cms.double(0.5),
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

#===============================================================================
# UFDiMuonAnalyzer

if thisIsData:
  process.load("UserArea.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_cff")
else:
  process.load("UserArea.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_MC_cff")

process.dimuons = process.DiMuons.clone()
process.dimuons.getFilename    = cms.untracked.string("DYJets_13TeV_Stage1.root")
process.dimuons.pfJetsTag = cms.InputTag("cleanJets")

#===============================================================================

process.p = cms.Path(#
                     process.cleanJets*
                     process.dimuons
                     )

#process.outpath = cms.EndPath()

## #Test to dump file content
## process.output = cms.OutputModule("PoolOutputModule",
##                                   outputCommands = cms.untracked.vstring("keep *"),
##                                   fileName = cms.untracked.string('dump.root')
##                                   )
## 
## process.out_step = cms.EndPath(process.output)

#===============================================================================

#process.source.fileNames.extend(
#[
##'file:/data/0b/digiovan/code/higgs/dev/addEle/CMSSW_5_3_3_patch3/src/UserArea/test/DYJetsToLL.root'
##"file:/data/uftrig01b/jhugon/hmumu/devNtupler/testFiles/VBFHToMM_M125_8TeV-powheg-pythia6-tauola-RECO_1.root"
##"file:/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/testPriVtxConstr/TTJetsSkims/TTJets_10_1_crI.root"
##"file:/home/jhugon/scratchRaid7/hmumu/recoData/VBFHToMM_M125_8TeV-powheg-pythia6-tauola-RECO_1.root"
#]
#)
#process.out.outputCommands = cms.untracked.vstring("keep *")
#process.outpath = cms.EndPath(process.out)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )


