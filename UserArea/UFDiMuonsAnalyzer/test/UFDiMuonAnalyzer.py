import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

thisIsData = False

if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration.EventContent.EventContent_cff')


# global tag
if thisIsData:
    print 'Loading Global Tag For Data'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "GR_R_53_V16::All"
else:
    print 'Loading Global Tag For MC'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "START53_V14::All"


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
# -------- PoolSource END -------------

#===============================================================================

## PF2PAT
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK5"

jetCorrections = ('AK5PF', ['L1FastJet','L2Relative','L3Absolute'])
if thisIsData:
  jetCorrections = ('AK5PF', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])

usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=(not thisIsData), postfix=postfix, jetCorrections=jetCorrections, typeIMetCorrections=True)

process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
#process.pfPileUpPFlow.Vertices = 'goodOfflinePrimaryVertices'
process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet = False
# Compute the mean pt per unit area (rho) from the
# PFchs inputs
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsPFlow = kt4PFJets.clone(
    rParam = cms.double(0.6),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )
process.patJetCorrFactorsPFlow.rho = cms.InputTag("kt6PFJetsPFlow", "rho")
# Add the PV selector and KT6 producer to the sequence
getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsPFlow )

# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True 
getattr(process,"pfNoMuon"+postfix).enable = True 
getattr(process,"pfNoElectron"+postfix).enable = True 
getattr(process,"pfNoTau"+postfix).enable = True 
getattr(process,"pfNoJet"+postfix).enable = True

# Clean the Jets from good muons, apply loose jet Id
ccMuPreSel = "pt > 20. && isGlobalMuon "
ccMuPreSel += " && globalTrack().normalizedChi2 < 10 "
ccMuPreSel += " && isPFMuon "
ccMuPreSel += " && innerTrack().hitPattern().trackerLayersWithMeasurement > 5 "
ccMuPreSel += " && innerTrack().hitPattern().numberOfValidPixelHits > 0 "
ccMuPreSel += " && globalTrack().hitPattern().numberOfValidMuonHits > 0 "
ccMuPreSel += " && numberOfMatchedStations > 1 && dB < 0.2 && abs(eta) < 2.4 "
ccMuPreSel += " && ( chargedHadronIso + neutralHadronIso + photonIso ) < 0.10 * pt"

jetSelection = 'neutralEmEnergy/energy < 0.99 '
jetSelection += ' && neutralHadronEnergy/energy < 0.99 '
jetSelection += ' && (chargedMultiplicity + neutralMultiplicity) > 1 '
jetSelection += ' && ((abs(eta)>2.4) || (chargedMultiplicity > 0 '
jetSelection += ' && chargedHadronEnergy/energy > 0.0'
jetSelection += ' && chargedEmEnergy/energy < 0.99))'

process.cleanPatJetsPFlow = cms.EDProducer("PATJetCleaner",
          src = cms.InputTag("selectedPatJetsPFlow"),
          preselection = cms.string(jetSelection),
          checkOverlaps = cms.PSet(
             muons = cms.PSet(
               src       = cms.InputTag("selectedPatMuonsPFlow"),
               algorithm = cms.string("byDeltaR"),
               preselection        = cms.string(ccMuPreSel),
               deltaR              = cms.double(0.5),
               checkRecoComponents = cms.bool(False),
               pairCut             = cms.string(""),
               requireNoOverlaps   = cms.bool(True),
             ),
         ),
         finalCut = cms.string('')
)

process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.jets = cms.InputTag("cleanPatJetsPFlow")
process.puJetMva.jets = cms.InputTag("cleanPatJetsPFlow")

#===============================================================================
# UFDiMuonAnalyzer
process.load("UserArea.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_nocuts_cff")
process.dimuons = process.DiMuons.clone()

if thisIsData:
  process.dimuons.isMonteCarlo   = cms.bool(False) 
else:
  process.dimuons.isMonteCarlo   = cms.bool(True) 
process.dimuons.checkTrigger   = cms.bool(False)
process.dimuons.processName    = cms.string("HLT")
process.dimuons.triggerNames   = cms.vstring("HLT_IsoMu24","HLT_Mu17_Mu8")
process.dimuons.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.dimuons.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")

process.dimuons.metTag         = cms.InputTag("patMETsPFlow")
process.dimuons.pfJetsTag      = cms.InputTag("cleanPatJetsPFlow")
process.dimuons.genJetsTag     = cms.InputTag("null")

process.dimuons.puJetMvaFullDiscTag = cms.InputTag("puJetMva","fullDiscriminant")
process.dimuons.puJetMvaFullIdTag = cms.InputTag("puJetMva","fullId")
process.dimuons.puJetMvaSimpleDiscTag = cms.InputTag("puJetMva","simpleDiscriminant")
process.dimuons.puJetMvaSimpleIdTag = cms.InputTag("puJetMva","simpleId")
process.dimuons.puJetMvaCutDiscTag = cms.InputTag("puJetMva","cutbasedDiscriminant")
process.dimuons.puJetMvaCutIdTag = cms.InputTag("puJetMva","cutbased")

process.dimuons.puJetMvaFullDiscTag = cms.InputTag("null")
process.dimuons.puJetMvaFullIdTag = cms.InputTag("null")
process.dimuons.puJetMvaSimpleDiscTag = cms.InputTag("null")
process.dimuons.puJetMvaSimpleIdTag = cms.InputTag("null")
process.dimuons.puJetMvaCutDiscTag = cms.InputTag("null")
process.dimuons.puJetMvaCutIdTag = cms.InputTag("null")


#===============================================================================

process.p = cms.Path( getattr(process,"patPF2PATSequence"+postfix)*
                      process.cleanPatJetsPFlow*
                      process.puJetId*
                      process.puJetMva*
                      process.dimuons
                    )

process.outpath = cms.EndPath()
#===============================================================================

process.dimuons.getFilename    = cms.untracked.string("vbfHmumu125.root")

process.source.fileNames.extend(
[
#"file:/data/uftrig01b/jhugon/hmumu/devNtupler/testFiles/VBFHToMM_M125_8TeV-powheg-pythia6-tauola-RECO_1.root"
#"file:/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/testPriVtxConstr/TTJetsSkims/TTJets_10_1_crI.root"
#"file:/home/jhugon/scratchRaid7/hmumu/recoData/VBFHToMM_M125_8TeV-powheg-pythia6-tauola-RECO_1.root"
]
)
#process.out.outputCommands = cms.untracked.vstring("keep *")
#process.outpath = cms.EndPath(process.out)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

