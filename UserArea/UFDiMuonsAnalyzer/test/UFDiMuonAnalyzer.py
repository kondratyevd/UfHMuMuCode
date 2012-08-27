import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

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
    process.GlobalTag.globaltag = "GR_R_52_V9::All"
else:
    print 'Loading Global Tag For MC'
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = "START53_V10::All"


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
# -------- PoolSource END -------------


#===============================================================================
# remove all the cuts
process.load("UserArea.UFDiMuonsAnalyzer.UFDiMuonAnalyzer_nocuts_cff")
process.dimuons = process.DiMuons.clone()

process.dimuons.getFilename    = cms.untracked.string("vbfHmumu150.root")
if thisIsData:
  process.dimuons.isMonteCarlo   = cms.bool(False) 
else:
  process.dimuons.isMonteCarlo   = cms.bool(True) 
process.dimuons.checkTrigger   = cms.bool(False)
process.dimuons.processName    = cms.string("HLT")
process.dimuons.triggerNames   = cms.vstring("HLT_Mu40_eta2p1","HLT_Mu22_TkMu22")
process.dimuons.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.dimuons.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")

process.dimuons.metTag         = cms.InputTag("pfMet")
process.dimuons.pfJetsTag      = cms.InputTag("cleanPatJets")
process.dimuons.genJetsTag     = cms.InputTag("null")

process.dimuons.puJetMvaFullDiscTag = cms.InputTag("puJetMva","fullDiscriminant")
process.dimuons.puJetMvaFullIdTag = cms.InputTag("puJetMva","fullId")
process.dimuons.puJetMvaSimpleDiscTag = cms.InputTag("puJetMva","simpleDiscriminant")
process.dimuons.puJetMvaSimpleIdTag = cms.InputTag("puJetMva","simpleId")
process.dimuons.puJetMvaCutDiscTag = cms.InputTag("puJetMva","cutbasedDiscriminant")
process.dimuons.puJetMvaCutIdTag = cms.InputTag("puJetMva","cutbased")


## ADDING PAT
removeMCMatching(process, ['All'])

if thisIsData:
    print "\nData Jet Corrections"
    switchJetCollection(process,cms.InputTag('ak5PFJets'),
                        doJTA        = True,
                        doBTagging   = True,
                        jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual'])),
                        doType1MET   = True,
                        #   genJetCollection=cms.InputTag("ak5GenJets"),
                        doJetID      = True
                        )

else:
    print "\nMC Jet Corrections"
    switchJetCollection(process,cms.InputTag('ak5PFJets'),
                        doJTA        = True,
                        doBTagging   = True,
                        jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute'])),
                        doType1MET   = True,
                        genJetCollection=cms.InputTag("ak5GenJets"),
                        doJetID      = True
                        )

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
process.cleanPatJets = cms.EDProducer("PATJetCleaner",
          src = cms.InputTag("patJets"),
          preselection = cms.string(jetSelection),
          checkOverlaps = cms.PSet(
             muons = cms.PSet(
               src       = cms.InputTag("selectedPatMuons"),
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
process.puJetId.jets = cms.InputTag("cleanPatJets")
process.puJetMva.jets = cms.InputTag("cleanPatJets")

#===============================================================================

process.p = cms.Path( process.patDefaultSequence*
                      process.cleanPatJets*
                      process.puJetId*
                      process.puJetMva*
                      process.dimuons
                    )

process.outpath = cms.EndPath()
#===============================================================================

process.source.fileNames.extend(
[
]
)

