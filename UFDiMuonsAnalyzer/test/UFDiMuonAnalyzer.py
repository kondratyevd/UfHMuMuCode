import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonAnalyzer")

thisIsData = True
thisIs2011 = False

if thisIsData:
    print 'Running over data sample'
else:
    print 'Running over MC sample'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

## Geometry and Detector Conditions (needed for a few patTuple production steps)

if thisIs2011:
    print "Analyzing a 7 TeV dataset"
    process.load("Configuration.StandardSequences.GeometryExtended_cff")
else:
    print "Analyzing a 8 TeV dataset"
    process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.EventContent.EventContent_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond


# global tag
print 'Loading Global Tag: FT_53_V21_AN4'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "FT_53_V21_AN4::All"


# ------------ PoolSource -------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:SingleMuRun2012A-22Jan2013-v1_1_3_W7l.root'))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
# -------- PoolSource END -------------

#===============================================================================

## Test JEC from test instances of the global DB
#process.load("PhysicsTools.PatAlgos.patTestJEC_cfi")

## Test JEC from local sqlite file
#process.load("PhysicsTools.PatAlgos.patTestJEC_local_cfi")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
			       dropMetaData = cms.untracked.string('ALL'),
	                       splitLevel = cms.untracked.int32(99),
                               SelectEvents = cms.untracked.PSet(      SelectEvents = cms.vstring('p')     ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *' )
                               )

process.outpath = cms.EndPath(process.out)


# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")


# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

# An empty postfix means that only PF2PAT is run,
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
# collections have standard names + postfix (e.g. patElectronPFlow)
#### Standard PAT Configuration File
##process.load("PhysicsTools.PatAlgos.patSequences_cff")
##
#### PF2PAT
###from PhysicsTools.PatAlgos.patTemplate_cfg import *
###from PhysicsTools.PatAlgos.tools.coreTools import *
##from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK5"

# this is to fix a problem with the standard PAT which was using L1Offset
# which is not anymore available for gt FT_53_V21_AN4
process.patJetCorrFactors.levels =  ['L1FastJet','L2Relative','L3Absolute']
process.patJetCorrFactors.useRho = True

jetCorrections = ('AK5PF', ['L1FastJet','L2Relative','L3Absolute'])
if thisIsData:
  jetCorrections = ('AK5PF', ['L1FastJet','L2Relative','L3Absolute','L2L3Residual'])

usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=(not thisIsData), postfix=postfix, jetCorrections=jetCorrections, typeIMetCorrections=True)
usePFIso( process ) 

if thisIsData:
    # removing MC matching for standard PAT sequence
    # for the PF2PAT+PAT sequence, it is done in the usePF2PAT function
    removeMCMatchingPF2PAT( process, '' )

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

# Taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/VHbbAnalysis/HbbAnalyzer/test/patMC.py?revision=1.9&view=markup
ccElePreSel = "pt > 15.0 && abs(eta) < 2.5 &&"
ccElePreSel += "(isEE || isEB) && !isEBEEGap &&"
ccElePreSel += " (chargedHadronIso + neutralHadronIso + photonIso)/pt <0.10 &&"
ccElePreSel += "dB < 0.02 && "  #dB is computed wrt PV but is transverse only, no info about dZ(vertex) 
ccElePreSel += "( "
ccElePreSel += "(isEE && ("
ccElePreSel += "abs(deltaEtaSuperClusterTrackAtVtx) < 0.005 &&  abs(deltaPhiSuperClusterTrackAtVtx) < 0.02 && sigmaIetaIeta < 0.03 && hadronicOverEm < 0.10 &&  abs(1./ecalEnergy*(1.-eSuperClusterOverP)) < 0.05 "
ccElePreSel += ")) || " 
ccElePreSel += "(isEB && (  "
ccElePreSel += "abs(deltaEtaSuperClusterTrackAtVtx) < 0.004 &&  abs(deltaPhiSuperClusterTrackAtVtx) < 0.03 && sigmaIetaIeta < 0.01 && hadronicOverEm < 0.12 && abs(1./ecalEnergy*(1.-eSuperClusterOverP)) < 0.05"
ccElePreSel += "))"
ccElePreSel += ")" 

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
             electrons = cms.PSet(
               src       = cms.InputTag("selectedPatElectronsPFlow"),
               algorithm = cms.string("byDeltaR"),
               preselection        = cms.string(ccElePreSel),
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
process.dimuons.triggerNames   = cms.vstring("HLT_IsoMu24_eta2p1","HLT_Mu40_eta2p1")
process.dimuons.triggerResults = cms.InputTag("TriggerResults","","HLT")
process.dimuons.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")

process.dimuons.eleColl = cms.InputTag("selectedPatElectronsPFlow")
process.dimuons.eleTriggerName = cms.string("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL");

process.dimuons.metTag         = cms.InputTag("patMETsPFlow")
process.dimuons.pfJetsTag      = cms.InputTag("cleanPatJetsPFlow")
process.dimuons.genJetsTag     = cms.InputTag("null")

process.dimuons.puJetMvaFullDiscTag = cms.InputTag("puJetMva","full53xDiscriminant")
process.dimuons.puJetMvaFullIdTag = cms.InputTag("puJetMva","full53xId")
process.dimuons.puJetMvaSimpleDiscTag = cms.InputTag("puJetMva","simpleDiscriminant")
process.dimuons.puJetMvaSimpleIdTag = cms.InputTag("puJetMva","simpleId")
process.dimuons.puJetMvaCutDiscTag = cms.InputTag("puJetMva","cutbasedDiscriminant")
process.dimuons.puJetMvaCutIdTag = cms.InputTag("puJetMva","cutbased")


########## Hbb specific stuff starts here ########################
# rho2.5 calculation
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets25 = process.kt4PFJets.clone( src = 'pfNoElectron'+postfix,rParam = 0.6,doRhoFastjet = True,Ghost_EtaMax = 2.5, Rho_EtaMax = 2.5 )

#Below: special case only for HZZ-4l in order to be synchronized with the muon rho computation in the same analysis:

process.kt6PFJets = process.kt4PFJets.clone(src = 'pfNoElectron'+postfix,
                                            rParam = cms.double(0.6),
                                            doAreaFastjet = cms.bool(True),
                                            doRhoFastjet = cms.bool(True)
                                            )

process.kt6PFJets25asHtoZZto4l = process.kt6PFJets.clone(rParam = cms.double(0.6),
                                                         Rho_EtaMax = cms.double(2.5),
                                                         Ghost_EtaMax = cms.double(2.5),
                                                         )

# end of rho2.5 calculation

#===============================================================================
# Electron Selection

#Electron ID
#process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )

# try
process.patElectronsPFlow.electronIDSources = cms.PSet(
#process.patElectrons.electronIDSources = cms.PSet(
    #MVA
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
)   

process.patElectronsPFlow.isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso")
    )
process.patElectronsPFlow.isolationValuesNoPFId = cms.PSet(
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso")
    )

#===============================================================================

process.p = cms.Path(#
                     process.mvaID*                  
                     process.patDefaultSequence*
                     getattr(process,"patPF2PATSequence"+postfix)*
                     #process.patPF2PATSequencePFlow *
                     process.cleanPatJetsPFlow*
                     process.puJetIdSqeuence*#process.puJetId*
                     #process.puJetMva*
                     process.kt6PFJets25*
                     process.kt6PFJets25asHtoZZto4l*
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

process.dimuons.getFilename    = cms.untracked.string("SingleMuRun2012A-22Jan2013-v1.root")

process.source.fileNames.extend(
[
#'file:/data/0b/digiovan/code/higgs/dev/addEle/CMSSW_5_3_3_patch3/src/UserArea/test/DYJetsToLL.root'
#"file:/data/uftrig01b/jhugon/hmumu/devNtupler/testFiles/VBFHToMM_M125_8TeV-powheg-pythia6-tauola-RECO_1.root"
#"file:/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/testPriVtxConstr/TTJetsSkims/TTJets_10_1_crI.root"
#"file:/home/jhugon/scratchRaid7/hmumu/recoData/VBFHToMM_M125_8TeV-powheg-pythia6-tauola-RECO_1.root"
]
)
#process.out.outputCommands = cms.untracked.vstring("keep *")
#process.outpath = cms.EndPath(process.out)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )


