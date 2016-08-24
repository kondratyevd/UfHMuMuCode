import FWCore.ParameterSet.Config as cms

## Assumes you are looking at miniAOD

DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         muonColl = cms.InputTag("slimmedMuons"),
                         
                         isVerbose = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(False),
                         nMuons = cms.int32(2),
                         
                         # muon kinematic cuts
                         isGlobal    =cms.int32(0),
                         isStandAlone=cms.int32(0),
                         isTracker   =cms.int32(1),
                         ptMin       = cms.double(10),
                         etaMax      = cms.double(2.4),
    
                         primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         beamSpotTag = cms.InputTag("offlineBeamSpot"),
                         prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                         packedGenParticleTag = cms.InputTag("packedGenParticles"),
                         
                         #triggerName = cms.string("@"),#wild card: all triggers
                         # HLT trigger info
                         checkTrigger = cms.bool(True),

                         processName = cms.string("HLT"),
                         triggerNames = cms.vstring("HLT_IsoMu20", "HLT_IsoTkMu20","HLT_IsoMu22", "HLT_IsoTkMu22","HLT_IsoMu24","HLT_IsoTkMu24"), # up to 6 trigger names
                         triggerResults = cms.InputTag("TriggerResults","","HLT"),
                         triggerObjs = cms.InputTag("selectedPatTrigger"),
                         
                         metTag         = cms.InputTag("slimmedMETs"),
                         pfJetsTag      = cms.InputTag("slimmedJets"),
                         genJetsTag     = cms.InputTag("slimmedGenJets")
                         )
