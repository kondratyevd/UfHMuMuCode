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
                         normChiSquareMax= cms.double(999),

                         # number of hits cuts
                         minMuonHits	  = cms.int32(-999),
                         minPixelHits	  = cms.int32(-999),
                         minStripHits	  = cms.int32(-999),
                         minTrackerHits	  = cms.int32(-999),
                         minSegmentMatches= cms.int32(-999),
                         minNumOfMatchedStations = cms.int32(-999),

                         minPixelLayers  = cms.int32(-999),
                         minTrackerLayers= cms.int32(-999),
                         minStripLayers  = cms.int32(-999),
                         validFracTrackerMin= cms.int32(-999),   
    
                         primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         beamSpotTag = cms.InputTag("offlineBeamSpot"),
                         prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                         packedGenParticleTag = cms.InputTag("packedGenParticles"),
                         d0Max = cms.double(999), 
                         
                         #track isolation
                         trackIsoMaxSumPt = cms.double(999),
                         relCombIsoMax    = cms.double(999),
                         
                         #triggerName = cms.string("@"),#wild card: all triggers
                         # HLT trigger info
                         checkTrigger = cms.bool(True),

                         processName = cms.string("HLT"),
                         triggerNames = cms.vstring("HLT_IsoMu20","HLT_IsoTkMu20","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"), # up to 3 trigger names
                         triggerResults = cms.InputTag("TriggerResults","","HLT"),
                         triggerObjs = cms.InputTag("selectedPatTrigger"),
                         
                         metTag         = cms.InputTag("slimmedMETs"),
                         pfJetsTag      = cms.InputTag("slimmedJets"),
                         genJetsTag     = cms.InputTag("slimmedGenJets")
                         
                         )
