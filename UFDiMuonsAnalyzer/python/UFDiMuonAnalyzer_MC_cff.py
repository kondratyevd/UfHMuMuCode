import FWCore.ParameterSet.Config as cms

## Assumes you are looking at miniAOD

DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         muonColl = cms.InputTag("slimmedMuons"),
                         
                         isVerbose = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(True),
                         nMuons = cms.int32(2),
                         
                         # muon kinematic cuts
                         isGlobal    =cms.int32(0),
                         isStandAlone=cms.int32(0),
                         isTracker   =cms.int32(0),
                         ptMin       = cms.double(20),
                         etaMax      = cms.double(2.4),
                         normChiSquareMax= cms.double(999),

                         # number of hits cuts
                         minMuonHits	  = cms.int32(1),
                         minPixelHits	  = cms.int32(1),
                         minStripHits	  = cms.int32(-999),
                         minTrackerHits	  = cms.int32(-999),
                         minSegmentMatches= cms.int32(2),
                         minNumOfMatchedStations = cms.int32(-999),

                         minPixelLayers  = cms.int32(-999),
                         minTrackerLayers= cms.int32(8),
                         minStripLayers  = cms.int32(-999),
                         validFracTrackerMin= cms.int32(-999),   
    
                         primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         beamSpotTag = cms.InputTag("offlineBeamSpot"),
                         prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                         packedGenParticleTag = cms.InputTag("packedGenParticles"),
                         d0Max = cms.double(0.2), 
                         
                         #track isolation
                         trackIsoMaxSumPt = cms.double(9999),
                         relCombIsoMax    = cms.double(0.1),
                         
                         #triggerName = cms.string("@"),#wild card: all triggers
                         # HLT trigger info
                         checkTrigger = cms.bool(False),

                         processName = cms.string("HLT"),
                         triggerNames = cms.vstring("HLT_IsoMu24_eta2p1"), # up to 3 trigger names
                         triggerResults = cms.InputTag("TriggerResults","","HLT"),
                         triggerObjs = cms.InputTag("selectedPatTrigger"),
                         
                         metTag         = cms.InputTag("slimmedMETs"),
                         pfJetsTag      = cms.InputTag("slimmedJets"),
                         genJetsTag     = cms.InputTag("slimmedGenJets")
                         
                         )
