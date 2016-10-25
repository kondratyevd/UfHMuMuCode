import FWCore.ParameterSet.Config as cms

## Assumes you are looking at miniAOD

DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         muonColl = cms.InputTag("slimmedMuons"),
                         
                         isVerbose    = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(True),
                         
                         # muon cuts for dimuon pair
                         nMuons = cms.int32(2),        # require minimum number of muons in the event to pass the following dimu pair criteria
                         isGlobal    =cms.int32(0),
                         isStandAlone=cms.int32(0),
                         isTracker   =cms.int32(1),
                         ptMin       = cms.double(10),
                         etaMax      = cms.double(2.4),
    
                         # HLT trigger info
                         checkTrigger = cms.bool(False),    # For MC we don't require that at least one HLT trigger in the list fires
                         processName  = cms.string("HLT"),
                         triggerNames = cms.vstring("HLT_IsoMu20", "HLT_IsoTkMu20","HLT_IsoMu22", "HLT_IsoTkMu22","HLT_IsoMu24","HLT_IsoTkMu24"), # up to 6 trigger names

                         triggerResults = cms.InputTag("TriggerResults","","HLT"),
                         triggerObjs    = cms.InputTag("selectedPatTrigger"),
                         
                         # Vertex and Beam Spot
                         primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         beamSpotTag      = cms.InputTag("offlineBeamSpot"),

                         # Gen particle collections besides jets
                         prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                         packedGenParticleTag = cms.InputTag("packedGenParticles"),
                         
                         # MET
                         metTag         = cms.InputTag("slimmedMETs"),
 
                         # Jets
                         pfJetsTag      = cms.InputTag("slimmedJets"),
                         genJetsTag     = cms.InputTag("slimmedGenJets"),
                         btagNames      = cms.untracked.vstring(["pfCombinedInclusiveSecondaryVertexV2BJetTags"]),

                         # Electrons
                         electronColl              = cms.InputTag("slimmedElectrons"),
                         electronCutBasedId_veto   = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                         electronCutBasedId_loose  = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                         electronCutBasedId_medium = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                         electronCutBasedId_tight  = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),

                         # Taus
                         tauColl    = cms.InputTag("slimmedTaus"),
                         tauIDNames = cms.untracked.vstring(["byCombinedIsolationDeltaBetaCorrRaw3Hits", "againstElectronLooseMVA6", "againstElectronVLooseMVA6",
                                                             "decayModeFindingNewDMs", "byLooseCombinedIsolationDeltaBetaCorr3Hits", 
                                                             "byMediumCombinedIsolationDeltaBetaCorr3Hits", "byTightCombinedIsolationDeltaBetaCorr3Hits", "againstMuonLoose3", 
                                                             "againstMuonTight3", "againstElectronVTightMVA6", "againstElectronTightMVA6", "againstElectronMediumMVA6"]),
                         )
