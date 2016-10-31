import FWCore.ParameterSet.Config as cms

## Assumes you are looking at miniAOD

DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         muonColl = cms.InputTag("slimmedMuons"),
                         
                         isVerbose    = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(False),
                         
                         # muon cuts for dimuon pair
                         nMuons      = cms.int32(2),        # require minimum number of muons in the event to pass the following dimu pair criteria
                         isGlobal    =cms.int32(0),
                         isStandAlone=cms.int32(0),
                         isTracker   =cms.int32(1),
                         ptMin       = cms.double(10),
                         etaMax      = cms.double(2.4),
    
                         # HLT trigger info
                         checkTrigger = cms.bool(True),    # For Data we require that at least one HLT trigger in the list fires
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
                         btagNames      = cms.vstring(["pfCombinedInclusiveSecondaryVertexV2BJetTags"]),

                         # Electrons
                         electronColl             = cms.InputTag("slimmedElectrons"),
                         electronCutBasedIdVeto   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                         electronCutBasedIdLoose  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                         electronCutBasedIdMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                         electronCutBasedIdTight  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),

                         # Taus
                         tauColl    = cms.InputTag("slimmedTaus"),
                         tauIDNames = cms.vstring(["decayModeFinding",
                                                   "byLooseCombinedIsolationDeltaBetaCorr3Hits",
                                                   "byMediumCombinedIsolationDeltaBetaCorr3Hits",
                                                   "byTightCombinedIsolationDeltaBetaCorr3Hits",
                                                   "byVLooseIsolationMVArun2v1DBoldDMwLT",
                                                   "byLooseIsolationMVArun2v1DBoldDMwLT",
                                                   "byMediumIsolationMVArun2v1DBoldDMwLT",
                                                   "byTightIsolationMVArun2v1DBoldDMwLT",
                                                   "byVTightIsolationMVArun2v1DBoldDMwLT",
                                                   "againstMuonLoose3",
                                                   "againstMuonTight3",
                                                   "againstElectronVLooseMVA6",
                                                   "againstElectronLooseMVA6",
                                                   "againstElectronMediumMVA6",
                                                   "againstElectronTightMVA6",
                                                   "againstElectronVTightMVA6"]),

                         )
