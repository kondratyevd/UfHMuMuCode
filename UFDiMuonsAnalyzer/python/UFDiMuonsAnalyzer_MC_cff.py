import FWCore.ParameterSet.Config as cms

## Assumes you are looking at miniAOD
DiMuons = cms.EDAnalyzer('UFDiMuonsAnalyzer',
                         
                         isVerbose    = cms.untracked.bool(False),
                         isMonteCarlo = cms.bool(True),
                         doSys        = cms.bool(False),
                         slimOut      = cms.bool(True),

                         ## Event selection cuts
                         skim_nMuons = cms.int32(0),
                         skim_trig   = cms.bool(False),
                         
                         ## HLT trigger info
                         processName  = cms.string("HLT"),
                         ## processName  = cms.string("HLT2"),  ## For reHLT MC samples?
                         ## Unprescaled triggers at the end of 2016
                         ## https://cmswbm.web.cern.ch/cmswbm/cmsdb/servlet/TriggerMode?KEY=l1_hlt_collisions2016/v450
                         trigNames = cms.vstring("HLT_IsoMu22_eta2p1", "HLT_IsoTkMu22_eta2p1", 
                                                    "HLT_IsoMu24", "HLT_IsoTkMu24", 
                                                    "HLT_Mu50", "HLT_TkMu50"),

                         ## trigResults = cms.InputTag("TriggerResults","","HLT"),
                         trigResults = cms.InputTag("TriggerResults","","HLT2"),  ## For reHLT MC samples
                         trigObjs    = cms.InputTag("selectedPatTrigger"),
                         
                         ## Vertex and Beam Spot
                         primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                         beamSpotTag      = cms.InputTag("offlineBeamSpot"),

                         ## Muons
                         muonColl   = cms.InputTag("slimmedMuons"),
                         doSys_KaMu = cms.bool(False),
                         doSys_Roch = cms.bool(True),

                         ## Electrons
                         eleColl     = cms.InputTag("slimmedElectrons"),
                         eleIdVeto   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                         eleIdLoose  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                         eleIdMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                         eleIdTight  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),

                         ## Taus
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
                         ## Jets
                         jetsTag  = cms.InputTag("slimmedJets"),
                         btagName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),

                         ## MET
                         metTag         = cms.InputTag("slimmedMETs"),

                         ## GEN particle collections
                         genJetsTag           = cms.InputTag("slimmedGenJets"),
                         prunedGenParticleTag = cms.InputTag("prunedGenParticles"),
                         packedGenParticleTag = cms.InputTag("packedGenParticles"),
                         
                         ## Object selection
                         vertex_ndof_min = cms.double( 4.0),
                         vertex_rho_max  = cms.double( 2.0),
                         vertex_z_max    = cms.double(24.0),

                         muon_ID        = cms.string("loose"),
                         muon_pT_min    = cms.double(10.0),
                         muon_eta_max   = cms.double( 2.4),
                         muon_trig_dR   = cms.double( 0.1),
                         muon_use_pfIso = cms.bool(True),
                         muon_iso_dR    = cms.double( 0.3),
                         muon_iso_max   = cms.double(999.),

                         ele_ID      = cms.string("loose"),
                         ele_pT_min  = cms.double(10.),
                         ele_eta_max = cms.double(2.5),

                         tau_pT_min  = cms.double(10.),
                         tau_eta_max = cms.double(2.5),
                         
                         jet_ID      = cms.string("loose"),
                         jet_pT_min  = cms.double(20.0),
                         jet_eta_max = cms.double(4.7),

                         ## Event weights
                         PU_wgt_file = cms.string("PU_wgt_2016_Summer16_v0.root"),
                         Trig_eff_3_file = cms.string("EfficienciesAndSF_Period3.root"),
                         Trig_eff_4_file = cms.string("EfficienciesAndSF_Period4.root"),

                         )
