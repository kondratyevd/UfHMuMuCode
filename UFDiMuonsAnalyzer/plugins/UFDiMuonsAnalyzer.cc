
#include "UfHMuMuCode/UFDiMuonsAnalyzer/plugins/UFDiMuonsAnalyzer.h"

// Constructor
UFDiMuonsAnalyzer::UFDiMuonsAnalyzer(const edm::ParameterSet& iConfig):
  _numEvents(0)
{
  // Initialize the weighted count and the trees.
  // Use the file service to make the trees so that it knows to save them.
  _sumEventWeights = 0;
  _outTree = fs->make<TTree>("tree", "myTree");
  _outTreeMetadata = fs->make<TTree>("metadata", "Metadata Tree");

  // Get the collections designated from the python config file and load them into the tokens

  // Boolean switches from config file
  _isVerbose	 = iConfig.getUntrackedParameter<bool>("isVerbose", false);
  _isMonteCarlo	 = iConfig.getParameter         <bool>("isMonteCarlo");
  _doSys         = iConfig.getParameter         <bool>("doSys");
  _slimOut       = iConfig.getParameter         <bool>("slimOut");
  
  // Event selection from config file
  _skim_nMuons = iConfig.getParameter<int>  ("skim_nMuons");
  _skim_trig   = iConfig.getParameter<bool> ("skim_trig");

  // Trigger info
  _processName  = iConfig.getParameter            <std::string> ("processName");
  _trigNames    = iConfig.getParameter<std::vector<std::string>>("trigNames");

  _trigResultsToken = consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("trigResults"));
  _trigObjsToken    = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("trigObjs"));

  // Underlying event
  _beamSpotToken      = consumes<reco::BeamSpot>        (iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  _primaryVertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"));
  _PupInfoToken       = consumes< std::vector<PileupSummaryInfo> >           (edm::InputTag ("slimmedAddPileupInfo"));

  // Muons
  _muonCollToken = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonColl"));

  // Electrons
  _eleCollToken     = consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>("eleColl"));
  _eleIdVetoToken   = consumes< edm::ValueMap<bool> >    (iConfig.getParameter<edm::InputTag>("eleIdVeto"));
  _eleIdLooseToken  = consumes< edm::ValueMap<bool> >    (iConfig.getParameter<edm::InputTag>("eleIdLoose"));
  _eleIdMediumToken = consumes< edm::ValueMap<bool> >    (iConfig.getParameter<edm::InputTag>("eleIdMedium"));
  _eleIdTightToken  = consumes< edm::ValueMap<bool> >    (iConfig.getParameter<edm::InputTag>("eleIdTight"));

  // Taus
  _tauCollToken = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauColl"));
  _tauIDNames   = iConfig.getParameter<std::vector<std::string> >("tauIDNames");

   // Jets / MET
  _jetsToken = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsTag"));
  _btagName  = iConfig.getParameter<std::string>("btagName");
  _metToken  = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metTag"));  // Correct(ed) MET? - AWB 08.11.16

  // GEN objects
  _genJetsToken           = consumes<reco::GenJetCollection>          (iConfig.getParameter<edm::InputTag>("genJetsTag"));
  _prunedGenParticleToken = consumes<reco::GenParticleCollection>     (iConfig.getParameter<edm::InputTag>("prunedGenParticleTag"));
  _packedGenParticleToken = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticleTag"));
  _genEvtInfoToken        = consumes<GenEventInfoProduct>                                  (edm::InputTag ("generator"));

  // Object selection from config file
  _vertex_ndof_min = iConfig.getParameter<double> ("vertex_ndof_min");
  _vertex_rho_max  = iConfig.getParameter<double> ("vertex_rho_max");
  _vertex_z_max    = iConfig.getParameter<double> ("vertex_z_max");

  _muon_ID        = iConfig.getParameter<std::string> ("muon_ID");
  _muon_pT_min    = iConfig.getParameter<double>      ("muon_pT_min");
  _muon_eta_max   = iConfig.getParameter<double>      ("muon_eta_max");
  _muon_trig_dR   = iConfig.getParameter<double>      ("muon_trig_dR");
  _muon_use_pfIso = iConfig.getParameter<bool>        ("muon_use_pfIso");
  _muon_iso_dR    = iConfig.getParameter<double>      ("muon_iso_dR");
  _muon_iso_max   = iConfig.getParameter<double>      ("muon_iso_max");

  _ele_ID      = iConfig.getParameter<std::string> ("ele_ID");
  _ele_pT_min  = iConfig.getParameter<double>      ("ele_pT_min");
  _ele_eta_max = iConfig.getParameter<double>      ("ele_eta_max");

  _tau_pT_min  = iConfig.getParameter<double>       ("tau_pT_min");
  _tau_eta_max = iConfig.getParameter<double>       ("tau_eta_max");

  _jet_ID      = iConfig.getParameter<std::string> ("jet_ID");
  _jet_pT_min  = iConfig.getParameter<double>      ("jet_pT_min");
  _jet_eta_max = iConfig.getParameter<double>      ("jet_eta_max");

  if (_isMonteCarlo) _KaMu_calib = KalmanMuonCalibrator("MC_80X_13TeV");
  else               _KaMu_calib = KalmanMuonCalibrator("DATA_80X_13TeV");
  _doSys_KaMu  = iConfig.getParameter<bool>("doSys_KaMu");

  // Jigger path name for crab
  edm::FileInPath cfg_RochCor("RochCor/Calibration/data/Feb06/config.txt");
  std::string path_RochCor = cfg_RochCor.fullPath().c_str();
  std::string file_RochCor = "/config.txt";
  std::string::size_type find_RochCor = path_RochCor.find(file_RochCor);
  if (find_RochCor != std::string::npos)
    path_RochCor.erase(find_RochCor, file_RochCor.length());

  std::cout << "Rochester correction files located in " << path_RochCor << std::endl;
  _Roch_calib.init(path_RochCor);
  _doSys_Roch = iConfig.getParameter<bool>("doSys_Roch");

  if (_isMonteCarlo) {
    edm::FileInPath path_PU_wgt("UfHMuMuCode/UFDiMuonsAnalyzer/data/Pileup/"+iConfig.getParameter<std::string>("PU_wgt_file"));
    _PU_wgt_file      = new TFile(path_PU_wgt.fullPath().c_str());
    _PU_wgt_hist      = (TH1D*) _PU_wgt_file->Get("PU_wgt");
    _PU_wgt_hist_up   = (TH1D*) _PU_wgt_file->Get("PU_wgt_up");
    _PU_wgt_hist_down = (TH1D*) _PU_wgt_file->Get("PU_wgt_down");
  }

  edm::FileInPath path_IsoMu_eff_3("UfHMuMuCode/UFDiMuonsAnalyzer/data/MuonTrig/"+iConfig.getParameter<std::string>("Trig_eff_3_file"));
  edm::FileInPath path_IsoMu_eff_4("UfHMuMuCode/UFDiMuonsAnalyzer/data/MuonTrig/"+iConfig.getParameter<std::string>("Trig_eff_4_file"));
  _IsoMu_eff_3_file = new TFile(path_IsoMu_eff_3.fullPath().c_str());
  _IsoMu_eff_4_file = new TFile(path_IsoMu_eff_4.fullPath().c_str());
  _IsoMu_eff_3_hist = (TH2F*) _IsoMu_eff_3_file->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
  _IsoMu_eff_4_hist = (TH2F*) _IsoMu_eff_4_file->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA");
  
} // End constructor: UFDiMuonsAnalyzer::UFDiMuonsAnalyzer

// Destructor
UFDiMuonsAnalyzer::~UFDiMuonsAnalyzer() {

  if (_isMonteCarlo)
    _PU_wgt_file->Close();

  _IsoMu_eff_3_file->Close();
  _IsoMu_eff_4_file->Close();

}

// Called once per event
void UFDiMuonsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // ------------------------
  // COUNT EVENTS AND WEIGHTS 
  // ------------------------
  _numEvents++;
  if (_isVerbose) std::cout << "\n\n A N A L Y Z I N G   E V E N T = " << _numEvents << std::endl << std::endl;
  
  if (!_isMonteCarlo)
    _sumEventWeights += 1;
  else {
    // The generated weight. Due to the interference of terms in QM in the
    // NLO simulations there are negative weights that need to be accounted for. 
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(_genEvtInfoToken, genEvtInfo );
    _GEN_wgt = (genEvtInfo->weight() > 0) ? 1 : -1;  // Why don't we use the decimal weight? - AWB 08.11.16
    _sumEventWeights += _GEN_wgt;
  }

  // -------------------
  // HLT TRIGGER HANDLES
  // -------------------
  if (_isVerbose) std::cout << "\nAccessing HLT info" << std::endl;
  edm::Handle<edm::TriggerResults> trigResultsHandle;
  iEvent.getByToken(_trigResultsToken, trigResultsHandle);
  if (!trigResultsHandle.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }

  edm::Handle<pat::TriggerObjectStandAloneCollection> trigObjsHandle;
  iEvent.getByToken(_trigObjsToken, trigObjsHandle);
  if (!trigObjsHandle.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerObjects product from Event!" << std::endl;
    return;
  }

  // If all the HLT paths are not fired, will discard the event immediately
  // Do we store trigger decisions somewhere? - AWB 12.11.16
  if ( _skim_trig && !isHltPassed(iEvent, iSetup, trigResultsHandle, _trigNames) ) {
    if (_isVerbose) std::cout << "None of the HLT paths fired -> discard the event\n";
    return;
  }

  // ----------------
  // RUN / EVENT INFO
  // ----------------
  if (_isVerbose) std::cout << "\nFilling EventInfo" << std::endl;
  FillEventInfo( _eventInfo, iEvent );

  // ----------------------------
  // BEAMSPOT / VERTICES / PILEUP
  // ----------------------------
  if (_isVerbose) std::cout << "\nFilling VertexInfo" << std::endl;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(_beamSpotToken, beamSpotHandle);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(_primaryVertexToken, vertices);

  reco::VertexCollection verticesSelected = SelectVertices( vertices, _vertex_ndof_min,
							    _vertex_rho_max, _vertex_z_max );
  if (verticesSelected.size() == 0) {
    std::cout << "BUGGY EVENT!  There are no good vertices!  Skipping ..." << std::endl;
    return;
  }
  reco::Vertex primaryVertex = verticesSelected.at(0);
  
  bool _onlyPV = true;  // Only fill primary vertex
  FillVertexInfos( _vertexInfos, _nVertices, verticesSelected, _onlyPV );
  
  // -----
  // MUONS
  // -----
  if (_isVerbose) std::cout << "\nFilling MuonInfo" << std::endl;
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(_muonCollToken, muons);

  pat::MuonCollection muonsSelected = SelectMuons( muons, primaryVertex, _muon_ID,
						   _muon_pT_min, _muon_eta_max, _muon_trig_dR,
						   _muon_use_pfIso, _muon_iso_dR, _muon_iso_max );
  // Throw away event if there are too few muons
  if ( muonsSelected.size() < (unsigned int) _skim_nMuons )
    return;

  // Sort the selected muons by pT
  sort(muonsSelected.begin(), muonsSelected.end(), sortMuonsByPt);
  
  // Get GEN muons for Rochester corrections
  edm::Handle<reco::GenParticleCollection> genPartons;
  iEvent.getByToken(_prunedGenParticleToken, genPartons);

  FillMuonInfos( _muonInfos, muonsSelected, primaryVertex, verticesSelected.size(), beamSpotHandle, 
		 iEvent, iSetup, trigObjsHandle, trigResultsHandle, _trigNames,
		 _muon_trig_dR, _muon_use_pfIso, _muon_iso_dR, !(_isMonteCarlo), 
		 _KaMu_calib, _doSys_KaMu, _Roch_calib, _doSys_Roch, genPartons ); 
  _nMuons = _muonInfos.size();

  CalcTrigEff( _IsoMu_eff_3, _IsoMu_eff_3_up, _IsoMu_eff_3_down, 
	       _IsoMu_eff_3_hist, _muonInfos, false );
  CalcTrigEff( _IsoMu_eff_4, _IsoMu_eff_4_up, _IsoMu_eff_4_down, 
	       _IsoMu_eff_4_hist, _muonInfos, false );
  CalcTrigEff( _IsoMu_eff_bug, _IsoMu_eff_bug_up, _IsoMu_eff_bug_down, 
	       _IsoMu_eff_3_hist, _muonInfos, true );


  // ------------
  // DIMUON PAIRS
  // ------------
  if (_isVerbose) std::cout << "\nFilling PairInfo" << std::endl;
  FillPairInfos( _pairInfos, _muonInfos );
  _nPairs = _pairInfos.size();
  // Throw away events with only low-mass pairs
  if ( _skim_nMuons == 2 && _nPairs == 1 )
    if ( _pairInfos.at(0).mass < 12 )
      return;

  // -----------------------
  // MC PILEUP / GEN WEIGHTS
  // -----------------------
  if (_isVerbose) std::cout << "\nAccessing piluep info" << std::endl;
  _nPU = -1;
  _PU_wgt      = 1.0;
  _PU_wgt_up   = 1.0;
  _PU_wgt_down = 1.0;
  if (_isMonteCarlo) {
    edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(_PupInfoToken, PupInfo);
    
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if (PVI->getBunchCrossing() == 0) { 
	_nPU = PVI->getTrueNumInteractions();
	_PU_wgt      = _PU_wgt_hist->GetBinContent(_nPU);
	_PU_wgt_up   = _PU_wgt_hist_up->GetBinContent(_nPU);
	_PU_wgt_down = _PU_wgt_hist_down->GetBinContent(_nPU);
	continue;
      }
    }
  }
  
  // ---------
  // ELECTRONS
  // ---------
  if (_isVerbose) std::cout << "\nFilling EleInfo" << std::endl;
  edm::Handle<edm::View<pat::Electron>>  eles;
  edm::Handle<edm::ValueMap<bool>>       ele_id_veto;
  edm::Handle<edm::ValueMap<bool>>       ele_id_loose;
  edm::Handle<edm::ValueMap<bool>>       ele_id_medium;
  edm::Handle<edm::ValueMap<bool>>       ele_id_tight;

  iEvent.getByToken(_eleCollToken,     eles);
  iEvent.getByToken(_eleIdVetoToken,   ele_id_veto); 
  iEvent.getByToken(_eleIdLooseToken,  ele_id_loose); 
  iEvent.getByToken(_eleIdMediumToken, ele_id_medium); 
  iEvent.getByToken(_eleIdTightToken,  ele_id_tight); 

  std::vector<std::array<bool, 4>> ele_ID_pass;
  pat::ElectronCollection elesSelected = SelectEles( eles, primaryVertex, ele_id_veto,
  						     ele_id_loose, ele_id_medium, ele_id_tight,
  						     _ele_ID, _ele_pT_min, _ele_eta_max,
						     ele_ID_pass );
  
  // Sort the selected electrons by pT
  sort(elesSelected.begin(), elesSelected.end(), sortElesByPt);
  
  FillEleInfos( _eleInfos, elesSelected, primaryVertex, iEvent, ele_ID_pass );
  _nEles = _eleInfos.size();

  // // ----
  // // TAUS
  // // ----
  // if (_isVerbose) std::cout << "\nFilling TauInfo" << std::endl;
  // edm::Handle<pat::TauCollection> taus;
  // iEvent.getByToken(_tauCollToken, taus);

  // pat::TauCollection tausSelected = SelectTaus( taus, _tau_pT_min, _tau_eta_max );

  // // Sort the selected taus by pT
  // sort(tausSelected.begin(), tausSelected.end(), sortTausByPt);

  // FillTauInfos( _tauInfos, tausSelected, _tauIDNames );  // Sort first? - AWB 09.11.16
  // _nTaus = _tauInfos.size();
    
  // ----
  // JETS
  // ----
  if (_isVerbose) std::cout << "\nFilling JetInfo" << std::endl;
  edm::Handle < pat::JetCollection > jets;
  if(!_jetsToken.isUninitialized()) 
    iEvent.getByToken(_jetsToken, jets);

  // Following https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  // Could use ESTransientHandle? https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHowToGetDataFromES#Optimizing_memory_usage
  iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  
  // Might we need to re-compute the central scale via this recipe? - AWB 14.11.16
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets

  pat::JetCollection jetsSelected = SelectJets( jets, JetCorPar, "none",
						_jet_ID, _jet_pT_min, _jet_eta_max );
  
  sort(jetsSelected.begin(), jetsSelected.end(), sortJetsByPt);

  FillJetInfos( _jetInfos, _nJetsFwd, 
		_nBLoose, _nBMed, _nBTight, 
		jetsSelected, _btagName );
  _nJets = _jetInfos.size();
  _nJetsCent = _nJets - _nJetsFwd;
  if (_slimOut) FillSlimJetInfos( _slimJetInfos, _jetInfos );
  
  // Alternate collections corresponding to jet energy scale up / down
  if (_doSys) {
    pat::JetCollection jets_JES_up   = SelectJets( jets, JetCorPar, "JES_up",
						   _jet_ID, _jet_pT_min, _jet_eta_max );
    pat::JetCollection jets_JES_down = SelectJets( jets, JetCorPar, "JES_down",
						   _jet_ID, _jet_pT_min, _jet_eta_max );
    // pat::JetCollection jets_JER_up   = SelectJets( jets, JetCorPar, "JER_up",
    // 						   _jet_ID, _jet_pT_min, _jet_eta_max );
    // pat::JetCollection jets_JER_down = SelectJets( jets, JetCorPar, "JER_down",
    // 						   _jet_ID, _jet_pT_min, _jet_eta_max );

    // delete * JetCorPar; // Avoid nasty memory leak

    sort(jets_JES_up.begin(),   jets_JES_up.end(),   sortJetsByPt);
    sort(jets_JES_down.begin(), jets_JES_down.end(), sortJetsByPt);
    // sort(jets_JER_up.begin(),   jets_JER_up.end(),   sortJetsByPt);
    // sort(jets_JER_down.begin(), jets_JER_down.end(), sortJetsByPt);
    
    FillJetInfos( _jetInfos_JES_up,   _nJetsFwd_JES_up,   
		  _nBLoose_JES_up, _nBMed_JES_up, _nBTight_JES_up, 
		  jets_JES_up,   _btagName );
    FillJetInfos( _jetInfos_JES_down, _nJetsFwd_JES_down, 
		  _nBLoose_JES_down, _nBMed_JES_down, _nBTight_JES_down, 
		  jets_JES_down, _btagName );
    _nJets_JES_up   = _jetInfos_JES_up.size();
    _nJets_JES_down = _jetInfos_JES_down.size();
    _nJetsCent_JES_up   = _nJets_JES_up   - _nJetsFwd_JES_up;
    _nJetsCent_JES_down = _nJets_JES_down - _nJetsFwd_JES_down;
    if (_slimOut) FillSlimJetInfos( _slimJetInfos_JES_up,   _jetInfos_JES_up   );
    if (_slimOut) FillSlimJetInfos( _slimJetInfos_JES_down, _jetInfos_JES_down );

    // FillJetInfos( _jetInfos_JER_up,   _nJetsFwd_JER_up,   
    // 		  _nBLoose_JER_up, _nBMed_JER_up, _nBTight_JER_up, 
    // 		  jets_JER_up,   _btagName );
    // FillJetInfos( _jetInfos_JER_down, _nJetsFwd_JER_down, 
    // 		  _nBLoose_JER_down, _nBMed_JER_down, _nBTight_JER_down, 
    // 		  jets_JER_down, _btagName );
    // _nJets_JER_up   = _jetInfos_JER_up.size();
    // _nJets_JER_down = _jetInfos_JER_down.size();
    // _nJetsCent_JER_up   = _nJets_JER_up   - _nJetsFwd_JER_up;
    // _nJetsCent_JER_down = _nJets_JER_down - _nJetsFwd_JER_down;
    // if (_slimOut) FillSlimJetInfos( _slimJetInfos_JER_up,   _jetInfos_JER_up   );
    // if (_slimOut) FillSlimJetInfos( _slimJetInfos_JER_down, _jetInfos_JER_down );
  }


  // As of January 12, not usable - AWB 16.01.17
  // https://indico.cern.ch/event/598195/contributions/2429716/attachments/1394701/2125692/met_tails_ppd.pdf
  // ---
  // MET
  // ---
  if (_isVerbose) std::cout << "\nFilling MetInfo" << std::endl;
  edm::Handle < pat::METCollection > mets;
  if (!_metToken.isUninitialized()) 
    iEvent.getByToken(_metToken, mets);
  
  FillMetInfo( _metInfo, mets, iEvent );

  // Not using all the correct filters - will need to if we use MET in the analysis - AWB 14.11.16
  //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
  // Should propagate JES uncertainties to MET, if we end up using MET - AWB 14.11.16
  //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#Type_I_Correction_Propagation_of
  //   https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/src/MiniAODHelper.cc#L1085


  // ---
  // MHT
  // ---
  if (_isVerbose) std::cout << "\nFilling MhtInfo" << std::endl;
  
  FillMhtInfo( _mhtInfo, _muonInfos, _eleInfos, _jetInfos ); 
  FillMhtInfo( _mhtInfo_JES_up, _muonInfos, _eleInfos, _jetInfos_JES_up ); 
  FillMhtInfo( _mhtInfo_JES_down, _muonInfos, _eleInfos, _jetInfos_JES_down ); 
  // FillMhtInfo( _mhtInfo_JER_up, _muonInfos, _eleInfos, _jetInfos_JER_up ); 
  // FillMhtInfo( _mhtInfo_JER_down, _muonInfos, _eleInfos, _jetInfos_JER_down ); 


  // -------------
  // GEN PARTICLES
  // -------------
  if (_isMonteCarlo) {
    // Parents
    if (_isVerbose) std::cout << "\nFilling GenParentInfo" << std::endl;

    FillGenParentInfos( _genParentInfos, genPartons, 
			std::vector<int> {6, 22, 23, 24, 25}, 
			_isMonteCarlo );
    _nGenParents = _genParentInfos.size();

    // Muons
    if (_isVerbose) std::cout << "\nFilling GenMuonInfo" << std::endl;
    FillGenMuonInfos( _genMuonInfos, _genParentInfos, genPartons, _isMonteCarlo );
    _nGenMuons = _genMuonInfos.size();
    _nGenParents = _genParentInfos.size();

    if (_isVerbose) std::cout << "\nFilling GenPairInfo" << std::endl;
    FillGenPairInfos( _genPairInfos, _genMuonInfos );
    _nGenPairs = _genPairInfos.size();

    // Jets
    if (_isVerbose) std::cout << "\nFilling GenJetInfo" << std::endl;
    edm::Handle < reco::GenJetCollection > genJets;
    if (!_genJetsToken.isUninitialized()) 
      iEvent.getByToken(_genJetsToken, genJets);

    FillGenJetInfos( _genJetInfos, genJets, _isMonteCarlo );
    _nGenJets = _genJetInfos.size();

  } // End conditional: if (_isMonteCarlo)


  // ============================
  // Store everything in an NTuple
  // ============================
  if (_isVerbose) std::cout << "\nFilling tree" << std::endl;
  _outTree->Fill();
  if (_isVerbose) std::cout << "\nD O N E !!! WITH EVENT!" << std::endl;
  return;  
} // End void UFDiMuonsAnalyzer::analyze

	
///////////////////////////////////////////////////////////
// BeginJob ==============================================
//////////////////////////////////////////////////////////

// Method called once each job just before starting event loop
// Set up TTrees where we save all of the info gathered in the analyzer
void UFDiMuonsAnalyzer::beginJob() {

  displaySelection();

  // // Include user-defined structs (or classes) in the output tree
  // gROOT->ProcessLine("#include <UfHMuMuCode/UFDiMuonsAnalyzer/interface/LinkDef.h>");

  if (_slimOut) {
    _outTree->Branch("jets",          (SlimJetInfos*)      &_slimJetInfos          );
    _outTree->Branch("jets_JES_up",   (SlimJetInfos*)      &_slimJetInfos_JES_up   );
    _outTree->Branch("jets_JES_down", (SlimJetInfos*)      &_slimJetInfos_JES_down );
    // _outTree->Branch("jets_JER_up",   (SlimJetInfos*)      &_slimJetInfos_JER_up   );
    // _outTree->Branch("jets_JER_down", (SlimJetInfos*)      &_slimJetInfos_JER_down );
  } else {
  _outTree->Branch("jets",          (JetInfos*)      &_jetInfos          );
  _outTree->Branch("jets_JES_up",   (JetInfos*)      &_jetInfos_JES_up   );
  _outTree->Branch("jets_JES_down", (JetInfos*)      &_jetInfos_JES_down );
  // _outTree->Branch("jets_JER_up",   (JetInfos*)      &_jetInfos_JER_up   );
  // _outTree->Branch("jets_JER_down", (JetInfos*)      &_jetInfos_JER_down );

  _outTree->Branch("vertices",      (VertexInfos*)   &_vertexInfos       );
  _outTree->Branch("taus",          (TauInfos*)      &_tauInfos          );
  }

  // Set up the _outTree branches
  _outTree->Branch("event",         (EventInfo*)     &_eventInfo         );
  _outTree->Branch("muons",         (MuonInfos*)     &_muonInfos         );
  _outTree->Branch("pairs",         (PairInfos*)     &_pairInfos         );
  _outTree->Branch("eles",          (EleInfos*)      &_eleInfos          );
  // _outTree->Branch("met",           (MetInfo*)       &_metInfo           );
  _outTree->Branch("mht",           (MhtInfo*)       &_mhtInfo           );
  _outTree->Branch("mht_JES_up",    (MhtInfo*)       &_mhtInfo_JES_up    );
  _outTree->Branch("mht_JES_down",  (MhtInfo*)       &_mhtInfo_JES_down  );
  // _outTree->Branch("mht_JER_up",    (MhtInfo*)       &_mhtInfo_JER_up    );
  // _outTree->Branch("mht_JER_down",  (MhtInfo*)       &_mhtInfo_JER_down  );

  _outTree->Branch("nVertices",          (int*) &_nVertices          );
  _outTree->Branch("nMuons",             (int*) &_nMuons             );
  _outTree->Branch("nPairs",             (int*) &_nPairs             );
  _outTree->Branch("nEles",              (int*) &_nEles              );
  _outTree->Branch("nTaus",              (int*) &_nTaus              );
  _outTree->Branch("nJets",              (int*) &_nJets              );
  _outTree->Branch("nJetsCent",          (int*) &_nJetsCent          );
  _outTree->Branch("nJetsFwd",           (int*) &_nJetsFwd           );
  _outTree->Branch("nBLoose",            (int*) &_nBLoose            );
  _outTree->Branch("nBMed",              (int*) &_nBMed              );
  _outTree->Branch("nBTight",            (int*) &_nBTight            );
  _outTree->Branch("nJets_JES_up",       (int*) &_nJets_JES_up       );
  _outTree->Branch("nJetsCent_JES_up",   (int*) &_nJetsCent_JES_up   );
  _outTree->Branch("nBLoose_JES_up",     (int*) &_nBLoose_JES_up     );
  _outTree->Branch("nBMed_JES_up",       (int*) &_nBMed_JES_up       );
  _outTree->Branch("nBTight_JES_up",     (int*) &_nBTight_JES_up     );
  _outTree->Branch("nJetsFwd_JES_up",    (int*) &_nJetsFwd_JES_up    );
  _outTree->Branch("nJets_JES_down",     (int*) &_nJets_JES_down     );
  _outTree->Branch("nJetsCent_JES_down", (int*) &_nJetsCent_JES_down );
  _outTree->Branch("nJetsFws_JES_down",  (int*) &_nJetsFwd_JES_down  );
  _outTree->Branch("nBLoose_JES_down",   (int*) &_nBLoose_JES_down   );
  _outTree->Branch("nBMed_JES_down",     (int*) &_nBMed_JES_down     );
  _outTree->Branch("nBTight_JES_down",   (int*) &_nBTight_JES_down   );
  // _outTree->Branch("nJets_JER_up",       (int*) &_nJets_JER_up       );
  // _outTree->Branch("nJetsCent_JER_up",   (int*) &_nJetsCent_JER_up   );
  // _outTree->Branch("nJetsFwd_JER_up",    (int*) &_nJetsFwd_JER_up    );
  // _outTree->Branch("nBLoose_JER_up",     (int*) &_nBLoose_JER_up     );
  // _outTree->Branch("nBMed_JER_up",       (int*) &_nBMed_JER_up       );
  // _outTree->Branch("nBTight_JER_up",     (int*) &_nBTight_JER_up     );
  // _outTree->Branch("nJets_JER_down",     (int*) &_nJets_JER_down     );
  // _outTree->Branch("nJetsCent_JER_down", (int*) &_nJetsCent_JER_down );
  // _outTree->Branch("nJetsFws_JER_down",  (int*) &_nJetsFwd_JER_down  );
  // _outTree->Branch("nBLoose_JER_down",   (int*) &_nBLoose_JER_down   );
  // _outTree->Branch("nBMed_JER_down",     (int*) &_nBMed_JER_down     );
  // _outTree->Branch("nBTight_JER_down",   (int*) &_nBTight_JER_down   );

  _outTree->Branch("hltPaths",      &_trigNames);
  _outTree->Branch("btagName",      &_btagName   );
  _outTree->Branch("tauIDNames",    &_tauIDNames  );

  // Weights and efficiencies
  _outTree->Branch("IsoMu_eff_3",        &_IsoMu_eff_3,        "IsoMu_eff_3/F"        );
  _outTree->Branch("IsoMu_eff_3_up",     &_IsoMu_eff_3_up,     "IsoMu_eff_3_up/F"     );
  _outTree->Branch("IsoMu_eff_3_down",   &_IsoMu_eff_3_down,   "IsoMu_eff_3_down/F"   );
  _outTree->Branch("IsoMu_eff_4",        &_IsoMu_eff_4,        "IsoMu_eff_4/F"        );
  _outTree->Branch("IsoMu_eff_4_up",     &_IsoMu_eff_4_up,     "IsoMu_eff_4_up/F"     );
  _outTree->Branch("IsoMu_eff_4_down",   &_IsoMu_eff_4_down,   "IsoMu_eff_4_down/F"   );
  _outTree->Branch("IsoMu_eff_bug",      &_IsoMu_eff_bug,      "IsoMu_eff_bug/F"      );
  _outTree->Branch("IsoMu_eff_bug_up",   &_IsoMu_eff_bug_up,   "IsoMu_eff_bug_up/F"   );
  _outTree->Branch("IsoMu_eff_bug_down", &_IsoMu_eff_bug_down, "IsoMu_eff_bug_down/F" );


  // MC information
  if (_isMonteCarlo) {

    _outTree->Branch("nPU", 	    &_nPU,         "nPU/I"         );
    _outTree->Branch("PU_wgt",      &_PU_wgt,      "PU_wgt/F"      );
    _outTree->Branch("PU_wgt_up",   &_PU_wgt_up,   "PU_wgt_up/F"   );
    _outTree->Branch("PU_wgt_down", &_PU_wgt_down, "PU_wgt_down/F" );
    _outTree->Branch("GEN_wgt",     &_GEN_wgt,     "GEN_wgt/I"     );
    
    _outTree->Branch("nGenParents", (int*) &_nGenParents );
    _outTree->Branch("nGenMuons",   (int*) &_nGenMuons   );
    _outTree->Branch("nGenPairs",   (int*) &_nGenPairs   );
    
    _outTree->Branch("genParents", (GenParentInfos*) &_genParentInfos );
    _outTree->Branch("genMuons",   (GenMuonInfos*)   &_genMuonInfos   );
    _outTree->Branch("genPairs",   (GenPairInfos*)   &_genPairInfos   );

    if (!_slimOut) {
      _outTree->Branch("nGenJets", (int*)         &_nGenJets    );
      _outTree->Branch("genJets",  (GenJetInfos*) &_genJetInfos );
    }  // End conditional: if (!_slimOut)

  }

}

///////////////////////////////////////////////////////////
// BeginJob ==============================================
//////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::endJob() 
{
// Method called once each job just after ending the event loop
// Set up the meta data ttree and save it.

  std::cout << "Total Number of Events Read: "<< _numEvents << std::endl <<std::endl;
  std::cout << "Number of events weighted: "  << _sumEventWeights << std::endl <<std::endl;

  std::cout<<"number of dimuon candidates: "
           <<_outTree->GetEntries()<<std::endl;

  // create the metadata tree branches
  _outTreeMetadata->Branch("originalNumEvents", &_numEvents,        "originalNumEvents/I");
  _outTreeMetadata->Branch("sumEventWeights",   &_sumEventWeights,  "sumEventWeights/I");
  _outTreeMetadata->Branch("isMonteCarlo",      &_isMonteCarlo,     "isMonteCarlo/O");

  _outTreeMetadata->Branch("trigNames",  "std::vector< std::string >", &_trigNames);
  _outTreeMetadata->Branch("btagName",   "std::string",                &_btagName);
  _outTreeMetadata->Branch("tauIDNames", "std::vector< std::string >", &_tauIDNames);

  _outTreeMetadata->Fill();

}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::isHltPassed(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				    const edm::Handle<edm::TriggerResults>& trigResultsHandle,
                                    const std::vector<std::string> desiredTrigNames) 
{
// this method will simply check if the selected HLT path (via triggerName)
// is run and accepted and no error are found
//
// bool true  if (run && accept && !error)
//      false if any other combination
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;


  const boost::regex re("_v[0-9]+");

  const TriggerNames &trigNames = iEvent.triggerNames(*trigResultsHandle);

  const unsigned nTrigs = trigResultsHandle->size();
  for (unsigned iTrig = 0; iTrig < nTrigs; ++iTrig)
  {
    const string trigName = trigNames.triggerName(iTrig);
    string trigNameStripped = boost::regex_replace(trigName,re,"",boost::match_default | boost::format_sed);
    for(std::vector<std::string>::const_iterator desiredTrigName=desiredTrigNames.begin();
            desiredTrigName!=desiredTrigNames.end();desiredTrigName++)
    {
      if (*desiredTrigName == trigNameStripped && trigResultsHandle->accept(iTrig))
      {
        stringstream debugString;
        debugString << "isHltPassed:" <<endl;
        debugString << "  Trigger "<<iTrig<<": "<< trigName << "("<<trigNameStripped<<") passed: "<<trigResultsHandle->accept(iTrig)<<endl;
        debugString << "    Desired Trigger Names: ";
        debugString <<"'"<< *desiredTrigName<<"' ";
        debugString << endl << "    Accept Trigger" << endl;
        LogVerbatim("UFHLTTests") << debugString.str();
        return true;
      }
    }
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::displaySelection() {

  std::cout << "\n\n*** UFDiMuonsAnalyzer Configuration ***\n";

  // cout object / event cuts - AWB 14.11.16
  std::cout << "_vertex_ndof_min = " << _vertex_ndof_min << std::endl;
  std::cout << "_vertex_rho_max  = " << _vertex_rho_max << std::endl;
  std::cout << "_vertex_z_max    = " << _vertex_z_max << std::endl;

  std::cout << "_muon_ID        = " << _muon_ID << std::endl;
  std::cout << "_muon_pT_min    = " << _muon_pT_min << std::endl;
  std::cout << "_muon_eta_max   = " << _muon_eta_max << std::endl;
  std::cout << "_muon_trig_dR   = " << _muon_trig_dR << std::endl;
  std::cout << "_muon_use_pfIso = " << _muon_use_pfIso << std::endl;
  std::cout << "_muon_iso_dR    = " << _muon_iso_dR << std::endl;
  std::cout << "_muon_iso_max   = " << _muon_iso_max << std::endl;

  std::cout << "_ele_ID      = " << _ele_ID << std::endl;
  std::cout << "_ele_pT_min  = " << _ele_pT_min << std::endl;
  std::cout << "_ele_eta_max = " << _ele_eta_max << std::endl;

  std::cout << "_tau_pT_min  = " << _tau_pT_min << std::endl;
  std::cout << "_tau_eta_max = " << _tau_eta_max << std::endl;

  std::cout << "_jet_ID      = " << _jet_ID << std::endl;
  std::cout << "_jet_pT_min  = " << _jet_pT_min << std::endl;
  std::cout << "_jet_eta_max = " << _jet_eta_max << std::endl;

  // module config parameters
  std::cout << " - _skim_trig: " << _skim_trig << std::endl;
  std::cout << " - Triggers To Probe:\n";
  unsigned int trigSize = _trigNames.size();
  for (unsigned int i=0; i < trigSize; i++) 
    std::cout << "    * trigNames["<<i<<"]: " << _trigNames[i] << std::endl;
  
  std::cout << std::endl << std::endl;

}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::sortMuonsByPt    (pat::Muon i,     pat::Muon j    ) { return (i.pt() > j.pt()); }
bool UFDiMuonsAnalyzer::sortElesByPt     (pat::Electron i, pat::Electron j) { return (i.pt() > j.pt()); }
bool UFDiMuonsAnalyzer::sortTausByPt     (pat::Tau i,      pat::Tau j     ) { return (i.pt() > j.pt()); }
bool UFDiMuonsAnalyzer::sortJetsByPt     (pat::Jet i,      pat::Jet j     ) { return (i.pt() > j.pt()); }

