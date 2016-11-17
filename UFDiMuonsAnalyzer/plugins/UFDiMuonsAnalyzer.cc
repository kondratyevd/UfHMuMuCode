
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
  _isVerbose	= iConfig.getUntrackedParameter<bool>("isVerbose", false);
  _isMonteCarlo	= iConfig.getParameter         <bool>("isMonteCarlo");
  _doSyst       = iConfig.getParameter         <bool>("doSyst");
  
  // Event selection from config file
  _skim_nMuons = iConfig.getParameter<int>  ("skim_nMuons");
  _skim_trigger = iConfig.getParameter<bool>("skim_trigger");

  // Trigger info
  _processName  = iConfig.getParameter            <std::string> ("processName");
  _triggerNames = iConfig.getParameter<std::vector<std::string>>("triggerNames");

  _triggerResultsToken = consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerResults"));
  _triggerObjsToken    = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerObjs"));

  // Underlying event
  _beamSpotToken      = consumes<reco::BeamSpot>        (iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  _primaryVertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"));
  _PupInfoToken       = consumes< std::vector<PileupSummaryInfo> >           (edm::InputTag ("slimmedAddPileupInfo"));

  // Muons
  _muonCollToken = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonColl"));

  // Electrons
  _electronCollToken     = consumes<edm::View<pat::Electron>> (iConfig.getParameter<edm::InputTag>("electronColl"));
  _electronIdVetoToken   = consumes< edm::ValueMap<bool> >   (iConfig.getParameter<edm::InputTag>("electronIdVeto"));
  _electronIdLooseToken  = consumes< edm::ValueMap<bool> >   (iConfig.getParameter<edm::InputTag>("electronIdLoose"));
  _electronIdMediumToken = consumes< edm::ValueMap<bool> >   (iConfig.getParameter<edm::InputTag>("electronIdMedium"));
  _electronIdTightToken  = consumes< edm::ValueMap<bool> >   (iConfig.getParameter<edm::InputTag>("electronIdTight"));

  // Taus
  _tauCollToken = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauColl"));
  _tauIDNames   = iConfig.getParameter<std::vector<std::string> >("tauIDNames");

   // Jets / MET
  _jetsToken = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsTag"));
  _btagNames   = iConfig.getParameter<std::vector<std::string>>("btagNames");
  _metToken    = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metTag"));  // Correct(ed) MET? - AWB 08.11.16

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

  _electron_ID      = iConfig.getParameter<std::string> ("electron_ID");
  _electron_pT_min  = iConfig.getParameter<double>      ("electron_pT_min");
  _electron_eta_max = iConfig.getParameter<double>      ("electron_eta_max");

  _tau_pT_min  = iConfig.getParameter<double>       ("tau_pT_min");
  _tau_eta_max = iConfig.getParameter<double>       ("tau_eta_max");

  _jet_ID      = iConfig.getParameter<std::string> ("jet_ID");
  _jet_pT_min  = iConfig.getParameter<double>      ("jet_pT_min");
  _jet_eta_max = iConfig.getParameter<double>      ("jet_eta_max");

} // End constructor: UFDiMuonsAnalyzer::UFDiMuonsAnalyzer

// Destructor
UFDiMuonsAnalyzer::~UFDiMuonsAnalyzer() {}

// Called once per event
void UFDiMuonsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // ------------------------
  // COUNT EVENTS AND WEIGHTS 
  // ------------------------
  _numEvents++;
  if (_isVerbose) 
    std::cout << "\n\n A N A L Y Z I N G   E V E N T = " 
	      << _numEvents << std::endl << std::endl;
  
  if (!_isMonteCarlo)
    _sumEventWeights += 1;
  else {
    // The generated weight. Due to the interference of terms in QM in the
    // NLO simulations there are negative weights that need to be accounted for. 
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(_genEvtInfoToken, genEvtInfo );
    _genWeight = (genEvtInfo->weight() > 0)? 1 : -1;  // Why don't we use the decimal weight? - AWB 08.11.16
    _sumEventWeights += _genWeight;
  }

  // -------------------
  // HLT TRIGGER HANDLES
  // -------------------
  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(_triggerResultsToken, triggerResultsHandle);
  if (!triggerResultsHandle.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjsHandle;
  iEvent.getByToken(_triggerObjsToken, triggerObjsHandle);
  if (!triggerObjsHandle.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerObjects product from Event!" << std::endl;
    return;
  }

  // If all the HLT paths are not fired, will discard the event immediately
  // Do we store trigger decisions somewhere? - AWB 12.11.16
  if ( _skim_trigger && !isHltPassed(iEvent, iSetup, triggerResultsHandle, _triggerNames) ) {
    if (_isVerbose) std::cout << "None of the HLT paths fired -> discard the event\n";
    return;
  }

  // ----------------
  // RUN / EVENT INFO
  // ----------------
  _eventInfo.init();
  _eventInfo.fill(iEvent);

  // ----------------------------
  // BEAMSPOT / VERTICES / PILEUP
  // ----------------------------
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(_beamSpotToken, beamSpotHandle);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(_primaryVertexToken, vertices);

  reco::VertexCollection verticesSelected = _vertexInfo.select( vertices, _vertex_ndof_min,
								_vertex_rho_max, _vertex_z_max );
  reco::Vertex primaryVertex = verticesSelected.at(0);
  
  _vertexInfo.init();
  _vertexInfo.fill(verticesSelected);
  
  // -----
  // MUONS
  // -----
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(_muonCollToken, muons);

  pat::MuonCollection muonsSelected = _muonInfo.select( muons, primaryVertex, _muon_ID,
							_muon_pT_min, _muon_eta_max, _muon_trig_dR,
							_muon_use_pfIso, _muon_iso_dR, _muon_iso_max );
  // Throw away event if there are too few muons
  if ( muonsSelected.size() < (unsigned int) _skim_nMuons )
    return;

  // Sort the selected muons by pT
  sort(muonsSelected.begin(), muonsSelected.end(), sortMuonsByPt);
  
  _muonInfo.init();
  _muonInfo.fill( muonsSelected, primaryVertex, verticesSelected.size(), beamSpotHandle, 
		  iEvent, iSetup, triggerObjsHandle, triggerResultsHandle, _triggerNames,
		  _muon_trig_dR, _muon_use_pfIso, _muon_iso_dR ); 

  // ------------
  // DIMUON PAIRS
  // ------------
  _diMuonInfo.init();
  _diMuonInfo.fill(_muonInfo);

  // -----------------------
  // MC PILEUP / GEN WEIGHTS
  // -----------------------
  _nPU = -1;
  if (_isMonteCarlo) {
    edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(_PupInfoToken, PupInfo);
    
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if (PVI->getBunchCrossing() == 0) { 
	_nPU = PVI->getTrueNumInteractions();
	continue;
      }
    }
  }
  
  // ---------
  // ELECTRONS
  // ---------
  edm::Handle<edm::View<pat::Electron>>  electrons;
  edm::Handle<edm::ValueMap<bool>>      ele_id_veto;
  edm::Handle<edm::ValueMap<bool>>      ele_id_loose;
  edm::Handle<edm::ValueMap<bool>>      ele_id_medium;
  edm::Handle<edm::ValueMap<bool>>      ele_id_tight;

  iEvent.getByToken(_electronCollToken,     electrons);
  iEvent.getByToken(_electronIdVetoToken,   ele_id_veto); 
  iEvent.getByToken(_electronIdLooseToken,  ele_id_loose); 
  iEvent.getByToken(_electronIdMediumToken, ele_id_medium); 
  iEvent.getByToken(_electronIdTightToken,  ele_id_tight); 

  pat::ElectronCollection electronsSelected = _electronInfo.select( electrons, primaryVertex, ele_id_veto,
								    ele_id_loose, ele_id_medium, ele_id_tight,
								    _electron_ID, _electron_pT_min, _electron_eta_max );
  
  // Sort the selected electrons by pT
  sort(electronsSelected.begin(), electronsSelected.end(), sortElectronsByPt);

  _electronInfo.init();
  _electronInfo.fill( electronsSelected, primaryVertex, iEvent, ele_id_veto,
		      ele_id_loose, ele_id_medium, ele_id_tight );

  // ----
  // TAUS
  // ----
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(_tauCollToken, taus);

  pat::TauCollection tausSelected = _tauInfo.select( taus, _tau_pT_min, _tau_eta_max );

  // Sort the selected taus by pT
  sort(tausSelected.begin(), tausSelected.end(), sortTausByPt);

  _tauInfo.init();
  _tauInfo.fill(tausSelected, _tauIDNames);  // Sort first? - AWB 09.11.16

  // ----
  // JETS
  // ----
  edm::Handle < pat::JetCollection > jets;
  if(!_jetsToken.isUninitialized()) 
    iEvent.getByToken(_jetsToken, jets);

  // Following https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  
  // Might we need to re-compute the central scale via this recipe? - AWB 14.11.16
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets

  pat::JetCollection jetsSelected = _jetInfo.select( jets, JetCorPar, "none",
						     _jet_ID, _jet_pT_min, _jet_eta_max );
  
  sort(jetsSelected.begin(), jetsSelected.end(), sortJetsByPt);

  _jetInfo.init();
  _jetInfo.fill(jetsSelected, _btagNames);

  // Alternate collections corresponding to jet energy scale up / down
  if (_doSyst) {
    pat::JetCollection jets_JES_up   = _jetInfo.select( jets, JetCorPar, "JES_up",
							_jet_ID, _jet_pT_min, _jet_eta_max );
    pat::JetCollection jets_JES_down = _jetInfo.select( jets, JetCorPar, "JES_down",
							_jet_ID, _jet_pT_min, _jet_eta_max );

    sort(jets_JES_up.begin(),   jets_JES_up.end(),   sortJetsByPt);
    sort(jets_JES_down.begin(), jets_JES_down.end(), sortJetsByPt);
    
    _jetInfo_JES_up.init();
    _jetInfo_JES_down.init();

    _jetInfo_JES_up.fill  (jets_JES_up,   _btagNames);
    _jetInfo_JES_down.fill(jets_JES_down, _btagNames);
  }


  // ---
  // MET
  // ---
  edm::Handle < pat::METCollection > mets;
  if (!_metToken.isUninitialized()) 
    iEvent.getByToken(_metToken, mets);
  
  _metInfo.init();
  _metInfo.fill(mets, iEvent);

  // Not using all the correct filters - will need to if we use MET in the analysis - AWB 14.11.16
  //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
  // Should propagate JES uncertainties to MET, if we end up using MET - AWB 14.11.16
  //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections#Type_I_Correction_Propagation_of
  //   https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/src/MiniAODHelper.cc#L1085


  // --------
  // GEN JETS
  // --------
  if (_isMonteCarlo) {
    edm::Handle < reco::GenJetCollection > genJets;
    if (!_genJetsToken.isUninitialized()) 
      iEvent.getByToken(_genJetsToken, genJets);

    _genJetInfo.init();
    _genJetInfo.fill(genJets, _isMonteCarlo);
  }

  // -------------------------------------------
  // MONTE CARLO GEN INFO: MUONS, Gamma, H, W, Z
  // -------------------------------------------
  
  // Make optional - AWB 08.11.16
  if (_isMonteCarlo) {

    // initialize Gamma to default values
    _genGpreFSR.init();  _genM1GpreFSR.init();  _genM2GpreFSR.init();
    _genGpostFSR.init(); _genM1GpostFSR.init(); _genM2GpostFSR.init();

    // initialize Z to default values
    _genZpreFSR.init();  _genM1ZpreFSR.init();  _genM2ZpreFSR.init();
    _genZpostFSR.init(); _genM1ZpostFSR.init(); _genM2ZpostFSR.init();

    // initialize H to default values
    _genHpreFSR.init();  _genM1HpreFSR.init();  _genM2HpreFSR.init();
    _genHpostFSR.init(); _genM1HpostFSR.init(); _genM2HpostFSR.init();

    // initialize W to default values
    _genWpreFSR.init();  _genMWpreFSR.init();
    _genWpostFSR.init(); _genMWpostFSR.init();

    edm::Handle<reco::GenParticleCollection> prunedGenParticles;
    iEvent.getByToken(_prunedGenParticleToken, prunedGenParticles);

    for (reco::GenParticle g: *prunedGenParticles) {
      // Looks for gamma, W, Z, H and the muons from them
      fillBosonAndMuDaughters( g, 
			       _genGpreFSR,  _genM1GpreFSR,  _genM2GpreFSR,
			       _genGpostFSR, _genM1GpostFSR, _genM2GpostFSR,
			       _genZpreFSR,  _genM1ZpreFSR,  _genM2ZpreFSR,
			       _genZpostFSR, _genM1ZpostFSR, _genM2ZpostFSR,
			       _genHpreFSR,  _genM1HpreFSR,  _genM2HpreFSR,
			       _genHpostFSR, _genM1HpostFSR, _genM2HpostFSR,
			       _genWpreFSR,  _genMWpreFSR,
			       _genWpostFSR, _genMWpostFSR );
    }    
  } // End conditional: if (_isMonteCarlo)



  // ============================
  // Store everything in an NTuple
  // ============================
  _outTree->Fill();
  return;  
} // End void UFDiMuonsAnalyzer::analyze

	
///////////////////////////////////////////////////////////
// BeginJob ==============================================
//////////////////////////////////////////////////////////

// Method called once each job just before starting event loop
// Set up TTrees where we save all of the info gathered in the analyzer
void UFDiMuonsAnalyzer::beginJob() {

  displaySelection();

  // Set up the _outTree branches
  _outTree->Branch("event",         &_eventInfo,        _eventInfo.getVarString()        );
  _outTree->Branch("vertices",      &_vertexInfo,       _vertexInfo.getVarString()       );
  _outTree->Branch("muons",         &_muonInfo,         _muonInfo.getVarString()         ); 
  _outTree->Branch("diMuons",       &_diMuonInfo,       _diMuonInfo.getVarString()       );
  _outTree->Branch("electrons",     &_electronInfo,     _electronInfo.getVarString()     ); 
  _outTree->Branch("taus",          &_tauInfo,          _tauInfo.getVarString()          ); 
  _outTree->Branch("met",           &_metInfo,          _metInfo.getVarString()          );
  _outTree->Branch("jets",          &_jetInfo,          _jetInfo.getVarString()          );
  _outTree->Branch("jets_JES_up",   &_jetInfo_JES_up,   _jetInfo_JES_up.getVarString()   );
  _outTree->Branch("jets_JES_down", &_jetInfo_JES_down, _jetInfo_JES_down.getVarString() );

  _outTree->Branch("hltPaths",      &_triggerNames);
  _outTree->Branch("btagNames",     &_btagNames   );
  _outTree->Branch("tauIDNames",    &_tauIDNames  );

  // MC information
  if (_isMonteCarlo) {

     // Off shell gamma block
    _outTree->Branch("genGpreFSR",   &_genGpreFSR,   _genGpreFSR.getVarString()   );
    _outTree->Branch("genM1GpreFSR", &_genM1GpreFSR, _genM1GpreFSR.getVarString() );
    _outTree->Branch("genM2GpreFSR", &_genM2GpreFSR, _genM2GpreFSR.getVarString() );

    _outTree->Branch("genGpostFSR",   &_genGpostFSR,   _genGpostFSR.getVarString()   );
    _outTree->Branch("genM1GpostFSR", &_genM1GpostFSR, _genM1GpostFSR.getVarString() );
    _outTree->Branch("genM2GpostFSR", &_genM2GpostFSR, _genM2GpostFSR.getVarString() );

    // Z block
    _outTree->Branch("genZpreFSR",   &_genZpreFSR,   _genZpreFSR.getVarString()   );
    _outTree->Branch("genM1ZpreFSR", &_genM1ZpreFSR, _genM1ZpreFSR.getVarString() );
    _outTree->Branch("genM2ZpreFSR", &_genM2ZpreFSR, _genM2ZpreFSR.getVarString() );

    _outTree->Branch("genZpostFSR",   &_genZpostFSR,   _genZpostFSR.getVarString()   );
    _outTree->Branch("genM1ZpostFSR", &_genM1ZpostFSR, _genM1ZpostFSR.getVarString() );
    _outTree->Branch("genM2ZpostFSR", &_genM2ZpostFSR, _genM2ZpostFSR.getVarString() );

    // W block
    _outTree->Branch("genWpreFSR",   &_genWpreFSR,   _genWpreFSR.getVarString()   );
    _outTree->Branch("genMWpreFSR",  &_genMWpreFSR,  _genMWpreFSR.getVarString()  );
    _outTree->Branch("genWpostFSR",  &_genWpostFSR,  _genWpostFSR.getVarString()  );
    _outTree->Branch("genMWpostFSR", &_genMWpostFSR, _genMWpostFSR.getVarString() );

    // H block
    _outTree->Branch("genHpreFSR",   &_genHpreFSR,   _genHpreFSR.getVarString()   );
    _outTree->Branch("genM1HpreFSR", &_genM1HpreFSR, _genM1HpreFSR.getVarString() );
    _outTree->Branch("genM2HpreFSR", &_genM2HpreFSR, _genM2HpreFSR.getVarString() );

    _outTree->Branch("genHpostFSR",   &_genHpostFSR,   _genHpostFSR.getVarString()   );
    _outTree->Branch("genM1HpostFSR", &_genM1HpostFSR, _genM1HpostFSR.getVarString() );
    _outTree->Branch("genM2HpostFSR", &_genM2HpostFSR, _genM2HpostFSR.getVarString() );

    _outTree->Branch("genJets", &_genJetInfo, _genJetInfo.getVarString() );

    _outTree->Branch("nPU", 	  &_nPU   	,"nPU/I"      );              
    _outTree->Branch("genWeight", &_genWeight   ,"genWeight/I");              
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

  _outTreeMetadata->Branch("triggerNames"  ,"std::vector< std::string > >", &_triggerNames);
  _outTreeMetadata->Branch("btagNames"     ,"std::vector< std::string > >", &_btagNames);
  _outTreeMetadata->Branch("tauIDNames"    ,"std::vector< std::string > >", &_tauIDNames);

  _outTreeMetadata->Fill();

}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::isHltPassed(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
				    const edm::Handle<edm::TriggerResults>& triggerResultsHandle,
                                    const std::vector<std::string> desiredTriggerNames) 
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

  const TriggerNames &triggerNames = iEvent.triggerNames(*triggerResultsHandle);

  const unsigned nTriggers = triggerResultsHandle->size();
  for (unsigned iTrigger = 0; iTrigger < nTriggers; ++iTrigger)
  {
    const string triggerName = triggerNames.triggerName(iTrigger);
    string triggerNameStripped = boost::regex_replace(triggerName,re,"",boost::match_default | boost::format_sed);
    for(std::vector<std::string>::const_iterator desiredTriggerName=desiredTriggerNames.begin();
            desiredTriggerName!=desiredTriggerNames.end();desiredTriggerName++)
    {
      if (*desiredTriggerName == triggerNameStripped && triggerResultsHandle->accept(iTrigger))
      {
        stringstream debugString;
        debugString << "isHltPassed:" <<endl;
        debugString << "  Trigger "<<iTrigger<<": "<< triggerName << "("<<triggerNameStripped<<") passed: "<<triggerResultsHandle->accept(iTrigger)<<endl;
        debugString << "    Desired Trigger Names: ";
        debugString <<"'"<< *desiredTriggerName<<"' ";
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

  std::cout << "_electron_ID      = " << _electron_ID << std::endl;
  std::cout << "_electron_pT_min  = " << _electron_pT_min << std::endl;
  std::cout << "_electron_eta_max = " << _electron_eta_max << std::endl;

  std::cout << "_tau_pT_min  = " << _tau_pT_min << std::endl;
  std::cout << "_tau_eta_max = " << _tau_eta_max << std::endl;

  std::cout << "_jet_ID      = " << _jet_ID << std::endl;
  std::cout << "_jet_pT_min  = " << _jet_pT_min << std::endl;
  std::cout << "_jet_eta_max = " << _jet_eta_max << std::endl;

  // module config parameters
  std::cout << " - _skim_trigger: " << _skim_trigger << std::endl;
  std::cout << " - Triggers To Probe:\n";
  unsigned int triggerSize = _triggerNames.size();
  for (unsigned int i=0; i < triggerSize; i++) 
    std::cout << "    * triggerNames["<<i<<"]: " << _triggerNames[i] << std::endl;
  
  std::cout << std::endl << std::endl;

}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::sortMuonsByPt    (pat::Muon i,     pat::Muon j    ) { return (i.pt() > j.pt()); }
bool UFDiMuonsAnalyzer::sortElectronsByPt(pat::Electron i, pat::Electron j) { return (i.pt() > j.pt()); }
bool UFDiMuonsAnalyzer::sortTausByPt     (pat::Tau i,      pat::Tau j     ) { return (i.pt() > j.pt()); }
bool UFDiMuonsAnalyzer::sortJetsByPt     (pat::Jet i,      pat::Jet j     ) { return (i.pt() > j.pt()); }
