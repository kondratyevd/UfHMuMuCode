///////////////////////////////////////////////////////////
//=========================================================
// UFDiMuonsAnalyzer.cxx                                 //
//=========================================================
//                                                       //
// v1 coded in 2010 by Gian Pierro.                      //
// Updated by various graduate students since.           //
// Analyzer to grab the info needed for HToMuMu          //
// and store it into a TTree.                            //
//                                                       //
//========================================================
///////////////////////////////////////////////////////////

#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/UFDiMuonsAnalyzer.h"


///////////////////////////////////////////////////////////
// Constructors/Destructors===============================
//////////////////////////////////////////////////////////

UFDiMuonsAnalyzer::UFDiMuonsAnalyzer(const edm::ParameterSet& iConfig):_numEvents(0)
{
  // Initialize the weighted count and the trees.
  // Use the file service to make the trees so that it knows to save them.
  _sumEventWeights = 0;
  _outTree = fs->make<TTree>("tree", "myTree");
  _outTreeMetadata = fs->make<TTree>("metadata", "Metadata Tree");

  // Get the collections designated from the python config file and load them into the tokens
  // muons
  _muonCollToken = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonColl"));

  // electrons
  _electronCollToken = consumes< edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronColl"));
  _electronCutBasedIdVetoToken = consumes< edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronCutBasedIdVeto"));
  _electronCutBasedIdLooseToken = consumes< edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronCutBasedIdLoose"));
  _electronCutBasedIdMediumToken = consumes< edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronCutBasedIdMedium"));
  _electronCutBasedIdTightToken = consumes< edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronCutBasedIdTight"));

  // taus
  _tauCollToken = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauColl"));
  _tauIDNames = iConfig.getParameter<std::vector<std::string> >("tauIDNames");

   // jets
  _metToken     = consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("metTag"));
  _pfJetsToken  = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("pfJetsTag"));
  _btagNames = iConfig.getParameter<std::vector<std::string> >("btagNames");

  // event stuff
  _beamSpotToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  _primaryVertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"));
  _PupInfoToken = consumes< std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));

  // gen info
  _genJetsToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsTag"));
  _prunedGenParticleToken = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag" ));
  _packedGenParticleToken = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticleTag" ));
  _genEvtInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

  // Get boolean switches from config file
  _isVerbose	= iConfig.getUntrackedParameter<bool>("isVerbose", false);
  _isMonteCarlo	= iConfig.getParameter<bool>("isMonteCarlo");
  _checkTrigger = iConfig.getParameter<bool>("checkTrigger");

  // Get selection criteria from config file
  _nMuons  = iConfig.getParameter<int>("nMuons");
  _isGlobal = iConfig.getParameter<int>("isGlobal");
  _isTracker = iConfig.getParameter<int>("isTracker");
  _ptMin  = iConfig.getParameter<double>("ptMin");
  _etaMax = iConfig.getParameter<double>("etaMax");

  if (!_isGlobal && !_isTracker) 
    std::cout << "\n\nWARNING: you are not requiring the muon to be not TRACKER nor GLOBAL\n" << "Be aware of the fact that StandAlone muon only are rejected anyway in the code\n";

  // Get trigger information from config file
  _processName = iConfig.getParameter<std::string>("processName");
  _triggerNames = iConfig.getParameter<std::vector <std::string> >("triggerNames");

  _triggerResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
  _triggerObjsToken = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjs"));

}

// Destructor
UFDiMuonsAnalyzer::~UFDiMuonsAnalyzer() {}

///////////////////////////////////////////////////////////
// Analyze ===============================================
//////////////////////////////////////////////////////////

// Need to break this up into functions for easier readability
// Muons, Jets, Vertices
// Could make some helper classes

void UFDiMuonsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
// function called every event

  // -----------------------------------------
  // COUNT EVENTS AND WEIGHTS 
  // -----------------------------------------
  _numEvents++;
  if (!_isMonteCarlo)
  {
      _sumEventWeights += 1;
  }
  else
  {
    // The generated weight. Due to the interference of terms in QM in the NLO simulations
    // there are negative weights that need to be accounted for. 
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(_genEvtInfoToken, genEvtInfo );
    _genWeight = (genEvtInfo->weight() > 0)? 1 : -1;
    _sumEventWeights += _genWeight;
  }

  if (_isVerbose) 
    std::cout << "\n\n A N A LI Z I N G   E V E N T = " 
	      << _numEvents << std::endl << std::endl;
 
  // -----------------------------------------
  // HLT HANDLES
  // -----------------------------------------
  
  iEvent.getByToken(_triggerResultsToken, _triggerResultsHandle);
  if (!_triggerResultsHandle.isValid()) 
  {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }

  iEvent.getByToken(_triggerObjsToken, _triggerObjsHandle);
  if (!_triggerObjsHandle.isValid()) 
  {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerObjects product from Event!" << std::endl;
    return;
  }

  // if all the HLT paths are not fired, will discard the event immediately
  if (_checkTrigger) 
  {
    if ( !isHltPassed(iEvent,iSetup,_triggerNames) )
    {
      if (_isVerbose) std::cout << "None of the HLT paths fired -> discard the event\n";
      return;
    }
     
  }

  // -----------------------------------------
  // RUN/EVENT INFO
  // -----------------------------------------
  
  _eventInfo.init();

  int theRun   = iEvent.id().run();
  int theLumi  = iEvent.luminosityBlock();
  long long int theEvent = iEvent.id().event();
  int theBx    = iEvent.bunchCrossing();
  int theOrbit = iEvent.orbitNumber();
  
  _eventInfo.run   = theRun;
  _eventInfo.lumi  = theLumi;
  _eventInfo.event = theEvent;
  _eventInfo.bx    = theBx;
  _eventInfo.orbit = theOrbit;


  // -----------------------------------------
  // VERTICES AND PILEUP
  // -----------------------------------------
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(_primaryVertexToken, vertices);

  // init vertices
  _vertexInfo.init();
  _vertexInfo.nVertices   = 0;

  if (vertices.isValid()) 
  {
    //std::cout << "vertex->size():"<< vertices->size() << std::endl;
    int iVertex=0;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) 
      {
        _vertexInfo.isValid[iVertex] = 0;
        continue;
      }
    
      _vertexInfo.isValid[iVertex] = 1;
      
      _vertexInfo.x[iVertex]        = vtx->position().X();	
      _vertexInfo.y[iVertex]        = vtx->position().Y();	
      _vertexInfo.z[iVertex]        = vtx->position().Z();	
      _vertexInfo.xErr[iVertex]     = vtx->xError();	
      _vertexInfo.yErr[iVertex]     = vtx->yError();	
      _vertexInfo.zErr[iVertex]     = vtx->zError();	
      _vertexInfo.chi2[iVertex]     = vtx->chi2();	
      _vertexInfo.ndf[iVertex]      = vtx->ndof();	
      _vertexInfo.normChi2[iVertex] = vtx->normalizedChi2();
      
      iVertex++;
      _vertexInfo.nVertices++;
    }
  }
  else std::cout << "VertexCollection is NOT valid -> vertex Info NOT filled!" << std::endl;
  
  // Get MC Truth Pileup and Generated Weights
  // addPileupInfo for miniAOD version 1 same for AOD
  // slimmedAddPileupInfo for miniAOD version 2
  _nPU = -1;
  if (_isMonteCarlo) 
  {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByToken(_PupInfoToken, PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) 
    {
       int BX = PVI->getBunchCrossing();

       if(BX == 0) 
       { 
          _nPU = PVI->getTrueNumInteractions();
          continue;
       }
    }

  }
  
  // B E A M S P O T
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(_beamSpotToken, beamSpotHandle);

  // -----------------------------------------
  // MONTE CARLO GEN INFO: MUONS, Gamma, H, W, Z
  // -----------------------------------------
  
  if (_isMonteCarlo) {

    // initialize Gamma to default values
    _genGpreFSR.init(); _genM1GpreFSR.init(); _genM2GpreFSR.init();
    _genGpostFSR.init();_genM1GpostFSR.init();_genM2GpostFSR.init();

    // initialize Z to default values
    _genZpreFSR.init(); _genM1ZpreFSR.init(); _genM2ZpreFSR.init();
    _genZpostFSR.init();_genM1ZpostFSR.init();_genM2ZpostFSR.init();

    // initialize H to default values
    _genHpreFSR.init(); _genM1HpreFSR.init(); _genM2HpreFSR.init();
    _genHpostFSR.init();_genM1HpostFSR.init();_genM2HpostFSR.init();

    // initialize W to default values
    _genWpreFSR.init(); _genMWpreFSR.init();
    _genWpostFSR.init();_genMWpostFSR.init();

    edm::Handle<reco::GenParticleCollection> prunedGenParticles;
    iEvent.getByToken(_prunedGenParticleToken, prunedGenParticles);

    //std::cout << "\n====================================\n"; 
    for (reco::GenParticle g: *prunedGenParticles) 
    {
        reco::GenParticle* gen = &g;
        fillBosonAndMuDaughters(gen); // looks for gamma, W, Z, H and the muons from them 
    } // loop over gen level

  }// end _isMonteCarlo

  // -----------------------------------------
  // ELECTRONS
  // -----------------------------------------
  _electronInfo.init();
  
  edm::Handle< edm::View<pat::Electron> > electrons;
  iEvent.getByToken(_electronCollToken, electrons);
  _electronInfo.nElectrons = 0;
  for (size_t i = 0; i < electrons->size(); ++i)
  {
    const edm::Ptr<pat::Electron> e = electrons->ptrAt(i);

    if( e->pt() < 10 ) // keep only electrons above 10 GeV
      continue;

    if((unsigned int)_electronInfo.nElectrons < _electronInfo.arraySize)
        fillElectron(_electronInfo.nElectrons, e, vertices, iEvent);
    
    _electronInfo.nElectrons++; 
  }
  // -----------------------------------------
  // TAUS
  // -----------------------------------------
  
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(_tauCollToken, taus);

  _tauInfo.init();

  _tauInfo.nTaus = 0;
  for (pat::TauCollection::const_iterator tau = taus->begin(), tausEnd = taus->end(); tau !=tausEnd; ++tau)
  {
    if( tau->pt() < 20) // keep only taus above 20 GeV
        continue;

    if((unsigned int)_tauInfo.nTaus < _tauInfo.arraySize)
        fillTau(_tauInfo.nTaus, *tau, vertices, iEvent);

    _tauInfo.nTaus++;
  }

  // -----------------------------------------
  // MET, JETS, GENJETS
  // -----------------------------------------
  
  edm::Handle < std::vector<pat::MET> > mets;
  if(!_metToken.isUninitialized()) iEvent.getByToken(_metToken, mets);

  // init met, set values to -999
  _metInfo.init();

  // fill actual met values
  if( mets.isValid() ){
    _metInfo.px = (*mets)[0].px();
    _metInfo.py = (*mets)[0].py();
    _metInfo.pt = (*mets)[0].pt();
    _metInfo.phi= (*mets)[0].phi();
    _metInfo.sumEt = (*mets)[0].sumEt();
  }

  edm::Handle < std::vector<pat::Jet> > jets;
  if(!_pfJetsToken.isUninitialized()) iEvent.getByToken(_pfJetsToken, jets);

  // init pfjets, set values to -999
  _pfJetInfo.init();

  //// Get JEC Uncertainty Calculator
  //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl); 
  //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  //JetCorrectionUncertainty *jecUncCalculator = new JetCorrectionUncertainty(JetCorPar);

  // fill actual jet values
  if( jets.isValid() )
  {
    _pfJetInfo.nJets = 0;

    for(unsigned int i=0; i<jets->size(); i++)
    {
      _pfJetInfo.nJets++;
      if( i<_pfJetInfo.arraySize )
      {
        const pat::Jet& jet = jets->at(i);
        _pfJetInfo.px[i] = jet.px();
        _pfJetInfo.py[i] = jet.py();
        _pfJetInfo.pz[i] = jet.pz();
        _pfJetInfo.pt[i] = jet.pt();
        _pfJetInfo.eta[i]= jet.eta();
        _pfJetInfo.phi[i]= jet.phi();
        _pfJetInfo.mass[i]  = jet.mass();
        _pfJetInfo.charge[i]  = jet.charge();
        _pfJetInfo.partonFlavour[i] = jet.partonFlavour();
        // Energy Fractions
        _pfJetInfo.chf[i]  = jet.chargedHadronEnergyFraction();
        _pfJetInfo.nhf[i]  = jet.neutralHadronEnergyFraction();
        _pfJetInfo.cef[i]  = jet.chargedEmEnergyFraction();
        _pfJetInfo.nef[i]  = jet.neutralEmEnergyFraction();
        _pfJetInfo.muf[i]  = jet.muonEnergyFraction();
        _pfJetInfo.hfhf[i]  = jet.HFHadronEnergyFraction();
        _pfJetInfo.hfef[i]  = jet.HFEMEnergyFraction();
        // Multiplicities
        _pfJetInfo.cm[i]  = jet.chargedMultiplicity();
        _pfJetInfo.chm[i]  = jet.chargedHadronMultiplicity();
        _pfJetInfo.nhm[i]  = jet.neutralHadronMultiplicity();
        _pfJetInfo.cem[i]  = jet.electronMultiplicity();
        _pfJetInfo.nem[i]  = jet.photonMultiplicity();
        _pfJetInfo.mum[i]  = jet.muonMultiplicity();
        _pfJetInfo.hfhm[i]  = jet.HFHadronMultiplicity();
        _pfJetInfo.hfem[i]  = jet.HFEMMultiplicity();
        //               _pfJetInfo.pfJetCh[i] = jet.jetCharge();
        // Get JEC Uncertainty
        //jecUncCalculator->setJetEta(jet.eta());
        //jecUncCalculator->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
        //_pfJetInfo.jecUnc[i] = jecUncCalculator->getUncertainty(true);
        _pfJetInfo.jecUnc[i] = -1.0;
        _pfJetInfo.jecFactor[i]  = jet.jecFactor("Uncorrected");
        // b-Tag
        _pfJetInfo.isB[i]  = jet.bDiscriminator(_btagNames[0]);
        _pfJetInfo.puid[i]  = jet.userFloat("pileupJetId:fullDiscriminant");
        //PAT matched Generator Jet
        const reco::GenJet* genJet = jet.genJet();
        if (genJet != NULL)
        {
            _pfJetInfo.genMatched[i] = true;
            _pfJetInfo.genPx[i] = genJet->px();
            _pfJetInfo.genPy[i] = genJet->py();
            _pfJetInfo.genPz[i] = genJet->pz();
            _pfJetInfo.genPt[i] = genJet->pt();
            _pfJetInfo.genEta[i]= genJet->eta();
            _pfJetInfo.genPhi[i]= genJet->phi();
            _pfJetInfo.genMass[i]  = genJet->mass();
            double genJetEnergy = genJet->energy();
            _pfJetInfo.genEMF[i]  = genJet->emEnergy()/genJetEnergy;
            _pfJetInfo.genHadF[i]  = genJet->hadEnergy()/genJetEnergy;
            _pfJetInfo.genInvF[i]  = genJet->invisibleEnergy()/genJetEnergy;
            _pfJetInfo.genAuxF[i]  = genJet->auxiliaryEnergy()/genJetEnergy;

        }
        else
        {
            _pfJetInfo.genMatched[i] = false;
            _pfJetInfo.genPx[i] =-1;
            _pfJetInfo.genPy[i] =-1;
            _pfJetInfo.genPz[i] =-1;
            _pfJetInfo.genPt[i] =-1;
            _pfJetInfo.genEta[i]=-1;
            _pfJetInfo.genPhi[i]=-1;
            _pfJetInfo.genMass[i]  =-1;
            _pfJetInfo.genEMF[i]  =-1;
            _pfJetInfo.genHadF[i]  =-1;
            _pfJetInfo.genInvF[i]  =-1;
            _pfJetInfo.genAuxF[i]  =-1;
        }

        //delete genJet;	  
      
      }
    }
  }
  
//  delete jecUncCalculator;

  edm::Handle < reco::GenJetCollection > genJets;
  if(!_genJetsToken.isUninitialized()) iEvent.getByToken(_genJetsToken, genJets);

  _genJetInfo.init();

  if( genJets.isValid() ){
    reco::GenJetCollection sortedGenJets = (*genJets);
    sort(sortedGenJets.begin(), sortedGenJets.end(), sortGenJetFunc);
    for(unsigned int i=0; i<sortedGenJets.size(); i++){
      _genJetInfo.nJets++;
      if( i<_genJetInfo.arraySize ){
        _genJetInfo.px[i] = sortedGenJets[i].px();
        _genJetInfo.py[i] = sortedGenJets[i].py();
        _genJetInfo.pz[i] = sortedGenJets[i].pz();
        _genJetInfo.pt[i] = sortedGenJets[i].pt();
        _genJetInfo.eta[i] = sortedGenJets[i].eta();
        _genJetInfo.phi[i] = sortedGenJets[i].phi();
        _genJetInfo.mass[i] = sortedGenJets[i].mass();
      }
    }
  }

  // -----------------------------------------
  // MUONS
  // -----------------------------------------
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(_muonCollToken, muons);

  // reco muons collection 
  pat::MuonCollection muonsSelected; 
  muonsSelected.clear(); 

  int countMuons = 0;

  // pre-selection: just check the muons are at least tracker muons.
  for (pat::MuonCollection::const_iterator muon = muons->begin(), muonsEnd = muons->end(); muon !=muonsEnd; ++muon)
  {
    countMuons++;

    // basic cuts
    if (!passBasicMuonSelection(*muon)) {
      if (_isVerbose) std::cout << "Muon NOT passing basic selections\n"; 
      continue;
    }
    if (_isVerbose) std::cout << "Muon passing the basic selections\n\n";     
    
    // put this muons in the collection
    muonsSelected.push_back(*muon);

  }
  // ===========================================================================
 
  if (_isVerbose) std::cout << " Found " << muonsSelected.size() 
			    << " muons candidates\n\n";

  // at least _nMuons muons
  if ( muonsSelected.size() < _nMuons  ) return;

  // -----------------------------------------
  // INIT MUONS & DIMUON CANDIDATE
  // -----------------------------------------
  
  _muonInfo.init();
  _dimuCandInfo.init();

  _muonInfo.nMuons = muonsSelected.size();
  _muonInfo.nMuonPairs = 0;

  // Zero Muons /////////////////////////////////////////////
  if (muonsSelected.size() == 0) {
    if (_isVerbose) std::cout << "0 reco'd muons...\n";
    _outTree->Fill();
    return;
  }

  // One Muon /////////////////////////////////////////////
  if (muonsSelected.size() == 1) {
    if (_isVerbose) std::cout << "Only 1 reco'd muon...\n";
    pat::Muon mu = muonsSelected.at(0);
    // store all the info
    // muon 1

    fillMuon(0, mu, vertices, beamSpotHandle, iEvent, iSetup);
    _outTree->Fill();
    
    // you can exit and pass to the new event
    return;
  }
  
  /////////////////////////////////////////////////////////////////// 
  // More than one muon /////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////// 
  
  // ... if we have two or more muons we construct dimuon candidates
  
  // grab the candidates: they are already sorted by
  // distance to the Z mass PDG value: the closer the better
  MuonPairs dimuons = GetMuonPairs(&muonsSelected);

  // sort the selected muons by pt
  pat::MuonCollection sortedMuons = muonsSelected;
  sort(sortedMuons.begin(), sortedMuons.end(), sortMuonFunc);

  // loop over the candidates
  for (MuonPairs::const_iterator pair= dimuons.begin(); pair != dimuons.end(); ++pair)
  {
    
    pat::Muon mu1  = pair->first;
    pat::Muon mu2  = pair->second;

    // only consider dimuon events with opposite charge
    if(mu1.charge() == mu2.charge()) continue;
 
    _muonInfo.nMuonPairs++;
    fillDimuonCandidate(&*pair, vertices, beamSpotHandle, iEvent, iSetup);
    // fill muons besides the reco candidate pair, looking at the other muons in sortedMuons = sort(muonsSelected, byPt)
    // muons are only considered if they pass the basic selection criteria: at least tracker, pt > ptMin, eta < etaMax
    fillOtherMuons(&*pair, sortedMuons, vertices, beamSpotHandle, iEvent, iSetup);

    // ===========================================================================
    // store everything in a ntuple
    _outTree->Fill();
  
  }
  
  return;
  
}
	
///////////////////////////////////////////////////////////
// BeginJob ==============================================
//////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::beginJob()
{
// Method called once each job just before starting event loop
// Set up TTrees where we save all of the info gathered in the analyzer

  displaySelection();

  // eventInfo;
  // set up the _outTree branches
  _outTree->Branch("eventInfo",     &_eventInfo,    _eventInfo.getVarString());
  _outTree->Branch("vertexInfo",    &_vertexInfo,   _vertexInfo.getVarString());
  _outTree->Branch("recoMuons",     &_muonInfo,     _muonInfo.getVarString()); 
  _outTree->Branch("dimuCand",      &_dimuCandInfo, _dimuCandInfo.getVarString());
  _outTree->Branch("recoElectrons", &_electronInfo, _electronInfo.getVarString()); 
  _outTree->Branch("recoTaus",      &_tauInfo,      _tauInfo.getVarString()); 
  _outTree->Branch("met",           &_metInfo,      _metInfo.getVarString());
  _outTree->Branch("pfJets",        &_pfJetInfo,    _pfJetInfo.getVarString());

  _outTree->Branch("hltPaths",      &_triggerNames);
  _outTree->Branch("btagNames",     &_btagNames);
  _outTree->Branch("tauIDNames",    &_tauIDNames);

  // MC information
  if (_isMonteCarlo) {

     // Off shell gamma block
    _outTree->Branch("genGpreFSR",  &_genGpreFSR,  _genGpreFSR.getVarString());
    _outTree->Branch("genM1GpreFSR",&_genM1GpreFSR,_genM1GpreFSR.getVarString());
    _outTree->Branch("genM2GpreFSR",&_genM2GpreFSR,_genM2GpreFSR.getVarString());

    _outTree->Branch("genGpostFSR",  &_genGpostFSR, _genGpostFSR.getVarString());
    _outTree->Branch("genM1GpostFSR",&_genM1GpostFSR, _genM1GpostFSR.getVarString());
    _outTree->Branch("genM2GpostFSR",&_genM2GpostFSR, _genM2GpostFSR.getVarString());

    // Z block
    _outTree->Branch("genZpreFSR",  &_genZpreFSR, _genZpreFSR.getVarString());
    _outTree->Branch("genM1ZpreFSR",&_genM1ZpreFSR, _genM1ZpreFSR.getVarString());
    _outTree->Branch("genM2ZpreFSR",&_genM2ZpreFSR, _genM2ZpreFSR.getVarString());

    _outTree->Branch("genZpostFSR",  &_genZpostFSR, _genZpostFSR.getVarString());
    _outTree->Branch("genM1ZpostFSR",&_genM1ZpostFSR, _genM1ZpostFSR.getVarString());
    _outTree->Branch("genM2ZpostFSR",&_genM2ZpostFSR, _genM2ZpostFSR.getVarString());

    // W block
    _outTree->Branch("genWpreFSR",  &_genWpreFSR, _genWpreFSR.getVarString());
    _outTree->Branch("genMWpreFSR", &_genMWpreFSR , _genMWpreFSR.getVarString());
    _outTree->Branch("genWpostFSR",  &_genWpostFSR, _genWpostFSR.getVarString());
    _outTree->Branch("genMWpostFSR",&_genMWpostFSR, _genMWpostFSR.getVarString());

    // H block
    _outTree->Branch("genHpreFSR",  &_genHpreFSR, _genHpreFSR.getVarString());
    _outTree->Branch("genM1HpreFSR",&_genM1HpreFSR, _genM1HpreFSR.getVarString());
    _outTree->Branch("genM2HpreFSR",&_genM2HpreFSR, _genM2HpreFSR.getVarString());

    _outTree->Branch("genHpostFSR",  &_genHpostFSR, _genHpostFSR.getVarString());
    _outTree->Branch("genM1HpostFSR",&_genM1HpostFSR, _genM1HpostFSR.getVarString());
    _outTree->Branch("genM2HpostFSR",&_genM2HpostFSR, _genM2HpostFSR.getVarString());

    _outTree->Branch("genJets", &_genJetInfo, _genJetInfo.getVarString());

    _outTree->Branch("nPU", 	  &_nPU   	,"nPU/I");              
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


TLorentzVector const UFDiMuonsAnalyzer::GetLorentzVector(UFDiMuonsAnalyzer::MuonPair const* pair)
{
  
  TLorentzVector muon1, muon2, sum;
  double const MASS_MUON = 0.105658367;    //GeV/c2

  //reco::TrackRef const muon1Track = pair->first . innerTrack();
  //reco::TrackRef const muon2Track = pair->second. innerTrack();

  //muon1.SetPtEtaPhiM(muon1Track->pt(), muon1Track->eta(), muon1Track->phi(), MASS_MUON);
  //muon2.SetPtEtaPhiM(muon2Track->pt(), muon2Track->eta(), muon2Track->phi(), MASS_MUON);

  reco::Track const muon1Track = *(pair->first . innerTrack());
  reco::Track const muon2Track = *(pair->second. innerTrack());

  muon1.SetPtEtaPhiM(muon1Track.pt(), muon1Track.eta(), muon1Track.phi(), MASS_MUON);
  muon2.SetPtEtaPhiM(muon2Track.pt(), muon2Track.eta(), muon2Track.phi(), MASS_MUON);

  sum = muon1+muon2;
  return sum;

}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

TLorentzVector const UFDiMuonsAnalyzer::GetLorentzVector(UFDiMuonsAnalyzer::TrackPair const* pair) {
  
  TLorentzVector muon1, muon2, sum;
  double const MASS_MUON = 0.105658367;    //GeV/c2
  
  reco::Track const muon1Track = pair->first;
  reco::Track const muon2Track = pair->second;

  muon1.SetPtEtaPhiM(muon1Track.pt(), muon1Track.eta(), muon1Track.phi(), MASS_MUON);
  muon2.SetPtEtaPhiM(muon2Track.pt(), muon2Track.eta(), muon2Track.phi(), MASS_MUON);

  sum = muon1+muon2;
  return sum;
  
}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::isHltPassed(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
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

  const TriggerNames &triggerNames = iEvent.triggerNames(*_triggerResultsHandle);

  const unsigned nTriggers = _triggerResultsHandle->size();
  for (unsigned iTrigger = 0; iTrigger < nTriggers; ++iTrigger)
  {
    const string triggerName = triggerNames.triggerName(iTrigger);
    string triggerNameStripped = boost::regex_replace(triggerName,re,"",boost::match_default | boost::format_sed);
    for(std::vector<std::string>::const_iterator desiredTriggerName=desiredTriggerNames.begin();
            desiredTriggerName!=desiredTriggerNames.end();desiredTriggerName++)
    {
      if (*desiredTriggerName == triggerNameStripped && _triggerResultsHandle->accept(iTrigger))
      {
        stringstream debugString;
        debugString << "isHltPassed:" <<endl;
        debugString << "  Trigger "<<iTrigger<<": "<< triggerName << "("<<triggerNameStripped<<") passed: "<<_triggerResultsHandle->accept(iTrigger)<<endl;
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
//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::isHltMatched(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& desiredTriggerName, 
                                     const pat::TriggerObjectStandAloneCollection& triggerObjects, const pat::Muon& mu)
{
// same check for isHltPassed +
// check if the muon is the one firing the HLT path
  using namespace std;
  using namespace edm;
  using namespace pat;
  using namespace reco;
  using namespace trigger;


  const boost::regex re("_v[0-9]+");

  const TriggerNames &triggerNames = iEvent.triggerNames(*_triggerResultsHandle);

  const unsigned nTriggers = _triggerResultsHandle->size();
  for (unsigned iTrigger = 0; iTrigger < nTriggers; ++iTrigger)
  {
    const string triggerName = triggerNames.triggerName(iTrigger);
    string triggerNameStripped = boost::regex_replace(triggerName,re,"",boost::match_default | boost::format_sed);
    if (desiredTriggerName == triggerNameStripped && _triggerResultsHandle->accept(iTrigger))
    {
      stringstream debugString;
      debugString << "isHltMatched: ";
      debugString << "'" << desiredTriggerName<<"'\n";
      debugString << "  Trigger "<<iTrigger<<": "<< triggerName << "("<<triggerNameStripped<<") passed: "<<_triggerResultsHandle->accept(iTrigger)<<endl;
      for(TriggerObjectStandAloneCollection::const_iterator trigObj=triggerObjects.begin(); 
          trigObj!=triggerObjects.end();trigObj++)
      {
        TriggerObjectStandAlone tmpTrigObj(*trigObj); // Get rid of const which messes up unpackPathNames
        tmpTrigObj.unpackPathNames(triggerNames);
        bool isRightObj = tmpTrigObj.hasPathName(triggerName,true,true); // name, check that is l3 filter accepted, check that is last filter
        if (isRightObj)
        {
          debugString << "    TriggerObject:  "<<tmpTrigObj.collection() << endl;
          bool isMatched = (deltaR(tmpTrigObj,mu) < 0.2);
          if (isMatched) 
          {
            debugString << "      is Matched*****"  <<endl;
            LogVerbatim("UFHLTTests") << debugString.str();
            return true;
          }
        } // if isRightObj
      } // trigObj lookp
      LogVerbatim("UFHLTTests") << debugString.str();
    }// if desiredTrigger
  }// iTrigger loop

  return false;
}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::passBasicMuonSelection(const pat::Muon& muon)
{
// Most basic cuts on the muons tested here: global mu, tracker mu, pt, eta
  
  bool passKinCuts=false;
  
  // =========================================================== //
  // What was corresponding to the old Loose VBTF
  // =========================================================== //
  if (_isVerbose) 
  {
    std::cout<< "is Global?"     << muon.isGlobalMuon()     << std::endl;
    std::cout<< "is Tracker?"    << muon.isTrackerMuon()    << std::endl;
    std::cout<< "is StandAlone?" << muon.isStandAloneMuon() << std::endl;
  }

  // reconstruction cuts
  if (!muon.isGlobalMuon() && _isGlobal ) return passKinCuts; // gbl muon
  if (!muon.isTrackerMuon() && _isTracker) return passKinCuts; // trk muon

  // do not accept muons which are standalone only
  if(!muon.isGlobalMuon() && !muon.isTrackerMuon()) return passKinCuts;

  reco::Track globalTrack;
  if (muon.isGlobalMuon())       globalTrack = *(muon.globalTrack());
  else if (muon.isTrackerMuon()) globalTrack = *(muon.innerTrack());
  else 
  {
    // redundant: just in case...
    std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
    return false;
  }

  if (_isVerbose) 
  {
    std::cout<< "muon.pt(): "  << muon.pt() << " [ptMin="  << _ptMin  <<"]" << std::endl;
    std::cout<< "fabs(muon.eta()): " << fabs(muon.eta())<< " [etaMax=" << _etaMax <<"]" << std::endl;
  }
    
  // kinematic cuts
  if (muon.pt()        < _ptMin ) return passKinCuts; // pt cut
  if (fabs(muon.eta()) > _etaMax) return passKinCuts; // eta cut

  if (_isVerbose) std::cout << "passing kinematic cuts\n"; 

  passKinCuts=true;
  return passKinCuts;
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

UFDiMuonsAnalyzer::MuonPairs const UFDiMuonsAnalyzer::GetMuonPairs(pat::MuonCollection const* muons) const 
{
                                                                           
  
  MuonPairs muonpairs; 
  muonpairs.clear();

  for (pat::MuonCollection::const_iterator muon1 = muons->begin(), 
         muonsEnd = muons->end(); muon1 != muonsEnd; ++muon1){

    for (pat::MuonCollection::const_iterator muon2 = muon1+1; 
         muon2 != muonsEnd; ++muon2){
      muonpairs.push_back( UFDiMuonsAnalyzer::MuonPair(*muon1,*muon2) );
    }
  }
  
  std::sort(muonpairs.begin(),muonpairs.end(),sortMuonObject);

  return muonpairs;
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillMuon(unsigned int i, const pat::Muon& mu, const edm::Handle<reco::VertexCollection>& vertices, const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                                 const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    _muonInfo.isGlobal[i]     = mu.isGlobalMuon(); 
    _muonInfo.isTracker[i]    = mu.isTrackerMuon(); 
    _muonInfo.isStandAlone[i] = mu.isStandAloneMuon(); 

    reco::Track track;
    if      (mu.isGlobalMuon())  track = *(mu.globalTrack());
    else if (mu.isTrackerMuon()) track = *(mu.innerTrack());
    else 
    {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    if(i >= _muonInfo.arraySize)
    {
      std::cout << "ERROR: Tried to add muon info to an index that is out of bounds.";
      return;
    }

    _muonInfo.charge[i] = mu.charge();
    _muonInfo.pt[i]     = mu.pt();  
    _muonInfo.ptErr[i]  = track.ptError(); 
    _muonInfo.eta[i]    = mu.eta(); 
    _muonInfo.phi[i]    = mu.phi();

    // redundant if the muon is tracker-only muon
    if (mu.isTrackerMuon()) {
      _muonInfo.trkPt[i]   = mu.innerTrack()->pt();                    
      _muonInfo.trkPtErr[i] = mu.innerTrack()->ptError();
      _muonInfo.trketa[i]  = mu.innerTrack()->eta();                   
      _muonInfo.trkPhi[i]  = mu.innerTrack()->phi();                   
    }

    _muonInfo.d0_BS[i]= mu.innerTrack()->dxy( beamSpotHandle->position() );
    _muonInfo.dz_BS[i]= mu.innerTrack()->dz ( beamSpotHandle->position() );

    reco::Vertex bestVtx1;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx!=vertices->end(); ++vtx)
    {
      if (!vtx->isValid()) continue;
  
      _muonInfo.d0_PV[i]= track.dxy( vtx->position() );
      _muonInfo.dz_PV[i]= track.dz ( vtx->position() );
      bestVtx1 = *vtx;
    
      //exit at the first available vertex
      break;
    }

    // type of muon
    _muonInfo.isTightMuon[i]  =  muon::isTightMuon(mu, bestVtx1);
    _muonInfo.isMediumMuon[i] =  muon::isMediumMuon(mu);
    _muonInfo.isLooseMuon[i]  =  muon::isLooseMuon(mu);

    //isolation
    _muonInfo.trackIsoSumPt[i]    = mu.isolationR03().sumPt;
    _muonInfo.trackIsoSumPtCorr[i]= mu.isolationR03().sumPt; // no correction with only 1 muon
    _muonInfo.ecalIso[i]          = mu.isolationR03().emEt;
    _muonInfo.hcalIso[i]          = mu.isolationR03().hadEt;

    double isovar = mu.isolationR03().sumPt;
    isovar += mu.isolationR03().hadEt; //tracker + HCAL 
    isovar /= mu.pt(); // relative combine isolation
    _muonInfo.relCombIso[i]=isovar;

  
    // PF Isolation
    _muonInfo.isPFMuon[i] = mu.isPFMuon();
    
    if ( mu.isPFMuon() ) {
      
      reco::Candidate::LorentzVector pfmuon = mu.pfP4();
      
      _muonInfo.pfPt[i]  = pfmuon.Pt();
      _muonInfo.pfEta[i] = pfmuon.Eta();
      _muonInfo.pfPhi[i] = pfmuon.Phi();

      _muonInfo.sumChargedHadronPtR03[i]   = mu.pfIsolationR03().sumChargedHadronPt  ;
      _muonInfo.sumChargedParticlePtR03[i] = mu.pfIsolationR03().sumChargedParticlePt;
      _muonInfo.sumNeutralHadronEtR03[i]   = mu.pfIsolationR03().sumNeutralHadronEt  ;
      _muonInfo.sumPhotonEtR03[i]          = mu.pfIsolationR03().sumPhotonEt         ;
      _muonInfo.sumPUPtR03[i]              = mu.pfIsolationR03().sumPUPt             ;
      
      _muonInfo.sumChargedHadronPtR04[i]   = mu.pfIsolationR04().sumChargedHadronPt  ;
      _muonInfo.sumChargedParticlePtR04[i] = mu.pfIsolationR04().sumChargedParticlePt;
      _muonInfo.sumNeutralHadronEtR04[i]   = mu.pfIsolationR04().sumNeutralHadronEt  ;
      _muonInfo.sumPhotonEtR04[i]          = mu.pfIsolationR04().sumPhotonEt         ;
      _muonInfo.sumPUPtR04[i]              = mu.pfIsolationR04().sumPUPt             ;
    }

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<_triggerNames.size();iTrigger++) 
      _muonInfo.isHltMatched[i][iTrigger] = isHltMatched(iEvent, iSetup, _triggerNames[iTrigger], *_triggerObjsHandle, mu);

    //std::cout << std::endl;
    //std::cout << "!!! muonInfo[" << i << "] ..." << std::endl;
    //std::cout << "charge: " << _muonInfo.charge[i] << std::endl;
    //std::cout << "pt: " << _muonInfo.pt[i] << std::endl;
    //std::cout << "eta: " << _muonInfo.eta[i] << std::endl;
    //std::cout << "phi: " << _muonInfo.phi[i] << std::endl;
    //std::cout << std::endl;
    //std::cout << "!!! pat::mu..." << std::endl;
    //std::cout << "charge: " << mu.charge() << std::endl;
    //std::cout << "pt: " << mu.pt() << std::endl;
    //std::cout << "eta: " << mu.eta() << std::endl;
    //std::cout << "phi: " << mu.phi() << std::endl;
    //std::cout << std::endl;
   
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillDimuonCandidate(const UFDiMuonsAnalyzer::MuonPair* pair, const edm::Handle<reco::VertexCollection>& vertices, 
                                            const edm::Handle<reco::BeamSpot>& beamSpotHandle, const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    pat::Muon mu1 = pair->first;
    pat::Muon mu2 = pair->second;

    fillMuon(0, mu1, vertices, beamSpotHandle, iEvent, iSetup);
    fillMuon(1, mu2, vertices, beamSpotHandle, iEvent, iSetup);

    if( mu1.innerTrack().isNonnull() ){
      double dEta =      mu1.track()->eta() - mu2.track()->eta();
      double dPhi = fabs(mu1.track()->phi() - mu2.track()->phi());
      if( dPhi>3.1415927 ) dPhi = 2.*3.1415927 - dPhi;
      if( sqrt(dEta*dEta+dPhi*dPhi)<0.3 && _muonInfo.trackIsoSumPt[1]>0.9*mu1.track()->pt() ) 
        _muonInfo.trackIsoSumPtCorr[0] = _muonInfo.trackIsoSumPt[1] - mu1.track()->pt();
    }

    if( mu2.innerTrack().isNonnull() ){
      double dEta =      mu2.track()->eta() - mu1.track()->eta();
      double dPhi = fabs(mu2.track()->phi() - mu1.track()->phi());
      if( dPhi>3.1415927 ) dPhi = 2.*3.1415927 - dPhi;
      if( sqrt(dEta*dEta+dPhi*dPhi)<0.3 && _muonInfo.trackIsoSumPt[0]>0.9*mu2.track()->pt() ) 
        _muonInfo.trackIsoSumPtCorr[1] = _muonInfo.trackIsoSumPt[0] - mu2.track()->pt();
    }

    // combine the info for the dimuon candidate
    TLorentzVector mother=GetLorentzVector(pair); 

    _dimuCandInfo.recoCandMass = mother.M();
    _dimuCandInfo.recoCandPt   = mother.Pt();
    _dimuCandInfo.recoCandEta  = mother.PseudoRapidity();
    _dimuCandInfo.recoCandY    = mother.Rapidity();
    _dimuCandInfo.recoCandPhi  = mother.Phi();

    _dimuCandInfo.angleDiMuons = acos(-mu1.track()->momentum().Dot(mu2.track()->momentum()/
                                                      mu1.track()->p()/mu2.track()->p())); 
    
    
    if ( mu1.isPFMuon() && mu2.isPFMuon() ) {

      TLorentzVector muon1_pf, muon2_pf, mother_pf;
      double const MASS_MUON = 0.105658367; //GeV/c2

      reco::Candidate::LorentzVector pfmuon1 = mu1.pfP4();
      reco::Candidate::LorentzVector pfmuon2 = mu2.pfP4();
      
      muon1_pf.SetPtEtaPhiM(pfmuon1.Pt(), pfmuon1.Eta(), pfmuon1.Phi(), MASS_MUON);
      muon2_pf.SetPtEtaPhiM(pfmuon2.Pt(), pfmuon2.Eta(), pfmuon2.Phi(), MASS_MUON);

      mother_pf = muon1_pf+muon2_pf;
      
      _dimuCandInfo.recoCandMassPF = mother_pf.M();
      _dimuCandInfo.recoCandPtPF   = mother_pf.Pt();
                              
      _dimuCandInfo.recoCandEtaPF  = mother_pf.PseudoRapidity();
      _dimuCandInfo.recoCandYPF    = mother_pf.Rapidity();
      _dimuCandInfo.recoCandPhiPF  = mother_pf.Phi();
    
    }
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillOtherMuons(const UFDiMuonsAnalyzer::MuonPair* pair, const pat::MuonCollection& sortedMuons, 
                                       const edm::Handle<reco::VertexCollection>& vertices, const edm::Handle<reco::BeamSpot>& beamSpotHandle, 
                                       const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
// Fill the muons that aren't the reco dimuon candidate into the muon info spots other than [0] and [1]
// Other muons are sorted by pt

  pat::Muon mu1 = pair->first;
  pat::Muon mu2 = pair->second;

  // dimuon candidate muons occupy muonInfo[0] and muonInfo[1]
  // start adding the other muons at muonInfo[2]
  unsigned int i=2;

  // pre-selection: just check the muons are at least tracker muons.
  for (pat::MuonCollection::const_iterator muon = sortedMuons.begin(), muonsEnd = sortedMuons.end(); muon !=muonsEnd; ++muon)
  {
    // muon is the same as mu1
    if (mu1.pt() == (*muon).pt() && mu1.phi() == (*muon).phi() && mu1.eta() == (*muon).eta())
        continue;

    // muon is the same as mu2
    if (mu2.pt() == (*muon).pt() && mu2.phi() == (*muon).phi() && mu2.eta() == (*muon).eta())
        continue;

    // put this muon in the collection if the collection is not full and the muon is at least a tracker muon
    if(i<_muonInfo.arraySize && (*muon).isTrackerMuon()) 
    {
        fillMuon(i, *muon, vertices, beamSpotHandle, iEvent, iSetup);
        i++;
    }
    else return;
  }
}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillBosonAndMuDaughters(const reco::Candidate* boson)
{
  // photon, Z, W, H
  if(boson->status() == 62)
  {
      if(!(abs(boson->pdgId()) == 22 || abs(boson->pdgId()) == 23 || abs(boson->pdgId()) == 24 || abs(boson->pdgId()) == 25)) return;
  }
  // technically 21 is an incoming particle from the feynman diagram
  // q, anti-q -> lept,lept but no intermediate gen particle. I think this is q,anti-q->gamma*->lept,lept
  // whenever there is no Z in the dyJetsToLL sample there is this q,antiq -> lept, lept where the quark 
  // and the antiquark have the same two leptons as daughters
  // we will have to reconstruct the values for the off shell gamma
  else if(boson->status() == 21) ;

  // Don't care about other situations
  else return;

  // initialize the temporary structs
  _genPartInfo bosonInfo;
  _TrackInfo mu1preFSR;
  _TrackInfo mu1postFSR;
  _TrackInfo mu2preFSR;
  _TrackInfo mu2postFSR;

  bosonInfo.init();
  mu1preFSR.init();
  mu1postFSR.init();
  mu2preFSR.init();
  mu2postFSR.init();

  bosonInfo.mass = boson->mass(); 
  bosonInfo.pt   = boson->pt();   
  bosonInfo.eta  = boson->eta();  
  bosonInfo.y    = boson->rapidity();    
  bosonInfo.phi  = boson->phi();  
  bosonInfo.charge = boson->charge();

  TLorentzVector l1, l2, mother;
  bool moreThanOneLeptPair = false;

  // check mass of particles for daughters to see if muon or not

  // Get the daughter muons for the boson 
  for(unsigned int i=0; i<boson->numberOfDaughters(); i++)
  {
      const reco::Candidate* daughter = boson->daughter(i);

      // get information about the lepton daughters to reconstruct the virtual photon later
      if(daughter->pdgId() == 11 || daughter->pdgId() == 13 || daughter->pdgId() == 15)
      {
          for(unsigned int j=0; j<boson->numberOfDaughters(); j++)
          {
              // we already know you can't make a lepton pair with yourself, so skip this one if it's the case
              if(i==j) continue;
              const reco::Candidate* daughter2 = boson->daughter(j);

              // found a pair of opposite signed lepton daughters
              if(daughter->pdgId() == -1*daughter2->pdgId()) 
              {
                  // we already found a pair, l1 AND l2 were initialized already
                  if(l1.M() !=0 && l2.M() != 0) moreThanOneLeptPair = true;
                  l1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
                  l2.SetPtEtaPhiM(daughter2->pt(), daughter2->eta(), daughter2->phi(), daughter2->mass());
              }
          }
      }

      // status 23 muon, intermediate particle from a decay
      if(daughter->pdgId() == 13 && daughter->status() == 23)
      {
          // we have an intermediate status 23 muon, save the intermediate values as preFSR
          mu1preFSR.pt = daughter->pt();
          mu1preFSR.eta = daughter->eta();
          mu1preFSR.phi = daughter->phi();
          mu1preFSR.charge = daughter->charge();

          // If it did not radiate then the post and pre are the same
          mu1postFSR = mu1preFSR;

          // if it did radiate, get the postFSR final state, status 1, version of this muon
          // and overwrite the postFSR quantities
          for(unsigned int i=0; i<daughter->numberOfDaughters(); i++)
          {
              const reco::Candidate* postFSRcand = daughter->daughter(i);
              if(postFSRcand->pdgId() == 13 && daughter->status() == 1)
              {
                  mu1postFSR.pt = postFSRcand->pt();
                  mu1postFSR.eta = postFSRcand->eta();
                  mu1postFSR.phi = postFSRcand->phi();
                  mu1postFSR.charge = postFSRcand->charge();
              }
          }
      }
      // status 23 antimuon, intermediate particle from a decay
      else if(daughter->pdgId() == -13 && daughter->status() == 23)
      {
          // we have an intermediate status 23 muon, save the intermediate values as preFSR
          mu2preFSR.pt = daughter->pt();
          mu2preFSR.eta = daughter->eta();
          mu2preFSR.phi = daughter->phi();
          mu2preFSR.charge = daughter->charge();

          mu2postFSR = mu2preFSR;

          // if it did radiate, get the postFSR final state, status 1, version of this muon
          // and overwrite the postFSR quantities
          for(unsigned int i=0; i<daughter->numberOfDaughters(); i++)
          {
              const reco::Candidate* postFSRcand = daughter->daughter(i);
              if(postFSRcand->pdgId() == -13 && daughter->status() == 1)
              {
                  mu2postFSR.pt = postFSRcand->pt();
                  mu2postFSR.eta = postFSRcand->eta();
                  mu2postFSR.phi = postFSRcand->phi();
                  mu2postFSR.charge = postFSRcand->charge();
              }
          }
      }
      // final state muon
      else if(daughter->pdgId() == 13 && daughter->status() == 1)
      {
          // no intermediate status 23 muon that radiated only final state status 1, so pre and post are the same
          mu1preFSR.pt = daughter->pt();
          mu1preFSR.eta = daughter->eta();
          mu1preFSR.phi = daughter->phi();
          mu1preFSR.charge = daughter->charge();

          // no radiation, post and pre are the same
          mu1postFSR = mu1preFSR;
      }
      // final state antimuon
      else if(daughter->pdgId() == -13 && daughter->status() == 1)
      {
          // no intermediate status 23 muon that radiated only final state status 1, so pre and post are the same
          mu2preFSR.pt = daughter->pt();
          mu2preFSR.eta = daughter->eta();
          mu2preFSR.phi = daughter->phi();
          mu2preFSR.charge = daughter->charge();

          // no radiation, post and pre are the same
          mu2postFSR = mu2preFSR;
      }
    
  }

  // no muons found, Z/gamma* decayed to other leptons
  if(mu1preFSR.pt < 0 && mu1postFSR.pt < 0 && mu2preFSR.pt < 0 && mu2postFSR.pt < 0 && l1.M()!=0 && l2.M()!=0)
  {
          // use the mass to let us know what the Z/gamma* decayed to
          mu1preFSR.pt = -l1.M();
          mu1preFSR.eta = -l1.M();
          mu1preFSR.phi = -l1.M();;
          mu1preFSR.charge = -l1.M();

          mu2preFSR.pt = -l2.M();
          mu2preFSR.eta = -l2.M();
          mu2preFSR.phi = -l2.M();;
          mu2preFSR.charge = -l2.M();

          mu1postFSR = mu1preFSR;
          mu2postFSR = mu2preFSR;
  }

  // fill the appropriate boson and daughters
  if(boson->status() == 21)
  {
  // in DY this is an incoming quark, anti-quark annihilation 
  // 21 could be any incoming particle in other samples though
  
      // if the virtual photon went to two leptons then reconstruct it 
      if(l1.M() != 0 && l2.M() != 0)
      {
          mother = l1 + l2;

          bosonInfo.mass = mother.M(); 
          bosonInfo.pt   = mother.Pt();   
          bosonInfo.eta  = -111;  
          bosonInfo.y    = mother.Rapidity();    
          bosonInfo.phi  = mother.Phi();  
          bosonInfo.charge = 0;

          // Not sure what to do if the virtual photon decayed to a bunch of leptons
          if(moreThanOneLeptPair)
          {
              bosonInfo.mass = -333; 
              bosonInfo.pt   = -333;
              bosonInfo.eta  = -333;
              bosonInfo.y    = -333;    
              bosonInfo.phi  = -333;
              bosonInfo.charge = -333;
          }
      }
      else
      {
          bosonInfo.mass = -999; 
          bosonInfo.pt   = -999;
          bosonInfo.eta  = -999;
          bosonInfo.y    = -999;    
          bosonInfo.phi  = -999;
          bosonInfo.charge = -999;
      }

      _genGpreFSR = bosonInfo;
      _genM1GpreFSR = mu1preFSR;
      _genM2GpreFSR = mu2preFSR;

      _genGpostFSR = bosonInfo;
      _genM1GpostFSR = mu1postFSR;
      _genM2GpostFSR = mu2postFSR;
  }
  // Z
  if(abs(boson->pdgId()) == 23)
  {
      _genZpreFSR = bosonInfo;
      _genM1ZpreFSR = mu1preFSR;
      _genM2ZpreFSR = mu2preFSR;

      _genZpostFSR = bosonInfo;
      _genM1ZpostFSR = mu1postFSR;
      _genM2ZpostFSR = mu2postFSR;
  }
  // W
  if(abs(boson->pdgId()) == 24)
  {
      _genWpreFSR = bosonInfo;
      if(bosonInfo.charge == mu1preFSR.charge) _genMWpreFSR = mu1preFSR;
      if(bosonInfo.charge == mu2preFSR.charge) _genMWpreFSR = mu2preFSR;

      _genWpostFSR = bosonInfo;
      if(bosonInfo.charge == mu1postFSR.charge) _genMWpreFSR = mu1postFSR;
      if(bosonInfo.charge == mu2postFSR.charge) _genMWpreFSR = mu2postFSR;
  }
  // H
  if(abs(boson->pdgId()) == 25)
  {
      _genHpreFSR = bosonInfo;
      _genM1HpreFSR = mu1preFSR;
      _genM2HpreFSR = mu2preFSR;

      _genHpostFSR = bosonInfo;
      _genM1HpostFSR = mu1postFSR;
      _genM2HpostFSR = mu2postFSR;
  }

}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillElectron(unsigned int i, const edm::Ptr<pat::Electron>& e, const edm::Handle<reco::VertexCollection>& vertices, const edm::Event& iEvent) 
{
    // Electron Cut Based ID
    edm::Handle< edm::ValueMap<bool> > eVetoHandle, eLooseHandle, eMediumHandle, eTightHandle;
    iEvent.getByToken(_electronCutBasedIdTightToken, eTightHandle);
    iEvent.getByToken(_electronCutBasedIdMediumToken, eMediumHandle);
    iEvent.getByToken(_electronCutBasedIdLooseToken, eLooseHandle);
    iEvent.getByToken(_electronCutBasedIdVetoToken, eVetoHandle);

    _electronInfo.charge[i] = e->charge();
    _electronInfo.pt[i]     = e->pt();  
    _electronInfo.eta[i]    = e->eta(); 
    _electronInfo.phi[i]    = e->phi();

    // type of electron
    _electronInfo.isTightElectron[i]  =  (*eTightHandle)[e];
    _electronInfo.isMediumElectron[i] =  (*eMediumHandle)[e];
    _electronInfo.isLooseElectron[i]  =  (*eLooseHandle)[e];
    _electronInfo.isVetoElectron[i]   =  (*eVetoHandle)[e];

    _electronInfo.passConversionVeto[i] = e->passConversionVeto();

    reco::Vertex bestVtx1;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx!=vertices->end(); ++vtx)
    {
      if (!vtx->isValid()) continue;
  
      _electronInfo.d0_PV[i]= e->gsfTrack()->dxy(vtx->position());
      _electronInfo.dz_PV[i]= e->gsfTrack()->dz(vtx->position());
      bestVtx1 = *vtx;
    
      //exit at the first available vertex
      break;
    }
   
    _electronInfo.missingInnerHits[i] = e->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    // PF Isolation
    _electronInfo.isPFElectron[i] = e->isPF();

    _electronInfo.sumChargedHadronPtR03[i]   = e->pfIsolationVariables().sumChargedHadronPt  ;
    _electronInfo.sumNeutralHadronEtR03[i]   = e->pfIsolationVariables().sumNeutralHadronEt  ;
    _electronInfo.sumPhotonEtR03[i]          = e->pfIsolationVariables().sumPhotonEt         ;
    _electronInfo.sumPUPtR03[i]              = e->pfIsolationVariables().sumPUPt             ;
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillTau(unsigned int i, const pat::Tau& tau, const edm::Handle<reco::VertexCollection>& vertices, const edm::Event& iEvent) 
{

    if(!tau.leadChargedHadrCand().isNull())
    {
        pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
 
        _tauInfo.dz_PV[i] = packedLeadTauCand->dz();
        _tauInfo.d0_PV[i] = packedLeadTauCand->dxy();
    }

    _tauInfo.charge[i]  = tau.charge();
    _tauInfo.pt[i]      = tau.pt();
    _tauInfo.eta[i]     = tau.eta();
    _tauInfo.phi[i]     = tau.phi();
    _tauInfo.isPFTau[i] = tau.isPFTau();

    for(unsigned int id=0; id<_tauInfo.idArraySize && id<_tauIDNames.size(); id++)
    {
        _tauInfo.tauID[i][id] = tau.tauID(_tauIDNames[id]);
    }
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::displaySelection() {

  std::cout << "\n\n*** UFDiMuonsAnalyzer Configuration ***\n";

  // variable to cuts over
  std::cout << " - _isGlobal:         " << _isGlobal << std::endl;
  std::cout << " - _isTracker:        " << _isTracker << std::endl;
  std::cout << " - _ptMin:            " << _ptMin << std::endl;
  std::cout << " - _etaMax:           " << _etaMax<< std::endl;

  // module config parameters
  std::cout << " - _checkTrigger: " << _checkTrigger << std::endl;
  std::cout << " - Triggers To Probe:\n";
  unsigned int triggerSize = _triggerNames.size();
  for (unsigned int i=0; i < triggerSize; i++) 
    std::cout << "    * triggerNames["<<i<<"]: " << _triggerNames[i] << std::endl;
  
  std::cout << std::endl << std::endl;

}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::sortGenJetFunc(reco::GenJet i, reco::GenJet j){ return (i.pt()>j.pt()); }
bool UFDiMuonsAnalyzer::sortMuonFunc(pat::Muon i, pat::Muon j){ return (i.pt()>j.pt()); }
bool UFDiMuonsAnalyzer::sortElectronFunc(pat::Electron i, pat::Electron j){ return (i.pt()>j.pt()); }
bool UFDiMuonsAnalyzer::sortTauFunc(pat::Tau i, pat::Tau j){ return (i.pt()>j.pt()); }

