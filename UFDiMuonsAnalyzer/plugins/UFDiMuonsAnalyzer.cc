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

  // Get the collections designated from the config file and load them into the tokens
  _muonCollToken = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonColl"));
  _beamSpotToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  _prunedGenParticleToken = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticleTag" ));
  _packedGenParticleToken = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticleTag" ));
  _primaryVertexToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"));
  _metToken     = consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("metTag"));
  _pfJetsToken  = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("pfJetsTag"));
  _genJetsToken = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsTag"));
  _genEvtInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  _PupInfoToken = consumes< std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));

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
  
  int theRun   = iEvent.id().run();
  int theLumi  = iEvent.luminosityBlock();
  long long int theEvent = iEvent.id().event();
  int theBx    = iEvent.bunchCrossing();
  int theOrbit = iEvent.orbitNumber();
  
  eventInfo.run   = theRun;
  eventInfo.lumi  = theLumi;
  eventInfo.event = theEvent;
  eventInfo.bx    = theBx;
  eventInfo.orbit = theOrbit;


  // -----------------------------------------
  // VERTICES AND PILEUP
  // -----------------------------------------
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(_primaryVertexToken, vertices);
 
  for (int i=0;i<20;i++) {
    vertexInfo.isValid[i]  = 0;
    vertexInfo.x[i]        = -999;     
    vertexInfo.y[i]        = -999;     
    vertexInfo.z[i]        = -999;     
    vertexInfo.xErr[i]     = -999;
    vertexInfo.yErr[i]     = -999;
    vertexInfo.zErr[i]     = -999;
    vertexInfo.chi2[i]     = -999;
    vertexInfo.ndf[i]      = -999;
    vertexInfo.normChi2[i] = -999;
  }      
  vertexInfo.nVertices   = 0;
  
  // primary vertex

  // init (vertices)
  if (vertices.isValid()) {
    
    //std::cout << "vertex->size():"<< vertices->size() << std::endl;
    int iVertex=0;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) {
        vertexInfo.isValid[iVertex] = 0;
        continue;
      }
    
      vertexInfo.isValid[iVertex] = 1;
      
      vertexInfo.x[iVertex]        = vtx->position().X();	
      vertexInfo.y[iVertex]        = vtx->position().Y();	
      vertexInfo.z[iVertex]        = vtx->position().Z();	
      vertexInfo.xErr[iVertex]     = vtx->xError();	
      vertexInfo.yErr[iVertex]     = vtx->yError();	
      vertexInfo.zErr[iVertex]     = vtx->zError();	
      vertexInfo.chi2[iVertex]     = vtx->chi2();	
      vertexInfo.ndf[iVertex]      = vtx->ndof();	
      vertexInfo.normChi2[iVertex] = vtx->normalizedChi2();
      
      iVertex++;
      vertexInfo.nVertices++;
    }
  }
  else std::cout << "VertexCollection is NOT valid -> vertex Info NOT filled!\n";
  
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
  // MONTE CARLO GEN INFO: MUONS, H, W, Z
  // -----------------------------------------
  
  if (_isMonteCarlo) {

    // initialize Z to default values
    initGenPart(_genZpreFSR); initTrack(_genM1ZpreFSR); initTrack(_genM2ZpreFSR);
    initGenPart(_genZpostFSR);initTrack(_genM1ZpostFSR);initTrack(_genM2ZpostFSR);

    // initialize H to default values
    initGenPart(_genHpreFSR); initTrack(_genM1HpreFSR); initTrack(_genM2HpreFSR);
    initGenPart(_genHpostFSR);initTrack(_genM1HpostFSR);initTrack(_genM2HpostFSR);

    // initialize W to default values
    initGenPart(_genWpreFSR); initTrack(_genMWpreFSR);
    initGenPart(_genWpostFSR); initTrack(_genMWpostFSR);

    edm::Handle<reco::GenParticleCollection> prunedGenParticles;
    iEvent.getByToken(_prunedGenParticleToken, prunedGenParticles);

    reco::GenParticleCollection hardProcessMuons;

    bool foundW = false;
    bool foundZ = false;
    bool foundH = false;
    //std::cout << "\n====================================\n"; 
    for (reco::GenParticleCollection::const_iterator gen = prunedGenParticles->begin(), genEnd = prunedGenParticles->end(); 
         gen != genEnd; ++gen) 
    {
      
        int id = gen->pdgId();
        int genstatus = gen->status();

        // In Pythia 6, take status code 3 to mean hard process particle
        // and status code 2 to mean later in the process

        // Status code 22 means hard process intermediate particle in Pythia 8
        if (abs(id) == 23 && (genstatus == 22 || genstatus == 3)) {
          foundZ = true;
          _genZpreFSR.mass = gen->mass(); 
          _genZpreFSR.pt   = gen->pt();   
          _genZpreFSR.eta  = gen->eta();  
          _genZpreFSR.y    = gen->rapidity();    
          _genZpreFSR.phi  = gen->phi();  
          continue;
        }
        if (abs(id) == 24 && (genstatus == 22 || genstatus == 3)) {
          foundW = true;
          _genWpreFSR.mass = gen->mass(); 
          _genWpreFSR.pt   = gen->pt();   
          _genWpreFSR.eta  = gen->eta();  
          _genWpreFSR.y    = gen->rapidity();    
          _genWpreFSR.phi  = gen->phi();  
          continue;
        }
        if (abs(id) == 25 && (genstatus == 22 || genstatus == 3)) {
          foundH = true;
          _genHpreFSR.mass = gen->mass(); 
          _genHpreFSR.pt   = gen->pt();   
          _genHpreFSR.eta  = gen->eta();  
          _genHpreFSR.y    = gen->rapidity();    
          _genHpreFSR.phi  = gen->phi();  
          continue;
        }
     
        // Status code 23 means hard process outgoing particle in Pythia 8
        if (abs(id) == 13 && (genstatus == 23 || genstatus == 3))  {
          hardProcessMuons.push_back(*gen);
          continue;
        } // muon 

        // Status code 62 means after ISR,FSR, and primoridial pt from UE
        if (abs(id) == 23 && (genstatus == 62 || genstatus == 2)) {
          _genZpostFSR.mass = gen->mass(); 
          _genZpostFSR.pt   = gen->pt();   
          _genZpostFSR.eta  = gen->eta();  
          _genZpostFSR.y    = gen->rapidity();    
          _genZpostFSR.phi  = gen->phi();  
          continue;
        }
        if (abs(id) == 25 && (genstatus == 62 || genstatus == 2)) {
          _genHpostFSR.mass = gen->mass(); 
          _genHpostFSR.pt   = gen->pt();   
          _genHpostFSR.eta  = gen->eta();  
          _genHpostFSR.y    = gen->rapidity();    
          _genHpostFSR.phi  = gen->phi();  
          continue;
        }

        if (abs(id) == 24 && (genstatus == 62 || genstatus == 2)) {
          foundW = true;
          _genWpostFSR.mass = gen->mass(); 
          _genWpostFSR.pt   = gen->pt();   
          _genWpostFSR.eta  = gen->eta();  
          _genWpostFSR.y    = gen->rapidity();    
          _genWpostFSR.phi  = gen->phi();  
          continue;
        }

    } // loop over gen level

    if (foundZ && hardProcessMuons.size()==2){
      reco::GenParticle& gen1 = hardProcessMuons[0];
      reco::GenParticle& gen2 = hardProcessMuons[1];
      _genM1ZpreFSR.pt = gen1.pt();
      _genM1ZpreFSR.eta = gen1.eta();
      _genM1ZpreFSR.phi = gen1.phi();
      _genM1ZpreFSR.charge = gen1.charge();
      _genM2ZpreFSR.pt = gen2.pt();
      _genM2ZpreFSR.eta = gen2.eta();
      _genM2ZpreFSR.phi = gen2.phi();
      _genM2ZpreFSR.charge = gen2.charge();
    }
    if (foundW && hardProcessMuons.size()==1){
      reco::GenParticle& gen1 = hardProcessMuons[0];
      _genMWpreFSR.pt = gen1.pt();
      _genMWpreFSR.eta = gen1.eta();
      _genMWpreFSR.phi = gen1.phi();
      _genMWpreFSR.charge = gen1.charge();
    }
    if (foundH && hardProcessMuons.size()==2){
      reco::GenParticle& gen1 = hardProcessMuons[0];
      reco::GenParticle& gen2 = hardProcessMuons[1];
      _genM1HpreFSR.pt = gen1.pt();
      _genM1HpreFSR.eta = gen1.eta();
      _genM1HpreFSR.phi = gen1.phi();
      _genM1HpreFSR.charge = gen1.charge();
      _genM2HpreFSR.pt = gen2.pt();
      _genM2HpreFSR.eta = gen2.eta();
      _genM2HpreFSR.phi = gen2.phi();
      _genM2HpreFSR.charge = gen2.charge();
    }

    edm::Handle<pat::PackedGenParticleCollection> packedGenParticles;
    iEvent.getByToken(_packedGenParticleToken, packedGenParticles);
    pat::PackedGenParticleCollection finalStateGenMuons;
  
    for (pat::PackedGenParticleCollection::const_iterator gen = packedGenParticles->begin(), 
            genEnd = packedGenParticles->end(); 
         gen != genEnd; ++gen) 
    {
        if (abs(gen->pdgId()) == 13)  {
          finalStateGenMuons.push_back(*gen);
        } // muon 
    }

    if (foundZ && finalStateGenMuons.size()==2){
      pat::PackedGenParticle& gen1 = finalStateGenMuons[0];
      pat::PackedGenParticle& gen2 = finalStateGenMuons[1];
      _genM1ZpostFSR.pt = gen1.pt();
      _genM1ZpostFSR.eta = gen1.eta();
      _genM1ZpostFSR.phi = gen1.phi();
      _genM1ZpostFSR.charge = gen1.charge();
      _genM2ZpostFSR.pt = gen2.pt();
      _genM2ZpostFSR.eta = gen2.eta();
      _genM2ZpostFSR.phi = gen2.phi();
      _genM2ZpostFSR.charge = gen2.charge();
    }
    if (foundW && finalStateGenMuons.size()==1){
      pat::PackedGenParticle& gen1 = finalStateGenMuons[0];
      _genMWpostFSR.pt = gen1.pt();
      _genMWpostFSR.eta = gen1.eta();
      _genMWpostFSR.phi = gen1.phi();
      _genMWpostFSR.charge = gen1.charge();
    }
    if (foundH && finalStateGenMuons.size()==2){
      pat::PackedGenParticle& gen1 = finalStateGenMuons[0];
      pat::PackedGenParticle& gen2 = finalStateGenMuons[1];
      _genM1HpostFSR.pt = gen1.pt();
      _genM1HpostFSR.eta = gen1.eta();
      _genM1HpostFSR.phi = gen1.phi();
      _genM1HpostFSR.charge = gen1.charge();
      _genM2HpostFSR.pt = gen2.pt();
      _genM2HpostFSR.eta = gen2.eta();
      _genM2HpostFSR.phi = gen2.phi();
      _genM2HpostFSR.charge = gen2.charge();
    }

  
  }// end _isMonteCarlo

  // -----------------------------------------
  // JETS, GENJETS
  // -----------------------------------------
  
  edm::Handle < std::vector<pat::MET> > mets;
  if(!_metToken.isUninitialized()) iEvent.getByToken(_metToken, mets);
  bzero(&_metInfo,sizeof(_MetInfo));

  if( mets.isValid() ){
    _metInfo.px = (*mets)[0].px();
    _metInfo.py = (*mets)[0].py();
    _metInfo.pt = (*mets)[0].pt();
    _metInfo.phi= (*mets)[0].phi();
    _metInfo.sumEt = (*mets)[0].sumEt();
  }

  edm::Handle < std::vector<pat::Jet> > jets;
  if(!_pfJetsToken.isUninitialized()) iEvent.getByToken(_pfJetsToken, jets);
  bzero(&_pfJetInfo,sizeof(_PFJetInfo));

  //// Get JEC Uncertainty Calculator
  //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl); 
  //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  //JetCorrectionUncertainty *jecUncCalculator = new JetCorrectionUncertainty(JetCorPar);

  if( jets.isValid() )
  {

    for(unsigned int i=0; i<jets->size(); i++)
    {
      _pfJetInfo.nJets++;
      if( i<10 )
      {
        const pat::Jet& jet = jets->at(i);
        _pfJetInfo.px[i] = jet.px();
        _pfJetInfo.py[i] = jet.py();
        _pfJetInfo.pz[i] = jet.pz();
        _pfJetInfo.pt[i] = jet.pt();
        _pfJetInfo.eta[i]= jet.eta();
        _pfJetInfo.phi[i]= jet.phi();
        _pfJetInfo.mass[i]  = jet.mass();
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
        _pfJetInfo.csv[i]  = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
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
  bzero(&_genJetInfo,sizeof(_GenJetInfo));

  if( genJets.isValid() ){
    reco::GenJetCollection sortedGenJets = (*genJets);
    sort(sortedGenJets.begin(), sortedGenJets.end(), sortGenJetFunc);
    for(unsigned int i=0; i<sortedGenJets.size(); i++){
      _genJetInfo.nJets++;
      if( i<10 ){
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

  // pre-selection: just check the muons are at least tracker muons.
  for (pat::MuonCollection::const_iterator muon = muons->begin(), 
         muonsEnd = muons->end(); muon !=muonsEnd; ++muon){
    
    // 
    if (!isPreselected(*muon, beamSpotHandle)) {
      if (_isVerbose) std::cout << "Muon NOT passing pre-selections\n"; 
      continue;
    }
    if (_isVerbose) std::cout << "Muon passing the defined pre-selections\n\n";     


    // all cuts but isolation
    if (!passKinCuts(*muon, beamSpotHandle)) {
      if (_isVerbose) std::cout << "Muon NOT passing kinematic selections\n"; 
      continue;
    }
    if (_isVerbose) std::cout << "Muon passing the defined kinematic selections\n\n";     
    
    // put this muons in the collection
    muonsSelected.push_back(*muon);

  }
  // ===========================================================================
 
  if (_isVerbose) std::cout << " Found " << muonsSelected.size() 
			    << " muons candidates\n\n";

  // at least _nMuons muons
  if ( muonsSelected.size() < _nMuons  ) return;

  // -----------------------------------------
  // INIT MUON VALUES
  // -----------------------------------------
  
  initMuon(_muon1);
  initMuon(_muon2);

  initTrack(_muon1vc);
  initTrack(_muon2vc);

  initTrack(_muon1pvc);
  initTrack(_muon2pvc);

  // -----------------------------------------
  // INIT DIMUON VALUES
  // -----------------------------------------
  
  _recoCandMass = -999;
  _recoCandPt   = -999;
  _recoCandEta  = -999;
  _recoCandY    = -999;
  _recoCandPhi  = -999;

  _recoCandMassPF = -999;
  _recoCandPtPF   = -999;
  _recoCandEtaPF  = -999;
  _recoCandYPF    = -999;
  _recoCandPhiPF  = -999;

  _recoCandMassVC = -999; 
  _recoCandMassResVC    = -999;
  _recoCandMassResCovVC = -999;
  _recoCandPtVC   = -999;
  _recoCandEtaVC  = -999;
  _recoCandYVC    = -999;
  _recoCandPhiVC  = -999;

  _recoCandMassPVC = -999; 
  _recoCandMassResPVC    = -999;
  _recoCandMassResCovPVC = -999;
  _recoCandPtPVC   = -999;
  _recoCandEtaPVC  = -999;
  _recoCandYPVC    = -999;
  _recoCandPhiPVC  = -999;

  _angleDiMuons = -999;
  _vertexIsValid = -999;
  _vertexNormChiSquare = -999;
  _vertexChiSquare = -999;
  _vertexNDF = -999;
  _vertexX = -999;
  _vertexY = -999;
  _vertexZ = -999;

 

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
    _muon1.isGlobal     = mu.isGlobalMuon(); 
    _muon1.isTracker    = mu.isTrackerMuon(); 
    _muon1.isStandAlone = mu.isStandAloneMuon(); 

    reco::Track track;
    if      (mu.isGlobalMuon())  track = *(mu.globalTrack());
    else if (mu.isTrackerMuon()) track = *(mu.innerTrack());
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon1.charge = mu.charge();
    _muon1.pt     = mu.pt();  
    _muon1.ptErr  = track.ptError(); 
    _muon1.eta    = mu.eta(); 
    _muon1.phi    = mu.phi();
   
    // redundant if the muon is tracker-only muon
    if (mu.isTrackerMuon()) {
      _muon1.trkPt   = mu.innerTrack()->pt();                    
      _muon1.trkPtErr= mu.innerTrack()->ptError();
      _muon1.trketa  = mu.innerTrack()->eta();                   
      _muon1.trkPhi  = mu.innerTrack()->phi();                   
    }

    _muon1.d0_BS= mu.innerTrack()->dxy( beamSpotHandle->position() );
    _muon1.dz_BS= mu.innerTrack()->dz ( beamSpotHandle->position() );

    reco::Vertex bestVtx1;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon1.d0_PV= track.dxy( vtx->position() );
      _muon1.dz_PV= track.dz ( vtx->position() );
      bestVtx1 = *vtx;
    
      //exit at the first available vertex
      break;
    }

    // type of muon
    _muon1.isTightMuon  =  muon::isTightMuon(mu, bestVtx1);
    _muon1.isMediumMuon =  muon::isMediumMuon(mu);
    _muon1.isLooseMuon  =  muon::isLooseMuon(mu);

    //isolation
    _muon1.trackIsoSumPt    = mu.isolationR03().sumPt;
    _muon1.trackIsoSumPtCorr= mu.isolationR03().sumPt; // no correction with only 1 muon
    _muon1.ecalIso          = mu.isolationR03().emEt;
    _muon1.hcalIso          = mu.isolationR03().hadEt;

    double isovar = mu.isolationR03().sumPt;
    isovar += mu.isolationR03().hadEt; //tracker + HCAL 
    isovar /= mu.pt(); // relative combine isolation
    _muon1.relCombIso=isovar;

  
    // PF Isolation
    _muon1.isPFMuon = mu.isPFMuon();
    
    if ( mu.isPFMuon() ) {
      
      reco::Candidate::LorentzVector pfmuon = mu.pfP4();
      
      _muon1.pfPt  = pfmuon.Pt();
      _muon1.pfEta = pfmuon.Eta();
      _muon1.pfPhi = pfmuon.Phi();

      _muon1.sumChargedHadronPtR03   = mu.pfIsolationR03().sumChargedHadronPt  ;
      _muon1.sumChargedParticlePtR03 = mu.pfIsolationR03().sumChargedParticlePt;
      _muon1.sumNeutralHadronEtR03   = mu.pfIsolationR03().sumNeutralHadronEt  ;
      _muon1.sumPhotonEtR03          = mu.pfIsolationR03().sumPhotonEt         ;
      _muon1.sumPUPtR03              = mu.pfIsolationR03().sumPUPt             ;
      
      _muon1.sumChargedHadronPtR04   = mu.pfIsolationR04().sumChargedHadronPt  ;
      _muon1.sumChargedParticlePtR04 = mu.pfIsolationR04().sumChargedParticlePt;
      _muon1.sumNeutralHadronEtR04   = mu.pfIsolationR04().sumNeutralHadronEt  ;
      _muon1.sumPhotonEtR04          = mu.pfIsolationR04().sumPhotonEt         ;
      _muon1.sumPUPtR04              = mu.pfIsolationR04().sumPUPt             ;
    }

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<_triggerNames.size();iTrigger++) 
      _muon1.isHltMatched[iTrigger] = isHltMatched(iEvent, iSetup, _triggerNames[iTrigger], *_triggerObjsHandle, mu);

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

  // loop over the candidates
  for (MuonPairs::const_iterator pair= dimuons.begin(); 
       pair != dimuons.end(); ++pair){
    
    pat::Muon mu1  = pair->first;
    pat::Muon mu2  = pair->second;

    // - relative combined isolation -
    double isovar1 = mu1.isolationR03().sumPt;
    isovar1 += mu1.isolationR03().hadEt; //tracker + HCAL 
    isovar1 /= mu1.pt(); // relative combine isolation

    double isovar2 = mu2.isolationR03().sumPt;
    isovar2 += mu2.isolationR03().hadEt; //tracker + HCAL
    isovar2 /= mu2.pt(); // relative combine isolation
    
    _muon1.relCombIso=isovar1;
    _muon2.relCombIso=isovar2;

    _muon1.trackIsoSumPt=mu1.isolationR03().sumPt ;
    _muon2.trackIsoSumPt=mu2.isolationR03().sumPt ;
  
    // correction for the high boosts (if needed)
    _muon1.trackIsoSumPtCorr=mu1.isolationR03().sumPt ;
    _muon2.trackIsoSumPtCorr=mu2.isolationR03().sumPt ;

    if( mu2.innerTrack().isNonnull() ){
      double dEta =      mu2.track()->eta() - mu1.track()->eta();
      double dPhi = fabs(mu2.track()->phi() - mu1.track()->phi());
      if( dPhi>3.1415927 ) dPhi = 2.*3.1415927 - dPhi;
      if( sqrt(dEta*dEta+dPhi*dPhi)<0.3 && _muon1.trackIsoSumPt>0.9*mu2.track()->pt() ) 
        _muon1.trackIsoSumPtCorr = _muon1.trackIsoSumPt - mu2.track()->pt();
    }

    if( mu1.innerTrack().isNonnull() ){
      double dEta =      mu1.track()->eta() - mu2.track()->eta();
      double dPhi = fabs(mu1.track()->phi() - mu2.track()->phi());
      if( dPhi>3.1415927 ) dPhi = 2.*3.1415927 - dPhi;
      if( sqrt(dEta*dEta+dPhi*dPhi)<0.3 && _muon2.trackIsoSumPt>0.9*mu1.track()->pt() ) 
        _muon2.trackIsoSumPtCorr = _muon2.trackIsoSumPt - mu1.track()->pt();
    }


    // PF Isolation
    _muon1.isPFMuon = mu1.isPFMuon();
    
    if ( mu1.isPFMuon() ) {
      
      reco::Candidate::LorentzVector pfmuon = mu1.pfP4();
      
      _muon1.pfPt  = pfmuon.Pt();
      _muon1.pfEta = pfmuon.Eta();
      _muon1.pfPhi = pfmuon.Phi();

      _muon1.sumChargedHadronPtR03   = mu1.pfIsolationR03().sumChargedHadronPt  ;
      _muon1.sumChargedParticlePtR03 = mu1.pfIsolationR03().sumChargedParticlePt;
      _muon1.sumNeutralHadronEtR03   = mu1.pfIsolationR03().sumNeutralHadronEt  ;
      _muon1.sumPhotonEtR03          = mu1.pfIsolationR03().sumPhotonEt         ;
      _muon1.sumPUPtR03              = mu1.pfIsolationR03().sumPUPt             ;
                                         
      _muon1.sumChargedHadronPtR04   = mu1.pfIsolationR04().sumChargedHadronPt  ;
      _muon1.sumChargedParticlePtR04 = mu1.pfIsolationR04().sumChargedParticlePt;
      _muon1.sumNeutralHadronEtR04   = mu1.pfIsolationR04().sumNeutralHadronEt  ;
      _muon1.sumPhotonEtR04          = mu1.pfIsolationR04().sumPhotonEt         ;
      _muon1.sumPUPtR04              = mu1.pfIsolationR04().sumPUPt             ;
    }


    _muon2.isPFMuon = mu2.isPFMuon();

    if ( mu2.isPFMuon() ) {
      
      reco::Candidate::LorentzVector pfmuon = mu2.pfP4();
      
      _muon2.pfPt  = pfmuon.Pt();
      _muon2.pfEta = pfmuon.Eta();
      _muon2.pfPhi = pfmuon.Phi();

      _muon2.sumChargedHadronPtR03   = mu2.pfIsolationR03().sumChargedHadronPt  ;
      _muon2.sumChargedParticlePtR03 = mu2.pfIsolationR03().sumChargedParticlePt;
      _muon2.sumNeutralHadronEtR03   = mu2.pfIsolationR03().sumNeutralHadronEt  ;
      _muon2.sumPhotonEtR03          = mu2.pfIsolationR03().sumPhotonEt         ;
      _muon2.sumPUPtR03              = mu2.pfIsolationR03().sumPUPt             ;
                                         
      _muon2.sumChargedHadronPtR04   = mu2.pfIsolationR04().sumChargedHadronPt  ;
      _muon2.sumChargedParticlePtR04 = mu2.pfIsolationR04().sumChargedParticlePt;
      _muon2.sumNeutralHadronEtR04   = mu2.pfIsolationR04().sumNeutralHadronEt  ;
      _muon2.sumPhotonEtR04          = mu2.pfIsolationR04().sumPhotonEt         ;
      _muon2.sumPUPtR04              = mu2.pfIsolationR04().sumPUPt             ;
    }

    bool passSelection=false;

    // pass kinematic tighter selections
    if ( passKinCuts(mu1, beamSpotHandle) && 
         passKinCuts(mu2, beamSpotHandle)  ) passSelection=true; 

    //// chek the trigger
    //if (passSelection) {
    //    
    //  // set the event as not passed  
    //  // and then check it passes the trigger
    //  passSelection=!true;

    //  if ( isHltMatched(iEvent,iSetup,_triggerNames, mu1, mu2) ) {
    //    if (_isVerbose) std::cout << "Both Muons TIGHT and At Least One Matches the Trigger\n";
    //    passSelection=true;
    //  }
    //}
    // ===== END Requirements ====
    
    // do we pass the selection? Yes -> store the event
    // else move to the next event
    if (!passSelection) return;

    // store all the info
    // muon 1
    _muon1.isGlobal     = mu1.isGlobalMuon(); 
    _muon1.isTracker    = mu1.isTrackerMuon(); 
    _muon1.isStandAlone = mu1.isStandAloneMuon(); 


    reco::Track track1;
    if      (mu1.isGlobalMuon())  track1 = *(mu1.globalTrack());
    else if (mu1.isTrackerMuon()) track1 = *(mu1.innerTrack());
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon1.charge = mu1.charge();
    _muon1.pt     = mu1.pt();  
    _muon1.ptErr  = track1.ptError(); 
    _muon1.eta    = mu1.eta(); 
    _muon1.phi    = mu1.phi();
   
    if (mu1.isTrackerMuon()) {
      _muon1.trkPt   = mu1.innerTrack()->pt();                    
      _muon1.trkPtErr= mu1.innerTrack()->ptError();
      _muon1.trketa  = mu1.innerTrack()->eta();                   
      _muon1.trkPhi  = mu1.innerTrack()->phi();                   
    }

    _muon1.d0_BS= track1.dxy( beamSpotHandle->position() );
    _muon1.dz_BS= track1.dz ( beamSpotHandle->position() );

    reco::Vertex bestVtx1;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon1.d0_PV= track1.dxy( vtx->position() );
      _muon1.dz_PV= track1.dz ( vtx->position() );
      bestVtx1 = *vtx;
    
      //exit at the first available vertex
      break;
    }

    // type of muon
    _muon1.isTightMuon  =  muon::isTightMuon(mu1, bestVtx1);
    _muon1.isMediumMuon =  muon::isMediumMuon(mu1);
    _muon1.isLooseMuon  =  muon::isLooseMuon(mu1);

    //tracker iso and rel comb iso already taken care of
    _muon1.ecalIso = mu1.isolationR03().emEt ;
    _muon1.hcalIso = mu1.isolationR03().hadEt ;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<_triggerNames.size();iTrigger++) 
      _muon1.isHltMatched[iTrigger] = isHltMatched(iEvent, iSetup, _triggerNames[iTrigger], *_triggerObjsHandle, mu1);
    
    // muon 2
    _muon2.isGlobal     = mu2.isGlobalMuon(); 
    _muon2.isTracker    = mu2.isTrackerMuon(); 

    reco::Track track2;
    if      (mu2.isGlobalMuon())  track2 = *(mu2.globalTrack());
    else if (mu2.isTrackerMuon()) track2 = *(mu2.innerTrack());
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon2.charge = mu2.charge(); 
    _muon2.pt     = mu2.pt(); 
    _muon2.ptErr  = track2.ptError(); 
    _muon2.eta    = mu2.eta(); 
    _muon2.phi    = mu2.phi();
   
    if (mu2.isTrackerMuon()) {
      _muon2.trkPt   = mu2.innerTrack()->pt();                    
      _muon2.trkPtErr= mu2.innerTrack()->ptError();
      _muon2.trketa  = mu2.innerTrack()->eta();                   
      _muon2.trkPhi  = mu2.innerTrack()->phi();                   
    }

    _muon2.d0_BS= track2.dxy( beamSpotHandle->position() );
    _muon2.dz_BS= track2.dz ( beamSpotHandle->position() );
   
    reco::Vertex bestVtx2;
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon2.d0_PV= track2.dxy( vtx->position() );
      _muon2.dz_PV= track2.dz ( vtx->position() );
      bestVtx2 = *vtx;
    
      //exit at the first available vertex
      break;
    }

    // type of muon
    _muon2.isTightMuon  =  muon::isTightMuon(mu2, bestVtx2);
    _muon2.isMediumMuon =  muon::isMediumMuon(mu2);
    _muon2.isLooseMuon  =  muon::isLooseMuon(mu2);

    //tracker iso and rel comb iso already taken care of
    _muon2.ecalIso = mu2.isolationR03().emEt ;
    _muon2.hcalIso = mu2.isolationR03().hadEt ;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<_triggerNames.size();iTrigger++) 
      _muon2.isHltMatched[iTrigger] = isHltMatched(iEvent, iSetup, _triggerNames[iTrigger], *_triggerObjsHandle, mu2);
    
    // combine the info
    // muons collection
    TLorentzVector mother=GetLorentzVector(&*pair); 

    _recoCandMass = mother.M();
    _recoCandPt   = mother.Pt();
    _recoCandEta  = mother.PseudoRapidity();
    _recoCandY    = mother.Rapidity();
    _recoCandPhi  = mother.Phi();

    _angleDiMuons = acos(-mu1.track()->momentum().Dot(mu2.track()->momentum()/
                                                      mu1.track()->p()/mu2.track()->p())); 
    
    
    if ( mu1.isPFMuon() && mu2.isPFMuon() ) {

      TLorentzVector muon1_pf, muon2_pf, mother_pf;
      double const MASS_MUON = 0.105658367; //GeV/c2

      reco::Candidate::LorentzVector pfmuon1 = mu1.pfP4();
      reco::Candidate::LorentzVector pfmuon2 = mu2.pfP4();
      
      muon1_pf.SetPtEtaPhiM(pfmuon1.Pt(), pfmuon1.Eta(), pfmuon1.Phi(), MASS_MUON);
      muon2_pf.SetPtEtaPhiM(pfmuon2.Pt(), pfmuon2.Eta(), pfmuon2.Phi(), MASS_MUON);

      mother_pf = muon1_pf+muon2_pf;
      
      _recoCandMassPF = mother_pf.M();
      _recoCandPtPF   = mother_pf.Pt();
                              
      _recoCandEtaPF  = mother_pf.PseudoRapidity();
      _recoCandYPF    = mother_pf.Rapidity();
      _recoCandPhiPF  = mother_pf.Phi();
    
    }

    // ===========================================================================
    // store everything in a ntuple
    _outTree->Fill();
  
    if (_isVerbose) {
      std::cout<<"\t"<<theRun    <<":"<<theLumi<<":"<<theEvent
               <<"\t"<<mother.M() <<":"<<mother.Pt()
               <<"\t"<<_muon1.eta<<":"<<_muon2.eta
               <<"\t"<<_muon1.pt <<":"<<_muon2.pt
               <<std::endl;
    }
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
  _outTree->Branch( "eventInfo",  &eventInfo, "run/I:lumi/I:event/L:bx/I:orbit/I");
  
  // rho
  _outTree->Branch("rho"  ,            &_rho,              "rho/F"             );
  _outTree->Branch("rho25",            &_rho25,            "rho25/F"           );
  _outTree->Branch("rho25asHtoZZto4l", &_rho25asHtoZZto4l, "rho25asHtoZZto4l/F");

  _outTree->Branch("vertexInfo", &vertexInfo, "nVertices/I:isValid[20]/I:"
		   "x[20]/F:y[20]/F:z[20]/F:xErr[20]/F:yErr[20]/F:zErr[20]/F:"
		   "chi2[20]/F:ndf[20]/I:normChi2[20]/F");

  _outTree->Branch("reco1", &_muon1, 
                   "isTracker/I:isStandAlone/I:isGlobal/I:"
                   "isTightMuon/I:isMediumMuon/I:isLooseMuon/I:"
                   "charge/I:pt/F:ptErr/F:eta/F:phi/F:"
                   "trkPt/F:trkPtErr/F:trkEta/F:trkPhi/F:"
                   "d0_BS/F:dz_BS/F:"
                   "d0_PV/F:dz_PV/F:"
                   "trackIsoSumPt/F:"
                   "trackIsoSumPtCorr/F:"      
                   "hcalIso/F:"
                   "ecalIso/F:"
                   "relCombIso/F:"
                   "isPFMuon/I:"
                   "pfPt/F:"
                   "pfEta/F:"
                   "pfPhi/F:"
                   "sumChargedHadronPtR03/F:"
                   "sumChargedParticlePtR03/F:"
                   "sumNeutralHadronEtR03/F:"
                   "sumPhotonEtR03/F:"
                   "sumPUPtR03/F:"
                   "sumChargedHadronPtR04/F:"
                   "sumChargedParticlePtR04/F:"
                   "sumNeutralHadronEtR04/F:"
                   "sumPhotonEtR04/F:"
                   "sumPUPtR04/F:"
                   "isHltMatched[3]/I:"
                   "hltPt[3]/F:"
                   "hltEta[3]/F:"
                   "hltPhi[3]/F");

  _outTree->Branch("reco2", &_muon2, 
                   "isTracker/I:isStandAlone/I:isGlobal/I:"
                   "isTightMuon/I:isMediumMuon/I:isLooseMuon/I:"
                   "charge/I:pt/F:ptErr/F:eta/F:phi/F:"
                   "trkPt/F:trkPtErr/F:trkEta/F:trkPhi/F:"
                   "d0_BS/F:dz_BS/F:"
                   "d0_PV/F:dz_PV/F:"
                   "trackIsoSumPt/F:"      
                   "trackIsoSumPtCorr/F:"      
                   "hcalIso/F:"
                   "ecalIso/F:"
                   "relCombIso/F:"      
                   "isPFMuon/I:"
                   "pfPt/F:"
                   "pfEta/F:"
                   "pfPhi/F:"
                   "sumChargedHadronPtR03/F:"
                   "sumChargedParticlePtR03/F:"
                   "sumNeutralHadronEtR03/F:"
                   "sumPhotonEtR03/F:"
                   "sumPUPtR03/F:"
                   "sumChargedHadronPtR04/F:"
                   "sumChargedParticlePtR04/F:"
                   "sumNeutralHadronEtR04/F:"
                   "sumPhotonEtR04/F:"
                   "sumPUPtR04/F:"
                   "isHltMatched[3]/I:"
                   "hltPt[3]/F:"
                   "hltEta[3]/F:"
                   "hltPhi[3]/F");

  _outTree->Branch("reco1vc", &_muon1vc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
  _outTree->Branch("reco2vc", &_muon2vc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

  _outTree->Branch("reco1pvc", &_muon1pvc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
  _outTree->Branch("reco2pvc", &_muon2pvc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

  _outTree->Branch("hltPaths",    &_triggerNames);

  // mass and pt for the fit
  _outTree->Branch("recoCandMass", &_recoCandMass, "recoCandMass/F");
  _outTree->Branch("recoCandPt"  , &_recoCandPt  , "recoCandPt/F");
  _outTree->Branch("recoCandEta" , &_recoCandEta , "recoCandEta/F");
  _outTree->Branch("recoCandY"   , &_recoCandY   , "recoCandY/F");
  _outTree->Branch("recoCandPhi" , &_recoCandPhi , "recoCandPhi/F");

  _outTree->Branch("recoCandMassPF", &_recoCandMassPF, "recoCandMassPF/F");
  _outTree->Branch("recoCandPtPF"  , &_recoCandPtPF  , "recoCandPtPF/F");
  _outTree->Branch("recoCandEtaPF" , &_recoCandEtaPF , "recoCandEtaPF/F");
  _outTree->Branch("recoCandYPF"   , &_recoCandYPF   , "recoCandYPF/F");
  _outTree->Branch("recoCandPhiPF" , &_recoCandPhiPF , "recoCandPhiPF/F");

  _outTree->Branch("recoCandMassVC", &_recoCandMassVC, "recoCandMassVC/F");
  _outTree->Branch("recoCandMassResVC",   &_recoCandMassResVC,   "recoCandMassResVC/F");
  _outTree->Branch("recoCandMassResCovVC",&_recoCandMassResCovVC,"recoCandMassResCovVC/F");
  _outTree->Branch("recoCandPtVC"  , &_recoCandPtVC  , "recoCandPtVC/F");
  _outTree->Branch("recoCandEtaVC" , &_recoCandEtaVC , "recoCandEtaVC/F");
  _outTree->Branch("recoCandYVC"   , &_recoCandYVC   , "recoCandYVC/F");
  _outTree->Branch("recoCandPhiVC" , &_recoCandPhiVC , "recoCandPhiVC/F");

  _outTree->Branch("recoCandMassPVC", &_recoCandMassPVC, "recoCandMassPVC/F");
  _outTree->Branch("recoCandMassResPVC",   &_recoCandMassResPVC,   "recoCandMassResPVC/F");
  _outTree->Branch("recoCandMassResCovPVC",&_recoCandMassResCovPVC,"recoCandMassResCovPVC/F");
  _outTree->Branch("recoCandPtPVC"  , &_recoCandPtPVC  , "recoCandPtPVC/F");
  _outTree->Branch("recoCandEtaPVC" , &_recoCandEtaPVC , "recoCandEtaPVC/F");
  _outTree->Branch("recoCandYPVC"   , &_recoCandYPVC   , "recoCandYPVC/F");
  _outTree->Branch("recoCandPhiPVC" , &_recoCandPhiPVC , "recoCandPhiPVC/F");

  _outTree->Branch("angleDiMuons",        &_angleDiMuons       ,"angleDiMuons/F");
  _outTree->Branch("vertexIsValid",       &_vertexIsValid      ,"vertexIsValid/I");          
  _outTree->Branch("vertexNormChiSquare", &_vertexNormChiSquare,"vertexNormChiSquare/F");  
  _outTree->Branch("vertexChiSquare",     &_vertexChiSquare    ,"vertexChiSquare/F");      
  _outTree->Branch("vertexNDF",           &_vertexNDF          ,"vertexNDF/F");            
  _outTree->Branch("vertexX",             &_vertexX            ,"vertexX/F");              
  _outTree->Branch("vertexY",             &_vertexY            ,"vertexY/F");              
  _outTree->Branch("vertexZ",             &_vertexZ            ,"vertexZ/F");              

  _outTree->Branch("met",    &_metInfo,   "px/F:py/F:pt/F:phi/F:sumEt/F");
  _outTree->Branch("pfJets", &_pfJetInfo, "nJets/I:px[10]/F:py[10]/F:pz[10]/F:pt[10]/F:eta[10]/F:phi[10]/F:mass[10]/F:charge[10]/I:partonFlavour[10]:chf[10]/F:nhf[10]/F:cef[10]/F:nef[10]/F:muf[10]/F:hfhf[10]/F:hfef[10]/F:cm[10]/I:chm[10]/I:nhm[10]/I:cem[10]/I:nem[10]/I:mum[10]/I:hfhm[10]/I:hfem[10]/I:jecFactor[10]/F:jecUnc[10]/F:csv[10]/F:genPx[10]/F:genPy[10]/F:genPz[10]/F:genPt[10]/F:genEta[10]/F:genPhi[10]/F:genMass[10]/F:genEMF[10]/F:genHadF[10]/F:genInvF[10]/F:genAux[10]/F");

  // MC information
  if (_isMonteCarlo) {

    _outTree->Branch("genZpreFSR",  &_genZpreFSR  ,"mass/F:pt/F:eta/F:y/F:phi/F");
    _outTree->Branch("genM1ZpreFSR",&_genM1ZpreFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("genM2ZpreFSR",&_genM2ZpreFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");

    _outTree->Branch("genZpostFSR",  &_genZpostFSR  ,"mass/F:pt/F:eta/F:y/F:phi/F");
    _outTree->Branch("genM1ZpostFSR",&_genM1ZpostFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("genM2ZpostFSR",&_genM2ZpostFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");

    // W block
    _outTree->Branch("genWpreFSR",  &_genWpreFSR  ,"mass/F:pt/F:eta/F:y/F:phi/F");
    _outTree->Branch("genMWpreFSR", &_genMWpreFSR ,"charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("genWpostFSR",  &_genWpostFSR  ,"mass/F:pt/F:eta/F:y/F:phi/F");
    _outTree->Branch("genMWpostFSR",&_genMWpostFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");

    // H block
    _outTree->Branch("genHpreFSR",  &_genHpreFSR  ,"mass/F:pt/F:eta/F:y/F:phi/F");
    _outTree->Branch("genM1HpreFSR",&_genM1HpreFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("genM2HpreFSR",&_genM2HpreFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");

    _outTree->Branch("genHpostFSR",  &_genHpostFSR  ,"mass/F:pt/F:eta/F:y/F:phi/F");
    _outTree->Branch("genM1HpostFSR",&_genM1HpostFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("genM2HpostFSR",&_genM2HpostFSR,"charge/I:pt/F:ptErr/F:eta/F:phi/F");

    _outTree->Branch("genJets", &_genJetInfo, "nJets/I:px[10]/F:py[10]/F:pz[10]/F:pt[10]/F:eta[10]/F:phi[10]/F:mass[10]/F:charge[10]/I");

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
  _outTreeMetadata->Branch("originalNumEvents"  ,            &_numEvents,            "originalNumEvents/I"             );
  _outTreeMetadata->Branch("sumEventWeights"  ,            &_sumEventWeights,      "sumEventWeights/I"             );
  _outTreeMetadata->Branch("isMonteCarlo"  ,            &_isMonteCarlo,              "isMonteCarlo/O"             );
  std::vector <std::string> * triggerNamesPointer = &_triggerNames;
  _outTreeMetadata->Branch("triggerNames"  ,"std::vector< std::string > >", &triggerNamesPointer);
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
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::isPreselected(const pat::Muon& muon,
                                      edm::Handle<reco::BeamSpot> beamSpotHandle)
{
// Most basic type cuts on the muons tested here
// Currently only checks if the muon is at least a tracker muon

  bool pass=false;

  if (_isVerbose) {
    std::cout<< "is Global?"     << muon.isGlobalMuon()     << std::endl;
    std::cout<< "is Tracker?"    << muon.isTrackerMuon()    << std::endl;
    std::cout<< "is StandAlone?" << muon.isStandAloneMuon() << std::endl;
  }

  // reconstruction cuts
  if (!muon.isGlobalMuon()  && _isGlobal ) return pass; // gbl muon
  if (!muon.isTrackerMuon() && _isTracker) return pass; // trk muon

  // do not accept muons which are standalone only
  // cannot get the cocktail fit which needs at least tracker...
  if(!muon.isGlobalMuon() && !muon.isTrackerMuon()) return pass;
  
  // if all the cuts are passed
  pass=true;
  return pass;
  
}

////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::passKinCuts(const pat::Muon& muon,
                                    edm::Handle<reco::BeamSpot> beamSpotHandle)
{
// Most basic kinematic cuts on the muons tested here
// Redundantly checks that the muon is at least a tracker muon
// Check that pt is greater than ptmin specified in the python config file 
// Check that the muon has eta less than eta max specified in python config file
  
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

void UFDiMuonsAnalyzer::initMuon(_MuonInfo& muon) 
{
// initialize values in the muon info data structure

  muon.isTracker    = -999;
  muon.isStandAlone = -999;
  muon.isGlobal     = -999;

  muon.isTightMuon    = -999;
  muon.isMediumMuon   = -999;
  muon.isLooseMuon    = -999;

  muon.charge = -999;
  muon.pt     = -999;
  muon.eta    = -999; 
  muon.phi    = -999;
  
  muon.d0_BS= -999;
  muon.dz_BS= -999;
  
  muon.d0_PV= -999;
  muon.dz_PV= -999;
  
  muon.trackIsoSumPt     = -999;
  muon.trackIsoSumPtCorr = -999;
  muon.hcalIso           = -999;
  muon.ecalIso           = -999;
  muon.relCombIso        = -999;

  muon.isPFMuon = -999;

  muon.pfPt  = -999;
  muon.pfEta = -999;
  muon.pfPhi = -999;

  muon.sumChargedHadronPtR03   = -999;
  muon.sumChargedParticlePtR03 = -999;
  muon.sumNeutralHadronEtR03   = -999;
  muon.sumPhotonEtR03          = -999;
  muon.sumPUPtR03              = -999;
  
  muon.sumChargedHadronPtR04   = -999;
  muon.sumChargedParticlePtR04 = -999;
  muon.sumNeutralHadronEtR04   = -999;
  muon.sumPhotonEtR04          = -999;
  muon.sumPUPtR04              = -999;

  for (unsigned int iTrigger=0;iTrigger<3;iTrigger++) {
    muon.isHltMatched[iTrigger] = -999;
    muon.hltPt[iTrigger] = -999;
    muon.hltEta[iTrigger] = -999;
    muon.hltPhi[iTrigger] = -999;
  }

}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::initTrack(_TrackInfo& track) 
{
// Initialize track info data structure
  track.charge = -999; 
  track.pt     = -999; 
  track.ptErr  = -999;
  track.eta    = -999; 
  track.phi    = -999;
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::initGenPart(_genPartInfo& part)
{
// Initialize gen info data structure
  part.mass = -999;
  part.pt   = -999;
  part.eta  = -999;
  part.y    = -999;
  part.phi  = -999;
}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

bool UFDiMuonsAnalyzer::checkMother(const reco::Candidate &part,
                                    int momPdgId){

  bool matchFound = false;

  // loop over all the mothers
  int nMothers = part.numberOfMothers();
  bool hasMother = nMothers ? true : false;
  const reco::Candidate *mom = NULL;
  if (hasMother) mom=part.mother();          

  while (hasMother && mom) {
    // exit if
    // 1. match is found
    //std::cout << "   --- momPdgId = " << mom->pdgId() << std::endl;
    if (abs(mom->pdgId()) == abs(momPdgId)) {
      matchFound = true;
      hasMother=false;
    }

    // 2. we are at the parton level
    if (abs(mom->pdgId()) < 10) hasMother=false;

    // 3. there is no other mom
    if (mom->numberOfMothers() == 0) mom=NULL;
    else                             mom=mom->mother();
 
  }

  mom = NULL;
  delete mom;

  return matchFound;

}

////////////////////////////////////////////////////////////////////////////
//-- ----------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

void UFDiMuonsAnalyzer::fillDiMuonGenPart(const reco::GenParticleCollection &genColl,
                                          _genPartInfo& part,
                                          _TrackInfo&  muon1,
                                          _TrackInfo&  muon2) {


  if (genColl.size() != 2) return;

   muon1.charge = genColl[0].charge(); 
   muon1.pt     = genColl[0].pt(); 
   muon1.eta    = genColl[0].eta(); 
   muon1.phi    = genColl[0].phi();	

   muon2.charge = genColl[1].charge(); 
   muon2.pt     = genColl[1].pt(); 
   muon2.eta    = genColl[1].eta(); 
   muon2.phi    = genColl[1].phi();	

   TLorentzVector vtrue1, vtrue2, vtrueMother;

   vtrue1.SetPtEtaPhiM(genColl[0].pt(), 
                       genColl[0].eta(), 
                       genColl[0].phi(), 
                       genColl[0].mass());

   vtrue2.SetPtEtaPhiM(genColl[1].pt(), 
                       genColl[1].eta(), 
                       genColl[1].phi(), 
                       genColl[1].mass());

   vtrueMother = vtrue1+vtrue2;	

   part.mass = vtrueMother.M();
   part.pt   = vtrueMother.Pt();
   part.eta  = vtrueMother.PseudoRapidity();
   part.y    = vtrueMother.Rapidity();
   part.phi  = vtrueMother.Phi();

   return;
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

