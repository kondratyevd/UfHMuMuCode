// -*- C++ -*-
//
// Package:    UFDiMuonsAnalyzer
// Class:      UFDiMuonsAnalyzer
// 
/**\class UFDiMuonsAnalyzer UserArea/UFDiMuonsAnalyzer/src/UFDiMuonsAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Gian Piero Di Giovanni,32 4-B08,+41227674961,
//         Created:  Thur Oct 21 10:44:13 CEST 2010
// $Id: UFDiMuonsAnalyzer.cc,v 1.15 2013/06/14 14:23:35 digiovan Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/regex.hpp>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBranch.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"


#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

//
// math classes
//
#include "DataFormats/Math/interface/deltaR.h"
//
// trigger
// 
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//
// vertexing
//
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//
// gen particles
//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// 2010.11.21 Adding the Muon Cocktail
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Math/interface/deltaR.h"

// pfJets and MET
//#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

// PU Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/View.h"

bool sortGenJetFunc(reco::GenJet i, reco::GenJet j){ return (i.pt()>j.pt()); }

// Add the data formats
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

// general 
double const PDG_MASS_Z  = 91.1876;//GeV/c2
double const PDG_WIDTH_Z = 2.4952; //GeV/c2
double const MASS_MUON = 0.105658367;    //GeV/c2

//
// class declaration
//
class UFDiMuonsAnalyzer : public edm::EDAnalyzer {

public:

  explicit UFDiMuonsAnalyzer(const edm::ParameterSet&);
  ~UFDiMuonsAnalyzer();
  
  int _numEvents;

  edm::ESHandle<GlobalTrackingGeometry> globalTrackingGeometry;
  //  MuonServiceProxy* theService;

  // muons pairs
  typedef std::pair<pat::Muon,pat::Muon> MuonPair;
  typedef std::vector<MuonPair>            MuonPairs;

  // sort the dimuons putting in front the candidates whose mass is 
  // the closest to the Z PDG value
  struct sortMuonsClass {
    bool operator() (MuonPair pair1, 
  		     MuonPair pair2) {

      TLorentzVector muon11, muon12, dimuon1;

      //reco::TrackRef const muon11Track = pair1.first . innerTrack();
      //reco::TrackRef const muon12Track = pair1.second. innerTrack();

      //muon11.SetPtEtaPhiM(muon11Track->pt(), muon11Track->eta(), muon11Track->phi(), MASS_MUON);
      //muon12.SetPtEtaPhiM(muon12Track->pt(), muon12Track->eta(), muon12Track->phi(), MASS_MUON);

      reco::Track const muon11Track = *(pair1.first . innerTrack());
      reco::Track const muon12Track = *(pair1.second. innerTrack());

      muon11.SetPtEtaPhiM(muon11Track.pt(), muon11Track.eta(), muon11Track.phi(), MASS_MUON);
      muon12.SetPtEtaPhiM(muon12Track.pt(), muon12Track.eta(), muon12Track.phi(), MASS_MUON);

      dimuon1 = muon11+muon12;


      TLorentzVector muon21, muon22, dimuon2;

      //reco::TrackRef const muon21Track = pair2.first . innerTrack();
      //reco::TrackRef const muon22Track = pair2.second. innerTrack();

      //muon21.SetPtEtaPhiM(muon21Track->pt(), muon21Track->eta(), muon21Track->phi(), MASS_MUON);
      //muon22.SetPtEtaPhiM(muon22Track->pt(), muon22Track->eta(), muon22Track->phi(), MASS_MUON);

      reco::Track const muon21Track = *(pair2.first . innerTrack());
      reco::Track const muon22Track = *(pair2.second. innerTrack());

      muon21.SetPtEtaPhiM(muon21Track.pt(), muon21Track.eta(), muon21Track.phi(), MASS_MUON);
      muon22.SetPtEtaPhiM(muon22Track.pt(), muon22Track.eta(), muon22Track.phi(), MASS_MUON);

      dimuon2 = muon21+muon22;

      return ( fabs(dimuon1.M()-PDG_MASS_Z) < fabs(dimuon2.M()-PDG_MASS_Z) );
    }
  } sortMuonObject;

  // tracks pairs, e.g. cocktail
  typedef std::pair<reco::Track,reco::Track> TrackPair;
  typedef std::vector<TrackPair> TrackPairs;

  // general event information	
  _EventInfo eventInfo;

  // rho
  float _rho;
  float _rho25;
  float _rho25asHtoZZto4l;

  // vertex information
  _VertexInfo vertexInfo;


  // combined muon-muon info
  float _recoCandMass;
  float _recoCandPt;
  float _recoCandEta; // pseudo rapidity
  float _recoCandY;   // rapidity
  float _recoCandPhi; // phi

  // combined muon-muon info
  float _recoCandMassPF;
  float _recoCandPtPF;
  float _recoCandEtaPF; 
  float _recoCandYPF;   
  float _recoCandPhiPF; 

  // Vertex Constrained muon-muon info
  float _recoCandMassVC;
  float _recoCandMassResVC;
  float _recoCandMassResCovVC;
  float _recoCandPtVC;
  float _recoCandEtaVC; 
  float _recoCandYVC;   
  float _recoCandPhiVC; 

  // Vertex Constrained from PV info
  float _recoCandMassPVC;
  float _recoCandMassResPVC;
  float _recoCandMassResCovPVC;
  float _recoCandPtPVC;
  float _recoCandEtaPVC; 
  float _recoCandYPVC;   
  float _recoCandPhiPVC; 


  float _angleDiMuons;
  int   _vertexIsValid;
  float _vertexNormChiSquare;
  float _vertexChiSquare;
  float _vertexNDF;
  float _vertexX;
  float _vertexY;
  float _vertexZ;

  _MuonInfo _muon1, _muon2;
  // dimuon candidates from refit with Vertex Constraint
  _TrackInfo _muon1vc,_muon2vc; //dimuon
  _TrackInfo _muon1pvc,_muon2pvc; //pv



  // generator level info Z pre-FSR
  _genPartInfo _genZpreFSR;
  _TrackInfo   _genM1ZpreFSR,_genM2ZpreFSR;

  // generator level info Z post-FSR
  _genPartInfo _genZpostFSR;
  _TrackInfo   _genM1ZpostFSR,_genM2ZpostFSR;


  // generator level info W pre-FSR
  _genPartInfo _genWpreFSR;
  _TrackInfo   _genMWpreFSR;

  // generator level info W post-FSR
  _genPartInfo _genWpostFSR;
  _TrackInfo   _genMWpostFSR;


  // generator level info H pre-FSR
  _genPartInfo _genHpreFSR;
  _TrackInfo   _genM1HpreFSR,_genM2HpreFSR;

  // generator level info H post-FSR
  _genPartInfo _genHpostFSR;
  _TrackInfo   _genM1HpostFSR,_genM2HpostFSR;

  // Jets and MET
  _MetInfo     _metInfo;
  _PFJetInfo   _pfJetInfo;
  _GenJetInfo  _genJetInfo;

  int _nPU;
  int _genWeight;
  int _sumEventWeights;

  // where to save all the info  
  TTree* _outTree;
  TTree* _outTreeMetadata;


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::Service<TFileService> fs;

  void initMuon(_MuonInfo& muon);
  void initTrack(_TrackInfo& track);
  void initGenPart(_genPartInfo& part);

  bool checkMother(const reco::Candidate &part, int momPdgId);
  void fillDiMuonGenPart(const reco::GenParticleCollection &genColl,
                         _genPartInfo& part,
                         _TrackInfo&  muon1,
                         _TrackInfo&  muon2); 

  // method to select the muons
  bool isHltPassed (const edm::Event&, const edm::EventSetup&, const std::vector<std::string> triggerNames);
  bool isHltMatched(const edm::Event&, const edm::EventSetup&, 
                    const std::string& triggerName, 
                    const pat::TriggerObjectStandAloneCollection& triggerObjects,
                    const pat::Muon&);

  bool isPreselected(const pat::Muon& muon,
                     edm::Handle<reco::BeamSpot> beamSpotHandle);
  
  bool passKinCuts(const pat::Muon& muon,
                   edm::Handle<reco::BeamSpot> beamSpotHandle);

  void displaySelection();


  // muons
  edm::InputTag _muonColl;
  edm::InputTag _beamSpotTag;		
  edm::InputTag _prunedGenParticleTag;		
  edm::InputTag _packedGenParticleTag;		
  edm::InputTag _primaryVertexTag;		

  // variable to cuts over
  double _isGlobal;
  double _isTracker;  

  double _ptMin;
  double _etaMax;
  double _normChiSquareMax;
  double _d0Max;

  int _numPixelLayersMin;   
  int _numTrackerLayersMin; 
  int _numStripLayersMin;   

  double _validFracTrackerMin; 

  int _numValidMuonHitsMin;
  int _numValidPixelHitsMin;
  int _numValidStripHitsMin;
  int _numValidTrackerHitsMin;
  int _numSegmentMatchesMin;
  int _numOfMatchedStationsMin;

  // track iso
  double _trackIsoMaxSumPt;
  // relative combines iso
  double _relCombIsoMax;

  unsigned int _nMuons;

  // module config parameters
  std::string   processName_;
  std::vector < std::string > triggerNames_;

  std::vector < int > l1Prescale_;
  std::vector < int > hltPrescale_;

  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerObjsTag_;

  edm::InputTag metTag;
  edm::InputTag pfJetsTag;
  edm::InputTag genJetsTag;

  edm::InputTag puJetMvaFullDiscTag;
  edm::InputTag puJetMvaFullIdTag;
  edm::InputTag puJetMvaSimpleDiscTag;
  edm::InputTag puJetMvaSimpleIdTag;
  edm::InputTag puJetMvaCutDiscTag;
  edm::InputTag puJetMvaCutIdTag;

  // additional class data memebers
  bool _checkTrigger; // activate or not the trigger checking   
  edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
  edm::Handle<pat::TriggerObjectStandAloneCollection>   triggerObjsHandle_;

  MuonPairs  const GetMuonPairs (pat::MuonCollection  const* muons ) const;
  TrackPairs const GetTrackPairs(reco::TrackCollection const* tracks) const;

  TLorentzVector const GetLorentzVector(UFDiMuonsAnalyzer::MuonPair  const* pair) ;//const; 
  TLorentzVector const GetLorentzVector(UFDiMuonsAnalyzer::TrackPair const* pair) ;//const; 

  // useful switches
  bool _isVerbose;    // debug mode
  bool _isMonteCarlo; // save MC truth
  bool _isTriggerEmulated;
  double _emulatedPt;

  std::vector<float> getPUJetIDDisc(edm::Handle<edm::View<pat::Jet> >  jets, const edm::Event& event, edm::InputTag tag);
  std::vector<int> getPUJetID(edm::Handle<edm::View<pat::Jet> >  jets, const edm::Event& event, edm::InputTag tag);
};


//
// constructors and destructor
//
UFDiMuonsAnalyzer::UFDiMuonsAnalyzer(const edm::ParameterSet& iConfig):
  _numEvents(0)
{
  _sumEventWeights = 0;

  _outTree = fs->make<TTree>("tree", "myTree");
  _outTreeMetadata = fs->make<TTree>("metadata", "Metadata Tree");

  //now do what ever initialization is needed
  _muonColl	= iConfig.getParameter<edm::InputTag>("muonColl");

  _beamSpotTag	= iConfig.getParameter<edm::InputTag>("beamSpotTag");
  _prunedGenParticleTag	= iConfig.getParameter<edm::InputTag>("prunedGenParticleTag" );
  _packedGenParticleTag	= iConfig.getParameter<edm::InputTag>("packedGenParticleTag" );
  _primaryVertexTag	= iConfig.getParameter<edm::InputTag>("primaryVertexTag");

  _isVerbose	= iConfig.getUntrackedParameter<bool>("isVerbose",   false);
  _isMonteCarlo	= iConfig.getParameter<bool>("isMonteCarlo");
  _isTriggerEmulated   = iConfig.getUntrackedParameter<bool>("isTriggerEmulated",   false);
  _emulatedPt   = iConfig.getUntrackedParameter<double>("emulatedPt",  -999);

  // selection cuts
  _isGlobal               = iConfig.getParameter<int>("isGlobal");
  _isTracker              = iConfig.getParameter<int>("isTracker");

  if (!_isGlobal && !_isTracker) {
    std::cout << "\n\nWARNING: you are not requiring the muon to be not TRACKER nor GLOBAL\n"
	      << "Be aware of the fact that StandAlone muon only are rejected anyway in the code\n";
  }
 
  _ptMin                  = iConfig.getParameter<double>("ptMin");
  _etaMax                 = iConfig.getParameter<double>("etaMax");

  _normChiSquareMax	  = iConfig.getParameter<double>("normChiSquareMax");

  _numPixelLayersMin      = iConfig.getParameter<int>("minPixelLayers");  
  _numTrackerLayersMin    = iConfig.getParameter<int>("minTrackerLayers");
  _numStripLayersMin      = iConfig.getParameter<int>("minStripLayers");  
                                                 
  _validFracTrackerMin    = iConfig.getParameter<int>("validFracTrackerMin");

  _numValidMuonHitsMin	  = iConfig.getParameter<int>("minMuonHits");
  _numValidPixelHitsMin	  = iConfig.getParameter<int>("minPixelHits");
  _numValidStripHitsMin	  = iConfig.getParameter<int>("minStripHits");
  _numValidTrackerHitsMin = iConfig.getParameter<int>("minTrackerHits");
  _numSegmentMatchesMin	  = iConfig.getParameter<int>("minSegmentMatches");
  _numOfMatchedStationsMin= iConfig.getParameter<int>("minNumOfMatchedStations");
  
  _d0Max                  = iConfig.getParameter<double>("d0Max");
  _trackIsoMaxSumPt       = iConfig.getParameter<double>("trackIsoMaxSumPt");
  _relCombIsoMax          = iConfig.getParameter<double>("relCombIsoMax");

  _nMuons                 = iConfig.getParameter<int>("nMuons");


  metTag     = iConfig.getParameter<edm::InputTag>("metTag");
  pfJetsTag  = iConfig.getParameter<edm::InputTag>("pfJetsTag");
  genJetsTag = iConfig.getParameter<edm::InputTag>("genJetsTag");

  //HLT trigger initialization
  _checkTrigger	     = iConfig.getParameter<bool>("checkTrigger");

  processName_       = iConfig.getParameter<std::string>("processName");
  triggerNames_  = iConfig.getParameter<std::vector <std::string> >("triggerNames");

  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  triggerObjsTag_ = iConfig.getParameter<edm::InputTag>("triggerObjs");
}


UFDiMuonsAnalyzer::~UFDiMuonsAnalyzer() {}


// ------------ method called to for each event  ------------
void UFDiMuonsAnalyzer::analyze(const edm::Event& iEvent, 
                                const edm::EventSetup& iSetup)
{
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
    iEvent.getByLabel( edm::InputTag("generator"), genEvtInfo );
    _genWeight = (genEvtInfo->weight() > 0)? 1 : -1;
    _sumEventWeights += _genWeight;
  }

  if (_isVerbose) 
    std::cout << "\n\n A N A L Y Z I N G   E V E N T = " 
	      << _numEvents << std::endl << std::endl;
 
  // -----------------------------------------
  // H L T   H A N D L E S
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }
  iEvent.getByLabel(triggerObjsTag_,triggerObjsHandle_);
  if (!triggerObjsHandle_.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerObjects product from Event!" << std::endl;
    return;
  }

  // if all the HLT paths are not fired, will discard the event immediately
  if (_checkTrigger) 
  {
    if ( !isHltPassed(iEvent,iSetup,triggerNames_) )
    {
      if (_isVerbose) std::cout << "None of the HLT paths fired -> discard the event\n";
      return;
    }
     
  }

  //// once the triggernames are cleared add the prescales
  //for (unsigned int iTrigger=0; iTrigger<triggerNames_.size(); iTrigger++) {
  //  std::pair<int,int> prescaleValues = hltConfig_.prescaleValues(iEvent, iSetup,triggerNames_[iTrigger]);
  //  l1Prescale_[iTrigger]  = prescaleValues.first;
  //  hltPrescale_[iTrigger] = prescaleValues.second;

  //}

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


  //vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(_primaryVertexTag, vertices);
 
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
    iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);

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
  iEvent.getByLabel(_beamSpotTag, beamSpotHandle);

  //
  // get the sim mass if it exists
  //
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
    iEvent.getByLabel(_prunedGenParticleTag, prunedGenParticles);

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
    iEvent.getByLabel(_packedGenParticleTag, packedGenParticles);
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

  // ===========================================================================
  // Jet Info
  edm::Handle < std::vector<pat::MET> > mets;
  if( metTag.label() != "null" ) iEvent.getByLabel(metTag, mets);
  bzero(&_metInfo,sizeof(_MetInfo));

  if( mets.isValid() ){
    _metInfo.px = (*mets)[0].px();
    _metInfo.py = (*mets)[0].py();
    _metInfo.pt = (*mets)[0].pt();
    _metInfo.phi= (*mets)[0].phi();
    _metInfo.sumEt = (*mets)[0].sumEt();
  }

  edm::Handle < std::vector<pat::Jet> > jets;
  if( pfJetsTag.label() != "null" ) iEvent.getByLabel(pfJetsTag, jets);
  bzero(&_pfJetInfo,sizeof(_PFJetInfo));

  //// Get JEC Uncertainty Calculator
  //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PF",JetCorParColl); 
  //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  //JetCorrectionUncertainty *jecUncCalculator = new JetCorrectionUncertainty(JetCorPar);

  if( jets.isValid() ){

    for(unsigned int i=0; i<jets->size(); i++){
      _pfJetInfo.nJets++;
      if( i<10 ){
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
  if( genJetsTag.label() != "null" ) iEvent.getByLabel(genJetsTag, genJets);
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

  // ===========================================================================
  // M U O N S
  //
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(_muonColl, muons);

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

  // initialize the muons to default values
  initMuon(_muon1);
  initMuon(_muon2);

  initTrack(_muon1vc);
  initTrack(_muon2vc);

  initTrack(_muon1pvc);
  initTrack(_muon2pvc);

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

 


  if (muonsSelected.size() == 0) {
    if (_isVerbose) std::cout << "0 reco'd muons...\n";
    _outTree->Fill();
    return;
  }

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

    _muon1.normChiSquare=track.normalizedChi2();
    _muon1.d0_BS= mu.innerTrack()->dxy( beamSpotHandle->position() );
    _muon1.dz_BS= mu.innerTrack()->dz ( beamSpotHandle->position() );

    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon1.d0_PV= track.dxy( vtx->position() );
      _muon1.dz_PV= track.dz ( vtx->position() );
    
      //exit at the first available vertex
      break;
    }

    _muon1.numPixelLayers   = track.hitPattern().pixelLayersWithMeasurement();   
    _muon1.numTrackerLayers = track.hitPattern().trackerLayersWithMeasurement(); 
    _muon1.numStripLayers   = track.hitPattern().stripLayersWithMeasurement();   
    
    _muon1.validFracTracker = track.validFraction(); 

    _muon1.numValidMuonHits    = track.hitPattern().numberOfValidMuonHits();
    _muon1.numValidPixelHits   = track.hitPattern().numberOfValidPixelHits();
    _muon1.numValidTrackerHits = track.hitPattern().numberOfValidTrackerHits();
    _muon1.numValidStripHits   = track.hitPattern().numberOfValidStripHits();
    _muon1.numSegmentMatches   = mu.numberOfMatches();
    _muon1.numOfMatchedStations= mu.numberOfMatchedStations();

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

    _muon1.segmentCompatibility = muon::segmentCompatibility(mu);
    _muon1.combinedQualityChi2LocalPosition = mu.combinedQuality().chi2LocalPosition;
    _muon1.combinedQualityTrkKink = mu.combinedQuality().trkKink;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++) 
      _muon1.isHltMatched[iTrigger] = isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], *triggerObjsHandle_, mu);

    _outTree->Fill();
    
    // you can exit and pass to the new event
    return;
  }
  

  // ... if we got two or more muons we construct dimuon candidates
  
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
     
      
    // both muons have to be isolated
    if (_muon1.trackIsoSumPt > _trackIsoMaxSumPt) return;
    if (_muon2.trackIsoSumPt > _trackIsoMaxSumPt) return;

    if (_muon1.relCombIso > _relCombIsoMax) return;
    if (_muon2.relCombIso > _relCombIsoMax) return;

    // pass kinematic tighter selections
    if ( passKinCuts(mu1, beamSpotHandle) && 
         passKinCuts(mu2, beamSpotHandle)  ) passSelection=true; 

    //// chek the trigger
    //if (passSelection) {
    //    
    //  // set the event as not passed  
    //  // and then check it passes the trigger
    //  passSelection=!true;

    //  if ( isHltMatched(iEvent,iSetup,triggerNames_, mu1, mu2) ) {
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

    _muon1.normChiSquare=track1.normalizedChi2();
    _muon1.d0_BS= track1.dxy( beamSpotHandle->position() );
    _muon1.dz_BS= track1.dz ( beamSpotHandle->position() );

    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon1.d0_PV= track1.dxy( vtx->position() );
      _muon1.dz_PV= track1.dz ( vtx->position() );
    
      //exit at the first available vertex
      break;
    }

    _muon1.numPixelLayers   = track1.hitPattern().pixelLayersWithMeasurement();   
    _muon1.numTrackerLayers = track1.hitPattern().trackerLayersWithMeasurement(); 
    _muon1.numStripLayers   = track1.hitPattern().stripLayersWithMeasurement();   
    
    _muon1.validFracTracker = track1.validFraction(); 

    _muon1.numValidMuonHits    = track1.hitPattern().numberOfValidMuonHits();
    _muon1.numValidPixelHits   = track1.hitPattern().numberOfValidPixelHits();
    _muon1.numValidTrackerHits = track1.hitPattern().numberOfValidTrackerHits();
    _muon1.numValidStripHits   = track1.hitPattern().numberOfValidStripHits();
    _muon1.numSegmentMatches   = mu1.numberOfMatches();
    _muon1.numOfMatchedStations= mu1.numberOfMatchedStations();

    //tracker iso and rel comb iso already taken care of
    _muon1.ecalIso = mu1.isolationR03().emEt ;
    _muon1.hcalIso = mu1.isolationR03().hadEt ;

    _muon1.segmentCompatibility = muon::segmentCompatibility(mu1);
    _muon1.combinedQualityChi2LocalPosition = mu1.combinedQuality().chi2LocalPosition;
    _muon1.combinedQualityTrkKink = mu1.combinedQuality().trkKink;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++) 
      _muon1.isHltMatched[iTrigger] = isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], *triggerObjsHandle_, mu1);
    
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

    _muon2.numPixelLayers   = track2.hitPattern().pixelLayersWithMeasurement();   
    _muon2.numTrackerLayers = track2.hitPattern().trackerLayersWithMeasurement(); 
    _muon2.numStripLayers   = track2.hitPattern().stripLayersWithMeasurement();   
    
    _muon2.validFracTracker = track2.validFraction(); 

    _muon2.normChiSquare=track2.normalizedChi2();
    _muon2.d0_BS= track2.dxy( beamSpotHandle->position() );
    _muon2.dz_BS= track2.dz ( beamSpotHandle->position() );
   
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon2.d0_PV= track2.dxy( vtx->position() );
      _muon2.dz_PV= track2.dz ( vtx->position() );
    
      //exit at the first available vertex
      break;
    }

    _muon2.numValidMuonHits    = track2.hitPattern().numberOfValidMuonHits();
    _muon2.numValidPixelHits   = track2.hitPattern().numberOfValidPixelHits();
    _muon2.numValidTrackerHits = track2.hitPattern().numberOfValidTrackerHits();
    _muon2.numValidStripHits   = track2.hitPattern().numberOfValidStripHits();
    _muon2.numSegmentMatches   = mu2.numberOfMatches();
    _muon2.numOfMatchedStations= mu2.numberOfMatchedStations();

    //tracker iso and rel comb iso already taken care of
    _muon2.ecalIso = mu2.isolationR03().emEt ;
    _muon2.hcalIso = mu2.isolationR03().hadEt ;

    _muon2.segmentCompatibility = muon::segmentCompatibility(mu2);
    _muon2.combinedQualityChi2LocalPosition = mu2.combinedQuality().chi2LocalPosition;
    _muon2.combinedQualityTrkKink = mu2.combinedQuality().trkKink;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++) 
      _muon2.isHltMatched[iTrigger] = isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], *triggerObjsHandle_, mu2);
    
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
	
// ------------ method called once each job just before starting event loop  ------------
void UFDiMuonsAnalyzer::beginJob()
{
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
                   "charge/I:pt/F:ptErr/F:eta/F:phi/F:"
                   "trkPt/F:trkPtErr/F:trkEta/F:trkPhi/F:"
                   "normChiSquare/F:"
                   "d0_BS/F:dz_BS/F:"
                   "d0_PV/F:dz_PV/F:"
                   "numPixelLayers/I:"
                   "numTrackerLayers/I:"
                   "numStripLayers/I:"
                   "validFracTracker/F:"
                   "numValidMuonHits/I:"
                   "numValidPixelHits/I:"    
                   "numValidTrackerHits/I:"  
                   "numValidStripHits/I:"    
                   "numSegmentMatches/I:"    
                   "numOfMatchedStations/I:"
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
                   "hltPhi[3]/F:"
                   "segmentCompatibility/F:"
                   "combinedQualityChi2LocalPosition/F:"
                   "combinedQualityTrkKink/F");

  _outTree->Branch("reco2", &_muon2, 
                   "isTracker/I:isStandAlone/I:isGlobal/I:"
                   "charge/I:pt/F:ptErr/F:eta/F:phi/F:"
                   "trkPt/F:trkPtErr/F:trkEta/F:trkPhi/F:"
                   "normChiSquare/F:"
                   "d0_BS/F:dz_BS/F:"
                   "d0_PV/F:dz_PV/F:"
                   "numPixelLayers/I:"
                   "numTrackerLayers/I:"
                   "numStripLayers/I:"
                   "validFracTracker/F:"
                   "numValidMuonHits/I:"
                   "numValidPixelHits/I:"    
                   "numValidTrackerHits/I:"  
                   "numValidStripHits/I:"    
                   "numSegmentMatches/I:"   
                   "numOfMatchedStations/I:" 
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
                   "hltPhi[3]/F:"
                   "segmentCompatibility/F:"
                   "combinedQualityChi2LocalPosition/F:"
                   "combinedQualityTrkKink/F");

  _outTree->Branch("reco1vc", &_muon1vc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
  _outTree->Branch("reco2vc", &_muon2vc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

  _outTree->Branch("reco1pvc", &_muon1pvc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
  _outTree->Branch("reco2pvc", &_muon2pvc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

  _outTree->Branch("hltPaths",    &triggerNames_);
  _outTree->Branch("l1Prescale",  &l1Prescale_  );
  _outTree->Branch("hltPrescale", &hltPrescale_ );
 

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


// ------------ method called once each job just after ending the event loop  ------------
void UFDiMuonsAnalyzer::endJob() {

  std::cout << "Total Number of Events Read: "<< _numEvents << std::endl <<std::endl;
  std::cout << "Number of events weighted: "  << _sumEventWeights << std::endl <<std::endl;

  std::cout<<"number of candidate dimuons candidates: "
           <<_outTree->GetEntries()<<std::endl;

  // create the metadata tree branches
  _outTreeMetadata->Branch("originalNumEvents"  ,            &_numEvents,            "originalNumEvents/I"             );
  _outTreeMetadata->Branch("sumEventWeights"  ,            &_sumEventWeights,      "sumEventWeights/I"             );
  _outTreeMetadata->Branch("isMonteCarlo"  ,            &_isMonteCarlo,              "isMonteCarlo/O"             );
  std::vector <std::string> * triggerNamesPointer = &triggerNames_;
  _outTreeMetadata->Branch("triggerNames"  ,"std::vector< std::string > >", &triggerNamesPointer);
  _outTreeMetadata->Fill();

}


TLorentzVector const UFDiMuonsAnalyzer::GetLorentzVector(UFDiMuonsAnalyzer::MuonPair const* pair) {
  
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

// this method will simply check is the selected HLT path (via triggerName)
// is run and accepted and no error are found
//
// bool true  if (run && accept && !error)
//      false if any other combination
bool UFDiMuonsAnalyzer::isHltPassed(const edm::Event& iEvent, 
                                    const edm::EventSetup& iSetup, 
                                    const std::vector<std::string> desiredTriggerNames) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;


  const boost::regex re("_v[0-9]+");

  const TriggerNames &triggerNames = iEvent.triggerNames(*triggerResultsHandle_);

  const unsigned nTriggers = triggerResultsHandle_->size();
  for (unsigned iTrigger = 0; iTrigger < nTriggers; ++iTrigger)
  {
    const string triggerName = triggerNames.triggerName(iTrigger);
    string triggerNameStripped = boost::regex_replace(triggerName,re,"",boost::match_default | boost::format_sed);
    for(std::vector<std::string>::const_iterator desiredTriggerName=desiredTriggerNames.begin();
            desiredTriggerName!=desiredTriggerNames.end();desiredTriggerName++)
    {
      if (*desiredTriggerName == triggerNameStripped && triggerResultsHandle_->accept(iTrigger))
      {
        stringstream debugString;
        debugString << "isHltPassed:" <<endl;
        debugString << "  Trigger "<<iTrigger<<": "<< triggerName << "("<<triggerNameStripped<<") passed: "<<triggerResultsHandle_->accept(iTrigger)<<endl;
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

// same check for isHltPassed +
// check if the muon is the one firing the HLT path
bool UFDiMuonsAnalyzer::isHltMatched(const edm::Event& iEvent, 
                                     const edm::EventSetup& iSetup, 
                                     const std::string& desiredTriggerName, 
                                     const pat::TriggerObjectStandAloneCollection& triggerObjects,
                                     const pat::Muon& mu)
{
  using namespace std;
  using namespace edm;
  using namespace pat;
  using namespace reco;
  using namespace trigger;


  const boost::regex re("_v[0-9]+");

  const TriggerNames &triggerNames = iEvent.triggerNames(*triggerResultsHandle_);

  const unsigned nTriggers = triggerResultsHandle_->size();
  for (unsigned iTrigger = 0; iTrigger < nTriggers; ++iTrigger)
  {
    const string triggerName = triggerNames.triggerName(iTrigger);
    string triggerNameStripped = boost::regex_replace(triggerName,re,"",boost::match_default | boost::format_sed);
    if (desiredTriggerName == triggerNameStripped && triggerResultsHandle_->accept(iTrigger))
    {
      stringstream debugString;
      debugString << "isHltMatched: ";
      debugString << "'" << desiredTriggerName<<"'\n";
      debugString << "  Trigger "<<iTrigger<<": "<< triggerName << "("<<triggerNameStripped<<") passed: "<<triggerResultsHandle_->accept(iTrigger)<<endl;
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

bool UFDiMuonsAnalyzer::isPreselected(const pat::Muon& muon,
                                      edm::Handle<reco::BeamSpot> beamSpotHandle){

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

bool UFDiMuonsAnalyzer::passKinCuts(const pat::Muon& muon,
                                    edm::Handle<reco::BeamSpot> beamSpotHandle){
  
  bool passKinCuts=false;
  
  // =========================================================== //
  // What was corresponding to the old Loose VBTF
  // =========================================================== //
  if (_isVerbose) {
    std::cout<< "is Global?"     << muon.isGlobalMuon()     << std::endl;
    std::cout<< "is Tracker?"    << muon.isTrackerMuon()    << std::endl;
    std::cout<< "is StandAlone?" << muon.isStandAloneMuon() << std::endl;
  }

  // reconstruction cuts
  if (!muon.isGlobalMuon()    && _isGlobal )    return passKinCuts; // gbl muon
  if (!muon.isTrackerMuon()   && _isTracker)    return passKinCuts; // trk muon

  // do not accept muons which are standalone only
  if(!muon.isGlobalMuon() && !muon.isTrackerMuon()) return passKinCuts;

  reco::Track globalTrack;
  if (muon.isGlobalMuon())       globalTrack = *(muon.globalTrack());
  else if (muon.isTrackerMuon()) globalTrack = *(muon.innerTrack());
  else {
    // redundant: just in case...
    std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
    return false;
  }

  if (_isVerbose) {
    std::cout<< "muon.pt(): "  << muon.pt() << " [ptMin="  << _ptMin  <<"]" << std::endl;
    std::cout<< "fabs(muon.eta()): " << fabs(muon.eta())<< " [etaMax=" << _etaMax <<"]" << std::endl;
    std::cout<< "numberOfValidTrackerHits(): " << globalTrack.hitPattern().numberOfValidTrackerHits() 
             << " [min=" << _numValidTrackerHitsMin << "]" << std::endl;
    
    std::cout<<"d0: " << fabs(globalTrack.dxy(beamSpotHandle->position()))
             << " [max=" << _d0Max << "]" << std::endl;
    
  }

    
  // kinematic cuts
  if (muon.pt()        < _ptMin ) return passKinCuts; // pt cut
  if (fabs(muon.eta()) > _etaMax) return passKinCuts; // eta cut
  if (globalTrack.hitPattern().numberOfValidTrackerHits() < _numValidTrackerHitsMin) return passKinCuts; // # hits in tracker

  // beam spot cut
  if (fabs(globalTrack.dxy(beamSpotHandle->position())) > _d0Max) return passKinCuts;


  // =========================================================== //
  // What was corresponding to the old Tight VBTF
  // + some additions
  // =========================================================== //
  if (_isVerbose) {
    std::cout << "numberOfValidMuonHits: " 
              << globalTrack.hitPattern().numberOfValidMuonHits() 
              << " [min=" << _numValidMuonHitsMin << "]" 
              << std::endl;
  
    std::cout << "numberOfValidPixelHits(): " 
              << globalTrack.hitPattern().numberOfValidPixelHits() 
              << " [min=" << _numValidPixelHitsMin << "]" 
              << std::endl;
     
    std::cout << "numberOfValidStripHits(): " 
              << globalTrack.hitPattern().numberOfValidStripHits() 
              << " [min=" << _numValidStripHitsMin << "]" 
              << std::endl;
    
    std::cout << "pixelLayersWithMeasurement(): " 
              << globalTrack.hitPattern().pixelLayersWithMeasurement() 
              << " [min=" << _numPixelLayersMin << "]" 
              << std::endl;

    std::cout << "trackerLayersWithMeasurement(): " 
              << globalTrack.hitPattern().trackerLayersWithMeasurement() 
              << " [min=" << _numTrackerLayersMin << "]" 
              << std::endl;

    std::cout << "stripLayersWithMeasurement(): " 
              << globalTrack.hitPattern().stripLayersWithMeasurement() 
              << " [min=" << _numStripLayersMin << "]" 
              << std::endl;

    std::cout << "validFraction(): " 
              << globalTrack. validFraction()
              << " [min=" << _validFracTrackerMin << "]" 
              << std::endl;

    std::cout << "muon.numberOfMatches(): " 
              <<  muon.numberOfMatches() 
              << " [min=" << _numSegmentMatchesMin << "]" 
              << std::endl;
  
    std::cout << "muon.numberOfMatchedStations(): " 
              <<  muon.numberOfMatchedStations() 
              << " [min=" << _numOfMatchedStationsMin << "]" 
              << std::endl;
  
    std::cout<<"globalTrack.normalizedChi2(): " 
             << globalTrack.normalizedChi2()
             << " [max=" << _normChiSquareMax <<"]" 
             << std::endl;
   
  }

  if (globalTrack.hitPattern().pixelLayersWithMeasurement()   < _numPixelLayersMin)   return passKinCuts;   
  if (globalTrack.hitPattern().trackerLayersWithMeasurement() < _numTrackerLayersMin) return passKinCuts; 
  if (globalTrack.hitPattern().stripLayersWithMeasurement()   < _numTrackerLayersMin) return passKinCuts;   
  		 
  if (globalTrack.validFraction() < _validFracTrackerMin) return passKinCuts; 

  if ( globalTrack.hitPattern().numberOfValidMuonHits() < _numValidMuonHitsMin  )  return passKinCuts;
  if ( globalTrack.hitPattern().numberOfValidPixelHits() < _numValidPixelHitsMin ) return passKinCuts;
  if ( globalTrack.hitPattern().numberOfValidStripHits() < _numValidStripHitsMin ) return passKinCuts;
  if ( muon.numberOfMatches() < _numSegmentMatchesMin   )         return passKinCuts;
  if ( muon.numberOfMatchedStations() < _numOfMatchedStationsMin) return passKinCuts;
  if ( globalTrack.normalizedChi2() > _normChiSquareMax)          return passKinCuts;

  if (_isVerbose) std::cout << "passing kinematic cuts\n"; 

  passKinCuts=true;
  return passKinCuts;
}

UFDiMuonsAnalyzer::MuonPairs const UFDiMuonsAnalyzer::GetMuonPairs(pat::MuonCollection const* muons) const {
                                                                           
  
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

void UFDiMuonsAnalyzer::initMuon(_MuonInfo& muon) {

  muon.isTracker    = -999;
  muon.isStandAlone = -999;
  muon.isGlobal     = -999;

  muon.charge = -999;
  muon.pt     = -999;
  muon.eta    = -999; 
  muon.phi    = -999;
  
  muon.normChiSquare=-999;
  muon.d0_BS= -999;
  muon.dz_BS= -999;
  
  muon.d0_PV= -999;
  muon.dz_PV= -999;

  muon.numPixelLayers = -999; 
  muon.numTrackerLayers = -999;
  muon.numStripLayers = -999;  
  
  muon.validFracTracker = -999;

  muon.numValidMuonHits    = -999;
  muon.numValidPixelHits   = -999;
  muon.numValidTrackerHits = -999;
  muon.numValidStripHits   = -999;
  muon.numSegmentMatches   = -999;
  muon.numOfMatchedStations= -999;
  
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

  muon.segmentCompatibility = -999;
  muon.combinedQualityChi2LocalPosition = -999;
  muon.combinedQualityTrkKink = -999;

  for (unsigned int iTrigger=0;iTrigger<3;iTrigger++) {
    muon.isHltMatched[iTrigger] = -999;
    muon.hltPt[iTrigger] = -999;
    muon.hltEta[iTrigger] = -999;
    muon.hltPhi[iTrigger] = -999;
  }

}

void UFDiMuonsAnalyzer::initTrack(_TrackInfo& track) {
  track.charge = -999; 
  track.pt     = -999; 
  track.ptErr  = -999;
  track.eta    = -999; 
  track.phi    = -999;
}

void UFDiMuonsAnalyzer::initGenPart(_genPartInfo& part){
  part.mass = -999;
  part.pt   = -999;
  part.eta  = -999;
  part.y    = -999;
  part.phi  = -999;
}

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


void UFDiMuonsAnalyzer::displaySelection() {

  std::cout << "\n\n*** UFDiMuonsAnalyzer Configuration ***\n";

  // variable to cuts over
  std::cout << " - _ptMin:            " << _ptMin << std::endl;
  std::cout << " - _etaMax:           " << _etaMax<< std::endl;
  std::cout << " - _normChiSquareMax: " << _normChiSquareMax << std::endl;
  std::cout << " - _d0Max:            " << _d0Max << std::endl;            

  std::cout << " - _numPixelLayersMin:   " << _numPixelLayersMin  << std::endl;   
  std::cout << " - _numTrackerLayersMin: " << _numTrackerLayersMin<< std::endl; 
  std::cout << " - _numStripLayersMin:   " << _numStripLayersMin  << std::endl;   

  std::cout << " - _validFracTrackerMin: " << _validFracTrackerMin<< std::endl; 

  std::cout << " - _numValidMuonHitsMin:    " << _numValidMuonHitsMin     << std::endl;
  std::cout << " - _numValidPixelHitsMin:   " << _numValidPixelHitsMin    << std::endl;
  std::cout << " - _numValidStripHitsMin:   " << _numValidStripHitsMin    << std::endl;
  std::cout << " - _numValidTrackerHitsMin: " << _numValidTrackerHitsMin  << std::endl;
  std::cout << " - _numSegmentMatchesMin:   " << _numSegmentMatchesMin    << std::endl;  
  std::cout << " - _numOfMatchedStationsMin:" << _numOfMatchedStationsMin << std::endl;  

  std::cout << " - _trackIsoMaxSumPt: " << _trackIsoMaxSumPt << std::endl;
  std::cout << " - _relCombIsoMax: "    << _relCombIsoMax    << std::endl;      

  // module config parameters
  std::cout << " - _checkTrigger: " << _checkTrigger << std::endl;
  std::cout << " - Additional Triggers To Probe:\n";
  unsigned int triggerSize = triggerNames_.size();
  for (unsigned int i=0; i < triggerSize; i++) 
    std::cout << "    * triggerNames["<<i<<"]: " << triggerNames_[i] << std::endl;
  
  std::cout << std::endl << std::endl;

}

std::vector<float> UFDiMuonsAnalyzer::getPUJetIDDisc(edm::Handle<edm::View<pat::Jet> >  jets, const edm::Event& event, edm::InputTag tag)
{
  if(!jets.isValid())
  {
    std::vector<float> result;
    return result;
  }
  unsigned nJets = jets->size();
  std::vector<float> result(nJets,-99999999.0);
  if(tag.label()=="null")
    return result;

  edm::Handle<edm::ValueMap<float> > puJetId;
  event.getByLabel(tag,puJetId);
  if (!puJetId.isValid())
    return result;

  for(unsigned i=0; i<nJets;++i)
  {
    float id = (*puJetId)[jets->refAt(i)];
    result[i] = id;
  }
  return result;
}
std::vector<int> UFDiMuonsAnalyzer::getPUJetID(edm::Handle<edm::View<pat::Jet> >  jets, const edm::Event& event, edm::InputTag tag)
{
  if(!jets.isValid())
  {
    std::vector<int> result;
    return result;
  }
  unsigned nJets = jets->size();
  std::vector<int> result(nJets,-99999999.0);
  if(tag.label()=="null")
    return result;
  edm::Handle<edm::ValueMap<int> > puJetId;
  event.getByLabel(tag,puJetId);
  if (!puJetId.isValid())
    return result;
  for(unsigned i=0; i<nJets;++i)
  {
    int id = (*puJetId)[jets->refAt(i)];
    result[i] = id;
  }
  return result;
}

bool isMediumMuon(const reco::Muon & recoMu) 
{
   bool goodGlob = recoMu.isGlobalMuon() && 
                   recoMu.globalTrack()->normalizedChi2() < 3 && 
                   recoMu.combinedQuality().chi2LocalPosition < 12 && 
                   recoMu.combinedQuality().trkKink < 20; 
   bool isMedium = muon::isLooseMuon(recoMu) && 
                   recoMu.innerTrack()->validFraction() > 0.8 && 
                   muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
   return isMedium; 
}

DEFINE_FWK_MODULE(UFDiMuonsAnalyzer);


