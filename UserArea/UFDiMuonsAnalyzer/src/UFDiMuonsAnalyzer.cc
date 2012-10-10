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
// $Id: UFDiMuonsAnalyzer.cc,v 1.4 2012/09/19 15:05:19 digiovan Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBranch.h"

// user include files
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
#include "DataFormats/Math/interface/deltaPhi.h"
//
// trigger
// 
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Common/interface/TriggerNames.h" 
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// vertexing
//
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//
// gen particles
//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// 2010.11.21 Adding the Muon Cocktail
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/Math/interface/deltaR.h"

// pfJets and MET
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include <algorithm>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

// PU Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

bool sortGenJetFunc(reco::GenJet i, reco::GenJet j){ return (i.pt()>j.pt()); }

// Add the data formats
#include "UserArea/UFDiMuonsAnalyzer/interface/DataFormats.h"

// Add the mass error calculation from WW analysis
#include "UserArea/UFDiMuonsAnalyzer/interface/CompositeCandMassResolution.h"

// general 
double const PDG_MASS_Z  = 91.1876;//GeV/c2
double const PDG_WIDTH_Z = 2.4952; //GeV/c2

//
// class declaration
//
class UFDiMuonsAnalyzer : public edm::EDAnalyzer {

public:

  explicit UFDiMuonsAnalyzer(const edm::ParameterSet&);
  ~UFDiMuonsAnalyzer();
  
  int _numEvents;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;

  edm::ESHandle<GlobalTrackingGeometry> globalTrackingGeometry;
  //  MuonServiceProxy* theService;

  // muons pairs
  typedef std::pair<reco::Muon,reco::Muon> MuonPair;
  typedef std::vector<MuonPair>            MuonPairs;

  // sort the dimuons putting in front the candidates whose mass is 
  // the closest to the Z PDG value
  struct sortMuonsClass {
    bool operator() (MuonPair pair1, 
  		     MuonPair pair2) {

      double const MASS_MUON = 0.105658367;    //GeV/c2
      TLorentzVector muon11, muon12, dimuon1;

      reco::TrackRef const muon11Track = pair1.first . innerTrack();
      reco::TrackRef const muon12Track = pair1.second. innerTrack();

      muon11.SetPtEtaPhiM(muon11Track->pt(), muon11Track->eta(), muon11Track->phi(), MASS_MUON);
      muon12.SetPtEtaPhiM(muon12Track->pt(), muon12Track->eta(), muon12Track->phi(), MASS_MUON);

      dimuon1 = muon11+muon12;


      TLorentzVector muon21, muon22, dimuon2;

      reco::TrackRef const muon21Track = pair2.first . innerTrack();
      reco::TrackRef const muon22Track = pair2.second. innerTrack();

      muon21.SetPtEtaPhiM(muon21Track->pt(), muon21Track->eta(), muon21Track->phi(), MASS_MUON);
      muon22.SetPtEtaPhiM(muon22Track->pt(), muon22Track->eta(), muon22Track->phi(), MASS_MUON);

      dimuon2 = muon21+muon22;

      return ( fabs(dimuon1.M()-PDG_MASS_Z) < fabs(dimuon2.M()-PDG_MASS_Z) );
    }
  } sortMuonObject;

  // tracks pairs, e.g. cocktail
  typedef std::pair<reco::Track,reco::Track> TrackPair;
  typedef std::vector<TrackPair> TrackPairs;

  // general event information	
  _EventInfo eventInfo;

  // vertex information
  _VertexInfo vertexInfo;


  // combined muon-muon info
  float _recoCandMass;
  float _recoCandMassRes;
  float _recoCandMassResCov;
  float _recoCandPt;
  float _recoCandEta; // pseudo rapidity
  float _recoCandY;   // rapidity
  float _recoCandPhi; // phi

  // combined muon-muon info
  float _recoCandMassPF;
  float _recoCandMassResPF;
  float _recoCandPtPF;
  float _recoCandEtaPF; 
  float _recoCandYPF;   
  float _recoCandPhiPF; 


  float _angleDiMuons;
  int   _vertexIsValid;
  float _vertexNormChiSquare;
  float _vertexChiSquare;
  float _vertexNDF;
  float _vertexX;
  float _vertexY;
  float _vertexZ;


  _MuonInfo _muon1, _muon2;

  // dimuon candidates from refit with vertex constraint
  _TrackInfo _muon1vc,_muon2vc;

  // Vertex Constrained muon-muon info
  float _recoCandMassVC;
  float _recoCandMassResVC;
  float _recoCandMassResCovVC;
  float _recoCandPtVC;
  float _recoCandEtaVC; 
  float _recoCandYVC;   
  float _recoCandPhiVC; 


  // generator level info post-FSR
  float _trueMass;
  float _truePt;
  float _trueEta; // pseudo rapidity
  float _trueY;   // rapidity
  float _truePhi; // phi

  _TrackInfo _true1,_true2;

  // generator level info pre-FSR
  float _trueMassPreFSR;
  float _truePtPreFSR;
  float _trueEtaPreFSR; // pseudo rapidity
  float _trueYPreFSR;   // rapidity
  float _truePhiPreFSR; // phi

  _TrackInfo _true1PreFSR,_true2PreFSR;

  // Jets and MET
  _MetInfo     _metInfo;
  _PFJetInfo   _pfJetInfo;
  _GenJetInfo  _genJetInfo;

  float _puJetFullDisc[10];
  float _puJetFullId[10];
  float _puJetSimpleDisc[10];
  float _puJetSimpleId[10];
  float _puJetCutDisc[10];
  float _puJetCutId[10];

  int _nPU;

  // where to save all the info  
  TFile* _outFile;	
  TTree* _outTree;


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void initMuon ( _MuonInfo& muon );
  void initTrack(_TrackInfo& track);

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);

  // method to select the muons
  bool isHltPassed (const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
  bool isHltMatched(const edm::Event&, const edm::EventSetup&, 
                    const std::string& triggerName, const std::string& filterName,
                    const reco::Muon&,
                    float& hltPt, float& hltEta, float& hltPhi);
  bool isHltMatched(const edm::Event&, const edm::EventSetup&, 
                    const std::vector < std::string>& triggerNames, 
                    const reco::Muon&, const reco::Muon&);
  bool isDoubleMu(const std::string& triggerName);    
  void addVersion(const HLTConfigProvider hltConfig_,
                  std::string& triggerBaseName,
                  std::string& triggerName);    
  void findLowestSingleMu(const HLTConfigProvider hltConfig_,
                          const edm::Event&, const edm::EventSetup&);
  std::string findFilterName(const HLTConfigProvider hltConfig_, const std::string& triggerName);
  double DR(double eta1, double eta2, double phi1, double phi2);

  bool isPreselected(const reco::Muon& muon,
                     edm::Handle<reco::BeamSpot> beamSpotHandle);
  
  bool passKinCuts(const reco::Muon& muon,
                   edm::Handle<reco::BeamSpot> beamSpotHandle);

  void displaySelection();

  // WW analysis mass resolution
  CompositeCandMassResolution  CCMassResolution;

  // muons
  edm::InputTag _muonColl;
  edm::InputTag _beamSpot;		
  std::string _getFilename;	

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
  std::vector < std::string > triggerBaseNames_;
  std::vector < std::string > filterNames_;

  std::vector < int > l1Prescale_;
  std::vector < int > hltPrescale_;

  edm::InputTag triggerResultsTag_;
  edm::InputTag triggerEventTag_;

  bool _selectLowestSingleMuTrigger;

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
  edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
  HLTConfigProvider hltConfig_;


  MuonPairs  const GetMuonPairs (reco::MuonCollection  const* muons ) const;
  TrackPairs const GetTrackPairs(reco::TrackCollection const* tracks) const;

  TLorentzVector const GetLorentzVector(UFDiMuonsAnalyzer::MuonPair  const* pair) ;//const; 
  TLorentzVector const GetLorentzVector(UFDiMuonsAnalyzer::TrackPair const* pair) ;//const; 

  double const GetMassRes(UFDiMuonsAnalyzer::MuonPair  const* pair);
  double const GetMassRes(UFDiMuonsAnalyzer::TrackPair const* pair);

  double const GetMassResCov(UFDiMuonsAnalyzer::MuonPair  const* pair);
  double const GetMassResCov(UFDiMuonsAnalyzer::TrackPair  const* pair);

  TransientVertex const GetVertexFromPair(UFDiMuonsAnalyzer::TrackPair const* muPair) const;
  TransientVertex const GetVertexFromTracks(reco::TrackRef trackref1, reco::TrackRef trackref2) const;


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

  //now do what ever initialization is needed
  _getFilename  = iConfig.getUntrackedParameter<std::string>("getFilename", "badger.root");
  _muonColl	= iConfig.getParameter<edm::InputTag>("muonColl");

  _beamSpot	= iConfig.getUntrackedParameter<edm::InputTag>("beamSpotTag",
                                                               edm::InputTag("offlineBeamSpot") );

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

  //HLT trigger initialization
  _checkTrigger	     = iConfig.getParameter<bool>("checkTrigger");

  processName_       = iConfig.getParameter<std::string>("processName");
  triggerBaseNames_  = iConfig.getParameter<std::vector <std::string> >("triggerNames");
  _selectLowestSingleMuTrigger = iConfig.getUntrackedParameter<bool>("selectLowestSingleMuTrigger", false);

  if (triggerBaseNames_[0]=="") {
    triggerBaseNames_.clear();
    triggerNames_.clear();
    l1Prescale_.clear();
    hltPrescale_.clear();
  }

  if (_selectLowestSingleMuTrigger) triggerBaseNames_.push_back("trigPlaceHolder");

  for (unsigned int i=0; i<triggerBaseNames_.size();i++) {
    triggerNames_.push_back("HLT_placeHolder");
    filterNames_.push_back("filter_placeHolder");
    l1Prescale_.push_back(-999);
    hltPrescale_.push_back(-999);
  }

  if (triggerNames_.size() > 3){
    std::cout << "ERROR: triggerNames has its maximum size to 3! -> change the code in case\n";
    assert(0);
  }

  if (triggerNames_.size() == 0){
    std::cout << "ERROR: You need to pass at least one valid HLT trigger\n";
    assert(0);
  }

  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResults");
  triggerEventTag_   = iConfig.getParameter<edm::InputTag>("triggerEvent");

  metTag     = iConfig.getParameter<edm::InputTag>("metTag");
  pfJetsTag  = iConfig.getParameter<edm::InputTag>("pfJetsTag");
  genJetsTag = iConfig.getParameter<edm::InputTag>("genJetsTag");

  puJetMvaFullDiscTag 	 = iConfig.getParameter<edm::InputTag>("puJetMvaFullDiscTag");
  puJetMvaFullIdTag 	 = iConfig.getParameter<edm::InputTag>("puJetMvaFullIdTag");
  puJetMvaSimpleDiscTag	 = iConfig.getParameter<edm::InputTag>("puJetMvaSimpleDiscTag");
  puJetMvaSimpleIdTag 	 = iConfig.getParameter<edm::InputTag>("puJetMvaSimpleIdTag");
  puJetMvaCutDiscTag 	 = iConfig.getParameter<edm::InputTag>("puJetMvaCutDiscTag");
  puJetMvaCutIdTag 	 = iConfig.getParameter<edm::InputTag>("puJetMvaCutIdTag");

}


UFDiMuonsAnalyzer::~UFDiMuonsAnalyzer() {}


// ------------ method called to for each event  ------------
void UFDiMuonsAnalyzer::analyze(const edm::Event& iEvent, 
                                const edm::EventSetup& iSetup)
{
  _numEvents++;

  if (_isVerbose) 
    std::cout << "\n\n A N A LI Z I N G   E V E N T = " 
	      << _numEvents << std::endl << std::endl;
 

  // -----------------------------------------
  // H L T   H A N D L E S
  iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerResults product from Event!" << std::endl;
    return;
  }
  iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    std::cout << "UFDiMuonsAnalyzer::analyze: Error in getting TriggerEvent product from Event!" << std::endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());
  
  if (_selectLowestSingleMuTrigger)
    findLowestSingleMu(hltConfig_,iEvent,iSetup);

  // if all the HLT paths are not fired, will discard the event immediately
  if (_checkTrigger) {
    bool isPassed=false;
     
    for (unsigned int iTrigger=0; iTrigger<triggerNames_.size(); iTrigger++) 
      if ( isHltPassed(iEvent,iSetup,triggerNames_[iTrigger]) ) isPassed=true;
     
    if (!isPassed) {
      if (_isVerbose) std::cout << "None of the HLT paths fired -> discard the event\n";
      return;
    }
     
  }


  // once the triggernames are cleared add the prescales
  for (unsigned int iTrigger=0; iTrigger<triggerNames_.size(); iTrigger++) {
    std::pair<int,int> prescaleValues = hltConfig_.prescaleValues(iEvent, iSetup,triggerNames_[iTrigger]);
    l1Prescale_[iTrigger]  = prescaleValues.first;
    hltPrescale_[iTrigger] = prescaleValues.second;

  }



  int theRun   = iEvent.id().run();
  int theLumi  = iEvent.luminosityBlock();
  int theEvent = iEvent.id().event();
  int theBx    = iEvent.bunchCrossing();
  int theOrbit = iEvent.orbitNumber();
  
  eventInfo.run   = theRun;
  eventInfo.lumi  = theLumi;
  eventInfo.event = theEvent;
  eventInfo.bx    = theBx;
  eventInfo.orbit = theOrbit;


  //vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlinePrimaryVertices", vertices);
 
  for (int i=0;i<80;i++) {
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
  
  // Get MC Truth Pileup
  _nPU = -1;
  if (_isMonteCarlo) {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

       int BX = PVI->getBunchCrossing();

       if(BX == 0) { 
          _nPU = PVI->getTrueNumInteractions();
          continue;
       }
    }
  }
  


  // G E O M E T R Y
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transientTrackBuilder);
  // Comp Cand Mass Res
  CCMassResolution.init(iSetup);

  //
  // get the sim mass if it exists
  //
  if (_isMonteCarlo) {

    edm::Handle<reco::GenParticleCollection> allGenParticles;
    iEvent.getByLabel("genParticles", allGenParticles);
    ////const reco::GenParticleCollection* genParticles  = allGenParticles.failedToGet () ? 0 : &*allGenParticles;
    
    
    reco::GenParticleCollection theGens;
    reco::GenParticleCollection theGensPreFSR;

    //std::cout << "event" << std::endl;
    for (reco::GenParticleCollection::const_iterator gen = allGenParticles->begin(), genEnd = allGenParticles->end(); 
         gen != genEnd; ++gen) 
      {
        ////for (size_t i=0; i<genParticles->size(); ++i){
        ////  
        ////  const reco::Candidate &gen = (*genParticles)[i];
        ////  int id       = gen.pdgId();
        ////  int parentId = -999;
        ////  if (fabs(id) != 13 && fabs(id) != 22) continue;
        ////  
        ////  std::cout << id << ", " << gen.status() << ", " << gen.pt();
        ////  if (gen.numberOfMothers())
        ////    std::cout << "has a mother?" << gen.mother()->pdgId() << std::endl;
        ////  else
        ////    std::cout << std::endl;
        
        int id       = gen->pdgId();
        int parentId = -999;
        
        if (abs(id) == 13) parentId = gen->mother()->pdgId();
        
        //if (abs(id) == 13)  {
        //  std::cout << "muon, " << gen->status() << ", " << gen->pt() << std::endl;
        //  std::cout << "has a mother?" << gen->mother()->pdgId() << std::endl;
        //}
        //if (abs(id) == 23)  std::cout << "Z " << gen->status() << ", " << gen->mass() << ", " << gen->pt() << std::endl;
        
        if (gen->status() == 3 && abs(id) == 13 && abs(parentId) == 23) theGensPreFSR.push_back(*gen);
        if (gen->status() == 1 && abs(id) == 13 && abs(parentId) == 13) theGens.push_back(*gen);
      
      }
    
    _trueMassPreFSR = -999;
    _truePtPreFSR   = -999;
    
    _trueEtaPreFSR = -999;
    _trueYPreFSR   = -999;
    _truePhiPreFSR = -999;
    
    _true1PreFSR.charge = -999; _true1PreFSR.pt = -999; _true1PreFSR.ptErr = -999;_true1PreFSR.eta = -999; _true1PreFSR.phi = -999;
    _true2PreFSR.charge = -999; _true2PreFSR.pt = -999; _true2PreFSR.ptErr = -999;_true2PreFSR.eta = -999; _true2PreFSR.phi = -999;

    if (theGensPreFSR.size() == 2) { 
      TLorentzVector vtrue1, vtrue2, vtrueMother;
      
      vtrue1.SetPtEtaPhiM(theGensPreFSR[0].pt(), theGensPreFSR[0].eta(), theGensPreFSR[0].phi(), theGensPreFSR[0].mass());
      vtrue2.SetPtEtaPhiM(theGensPreFSR[1].pt(), theGensPreFSR[1].eta(), theGensPreFSR[1].phi(), theGensPreFSR[1].mass());

      _true1PreFSR.charge = theGensPreFSR[0].charge(); _true1PreFSR.pt = theGensPreFSR[0].pt(); 
      _true2PreFSR.charge = theGensPreFSR[1].charge(); _true2PreFSR.pt = theGensPreFSR[1].pt(); 
      _true1PreFSR.eta = theGensPreFSR[0].eta(); _true1PreFSR.phi=theGensPreFSR[0].phi();	
      _true2PreFSR.eta = theGensPreFSR[1].eta(); _true2PreFSR.phi=theGensPreFSR[1].phi();
      
      vtrueMother = vtrue1+vtrue2;	

      _trueMassPreFSR = vtrueMother.M();
      _truePtPreFSR   = vtrueMother.Pt();

      _trueEtaPreFSR = vtrueMother.PseudoRapidity();
      _trueYPreFSR   = vtrueMother.Rapidity();
      _truePhiPreFSR = vtrueMother.Phi();

    }


    _trueMass = -999;
    _truePt   = -999;

    _trueEta = -999;
    _trueY   = -999;
    _truePhi = -999;

    _true1.charge = -999; _true1.pt = -999; _true1.ptErr = -999;_true1.eta = -999; _true1.phi = -999;
    _true2.charge = -999; _true2.pt = -999; _true2.ptErr = -999;_true2.eta = -999; _true2.phi = -999;

    if (theGens.size() == 2) { 
      TLorentzVector vtrue1, vtrue2, vtrueMother;
      vtrue1.SetPtEtaPhiM(theGens[0].pt(), theGens[0].eta(), theGens[0].phi(), theGens[0].mass());
      vtrue2.SetPtEtaPhiM(theGens[1].pt(), theGens[1].eta(), theGens[1].phi(), theGens[1].mass());
      _true1.charge = theGens[0].charge(); _true1.pt = theGens[0].pt(); _true1.eta = theGens[0].eta(); _true1.phi=theGens[0].phi();
      _true2.charge = theGens[1].charge(); _true2.pt = theGens[1].pt(); _true2.eta = theGens[1].eta(); _true2.phi=theGens[1].phi();
	
      vtrueMother = vtrue1+vtrue2;	

      _trueMass = vtrueMother.M();
      _truePt   = vtrueMother.Pt();

      _trueEta = vtrueMother.PseudoRapidity();
      _trueY   = vtrueMother.Rapidity();
      _truePhi = vtrueMother.Phi();

    }

  }

  // ===========================================================================
  // M U O N S
  //
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(_muonColl, muons);

  // B E A M S P O T
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(_beamSpot, beamSpotHandle);

  // reco muons collection 
  reco::MuonCollection muonsSelected; 
  muonsSelected.clear(); 


  // pre-selection: just check the muons are at least tracker muons.
  for (reco::MuonCollection::const_iterator muon = muons->begin(), 
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

  _recoCandMass = -999; 
  _recoCandMassRes = -999;
  _recoCandMassResCov = -999;
  _recoCandPt   = -999;
  _recoCandEta  = -999;
  _recoCandY    = -999;
  _recoCandPhi  = -999;

  _recoCandMassPF = -999; 
  _recoCandMassResPF = -999;
  _recoCandPtPF   = -999;
  _recoCandEtaPF  = -999;
  _recoCandYPF    = -999;
  _recoCandPhiPF  = -999;

  _recoCandMassVC = -999; 
  _recoCandMassResVC = -999;
  _recoCandMassResCovVC = -999;
  _recoCandPtVC   = -999;
  _recoCandEtaVC  = -999;
  _recoCandYVC    = -999;
  _recoCandPhiVC  = -999;


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
    reco::Muon mu = muonsSelected.at(0);
    // store all the info
    // muon 1
    _muon1.isGlobal     = mu.isGlobalMuon(); 
    _muon1.isTracker    = mu.isTrackerMuon(); 
    _muon1.isStandAlone = mu.isStandAloneMuon(); 

    reco::TrackRef trackref;
    if      (mu.isGlobalMuon())  trackref = mu.globalTrack();
    else if (mu.isTrackerMuon()) trackref = mu.innerTrack();
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon1.charge = mu.charge();
    _muon1.pt     = mu.pt();  
    _muon1.ptErr  = trackref->ptError(); 
    _muon1.eta    = mu.eta(); 
    _muon1.phi    = mu.phi();
   
    // redundant if the muon is tracker-only muon
    if (mu.isTrackerMuon()) {
      _muon1.trkPt   = mu.innerTrack()->pt();                    
      _muon1.trkPtErr= mu.innerTrack()->ptError();
      _muon1.trketa  = mu.innerTrack()->eta();                   
      _muon1.trkPhi  = mu.innerTrack()->phi();                   
    }

    _muon1.normChiSquare=trackref->normalizedChi2();
    _muon1.d0_BS= mu.innerTrack()->dxy( beamSpotHandle->position() );
    _muon1.dz_BS= mu.innerTrack()->dz ( beamSpotHandle->position() );

    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon1.d0_PV= trackref->dxy( vtx->position() );
      _muon1.dz_PV= trackref->dz ( vtx->position() );
    
      //exit at the first available vertex
      break;
    }

    _muon1.numPixelLayers   = trackref->hitPattern().pixelLayersWithMeasurement();   
    _muon1.numTrackerLayers = trackref->hitPattern().trackerLayersWithMeasurement(); 
    _muon1.numStripLayers   = trackref->hitPattern().stripLayersWithMeasurement();   
    
    _muon1.validFracTracker = trackref->validFraction(); 

    _muon1.numValidMuonHits    = trackref->hitPattern().numberOfValidMuonHits();
    _muon1.numValidPixelHits   = trackref->hitPattern().numberOfValidPixelHits();
    _muon1.numValidTrackerHits = trackref->hitPattern().numberOfValidTrackerHits();
    _muon1.numValidStripHits   = trackref->hitPattern().numberOfValidStripHits();
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


    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++) 
      _muon1.isHltMatched[iTrigger] = ( isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], filterNames_[iTrigger], mu,
                                                     _muon1.hltPt[iTrigger], _muon1.hltEta[iTrigger], _muon1.hltPhi[iTrigger]) );
    
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
    
    reco::Muon mu1  = pair->first;
    reco::Muon mu2  = pair->second;

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

    // chek the trigger
    if (_checkTrigger && passSelection) {
        
      // set the event as not passed  
      // and then check it passes the trigger
      passSelection=!true;

      if ( isHltMatched(iEvent,iSetup,triggerNames_, mu1, mu2) ) {
        if (_isVerbose) std::cout << "Both Muons TIGHT and At Least One Matches the Trigger\n";
        passSelection=true;
      }
    }
    // ===== END Requirements ====
    
    
    // do we pass the selection? Yes -> store the event
    // else move to the next event
    if (!passSelection) return;

    // store all the info
    // muon 1
    _muon1.isGlobal     = mu1.isGlobalMuon(); 
    _muon1.isTracker    = mu1.isTrackerMuon(); 
    _muon1.isStandAlone = mu1.isStandAloneMuon(); 


    reco::TrackRef trackref1;
    if      (mu1.isGlobalMuon())     trackref1 = mu1.globalTrack();
    else if (mu1.isTrackerMuon())    trackref1 = mu1.innerTrack();
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon1.charge = mu1.charge();
    _muon1.pt     = mu1.pt();  
    _muon1.ptErr  = trackref1->ptError(); 
    _muon1.eta    = mu1.eta(); 
    _muon1.phi    = mu1.phi();
   
    if (mu1.isTrackerMuon()) {
      _muon1.trkPt   = mu1.innerTrack()->pt();                    
      _muon1.trkPtErr= mu1.innerTrack()->ptError();
      _muon1.trketa  = mu1.innerTrack()->eta();                   
      _muon1.trkPhi  = mu1.innerTrack()->phi();                   
    }

    _muon1.normChiSquare=trackref1->normalizedChi2();
    _muon1.d0_BS= trackref1->dxy( beamSpotHandle->position() );
    _muon1.dz_BS= trackref1->dz ( beamSpotHandle->position() );

    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon1.d0_PV= trackref1->dxy( vtx->position() );
      _muon1.dz_PV= trackref1->dz ( vtx->position() );
    
      //exit at the first available vertex
      break;
    }

    _muon1.numPixelLayers   = trackref1->hitPattern().pixelLayersWithMeasurement();   
    _muon1.numTrackerLayers = trackref1->hitPattern().trackerLayersWithMeasurement(); 
    _muon1.numStripLayers   = trackref1->hitPattern().stripLayersWithMeasurement();   
    
    _muon1.validFracTracker = trackref1->validFraction(); 

    _muon1.numValidMuonHits    = trackref1->hitPattern().numberOfValidMuonHits();
    _muon1.numValidPixelHits   = trackref1->hitPattern().numberOfValidPixelHits();
    _muon1.numValidTrackerHits = trackref1->hitPattern().numberOfValidTrackerHits();
    _muon1.numValidStripHits   = trackref1->hitPattern().numberOfValidStripHits();
    _muon1.numSegmentMatches   = mu1.numberOfMatches();
    _muon1.numOfMatchedStations= mu1.numberOfMatchedStations();

    //tracker iso and rel comb iso already taken care of
    _muon1.ecalIso = mu1.isolationR03().emEt ;
    _muon1.hcalIso = mu1.isolationR03().hadEt ;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++) 
      _muon1.isHltMatched[iTrigger] = ( isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], filterNames_[iTrigger], mu1,
                                                     _muon1.hltPt[iTrigger], _muon1.hltEta[iTrigger], _muon1.hltPhi[iTrigger]) );
    
    // muon 2
    _muon2.isGlobal     = mu2.isGlobalMuon(); 
    _muon2.isTracker    = mu2.isTrackerMuon(); 

    reco::TrackRef trackref2;
    if      (mu2.isGlobalMuon())     trackref2 = mu2.globalTrack();
    else if (mu2.isTrackerMuon())    trackref2 = mu2.innerTrack();
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return ;
    }

    _muon2.charge = mu2.charge(); 
    _muon2.pt     = mu2.pt(); 
    _muon2.ptErr  = trackref2->ptError(); 
    _muon2.eta    = mu2.eta(); 
    _muon2.phi    = mu2.phi();
   
    if (mu2.isTrackerMuon()) {
      _muon2.trkPt   = mu2.innerTrack()->pt();                    
      _muon2.trkPtErr= mu2.innerTrack()->ptError();
      _muon2.trketa  = mu2.innerTrack()->eta();                   
      _muon2.trkPhi  = mu2.innerTrack()->phi();                   
    }

    _muon2.numPixelLayers   = trackref2->hitPattern().pixelLayersWithMeasurement();   
    _muon2.numTrackerLayers = trackref2->hitPattern().trackerLayersWithMeasurement(); 
    _muon2.numStripLayers   = trackref2->hitPattern().stripLayersWithMeasurement();   
    
    _muon2.validFracTracker = trackref2->validFraction(); 

    _muon2.normChiSquare=trackref2->normalizedChi2();
    _muon2.d0_BS= trackref2->dxy( beamSpotHandle->position() );
    _muon2.dz_BS= trackref2->dz ( beamSpotHandle->position() );
   
    for (reco::VertexCollection::const_iterator vtx = vertices->begin(); 
         vtx!=vertices->end(); ++vtx){

      if (!vtx->isValid()) continue;
  
      _muon2.d0_PV= trackref2->dxy( vtx->position() );
      _muon2.dz_PV= trackref2->dz ( vtx->position() );
    
      //exit at the first available vertex
      break;
    }

    _muon2.numValidMuonHits    = trackref2->hitPattern().numberOfValidMuonHits();
    _muon2.numValidPixelHits   = trackref2->hitPattern().numberOfValidPixelHits();
    _muon2.numValidTrackerHits = trackref2->hitPattern().numberOfValidTrackerHits();
    _muon2.numValidStripHits   = trackref2->hitPattern().numberOfValidStripHits();
    _muon2.numSegmentMatches   = mu2.numberOfMatches();
    _muon2.numOfMatchedStations= mu2.numberOfMatchedStations();

    //tracker iso and rel comb iso already taken care of
    _muon2.ecalIso = mu2.isolationR03().emEt ;
    _muon2.hcalIso = mu2.isolationR03().hadEt ;

    // save trigger informations
    for (unsigned int iTrigger=0;iTrigger<triggerNames_.size();iTrigger++) 
      _muon2.isHltMatched[iTrigger] = ( isHltMatched(iEvent, iSetup, triggerNames_[iTrigger], filterNames_[iTrigger], mu2,
                                                     _muon2.hltPt[iTrigger], _muon2.hltEta[iTrigger], _muon2.hltPhi[iTrigger]) );
    
    // combine the info
    // muons collection
    TLorentzVector mother=GetLorentzVector(&*pair); 

    _recoCandMass = mother.M(); 
    _recoCandMassRes    = GetMassRes(&*pair);
    _recoCandMassResCov = GetMassResCov(&*pair);
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

    // calculate the vtx from the dimuon candidates
    TransientVertex vtx	= GetVertexFromTracks(mu1.innerTrack(),mu2.innerTrack());
    reco::Track refitTrack1 = *mu1.innerTrack();
    reco::Track refitTrack2 = *mu2.innerTrack();

    _vertexIsValid = vtx.isValid();
    if (vtx.isValid()) {
      
      _vertexNormChiSquare = vtx.normalisedChiSquared();
      _vertexChiSquare = vtx.totalChiSquared();
      _vertexNDF = vtx.degreesOfFreedom();
      _vertexX = vtx.position().x();
      _vertexY = vtx.position().y();
      _vertexZ = vtx.position().z();

      // refit the candidates (only if the vertex is valid ;-) )
      refitTrack1 = vtx.refittedTracks().front().track();
      refitTrack2 = vtx.refittedTracks().back().track();
    }

    
    // combine the refit info into a dimuon candidate
    TrackPair refitTrackPair(refitTrack1,refitTrack2);
    
    TLorentzVector motherRefit=GetLorentzVector(&refitTrackPair); 
    
    _recoCandMassVC = motherRefit.M(); 
    _recoCandMassResVC    = GetMassRes   (&refitTrackPair);
    _recoCandMassResCovVC = GetMassResCov(&refitTrackPair);
    _recoCandPtVC   = motherRefit.Pt();
    
    _recoCandEtaVC  = motherRefit.PseudoRapidity();
    _recoCandYVC    = motherRefit.Rapidity();
    _recoCandPhiVC  = motherRefit.Phi();

    // muon 1 refit info
    _muon1vc.charge = refitTrack1.charge(); 
    _muon1vc.pt     = refitTrack1.pt(); 
    _muon1vc.ptErr  = refitTrack1.ptError(); 
    _muon1vc.eta    = refitTrack1.eta(); 
    _muon1vc.phi    = refitTrack1.phi();

    _muon2vc.charge = refitTrack2.charge(); 
    _muon2vc.pt     = refitTrack2.pt(); 
    _muon2vc.ptErr  = refitTrack2.ptError(); 
    _muon2vc.eta    = refitTrack2.eta(); 
    _muon2vc.phi    = refitTrack2.phi();


//     if (_recoCandMass > 60 && _recoCandMass < 160) {

//       std::cout << "\n\n ======= dimuon candidate ======= \n\n";
//       std::cout <<    " Mass(mm): " << _recoCandMass   
//                 << "+/-" << _recoCandMassRes
//                 << " or +/-" << _recoCandMassResCov
//                 << ", MassVC(mm): " << _recoCandMassVC 
//                 << "+/-" << _recoCandMassResVC  
//                 << " or +/-" << _recoCandMassResCovVC
//                 << std::endl; 
      
//       std::cout <<    " pT(mm): " << _recoCandPt 
//                 << ", pTVC(mm): " << _recoCandPtVC 
//                 << std::endl; 

//       std::cout << " * 1st MUON:\n";
//       std::cout << "  charge: "<<_muon1.charge<< ", chargeVC: "<<_muon1vc.charge 
//                 << "\n    pt: "<<_muon1.pt    << ",     ptVC: "<<_muon1vc.pt     
//                 << "\n ptErr: "<<_muon1.ptErr << ",  ptErrVC: "<<_muon1vc.ptErr  
//                 << "\n   eta: "<<_muon1.eta   << ",    etaVC: "<<_muon1vc.eta    
//                 << "\n   phi: "<<_muon1.phi   << ",    phiVC: "<<_muon1vc.phi    
//                 << std::endl;

//       std::cout << " * 2nd MUON:\n";
//       std::cout << "  charge: "<<_muon2.charge<< ", chargeVC: "<<_muon2vc.charge 
//                 << "\n    pt: "<<_muon2.pt    << ",     ptVC: "<<_muon2vc.pt     
//                 << "\n ptErr: "<<_muon2.ptErr << ",  ptErrVC: "<<_muon2vc.ptErr  
//                 << "\n   eta: "<<_muon2.eta   << ",    etaVC: "<<_muon2vc.eta    
//                 << "\n   phi: "<<_muon2.phi   << ",    phiVC: "<<_muon2vc.phi    
//                 << std::endl;

//     }

    // ===========================================================================
    // Jet Info
    edm::Handle < std::vector<reco::PFMET> > mets;
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

    // Get JEC Uncertainty Calculator
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    JetCorrectionUncertainty *jecUncCalculator = new JetCorrectionUncertainty(JetCorPar);

    if( jets.isValid() ){
      //First Get PU Id's

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
          jecUncCalculator->setJetEta(jet.eta());
          jecUncCalculator->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
          _pfJetInfo.jecUnc[i] = jecUncCalculator->getUncertainty(true);
          _pfJetInfo.jecFactor[i]  = jet.jecFactor("Uncorrected");
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
	  
        }
      }
    }

    edm::Handle < edm::View<pat::Jet> > jetsForPUId;
    if( pfJetsTag.label() != "null" ) iEvent.getByLabel(pfJetsTag, jetsForPUId);
    std::vector<float> puIdFullDisc = getPUJetIDDisc(jetsForPUId,iEvent,puJetMvaFullDiscTag);
    std::vector<float> puIdSimpleDisc = getPUJetIDDisc(jetsForPUId,iEvent,puJetMvaSimpleDiscTag);
    std::vector<float> puIdCutDisc = getPUJetIDDisc(jetsForPUId,iEvent,puJetMvaCutDiscTag);
    std::vector<int> puIdFullId = getPUJetID(jetsForPUId,iEvent,puJetMvaFullIdTag);
    std::vector<int> puIdSimpleId = getPUJetID(jetsForPUId,iEvent,puJetMvaSimpleIdTag);
    std::vector<int> puIdCutId = getPUJetID(jetsForPUId,iEvent,puJetMvaCutIdTag);
    for(unsigned i=0; i<10;i++)
    {
      if(i<puIdFullDisc.size())
	_puJetFullDisc[i] = puIdFullDisc[i];
      else
	_puJetFullDisc[i] = -1000.0;

      if(i<puIdFullId.size())
	_puJetFullId[i] = puIdFullId[i];
      else
	_puJetFullId[i] = -1000.0;

      if(i<puIdSimpleDisc.size())
	_puJetSimpleDisc[i] = puIdSimpleDisc[i];
      else
	_puJetSimpleDisc[i] = -1000.0;

      if(i<puIdSimpleId.size())
	_puJetSimpleId[i] = puIdSimpleId[i];
      else
	_puJetSimpleId[i] = -1000.0;

      if(i<puIdCutDisc.size())
	_puJetCutDisc[i] = puIdCutDisc[i];
      else
	_puJetCutDisc[i] = -1000.0;

      if(i<puIdCutId.size())
	_puJetCutId[i] = puIdCutId[i];
      else
	_puJetCutId[i] = -1000.0;
    }

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

  // create the file
  _outFile	= new TFile(_getFilename.c_str(), "recreate");	
  _outFile->cd();

  // create the tree
  _outTree	= new TTree("tree", "myTree");


  // eventInfo;
  _outTree->Branch( "eventInfo",  &eventInfo, "run/I:lumi/I:event/I:bx/I:orbit/I");
  
  std::cout << "beginJob" << std::endl;
  std::cout << "triggerBaseNames_.size()= " << triggerBaseNames_.size() << std::endl;


  _outTree->Branch("vertexInfo", &vertexInfo, "nVertices/I:isValid[80]/I:"
		   "x[80]/F:y[80]/F:z[80]/F:xErr[80]/F:yErr[80]/F:zErr[80]/F:"
		   "chi2[80]/F:ndf[80]/I:normChi2[80]/F");

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
                   "hltPhi[3]/F");

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
                   "hltPhi[3]/F");

  _outTree->Branch("hltPaths",    &triggerNames_);
  _outTree->Branch("filterNames", &filterNames_ );
  _outTree->Branch("l1Prescale",  &l1Prescale_  );
  _outTree->Branch("hltPrescale", &hltPrescale_ );
 

  // mass and pt for the fit
  _outTree->Branch("recoCandMass",   &_recoCandMass,   "recoCandMass/F");
  _outTree->Branch("recoCandMassRes",&_recoCandMassRes,"recoCandMassRes/F");
  _outTree->Branch("recoCandMassResCov",&_recoCandMassResCov,"recoCandMassResCov/F");
  _outTree->Branch("recoCandPt"  , &_recoCandPt  , "recoCandPt/F");
  _outTree->Branch("recoCandEta" , &_recoCandEta , "recoCandEta/F");
  _outTree->Branch("recoCandY"   , &_recoCandY   , "recoCandY/F");
  _outTree->Branch("recoCandPhi" , &_recoCandPhi , "recoCandPhi/F");

  _outTree->Branch("recoCandMassPF", &_recoCandMassPF,   "recoCandMassPF/F");
  _outTree->Branch("recoCandMassResPF", &_recoCandMassResPF,"recoCandMassResPF/F");
  _outTree->Branch("recoCandPtPF"  , &_recoCandPtPF  , "recoCandPtPF/F");
  _outTree->Branch("recoCandEtaPF" , &_recoCandEtaPF , "recoCandEtaPF/F");
  _outTree->Branch("recoCandYPF"   , &_recoCandYPF   , "recoCandYPF/F");
  _outTree->Branch("recoCandPhiPF" , &_recoCandPhiPF , "recoCandPhiPF/F");

  _outTree->Branch("reco1vc", &_muon1vc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
  _outTree->Branch("reco2vc", &_muon2vc, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

  _outTree->Branch("recoCandMassVC",      &_recoCandMassVC,      "recoCandMassVC/F");
  _outTree->Branch("recoCandMassResVC",   &_recoCandMassResVC,   "recoCandMassResVC/F");
  _outTree->Branch("recoCandMassResCovVC",&_recoCandMassResCovVC,"recoCandMassResCovVC/F");
  _outTree->Branch("recoCandPtVC"  , &_recoCandPtVC  , "recoCandPtVC/F");
  _outTree->Branch("recoCandEtaVC" , &_recoCandEtaVC , "recoCandEtaVC/F");
  _outTree->Branch("recoCandYVC"   , &_recoCandYVC   , "recoCandYVC/F");
  _outTree->Branch("recoCandPhiVC" , &_recoCandPhiVC , "recoCandPhiVC/F");

  _outTree->Branch("angleDiMuons",        &_angleDiMuons       ,"angleDiMuons/F");
  _outTree->Branch("vertexIsValid",       &_vertexIsValid      ,"vertexIsValid/I");          
  _outTree->Branch("vertexNormChiSquare", &_vertexNormChiSquare,"vertexNormChiSquare/F");  
  _outTree->Branch("vertexChiSquare",     &_vertexChiSquare    ,"vertexChiSquare/F");      
  _outTree->Branch("vertexNDF",           &_vertexNDF          ,"vertexNDF/F");            
  _outTree->Branch("vertexX",             &_vertexX            ,"vertexX/F");              
  _outTree->Branch("vertexY",             &_vertexY            ,"vertexY/F");              
  _outTree->Branch("vertexZ",             &_vertexZ            ,"vertexZ/F");              

  _outTree->Branch("met",    &_metInfo,   "px/F:py/F:pt/F:phi/F:sumEt/F");
  _outTree->Branch("pfJets", &_pfJetInfo, "nJets/I:px[10]/F:py[10]/F:pz[10]/F:pt[10]/F:eta[10]/F:phi[10]/F:mass[10]/F:charge[10]/I:partonFlavour[10]:chf[10]/F:nhf[10]/F:cef[10]/F:nef[10]/F:muf[10]/F:hfhf[10]/F:hfef[10]/F:cm[10]/I:chm[10]/I:nhm[10]/I:cem[10]/I:nem[10]/I:mum[10]/I:hfhm[10]/I:hfem[10]/I:jecFactor[10]/F:jecUnc[10]/F:genPx[10]/F:genPy[10]/F:genPz[10]/F:genPt[10]/F:genEta[10]/F:genPhi[10]/F:genMass[10]/F:genEMF[10]/F:genHadF[10]/F:genInvF[10]/F:genAux[10]/F");

  _outTree->Branch("puJetFullDisc", 	&_puJetFullDisc 	 ,"puJetFullDisc[10]/F");              
  _outTree->Branch("puJetFullId", 	&_puJetFullId   	,"puJetFullId[10]/F");              
  _outTree->Branch("puJetSimpleDisc", 	&_puJetSimpleDisc 	 ,"puJetSimpleDisc[10]/F");              
  _outTree->Branch("puJetSimpleId", 	&_puJetSimpleId   	,"puJetSimpleId[10]/F");              
  _outTree->Branch("puJetCutDisc", 	&_puJetCutDisc 	 ,"puJetCutDisc[10]/F");              
  _outTree->Branch("puJetCutId", 	&_puJetCutId   	,"puJetCutId[10]/F");              

  // MC information
  if (_isMonteCarlo) {
    _outTree->Branch("trueMass"	, &_trueMass	, "trueMass/F");
    _outTree->Branch("truePt"	, &_truePt	, "truePt/F");
    _outTree->Branch("trueEta"  , &_trueEta     , "trueEta/F");
    _outTree->Branch("trueY"    , &_trueY       , "trueY/F");
    _outTree->Branch("truePhi"  , &_truePhi     , "truePhi/F");
    _outTree->Branch("true1"	, &_true1	, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("true2"	, &_true2	, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

    _outTree->Branch("trueMassPreFSR", &_trueMassPreFSR	, "trueMassPreFSR/F");
    _outTree->Branch("truePtPreFSR"  , &_truePtPreFSR	, "truePtPreFSR/F");
    _outTree->Branch("trueEtaPreFSR" , &_trueEtaPreFSR  , "trueEtaPreFSR/F");
    _outTree->Branch("trueYPreFSR"   , &_trueYPreFSR    , "trueYPreFSR/F");
    _outTree->Branch("truePhiPreFSR" , &_truePhiPreFSR  , "truePhiPreFSR/F");
    _outTree->Branch("true1PreFSR"   , &_true1PreFSR	, "charge/I:pt/F:ptErr/F:eta/F:phi/F");
    _outTree->Branch("true2PreFSR"   , &_true2PreFSR	, "charge/I:pt/F:ptErr/F:eta/F:phi/F");

    _outTree->Branch("genJets", &_genJetInfo, "nJets/I:px[10]/F:py[10]/F:pz[10]/F:pt[10]/F:eta[10]/F:phi[10]/F:mass[10]/F:charge[10]/I");

    _outTree->Branch("nPU", 	&_nPU   	,"nPU/I");              
  }

}


// ------------ method called once each job just after ending the event loop  ------------
void UFDiMuonsAnalyzer::endJob() {

  std::cout << "Total Number of Events Read: "<< _numEvents << std::endl <<std::endl;

  std::cout<<"number of candidate dimuons candidates: "
           <<_outTree->GetEntries()<<std::endl;
  _outFile->cd();
  _outTree->Write();

}


TLorentzVector const UFDiMuonsAnalyzer::GetLorentzVector(UFDiMuonsAnalyzer::MuonPair const* pair) {
  
  TLorentzVector muon1, muon2, sum;
  double const MASS_MUON = 0.105658367;    //GeV/c2

  reco::TrackRef const muon1Track = pair->first . innerTrack();
  reco::TrackRef const muon2Track = pair->second. innerTrack();

  muon1.SetPtEtaPhiM(muon1Track->pt(), muon1Track->eta(), muon1Track->phi(), MASS_MUON);
  muon2.SetPtEtaPhiM(muon2Track->pt(), muon2Track->eta(), muon2Track->phi(), MASS_MUON);

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


double const UFDiMuonsAnalyzer::GetMassRes(UFDiMuonsAnalyzer::MuonPair const* pair) {

  double const MASS_MUON = 0.105658367;    //GeV/c2
  
  // get the dimuon candidate
  TLorentzVector dimuon = GetLorentzVector(pair);

  // get the dimuon mass
  double dimuonMass = dimuon.M();

  // get the reference to the 2 single muons
  reco::TrackRef const track1 = pair->first .innerTrack();
  reco::TrackRef const track2 = pair->second.innerTrack();


  TLorentzVector muon1withError, muon2, dimuon_var1; 
  muon1withError.SetPtEtaPhiM(track1->pt()+track1->ptError(), 
                              track1->eta(), track1->phi(), MASS_MUON);
  muon2.SetPtEtaPhiM(track2->pt(), track2->eta(), track2->phi(), MASS_MUON);
  dimuon_var1 = muon1withError+muon2;
  double deltaM1 = fabs(dimuon_var1.M() - dimuonMass );
  
  TLorentzVector muon1, muon2withError, dimuon_var2; 
  muon1.SetPtEtaPhiM(track1->pt(), track1->eta(), track1->phi(), MASS_MUON);
  muon2withError.SetPtEtaPhiM(track2->pt()+track2->ptError(), 
                              track2->eta(), track2->phi(), MASS_MUON);
  dimuon_var2 = muon1+muon2withError;
  double deltaM2 = fabs(dimuon_var2.M() - dimuonMass );

  return std::sqrt( (deltaM1*deltaM1) + (deltaM2*deltaM2) );

}



double const UFDiMuonsAnalyzer::GetMassResCov(UFDiMuonsAnalyzer::MuonPair const* pair) {

  std::vector<reco::Track> leaves;
  
  reco::Track mu1  = *pair->first .innerTrack();
  reco::Track mu2  = *pair->second.innerTrack();

  leaves.push_back(mu1);
  leaves.push_back(mu2);

  TLorentzVector mother=GetLorentzVector(pair);
  
  double massRes = CCMassResolution.getMassResolution(leaves,mother);
  return massRes;
}

double const UFDiMuonsAnalyzer::GetMassResCov(UFDiMuonsAnalyzer::TrackPair const* pair) {

  std::vector<reco::Track> leaves;
  
  reco::Track mu1  = pair->first;
  reco::Track mu2  = pair->second;

  leaves.push_back(mu1);
  leaves.push_back(mu2);

  TLorentzVector mother=GetLorentzVector(pair);
  
  double massRes = CCMassResolution.getMassResolution(leaves,mother);
  return massRes;
}


double const UFDiMuonsAnalyzer::GetMassRes(UFDiMuonsAnalyzer::TrackPair const* pair) {

  double const MASS_MUON = 0.105658367;    //GeV/c2
  
  // get the dimuon candidate
  TLorentzVector dimuon = GetLorentzVector(pair);

  // get the dimuon mass
  double dimuonMass = dimuon.M();

  // get the reference to the 2 single muons
  reco::Track const track1 = pair->first ;
  reco::Track const track2 = pair->second;


  TLorentzVector muon1withError, muon2, dimuon_var1; 
  muon1withError.SetPtEtaPhiM(track1.pt()+track1.ptError(), 
                              track1.eta(), track1.phi(), MASS_MUON);
  muon2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), MASS_MUON);
  dimuon_var1 = muon1withError+muon2;
  double deltaM1 = fabs(dimuon_var1.M() - dimuonMass );
  
  TLorentzVector muon1, muon2withError, dimuon_var2; 
  muon1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), MASS_MUON);
  muon2withError.SetPtEtaPhiM(track2.pt()+track2.ptError(), 
                              track2.eta(), track2.phi(), MASS_MUON);
  dimuon_var2 = muon1+muon2withError;
  double deltaM2 = fabs(dimuon_var2.M() - dimuonMass );

  return sqrt( (deltaM1*deltaM1) + (deltaM2*deltaM2) );

}

// 
// get vertex fit from two muons
// 
TransientVertex const UFDiMuonsAnalyzer::GetVertexFromPair(UFDiMuonsAnalyzer::TrackPair const* muPair) const {
  //  reco::TransientTrack ttrk1 = (*transientTrackBuilder).build(*muPair.first.globalTrack());
  //  reco::TransientTrack ttrk2 = (*transientTrackBuilder).build(*muPair.second.globalTrack());
  //  reco::TransientTrack ttrk1 = (*transientTrackBuilder).build(*muPair.first.innerTrack());
  //  reco::TransientTrack ttrk2 = (*transientTrackBuilder).build(*muPair.second.innerTrack());
  reco::TransientTrack ttrk1 = (*transientTrackBuilder).build(muPair->first);
  reco::TransientTrack ttrk2 = (*transientTrackBuilder).build(muPair->second);
		
  std::vector<reco::TransientTrack> t_trks; t_trks.clear();
  t_trks.push_back(ttrk1);
  t_trks.push_back(ttrk2);


  KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
  TransientVertex tv = kvf.vertex(t_trks);
  return tv;

}

// 
// get vertex fit from two tracks
// 
TransientVertex const UFDiMuonsAnalyzer::GetVertexFromTracks(reco::TrackRef trackref1, reco::TrackRef trackref2) const{

  reco::TransientTrack ttrk1 = (*transientTrackBuilder).build(*trackref1);
  reco::TransientTrack ttrk2 = (*transientTrackBuilder).build(*trackref2);
	
  std::vector<reco::TransientTrack> t_trks; t_trks.clear();
  t_trks.push_back(ttrk1);
  t_trks.push_back(ttrk2);


  KalmanVertexFitter kvf(true); // false means no smoothing which means no track re-fit
  TransientVertex tv = kvf.vertex(t_trks);
  return tv;

}
 

// check the HLT configuration for each run. It may change you know ;-)
void UFDiMuonsAnalyzer::beginRun(edm::Run const& iRun, 
                                 edm::EventSetup const& iSetup)
{

  if (_isVerbose)
    std::cout << "\n UFDiMuonsAnalyzer::beginRun \n";
 
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    
    // if you set the trigger to look for the lowest unprescaled single mu trigger 
    // only there is no need to check its existence as it automatically will be
    // found later in the code
    if (changed) {
      
      // check if trigger name in (new) config
      const unsigned int n(hltConfig_.size());
      
      // loop over all the trigger names provided
      unsigned int triggerSize = _selectLowestSingleMuTrigger ? triggerNames_.size()-1 : triggerNames_.size();
      for (unsigned int iTrigger=0; iTrigger<triggerSize; iTrigger++) {

        addVersion(hltConfig_, triggerBaseNames_[iTrigger], triggerNames_[iTrigger]);
        if (_isVerbose)
          std::cout << "The trigger after is " << triggerNames_[iTrigger]  << std::endl;

	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerNames_[iTrigger]));
	if (triggerIndex>=n) {
          std::cout << "\n\nHLTEventAnalyzerAOD::analyze:"
                    << " TriggerName \"" << triggerNames_[iTrigger] 
                    << "\" is NOT available in (new) config!" << std::endl << std::endl;
          std::cout << " The available TriggerNames are: " << std::endl;
	  hltConfig_.dump("Triggers");
          
          throw cms::Exception("UFDiMuonsAnalyzer")<< "Throwing an exception because "
                                                   << "the trigger path name you want to check DOES NOT EXIST";
	}
        else filterNames_[iTrigger] = findFilterName ( hltConfig_, triggerNames_[iTrigger] ) ;
        
      }
      // dear god you do not want to uncomment them... but one day you could be
      // interested so I leave them there as a potential reference.
      //hltConfig_.dump("Streams");
      //hltConfig_.dump("Datasets");
      //hltConfig_.dump("PrescaleTable");
      //hltConfig_.dump("ProcessPSet");
    }
  } else {
    cout << "HLTEventAnalyzerAOD::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
    
    throw cms::Exception("UFDiMuonsAnalyzer")<<"Wrong processName_(\""<<processName_
                                             <<"\"): please double check what you passed "
                                             <<"in the python file... ";
    
  }


}


// this method will simply check is the selected HLT path (via triggerName)
// is run and accepted and no error are found
//
// bool true  if (run && accept && !error)
//      false if any other combination
bool UFDiMuonsAnalyzer::isHltPassed(const edm::Event& iEvent, 
                                    const edm::EventSetup& iSetup, 
                                    const std::string& triggerName) {
  
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (_isVerbose)
    std::cout << "\nisHltPassed::Analyzing The Trigger "<< triggerName << "...\n";
      
  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "UFDiMuonsAnalyzer::isHltPassed: path "
	 << triggerName << " - not found!" << endl;
    return false;
  }

  // Results from TriggerResults product
  if (_isVerbose) {
    std::cout << " Trigger path status:"
              << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
              << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
              << " Error =" << triggerResultsHandle_->error(triggerIndex)
              << std::endl;
  }

  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error  = triggerResultsHandle_->error (triggerIndex);

  bool isPassed=true;  

  if (!wasRun) isPassed=false;
  if (!accept) isPassed=false;
  if ( error ) isPassed=false;

  return isPassed;
}

// same check for isHltPassed +
// check if the muon is the one firing the HLT path
bool UFDiMuonsAnalyzer::isHltMatched(const edm::Event& iEvent, 
                                     const edm::EventSetup& iSetup, 
                                     const std::string& triggerName, 
                                     const std::string& filterName, 
                                     const reco::Muon& muon,
                                     float& hltPt, float& hltEta, float& hltPhi){

  double DRmin = 999;

  if (_isVerbose)  std::cout << "isHltMatched::analyzing hlt matching" << std::endl;

  bool isMatched=false;

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (_isVerbose) std::cout << "ANALYZING THE TRIGGER "<< triggerName << "..."<<std::endl;

  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "UFDiMuonsAnalyzer::isHltMatched: path "
	 << triggerName << " - not found!" << endl;
    return isMatched;
  }

  // Results from TriggerResults product
  //std::cout << " Trigger path status:"
  //          << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
  //          << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
  //          << " Error =" << triggerResultsHandle_->error(triggerIndex)
  //          << std::endl;
  
  bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
  bool accept = triggerResultsHandle_->accept(triggerIndex);
  bool error  = triggerResultsHandle_->error (triggerIndex);

  // Results from TriggerResults product
  if (_isVerbose)
    std::cout << " Trigger path status:"
              << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
              << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
              << " Error =" << triggerResultsHandle_->error(triggerIndex)
              << std::endl;

  if (!wasRun) return isMatched;
  if (!accept) return isMatched;
  if ( error ) return isMatched;
  
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  if (_isVerbose)
    std::cout << " Last active module - label/type: "
              << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
              << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
              << std::endl;
  assert (moduleIndex<m);

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
  
    // check only the last module
    if (moduleLabel != filterName) continue;

    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
  
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
      if (_isVerbose)
        std::cout << " 'L3' filter in slot " << j 
                  << " - label/type "        << moduleLabel 
                  << "/" << moduleType << std::endl;
      
      const Vids& VIDS (triggerEventHandle_->filterIds (filterIndex));
      const Keys& KEYS (triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);

      const size_type n(max(nI,nK));
      if (_isVerbose)
        std::cout << "   " << n  << " accepted 'L3' objects found: " << std::endl;

      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
    
      for (size_type i=0; i!=n; ++i) {
        const TriggerObject& TO(TOC[KEYS[i]]);
        
   
        //double DPT=fabs(TO.pt()-muon.pt());

        if (_isVerbose) {
          std::cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
                    << TO.id()  << " " << TO.pt() << " " << TO.eta() << " " 
                    << TO.phi() << " " << TO.mass()
                    << std::endl;

          std::cout << " muon: eta=" << muon.eta()
                    << ", phi=" << muon.phi()
                    << std::endl;

          std::cout << " DR=" << DR( TO.eta(),muon.eta(), TO.phi(),muon.phi() )
                    << " [max is 0.2]"
                    << std::endl;
        }
        
        if (_isTriggerEmulated && TO.pt() < _emulatedPt) continue;

        if ( DR( TO.eta(),muon.eta(), TO.phi(),muon.phi() ) < DRmin ) {
          hltPt=TO.pt();
          hltEta=TO.eta();
          hltPhi=TO.phi();
          DRmin=DR( TO.eta(),muon.eta(), TO.phi(),muon.phi() );
          if ( DRmin < 0.2 ) isMatched=true;      
        }
      }
    }
    
  }
  
  return isMatched;
}


// same check for isHltPassed +
// check if the muon is the one firing the HLT path
bool UFDiMuonsAnalyzer::isHltMatched(const edm::Event& iEvent, 
                                     const edm::EventSetup& iSetup, 
                                     const std::vector<std::string>& triggerNames, 
                                     const reco::Muon& muon1,
                                     const reco::Muon& muon2){
  bool isMatched = false;

  for (unsigned int iTrigger=0; iTrigger<triggerNames.size(); iTrigger++) {

    bool isDoubleMuTrigger=isDoubleMu(triggerNames[iTrigger]);
    
    float hltPtTmp=-999;
    float hltEtaTmp=-999;
    float hltPhiTmp=-999;

    if ( isDoubleMuTrigger ) {
      if( isHltMatched(iEvent, iSetup, triggerNames[iTrigger], filterNames_[iTrigger], muon1, hltPtTmp,hltEtaTmp,hltPhiTmp) && 
          isHltMatched(iEvent, iSetup, triggerNames[iTrigger], filterNames_[iTrigger], muon2, hltPtTmp,hltEtaTmp,hltPhiTmp)  ) {
        isMatched=true;
        break;
      }
    }

    else{
      if( isHltMatched(iEvent, iSetup, triggerNames[iTrigger], filterNames_[iTrigger], muon1, hltPtTmp,hltEtaTmp,hltPhiTmp) || 
          isHltMatched(iEvent, iSetup, triggerNames[iTrigger], filterNames_[iTrigger], muon2, hltPtTmp,hltEtaTmp,hltPhiTmp)  ) {
        isMatched=true;
        break;
      }
    }
    
  }
  
  return isMatched;
}


double 
UFDiMuonsAnalyzer::DR(double eta1, double eta2,
                      double phi1, double phi2){
  
  double diffEta = eta1 - eta2;
  double diffPhi = phi1 - phi2;
  double dr = sqrt(diffEta*diffEta + diffPhi*diffPhi);

  return dr;
  
}

bool UFDiMuonsAnalyzer::isPreselected(const reco::Muon& muon,
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
  if(!muon.isGlobalMuon() && !muon.isTrackerMuon()) return pass;
  
  // if all the cuts are passed
  pass=true;
  return pass;
  
}

bool UFDiMuonsAnalyzer::passKinCuts(const reco::Muon& muon,
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

  reco::TrackRef globalTrack;
  if (muon.isGlobalMuon()) globalTrack = muon.globalTrack();
  else if (muon.isTrackerMuon()) globalTrack = muon.innerTrack();
  else {
    // redundant: just in case...
    std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
    return false;
  }

  if (_isVerbose) {
    std::cout<< "muon.pt(): "  << muon.pt() << " [ptMin="  << _ptMin  <<"]" << std::endl;
    std::cout<< "fabs(muon.eta()): " << fabs(muon.eta())<< " [etaMax=" << _etaMax <<"]" << std::endl;
    std::cout<< "numberOfValidTrackerHits(): " << globalTrack->hitPattern().numberOfValidTrackerHits() 
             << " [min=" << _numValidTrackerHitsMin << "]" << std::endl;
    
    std::cout<<"d0: " << fabs(globalTrack->dxy(beamSpotHandle->position()))
             << " [max=" << _d0Max << "]" << std::endl;
    
  }

    
  // kinematic cuts
  if (muon.pt()        < _ptMin ) return passKinCuts; // pt cut
  if (fabs(muon.eta()) > _etaMax) return passKinCuts; // eta cut
  if (globalTrack->hitPattern().numberOfValidTrackerHits() < _numValidTrackerHitsMin) return passKinCuts; // # hits in tracker

  // beam spot cut
  if (fabs(globalTrack->dxy(beamSpotHandle->position())) > _d0Max) return passKinCuts;


  // =========================================================== //
  // What was corresponding to the old Tight VBTF
  // + some additions
  // =========================================================== //
  if (_isVerbose) {
    std::cout << "numberOfValidMuonHits: " 
              << globalTrack->hitPattern().numberOfValidMuonHits() 
              << " [min=" << _numValidMuonHitsMin << "]" 
              << std::endl;
  
    std::cout << "numberOfValidPixelHits(): " 
              << globalTrack->hitPattern().numberOfValidPixelHits() 
              << " [min=" << _numValidPixelHitsMin << "]" 
              << std::endl;
     
    std::cout << "numberOfValidStripHits(): " 
              << globalTrack->hitPattern().numberOfValidStripHits() 
              << " [min=" << _numValidStripHitsMin << "]" 
              << std::endl;
    
    std::cout << "pixelLayersWithMeasurement(): " 
              << globalTrack->hitPattern().pixelLayersWithMeasurement() 
              << " [min=" << _numPixelLayersMin << "]" 
              << std::endl;

    std::cout << "trackerLayersWithMeasurement(): " 
              << globalTrack->hitPattern().trackerLayersWithMeasurement() 
              << " [min=" << _numTrackerLayersMin << "]" 
              << std::endl;

    std::cout << "stripLayersWithMeasurement(): " 
              << globalTrack->hitPattern().stripLayersWithMeasurement() 
              << " [min=" << _numStripLayersMin << "]" 
              << std::endl;

    std::cout << "validFraction(): " 
              << globalTrack-> validFraction()
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
  
    std::cout<<"globalTrack->normalizedChi2(): " 
             << globalTrack->normalizedChi2()
             << " [max=" << _normChiSquareMax <<"]" 
             << std::endl;
   
  }

  if (globalTrack->hitPattern().pixelLayersWithMeasurement()   < _numPixelLayersMin)   return passKinCuts;   
  if (globalTrack->hitPattern().trackerLayersWithMeasurement() < _numTrackerLayersMin) return passKinCuts; 
  if (globalTrack->hitPattern().stripLayersWithMeasurement()   < _numTrackerLayersMin) return passKinCuts;   
  		 
  if (globalTrack->validFraction() < _validFracTrackerMin) return passKinCuts; 

  if ( globalTrack->hitPattern().numberOfValidMuonHits() < _numValidMuonHitsMin  )  return passKinCuts;
  if ( globalTrack->hitPattern().numberOfValidPixelHits() < _numValidPixelHitsMin ) return passKinCuts;
  if ( globalTrack->hitPattern().numberOfValidStripHits() < _numValidStripHitsMin ) return passKinCuts;
  if ( muon.numberOfMatches() < _numSegmentMatchesMin   )         return passKinCuts;
  if ( muon.numberOfMatchedStations() < _numOfMatchedStationsMin) return passKinCuts;
  if ( globalTrack->normalizedChi2() > _normChiSquareMax)         return passKinCuts;

  if (_isVerbose) std::cout << "passing kinematic cuts\n"; 

  passKinCuts=true;
  return passKinCuts;
}

UFDiMuonsAnalyzer::MuonPairs const UFDiMuonsAnalyzer::GetMuonPairs(reco::MuonCollection const* muons) const {
                                                                           
  
  MuonPairs muonpairs; 
  muonpairs.clear();

  for (reco::MuonCollection::const_iterator muon1 = muons->begin(), 
         muonsEnd = muons->end(); muon1 != muonsEnd; ++muon1){

    for (reco::MuonCollection::const_iterator muon2 = muon1+1; 
         muon2 != muonsEnd; ++muon2){
      muonpairs.push_back( UFDiMuonsAnalyzer::MuonPair(*muon1,*muon2) );
    }
  }
  
  std::sort(muonpairs.begin(),muonpairs.end(),sortMuonObject);

  return muonpairs;
}

void UFDiMuonsAnalyzer::initTrack(_TrackInfo& track) {

  track.charge = -999; 
  track.pt     = -999; 
  track.ptErr  = -999;
  track.eta    = -999; 
  track.phi    = -999;

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

  for (unsigned int iTrigger=0;iTrigger<3;iTrigger++) {
    muon.isHltMatched[iTrigger] = -999;
    muon.hltPt[iTrigger] = -999;
    muon.hltEta[iTrigger] = -999;
    muon.hltPhi[iTrigger] = -999;
  }

}


void UFDiMuonsAnalyzer::displaySelection() {

  std::cout << "\n\n*** UFDiMuonsAnalyzer Configuration ***\n";
  std::cout << " - Events saved in file: " << _getFilename << std::endl;	

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
  std::cout << " - _selectLowestSingleMuTrigger: " << _selectLowestSingleMuTrigger << std::endl;
  std::cout << " - Additional Triggers To Probe:\n";
  unsigned int triggerSize = _selectLowestSingleMuTrigger ? triggerNames_.size()-1 : triggerNames_.size();
  for (unsigned int i=0; i < triggerSize; i++) 
    std::cout << "    * triggerBaseNames["<<i<<"]: " << triggerBaseNames_[i] << std::endl;
  
  std::cout << std::endl << std::endl;

}

void UFDiMuonsAnalyzer::findLowestSingleMu(const HLTConfigProvider hltConfig_,
                                           const edm::Event& iEvent,    
                                           const edm::EventSetup& iSetup) {

  // originally from Michele's code at DQM/Physics/src/EwkMuLumiMonitorDQM.cc
  
  // Check the prescales
  // see the trigger single muon which are present
  std::string lowestMuonUnprescaledTrig = "";
  bool lowestMuonUnprescaledTrigFound = false;
  const std::vector<std::string>& triggerNames = hltConfig_.triggerNames();
  
  for (size_t ts = 0; ts< triggerNames.size() ; ts++){
    std::string trig = triggerNames[ts];
    size_t f = trig.find("HLT_Mu");

    size_t fmr = trig.find("MR");

    size_t weirdName = trig.find("L1Mu10erJetC12WdEtaPhi1DiJetsC");

    if ( (f != std::string::npos) && 
         (fmr == std::string::npos) && 
         (weirdName == std::string::npos) ) {

      //std::cout << "single muon trigger present: " << trig << std::endl; 

      int prescaleSet = hltConfig_.prescaleSet(iEvent, iSetup);
      //std::cout << "prescaleSet=" << prescaleSet << std::endl;

      if (prescaleSet == -1) continue;
      
      std::pair<int,int> prescaleValues = hltConfig_.prescaleValues(iEvent, iSetup,trig);
      int l1prescale  = prescaleValues.first;
      int hltprescale = prescaleValues.second;

      //std::cout << "L1 prescale="<< l1prescale //prescaleValues.first
      //         << " HLT prescale=" << hltprescale << std::endl;//prescaleValues.second << std::endl;

      if (l1prescale != 1 || hltprescale != 1) continue;

      for (unsigned int n=9; n<100 ; n++ ){
        std::string lowestTrig= "finalVersion";
        std::string lowestTrigv0 = "HLT_Mu";
        std::stringstream out;
        out << n;
        std::string s = out.str(); 
        lowestTrigv0.append(s);
        lowestTrig = lowestTrigv0;  
        
        //std::cout << " lowestTrig: " << lowestTrig << std::endl; 
        if (trig==lowestTrig) lowestMuonUnprescaledTrig = trig ;
        if (trig==lowestTrig) {std::cout << " before loop, lowestTrig lowest single muon trigger present unprescaled: " << lowestTrig << std::endl; }
        if (trig==lowestTrig) lowestMuonUnprescaledTrigFound = true ;
        if (trig==lowestTrig) break ;
          
        for (unsigned int v = 1; v<10 ; v++ ){
          lowestTrig = lowestTrigv0;
          
          size_t eta2p1 = trig.find("eta2p1");
          if (eta2p1 != std::string::npos) lowestTrig.append("_eta2p1");

          lowestTrig.append("_v");
          std::stringstream oout;
          oout << v;
          std::string ss = oout.str(); 
          lowestTrig.append(ss);
          if (trig==lowestTrig) lowestMuonUnprescaledTrig = trig ;
          if (trig==lowestTrig) lowestMuonUnprescaledTrigFound = true ;
          //std::cout << "lowestTrig=" << lowestTrig << "found?" << lowestMuonUnprescaledTrigFound << std::endl;
          if (trig==lowestTrig) break ;
        }	
        if (lowestMuonUnprescaledTrigFound) break; 
	
      }
      if (lowestMuonUnprescaledTrigFound) break; 
    }
  }

  if (_isVerbose) {
    std::cout << "after break, lowest single muon trigger present unprescaled: " << lowestMuonUnprescaledTrig << std::endl; 
    if (lowestMuonUnprescaledTrig == "") 
      edm::LogError("UFDiMuonsAnalyzer") << "Lowest Unprescaled Single Muon Trigger NOT FOUND!\n";
    std::cout << "Filter Name is " << findFilterName ( hltConfig_, lowestMuonUnprescaledTrig ) << std::endl;
  }
  
  // *********************************** //
  triggerNames_.back() = lowestMuonUnprescaledTrig ;
  filterNames_.back() = findFilterName ( hltConfig_, lowestMuonUnprescaledTrig );
  // *********************************** //

}


void UFDiMuonsAnalyzer::addVersion(const HLTConfigProvider hltConfig_,
                                   std::string& triggerBaseName,
                                   std::string& triggerName) {

  //std::cout << "The trigger is " << triggerBaseName  << std::endl;
  
  // This function looks gets a trigger name, e.g. HLT_Mu15 and add to it
  // its version number, e.g. HLT_Mu15_v2, if the trigger is present in 
  // the hltConfig  
  const std::vector<std::string>& triggerNames = hltConfig_.triggerNames();
    
  for (size_t ts = 0; ts< triggerNames.size() ; ts++){
    std::string trig = triggerNames[ts];
    size_t f = trig.find(triggerBaseName);

    if (f != std::string::npos)  {

      //adding version extension
      for (unsigned int iVersion = 1; iVersion<60; iVersion++ ){
        std::string trigWithVersion= triggerBaseName;
        trigWithVersion.append("_v");
        std::stringstream ss;
        ss << iVersion;
        std::string version = ss.str(); 
        trigWithVersion.append(version);
        
        if (trig==trigWithVersion) {
          triggerName.replace(triggerName.begin(),    triggerName.end(),
                              trigWithVersion.begin(),trigWithVersion.end());
          return;
        }
        
      }
    }
    
  }
  
  std::cout << "TriggerBaseName \"" << triggerBaseName 
            << "\" was already \"versionized\" OR "
            << " the root is not found in the list of trigger which"
            << " means cmssw is going to crash providing the list"
            << " of available triggers\n";

  triggerName.replace(triggerName.begin(),    triggerName.end(),
                      triggerBaseName.begin(),triggerBaseName.end());

  return;

}


std::string UFDiMuonsAnalyzer::findFilterName(const HLTConfigProvider hltConfig_, 
                                              const std::string& triggerName){
  
  std::string L3FilterName_;
  
  const std::vector<std::string>& moduleLabs = hltConfig_.moduleLabels(triggerName); 
    
  // the l3 filter name is just the last module.... 
  size_t moduleLabsSizeMinus2 = moduleLabs.size() - 2 ;
  
 
  L3FilterName_ = moduleLabs[moduleLabsSizeMinus2];

  if (_isVerbose) std::cout<<"triggerName="    << triggerName
                           <<" -> filterName=" << L3FilterName_ << std::endl;        
  
  // *********************************** //
  return L3FilterName_ ;
  //return moduleLabs[moduleLabsSizeMinus2];
  // *********************************** //
}


bool UFDiMuonsAnalyzer::isDoubleMu(const std::string& triggerName) {

  bool isDoubleMuTrigger=false;
  
  size_t found;
  found=triggerName.find("DoubleMu");

  if (found!=std::string::npos) isDoubleMuTrigger=true;

  return isDoubleMuTrigger;
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

DEFINE_FWK_MODULE(UFDiMuonsAnalyzer);


