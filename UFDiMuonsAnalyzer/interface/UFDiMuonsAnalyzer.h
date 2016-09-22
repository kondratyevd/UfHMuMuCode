///////////////////////////////////////////////////////////
//=========================================================
// UFDiMuonsAnalyzer.h                                   //
//=========================================================
//                                                       //
// v1 coded in 2010 by Gian Pierro. v2 revamps that code.//
// v2 coded by Andrew Carnes in 2015.                    //
//                                                       //
//                                                       //
//========================================================
///////////////////////////////////////////////////////////

#ifndef ADD_UDAV2
#define ADD_UDAV2

///////////////////////////////////////////////////////////
// Includes ==============================================
//////////////////////////////////////////////////////////
// Tons and tons of includes or it wouldn't be CMSSW.

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

// CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// FWCore General
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"


// DataFormats General
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DetId/interface/DetId.h"

// Geometry
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// math classes
#include "DataFormats/Math/interface/deltaR.h"

// trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// vertexing
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// gen particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// Adding the Muon Cocktail
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

// PAT objects
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"


// Add the data formats that we store most of the info into
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

///////////////////////////////////////////////////////////
// Class Definition ======================================
//////////////////////////////////////////////////////////

class UFDiMuonsAnalyzer : public edm::EDAnalyzer {

public:

  ///////////////////////////////////////////////////////////
  // Constructors/Destructors===============================
  //////////////////////////////////////////////////////////

  explicit UFDiMuonsAnalyzer(const edm::ParameterSet&);
  ~UFDiMuonsAnalyzer();

  ///////////////////////////////////////////////////////////
  // Basic types ===========================================
  ///////////////////////////////////////////////////////////
  
  //Physical Constants
  double static constexpr PDG_MASS_Z  = 91.1876;      //GeV/c2
  double static constexpr PDG_WIDTH_Z = 2.4952;       //GeV/c2
  double static constexpr MASS_MUON = 0.105658367;    //GeV/c2

  int _numEvents;

  // tracks pairs, e.g. cocktail
  typedef std::pair<reco::Track,reco::Track> TrackPair;
  typedef std::vector<TrackPair> TrackPairs;

  // rho
  float _rho;
  float _rho25;
  float _rho25asHtoZZto4l;

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

  int _nPU;
  int _genWeight;
  int _sumEventWeights;


  ///////////////////////////////////////////////////////////
  // Structs  ==============================================
  //////////////////////////////////////////////////////////
  
  // general event information	
  _EventInfo eventInfo;

  // vertex information
  _VertexInfo vertexInfo;

   // muon info
  _MuonInfo _muon1, _muon2;

  // dimuon candidates from refit with Vertex Constraint
  _TrackInfo _muon1vc,_muon2vc; //dimuon
  _TrackInfo _muon1pvc,_muon2pvc; //pv

  // generator level info Gamma pre-FSR
  _genPartInfo _genGpreFSR;
  _TrackInfo   _genM1GpreFSR,_genM2GpreFSR;

  // generator level info Gamma post-FSR
  _genPartInfo _genGpostFSR;
  _TrackInfo   _genM1GpostFSR,_genM2GpostFSR;

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

  // muons pairs
  typedef std::pair<pat::Muon,pat::Muon> MuonPair;
  typedef std::vector<MuonPair>          MuonPairs;

  // sort the dimuons putting in front the candidates whose mass is 
  // the closest to the Z PDG value
  struct sortMuonsClass {
    bool operator() (MuonPair pair1, 
  		     MuonPair pair2) {

      TLorentzVector muon11, muon12, dimuon1;

      reco::Track const muon11Track = *(pair1.first . innerTrack());
      reco::Track const muon12Track = *(pair1.second. innerTrack());

      muon11.SetPtEtaPhiM(muon11Track.pt(), muon11Track.eta(), muon11Track.phi(), MASS_MUON);
      muon12.SetPtEtaPhiM(muon12Track.pt(), muon12Track.eta(), muon12Track.phi(), MASS_MUON);

      dimuon1 = muon11+muon12;

      TLorentzVector muon21, muon22, dimuon2;

      reco::Track const muon21Track = *(pair2.first . innerTrack());
      reco::Track const muon22Track = *(pair2.second. innerTrack());

      muon21.SetPtEtaPhiM(muon21Track.pt(), muon21Track.eta(), muon21Track.phi(), MASS_MUON);
      muon22.SetPtEtaPhiM(muon22Track.pt(), muon22Track.eta(), muon22Track.phi(), MASS_MUON);

      dimuon2 = muon21+muon22;

      return ( fabs(dimuon1.M()-PDG_MASS_Z) < fabs(dimuon2.M()-PDG_MASS_Z) );
    }
  } sortMuonObject;

  ///////////////////////////////////////////////////////////
  // Trees  ================================================
  //////////////////////////////////////////////////////////

  // where to save all the info  
  TTree* _outTree;
  TTree* _outTreeMetadata;


private:

  ///////////////////////////////////////////////////////////
  // Inheritance from EDAnalyzer ===========================
  //////////////////////////////////////////////////////////

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  ///////////////////////////////////////////////////////////
  // Useful Functions ======================================
  //////////////////////////////////////////////////////////

  edm::Service<TFileService> fs;

  void initMuon(_MuonInfo& muon);
  void initTrack(_TrackInfo& track);
  void initGenPart(_genPartInfo& part);

  void fillBosonAndMuDaughters(const reco::Candidate* boson);

  // methods to select the muons
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

  MuonPairs  const GetMuonPairs(pat::MuonCollection  const* muons ) const;
  TrackPairs const GetTrackPairs(reco::TrackCollection const* tracks) const;

  TLorentzVector const GetLorentzVector(MuonPair  const* pair) ;//const; 
  TLorentzVector const GetLorentzVector(TrackPair const* pair) ;//const; 

  static bool sortGenJetFunc(reco::GenJet i, reco::GenJet j);

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Gather info from the python config file ==============
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////
  // Handles/Tokens  =======================================
  //////////////////////////////////////////////////////////
  // you get access to the information from the collections designated by the python config file
  // using the handle-token infrastructure. Get the name of the collection from the config
  // and load it into the token, then pass the token and the handle into the Event.

  // Not sure why these are defined here rather than in the analyze function like all the other handles
  // Should get rid of these handles and put them in the correct place
  edm::Handle<edm::TriggerResults> _triggerResultsHandle;
  edm::Handle<pat::TriggerObjectStandAloneCollection> _triggerObjsHandle;

  // trigger
  edm::EDGetTokenT<edm::TriggerResults> _triggerResultsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _triggerObjsToken;
  
  // muons
  edm::EDGetTokenT<pat::MuonCollection> _muonCollToken;
  edm::EDGetTokenT<reco::BeamSpot> _beamSpotToken;		
  edm::EDGetTokenT<reco::GenParticleCollection> _prunedGenParticleToken;		
  edm::EDGetTokenT<pat::PackedGenParticleCollection> _packedGenParticleToken;		
  edm::EDGetTokenT<reco::VertexCollection> _primaryVertexToken;		
  edm::EDGetTokenT<GenEventInfoProduct> _genEvtInfoToken;		
  edm::EDGetTokenT< std::vector< PileupSummaryInfo > > _PupInfoToken;		

  // jets
  edm::EDGetTokenT<std::vector<pat::MET>> _metToken;
  edm::EDGetTokenT<std::vector<pat::Jet>> _pfJetsToken;
  edm::EDGetTokenT<reco::GenJetCollection> _genJetsToken;

  ///////////////////////////////////////////////////////////
  // Basic types  ===========================================
  //////////////////////////////////////////////////////////

  // Selection Criteria given in config file
  double _isGlobal;
  double _isTracker;
  double _ptMin;
  double _etaMax;
  unsigned int _nMuons;
  bool _checkTrigger; 

  // useful switches given in config file
  bool _isVerbose;         // debug mode
  bool _isMonteCarlo;      // save MC truth

  // More info loaded from the config file
  std::string   _processName;
  std::vector < std::string > _triggerNames;

  // Not currently used, since we don't use prescaled triggers with this analysis
  // Can probably delete these vars from the class
  std::vector < int > _l1Prescale;
  std::vector < int > _hltPrescale;
};

DEFINE_FWK_MODULE(UFDiMuonsAnalyzer);
#endif


