///////////////////////////////////////////////////////////
//=========================================================
// UFDiMuonsAnalyzer.h                                   //
//=========================================================
//                                                       //
// H->MuMu Analyzer for UF. Handed down and revamped     //
// too many times.                                       //
//                                                       //
//========================================================
///////////////////////////////////////////////////////////

#ifndef ADD_UDAV2
#define ADD_UDAV2

// All the possible includes in CMSSW
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"

// Add the data formats that we store most of the info into
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

// Code to fill GEN branches
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/FillGenHelper.h"

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

  // Needed here? - AWB 09.11.16
  /* //Physical Constants */
  /* double static constexpr PDG_MASS_Z  = 91.1876;      //GeV/c2 */
  /* double static constexpr PDG_WIDTH_Z = 2.4952;       //GeV/c2 */
  /* double static constexpr MASS_MUON = 0.105658367;    //GeV/c2 */

  // meta-data not given in python config file
  // info gathered from py cfg defined later (py-cfg meta data: isMonteCarlo, triggerNames, tauIDNames, bTagNames)
  int _numEvents;
  int _sumEventWeights;

  // tracks pairs, e.g. cocktail
  typedef std::pair<reco::Track,reco::Track> TrackPair;
  typedef std::vector<TrackPair> TrackPairs;

  // gen info
  int _nPU;        // true pileup
  int _genWeight;  // +1 or -1 weight for nlo samples, -1 simulates interference when filling histos


  ///////////////////////////////////////////////////////////
  // Structs  ==============================================
  //////////////////////////////////////////////////////////
  
  // general event information	
  _EventInfo _eventInfo;

  // array of vertex information
  _VertexInfo _vertexInfo;

  // array of muons, 0,1 locations are the dimuon candidate
  _MuonInfo _muonInfo;

  // info about the dimu candidate formed from the 0,1 muons in _muonInfo
  _DiMuonInfo _diMuonInfo;

  // array of electrons
  _ElectronInfo  _electronInfo;

  // array of taus
  _TauInfo  _tauInfo;

  // generator level info Gamma pre-FSR
  _GenPartInfo _genGpreFSR, _genM1GpreFSR, _genM2GpreFSR;

  // generator level info Gamma post-FSR
  _GenPartInfo _genGpostFSR, _genM1GpostFSR, _genM2GpostFSR;

  // generator level info Z pre-FSR
  _GenPartInfo _genZpreFSR, _genM1ZpreFSR, _genM2ZpreFSR;

  // generator level info Z post-FSR
  _GenPartInfo _genZpostFSR, _genM1ZpostFSR, _genM2ZpostFSR;

  // generator level info W pre-FSR
  _GenPartInfo _genWpreFSR, _genMWpreFSR;

  // generator level info W post-FSR
  _GenPartInfo _genWpostFSR, _genMWpostFSR;

  // generator level info H pre-FSR
  _GenPartInfo _genHpreFSR, _genM1HpreFSR,_genM2HpreFSR;

  // generator level info H post-FSR
  _GenPartInfo _genHpostFSR, _genM1HpostFSR,_genM2HpostFSR;

  // Jets and MET
  _JetInfo    _jetInfo;
  _JetInfo    _jetInfo_JES_up;
  _JetInfo    _jetInfo_JES_down;
  _GenJetInfo _genJetInfo;
  _MetInfo    _metInfo;

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

  // methods for selection
  bool isHltPassed (const edm::Event&, const edm::EventSetup&, 
		    const edm::Handle<edm::TriggerResults>& triggerResultsHandle,
		    const std::vector<std::string> triggerNames);

  void displaySelection();

  static bool sortMuonsByPt    (pat::Muon i,     pat::Muon j    );
  static bool sortElectronsByPt(pat::Electron i, pat::Electron j);
  static bool sortTausByPt     (pat::Tau i,      pat::Tau j     );
  static bool sortJetsByPt     (pat::Jet i,      pat::Jet j     );


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

  // Trigger
  std::string _processName;
  std::vector<std::string> _triggerNames;

  edm::EDGetTokenT<edm::TriggerResults> _triggerResultsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _triggerObjsToken;
  
  // Muons
  edm::EDGetTokenT<pat::MuonCollection> _muonCollToken;

  // Electrons
  edm::EDGetTokenT<edm::View<pat::Electron>> _electronCollToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _electronIdVetoToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _electronIdLooseToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _electronIdMediumToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _electronIdTightToken;

  // Taus
  edm::EDGetTokenT<pat::TauCollection> _tauCollToken;
  std::vector<std::string> _tauIDNames;

  // Jets / MET
  edm::EDGetTokenT<pat::METCollection> _metToken;
  edm::EDGetTokenT<pat::JetCollection> _jetsToken;
  std::vector<std::string> _btagNames;

  // Event info
  edm::EDGetTokenT<reco::BeamSpot> _beamSpotToken;		
  edm::EDGetTokenT<reco::VertexCollection> _primaryVertexToken;		
  edm::EDGetTokenT< std::vector< PileupSummaryInfo > > _PupInfoToken;		

  // Gen Info
  edm::EDGetTokenT<reco::GenJetCollection> _genJetsToken;
  edm::EDGetTokenT<reco::GenParticleCollection> _prunedGenParticleToken;		
  edm::EDGetTokenT<pat::PackedGenParticleCollection> _packedGenParticleToken;		
  edm::EDGetTokenT<GenEventInfoProduct> _genEvtInfoToken;		

  ///////////////////////////////////////////////////////////
  // Basic types  ===========================================
  //////////////////////////////////////////////////////////

  // Selection Criteria for event and objects in config file
  bool _isVerbose;   
  bool _isMonteCarlo;
  bool _doSyst;

  int  _skim_nMuons;
  bool _skim_trigger;

  double _vertex_ndof_min;
  double _vertex_rho_max;
  double _vertex_z_max;

  std::string _muon_ID;
  double _muon_pT_min;
  double _muon_eta_max;
  double _muon_trig_dR;
  bool   _muon_use_pfIso;
  double _muon_iso_dR;
  double _muon_iso_max;

  std::string _electron_ID;
  double _electron_pT_min;
  double _electron_eta_max;

  double _tau_pT_min;
  double _tau_eta_max;

  std::string _jet_ID;
  double _jet_pT_min;
  double _jet_eta_max;

  // Not currently used, since we don't use prescaled triggers with this analysis
  // Can probably delete these vars from the class
  std::vector < int > _l1Prescale;
  std::vector < int > _hltPrescale;
};

DEFINE_FWK_MODULE(UFDiMuonsAnalyzer);
#endif

