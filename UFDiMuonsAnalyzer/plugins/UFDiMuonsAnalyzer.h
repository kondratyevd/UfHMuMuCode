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
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/NTupleBranches.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/NTupleHelpers.h"

// Special calibration classes
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"
#include "RochCor/Calibration/interface/rochcor2016.h"
#include "RochCor/Calibration/interface/RoccoR.h"
#include <time.h>

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

  // meta-data not given in python config file
  // info gathered from py cfg defined later (py-cfg meta data: isMonteCarlo, trigNames, tauIDNames, bTagName)
  int _numEvents;
  int _sumEventWeights;

  // tracks pairs, e.g. cocktail
  typedef std::pair<reco::Track,reco::Track> TrackPair;
  typedef std::vector<TrackPair> TrackPairs;

  // GEN info
  int    _nPU;     // True number of vertices
  float  _PU_wgt;  // Pileup weight
  float  _PU_wgt_up;
  float  _PU_wgt_down;
  TH1D*  _PU_wgt_hist;
  TH1D*  _PU_wgt_hist_up;
  TH1D*  _PU_wgt_hist_down;
  TFile* _PU_wgt_file;
  int _GEN_wgt;    // +1 or -1 weight for nlo samples, -1 simulates interference when filling histos


  ///////////////////////////////////////////////////////////
  // Structs  ==============================================
  //////////////////////////////////////////////////////////
  
  // general event information	
  EventInfo _eventInfo;

  // vector of vertex information
  VertexInfos _vertexInfos;
  int _nVertices;

  // vector of muons
  MuonInfos _muonInfos;
  int _nMuons;

  // info about the dimuon candidates in _muonInfos
  PairInfos _pairInfos;
  int _nPairs;

  // vector of electrons
  EleInfos  _eleInfos;
  int _nEles;

  // vector of taus
  TauInfos  _tauInfos;
  int _nTaus;

  // generator level info Gamma pre-FSR
  GenPartInfo _genGpreFSR, _genM1GpreFSR, _genM2GpreFSR;

  // generator level info Gamma post-FSR
  GenPartInfo _genGpostFSR, _genM1GpostFSR, _genM2GpostFSR;

  // generator level info Z pre-FSR
  GenPartInfo _genZpreFSR, _genM1ZpreFSR, _genM2ZpreFSR;

  // generator level info Z post-FSR
  GenPartInfo _genZpostFSR, _genM1ZpostFSR, _genM2ZpostFSR;

  // generator level info W pre-FSR
  GenPartInfo _genWpreFSR, _genMWpreFSR;

  // generator level info W post-FSR
  GenPartInfo _genWpostFSR, _genMWpostFSR;

  // generator level info H pre-FSR
  GenPartInfo _genHpreFSR, _genM1HpreFSR,_genM2HpreFSR;

  // generator level info H post-FSR
  GenPartInfo _genHpostFSR, _genM1HpostFSR,_genM2HpostFSR;

  // Jets and MET
  JetInfos     _jetInfos;
  SlimJetInfos _slimJetInfos;
  int _nJets, _nJetsCent, _nJetsFwd;
  int _nBLoose, _nBMed, _nBTight;
  JetInfos     _jetInfos_JES_up;
  SlimJetInfos _slimJetInfos_JES_up;
  int _nJets_JES_up, _nJetsCent_JES_up, _nJetsFwd_JES_up;
  int _nBLoose_JES_up, _nBMed_JES_up, _nBTight_JES_up;
  JetInfos     _jetInfos_JES_down;
  SlimJetInfos _slimJetInfos_JES_down;
  int _nJets_JES_down, _nJetsCent_JES_down, _nJetsFwd_JES_down;
  int _nBLoose_JES_down, _nBMed_JES_down, _nBTight_JES_down;
  JetInfos     _jetInfos_JER_up;
  SlimJetInfos _slimJetInfos_JER_up;
  int _nJets_JER_up, _nJetsCent_JER_up, _nJetsFwd_JER_up;
  int _nBLoose_JER_up, _nBMed_JER_up, _nBTight_JER_up;
  JetInfos     _jetInfos_JER_down;
  SlimJetInfos _slimJetInfos_JER_down;
  int _nJets_JER_down, _nJetsCent_JER_down, _nJetsFwd_JER_down;
  int _nBLoose_JER_down, _nBMed_JER_down, _nBTight_JER_down;

  GenJetInfos _genJetInfos;
  int _nGenJets;

  MetInfo     _metInfo;
  MhtInfo     _mhtInfo;
  MhtInfo     _mhtInfo_JES_up;
  MhtInfo     _mhtInfo_JES_down;
  MhtInfo     _mhtInfo_JER_up;
  MhtInfo     _mhtInfo_JER_down;

  // Weights and efficiencies
  float  _IsoMu_eff_3;
  float  _IsoMu_eff_3_up;
  float  _IsoMu_eff_3_down;
  TH2F*  _IsoMu_eff_3_hist;
  TFile* _IsoMu_eff_3_file;
  float  _IsoMu_eff_4;
  float  _IsoMu_eff_4_up;
  float  _IsoMu_eff_4_down;
  TH2F*  _IsoMu_eff_4_hist;
  TFile* _IsoMu_eff_4_file;
  float  _IsoMu_eff_bug;
  float  _IsoMu_eff_bug_up;
  float  _IsoMu_eff_bug_down;

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
		    const edm::Handle<edm::TriggerResults>& trigResultsHandle,
		    const std::vector<std::string> trigNames);

  void displaySelection();

  static bool sortMuonsByPt (pat::Muon i,     pat::Muon j    );
  static bool sortElesByPt  (pat::Electron i, pat::Electron j);
  static bool sortTausByPt  (pat::Tau i,      pat::Tau j     );
  static bool sortJetsByPt  (pat::Jet i,      pat::Jet j     );

  KalmanMuonCalibrator _KaMu_calib;
  bool _doSys_KaMu;
  rochcor2016* _Roch_calib[201];
  bool _doSys_Roch;

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
  std::vector<std::string> _trigNames;

  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _trigObjsToken;
  
  // Muons
  edm::EDGetTokenT<pat::MuonCollection> _muonCollToken;

  // Electrons
  edm::EDGetTokenT<edm::View<pat::Electron>> _eleCollToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _eleIdVetoToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _eleIdLooseToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _eleIdMediumToken;
  edm::EDGetTokenT< edm::ValueMap<bool> >    _eleIdTightToken;

  // Taus
  edm::EDGetTokenT<pat::TauCollection> _tauCollToken;
  std::vector<std::string> _tauIDNames;

  // Jets / MET
  edm::EDGetTokenT<pat::METCollection> _metToken;
  edm::EDGetTokenT<pat::JetCollection> _jetsToken;
  std::string _btagName;

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
  bool _doSys;
  bool _slimOut;

  int  _skim_nMuons;
  bool _skim_trig;

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

  std::string _ele_ID;
  double _ele_pT_min;
  double _ele_eta_max;

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

