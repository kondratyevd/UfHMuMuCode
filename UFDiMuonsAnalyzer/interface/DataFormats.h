
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"

#ifndef DATA_FORMATS
#define DATA_FORMATS

////////////////////////////////////////////////
////   WARNING!!! CONSTRUCTION OF STRUCTS   ////
////////////////////////////////////////////////
// All variables in the struct must have the same length
// e.g. Int_t and Float_t have 4 bytes, Long64_t and Double_t have 8 bytes
// https://twiki.cern.ch/twiki/bin/view/Main/RootNotes#Conventions_and_Types
// https://root.cern.ch/root/html534/guides/users-guide/Trees.html#adding-a-branch-to-hold-a-list-of-variables

// Event info
struct _EventInfo {
  
  // Data fields
  Long64_t run;
  Long64_t lumi;
  Long64_t event;
  Long64_t bx;
  Long64_t orbit;
  
  // Tell the ttree how to save the struct
  static TString getVarString() {
    return TString("run/L:lumi/L:event/L:bx/L:orbit/L");
  };
  
  // Fill the data fields
  void init();
  void fill(const edm::Event& iEvent);

};

// Vertex info
struct _VertexInfo {

  const static unsigned int arraySize = 1; // Make variable - AWB 08.11.16
  Long64_t nVertices;

  Long64_t isValid [arraySize];
  Double_t x       [arraySize];	
  Double_t y       [arraySize];	
  Double_t z       [arraySize];	
  Double_t xErr    [arraySize];	
  Double_t yErr    [arraySize];	
  Double_t zErr    [arraySize];	
  Double_t chi2    [arraySize];
  Long64_t ndof    [arraySize];
  Double_t normChi2[arraySize];

  static TString getVarString() {
    TString r = 
      TString("nVertices/L:isValid[N]/L:") +
      TString("x[N]/D:y[N]/D:z[N]/D:xErr[N]/D:yErr[N]/D:zErr[N]/D:") +
      TString("chi2[N]/D:ndof[N]/L:normChi2[N]/D");
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  };
  
  void init();
  void fill( const reco::VertexCollection verticesSelected );
  reco::VertexCollection select( const edm::Handle<reco::VertexCollection>& vertices, const double _vertex_ndof_min,
				 const double _vertex_rho_max, const double _vertex_z_max );
  
};


// Muon info
struct _MuonInfo {

  const static unsigned int arraySize        = 10;
  const static unsigned int triggerArraySize =  6;
  Long64_t nMuons;

  Long64_t isTracker   [arraySize];
  Long64_t isStandAlone[arraySize];
  Long64_t isGlobal    [arraySize];

  Long64_t isTightID [arraySize];
  Long64_t isMediumID[arraySize];
  Long64_t isLooseID [arraySize];

  Long64_t charge[arraySize];
  Double_t pt    [arraySize];
  Double_t ptErr [arraySize];
  Double_t eta   [arraySize];
  Double_t phi   [arraySize];

  Double_t trkPt   [arraySize];
  Double_t trkPtErr[arraySize];
  Double_t trketa  [arraySize];
  Double_t trkPhi  [arraySize];

  Double_t d0_BS[arraySize];
  Double_t dz_BS[arraySize];

  Double_t d0_PV[arraySize];
  Double_t dz_PV[arraySize];

  Double_t relIso           [arraySize];
  Double_t relCombIso       [arraySize];
  Double_t trackIsoSumPt    [arraySize];
  Double_t trackIsoSumPtCorr[arraySize];
  Double_t hcalIso          [arraySize];
  Double_t ecalIso          [arraySize];

  // PF information
  Long64_t isPF[arraySize];

  Double_t pfPt [arraySize];
  Double_t pfEta[arraySize];
  Double_t pfPhi[arraySize];
  
  Double_t sumChargedHadronPtR03  [arraySize];  // sum-pt of charged Hadron
  Double_t sumChargedParticlePtR03[arraySize];  // sum-pt of charged Particles(inludes e/mu)
  Double_t sumNeutralHadronEtR03  [arraySize];  // sum pt of neutral hadrons
  Double_t sumPhotonEtR03         [arraySize];  // sum pt of PF photons
  Double_t sumPUPtR03             [arraySize];  // sum pt of charged Particles not from PV  (for Pu corrections)

  Double_t sumChargedHadronPtR04  [arraySize];
  Double_t sumChargedParticlePtR04[arraySize];
  Double_t sumNeutralHadronEtR04  [arraySize];
  Double_t sumPhotonEtR04         [arraySize];
  Double_t sumPUPtR04             [arraySize];

  Long64_t isHltMatched[arraySize][triggerArraySize];
  Double_t hltEff      [arraySize][triggerArraySize];
  Double_t hltPt       [arraySize][triggerArraySize];
  Double_t hltEta      [arraySize][triggerArraySize];
  Double_t hltPhi      [arraySize][triggerArraySize];

  static TString getVarString() {
    TString r =
      TString("nMuons/L:") +
      TString("isTracker[N]/L:isStandAlone[N]/L:isGlobal[N]/L:") +
      TString("isTightID[N]/L:isMediumID[N]/L:isLooseID[N]/L:") +
      TString("charge[N]/L:pt[N]/D:ptErr[N]/D:eta[N]/D:phi[N]/D:") +
      TString("trkPt[N]/D:trkPtErr[N]/D:trkEta[N]/D:trkPhi[N]/D:") +
      TString("d0_BS[N]/D:dz_BS[N]/D:") +
      TString("d0_PV[N]/D:dz_PV[N]/D:") +
      TString("relIso[N]/D:") +
      TString("relCombIso[N]/D:") +
      TString("trackIsoSumPt[N]/D:") +
      TString("trackIsoSumPtCorr[N]/D:") +
      TString("hcalIso[N]/D:") +
      TString("ecalIso[N]/D:") +
      TString("isPF[N]/L:") +
      TString("pfPt[N]/D:") +
      TString("pfEta[N]/D:") +
      TString("pfPhi[N]/D:") +
      TString("sumChargedHadronPtR03[N]/D:") +
      TString("sumChargedParticlePtR03[N]/D:") +
      TString("sumNeutralHadronEtR03[N]/D:") +
      TString("sumPhotonEtR03[N]/D:") +
      TString("sumPUPtR03[N]/D:") +
      TString("sumChargedHadronPtR04[N]/D:") +
      TString("sumChargedParticlePtR04[N]/D:") +
      TString("sumNeutralHadronEtR04[N]/D:") +
      TString("sumPhotonEtR04[N]/D:") +
      TString("sumPUPtR04[N]/D:") +
      TString("isHltMatched[N][T]/L:") +
      TString("hltEff[N][T]/D:") +
      TString("hltPt[N][T]/D:") +
      TString("hltEta[N][T]/D:") +
      TString("hltPhi[N][T]/D");
    
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    r.ReplaceAll("[T]",Form("[%d]", triggerArraySize));
    return r;
  };

  void init();
  void fill( const pat::MuonCollection muonsSelected,
	     const reco::Vertex primaryVertex, const int _nPV,
	     const edm::Handle<reco::BeamSpot>& beamSpotHandle,
	     const edm::Event& iEvent, const edm::EventSetup& iSetup,
	     const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
	     const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
	     const std::vector<std::string> _triggerNames, const double _muon_trig_dR,
	     const bool _muon_use_pfIso, const double _muon_iso_dR );

  pat::MuonCollection select( const edm::Handle<pat::MuonCollection>& muons,
			      const reco::Vertex primaryVertex, const std::string _muon_ID,
			      const double _muon_pT_min, const double _muon_eta_max, const double _muon_trig_dR,
			      const bool _muon_use_pfIso, const double _muon_iso_dR, const double _muon_iso_max );

  bool IsLoose ( const pat::Muon muon );
  bool IsMedium( const pat::Muon muon );
  bool IsTight ( const pat::Muon muon, const reco::Vertex primaryVertex );

  double CalcRelIsoPF ( const pat::Muon muon, const double _muon_iso_dR );
  double CalcRelIsoTrk( const pat::Muon muon, const double _muon_iso_dR );
  double CalcTrigEff  ( const pat::Muon muon, const int _nPV, const std::string _triggerName );
  
  bool IsHltMatched( const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup,
		     const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
		     const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
		     const std::string _desiredTriggerName, const double _muon_trig_dR );
  
};


struct _DiMuonInfo {

  const static unsigned int arraySize = 6; // With 4 muons
  Long64_t nPairs;

  Long64_t iMu1[arraySize];
  Long64_t iMu2[arraySize];

  Double_t mass[arraySize];
  Double_t pt  [arraySize];
  Double_t eta [arraySize];
  Double_t y   [arraySize];
  Double_t phi [arraySize];

  Double_t pfMass[arraySize];
  Double_t pfPt  [arraySize];
  Double_t pfEta [arraySize];
  Double_t pfY   [arraySize];
  Double_t pfPhi [arraySize];

  Double_t angle [arraySize];

  static TString getVarString() {
    TString r = TString("nPairs/L:iMu1[N]/L:iMu2[N]/L:") + 
      TString("mass[N]/D:pt[N]/D:eta[N]/D:y[N]/D:phi[N]/D:") +
      TString("pfMass[N]/D:pfPt[N]/D:pfEta[N]/D:pfY[N]/D:pfPhi[N]/D:angle[N]/D");
    
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  };
  
  void init();
  void fill(const _MuonInfo _muonInfo);
  static bool smaller_dMass(std::pair< double, std::pair<int, int> > i,
			    std::pair< double, std::pair<int, int> > j);

};


// Electron info
struct _ElectronInfo {

  const static unsigned int arraySize = 10;
  Long64_t nElectrons;
  
  // Electron cut based IDs
  Long64_t isTightID         [arraySize];
  Long64_t isMediumID        [arraySize];
  Long64_t isLooseID         [arraySize];
  Long64_t isVetoID          [arraySize];
  Long64_t passConversionVeto[arraySize];

  Long64_t charge[arraySize];
  Double_t pt    [arraySize];
  Double_t eta   [arraySize];
  Double_t phi   [arraySize];

  Double_t d0_PV           [arraySize];
  Double_t dz_PV           [arraySize];
  Double_t missingInnerHits[arraySize];
 
  Long64_t isPF[arraySize]; 

  Double_t relIso               [arraySize];
  Double_t sumChargedHadronPtR03[arraySize];  // sum-pt of charged Hadron 
  Double_t sumNeutralHadronEtR03[arraySize];  // sum pt of neutral hadrons
  Double_t sumPhotonEtR03       [arraySize];  // sum pt of PF photons
  Double_t sumPUPtR03           [arraySize];  // sum pt of charged Particles not from PV (for Pu corrections)

  static TString getVarString() {
    TString r = 
      TString("nElectrons/L:") +
      TString("isTightID[N]/L:isMediumID[N]/L:isLooseID[N]/L:isVetoID[N]/L:") +
      TString("passConversionVeto[N]/L:charge[N]/L:pt[N]/D:eta[N]/D:phi[N]/D:") +
      TString("d0_PV[N]/D:dz_PV[N]/D:missingInnerHits[N]/D:") +
      TString("isPF[N]/L:") +
      TString("relIso[N]/D:") + 
      TString("sumChargedHadronPtR03[N]/D:") +
      TString("sumNeutralHadronEtR03[N]/D:") +
      TString("sumPhotonEtR03[N]/D:") +
      TString("sumPUPtR03[N]/D");
    
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  };

  void init();
  void fill( const pat::ElectronCollection electronsSelected,
	     const reco::Vertex primaryVertex, const edm::Event& iEvent,
	     const edm::Handle< edm::ValueMap<bool> >& ele_id_veto, const edm::Handle< edm::ValueMap<bool> >& ele_id_loose,
	     const edm::Handle< edm::ValueMap<bool> >& ele_id_medium, const edm::Handle< edm::ValueMap<bool> >& ele_id_tight );

  pat::ElectronCollection select( const edm::Handle<edm::View<pat::Electron>>& electrons, const reco::Vertex primaryVertex,
				  const edm::Handle< edm::ValueMap<bool> >& ele_id_veto, const edm::Handle< edm::ValueMap<bool> >& ele_id_loose,
				  const edm::Handle< edm::ValueMap<bool> >& ele_id_medium, const edm::Handle< edm::ValueMap<bool> >& ele_id_tight,
				  const std::string _electron_ID, const double _electron_pT_min, const double _electron_eta_max );

  bool   PassKinematics( const pat::Electron ele, const reco::Vertex primaryVertex );
  double CalcRelIsoPF_DeltaBeta( const pat::Electron ele );

};

// Tau info
struct _TauInfo {

  const static unsigned int arraySize   = 10;
  const static unsigned int idArraySize = 16;
  Long64_t nTaus;

  // Tau IDs
  Double_t tauID[arraySize][idArraySize];
  
  Long64_t charge[arraySize];
  Double_t pt    [arraySize];
  Double_t eta   [arraySize];
  Double_t phi   [arraySize];

  Double_t d0_PV[arraySize];
  Double_t dz_PV[arraySize];
 
  Long64_t isPF[arraySize]; 

  static TString getVarString() {
    TString r = 
      TString("nTaus/L:tauID[N][I]/D:") +
      TString("charge[N]/L:pt[N]/D:eta[N]/D:phi[N]/D:") +
      TString("d0_PV[N]/D:dz_PV[N]/D:") +
      TString("isPF[N]/L");
    
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    r.ReplaceAll("[I]",Form("[%d]", idArraySize));
    return r;
  };

  void init();
  void fill( const pat::TauCollection tausSelected,
	     const std::vector<std::string> _tauIDNames);
  
  pat::TauCollection select( const edm::Handle<pat::TauCollection>& taus,
			     const double _tau_pT_min, const double _tau_eta_max );
  
  
};

// MET
struct _MetInfo {
  Double_t px;
  Double_t py;
  Double_t pt;
  Double_t phi;
  Double_t sumEt;

  static TString getVarString() {
    return TString("px/D:py/D:pt/D:phi/D:sumEt/D");
  };

  void init();
  void fill( const edm::Handle < pat::METCollection >& mets, const edm::Event& iEvent );
  
};


// Jets
struct _JetInfo {

  const static unsigned int arraySize = 10;
  Long64_t nJets;
  Long64_t nJetsCent;
  Long64_t nJetsFwd;
  Double_t px           [arraySize];
  Double_t py           [arraySize];
  Double_t pz           [arraySize];
  Double_t pt           [arraySize];
  Double_t eta          [arraySize];
  Double_t phi          [arraySize];
  Double_t mass         [arraySize];
  Double_t charge       [arraySize];
  Double_t isB          [arraySize];
  Long64_t partonFlavour[arraySize];

  /////// Energy Fractions //////
  Double_t chf [arraySize];  // Charged Hadron
  Double_t nhf [arraySize];  // NeutralHadron
  Double_t cef [arraySize];  // Charged EM
  Double_t nef [arraySize];  // Neutral EM
  Double_t muf [arraySize];  // Mu
  Double_t hfhf[arraySize];  // HF Hadron Fraction
  Double_t hfef[arraySize];  // HF EM Fraction

  /////// Multiplicities //////
  Long64_t cm  [arraySize];  // Total Charged Mult
  Long64_t chm [arraySize];  // Charged Hadron Mult
  Long64_t nhm [arraySize];  // NeutralHadron Mult
  Long64_t cem [arraySize];  // Charged EM Mult
  Long64_t nem [arraySize];  // Neutral EM Mult
  Long64_t mum [arraySize];  // Mu Mult
  Long64_t hfhm[arraySize];  // HF Hadron Mult
  Long64_t hfem[arraySize];  // HF EM Mult

  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  Double_t jecFactor[arraySize];
  Double_t jecUnc   [arraySize];  // Jet Energy Correction Uncertainty

  // Gen Jet Values
  Long64_t genMatched[arraySize];
  Double_t genPx     [arraySize];
  Double_t genPy     [arraySize];
  Double_t genPz     [arraySize];
  Double_t genPt     [arraySize];
  Double_t genEta    [arraySize];
  Double_t genPhi    [arraySize];
  Double_t genMass   [arraySize];

  ///// Gen Jet Energy Fractions ///////
  Double_t genEMF [arraySize];  // EM Fraction
  Double_t genHadF[arraySize];  // Had Fraction
  Double_t genInvF[arraySize];  // Invisible Fraction
  Double_t genAuxF[arraySize];  // Auxiliary Fraction (Undecayed Sigmas, etc.)

  Double_t puid[arraySize];  // PUID

  static TString getVarString() {
    TString r = 
      TString("nJets/L:nJetsCent/L:nJetsFwd/L:px[N]/D:py[N]/D:pz[N]/D:pt[N]/D:eta[N]/D:") +
      TString("phi[N]/D:mass[N]/D:charge[N]/D:isB[N]/D:partonFlavour[N]/L:chf[N]/D:") +
      TString("nhf[N]/D:cef[N]/D:nef[N]/D:muf[N]/D:hfhf[N]/D:hfef[N]/D:") +
      TString("cm[N]/L:chm[N]/L:nhm[N]/L:cem[N]/L:nem[N]/L:mum[N]/L:hfhm[N]/L:") +
      TString("hfem[N]/L:jecFactor[N]/D:jecUnc[N]/D:genMatched[N]/L:genPx[N]/D:genPy[N]/D:") +
      TString("genPz[N]/D:genPt[N]/D:genEta[N]/D:genPhi[N]/D:genMass[N]/D:genEMF[N]/D:") +
      TString("genHadF[N]/D:genInvF[N]/D:genAux[N]/D:puid[N]/D");
    
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  };

  void init();
  void fill(const pat::JetCollection jetsSelected, const std::vector<std::string> _btagNames);

  pat::JetCollection select( const edm::Handle<pat::JetCollection>& jets,
			     const JetCorrectorParameters& JetCorPar, const std::string _JES_syst,
			     const std::string _jet_ID, const double _jet_pT_min, const double _jet_eta_max );
    
};

// Generator level jets
struct _GenJetInfo {
  const static unsigned int arraySize = 10;
  Long64_t nJets;
  Double_t px    [arraySize];
  Double_t py    [arraySize];
  Double_t pz    [arraySize];
  Double_t pt    [arraySize];
  Double_t eta   [arraySize];
  Double_t phi   [arraySize];
  Double_t mass  [arraySize];
  Double_t charge[arraySize];

  static TString getVarString() {
    TString r = TString("nJets/L:px[N]/D:py[N]/D:pz[N]/D:pt[N]/D:eta[N]/D:phi[N]/D:mass[N]/D:charge[N]/D");
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  };

  void init();
  void fill(const edm::Handle < reco::GenJetCollection >& genJets, bool isMC);
  static bool sortByPt(reco::GenJet i, reco::GenJet j);

};

// Generator level candidate
struct _GenPartInfo {
  Double_t charge;
  Double_t mass;
  Double_t pt;
  Double_t ptErr;
  Double_t eta;  // pseudo rapidity
  Double_t y;    // rapidity
  Double_t phi;  // phi

  static TString getVarString() {
    return TString("charge/D:mass/D:pt/D:ptErr/D:eta/D:y/D:phi/D");
  };

  void init();
  void fill(const reco::Candidate& genPart);

};

#endif  // #ifndef DATA_FORMATS
