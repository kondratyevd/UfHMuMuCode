// muon info
typedef struct {

  int isTracker;
  int isStandAlone;
  int isGlobal;

  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

  float trkPt;
  float trkPtErr;
  float trketa;
  float trkPhi;
  
  float normChiSquare;
  float d0;
  float dz;
  
  int numPixelLayers;   //number of pixel layers with valid hits
  int numTrackerLayers; //number of tracker layers with valid hits 
  int numStripLayers;   //number of strip layers with valid hits

  float validFracTracker; //valid fraction of tracker hits

  int numValidMuonHits;
  int numValidPixelHits;
  int numValidTrackerHits;
  int numValidStripHits;
  int numSegmentMatches;
  int numOfMatchedStations;

  float trackIsoSumPt;
  float trackIsoSumPtCorr;
  float hcalIso;
  float ecalIso;
  float relCombIso;

  // PF information
  int isPFMuon;

  float pfPt;
  float pfEta;
  float pfPhi;
  
  float sumChargedHadronPtR03; // sum-pt of charged Hadron 
  float sumChargedParticlePtR03; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03;  // sum pt of neutral hadrons
  float sumPhotonEtR03;  // sum pt of PF photons
  float sumPUPtR03;  // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04; 
  float sumChargedParticlePtR04;
  float sumNeutralHadronEtR04;  
  float sumPhotonEtR04;
  float sumPUPtR04;

  int isHltMatched[3];
  float hltPt[3];
  float hltEta[3];
  float hltPhi[3];
  
} _MuonInfo;

// basic track info for true and reco 
typedef struct {
  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;
} _TrackInfo;

