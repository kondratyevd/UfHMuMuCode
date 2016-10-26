// event info

struct _EventInfo{

  // data fields
  int run;
  int lumi;
  long long int event;
  int bx;
  int orbit;

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      return TString("run/I:lumi/I:event/L:bx/I:orbit/I");
  };

};


// vertex info
struct _VertexInfo{

  // data fields
  const static unsigned int arraySize = 20;
  int nVertices;
  int isValid[arraySize];
  float x[arraySize];	
  float y[arraySize];	
  float z[arraySize];	
  float xErr[arraySize];	
  float yErr[arraySize];	
  float zErr[arraySize];	
  float chi2[arraySize];
  int ndf[arraySize];
  float normChi2[arraySize];

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r =   TString("nVertices/I:isValid[N]/I:")+
                    TString("x[N]/F:y[N]/F:z[N]/F:xErr[N]/F:yErr[N]/F:zErr[N]/F:")+
                    TString("chi2[N]/F:ndf[N]/I:normChi2[N]/F");
      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      return r;
  };

};


// basic track info
struct _TrackInfo{
  int charge;
  float pt;
  float ptErr;
  float eta;
  float phi;

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      return TString("charge/I:pt/F:ptErr/F:eta/F:phi/F");
  };

};

// muon info
struct _MuonInfo{

  const static unsigned int arraySize = 10;
  const static unsigned int triggerArraySize = 6;
  int nMuons;
  int nMuonPairs;

  int isTracker[arraySize];
  int isStandAlone[arraySize];
  int isGlobal[arraySize];

  int isTightMuon[arraySize];
  int isMediumMuon[arraySize];
  int isLooseMuon[arraySize];

  int charge[arraySize];
  float pt[arraySize];
  float ptErr[arraySize];
  float eta[arraySize];
  float phi[arraySize];

  float trkPt[arraySize];
  float trkPtErr[arraySize];
  float trketa[arraySize];
  float trkPhi[arraySize];

  float d0_BS[arraySize];
  float dz_BS[arraySize];

  float d0_PV[arraySize];
  float dz_PV[arraySize];

  float trackIsoSumPt[arraySize];
  float trackIsoSumPtCorr[arraySize];
  float hcalIso[arraySize];
  float ecalIso[arraySize];
  float relCombIso[arraySize];

  // PF information
  int isPFMuon[arraySize];

  float pfPt[arraySize];
  float pfEta[arraySize];
  float pfPhi[arraySize];
  
  float sumChargedHadronPtR03[arraySize];   // sum-pt of charged Hadron 
  float sumChargedParticlePtR03[arraySize]; // sum-pt of charged Particles(inludes e/mu) 
  float sumNeutralHadronEtR03[arraySize];   // sum pt of neutral hadrons
  float sumPhotonEtR03[arraySize];          // sum pt of PF photons
  float sumPUPtR03[arraySize];              // sum pt of charged Particles not from PV  (for Pu corrections)

  float sumChargedHadronPtR04[arraySize]; 
  float sumChargedParticlePtR04[arraySize];
  float sumNeutralHadronEtR04[arraySize];  
  float sumPhotonEtR04[arraySize];
  float sumPUPtR04[arraySize];

  int isHltMatched[arraySize][triggerArraySize];
  float hltPt[arraySize][triggerArraySize];
  float hltEta[arraySize][triggerArraySize];
  float hltPhi[arraySize][triggerArraySize];

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r =  TString("nMuons/I:nMuonPairs/I:")+
                   TString("isTracker[N]/I:isStandAlone[N]/I:isGlobal[N]/I:")+
                   TString("isTightMuon[N]/I:isMediumMuon[N]/I:isLooseMuon[N]/I:")+
                   TString("charge[N]/I:pt[N]/F:ptErr[N]/F:eta[N]/F:phi[N]/F:")+
                   TString("trkPt[N]/F:trkPtErr[N]/F:trkEta[N]/F:trkPhi[N]/F:")+
                   TString("d0_BS[N]/F:dz_BS[N]/F:")+
                   TString("d0_PV[N]/F:dz_PV[N]/F:")+
                   TString("trackIsoSumPt[N]/F:")+
                   TString("trackIsoSumPtCorr[N]/F:")+
                   TString("hcalIso[N]/F:")+
                   TString("ecalIso[N]/F:")+
                   TString("relCombIso[N]/F:")+
                   TString("isPFMuon[N]/I:")+
                   TString("pfPt[N]/F:")+
                   TString("pfEta[N]/F:")+
                   TString("pfPhi[N]/F:")+
                   TString("sumChargedHadronPtR03[N]/F:")+
                   TString("sumChargedParticlePtR03[N]/F:")+
                   TString("sumNeutralHadronEtR03[N]/F:")+
                   TString("sumPhotonEtR03[N]/F:")+
                   TString("sumPUPtR03[N]/F:")+
                   TString("sumChargedHadronPtR04[N]/F:")+
                   TString("sumChargedParticlePtR04[N]/F:")+
                   TString("sumNeutralHadronEtR04[N]/F:")+
                   TString("sumPhotonEtR04[N]/F:")+
                   TString("sumPUPtR04[N]/F:")+
                   TString("isHltMatched[N][T]/I:")+
                   TString("hltPt[N][T]/F:")+
                   TString("hltEta[N][T]/F:")+
                   TString("hltPhi[N][T]/F");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      r.ReplaceAll("[T]",Form("[%d]", triggerArraySize));
      return r;
  };

};

// electron info
struct _ElectronInfo{

  const static unsigned int arraySize = 10;
  int nElectrons;

  // electron cut based IDs
  int isTightElectron[arraySize];
  int isMediumElectron[arraySize];
  int isLooseElectron[arraySize];
  int isVetoElectron[arraySize];
  int passConversionVeto[arraySize];

  int charge[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];

  float d0_PV[arraySize];
  float dz_PV[arraySize];
  float missingInnerHits[arraySize];
 
  int isPFElectron[arraySize]; 

  float sumChargedHadronPtR03[arraySize];   // sum-pt of charged Hadron 
  float sumNeutralHadronEtR03[arraySize];   // sum pt of neutral hadrons
  float sumPhotonEtR03[arraySize];          // sum pt of PF photons
  float sumPUPtR03[arraySize];              // sum pt of charged Particles not from PV  (for Pu corrections)

  // tell the ttree how to save the struct
  static TString getVarString()
  {
      TString r = TString("nElectrons/I:")+
                  TString("isTightElectron[N]/I:isMediumElectron[N]/I:isLooseElectron[N]/I:isVetoElectron[N]/I:")+
                  TString("passConversionVeto[N]/I:charge[N]/I:pt[N]/F:eta[N]/F:phi[N]/F:")+
                  TString("d0_PV[N]/F:dz_PV[N]/F:missingInnerHits[N]/F:")+
                  TString("isPFElectron[N]/I:")+
                  TString("sumChargedHadronPtR03[N]/F:")+
                  TString("sumNeutralHadronEtR03[N]/F:")+
                  TString("sumPhotonEtR03[N]/F:")+
                  TString("sumPUPtR03[N]/F");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      return r;
  };
};

//MET
struct _MetInfo{
  float px;
  float py;
  float pt;
  float phi;
  float sumEt;

  static TString getVarString()
  {
      return TString("px/F:py/F:pt/F:phi/F:sumEt/F");
  };

};


// pf Jets
struct _PFJetInfo{

  const static unsigned int arraySize = 10;
  int   nJets;
  float px[arraySize];
  float py[arraySize];
  float pz[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];
  float mass[arraySize];
  int   charge[arraySize];
  float isB[arraySize];
  int   partonFlavour[arraySize];
  /////// Energy Fractions //////
  //Charged Hadron
  float chf[arraySize];
  //NeutralHadron
  float nhf[arraySize];
  //Charged EM
  float cef[arraySize];
  //Neutral EM
  float nef[arraySize];
  //Mu
  float muf[arraySize];
  // HF Hadron Fraction
  float hfhf[arraySize];
  // HF EM Fraction
  float hfef[arraySize];
  /////// Multiplicities //////
  // Total Charged Mult
  int cm[arraySize];
  //Charged Hadron Mult
  int chm[arraySize];
  //NeutralHadron Mult
  int nhm[arraySize];
  //Charged EM Mult
  int cem[arraySize];
  //Neutral EM Mult
  int nem[arraySize];
  //Mu Mult
  int mum[arraySize];
  // HF Hadron Mult
  int hfhm[arraySize];
  // HF EM Mult
  int hfem[arraySize];
  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  float jecFactor[arraySize];
  // Jet Energy Correction Uncertainty
  float jecUnc[arraySize];
  // b-Tag
  float csv[arraySize];
  // Gen Jet Values
  bool genMatched[arraySize];
  float genPx[arraySize];
  float genPy[arraySize];
  float genPz[arraySize];
  float genPt[arraySize];
  float genEta[arraySize];
  float genPhi[arraySize];
  float genMass[arraySize];
  ///// Gen Jet Energy Fractions ///////
  // EM Fraction
  float genEMF[arraySize];
  // Had Fraction
  float genHadF[arraySize];
  // Invisible Fraction
  float genInvF[arraySize];
  // Auxiliary Fraction (Undecayed Sigmas, etc.)
  float genAuxF[arraySize];
  // PUID
  float puid[arraySize];

  static TString getVarString()
  {
       TString r = TString("nJets/I:px[N]/F:py[N]/F:pz[N]/F:pt[N]/F:eta[N]/F:")+
                   TString("phi[N]/F:mass[N]/F:charge[N]/I:isB[N]/F:partonFlavour[N]/I:chf[N]/F:")+
                   TString("nhf[N]/F:cef[N]/F:nef[N]/F:muf[N]/F:hfhf[N]/F:hfef[N]/F:")+
                   TString("cm[N]/I:chm[N]/I:nhm[N]/I:cem[N]/I:nem[N]/I:mum[N]/I:hfhm[N]/I:")+
                   TString("hfem[N]/I:jecFactor[N]/F:jecUnc[N]/F:csv[N]/F:genMatched[N]/F:genPx[N]/F:genPy[N]/F:")+
                   TString("genPz[N]/F:genPt[N]/F:genEta[N]/F:genPhi[N]/F:genMass[N]/F:genEMF[N]/F:")+
                   TString("genHadF[N]/F:genInvF[N]/F:genAux[N]/F:puid[N]:F");

      r.ReplaceAll("[N]",Form("[%d]", arraySize));
      return r;
  };

};

// generator level jets
struct _GenJetInfo{
  const static unsigned int arraySize = 10;
  int nJets;
  float px[arraySize];
  float py[arraySize];
  float pz[arraySize];
  float pt[arraySize];
  float eta[arraySize];
  float phi[arraySize];
  float mass[arraySize];
  int   charge[arraySize];

  static TString getVarString()
  {
    TString r = TString("nJets/I:px[N]/F:py[N]/F:pz[N]/F:pt[N]/F:eta[N]/F:phi[N]/F:mass[N]/F:charge[N]/I");
    r.ReplaceAll("[N]",Form("[%d]", arraySize));
    return r;
  }
};

// generator level composite Candidate
struct _genPartInfo{
  int charge;
  float mass;
  float pt;
  float eta;  // pseudo rapidity
  float y;    // rapidity
  float phi;  // phi

  static TString getVarString()
  {
      return TString("charge/I:mass/F:pt/F:eta/F:y/F:phi/F");
  };
};
