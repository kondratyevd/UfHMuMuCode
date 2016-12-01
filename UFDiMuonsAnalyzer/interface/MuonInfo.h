
#ifndef MUON_INFO
#define MUON_INFO

#include <vector>
#include "TMath.h"

struct MuonInfo {

  Bool_t isTracker   ;
  Bool_t isStandAlone;
  Bool_t isGlobal    ;

  Bool_t isTightID ;
  Bool_t isMediumID;
  Bool_t isLooseID ;

  Int_t    charge;
  Double_t pt    ;
  Double_t ptErr ;
  Double_t eta   ;
  Double_t phi   ;

  Double_t trkPt   ;
  Double_t trkPtErr;
  Double_t trketa  ;
  Double_t trkPhi  ;

  Float_t d0_BS;
  Float_t dz_BS;

  Float_t d0_PV;
  Float_t dz_PV;

  Float_t relIso           ;
  Float_t relCombIso       ;
  Float_t trackIsoSumPt    ;
  Float_t trackIsoSumPtCorr;
  Float_t hcalIso          ;
  Float_t ecalIso          ;

  // PF information
  Bool_t isPF;

  Double_t pfPt ;
  Double_t pfEta;
  Double_t pfPhi;
  
  Float_t sumChargedHadronPtR03  ;  // sum-pt of charged Hadron
  Float_t sumChargedParticlePtR03;  // sum-pt of charged Particles(inludes e/mu)
  Float_t sumNeutralHadronEtR03  ;  // sum pt of neutral hadrons
  Float_t sumPhotonEtR03         ;  // sum pt of PF photons
  Float_t sumPUPtR03             ;  // sum pt of charged Particles not from PV  (for Pu corrections)

  Float_t sumChargedHadronPtR04  ;
  Float_t sumChargedParticlePtR04;
  Float_t sumNeutralHadronEtR04  ;
  Float_t sumPhotonEtR04         ;
  Float_t sumPUPtR04             ;

  /* Bool_t  isHltMatched[triggerArraySize]; */
  /* Float_t hltEff      [triggerArraySize]; */
  /* Float_t hltPt       [triggerArraySize]; */
  /* Float_t hltEta      [triggerArraySize]; */
  /* Float_t hltPhi      [triggerArraySize]; */

  void init();

};

struct MuonInfos {

  Int_t nMuons;
  std::vector<MuonInfo> muons;

  void init();

};

#endif  // #ifndef MUON_INFO
