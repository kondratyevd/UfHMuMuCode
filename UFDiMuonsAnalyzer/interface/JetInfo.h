
#ifndef JET_INFO
#define JET_INFO

#include <vector>
#include "TMath.h"

struct JetInfo {

  Float_t px           ;
  Float_t py           ;
  Float_t pz           ;
  Float_t pt           ;
  Float_t eta          ;
  Float_t phi          ;
  Float_t mass         ;
  Float_t charge       ;
  Bool_t  isB          ;
  Int_t   partonFlavour;

  /////// Energy Fractions //////
  Float_t chf ;  // Charged Hadron
  Float_t nhf ;  // NeutralHadron
  Float_t cef ;  // Charged EM
  Float_t nef ;  // Neutral EM
  Float_t muf ;  // Mu
  Float_t hfhf;  // HF Hadron Fraction
  Float_t hfef;  // HF EM Fraction

  /////// Multiplicities //////
  Int_t cm  ;  // Total Charged Mult
  Int_t chm ;  // Charged Hadron Mult
  Int_t nhm ;  // NeutralHadron Mult
  Int_t cem ;  // Charged EM Mult
  Int_t nem ;  // Neutral EM Mult
  Int_t mum ;  // Mu Mult
  Int_t hfhm;  // HF Hadron Mult
  Int_t hfem;  // HF EM Mult

  // Jet Correction Factor--Above momentum is already corrected!!
  // This factor will return momentum to uncorrected value!!
  Float_t jecFactor;
  Float_t jecUnc   ;  // Jet Energy Correction Uncertainty

  // Gen Jet Values
  Bool_t  genMatched;
  Float_t genPx     ;
  Float_t genPy     ;
  Float_t genPz     ;
  Float_t genPt     ;
  Float_t genEta    ;
  Float_t genPhi    ;
  Float_t genMass   ;

  ///// Gen Jet Energy Fractions ///////
  Float_t genEMF ;  // EM Fraction
  Float_t genHadF;  // Had Fraction
  Float_t genInvF;  // Invisible Fraction
  Float_t genAuxF;  // Auxiliary Fraction (Undecayed Sigmas, etc.)

  Float_t puid;  // PUID

  void init();

};

struct JetInfos {

  Int_t nJets;
  Int_t nJetsCent;
  Int_t nJetsFwd;
  std::vector<JetInfo> jets;

  void init();

};

#endif  // #ifndef JET_INFO
