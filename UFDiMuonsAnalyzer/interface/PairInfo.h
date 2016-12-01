
#ifndef PAIR_INFO
#define PAIR_INFO

#include <vector>
#include "TMath.h"

struct PairInfo {

  Int_t iMu1;
  Int_t iMu2;

  Double_t mass;
  Double_t pt  ;
  Double_t eta ;
  Double_t y   ;
  Double_t phi ;

  Double_t pfMass;
  Double_t pfPt  ;
  Double_t pfEta ;
  Double_t pfY   ;
  Double_t pfPhi ;

  Double_t angle;

  void init();

};

struct PairInfos {
  
  Int_t nPairs;
  std::vector<PairInfo> pairs;

  void init();

};

#endif  // #ifndef PAIR_INFO
