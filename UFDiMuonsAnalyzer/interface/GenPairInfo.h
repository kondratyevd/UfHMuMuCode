
#ifndef GENPAIR_INFO
#define GENPAIR_INFO

#include <vector>
#include "TMath.h"

struct GenPairInfo {

  Int_t iMu1;
  Int_t iMu2;
  Int_t mother_ID;
  Int_t postFSR;

  Double_t mass;
  Double_t pt  ;
  Double_t eta ;
  Double_t y   ;
  Double_t phi ;
  Double_t angle;

  void init();

};

typedef std::vector<GenPairInfo> GenPairInfos;

#endif  // #ifndef GENPAIR_INFO
