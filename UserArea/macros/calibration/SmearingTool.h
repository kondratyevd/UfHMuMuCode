#include <iostream>
#include <TChain.h>
#include <TClonesArray.h>
#include <TString.h>
#include <map>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include "TF1.h"
#include <stdio.h>
#include <math.h>
// class init
class SmearingTool {
 public:
  //SmearingTool();
  //SmearingTool(int seed);
  SmearingTool();
  float PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen, TString ParVar = "null",float ParSig = 0);

 private:
  static const float mean[72];
  static const float sig1[72];
  static const float sig2[72];
  static const float Asig2[72];
  static const float ERRmean[72];
  static const float ERRsig1[72];
  static const float ERRsig2[72];
  static const float ERRAsig2[72];
  static const float ResRMS[72];
  static const float ErrResRMS[72];

  static const float PTbin[10];
  static const float ETAbin[9];

};
// end: class init
