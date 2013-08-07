#ifndef helpers_h
#define helpers_h
#include "DataFormats.h"
#include <algorithm>
#include "TRandom.h"
#include "TLorentzVector.h"

enum PUJetID
{
  puJetTight  = 0,
  puJetMedium = 1,
  puJetLoose  = 2
};

enum SelectionCodes
{
  notSelected = 0,

  vbfPresel = 1,
  vbfPresel_passVBFBDTCut = 2,
  vbfPresel_isBB_passVBFBDTCut = 4,
  vbfPresel_isNotBB_passVBFBDTCut = 8,

  incPresel = 16,
  incPresel_passIncBDTCut = 32,
  incPresel_isBB_passIncBDTCut = 64,
  incPresel_isBO_passIncBDTCut = 128,
  incPresel_isBE_passIncBDTCut = 256,
  incPresel_isOO_passIncBDTCut = 512,
  incPresel_isOE_passIncBDTCut = 1024,
  incPresel_isEE_passIncBDTCut = 2048,
  incPresel_isNotBB_passIncBDTCut = 4096
};



float getRelIso(_MuonInfo& muon);

bool isKinTight_2012(_MuonInfo& muon);

bool isKinTight_2011(_MuonInfo& muon);

bool isKinTight_2012_noAcc_noIso(_MuonInfo& muon);

bool isKinTight_2011_noAcc_noIso(_MuonInfo& muon);

bool passPUJetID(int flag, PUJetID desiredLevel);

float smearMC(float trueVal, float recoVal, float calib, float smearRatio,TRandom random, bool debug=false);

bool isHltMatched(_MuonInfo& muon1, _MuonInfo& muon2, std::vector<int> allowedPaths);
bool isHltMatched(_MuonInfo& muon1, std::vector<int> allowedPaths);

float resolutionBias(float eta);
float jerCorr   (float ptold, float oldgenpt, float etaold);
float corrPtUp  (float ptold, float oldgenpt, float etaold);
float corrPtDown(float ptold, float oldgenpt, float etaold);

int whichSelection(_MuonInfo& mu1, _MuonInfo& mu2,
                   std::vector<int> allowedHLTPaths,
                   std::string& runPeriod,
                   _PFJetInfo jets,
                   bool passIncBDTCut,
                   bool passVBFBDTCut,
                   double sigmasJEC=0,
                   double jerUncertainty=0);

void setStyle();

#endif
