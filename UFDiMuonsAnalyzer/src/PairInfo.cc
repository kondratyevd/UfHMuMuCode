
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/PairInfo.h"

void PairInfo::init() {

  iMu1 = -999;
  iMu2 = -999;
  
  mass = -999;
  pt   = -999;
  eta  = -999;
  y    = -999;
  phi  = -999;
  
  pfMass = -999;
  pfPt   = -999;
  pfEta  = -999;
  pfY    = -999;
  pfPhi  = -999;
  
  angle  = -999;

} // End void PairInfo::init()

void PairInfos::init() {

  nPairs = 0;
  pairs.clear();

}
