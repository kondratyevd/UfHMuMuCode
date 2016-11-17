
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _GenPartInfo::init() {
  charge = -999;
  mass   = -999;
  pt     = -999;
  ptErr  = -999;
  eta    = -999;
  y      = -999;
  phi    = -999;
}

void _GenPartInfo::fill(const reco::Candidate& genPart) {
  charge = genPart.charge();
  mass   = genPart.mass();
  pt     = genPart.pt();
  // ptErr would need to be filled with pat::Track
  eta    = genPart.eta();
  y      = genPart.rapidity();
  phi    = genPart.phi();
}

