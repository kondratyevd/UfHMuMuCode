
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _EventInfo::init() {
  run   = -999;
  lumi  = -999;
  event = -999;
  bx    = -999;
  orbit = -999;
}

void _EventInfo::fill(const edm::Event& iEvent) {
  run   = iEvent.id().run();
  lumi  = iEvent.luminosityBlock();
  event = iEvent.id().event();
  bx    = iEvent.bunchCrossing();
  orbit = iEvent.orbitNumber();
}
