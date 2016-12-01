
#ifndef MET_HELPER
#define MET_HELPER

#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/MetInfo.h"


void FillMetInfo( MetInfo& _metInfo,
		  const edm::Handle < pat::METCollection >& mets, 
		  const edm::Event& iEvent );

#endif  // #ifndef MET_HELPER
