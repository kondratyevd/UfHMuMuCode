
#ifndef GENPAIR_HELPER
#define GENPAIR_HELPER

#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/GenPairInfo.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/GenMuonInfo.h"

void FillGenPairInfos( GenPairInfos& _genPairInfos, const GenMuonInfos _muonInfos );

// static
bool genPair_smaller_dMass(std::pair< double, std::pair<int, int> > i,
			   std::pair< double, std::pair<int, int> > j);

#endif  // #ifndef GENPAIR_HELPER
