
#ifndef PAIR_HELPER
#define PAIR_HELPER

#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/PairInfo.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/MuonInfo.h"

void FillPairInfos( PairInfos& _pairInfos, const MuonInfos _muonInfos );

// static
bool pair_smaller_dMass(std::pair< double, std::pair<int, int> > i,
			std::pair< double, std::pair<int, int> > j);

#endif  // #ifndef PAIR_HELPER
