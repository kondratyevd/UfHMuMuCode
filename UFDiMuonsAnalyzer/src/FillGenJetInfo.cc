
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _GenJetInfo::init() {
  nJets = 0;
  for(unsigned int i = 0; i < arraySize; i++) {
    px[i]     = -999;
    py[i]     = -999;
    pz[i]     = -999;
    pt[i]     = -999;
    eta[i]    = -999;
    phi[i]    = -999;
    mass[i]   = -999;
    charge[i] = -999;
  }
}

void _GenJetInfo::fill(const edm::Handle < reco::GenJetCollection >& genJets, bool isMC) { 

  if ( !genJets.isValid() ) {
    if (isMC) std::cout << "No valid collection of GEN jets" << std::endl;
    return;
  }

  reco::GenJetCollection sortedGenJets = (*genJets);
  sort(sortedGenJets.begin(), sortedGenJets.end(), sortByPt);
  for(unsigned int i = 0; i < sortedGenJets.size(); i++) {
    
    nJets++;
    if (i >= arraySize) {
      std::cout << "Found " << i+1 << "th GEN jet; only " << arraySize << " allowed in array" << std::endl;
      return;
    }
    
    px[i]     = sortedGenJets[i].px();
    py[i]     = sortedGenJets[i].py();
    pz[i]     = sortedGenJets[i].pz();
    pt[i]     = sortedGenJets[i].pt();
    eta[i]    = sortedGenJets[i].eta();
    phi[i]    = sortedGenJets[i].phi();
    mass[i]   = sortedGenJets[i].mass();
    charge[i] = sortedGenJets[i].charge();
  }
  
}

bool _GenJetInfo::sortByPt (reco::GenJet i, reco::GenJet j) { 
  return ( i.pt() > j.pt() ); 
}
