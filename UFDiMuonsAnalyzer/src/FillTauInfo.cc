
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _TauInfo::init() {
  nTaus = 0;
  for (unsigned int i = 0; i < arraySize; i++) {

    for(unsigned int id = 0; id < idArraySize; id++)
      tauID[i][id] = -999;
    
    charge[i] = -999;
    pt[i] = -999;
    eta[i] = -999;
    phi[i] = -999;
    
    d0_PV[i] = -999;
    dz_PV[i] = -999;
    
    isPF[i] = -999;
  }
}

void _TauInfo::fill( const pat::TauCollection tausSelected, 
		     const std::vector<std::string> _tauIDNames) {
  
  nTaus = tausSelected.size();
  
  for (int i = 0; i < nTaus; i++) {
    
    if ( i >= int(arraySize) ) {
      std::cout << "Found " << i+1 << "th tau; only " << arraySize << " allowed in array" << std::endl;
      return;
    }
    
    pat::Tau tau = tausSelected.at(i);

    // Basic kinematics
    charge[i]  = tau.charge();
    pt[i]      = tau.pt();
    eta[i]     = tau.eta();
    phi[i]     = tau.phi();

    // Basic quality
    isPF[i] = tau.isPFTau();

    // Basic vertexing?

    // // -------------------------
    // // ONLY SAVE BASIC VARIABLES - AWB 12.11.16
    // // -------------------------
    // if ( _outputLevel == "Slim")
    //   continue;

    // Standard quality
    for (unsigned int id = 0; id < _tauIDNames.size(); id++) {
      if (id < idArraySize)
	tauID[i][id] = tau.tauID(_tauIDNames[id]);
      else
	std::cout << "Found " << id+1 << "th tau ID name; only " << idArraySize << " allowed in array" << std::endl;
    }

    // Standard vertexing
    if ( !(tau.leadChargedHadrCand().isNull()) ) {
      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      
      dz_PV[i] = packedLeadTauCand->dz();
      d0_PV[i] = packedLeadTauCand->dxy();
    }

  } // End loop: for (int i = 0; i < nTaus; i++)
} // End void _TauInfo::fill


pat::TauCollection _TauInfo::select( const edm::Handle<pat::TauCollection>& taus,
                                     const double _tau_pT_min, const double _tau_eta_max ) {
  // Presumably should follow something from here:
  // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV

  pat::TauCollection tausSelected;
  tausSelected.clear();

  if ( !taus.isValid() ) {
    std::cout << "No valid tau collection" << std::endl;
    return tausSelected;
  }

  for ( pat::TauCollection::const_iterator tau = taus->begin(), tausEnd = taus->end(); tau != tausEnd; ++tau) {
    
    if ( tau->pt()          < _tau_pT_min  ) continue;
    if ( fabs( tau->eta() ) > _tau_eta_max ) continue;

    tausSelected.push_back(*tau);
  }

  return tausSelected;
}
