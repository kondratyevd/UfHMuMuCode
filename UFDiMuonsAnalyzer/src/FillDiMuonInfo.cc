
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _DiMuonInfo::init() {
  nPairs = 0;

  for (unsigned int i = 0; i < arraySize; i++) {

    iMu1[i] = -999;
    iMu2[i] = -999;

    mass[i] = -999;
    pt[i]   = -999;
    eta[i]  = -999;
    y[i]    = -999;
    phi[i]  = -999;
    
    pfMass[i] = -999;
    pfPt[i]   = -999;
    pfEta[i]  = -999;
    pfY[i]    = -999;
    pfPhi[i]  = -999;
    
    angle[i]  = -999;
  } // End loop: for (unsigned int i = 0; i < arraySize; i++)

} // End void _DiMuonInfo::init()

void _DiMuonInfo::fill(const _MuonInfo _muonInfo) {

  if (_muonInfo.nMuons < 2)
    return;

  double const MASS_MUON  = 0.105658367; // GeV/c^2
  double const PDG_MASS_Z = 91.1876;     // GeV/c^2

  TLorentzVector mu1_vec;
  TLorentzVector mu2_vec;
  TLorentzVector pair_vec;

  std::vector< std::pair< double, std::pair<int, int> > > dMass;
  dMass.clear();

  // Sort pairs by mass difference from Z boson (smallest to largest)
  for (int i = 0; i < _muonInfo.nMuons; i++) {
    for (int j = i+1; j < _muonInfo.nMuons; j++) {

      mu1_vec.SetPtEtaPhiM(_muonInfo.pt[i], _muonInfo.eta[i], _muonInfo.phi[i], MASS_MUON);
      mu2_vec.SetPtEtaPhiM(_muonInfo.pt[j], _muonInfo.eta[j], _muonInfo.phi[j], MASS_MUON);
      pair_vec = mu1_vec + mu2_vec;

      double mass_diff = fabs(pair_vec.M() - PDG_MASS_Z);
      if (_muonInfo.charge[i] == _muonInfo.charge[j])
	mass_diff += 100000.;  // Push same-sign pairs to the back of the line
      dMass.push_back(std::make_pair(mass_diff, std::make_pair(i, j)));
    }
  }

  nPairs = dMass.size();
  std::stable_sort( dMass.begin(), dMass.end(), smaller_dMass );

  for (int i = 0; i < int(dMass.size()); i++) {
    
    iMu1[i] = dMass.at(i).second.first;
    iMu2[i] = dMass.at(i).second.second;

    mu1_vec.SetPtEtaPhiM(_muonInfo.pt[iMu1[i]], _muonInfo.eta[iMu1[i]], _muonInfo.phi[iMu1[i]], MASS_MUON);
    mu2_vec.SetPtEtaPhiM(_muonInfo.pt[iMu2[i]], _muonInfo.eta[iMu2[i]], _muonInfo.phi[iMu2[i]], MASS_MUON);
    pair_vec = mu1_vec + mu2_vec;
    
    // Correct trackIsoSumPtCorr for other muon? - AWB 09.11.16

    mass[i] = pair_vec.M();
    pt[i]   = pair_vec.Pt();
    eta[i]  = pair_vec.PseudoRapidity();
    y[i]    = pair_vec.Rapidity();
    phi[i]  = pair_vec.Phi();

    angle[i] = mu1_vec.DeltaR(mu2_vec); // Need to check that this is correct - AWB 10.11.16
    // angle[i] = acos( -mu1.track()->momentum().Dot(mu2.track()->momentum() /
    // 						  mu1.track()->p()/mu2.track()->p()) );
    
    
    if ( _muonInfo.isPF[iMu1[i]] && _muonInfo.isPF[iMu2[i]] ) {      
      mu1_vec.SetPtEtaPhiM(_muonInfo.pfPt[iMu1[i]], _muonInfo.pfEta[iMu1[i]], _muonInfo.pfPhi[iMu1[i]], MASS_MUON);
      mu2_vec.SetPtEtaPhiM(_muonInfo.pfPt[iMu2[i]], _muonInfo.pfEta[iMu2[i]], _muonInfo.pfPhi[iMu2[i]], MASS_MUON);
      pair_vec = mu1_vec + mu2_vec;
      
      pfMass[i] = pair_vec.M();
      pfPt  [i] = pair_vec.Pt();
      
      pfEta [i] = pair_vec.PseudoRapidity();
      pfY   [i] = pair_vec.Rapidity();
      pfPhi [i] = pair_vec.Phi();
    }

  } // End loop: for (int i = 0; i < dMass.size(); i++)
} // End void _DiMuonInfo::fill(const _MuonInfo _muonInfo)

bool _DiMuonInfo::smaller_dMass(std::pair< double, std::pair<int, int> > i, 
				std::pair< double, std::pair<int, int> > j) {
  return (i.first < j.first);
}
