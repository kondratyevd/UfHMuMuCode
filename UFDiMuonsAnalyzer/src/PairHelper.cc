
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/PairHelper.h"

void FillPairInfos( PairInfos& _pairInfos, const MuonInfos _muonInfos ) {

  _pairInfos.clear();
  if (_muonInfos.size() < 2)
    return;

  double const MASS_MUON  = 0.105658367; // GeV/c^2
  double const PDG_MASS_Z = 91.1876;     // GeV/c^2

  TLorentzVector mu1_vec;
  TLorentzVector mu2_vec;
  TLorentzVector pair_vec;

  std::vector< std::pair< double, std::pair<int, int> > > dMass;
  dMass.clear();

  // Sort pairs by mass difference from Z boson (smallest to largest)
  for (int i = 0; i < int(_muonInfos.size()); i++) {
    for (int j = i+1; j < int(_muonInfos.size()); j++) {

      mu1_vec.SetPtEtaPhiM(_muonInfos.at(i).pt, _muonInfos.at(i).eta, _muonInfos.at(i).phi, MASS_MUON);
      mu2_vec.SetPtEtaPhiM(_muonInfos.at(j).pt, _muonInfos.at(j).eta, _muonInfos.at(j).phi, MASS_MUON);
      pair_vec = mu1_vec + mu2_vec;

      double mass_diff = fabs(pair_vec.M() - PDG_MASS_Z);
      if (_muonInfos.at(i).charge == _muonInfos.at(j).charge)
	mass_diff += 100000.;  // Push same-sign pairs to the back of the line
      dMass.push_back(std::make_pair(mass_diff, std::make_pair(i, j)));
    }
  }

  int nPairs = dMass.size();
  std::stable_sort( dMass.begin(), dMass.end(), pair_smaller_dMass );

  for (int i = 0; i < int(dMass.size()); i++) {

    PairInfo _pairInfo;
    _pairInfo.init();
    
    int iMu1 = dMass.at(i).second.first;
    int iMu2 = dMass.at(i).second.second;

    _pairInfo.iMu1 = iMu1;
    _pairInfo.iMu2 = iMu2;

    mu1_vec.SetPtEtaPhiM(_muonInfos.at(iMu1).pt, _muonInfos.at(iMu1).eta, _muonInfos.at(iMu1).phi, MASS_MUON);
    mu2_vec.SetPtEtaPhiM(_muonInfos.at(iMu2).pt, _muonInfos.at(iMu2).eta, _muonInfos.at(iMu2).phi, MASS_MUON);
    pair_vec = mu1_vec + mu2_vec;
    
    // Correct trackIsoSumPtCorr for other muon? - AWB 09.11.16

    _pairInfo.mass = pair_vec.M();
    _pairInfo.pt   = pair_vec.Pt();
    _pairInfo.eta  = pair_vec.PseudoRapidity();
    _pairInfo.y    = pair_vec.Rapidity();
    _pairInfo.phi  = pair_vec.Phi();

    _pairInfo.angle = mu1_vec.DeltaR(mu2_vec); // Need to check that this is correct - AWB 10.11.16
    // _pairInfo.angle = acos( -mu1.track()->momentum().Dot(mu2.track()->momentum() /
    // 						  mu1.track()->p()/mu2.track()->p()) );
    
    
    if ( _muonInfos.at(iMu1).isPF && _muonInfos.at(iMu2).isPF ) {      
      mu1_vec.SetPtEtaPhiM(_muonInfos.at(iMu1).pfPt, _muonInfos.at(iMu1).pfEta, _muonInfos.at(iMu2).pfPhi, MASS_MUON);
      mu2_vec.SetPtEtaPhiM(_muonInfos.at(iMu2).pfPt, _muonInfos.at(iMu2).pfEta, _muonInfos.at(iMu2).pfPhi, MASS_MUON);
      pair_vec = mu1_vec + mu2_vec;
      
      _pairInfo.pfMass = pair_vec.M();
      _pairInfo.pfPt   = pair_vec.Pt();
      
      _pairInfo.pfEta  = pair_vec.PseudoRapidity();
      _pairInfo.pfY    = pair_vec.Rapidity();
      _pairInfo.pfPhi  = pair_vec.Phi();
    }

    _pairInfos.push_back( _pairInfo );
  } // End loop: for (int i = 0; i < dMass.size(); i++)
  
  if ( int(_pairInfos.size()) != nPairs )
    std::cout << "Bizzare error: _pairInfos.size() = " << _pairInfos.size()
	      << ", nPairs = " << nPairs << std::endl;
  
} // End void FillPairInfos( PairInfos& _pairInfos, const MuonInfos _muonInfos )

bool pair_smaller_dMass( std::pair< double, std::pair<int, int> > i, 
				std::pair< double, std::pair<int, int> > j) {
  return (i.first < j.first);
}