#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/ElePairHelper.h"

void FillElePairInfos( ElePairInfos& _pairInfos, const EleInfos _eleInfos ) {

  _pairInfos.clear();
  if (_eleInfos.size() < 2)
    return;

  double const MASS_ELE  = 0.0005109989461; // GeV/c^2
  double const PI         = 3.14159265359;

  // 4-vectors: nominal, mu1_up, mu1_down, mu2_up, mu2_down
  eleVecSys ele1_vec;
  eleVecSys ele2_vec;
  pairVecSys_ele pair_vec;
  // double massErr;

  std::vector< std::pair< bool, std::pair<int, int> > > isOS;
  isOS.clear();

  // Sort pairs by OS/SS, then highest mu1 pT, then highest mu2 pT
  // Muons come sorted by pT, so only need to stable sort by OS/SS 
  for (int i = 0; i < int(_eleInfos.size()); i++) {
    for (int j = i+1; j < int(_eleInfos.size()); j++) {
      bool osPair = (_eleInfos.at(i).charge + _eleInfos.at(j).charge == 0);
      isOS.push_back(std::make_pair(osPair, std::make_pair(i, j)));
    }
  }

  int nPairs = isOS.size();
  std::stable_sort( isOS.begin(), isOS.end(), pair_is_OS_ele );
  
  for (int i = 0; i < int(isOS.size()); i++) {
    
    ElePairInfo _pairInfo;
    _pairInfo.init();
    
    int iEle1 = isOS.at(i).second.first;
    int iEle2 = isOS.at(i).second.second;
    
    _pairInfo.iEle1 = iEle1;
    _pairInfo.iEle2 = iEle2;
    _pairInfo.charge = _eleInfos.at(iEle1).charge + _eleInfos.at(iEle2).charge; 

    FillElePairMasses( ele1_vec, ele2_vec, pair_vec, 
      // massErr,
       MASS_ELE, 
        _eleInfos.at(iEle1), _eleInfos.at(iEle2),
        _eleInfos.at(iEle1).pt, _eleInfos.at(iEle2).pt
        // ,
        // _eleInfos.at(iEle1).ptErr, _eleInfos.at(iEle2).ptErr 
        );

    _pairInfo.mass    = pair_vec.nom.M();
    // _pairInfo.massErr = massErr;
    _pairInfo.pt      = pair_vec.nom.Pt();
    _pairInfo.eta     = pair_vec.nom.PseudoRapidity();
    _pairInfo.rapid   = pair_vec.nom.Rapidity();
    _pairInfo.phi     = pair_vec.nom.Phi();

    double _dR   = ele1_vec.nom.DeltaR(ele2_vec.nom);
    double _dEta = ele1_vec.nom.PseudoRapidity() - ele2_vec.nom.PseudoRapidity();
    double _dPhi = ele1_vec.nom.DeltaPhi(ele2_vec.nom);

    double _dThetaStarEta = acos( tanh(_dEta/2) );
    double _dPhiStar      = tan( (PI - fabs(_dPhi)) / 2) * sin(_dThetaStarEta);

    _pairInfo.dR   = _dR;
    _pairInfo.dEta = _dEta;
    _pairInfo.dPhi = _dPhi;
    _pairInfo.dPhiStar = _dPhiStar;

    _pairInfos.push_back( _pairInfo );
  } // End loop: for (int i = 0; i < isOS.size(); i++)
  
  if ( int(_pairInfos.size()) != nPairs )
    std::cout << "Bizzare error: electron _pairInfos.size() = " << _pairInfos.size()
        << ", nPairs = " << nPairs << std::endl;
  
} // End void FillElePairInfos( ElePairInfos& _pairInfos, const EleInfos _eleInfos )

bool pair_is_OS_ele( std::pair< bool, std::pair<int, int> > i, 
     std::pair< bool, std::pair<int, int> > j) {
  return (i.first || !j.first);
}

void FillElePairMasses( eleVecSys& ele1_vec, eleVecSys& ele2_vec, pairVecSys_ele& pair_vec, 
         // double& massErr, 
         const double MASS_ELE,
         const EleInfo _ele1, const EleInfo _ele2, 
         const double _ele1_pt, const double _ele2_pt
         // ,
         // const double _ele1_ptErr, const double _ele2_ptErr 
         ) {

  ele1_vec.nom.SetPtEtaPhiM(_ele1_pt, _ele1.eta, _ele1.phi, MASS_ELE);
  ele2_vec.nom.SetPtEtaPhiM(_ele2_pt, _ele2.eta, _ele2.phi, MASS_ELE);
  
  // ele1_vec.up.SetPtEtaPhiM(_ele1_pt + _ele1_ptErr, _ele1.eta, _ele1.phi, MASS_ELE);
  // ele2_vec.up.SetPtEtaPhiM(_ele2_pt + _ele2_ptErr, _ele2.eta, _ele2.phi, MASS_ELE);
  
  // ele1_vec.down.SetPtEtaPhiM(_ele1_pt - _ele1_ptErr, _ele1.eta, _ele1.phi, MASS_ELE);
  // ele2_vec.down.SetPtEtaPhiM(_ele2_pt - _ele2_ptErr, _ele2.eta, _ele2.phi, MASS_ELE);
  
  pair_vec.nom   = ele1_vec.nom  + ele2_vec.nom;
  // pair_vec.up1   = ele1_vec.up   + ele2_vec.nom;
  // pair_vec.down1 = ele1_vec.down + ele2_vec.nom;
  // pair_vec.up2   = ele1_vec.nom  + ele2_vec.up;
  // pair_vec.down2 = ele1_vec.nom  + ele2_vec.down;
  
  // massErr = sqrt( pow(pair_vec.nom.M() - pair_vec.up1.M(),   2) +
  //     pow(pair_vec.nom.M() - pair_vec.down1.M(), 2) + 
  //     pow(pair_vec.nom.M() - pair_vec.up2.M(),   2) + 
  //     pow(pair_vec.nom.M() - pair_vec.down2.M(), 2) ) / 2.;
  
} // End void FillElePairMasses()