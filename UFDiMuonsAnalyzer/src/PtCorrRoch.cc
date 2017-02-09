
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/PtCorrRoch.h"

void CorrectPtRoch( const RoccoR _calib, const bool _doSys, const TLorentzVector _mu_vec, 
		    float& _pt, float& _pt_sys_up, float& _pt_sys_down, const int _charge, 
		    const int _trk_layers, const float _GEN_pt, const bool _isData ) {
  
  _pt = _mu_vec.Pt();
  _pt_sys_up = -999;
  _pt_sys_down = -999;
  float q_term     = 1.0;
  float q_term_sys = -99;
  
  srand( time(NULL) );
  int iRand_1 = rand();
  int iRand_2 = (iRand_1 % 1000)*(iRand_1 % 1000);
  float fRand_1 = (iRand_1 % 1000) * 0.001;
  float fRand_2 = (iRand_2 % 1000) * 0.001;
  
  // For default computation, error set and error member are 0
  if (_isData)          q_term = _calib.kScaleDT( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), 0, 0 );
  else if (_GEN_pt > 0) q_term = _calib.kScaleFromGenMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),
							 _trk_layers, _GEN_pt, fRand_1, 0, 0 );
  else                  q_term = _calib.kScaleAndSmearMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),
							  _trk_layers, fRand_1, fRand_2, 0, 0 );
  int nUp   = 0;
  int nDown = 0;
  double sum_sq_up   = 0;
  double sum_sq_down = 0;
  
  // Throw 100 toys to generate uncertainty on correction
  // Even with 100 toys, only increases UFDiMuonsAnalyzer NTuple-maker time from 2 min for 10k events to 2:20
  for (int i = 0; i < 100; i++) {
    if (!_doSys) break;
    
    if (_isData)          q_term_sys = _calib.kScaleDT( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(), 1, i );
    else if (_GEN_pt > 0) q_term_sys = _calib.kScaleFromGenMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),
							       _trk_layers, _GEN_pt, fRand_1, 1, i );
    else                  q_term_sys = _calib.kScaleAndSmearMC( _charge, _mu_vec.Pt(), _mu_vec.Eta(), _mu_vec.Phi(),
								_trk_layers, fRand_1, fRand_2, 1, i );
    if ( q_term_sys >= q_term ) {
      nUp   += 1;
      sum_sq_up   += pow( q_term_sys - q_term, 2 );
    } else {
      nDown += 1;
      sum_sq_down += pow( q_term_sys - q_term, 2 );
    }
  }
  
  _pt          = _mu_vec.Pt() * q_term;
  _pt_sys_up   = ( _doSys ? _pt * sqrt(sum_sq_up   / nUp)   : -999 );
  _pt_sys_down = ( _doSys ? _pt * sqrt(sum_sq_down / nDown) : -999 ); 

  // std::cout << "\n*******   RESULT   *******" << std::endl;
  // std::cout << "_mu_vec.Pt() = " << _mu_vec.Pt() << ", corrected = " << _pt << std::endl;
  // std::cout << "nUp = " << nUp << ", _pt_sys_up = " << _pt_sys_up << std::endl;
  // std::cout << "nDown = " << nDown << ", _pt_sys_down = " << _pt_sys_down << std::endl;
  
}
