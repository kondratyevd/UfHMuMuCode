
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/PtCorrRoch.h"

void CorrectPtRoch( rochcor2016* _calib[201], const bool _doSys,
		    TLorentzVector& _mu_vec, float& _q_term, 
		    float& _pt_sys_up, float& _pt_sys_down,
		    const int _charge, const int _trk_layers, const bool _isData ) {
  
  _q_term = 1.0;
  _pt_sys_up = -999;
  _pt_sys_down = -999;
  int runopt = 0;   // No run dependence for 2016 prompt, so default runopt = 0 

  TLorentzVector mu_vec_orig = _mu_vec;
  TLorentzVector mu_vec_sys = _mu_vec;
  float q_term_orig = _q_term;
  float q_term_sys = _q_term;

  if (_isData) _calib[0]->momcor_data( _mu_vec, _charge*1.0, runopt, _q_term );
  else         _calib[0]->momcor_mc  ( _mu_vec, _charge*1.0, _trk_layers, _q_term );

  int nUp   = 0;
  int nDown = 0;
  double sum_sq_up   = 0;
  double sum_sq_down = 0;
 
  // Throw 200 toys to generate uncertainty on correction
  // Even with 200 toys, does not add noticeable time to the UFDiMuonsAnalyzer NTuple-maker
  for (int i = 1; i < 201; i++) {
    if (!_doSys) break;

    if (_isData) _calib[i]->momcor_data( mu_vec_sys, _charge*1.0, runopt, q_term_sys );
    else         _calib[i]->momcor_mc  ( mu_vec_sys, _charge*1.0, _trk_layers, q_term_sys );

    // std::cout << "pt " << mu_vec_sys.Pt() << " returned by index " << i << std::endl;
    
    if ( mu_vec_sys.Pt() >= _mu_vec.Pt() ) {
      nUp   += 1;
      sum_sq_up   += pow( mu_vec_sys.Pt() - _mu_vec.Pt(), 2 );
    } else {
      nDown += 1;
      sum_sq_down += pow( mu_vec_sys.Pt() - _mu_vec.Pt(), 2 );
    }
    
    mu_vec_sys = mu_vec_orig;
    q_term_sys = q_term_orig;
  }

  _pt_sys_up   = ( _doSys ? _mu_vec.Pt() + sqrt(sum_sq_up   / nUp)   : -999 );
  _pt_sys_down = ( _doSys ? _mu_vec.Pt() + sqrt(sum_sq_down / nDown) : -999 ); 

  // std::cout << "\n*******   RESULT   *******" << std::endl;
  // std::cout << "mu_vec_orig.Pt() = " << mu_vec_orig.Pt() << ", q_term_orig = " << q_term_orig << std::endl;
  // std::cout << "_mu_vec.Pt() = " << _mu_vec.Pt() << ", _q_term = " << _q_term << std::endl;
  // std::cout << "nUp = " << nUp << ", _pt_sys_up = " << _pt_sys_up << std::endl;
  // std::cout << "nDown = " << nDown << ", _pt_sys_down = " << _pt_sys_down << std::endl;

}

