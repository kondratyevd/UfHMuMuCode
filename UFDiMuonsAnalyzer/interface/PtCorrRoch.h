
#ifndef PT_CORR_ROCH
#define PT_CORR_ROCH

#include "RochCor/Calibration/interface/rochcor2016.h"
#include "RochCor/Calibration/interface/RoccoR.h"
#include <math.h>
#include "TLorentzVector.h"

void CorrectPtRoch( rochcor2016* _calib[201], const bool _doSys,
		    TLorentzVector& _mu_vec, float& _q_term, 
		    float& _pt_sys_up, float& _pt_sys_down,
                    const int _charge, const int _trk_layers, const bool _isData );

#endif  // #ifndef PT_CORR_ROCH
