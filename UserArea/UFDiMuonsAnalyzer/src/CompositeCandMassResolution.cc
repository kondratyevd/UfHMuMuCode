#include "UserArea/UFDiMuonsAnalyzer/interface/CompositeCandMassResolution.h"

#include <cmath>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

#include <TMatrixD.h>
#include "TLorentzVector.h"

void CompositeCandMassResolution::init(const edm::EventSetup &iSetup) {
    iSetup.get<IdealMagneticFieldRecord>().get(magfield_);
}

double CompositeCandMassResolution::getMassResolution(const std::vector<reco::Track> &leaves, 
                                                      const TLorentzVector &c) const  {
  
    double const MASS_MUON = 0.105658367;    //GeV/c2
    int n = leaves.size(), ndim = n*3;

    TMatrixDSym bigCov(ndim);
    TMatrixD jacobian(1,ndim);
    for (int i = 0, o = 0; i < n; ++i, o += 2) {
      const reco::Track &ci = leaves[i];
      fillP3Covariance(ci, bigCov, o);
      TLorentzVector tl;
      tl.SetPtEtaPhiM(ci.pt(), ci.eta(), ci.phi(), MASS_MUON);
 
      jacobian(0, o+0) = (c.Energy()*(tl.Px()/tl.Energy()) - c.Px())/c.M();
      jacobian(0, o+1) = (c.Energy()*(tl.Py()/tl.Energy()) - c.Py())/c.M();
      jacobian(0, o+2) = (c.Energy()*(tl.Pz()/tl.Energy()) - c.Pz())/c.M();
    }

    /*static int debug_ = 0;
      if (++debug_ < 20) {
      std::cout << "Big matrix:   " << std::endl; bigCov.Print();
      std::cout << "Jacobian:     " << std::endl; jacobian.Print();
      }*/
    TMatrixDSym massCov = bigCov.Similarity(jacobian);
    /*if (debug_ < 20) {
      std::cout << "Final matrix: " << std::endl; massCov.Print();
      }*/
    double dm2 = massCov(0,0);
    return (dm2 > 0 ? std::sqrt(dm2) : 0.0);
}

void CompositeCandMassResolution::fillP3Covariance(const reco::Track &t, TMatrixDSym &bigCov, int offset) const {
      GlobalTrajectoryParameters gp(GlobalPoint(t.vx(), t.vy(),  t.vz()),
                    GlobalVector(t.px(),t.py(),t.pz()),
                    t.charge(),
                    magfield_.product());
      JacobianCurvilinearToCartesian curv2cart(gp);
      CartesianTrajectoryError cartErr= ROOT::Math::Similarity(curv2cart.jacobian(), t.covariance());
      const AlgebraicSymMatrix66 mat = cartErr.matrix();
      for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) {
            bigCov(offset+i,offset+j) = mat(i+3,j+3);
      } } 
}

