#ifndef CompositeCandMassResolution_h
#define CompositeCandMassResolution_h

namespace reco { class Candidate; class Muon; class GsfElectron; class Track; class PFCandidate; }
namespace edm { class EventSetup; }

#include <vector>
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include <TMatrixDSym.h>
#include "TLorentzVector.h"

class CompositeCandMassResolution  {
    public:
        CompositeCandMassResolution() {}
        ~CompositeCandMassResolution() {}
        void init(const edm::EventSetup &iSetup);
        double getMassResolution(const std::vector<reco::Track> &leaves,
                                 const TLorentzVector &c) const ;
    private:
        void   fillP3Covariance(const reco::Track &t, TMatrixDSym &bigCov, int offset) const ;
        edm::ESHandle<MagneticField> magfield_;
};

#endif

