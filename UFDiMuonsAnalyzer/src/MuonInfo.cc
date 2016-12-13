
#import "UfHMuMuCode/UFDiMuonsAnalyzer/interface/MuonInfo.h"

void MuonInfo::init() {

    isTracker    = -999;
    isStandAlone = -999;
    isGlobal     = -999;
    
    isTightID    = -999;
    isMediumID   = -999;
    isLooseID    = -999;
    
    charge = -999;
    pt     = -999;
    eta    = -999;
    phi    = -999;

    trkPt    = -999;
    trkPtErr = -999;
    trkEta   = -999;
    trkPhi   = -999;
    
    pt_KaMu           = -999;
    ptErr_KaMu        = -999;
    pt_sys_up_KaMu    = -999;
    pt_sys_down_KaMu  = -999;
    pt_clos_up_KaMu   = -999;
    pt_clos_down_KaMu = -999;

    pt_Roch          = -999;
    q_term_Roch      = -999;
    pt_sys_up_Roch   = -999;
    pt_sys_down_Roch = -999;
     
    d0_BS = -999;
    dz_BS = -999;
    
    d0_PV = -999;
    dz_PV = -999;

    relIso = -999;
    
    trackIsoSumPt     = -999;
    trackIsoSumPtCorr = -999;
    hcalIso           = -999;
    ecalIso           = -999;
    relCombIso        = -999;
    
    isPF = -999;
    
    pfPt  = -999;
    pfEta = -999;
    pfPhi = -999;
    
    sumChargedHadronPtR03   = -999;
    sumChargedParticlePtR03 = -999;
    sumNeutralHadronEtR03   = -999;
    sumPhotonEtR03          = -999;
    sumPUPtR03              = -999;
    
    sumChargedHadronPtR04   = -999;
    sumChargedParticlePtR04 = -999;
    sumNeutralHadronEtR04   = -999;
    sumPhotonEtR04          = -999;
    sumPUPtR04              = -999;
    
    // for (unsigned int iTrigger = 0; iTrigger < triggerArraySize; iTrigger++) {
    //   isHltMatched[iTrigger] = -999;
    //   hltPt[iTrigger]        = -999;
    //   hltEta[iTrigger]       = -999;
    //   hltPhi[iTrigger]       = -999;
    // }

} // End void MuonInfo::init()

