#!/bin/bash

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/

#MAXEVENTS="-m 10000"

#Muon Selections
CUT="reco1.isGlobal && reco1.isPFMuon && reco1.pt>25.0 && abs(reco1.eta)<2.1 && reco1.numTrackerLayers>5 && abs(reco1.d0_PV)<0.2 && abs(reco1.dz_PV)<0.5 && reco1.numValidMuonHits>0 && reco1.numValidPixelHits>0 && reco1.numOfMatchedStations>1 && reco1.normChiSquare<10.0 && (reco1.sumChargedHadronPtR04+max(0.0,reco1.sumNeutralHadronEtR04+ reco1.sumPhotonEtR04 - 0.5*reco1.sumPUPtR04))/reco1.pt<0.12 && reco2.isGlobal && reco2.isPFMuon && reco2.pt>25.0 && abs(reco2.eta)<2.1 && reco2.numTrackerLayers>5 && abs(reco2.d0_PV)<0.2 && abs(reco2.dz_PV)<0.5 && reco2.numValidMuonHits>0 && reco2.numValidPixelHits>0 && reco2.numOfMatchedStations>1 && reco2.normChiSquare<10.0 && (reco2.sumChargedHadronPtR04+max(0.0,reco2.sumNeutralHadronEtR04+ reco2.sumPhotonEtR04 - 0.5*reco2.sumPUPtR04))/reco2.pt<0.12 && reco1.charge != reco2.charge"

#Just Iso
#CUT="(reco2.sumChargedHadronPtR04+max(0.0,reco2.sumNeutralHadronEtR04+ reco2.sumPhotonEtR04 - 0.5*reco2.sumPUPtR04))/reco2.pt<0.12 && (reco1.sumChargedHadronPtR04+max(0.0,reco1.sumNeutralHadronEtR04+ reco1.sumPhotonEtR04 - 0.5*reco1.sumPUPtR04))/reco1.pt<0.12"
#Just Pt/Eta
#CUT="reco1.pt>25.0 && abs(reco1.eta)<2.1 && reco2.pt>25.0 && abs(reco2.eta)<2.1"

PUOPTIONS="--dataPUHist pileupDists/PileUpHist2012A.root --mcPUHist pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root"

nice ./skim NtupleDYJetsToLL_muSelSkim.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root -c "$CUT" $PUOPTIONS $MAXEVENTS
nice ./skim NtupleRun2012A_muSelSkim.root $DIR/NtuplesDataSingleMuRun2012A-13Jul2012-v1/*.root -c "$CUT" $MAXEVENTS
