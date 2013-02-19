#!/bin/bash

nice scons -j4

TRAININGTREES="true"
TRAIN="true"
#OPTIONS=" -m 1000"

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-10/

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

if [ "$TRAININGTREES" = "true" ]; then
echo "#######################"
echo "Creating Training Trees"
echo "#######################"

nice ./analyzer DYJetsToLL_7TeV.root $DIR/NtuplesMCDYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYJetsToLL_minimal.root --trainingTree backgroundTreeDY_7TeV.root -r 7TeV $OPTIONS >& log7TeV2 &
#nice ./analyzer DYToMuMu_7TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1/minimal/DYToMuMu_minimal.root --trainingTree backgroundTreeDY_7TeV.root -r 7TeV $OPTIONS  >& log7TeV2 &
nice ./analyzer ttbar_7TeV.root $DIR/NtuplesMCTTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/TTJets_minimal.root --trainingTree backgroundTreeTT_7TeV.root -r 7TeV $OPTIONS  >& log7TeV2 &

nice ./analyzer ggHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu7TeV125.root --trainingTree signalTreeGG_7TeV.root -r 7TeV $OPTIONS 
nice ./analyzer vbfHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu7TeV125.root --trainingTree signalTreeVBF_7TeV.root -r 7TeV $OPTIONS 
nice ./analyzer zHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/zHmumu7TeV125.root --trainingTree signalTreeZH_7TeV.root -r 7TeV $OPTIONS 
nice ./analyzer wHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/wHmumu7TeV125.root --trainingTree signalTreeWH_7TeV.root -r 7TeV $OPTIONS 

nice ./analyzer WW_7TeV.root $DIR/NtuplesMCWW_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WW_minimal.root --trainingTree backgroundTreeWW_7TeV.root -r 7TeV $OPTIONS 
nice ./analyzer WZ_7TeV.root $DIR/NtuplesMCWZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WZ_minimal.root --trainingTree backgroundTreeWZ_7TeV.root -r 7TeV $OPTIONS 
nice ./analyzer ZZ_7TeV.root $DIR/NtuplesMCZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/ZZ_minimal.root --trainingTree backgroundTreeZZ_7TeV.root -r 7TeV $OPTIONS 

#nice ./analyzer DYToTauTau_7TeV.root $DIR/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYToTauTau_minimal.root --trainingTree backgroundTreeDYToTauTau_7TeV.root -r 7TeV $OPTIONS

wait

echo "#######################"
echo "#######################"
fi

if [ "$TRAIN" = "true" ]; then
echo "#######################"
echo "    Training MVAs"
echo "#######################"
echo "#######################" >& log7TeV2
echo "    Training MVAs" >& log7TeV2
echo "#######################" >& log7TeV2
echo "training Inclusive..."
nice ./mvaTrain inclusive_7TeV.cfg >& logMVAInc7TeV
echo "training VBF..."
nice ./mvaTrain vbf_7TeV.cfg >& logMVAVBF7TeV
echo "done training."
echo "#######################"
echo "#######################"
fi

# Run with full MVA
nice ./analyzer DYJetsToLL_7TeV.root $DIR/NtuplesMCDYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYJetsToLL_minimal.root -r 7TeV $OPTIONS  >& log7TeV2 &
#nice ./analyzer DYToMuMu_7TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1/minimal/DYToMuMu_minimal.root -r 7TeV $OPTIONS  >& log7TeV2 &
nice ./analyzer ttbar_7TeV.root $DIR/NtuplesMCTTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/TTJets_minimal.root -r 7TeV $OPTIONS  >& log7TeV2 &

nice ./analyzer ggHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu7TeV125.root -r 7TeV $OPTIONS 
nice ./analyzer vbfHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu7TeV125.root -r 7TeV $OPTIONS 
nice ./analyzer zHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/zHmumu7TeV125.root -r 7TeV $OPTIONS 
nice ./analyzer wHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/wHmumu7TeV125.root -r 7TeV $OPTIONS 

nice ./analyzer DYToTauTau_7TeV.root $DIR/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYToTauTau_minimal.root -r 7TeV $OPTIONS 
nice ./analyzer WW_7TeV.root $DIR/NtuplesMCWW_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WW_minimal.root -r 7TeV $OPTIONS 
nice ./analyzer WZ_7TeV.root $DIR/NtuplesMCWZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WZ_minimal.root -r 7TeV $OPTIONS  
nice ./analyzer ZZ_7TeV.root $DIR/NtuplesMCZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/ZZ_minimal.root -r 7TeV $OPTIONS  
#nice ./analyzer WJetsToLNu_7TeV.root $DIR/NtuplesMCWJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WJetsToLNu_minimal.root -r 7TeV $OPTIONS >& log7TeV2 &
#nice ./analyzer QCD_7TeV.root $DIR/NtuplesMCQCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_Fall11-PU_S6_START44_V9B-v1/minimal/QCD_Pt_20_MuEnrichedPt_15_minimal.root -r 7TeV $OPTIONS >& log7TeV2 &

nice ./analyzer SingleMuRun2011Av1.root $DIR/NtuplesDataSingleMuRun2011A-08Nov2011-v1/minimal/SingleMuRun2011A-08Nov2011-v1_minimal.root -r 7TeV $OPTIONS >& log7TeV2 &
nice ./analyzer SingleMuRun2011Bv1.root $DIR/NtuplesDataSingleMuRun2011B-19Nov2011-v1/minimal/SingleMuRun2011B-19Nov2011-v1_minimal.root -r 7TeV $OPTIONS 

wait

tar czf result.tgz ggHmumu*.root vbfHmumu*.root zHmumu*.root wHmumu*.root ttbar*.root DY*.root WW*.root WZ*.root ZZ*.root SingleMu*.root

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
