#!/bin/bash

nice scons -j4

#TRAININGTREES="true"
#TRAIN="true"
#OPTIONS=" -m 1000"

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

if [ "$TRAININGTREES" = "true" ]; then
echo "#######################"
echo "Creating Training Trees"
echo "#######################"

nice ./analyzer DYJetsToLL_8TeV.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYJetsToLL_minimal.root --trainingTree backgroundTreeDY_8TeV.root -r 8TeV $OPTIONS >& log_DYJets < /dev/null &
#nice ./analyzer DYToMuMu_8TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root --trainingTree backgroundTreeDY_8TeV.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer ttbar_8TeV.root $DIR/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/minimal/TTJets_minimal.root --trainingTree backgroundTreeTT_8TeV.root -r 8TeV $OPTIONS  >& log_TTJets < /dev/null&

nice ./analyzer ggHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV125.root --trainingTree signalTreeGG_8TeV.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer vbfHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV125.root --trainingTree signalTreeVBF_8TeV.root -r 8TeV $OPTIONS 
nice ./analyzer zHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/zHmumu8TeV125.root --trainingTree signalTreeZH_8TeV.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer wHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/wHmumu8TeV125.root --trainingTree signalTreeWH_8TeV.root -r 8TeV $OPTIONS 

nice ./analyzer WW_8TeV.root $DIR/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WW_minimal.root --trainingTree backgroundTreeWW_8TeV.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer WZ_8TeV.root $DIR/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WZ_minimal.root --trainingTree backgroundTreeWZ_8TeV.root -r 8TeV $OPTIONS 
nice ./analyzer ZZ_8TeV.root $DIR/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/ZZ_minimal.root --trainingTree backgroundTreeZZ_8TeV.root -r 8TeV $OPTIONS 

#nice ./analyzer DYToTauTau_8TeV.root $DIR/NtuplesMCDYToTauTau_M_20_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToTauTau_minimal.root --trainingTree backgroundTreeDYToTauTau_8TeV.root -r 8TeV $OPTIONS

wait

echo "#######################"
echo "#######################"
fi

if [ "$TRAIN" = "true" ]; then
echo "#######################"
echo "    Training MVAs"
echo "#######################"
echo "#######################" >& log2
echo "    Training MVAs" >& log2
echo "#######################" >& log2
echo "training Inclusive..."
nice ./mvaTrain inclusive_8TeV.cfg >& logMVAInc_8TeV
echo "training VBF..."
nice ./mvaTrain vbf_8TeV.cfg >& logMVAVBF_8TeV
echo "done training."
echo "#######################"
echo "#######################"
fi

# Run with full MVA
nice ./analyzer DYJetsToLL_8TeV.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYJetsToLL_minimal.root -r 8TeV $OPTIONS  >& log2 &
#nice ./analyzer DYToMuMu_8TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer ttbar_8TeV.root $DIR/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/minimal/TTJets_minimal.root -r 8TeV $OPTIONS  >& log2 &

nice ./analyzer ggHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV125.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer vbfHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV125.root -r 8TeV $OPTIONS 
nice ./analyzer zHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/zHmumu8TeV125.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer wHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/wHmumu8TeV125.root -r 8TeV $OPTIONS 

nice ./analyzer ggHmumu123_8TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV123.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer vbfHmumu123_8TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV123.root -r 8TeV $OPTIONS 
nice ./analyzer ggHmumu127_8TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV127.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer vbfHmumu127_8TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV127.root -r 8TeV $OPTIONS 

#nice ./analyzer DYToTauTau_8TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer WW_8TeV.root $DIR/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WW_minimal.root -r 8TeV $OPTIONS 
nice ./analyzer WZ_8TeV.root $DIR/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WZ_minimal.root -r 8TeV $OPTIONS   >& log2 &
nice ./analyzer ZZ_8TeV.root $DIR/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/ZZ_minimal.root -r 8TeV $OPTIONS  
#nice ./analyzer WJetsToLNu_8TeV.root $DIR/NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/*.root -r 8TeV $OPTIONS >& log2 &
#nice ./analyzer QCD_8TeV.root $DIR/NtuplesMCQCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/minimal/QCD_Pt_20_MuEnrichedPt_15_minimal.root -r 8TeV $OPTIONS >& log2 &

nice ./analyzer SingleMuRun2012Av1.root $DIR/NtuplesDataSingleMuRun2012A-13Jul2012-v1/minimal/SingleMuRun2012A-13Jul2012-v1_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer SingleMuRun2012Av1Recover.root $DIR/NtuplesDataSingleMuRun2012A-recover-06Aug2012-v1/minimal/SingleMuRun2012A-recover-06Aug2012-v1_minimal.root -r 8TeV $OPTIONS 
nice ./analyzer SingleMuRun2012Bv1.root $DIR/NtuplesDataSingleMuRun2012B-13Jul2012-v1/minimal/SingleMuRun2012B-13Jul2012-v1_minimal.root -r 8TeV $OPTIONS 
nice ./analyzer SingleMuRun2012Cv1.root $DIR/NtuplesDataSingleMuRun2012C-24Aug2012-v1/minimal/SingleMuRun2012C-24Aug2012-v1_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer SingleMuRun2012Cv2.root $DIR/NtuplesDataSingleMuRun2012C-PromptReco-v2/minimal/SingleMuRun2012C-PromptReco-v2_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer SingleMuRun2012D.root $DIR/NtuplesDataSingleMuRun2012D-PromptReco-v1/minimal/SingleMuRun2012D-PromptReco-v1_minimal.root -r 8TeV $OPTIONS 

wait

tar czf result.tgz ggHmumu*.root vbfHmumu*.root zHmumu*.root wHmumu*.root ttbar*.root DY*.root WW*.root WZ*.root ZZ*.root SingleMu*.root

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
