#!/bin/bash

nice scons -j4

#TRAININGTREES="true"
#TRAIN="true"
#OPTIONS=" -m 10000"
OPTIONS=" -m -1"

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

# Run with full MVA

#nice ./analyzer ggHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV125.root -r 8TeV $OPTIONS  
nice ./analyzer vbfHmumu125_8TeV.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV125.root -r 8TeV $OPTIONS 

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
