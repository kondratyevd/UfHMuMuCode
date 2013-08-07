#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<THStack.h>
#include<TGraphErrors.h>
#include<TGraphAsymmErrors.h>
#include<TCanvas.h>
#include<TFrame.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
using namespace std;

#include "/afs/cern.ch/user/d/digiovan/scripts/Init/modifiedStyle.C"

#include "scripts/htomm/includes/selections.h"

// DATASET
//    0: DYJetsToLL
//    1: DYToMuMu
//    2: TTJets
//    3: DYToTauTau
//    4: WJetsToLNu 
//    5: WW
//    6: WZ
//    7: ZZ
//    8: QCD
//    9: DY2JetsToLL

// 9990: DYJetsToLL ptZ_50to70
// 9991: DYJetsToLL ptZ_70to100
// 9992: DYJetsToLL ptZ>100

//   10: 2012A-13Jul2012-v1
//   11: 2012A-recover-06Aug2012-v1
//   12: 2012B-13Jul2012-v1
//   13: 2012C-24Aug2012-v1
//   14: 2012C-PromptReco-v2
//   15: 2012D-PromptReco-v1


//   20: 2012A-13Jul2012-v1
//   21: 2012B-13Jul2012-v4
//   22: 2012C-24Aug2012-v1
//   23: 2012C-PromptReco-v2

// WHICH SKIMS
//    0: no cut
//    1: minimal cuts:
//     * M_mm in [60,120] GeV, 
//     * angle_dimuon < pi - 0.02
//     * dimuon vertex chisquare cut
//    2: 2 tight muons + minimal cuts, NO acceptance/trigger cuts

// default is Drell-Yan MC with minimal cuts
void skims(int dataset = 0,
           int whichSkims = 1) {

  if      ( whichSkims == 0 ) std::cout << " \n*** MERGING ALL THE FILES ***\n\n";
  else if ( whichSkims == 1 ) std::cout << " \n*** APPLYING MINIMAL CUTS ***\n\n";
  else if ( whichSkims == 2 ) std::cout << " \n*** APPLYING TIGHT CUTS, BUT ACCEPTANCE AND TRIGGER SELECTIONS ***\n\n";
  else {
    std::cout << "NO SKIMS TYPE DEFINED!\n";
    return;
  }


  TChain* tree = new TChain("tree");

  // ======== Monte Carlo ==========
  if ( dataset == 0 ) {
    std::cout << "DYJetsToLL_M50 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDYJetsToLL
    #include "scripts/htomm/skims/V00-01-10/chainDYJetsToLL"
  }
  
  else if ( dataset == 1 ) {
    std::cout << "DYToMuMu DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDYToMuMu
    #include "scripts/htomm/skims/V00-01-10/chainDYToMuMu"
  }
  
  else if ( dataset == 2 ) {
    std::cout << "TTJets DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainTTJets
    #include "scripts/htomm/skims/V00-01-10/chainTTJets"
  }

  else if ( dataset == 3 ) {
    std::cout << "DYToTauTau DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDYToTauTau
    #include "scripts/htomm/skims/V00-01-10/chainDYToTauTau"
  }

  else if ( dataset == 4 ) {
    std::cout << "WJetsToLNu DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainWJetsToLNu
    #include "scripts/htomm/skims/V00-01-10/chainWJetsToLNu"
  }

  else if ( dataset == 5 ) {
    std::cout << "WW DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainWW
    #include "scripts/htomm/skims/V00-01-10/chainWW"
  }

  else if ( dataset == 6 ) {
    std::cout << "WZ DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainWZ
    #include "scripts/htomm/skims/V00-01-10/chainWZ"
  }

  else if ( dataset == 7 ) {
    std::cout << "ZZ DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainZZ
    #include "scripts/htomm/skims/V00-01-10/chainZZ"
  }

  else if ( dataset == 8 ) {
    std::cout << "QCD DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCQCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCQCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainQCD_Pt_20_MuEnrichedPt_15
    #include "scripts/htomm/skims/V00-01-10/chainQCD_Pt_20_MuEnrichedPt_15"
  }
  else if ( dataset == 9 ) {
    std::cout << "DY2JetsToLL_M50 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDY2JetsToLL
    #include "scripts/htomm/skims/V00-01-10/chainDY2JetsToLL"
  }


  // ======== Special Monte Carlo =========
  else if ( dataset == 9990 ) {
    std::cout << "DYJetsToLL_PtZ-50To70 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDYJetsToLL_PtZ-50To70
    #include "scripts/htomm/skims/V00-01-10/chainDYJetsToLL_PtZ-50To70"
  }

  else if ( dataset == 9991 ) {
    std::cout << "DYJetsToLL_PtZ-70To100 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDYJetsToLL_PtZ-70To100
    #include "scripts/htomm/skims/V00-01-10/chainDYJetsToLL_PtZ-70To100"
  }

  else if ( dataset == 9992 ) {
    std::cout << "DYJetsToLL_PtZ-100 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCDYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDYJetsToLL_PtZ-100
    #include "scripts/htomm/skims/V00-01-10/chainDYJetsToLL_PtZ-100"
  }



  // ======== Data ==========
  // Single Muon Dataset
  else if ( dataset == 10 ) {
    std::cout << "SingleMu 2012A-13Jul2012-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012A-13Jul2012-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012A-13Jul2012-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainSingleMuRun2012A-13Jul2012-v1
    #include "scripts/htomm/skims/V00-01-10/chainSingleMuRun2012A-13Jul2012-v1"
  }

  else if ( dataset == 11 ) {
    std::cout << "SingleMu 2012A-recover-06Aug2012-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012A-recover-06Aug2012-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012A-recover-06Aug2012-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainSingleMuRun2012A-recover-06Aug2012-v1
    #include "scripts/htomm/skims/V00-01-10/chainSingleMuRun2012A-recover-06Aug2012-v1"
  }

  else if ( dataset == 12 ) {
    std::cout << "SingleMu 2012B-13Jul2012-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012B-13Jul2012-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012B-13Jul2012-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainSingleMuRun2012B-13Jul2012-v1
    #include "scripts/htomm/skims/V00-01-10/chainSingleMuRun2012B-13Jul2012-v1"
  }

  else if ( dataset == 13 ) {
    std::cout << "SingleMu 2012C-24Aug2012-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012C-24Aug2012-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012C-24Aug2012-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainSingleMuRun2012C-24Aug2012-v1
    #include "scripts/htomm/skims/V00-01-10/chainSingleMuRun2012C-24Aug2012-v1"
  }

  else if ( dataset == 14 ) {
    std::cout << "SingleMu 2012C-PromptReco-v2 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012C-PromptReco-v2/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012C-PromptReco-v2/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainSingleMuRun2012C-PromptReco-v2
    #include "scripts/htomm/skims/V00-01-10/chainSingleMuRun2012C-PromptReco-v2"
  }

  else if ( dataset == 15 ) {
    std::cout << "SingleMu 2012D-PromptReco-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012D-PromptReco-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataSingleMuRun2012D-PromptReco-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainSingleMuRun2012D-PromptReco-v1
    #include "scripts/htomm/skims/V00-01-10/chainSingleMuRun2012D-PromptReco-v1"
  }


  // Double Muon Dataset
  else if ( dataset == 20 ) {
    std::cout << "DoubleMu 2012A-13Jul2012-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012A-13Jul2012-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012A-13Jul2012-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012A-13Jul2012-v1
    #include "scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012A-13Jul2012-v1"
  }

  else if ( dataset == 21 ) {
    std::cout << "DoubleMu 2012B-13Jul2012-v4 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012B-13Jul2012-v4/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012B-13Jul2012-v4/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012B-13Jul2012-v4
    #include "scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012B-13Jul2012-v4"
  }

  else if ( dataset == 22 ) {
    std::cout << "DoubleMu 2012C-24Aug2012-v1 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012C-24Aug2012-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012C-24Aug2012-v1/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012C-24Aug2012-v1
    #include "scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012C-24Aug2012-v1"
  }

  else if ( dataset == 23 ) {
    std::cout << "DoubleMu 2012C-PromptReco-v2 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012C-PromptReco-v2/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesDataDoubleMuRun2012C-PromptReco-v2/"$1"\");"}' >> scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012C-PromptReco-v2
    #include "scripts/htomm/skims/V00-01-10/chainDoubleMuRun2012C-PromptReco-v2"
  }

  else {
    std::cout << "NO DATASET DEFINED!\n";
    return;
  }


  // ======== DEFINE THE FILENAME ========
  TString path = "/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/";

  TFile *newfile;

  if      ( dataset == 0 ) path+= "NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 1 ) path+= "NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 2 ) path+= "NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/";
  else if ( dataset == 3 ) path+= "NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 4 ) path+= "NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/";
  else if ( dataset == 5 ) path+= "NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 6 ) path+= "NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 7 ) path+= "NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 8 ) path+= "NtuplesMCQCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/";
  else if ( dataset == 9 ) path+= "NtuplesMCDY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/";

  else if ( dataset == 9990 ) path+="NtuplesMCDYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/";
  else if ( dataset == 9990 ) path+="NtuplesMCDYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/";
  else if ( dataset == 9992 ) path+="NtuplesMCDYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/";

  else if ( dataset == 10 ) path+="NtuplesDataSingleMuRun2012A-13Jul2012-v1/";
  else if ( dataset == 11 ) path+="NtuplesDataSingleMuRun2012A-recover-06Aug2012-v1/";
  else if ( dataset == 12 ) path+="NtuplesDataSingleMuRun2012B-13Jul2012-v1/";
  else if ( dataset == 13 ) path+="NtuplesDataSingleMuRun2012C-24Aug2012-v1/";
  else if ( dataset == 14 ) path+="NtuplesDataSingleMuRun2012C-PromptReco-v2/";
  else if ( dataset == 15 ) path+="NtuplesDataSingleMuRun2012D-PromptReco-v1/";
                                 
  else if ( dataset == 20 ) path+="NtuplesDataDoubleMuRun2012A-13Jul2012-v1/";
  else if ( dataset == 21 ) path+="NtuplesDataDoubleMuRun2012B-13Jul2012-v4/";
  else if ( dataset == 22 ) path+="NtuplesDataDoubleMuRun2012C-24Aug2012-v1/";
  else if ( dataset == 23 ) path+="NtuplesDataDoubleMuRun2012C-PromptReco-v2/";

  else {
    std::cout << "NO DATASET DEFINED!\n";
    return;
  }



  if      (whichSkims == 0) path += "merged/";
  else if (whichSkims == 1) path += "minimal/";
  else if (whichSkims == 2) path += "complete/";
  else {
    std::cout << "NO SKIMS TYPE DEFINED!\n";
    return;
  }



  if      ( dataset == 0 ) path += "DYJetsToLL_";
  else if ( dataset == 1 ) path += "DYToMuMu_";
  else if ( dataset == 2 ) path += "TTJets_";
  else if ( dataset == 3 ) path += "DYToTauTau_";
  else if ( dataset == 4 ) path += "WJetsToLNu_";
  else if ( dataset == 5 ) path += "WW_";
  else if ( dataset == 6 ) path += "WZ_";
  else if ( dataset == 7 ) path += "ZZ_";
  else if ( dataset == 8 ) path += "QCD_Pt_20_MuEnrichedPt_15_";
  else if ( dataset == 9 ) path += "DY2JetsToLL_";


  else if ( dataset == 9990 ) path += "DYJetsToLL_PtZ-50To70_";
  else if ( dataset == 9990 ) path += "DYJetsToLL_PtZ-70To100_";
  else if ( dataset == 9992 ) path += "DYJetsToLL_PtZ-100_";


  else if ( dataset == 10 ) path += "SingleMuRun2012A-13Jul2012-v1_";
  else if ( dataset == 11 ) path += "SingleMuRun2012A-recover-06Aug2012-v1_";
  else if ( dataset == 12 ) path += "SingleMuRun2012B-13Jul2012-v1_";
  else if ( dataset == 13 ) path += "SingleMuRun2012C-24Aug2012-v1_";
  else if ( dataset == 14 ) path += "SingleMuRun2012C-PromptReco-v2_";
  else if ( dataset == 15 ) path += "SingleMuRun2012D-PromptReco-v1_";


  else if ( dataset == 20 ) path += "DoubleMuRun2012A-13Jul2012-v1_";
  else if ( dataset == 21 ) path += "DoubleMuRun2012B-13Jul2012-v4_";
  else if ( dataset == 22 ) path += "DoubleMuRun2012C-24Aug2012-v1_";
  else if ( dataset == 23 ) path += "DoubleMuRun2012C-PromptReco-v2_";

  else {
    std::cout << "NO DATASET DEFINED!\n";
    return;
  }


  if      (whichSkims == 0) newfile = new TFile(path+"merged.root","recreate");
  else if (whichSkims == 1) newfile = new TFile(path+"minimal.root","recreate");
  else if (whichSkims == 2) newfile = new TFile(path+"complete.root","recreate");
  else {
    std::cout << "NO SKIMS TYPE DEFINED!\n";
    return;
  }


  std::cout << "saving in file " << newfile -> GetName() << std::endl;
  TTree *newtree = tree->CloneTree(0);


  // ======== GET THE HANDLES FOR SKIMMING ========
  float recoCandMass;
  tree->SetBranchAddress("recoCandMass", &recoCandMass);

  float vertexNormChiSquare, angleDiMuons;

  tree->SetBranchAddress("vertexNormChiSquare", &vertexNormChiSquare);
  tree->SetBranchAddress("angleDiMuons"       , &angleDiMuons );
  
  _MuonInfo reco1, reco2;
  
  tree->SetBranchAddress("reco1", &reco1);
  tree->SetBranchAddress("reco2", &reco2);


  // ======== PERFORM THE SKIMMING ========
  cout<<"Loop over the " << tree->GetEntries() << " entries ...\n";
  for (int iEvt=0; iEvt < tree->GetEntries(); iEvt++) {
    
    if ( (iEvt % 500000)==0 ) cout << "event " << iEvt << endl;
    tree -> GetEntry(iEvt);

    // additional selection cuts
    if (recoCandMass <  60 && whichSkims > 0) continue;
    //if (recoCandMass > 160 && whichSkims > 0) continue;
    if (vertexNormChiSquare > 10 && whichSkims > 0) continue; 
    if (angleDiMuons > TMath::Pi()-0.02 && whichSkims > 0) continue;

    if (!isKinTight_2012_noAcc(reco1) && whichSkims > 1) continue;
    if (!isKinTight_2012_noAcc(reco2) && whichSkims > 1) continue;

    newtree->Fill(); 
  }

  newtree->Print();
  newtree->AutoSave();
  
  std::cout << "new tree has " << newtree -> GetEntries() << std::endl;
  delete newfile;
}


