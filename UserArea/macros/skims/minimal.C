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

#include "/data/0b/digiovan/code/higgs/CMSSW_5_3_3_patch3/src/UserArea/UFDiMuonsAnalyzer/interface/DataFormats.h"


void minimal(int dataset = 0) {

  //dataset
  //    0: DYJetsToLL
  //    1: DYToMuMu
  //    2: TTJets
  //    3: DYToTauTau
  //    4: WJetsToLNu 
  //    5: WW
  //    6: WZ
  //    7: ZZ
  //    8: QCD
          

  TChain* tree = new TChain("tree");

  // ======== Monte Carlo ==========
  if ( dataset == 0 ) {
    std::cout << "DYJetsToLL_M50 DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/chainDYJetsToLL
    #include "scripts/htomm/skims/chainDYJetsToLL"
  }
  
  else if ( dataset == 1 ) {
    std::cout << "DYToMuMu DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/chainDYToMuMu
    #include "scripts/htomm/skims/chainDYToMuMu"
  }
  
  else if ( dataset == 2 ) {
    std::cout << "TTJets DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/"$1"\");"}' >> scripts/htomm/skims/chainTTJets
    #include "scripts/htomm/skims/chainTTJets"
  }

  else if ( dataset == 3 ) {
    std::cout << "DYToTauTau DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/chainDYToTauTau
    #include "scripts/htomm/skims/chainDYToTauTau"
  }

  else if ( dataset == 4 ) {
    std::cout << "WJetsToLNu DATASET\n";
    #include "scripts/htomm/skims/chainWJetsToLNu"
  }

  else if ( dataset == 5 ) {
    std::cout << "WW DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/chainWW
    #include "scripts/htomm/skims/chainWW"
  }

  else if ( dataset == 6 ) {
    std::cout << "WZ DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/chainWZ
    #include "scripts/htomm/skims/chainWZ"
  }

  else if ( dataset == 7 ) {
    std::cout << "ZZ DATASET\n";
    //ls /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/ | grep root | awk '{print "tree->AddFile(\"/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/"$1"\");"}' >> scripts/htomm/skims/chainZZ
    #include "scripts/htomm/skims/chainZZ"
  }

  else if ( dataset == 8 ) {
    std::cout << "QCD DATASET: ";
    #include "scripts/htomm/skims/chainQCD_Pt-20_MuEnriched-Pt15"
  }


  else {
    std::cout << "NO DATASET DEFINED!\n";
    return;
  }


  // ======== DEFINE THE FILENAME ========
  TString path = "/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-00-03/";

  TFile *newfile;

  if ( dataset == 0 ) 
    newfile = new TFile(path+"NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYJetsToLL_minimal.root","recreate");

  else if ( dataset == 1 ) 
    newfile = new TFile(path+"NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root","recreate");

  else if ( dataset == 2 ) 
    newfile = new TFile(path+"NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3/minimal/TTJets_minimal.root","recreate");

  else if ( dataset == 3 ) 
    newfile = new TFile(path+"NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToTauTau_minimal.root","recreate");

  else if ( dataset == 4 ) 
    newfile = new TFile(path+"NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/minimal/WJetsToLNu_minimal.root","recreate");

  else if ( dataset == 5 ) 
    newfile = new TFile(path+"NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WW_minimal.root","recreate");

  else if ( dataset == 6 ) 
    newfile = new TFile(path+"NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WZ_minimal.root","recreate");

  else if ( dataset == 7 ) 
    newfile = new TFile(path+"NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/ZZ_minimal.root","recreate");

  else if ( dataset == 8 ) 
    newfile = new TFile(path+"/minimal/QCD_Pt-20_MuEnrichedPt-15_minimal.root","recreate");

  else {
    std::cout << "NO DATASET DEFINED!\n";
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
  

  // ======== PERFORM THE SKIMMING ========
  cout<<"Loop over the " << tree->GetEntries() << " entries ...\n";
  for (int iEvt=0; iEvt < tree->GetEntries(); iEvt++) {
    
    if ( (iEvt % 500000)==0 ) cout << "event " << iEvt << endl;
    tree -> GetEntry(iEvt);

    // additional selection cuts
    if (recoCandMass <  60) continue;
    if (recoCandMass > 160) continue;
    if (vertexNormChiSquare > 10) continue; 
    if (angleDiMuons > TMath::Pi()-0.02) continue;

    newtree->Fill(); 
  }

  newtree->Print();
  newtree->AutoSave();
  
  std::cout << "new tree has " << newtree -> GetEntries() << std::endl;
  delete newfile;
}


