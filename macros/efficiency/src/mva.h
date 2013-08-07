#ifndef hmumu_mva_h
#define hmumu_mva_h

#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>

#include "DataFormats.h"
#include "helpers.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace std;

class MVA
{
  public:


  MVA(const std::vector<std::string> configFileNames, const std::string outFileName);
  ~MVA();
  float getMVA(const std::string configFileName, const std::string mvaName);
  bool getMVAPassBDTCut(const std::string configFileName);
  void resetValues();
  float getSigEffCut(const std::string configFileName, float eff);

  /////////////////////

  float mDiMu;
  float ptDiMu;
  float yDiMu;
  float mDiJet;
  float ptDiJet;
  float yDiJet;
  float ptMu1;
  float ptMu2;
  float etaMu1;
  float etaMu2;
  float ptJet1;
  float ptJet2;
  float etaJet1;
  float etaJet2;
  float cosThetaStar;
  float cosThetaStarCS;
  float deltaEtaJets;
  float productEtaJets;
  float nJetsInRapidityGap;

  float deltaEtaMuons;
  float deltaPhiMuons;
  float deltaRMuons;
  float deltaPhiJets;
  float deltaRJets;
  float deltaPhiHJ1;

  float relIsoMu1;
  float relIsoMu2;
  float ht;
  float nJets;
  float htInRapidityGap;

  float nVtx;

  float weight;

  float met;
  float ptmiss;
  float nPU;

  int vbfPreselection;

  float mDiMuResSigUp;
  float mDiMuResSigDown;
  float mDiMuResASigUp;
  float mDiMuResASigDown;

  float puJetIDSimpleDiscJet1;
  float puJetIDSimpleDiscJet2;
  float puJetIDSimpleDiscJet3;

  int puJetIDSimpleJet1;
  int puJetIDSimpleJet2;
  int puJetIDSimpleJet3;

  float bdtValInc;
  float bdtValVBF;

  //////////////////////
  
  std::map<std::string,TMVA::Reader*> readers_;
  std::map<std::string,float> mvaCuts_;
  TFile* outFile_;
  TTree* outTree_;

  void writeEvent(){if(outTree_ != NULL) outTree_->Fill();};

};

#endif
