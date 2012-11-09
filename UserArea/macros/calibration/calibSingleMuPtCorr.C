/*
 *  unfolding.C
 *  Created by Joseph Gartner on 12/2/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 */

#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
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
#include<TMatrixD.h>

#include "Math/VectorUtil_Cint.h"
#include "Math/GenVector/LorentzVector.h"
#include <TLorentzVector.h>

#include "TF1.h"
#include <TGraphErrors.h>
#include<modifiedStyle.C>
#include<userStyle.C>
//include Rochechester correction
#include<rochcor.h>
#include "MuScleFitCorrector.h"
// include smearing tool
#include <SmearingTool.h>
/*
// ROOFIT
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooPlot.h"

#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

using namespace RooFit ;
using namespace std;
*/

#include "./DataFormat.h"
//#include "./DataFormatBoostedZ.h"

  int HiggsNtuple = 1; // 1 - Higgs to MuMu ntuple is used wiht DataFormat.h <- check it and isKinTight_2012 (d0 and Z) 
                       // 0 - Boosted Z ntuoles are used with DataFormatBoostedZ.h <- check it and isKinTight_2012 (d0 and Z) 

  int MuCorr = 0; // 0 - no muon correction, 
                  // 1 - Rochester Correction,
                  // 2 - MuscleFit correcton 
  int Ismear = 1; // 0 - take pt reco from MC
                  // 1 - smear pt with own tools using pt gen post FSR  

  TString RunYear = "2012ABCsmall"; // 2012; 2011A; 2011B; 2012ABCsmall 
  //const float PTbin[] = {20., 30., 40., 50., 60., 70., 100., 150., 300.};
  //const float PTbin[] = {20., 30., 40., 45., 50., 60., 70., 100.};
  const float PTbin[] = {20., 30., 35., 40., 45., 50., 60., 70., 100.}; //default
  const float ETAbin[] = {0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
  const int NPThist = (sizeof(PTbin)/sizeof(float)-1);
  const int NETAhist = (sizeof(ETAbin)/sizeof(float)-1);
  const int Nhist = NPThist*NETAhist;

  double Pi  = acos(-1.);
  //const float PTbin_inv[] = {20., 30., 50., 100.};
  //const float ETAbin_inv[] = {0., 0.9, 2.1};
  //const float PTbin_inv[] = {0., 20., 30., 40., 50., 70., 100., 300.};
  //const float PTbin_inv[] = {25., 30., 35., 40., 50., 70., 100., 300.};
  const float PTbin_inv[] = {25., 30., 35., 40., 45., 50., 70., 150.};
  //const float ETAbin_inv[] = {0., 0.3, 0.8, 1.2, 1.6, 2.1};
  const float ETAbin_inv[] = {-2.1, -1.6, -1.2, -0.8, 0., 0.8, 1.2, 1.6, 2.1};
  //const float ETAbin_tag[] = {0., 0.8, 1.2, 2.1};
  const float ETAbin_tag[] = {0., 0.8, 1.2, 2.1};
  const int NPHIbin_inv = 10; // 10 bins from 0 to 2Pi 
  const int NPThist_inv = (sizeof(PTbin_inv)/sizeof(float)-1);
  const int NETAhist_inv = (sizeof(ETAbin_inv)/sizeof(float)-1);
  const int NETAhist_tag = (sizeof(ETAbin_tag)/sizeof(float)-1);
  const int Nhist_inv = NPThist_inv*NETAhist_inv*NETAhist_tag;
  const int Nhist_inv1 = NPThist_inv*NETAhist_inv;


  const float ScaleFactor = 1.;
  // variation Asig2 by 3 sigma;
  const float ScaleAsig2 = 0;// by default + 0 sigma of Asig2;
  //const float ScaleAsig2 = 3.;// + 3 sigma of Asig2;
  //const float ScaleAsig2 = -3.;// + 3 sigma of Asig2;
   //TString Extra = "";
   //TString Extra = "_AsigPlus3sigma";
   //TString Extra = "_AsigMinus3sigma";
   //TString Extra = "_MATGRAPH";

  Double_t MASS_MUON = 0.105658367;    //GeV/c2
  Double_t mean[Nhist];
  Double_t sig1[Nhist];
  Double_t sig2[Nhist];
  Double_t Asig2[Nhist]; 

//////////// functions:

  void fitDiMuon(TH1F* histoMC[], TH1F* histoDATA[], TString Extra);
// ---- CP error ----

TH1F* hdthetav1v2;
TH1F* hmuonRes[Nhist];
TH1F* hDimuonRes[Nhist];

TH1F* hDimuonMass[Nhist_inv];
TH1F* hDimuonMassGEN[Nhist_inv];
TH1F* hDimuonMassDATA[Nhist_inv];
TH1F* hDimuonMassTagMuPlus[Nhist_inv];
TH1F* hDimuonMassTagMuPlusGEN[Nhist_inv];
TH1F* hDimuonMassTagMuPlusDATA[Nhist_inv];

TH1F* hPhiDimuonMass[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassGEN[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassDATA[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassTagMuPlus[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassTagMuPlusGEN[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassTagMuPlusDATA[NPHIbin_inv*NETAhist_inv*NETAhist_tag];


TH1F* hdeltaR80[Nhist_inv];
TH1F* hdeltaR85_95[Nhist_inv];
TH1F* hdeltaR105[Nhist_inv];
TH1F* hdeltaRDATA80[Nhist_inv];
TH1F* hdeltaRDATA85_95[Nhist_inv];
TH1F* hdeltaRDATA105[Nhist_inv];

// test of problematic regions:
TH1F* hTestMC_ZPt; 
TH1F* hTestMC_ptTag; 
TH1F* hTestMC_ZMass; 
TH1F* hTestMC_ZMassGen; 
TH1F* hTestMC_DeltaR80; 
TH1F* hTestDeltaRMatch; 

// unfolding matrix

TH1F *raw1;
TH1F *test1;
TH1F *key1;
TH1F *key1_MC;
TH1F *ratTest;
TH1F *ratRaw;


TH2F *rMatrixVisualization;
TH2F *rInvMatrixVisualization;
TH2F *rPriMatrixVisualization;

TCanvas *_c1;
TCanvas *_c2;
TCanvas *_c3;
TCanvas *_c4;
TCanvas *_c5;
TLegend *_leg1;
TLegend *_leg2;
        float binning[19] = {0.,2.5,5,7.5,10,12.5,15,17.5,20,30,40,50,70,90,110,150,190,250,600};
        // to plot from 1. but fill from 0.:
        float binning_histo[19] = {1.,2.5,5,7.5,10,12.5,15,17.5,20,30,40,50,70,90,110,150,190,250,600};
        //float binning_histo[19];
        //binning_histo[0] = 1.;
        //for(int j=1; j<19; j++)
        //{
        //    binning_histo[j] = binning[j];
        //}

        float mapping[18][18];
        float primeMapping[18][18];
        float recoCounter[18];
        float realCounter[18];

using namespace std;

class calibSingleMuPtCorr
{
public :
   //constructor
   calibSingleMuPtCorr(){}
   void main();
};
void bookhistos();     // to book histograms
void printhistos();
Double_t BWnonrel(Double_t*, Double_t* );
Double_t Gauss(Double_t*, Double_t* );
Double_t FuncVoigtian(Double_t*, Double_t* );
Double_t FuncVoigtianBG(Double_t*, Double_t* );
Double_t DoubleGauss(Double_t*, Double_t* );
bool isKinTight_2012(_MuonInfo& muon, float pTcorr);
float getRelIso(_MuonInfo& muon, float pTcorr);

void DrawWithRatio(TCanvas *canvas, char *cTitle,
                  TH1F *hNum, TH1F *hDen);

void axis1F(TH1F  *histo,
           TAxis *xaxis,
           TAxis *yaxis,
           char  *xtitle,
           char  *ytitle);


void calibSingleMuPtCorr::main(){


	gROOT->Clear();
  	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	userStyle();
	//modifiedStyle();
  	// ---- open the MC files ----
  	TChain* treeMC = new TChain("tree");
        // min cuts: vertex, cosmics, mass
        if(HiggsNtuple == 1){
           if(RunYear == 2012 || RunYear == "2012ABCsmall") treeMC -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root");
           if(RunYear == "2011A" || RunYear == "2011B" || RunYear == "2011")treeMC -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-01/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1/minimal/DYToMuMu_minimal.root");
        }
        if(HiggsNtuple == 0)
        treeMC -> AddFile("/data/uftrig01b/digiovan/root/CMSSW_4_2_8_patch6/Ntuples/MC/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythiaFall11-PU_S6_START42_V14B-v1/minimal/boostedZ_DYToMuMu_minimal.root");
        // output file
        TFile *theFile    =new TFile(Form("calibSingleMu"+RunYear+"PtCorr%dSmearPt%d.root", MuCorr, Ismear), "RECREATE");

        bookhistos();// histo init 
        
        float rMass;
        float rPt;
        float rEta;
        _MuonInfo reco1, reco2;
        _TrackInfo true1mu, true2mu;
        _genPartInfo genZpostFSR;
        treeMC->SetBranchAddress("reco1",&reco1);
        treeMC->SetBranchAddress("reco2",&reco2);
        treeMC->SetBranchAddress("genZpostFSR",&genZpostFSR);
        treeMC->SetBranchAddress("genM1ZpostFSR",&true1mu);
        treeMC->SetBranchAddress("genM2ZpostFSR",&true2mu);
        treeMC->SetBranchAddress("recoCandMass",&rMass);
        treeMC->SetBranchAddress("recoCandPt",&rPt);
        treeMC->SetBranchAddress("recoCandEta",&rEta);
        cout << "Nevent to process = " << treeMC->GetEntries() << endl;

        //Rochester correction 
        //To get the central value of the momentum correction 
        rochcor *rmcor = new rochcor(); // make the pointer of rochcor class
        //REMARK : Need to call "rochcor(seed)" to assign the systematic error 
        //rochcor *rmcor = new rochcor(seed); //where "seed" is the random seed number

        ///
/*
        TString fitParFileMC, fitParFileData;
        if(RunYear == "2011A" || RunYear == "2011B" || RunYear == "2011"){
           fitParFileMC = "/data/uftrig01b/kropiv/HiggsMuMu/CMSSW_5_3_3_patch3/src/Calibration/MuScleFitCorrector_v1/MuScleFit_2011_MC_42X.txt";// no MC_44 :(
           fitParFileData = "/data/uftrig01b/kropiv/HiggsMuMu/CMSSW_5_3_3_patch3/src/Calibration/MuScleFitCorrector_v1/MuScleFit_2011_DATA_44X.txt"; 
        } 
        if(RunYear == "2012" || RunYear == "2012ABCsmall" || RunYear == "2012B" || RunYear == "2012C" || RunYear == "2012A" ){
           fitParFileMC = "/data/uftrig01b/kropiv/HiggsMuMu/CMSSW_5_3_3_patch3/src/Calibration/MuScleFitCorrector_v1/MuScleFit_2012_MC_52X.txt";// no MC_53 :(
           fitParFileData = "/data/uftrig01b/kropiv/HiggsMuMu/CMSSW_5_3_3_patch3/src/Calibration/MuScleFitCorrector_v1/MuScleFit_2012_DATA_53X.txt"; 
        }
        MuScleFitCorrector* correctorMC_ = new MuScleFitCorrector(fitParFileMC); 
        MuScleFitCorrector* correctorData_ = new MuScleFitCorrector(fitParFileData); 
*/
        // pointer to Smearing Tool:
        SmearingTool *smearPT = new SmearingTool();  
        int nbad_tpt = 0;
        int nbad_tMass = 0;
        //for( int k=0; k<100000; k++)
        for( int k=0; k<treeMC->GetEntries(); k++)
        {
                //process progress
                if(k!=0 && (k%10000)==0)
                std::cout << "- processing event " << k << "\r" << std::flush;

                treeMC->GetEntry(k);
                //cout << "any selection event = " << k << " reco1.charge = " << reco1.charge << " rMass = " << rMass << endl;

                if (reco1.charge == reco2.charge) continue; // opposit charged muons

                // check that dtheata between 2 muons less then Pi-0.02
                TLorentzVector MuReco1, MuReco2, MuTrue1, MuTrue2;
                float MuTrue1charge = true1mu.charge; 
                float MuTrue2charge = true2mu.charge; 
                MuReco1.SetPtEtaPhiM(reco1.pt, reco1.eta, reco1.phi, MASS_MUON);
                MuReco2.SetPtEtaPhiM(reco2.pt, reco2.eta, reco2.phi, MASS_MUON);
                float theta12 = MuReco1.Angle(MuReco2.Vect()); // angle between MuReco & MuReco
                hdthetav1v2 -> Fill (theta12);
                //if( theta12 >= (Pi-0.02) ) continue; 

                ////////////////////////
                //make matching between MuTrue1,2 and MuReco1,2:
                ROOT::Math::PtEtaPhiEVector MuTrue1h(true1mu.pt, true1mu.eta, true1mu.phi, true1mu.pt); 
                ROOT::Math::PtEtaPhiEVector MuTrue2h(true2mu.pt, true2mu.eta, true2mu.phi, true2mu.pt); 
                ROOT::Math::PtEtaPhiEVector MuReco1h(MuReco1.Pt(), reco1.eta, reco1.phi, MuReco1.Pt()); 
                ROOT::Math::PtEtaPhiEVector MuReco2h(MuReco2.Pt(), reco2.eta, reco2.phi, MuReco2.Pt()); 
                float deltaR_r1r2 = ROOT::Math::VectorUtil::DeltaR(MuReco1h, MuReco2h);  
                Float_t muonRes_t1r1 = -10.; 
                Float_t muonRes_t1r2 = -10.;
                Float_t muonRes_t2r1 = -10.; 
                Float_t muonRes_t2r2 = -10.;
                if(true1mu.pt > 0.) muonRes_t1r1 = (MuReco1.Pt()-true1mu.pt)/true1mu.pt;  
                if(true1mu.pt > 0.) muonRes_t1r2 = (MuReco2.Pt()-true1mu.pt)/true1mu.pt;  
                if(true2mu.pt > 0.) muonRes_t2r2 = (MuReco2.Pt()-true2mu.pt)/true2mu.pt;
                if(true2mu.pt > 0.) muonRes_t2r1 = (MuReco1.Pt()-true2mu.pt)/true2mu.pt;
                float deltaR_t1r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue1h, MuReco1h);  
                float deltaR_t1r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue1h, MuReco2h);  
                float deltaR_t2r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue2h, MuReco1h);  
                float deltaR_t2r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue2h, MuReco2h);  
             
                Float_t muonResCorr_t1r1 = muonRes_t1r1;
                Float_t muonResCorr_t2r2 = muonRes_t2r2;
                if(deltaR_t1r1 > deltaR_t1r2) muonResCorr_t1r1 = muonRes_t1r2;
                if(deltaR_t2r2 > deltaR_t2r1) muonResCorr_t2r2 = muonRes_t2r1;
                MuTrue1.SetPtEtaPhiM(true1mu.pt, true1mu.eta, true1mu.phi, MASS_MUON);
                MuTrue2.SetPtEtaPhiM(true2mu.pt, true2mu.eta, true2mu.phi, MASS_MUON);
                if(deltaR_t1r1 > deltaR_t1r2){
                     MuTrue1.SetPtEtaPhiM(true2mu.pt, true2mu.eta, true2mu.phi, MASS_MUON);
                     MuTrue1charge = true2mu.charge; 
                     deltaR_t1r1 = deltaR_t1r2;
                }
                if(deltaR_t2r2 > deltaR_t2r1){
                     MuTrue2.SetPtEtaPhiM(true1mu.pt, true1mu.eta, true1mu.phi, MASS_MUON);
                     MuTrue2charge = true1mu.charge; 
                     deltaR_t2r2 = deltaR_t2r1; 
                }
                ////////////////////////
                // Muon Correction
                //If you run MC, apply the muon momentum correction,"momcor_mc()" function (only for MC) 
                // This is for MC,
                //                 TLorentzVector, mu, corresponds to mu- for charge = -1 or mu+ for charge = +1 
                //                 sysdev == 0 returns the central value of muon momentum correction 
                //                 sysdev == 1 returns the muon momentum correction smeared 1 standard deviation
                //                             in <1/pt> correction 
                //                             (global factor is also changed by 1 sigma total error ( stat 0X syst)
                //                 If you want to assign the systematic error using x standard (statistical) deviation, 
                //                 you can assign sysdev == x.
                //                 runopt == 0 for 2011A correction 
                //                 runopt == 1 for 2011B correction
                //              TLor.  float   float   int             
                //rmcor->momcor_mc(mu, charge, sysdev, runopt); 
  
                if (MuCorr == 1){
                  if(RunYear == "2011B"){
                     rmcor->momcor_mc(MuReco1, float(reco1.charge), 0, 1); 
                     rmcor->momcor_mc(MuReco2, float(reco2.charge), 0, 1);
                  }
                  if(RunYear == "2011A"){
                     rmcor->momcor_mc(MuReco1, float(reco1.charge), 0, 0); 
                     rmcor->momcor_mc(MuReco2, float(reco2.charge), 0, 0);
                  }
                } 
/*                if (MuCorr == 2){
                   TLorentzVector* MuReco1Pointer = &MuReco1;//make pointer to MuReco1
                   TLorentzVector* MuReco2Pointer = &MuReco2;
                   if(reco1.charge < 0){
                        correctorMC_->applyPtCorrection(*MuReco1Pointer,-1);
                        //correctorMC_->applyPtSmearing(*MuReco1Pointer,-1);
                   } 
                   else{
                        correctorMC_->applyPtCorrection(*MuReco1Pointer,1);
                        //correctorMC_->applyPtSmearing(*MuReco1Pointer,1);
                   }
                   if(reco2.charge < 0){
                        correctorMC_->applyPtCorrection(*MuReco2Pointer,-1);
                        //correctorMC_->applyPtSmearing(*MuReco2Pointer,-1);
                   } 
                   else{
                        correctorMC_->applyPtCorrection(*MuReco2Pointer,1);
                        //correctorMC_->applyPtSmearing(*MuReco2Pointer,1);
                   }
                }
*/
                ////////////////////////
                /// smear PT after Correction if it took place (correction doesn't work any more, use correction function for smering):
                /// float PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen, TString ParVar = "null",float ParSig = 0);
                if (Ismear == 1){
                   float pt1Smear = smearPT -> PTsmear(MuTrue1.Pt(), MuTrue1.Eta(), MuTrue1charge);
                   float pt2Smear = smearPT -> PTsmear(MuTrue2.Pt(), MuTrue2.Eta(), MuTrue2charge);
                   // fill Smear values:
                   MuReco1.SetPtEtaPhiM(pt1Smear, reco1.eta, reco1.phi, MASS_MUON);
                   MuReco2.SetPtEtaPhiM(pt2Smear, reco2.eta, reco2.phi, MASS_MUON);
                   //TEST 
                   //cout << "ptgen1 = " << MuTrue1.Pt() << " ptSmear1 = " << pt1Smear << " ptreco1 = " << reco1.pt << endl; 
                   //cout << "ptgen2 = " << MuTrue2.Pt() << " ptSmear2 = " << pt2Smear << " ptreco2 = " << reco2.pt << endl; 
		   //cout << "=====================" << endl; 
                }
                ////////////////////////
                TLorentzVector MuRecoCand = MuReco1 + MuReco2; 
                float rMassCorr = MuRecoCand.M(); 
                if (rMassCorr <  60) continue;// mass restriction
                if (rMassCorr > 120) continue;

  
                //check that generator part is filled 
                if(nbad_tpt < 10 && (true1mu.pt < -900. || true2mu.pt < -900)){//rejection for MATGRAPH
                   cout << "CHECK ntuple, possible problem with fill gen level while RECO level is fine: true1mu.pt = " << true1mu.pt << " tMass = " << genZpostFSR.mass << endl; 
                   nbad_tpt++;
                }
                if(true1mu.pt < -900. || true2mu.pt < -900) continue;//rejection for MATGRAPH
                if(nbad_tMass < 10 && genZpostFSR.mass < -900.){//rejection for MATGRAPH
                   cout << "CHECK ntuple, possible EXTRA problem with fill gen level while RECO level is fine: tMass = " << genZpostFSR.mass << endl; 
                   nbad_tMass++;
                }

                float pTcorr1 = MuReco1.Pt(); 
                float pTcorr2 = MuReco2.Pt(); 
                if( !isKinTight_2012(reco1, pTcorr1) || !isKinTight_2012(reco2, pTcorr2) ) continue;


                  
                if(  (reco1.charge > 0 && MuReco1.Pt() >= 25 && MuReco1.Pt() < 30 && fabs(reco1.eta) < 0.3 && fabs(reco2.eta) < 1.)
                  || (reco2.charge > 0 && MuReco2.Pt() >= 25 && MuReco2.Pt() < 30 && fabs(reco2.eta) < 0.3 && fabs(reco1.eta) < 1.) ) {
                     hTestMC_ZMass -> Fill(rMassCorr);
                     hTestMC_ZMassGen -> Fill(rMassCorr);
                     if(rMassCorr < 80){
                       hTestMC_DeltaR80 -> Fill(deltaR_r1r2);
                       hTestDeltaRMatch -> Fill(deltaR_t1r1);
                       hTestDeltaRMatch -> Fill(deltaR_t2r2);
                       if(reco1.charge < 0) hTestMC_ptTag -> Fill(MuReco1.Pt());
                       if(reco2.charge < 0) hTestMC_ptTag -> Fill(MuReco2.Pt());
                     }
                     hTestMC_ZPt -> Fill(rPt);
                }
                //if(muonResCorr_t1r1 < -0.5 && true1mu.pt >= 70.) cout << "muonResCorr_t1r1 = " << muonResCorr_t1r1 << " true pt = " << true1mu.pt << " reco pt = " << MuReco1.Pt() 
                //     << " true pt2 = " << true2mu.pt << " reco pt2 = " << MuReco2.Pt() << " deltaR_t1r1 = " << deltaR_t1r1 << " deltaR_t1r2 = " << deltaR_t1r2 << endl;  

                for(int iPT = 0; iPT < NPThist; iPT++){
                for(int iETA = 0; iETA < NETAhist; iETA++){
                   int iK = iPT + iETA*NPThist_inv;
                   if(true1mu.pt >= PTbin[iPT] && true1mu.pt < PTbin[iPT+1] 
                      && fabs(true1mu.eta) >= ETAbin[iETA] && fabs(true1mu.eta) < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t1r1);  
                   if(true2mu.pt >= PTbin[iPT] && true2mu.pt < PTbin[iPT+1] 
                      && fabs(true2mu.eta) >= ETAbin[iETA] && fabs(true2mu.eta) < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t2r2);  
                   if(genZpostFSR.pt >= PTbin[iPT] && genZpostFSR.pt < PTbin[iPT+1] 
                      && fabs(genZpostFSR.eta) >= ETAbin[iETA] && fabs(genZpostFSR.eta) < ETAbin[iETA+1]) hDimuonRes[iK] -> Fill( (rPt-genZpostFSR.pt)/genZpostFSR.pt);  
                }} 
 
                // Z shape on GEN level as pt and eta for muon
                // Z mass and resolution as pt and eta for muon
                for(int iPT = 0; iPT < NPThist_inv; iPT++){
                for(int iETA = 0; iETA < NETAhist_inv; iETA++){
                for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
                   int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;

                   // muon tag minus and in tag eta region
                   if( 
                         ( MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]
                           //&& MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1] 
                           //&& fabs(reco1.eta) >= ETAbin_inv[iETA] && fabs(reco1.eta) < ETAbin_inv[iETA+1]
                           && reco1.eta >= ETAbin_inv[iETA] && reco1.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2.charge < 0.
                         )||
                         ( MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1]
                           //&& MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1] 
                           //&& fabs(reco2.eta) >= ETAbin_inv[iETA] && fabs(reco2.eta) < ETAbin_inv[iETA+1]
                           && reco2.eta >= ETAbin_inv[iETA] && reco2.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1.charge < 0.
                         ) ){
                             //if (reco2.charge < 0.) cout << " mu2- iPT = " << iPT << " pt1 = " << MuReco1.Pt() << " iETA = " << iETA << " eta1 = " << reco1.eta
                             //<< " iETA_tag = " << iETA_tag << " eta2 = " << reco2.eta <<  " charge mu2 = " << reco2.charge << endl;
                             //if (reco1.charge < 0.) cout << " mu2- invert iPT = " << iPT << " pt1 = " << MuReco2.Pt() << " iETA = " << iETA << " eta1 = " << reco2.eta
                             //<< " iETA_tag = " << iETA_tag << " eta2 = " << reco1.eta << " charge mu2 = " << reco1.charge << endl;
                             hDimuonMass[iK] -> Fill(rMassCorr);
                             if(rMassCorr < 80)hdeltaR80[iK] -> Fill(deltaR_r1r2);
                             if(rMassCorr > 85 && rMassCorr < 95) hdeltaR85_95[iK] -> Fill(deltaR_r1r2);
                             if(rMassCorr > 105) hdeltaR105[iK] -> Fill(deltaR_r1r2);
                            } 
                   // muon tag plus charge and in tag eta region
                   if( 
                         ( MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]
                           //&& MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1] 
                           //&& fabs(reco1.eta) >= ETAbin_inv[iETA] && fabs(reco1.eta) < ETAbin_inv[iETA+1]
                           && reco1.eta >= ETAbin_inv[iETA] && reco1.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2.charge > 0.
                         )||
                         ( MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1]
                           //&& MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]    
                           //&& fabs(reco2.eta) >= ETAbin_inv[iETA] && fabs(reco2.eta) < ETAbin_inv[iETA+1]
                           && reco2.eta >= ETAbin_inv[iETA] && reco2.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1.charge > 0.
                         ) )hDimuonMassTagMuPlus[iK] -> Fill(rMassCorr);

                   // GEN level
                   // muon tag minus and in tag eta region
                   if( 
                         ( true1mu.pt >= PTbin_inv[iPT] && true1mu.pt < PTbin_inv[iPT+1]
                           //&& true2mu.pt >= PTbin_inv[iPT] && true2mu.pt < PTbin_inv[iPT+1] 
                           //&& fabs(true1mu.eta) >= ETAbin_inv[iETA] && fabs(true1mu.eta) < ETAbin_inv[iETA+1]
                           && true1mu.eta >= ETAbin_inv[iETA] && true1mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true2mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true2mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true2mu.charge < 0.
                         )||
                         ( true2mu.pt >= PTbin_inv[iPT] && true2mu.pt < PTbin_inv[iPT+1] 
                           //&& true1mu.pt >= PTbin_inv[iPT] && true1mu.pt < PTbin_inv[iPT+1]
                           //&& fabs(true2mu.eta) >= ETAbin_inv[iETA] && fabs(true2mu.eta) < ETAbin_inv[iETA+1]
                           && true2mu.eta >= ETAbin_inv[iETA] && true2mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true1mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true1mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true1mu.charge < 0.
                         ) ){
                             hDimuonMassGEN[iK] -> Fill(genZpostFSR.mass);
                            } 

                   // GEN level
                   // muon tag plus charge and in tag eta region
                   if( 
                         ( true1mu.pt >= PTbin_inv[iPT] && true1mu.pt < PTbin_inv[iPT+1]
                           //&& true2mu.pt >= PTbin_inv[iPT] && true2mu.pt < PTbin_inv[iPT+1] 
                           //&& fabs(true1mu.eta) >= ETAbin_inv[iETA] && fabs(true1mu.eta) < ETAbin_inv[iETA+1]
                           && true1mu.eta >= ETAbin_inv[iETA] && true1mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true2mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true2mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true2mu.charge > 0.
                         )||
                         ( true2mu.pt >= PTbin_inv[iPT] && true2mu.pt < PTbin_inv[iPT+1] 
                           //&& true1mu.pt >= PTbin_inv[iPT] && true1mu.pt < PTbin_inv[iPT+1]
                           //&& fabs(true2mu.eta) >= ETAbin_inv[iETA] && fabs(true2mu.eta) < ETAbin_inv[iETA+1]
                           && true2mu.eta >= ETAbin_inv[iETA] && true2mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true1mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true1mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true1mu.charge > 0.
                         ) )hDimuonMassTagMuPlusGEN[iK] -> Fill(genZpostFSR.mass);


                }}}

// Phi binning
                for(int iPHI = 0; iPHI < NPHIbin_inv; iPHI++){
                for(int iETA = 0; iETA < NETAhist_inv; iETA++){
                for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
                   int iK = iPHI + iETA*NPHIbin_inv + iETA_tag*NETAhist_inv*NPHIbin_inv;
                   float phibin = 2./NPHIbin_inv;
                   float phi1 = reco1.phi;
                   float phi2 = reco2.phi;
                   if(phi1 < 0.) phi1 = 2*Pi + phi1; 
                   if(phi2 < 0.) phi2 = 2*Pi + phi2; 

                   // muon tag minus and in tag eta region
                   if( 
                         ( phi1 >= phibin*Pi*iPHI && phi1 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco1.eta) >= ETAbin_inv[iETA] && fabs(reco1.eta) < ETAbin_inv[iETA+1]
                           && reco1.eta >= ETAbin_inv[iETA] && reco1.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2.charge < 0.
                         )||
                         ( phi2 >= phibin*Pi*iPHI && phi2 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco2.eta) >= ETAbin_inv[iETA] && fabs(reco2.eta) < ETAbin_inv[iETA+1]
                           && reco2.eta >= ETAbin_inv[iETA] && reco2.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1.charge < 0.
                         ) ){
                             hPhiDimuonMass[iK] -> Fill(rMassCorr);
                            } 
                   // muon tag plus charge and in tag eta region
                   if( 
                         ( phi1 >= phibin*Pi*iPHI && phi1 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco1.eta) >= ETAbin_inv[iETA] && fabs(reco1.eta) < ETAbin_inv[iETA+1]
                           && reco1.eta >= ETAbin_inv[iETA] && reco1.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2.charge > 0.
                         )||
                         ( phi2 >= phibin*Pi*iPHI && phi2 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco2.eta) >= ETAbin_inv[iETA] && fabs(reco2.eta) < ETAbin_inv[iETA+1]
                           && reco2.eta >= ETAbin_inv[iETA] && reco2.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1.charge > 0.
                         ) ){
                             hPhiDimuonMassTagMuPlus[iK] -> Fill(rMassCorr);
                            }
 
                   // GEN level
                   phi1 = true1mu.phi;
                   phi2 = true2mu.phi; 
                   if(phi1 < 0.) phi1 = 2*Pi + phi1; 
                   if(phi2 < 0.) phi2 = 2*Pi + phi2; 
                   // GEN: muon tag minus and in tag eta region
                   if(
                         ( phi1 >= phibin*Pi*iPHI && phi1 < phibin*Pi*(iPHI+1)
                           //&& fabs(true1mu.eta) >= ETAbin_inv[iETA] && fabs(true1mu.eta) < ETAbin_inv[iETA+1]
                           && true1mu.eta >= ETAbin_inv[iETA] && true1mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true2mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true2mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true2mu.charge < 0.
                         )||
                         ( phi2 >= phibin*Pi*iPHI && phi2 < phibin*Pi*(iPHI+1)
                           //&& fabs(true2mu.eta) >= ETAbin_inv[iETA] && fabs(true2mu.eta) < ETAbin_inv[iETA+1]
                           && true2mu.eta >= ETAbin_inv[iETA] && true2mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true1mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true1mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true1mu.charge < 0.
                         ) ){
                             hPhiDimuonMassGEN[iK] -> Fill(genZpostFSR.mass);
                            }
                   // GEN: muon tag plus and in tag eta region
                   if(
                         ( phi1 >= phibin*Pi*iPHI && phi1 < phibin*Pi*(iPHI+1)
                           //&& fabs(true1mu.eta) >= ETAbin_inv[iETA] && fabs(true1mu.eta) < ETAbin_inv[iETA+1]
                           && true1mu.eta >= ETAbin_inv[iETA] && true1mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true2mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true2mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true2mu.charge > 0.
                         )||
                         ( phi2 >= phibin*Pi*iPHI && phi2 < phibin*Pi*(iPHI+1)
                           //&& fabs(true2mu.eta) >= ETAbin_inv[iETA] && fabs(true2mu.eta) < ETAbin_inv[iETA+1]
                           && true2mu.eta >= ETAbin_inv[iETA] && true2mu.eta < ETAbin_inv[iETA+1]
                           && fabs(true1mu.eta) >= ETAbin_tag[iETA_tag] && fabs(true1mu.eta) < ETAbin_tag[iETA_tag+1]
                           && true1mu.charge > 0.
                         ) ){
                             hPhiDimuonMassTagMuPlusGEN[iK] -> Fill(genZpostFSR.mass);
                            }

                }}}


        }


  TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 5);
  TF1* fitDoubleGauss2 = new TF1("fitDoubleGauss2", DoubleGauss, -0.1, 0.1, 5);
  for(int iPT = 0; iPT < NPThist; iPT++){
  for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
       // 
       fitDoubleGauss->SetParameters(0., 0.011, 0.041, 0.25, 2300.);
       if(PTbin[iPT] > 45 && ETAbin[iETA] > 0.9) fitDoubleGauss->SetParameters(0., 0.019, 0.16, 0.03, 320.); 
       fitDoubleGauss->FixParameter(0,0.);
       fitDoubleGauss->SetParLimits(1, 0.005, 0.05);//restrict sigma1
       fitDoubleGauss->SetParName(0,"mean");
       fitDoubleGauss->SetParName(1,"sig1");
       fitDoubleGauss->SetParName(2,"Asig2");
       fitDoubleGauss->SetParName(3,"sig2");
       fitDoubleGauss->SetParName(4,"Norm");
        
       //hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");
       hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");

       fitDoubleGauss2->SetParameters(fitDoubleGauss->GetParameter(0),fitDoubleGauss->GetParameter(1),fitDoubleGauss->GetParameter(2),fitDoubleGauss->GetParameter(3),fitDoubleGauss->GetParameter(4));
       if(PTbin[iPT] > 45 && ETAbin[iETA] > 0.9) fitDoubleGauss->SetParameters(0., 0.019, 0.16, 0.03, 320.);
       fitDoubleGauss2->FixParameter(0,0.);
       fitDoubleGauss2->SetParLimits(1, 0.005, 0.05);//restrict sigma1
       fitDoubleGauss2->SetParName(0,"mean");
       fitDoubleGauss2->SetParName(1,"sig1");
       fitDoubleGauss2->SetParName(2,"Asig2");
       fitDoubleGauss2->SetParName(3,"sig2");
       fitDoubleGauss2->SetParName(4,"Norm");


       hmuonRes[iK] -> Fit(fitDoubleGauss2,"RLE");
       mean[iK] = fitDoubleGauss2->GetParameter(0);  
       sig1[iK] = fitDoubleGauss2->GetParameter(1);  
       sig2[iK] = fitDoubleGauss2->GetParameter(3);  
       Asig2[iK] = fitDoubleGauss2->GetParameter(2) + ScaleAsig2*fitDoubleGauss2->GetParError(2);  
       if(Asig2[iK] < 0.) Asig2[iK] = 0.;  
       //
  }}        

  cout << "********************* MC DONE successfull ******************** " << endl;
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
        TChain* treeData = new TChain("tree");
        //DoubleMu DataSets:
        //treeData -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataDoubleMuRun2012A-13Jul2012-v1/merged/DoubleMuRun2012A-13July2012-v1_merged.root");//old
        // min set of cuts, vertex, cosmics
        if(HiggsNtuple == 1) { 
           /// 2012 year
           if(RunYear == "2012ABCsmall"){
              treeData -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/NtuplesDataSingleMuRun2012A-13Jul2012-v1/minimal/SingleMuRun2012A-13Jul2012-v1_minimal.root");
              treeData -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/NtuplesDataSingleMuRun2012B-13Jul2012-v1/minimal/SingleMuRun2012B-13Jul2012-v1_minimal.root");
              treeData -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/NtuplesDataSingleMuRun2012C-24Aug2012-v1/minimal/SingleMuRun2012C-24Aug2012-v1_minimal.root");
           }
        }
        // ***** 2011A
        if(RunYear == "2011A")treeData -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-01/NtuplesDataSingleMuRun2011A-08Nov2011-v1/minimal/SingleMuRun2011A-08Nov2011-v1_minimal.root");
        // ***** 2011B
        if(RunYear == "2011B")treeData -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-01/NtuplesDataSingleMuRun2011B-19Nov2011-v1/minimal/SingleMuRun2011B-19Nov2011-v1_minimal.root");
        //treeData -> AddFile("");
        if(HiggsNtuple == 0)
        treeData -> AddFile("/data/uftrig01b/digiovan/root/CMSSW_4_2_8_patch6/Ntuples/Data/SingleMuRun2011B-PromptReco-v1_AOD/minimal/boostedZ_2011B_PromptRecov1_allgood_175832to180252_minimal_m20to160.root"); 


        float rMass_data;
        float rPt_data;
        float rEta_data;
        _MuonInfo reco1_data, reco2_data;
        treeData->SetBranchAddress("reco1",&reco1_data);
        treeData->SetBranchAddress("reco2",&reco2_data);
        treeData->SetBranchAddress("recoCandMass",&rMass_data);
        treeData->SetBranchAddress("recoCandPt",&rPt_data);
        treeData->SetBranchAddress("recoCandEta",&rEta_data);
        cout << "Nevent DATA to process = " << treeData->GetEntries() << endl;

        //for( int k=0; k<10000; k++)
        for( int k=0; k<treeData->GetEntries(); k++)
        {
                //process progress
                if(k!=0 && (k%10000)==0)
                std::cout << "- processing DATA event " << k << "\r" << std::flush;

                treeData->GetEntry(k);

                if (reco1_data.charge == reco2_data.charge) continue;

                TLorentzVector MuReco1, MuReco2;
                MuReco1.SetPtEtaPhiM(reco1_data.pt, reco1_data.eta, reco1_data.phi, MASS_MUON);
                MuReco2.SetPtEtaPhiM(reco2_data.pt, reco2_data.eta, reco2_data.phi, MASS_MUON);
                float theta12 = MuReco1.Angle(MuReco2.Vect()); // angle between MuReco1 & MuReco2
                //if( theta12 >= (Pi-0.02) ) continue; 
                //cout << "Loose selection event = " << k << " MuReco1.Pt() = " << MuReco1.Pt() << " reco1.eta = " << reco1.eta << endl;
                //     << " reco1.numValidTrackerHits = " reco1.numValidTrackerHits << " reco1.trackIsoSumPt = " << reco1.trackIsoSumPt 
                //     << "reco1.d0 = " << reco1.d0 << endl;
                //if( !isLoose(reco1)    || !isLoose(reco2)    ) continue; // both loose

                ////////////////////////
                // Muon Correction
                //If you run MC, apply the muon momentum correction,"momcor_mc()" function (only for MC) 
                // This is for MC,
                //                 TLorentzVector, mu, corresponds to mu- for charge = -1 or mu+ for charge = +1 
                //                 sysdev == 0 returns the central value of muon momentum correction 
                //                 sysdev == 1 returns the muon momentum correction smeared 1 standard deviation
                //                             in <1/pt> correction 
                //                             (global factor is also changed by 1 sigma total error ( stat 0X syst)
                //                 If you want to assign the systematic error using x standard (statistical) deviation, 
                //                 you can assign sysdev == x.
                //                 runopt == 0 for 2011A correction 
                //                 runopt == 1 for 2011B correction
                //                TLor.  float   float   int             
                //rmcor->momcor_data(mu, charge, sysdev, runopt); 
  
                if (MuCorr == 1){ // Rochester Correction
                  if(RunYear == "2011B"){
                     rmcor->momcor_data(MuReco1, float(reco1_data.charge), 0, 1); 
                     rmcor->momcor_data(MuReco2, float(reco2_data.charge), 0, 1);
                  }
                  if(RunYear == "2011A"){
                     rmcor->momcor_data(MuReco1, float(reco1_data.charge), 0, 0); 
                     rmcor->momcor_data(MuReco2, float(reco2_data.charge), 0, 0);
                  }
                } 
/*                if (MuCorr == 2){
                   TLorentzVector* MuReco1Pointer = &MuReco1;//make pointer to MuReco1
                   TLorentzVector* MuReco2Pointer = &MuReco2;
                   if(reco1.charge < 0){
                        correctorData_->applyPtCorrection(*MuReco1Pointer,-1);
                   } 
                   else{
                        correctorData_->applyPtCorrection(*MuReco1Pointer,1);
                   }
                   if(reco2.charge < 0){
                        correctorData_->applyPtCorrection(*MuReco2Pointer,-1);
                   } 
                   else{
                        correctorData_->applyPtCorrection(*MuReco2Pointer,1);
                   }
                }
*/
                ////////////////////////
                ////////////////////////

                TLorentzVector MuRecoCand = MuReco1 + MuReco2; 
                float rMassCorr = MuRecoCand.M(); 
                if (rMassCorr <  60) continue;// mass restriction
                if (rMassCorr > 120) continue;

                //cout << "any selection event = " << k << " reco1_data.charge = " << reco1_data.charge << " rMassCorr = " << rMassCorr << endl;
                float pTcorr1 = MuReco1.Pt(); 
                float pTcorr2 = MuReco2.Pt(); 
                if( !isKinTight_2012(reco1_data, pTcorr1) || !isKinTight_2012(reco2_data, pTcorr2) ) continue;


                ROOT::Math::PtEtaPhiEVector MuReco1h(MuReco1.Pt(), reco1_data.eta, reco1_data.phi, MuReco1.Pt()); 
                ROOT::Math::PtEtaPhiEVector MuReco2h(MuReco2.Pt(), reco2_data.eta, reco2_data.phi, MuReco2.Pt()); 
                float deltaR_r1r2 = ROOT::Math::VectorUtil::DeltaR(MuReco1h, MuReco2h);  

                for(int iPT = 0; iPT < NPThist_inv; iPT++){
                for(int iETA = 0; iETA < NETAhist_inv; iETA++){
                for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
                   int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;

                   // muon tag minus and in tag eta region
                   if( 
                         ( MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]
                           //&& MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1] 
                           //&& fabs(reco1_data.eta) >= ETAbin_inv[iETA] && fabs(reco1_data.eta) < ETAbin_inv[iETA+1]
                           && reco1_data.eta >= ETAbin_inv[iETA] && reco1_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2_data.charge < 0.
                         )||
                         ( MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1] 
                           //&& MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]
                           //&& fabs(reco2_data.eta) >= ETAbin_inv[iETA] && fabs(reco2_data.eta) < ETAbin_inv[iETA+1]
                           && reco2_data.eta >= ETAbin_inv[iETA] && reco2_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1_data.charge < 0.
                         ) ){
                             hDimuonMassDATA[iK] -> Fill(rMassCorr);
                             if(rMassCorr < 80)hdeltaRDATA80[iK] -> Fill(deltaR_r1r2);
                             if(rMassCorr > 85 && rMassCorr < 95) hdeltaRDATA85_95[iK] -> Fill(deltaR_r1r2);
                             if(rMassCorr > 105) hdeltaRDATA105[iK] -> Fill(deltaR_r1r2);
                             }
                   // muon tag plus charge and in tag eta region
                   if( 
                         ( MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]
                           //&& MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1]
                           //&& fabs(reco1_data.eta) >= ETAbin_inv[iETA] && fabs(reco1_data.eta) < ETAbin_inv[iETA+1]
                           && reco1_data.eta >= ETAbin_inv[iETA] && reco1_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2_data.charge > 0.
                         )||
                         ( MuReco2.Pt() >= PTbin_inv[iPT] && MuReco2.Pt() < PTbin_inv[iPT+1] 
                           //&& MuReco1.Pt() >= PTbin_inv[iPT] && MuReco1.Pt() < PTbin_inv[iPT+1]
                           //&& fabs(reco2_data.eta) >= ETAbin_inv[iETA] && fabs(reco2_data.eta) < ETAbin_inv[iETA+1]
                           && reco2_data.eta >= ETAbin_inv[iETA] && reco2_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1_data.charge > 0.
                         ) )hDimuonMassTagMuPlusDATA[iK] -> Fill(rMassCorr);

                }}}
// Phi binning
                for(int iPHI = 0; iPHI < NPHIbin_inv; iPHI++){
                for(int iETA = 0; iETA < NETAhist_inv; iETA++){
                for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
                   int iK = iPHI + iETA*NPHIbin_inv + iETA_tag*NETAhist_inv*NPHIbin_inv;
                   float phibin = 2./NPHIbin_inv;
                   float phi1 = reco1_data.phi;
                   float phi2 = reco2_data.phi;
                   if(phi1 < 0.) phi1 = 2*Pi + phi1; 
                   if(phi2 < 0.) phi2 = 2*Pi + phi2; 

                   // muon tag minus and in tag eta region
                   if( 
                         ( phi1 >= phibin*Pi*iPHI && phi1 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco1_data.eta) >= ETAbin_inv[iETA] && fabs(reco1_data.eta) < ETAbin_inv[iETA+1]
                           && reco1_data.eta >= ETAbin_inv[iETA] && reco1_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2_data.charge < 0.
                         )||
                         ( phi2 >= phibin*Pi*iPHI && phi2 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco2_data.eta) >= ETAbin_inv[iETA] && fabs(reco2_data.eta) < ETAbin_inv[iETA+1]
                           && reco2_data.eta >= ETAbin_inv[iETA] && reco2_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1_data.charge < 0.
                         ) ){
                             hPhiDimuonMassDATA[iK] -> Fill(rMassCorr);
                            } 
                   // muon tag plus charge and in tag eta region
                   if( 
                         ( phi1 >= phibin*Pi*iPHI && phi1 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco1_data.eta) >= ETAbin_inv[iETA] && fabs(reco1_data.eta) < ETAbin_inv[iETA+1]
                           && reco1_data.eta >= ETAbin_inv[iETA] && reco1_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco2_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco2_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco2_data.charge > 0.
                         )||
                         ( phi2 >= phibin*Pi*iPHI && phi2 < phibin*Pi*(iPHI+1) 
                           //&& fabs(reco2_data.eta) >= ETAbin_inv[iETA] && fabs(reco2_data.eta) < ETAbin_inv[iETA+1]
                           && reco2_data.eta >= ETAbin_inv[iETA] && reco2_data.eta < ETAbin_inv[iETA+1]
                           && fabs(reco1_data.eta) >= ETAbin_tag[iETA_tag] && fabs(reco1_data.eta) < ETAbin_tag[iETA_tag+1]
                           && reco1_data.charge > 0.
                         ) ){
                             hPhiDimuonMassTagMuPlusDATA[iK] -> Fill(rMassCorr);
                            }
 
                }}}

        }




  //printhistos();

   theFile->cd();
   theFile->Write();
   theFile->Close();
}


///////////////////////////////////////////
void bookhistos(){
//	key1->GetXaxis()->SetTitle("P_{T}(#mu#mu) [GeV]");
//	key1->GetYaxis()->SetTitle("# Events per GeV");
//	test1->SetLineColor(2);
//	test1->SetMarkerStyle(2);
//	test1->SetMarkerColor(2);
//	raw1->SetLineColor(3);
   //int NbinMass = 60;
   int NbinMass = 120;
   hTestMC_ZPt = new TH1F("hTestMC_ZPt", "Z pT", 20, 0., 100.);
   hTestMC_ptTag = new TH1F("hTestMC_ptTag", "Z pT", 20, 0., 100.);
   hTestMC_ZMass = new TH1F("hTestMC_ZMass", "RECO", NbinMass, 60., 120.);
   hTestMC_ZMassGen = new TH1F("hTestMC_ZMassGen", "GEN", NbinMass, 60., 120.);
   hTestMC_DeltaR80 = new TH1F("hTestDeltaRDATA80","deltaR", 200, 0., 5.);
   hTestDeltaRMatch = new TH1F("hTestDeltaRMatch","deltaR", 200, 0., 5.);

   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      //hmuonRes[iK] = new TH1F(Form("hmuonRes%d",iK), Form("#mu resolution, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 200, -0.2, 0.2);
      hmuonRes[iK] = new TH1F(Form("hmuonRes%d",iK), Form("MC Drell-Yan, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 80, -0.1, 0.1);
      hDimuonRes[iK] = new TH1F(Form("hDimuonRes%d",iK), Form("#mu#mu resolution, PT: %4.1f-%4.1f GeV, ETA: %4.1f-%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 100, -0.2, 0.2);
   }}
   for(int iPT = 0; iPT < NPThist_inv; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
   for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
      int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;

      hDimuonMass[iK] = new TH1F(Form("hDimuonMass%d",iK), Form("MC DY, #mu^{+} p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f, #mu^{-} tag #eta: %4.1f#divide%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hDimuonMassGEN[iK] = new TH1F(Form("hDimuonMassGEN%d",iK), Form("DY GEN, #mu^{+} p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f, #mu^{-} tag #eta: %4.1f#divide%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hDimuonMassDATA[iK] = new TH1F(Form("hDimuonMassDATA%d",iK), Form("DATA, #mu^{+} p_{T}: %4.1f-%4.1f #eta: %4.1f#divide%4.1f, #mu^{-} tag #eta: %4.1f#divide%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hDimuonMassTagMuPlus[iK] = new TH1F(Form("hDimuonMassTagMuPlus%d",iK), Form("MC DY, #mu^{-} p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f, #mu^{+} tag #eta: %4.1f#divide%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hDimuonMassTagMuPlusGEN[iK] = new TH1F(Form("hDimuonMassTagMuPlusGEN%d",iK), Form("DY GEN, #mu^{-} p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f, #mu^{+} tag #eta: %4.1f#divide%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hDimuonMassTagMuPlusDATA[iK] = new TH1F(Form("hDimuonMassTagMuPlusDATA%d",iK), Form("DATA, #mu^{-} p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f, #mu^{+} tag #eta: %4.1f#divide%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);

      hdeltaR80[iK] = new TH1F(Form("hdeltaR80%d",iK), Form("MC DY, PT: %4.1f-%4.1f GeV, #mu^{+} ETA: %4.1f-%4.1f, #mu^{-} tag ETA: %4.1f-%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), 200, 0., 5.);
      hdeltaR85_95[iK] = new TH1F(Form("hdeltaR85_95%d",iK), Form("MC DY, PT: %4.1f-%4.1f GeV, #mu^{+} ETA: %4.1f-%4.1f, #mu^{-} tag ETA: %4.1f-%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), 200, 0., 5.);
      hdeltaR105[iK] = new TH1F(Form("hdeltaR105%d",iK), Form("MC DY, PT: %4.1f-%4.1f GeV, #mu^{+} ETA: %4.1f-%4.1f, #mu^{-} tag ETA: %4.1f-%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), 200, 0., 5.);
      hdeltaRDATA80[iK] = new TH1F(Form("hdeltaRDATA80%d",iK), Form("DATA, PT: %4.1f-%4.1f GeV, #mu^{+} ETA: %4.1f-%4.1f, #mu^{-} tag ETA: %4.1f-%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), 200, 0., 5.);
      hdeltaRDATA85_95[iK] = new TH1F(Form("hdeltaRDATA85_95%d",iK), Form("DATA, PT: %4.1f-%4.1f GeV, #mu^{+} ETA: %4.1f-%4.1f, #mu^{-} tag ETA: %4.1f-%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), 200, 0., 5.);
      hdeltaRDATA105[iK] = new TH1F(Form("hdeltaRDATA105%d",iK), Form("DATA, PT: %4.1f-%4.1f GeV, #mu^{+} ETA: %4.1f-%4.1f, #mu^{-} tag ETA: %4.1f-%4.1f ", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), 200, 0., 5.);
   }}}

// Phi binning
   for(int iPHI = 0; iPHI < NPHIbin_inv; iPHI++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
   for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
      int iK = iPHI + iETA*NPHIbin_inv + iETA_tag*NETAhist_inv*NPHIbin_inv;

      float phi = 2./NPHIbin_inv;
      hPhiDimuonMass[iK] = new TH1F(Form("hPhiDimuonMass%d",iK), Form("MC DY, #mu^{+} #phi: %4.1f#pi-%4.1f#pi GeV, |#eta|: %4.1f-%4.1f, #mu^{-} tag |#eta|: %4.1f-%4.1f ", iPHI*phi, (iPHI+1)*phi, ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hPhiDimuonMassGEN[iK] = new TH1F(Form("hPhiDimuonMassGEN%d",iK), Form("DY GEN, #mu^{+} #phi: %4.1f#pi-%4.1f#pi GeV, |#eta|: %4.1f-%4.1f, #mu^{-} tag |#eta|: %4.1f-%4.1f ", iPHI*phi, (iPHI+1)*phi, ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hPhiDimuonMassDATA[iK] = new TH1F(Form("hPhiDimuonMassDATA%d",iK), Form("DATA, #mu^{+} #phi: %4.1f#pi-%4.1f#pi GeV, |#eta|: %4.1f-%4.1f, #mu^{-} tag |#eta|: %4.1f-%4.1f ", iPHI*phi, (iPHI+1)*phi, ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hPhiDimuonMassTagMuPlus[iK] = new TH1F(Form("hPhiDimuonMassTagMuPlus%d",iK), Form("MC DY, #mu^{-} #phi: %4.1f#pi-%4.1f#pi GeV, |#eta|: %4.1f-%4.1f, #mu^{+} tag |#eta|: %4.1f-%4.1f ", iPHI*phi, (iPHI+1)*phi, ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hPhiDimuonMassTagMuPlusGEN[iK] = new TH1F(Form("hPhiDimuonMassTagMuPlusGEN%d",iK), Form("DY GEN, #mu^{-} #phi: %4.1f#pi-%4.1f#pi GeV, |#eta|: %4.1f-%4.1f, #mu^{+} tag |#eta|: %4.1f-%4.1f ", iPHI*phi, (iPHI+1)*phi, ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);
      hPhiDimuonMassTagMuPlusDATA[iK] = new TH1F(Form("hPhiDimuonMassTagMuPlusDATA%d",iK), Form("DATA, #mu^{-} #phi: %4.1f#pi-%4.1f#pi GeV, |#eta|: %4.1f-%4.1f, #mu^{+} tag |#eta|: %4.1f-%4.1f ", iPHI*phi, (iPHI+1)*phi, ETAbin_inv[iETA], ETAbin_inv[iETA+1], ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]), NbinMass, 60., 120.);

   }}}

        hdthetav1v2 = new TH1F("hdthetav1v2","",160,0., Pi);

// unfolding matrix:
        key1_MC = new TH1F("key1_MC","",18,binning_histo);
        test1 = new TH1F("test1","",18,binning_histo);
        key1 = new TH1F("key1","",18,binning_histo);
        key1->GetXaxis()->SetTitle("P_{T}(#mu#mu) [GeV]");
        key1->GetYaxis()->SetTitle("# Events per GeV");
        raw1 = new TH1F("raw1","",18,binning_histo);
        test1->SetLineColor(2);
        test1->SetMarkerStyle(2);
        test1->SetMarkerColor(2);
        raw1->SetLineColor(3);
 
}
///////////////////////////////////////////
void printhistos(){
   /// print muon resolution data, MC and MC sim
   gStyle->SetOptFit(1111);
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1); 

  //TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 5);
  //TF1* fitDoubleGauss2 = new TF1("fitDoubleGauss2", DoubleGauss, -0.1, 0.1, 5);

 /////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////
 TString Extra = "TagMuMinus";
 //fitDiMuon(hDimuonMass, hDimuonMassDATA, Extra);
 Extra = "TagMuPlus";
 //fitDiMuon(hDimuonMassTagMuPlus, hDimuonMassTagMuPlusDATA, Extra);

 /////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////

}
///////////////////////////////////////////

// for 2012
float getRelIso(_MuonInfo& muon, float pTcorr)
{
  // tracker iso cut
  //if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_2012;  

  float result = muon.sumChargedHadronPtR04 +
    std::max(0.0,muon.sumNeutralHadronEtR04 + muon.sumPhotonEtR04 - 0.5*muon.sumPUPtR04);

  return result/pTcorr;
} 

bool isKinTight_2012(_MuonInfo& muon, float pTcorr)
{ 
  
  bool isKinTight_2012=false;
  
  if (!muon.isGlobal)            return isKinTight_2012;
  if (!muon.isPFMuon)            return isKinTight_2012;
  
  // acceptance cuts
  if (pTcorr < 25)              return isKinTight_2012; // pt cut
  if (fabs(muon.eta) > 2.1)      return isKinTight_2012; // eta cut
  
  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_2012; // # hits in tracker
  
  if(getRelIso(muon, pTcorr) > 0.12)     return isKinTight_2012;
  
  // beam spot cut
  //if(HiggsNtuple == 1){
     if (fabs(muon.d0_PV) > 0.2)       return isKinTight_2012;
     if (fabs(muon.dz_PV) > 5)         return isKinTight_2012;
  //if(HiggsNtuple == 0)
  //   if (fabs(muon.d0) > 0.2)       return isKinTight_2012;
  //   if (fabs(muon.dz) > 5)         return isKinTight_2012;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2012;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2012;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012;
  if ( muon.normChiSquare > 10)     return isKinTight_2012;

  isKinTight_2012=true;
  return isKinTight_2012;
}


///////////////////////// Fit functions ////////////////////

Double_t BWnonrel(Double_t* x, Double_t* par){


   //Double_t eff = par[0]*(TMath::Erf(par[1]*x[0]-par[2]) - 1.) + par[3];
   //Non rel. Breit-Wigner                               mass     width 
   Double_t bw_nonrel = par[0]*(TMath::BreitWigner(x[0], par[1], par[2])); 
   return bw_nonrel;
}

Double_t Gauss(Double_t* x, Double_t* par){


   //Double_t eff = par[0]*(TMath::Erf(par[1]*x[0]-par[2]) - 1.) + par[3];
   //Non rel. Breit-Wigner                               mass     width 
   //Double_t bw_nonrel = par[0]*(TMath::BreitWigner(x[0], par[1], par[2])); 
   Double_t gauss = par[0]*(TMath::Gaus(x[0], par[1], par[2])); 
   return gauss;
}
Double_t FuncVoigtian(Double_t* x, Double_t* par){

 //--------------------------------------------------------------------//
    //Fit parameters:
    //par[0]=Width (scale) Breit-Wigner
    //par[1]=Most Probable (MP, location) Breit mean
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

       // Numeric constants
       Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
       Double_t twoPi = 6.2831853071795;//2Pi

       // Control constants
       Double_t np = 100.0;      // number of convolution steps
       Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
       //Double_t np = 60.0;      // number of convolution steps
       //Double_t sc =   3.0;      // convolution extends to +-sc Gaussian sigmas

       // Variables
       Double_t xx;
       Double_t fland;
       Double_t sum = 0.0;
       Double_t xlow,xupp;
       Double_t step;
       Double_t i;


       // Range of convolution integral
       xlow = x[0] - sc * par[3];
       xupp = x[0] + sc * par[3];

       step = (xupp-xlow) / np;

       // Convolution integral of Breit and Gaussian by sum
       for(i=1.0; i<=np/2; i++) {
          xx = xlow + (i-.5) * step;
          //fland = TMath::BreitWigner(xx,par[1],par[0]);
          fland = TMath::BreitWigner(xx,par[1],par[2]);
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

          xx = xupp - (i-.5) * step;
         //fland = TMath::BreitWigner(xx,par[1],par[0]);
         fland = TMath::BreitWigner(xx,par[1],par[2]);
          sum += fland * TMath::Gaus(x[0],xx,par[3]);
       }

       return (par[0] * step * sum * invsq2pi / par[3]);
}

Double_t FuncVoigtianBG(Double_t* x, Double_t* par){

 //--------------------------------------------------------------------//
    //Fit parameters:
    //par[0]=Width (scale) Breit-Wigner
    //par[1]=Most Probable (MP, location) Breit mean
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

       // Numeric constants
       Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
       Double_t twoPi = 6.2831853071795;//2Pi

       // Control constants
       Double_t np = 100.0;      // number of convolution steps
       Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
       //Double_t np = 60.0;      // number of convolution steps
       //Double_t sc =   3.0;      // convolution extends to +-sc Gaussian sigmas

       // Variables
       Double_t xx;
       Double_t fland;
       Double_t sum = 0.0;
       Double_t xlow,xupp;
       Double_t step;
       Double_t i;


       // Range of convolution integral
       xlow = x[0] - sc * par[3];
       xupp = x[0] + sc * par[3];

       step = (xupp-xlow) / np;

       // Convolution integral of Breit and Gaussian by sum
       for(i=1.0; i<=np/2; i++) {
          xx = xlow + (i-.5) * step;
          //fland = TMath::BreitWigner(xx,par[1],par[0]);
          fland = TMath::BreitWigner(xx,par[1],par[2]);
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

          xx = xupp - (i-.5) * step;
         //fland = TMath::BreitWigner(xx,par[1],par[0]);
         fland = TMath::BreitWigner(xx,par[1],par[2]);
          sum += fland * TMath::Gaus(x[0],xx,par[3]);
       }

       //return (par[0] * step * sum * invsq2pi / par[3] + exp(par[4]-par[5]*x[0]+par[6]*x[0]*x[0]));
       return (par[0] * step * sum * invsq2pi / par[3] + exp(par[4]-par[5]*x[0]+par[6]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0]));
}

// Sum of background and peak function
Double_t DoubleGauss(Double_t *x, Double_t *par) 
{ 

    //return BWnonrel(x,par) + Gauss(x,&par[3]);
    Double_t dgauss = 0.;
    //dgauss = Gauss(x,par) + Gauss(x,&par[3]);

     if(par[1] < par[3]){//not normalized gauss both gauss are  = 1 at x[0]=par[0]  
               dgauss =  exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]);
               dgauss = dgauss + par[2]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[3]/par[3]);
               dgauss = par[4]*dgauss;
     }
     //if (par[2]> 0.) dgauss = dgauss/(1+par[2]);  

    return dgauss;
}

//------------------------------------------------------------------------------
// Draw projections and residuals
//------------------------------------------------------------------------------
void DrawWithRatio(TCanvas *canvas, char *cTitle,
                  TH1F *hNum, TH1F *hDen)
{
 
 // sanity check
 if (hNum->GetNbinsX() != hDen->GetNbinsX()){
   std::cout<< " *** Error: binning not consistent between data"
            << " and MC -> Exit!\n";
   return;
 }

 
 hNum->Sumw2();
 hDen->Sumw2();

 
 TH1F *hPull = (TH1F*)hNum ->Clone("hPull");
 hPull->Sumw2();
 hPull->Divide(hDen);
 
 //----------------------------------------------------------------------------
 // Create the pads
 //----------------------------------------------------------------------------
 TPad* pad1;
 TPad* pad2;
 
 pad1 = new TPad("pad1","This is pad1",0.02,0.30,0.98,0.98,0);
 pad2 = new TPad("pad2","This is pad2",0.02,0.01,0.98,0.29,0);
 
 pad1->SetLogx();
 pad2->SetLogx();
 pad1->SetBottomMargin(0.01);
 pad2->SetBottomMargin(0.33);
 pad2->SetTopMargin   (0.10);
 
 pad1->Draw(); // Projections pad
 pad2->Draw(); // Residuals   pad
        
        _leg2 = new TLegend(.68,.77,.98,.93);
        _leg2->AddEntry(hDen,"Gen ","l");
        _leg2->AddEntry(hNum,"Unfolded Reco","p");
 pad1->cd();
 hDen->Draw("histo");
 hNum->Draw("pe same");
 _leg2->Draw(); 
 PrintItLog(pad1,cTitle);

//   TLegend* leg = SetLegend(0.73, 0.7, 0.92, 0.89);
//   leg -> AddEntry(hDen," no mass cut","f");
//   leg -> AddEntry(hNum," 60 < M^{#mu #mu} < 120","p");
//   leg ->Draw("same");
 //----------------------------------------------------------------------------
 // Residuals pad
 //----------------------------------------------------------------------------
 pad2->cd();
 
 TAxis *xPull = NULL;
 TAxis *yPull = NULL;
 char xAxisName[200];
 sprintf(xAxisName,"%s",hDen->GetXaxis()->GetTitle());
 axis1F(hPull,xPull,yPull,xAxisName,"ratio");
 
 if (hPull->GetMaximum() > 100) {
   hPull->SetMinimum(-100);
   hPull->SetMaximum( 100);
 }
 
 hPull->GetXaxis()->SetLabelOffset(0.005);
 hPull->GetXaxis()->SetLabelSize  (0.11);
 hPull->GetXaxis()->CenterTitle(1);
 hPull->GetXaxis()->SetTitleOffset(1.10);
 hPull->GetXaxis()->SetTitleSize  (0.12);
 hPull->GetXaxis()->SetNdivisions(7);

 hPull->GetYaxis()->SetLabelSize  (0.09);
 hPull->GetYaxis()->CenterTitle(1);
 hPull->GetYaxis()->SetTitleOffset(0.5);
 hPull->GetYaxis()->SetTitleSize  (0.12);

 hPull->SetMaximum(1.5);
 hPull->SetMinimum(0.5);
 hPull->Draw("pe");

 pad2->Update();
 pad2->GetFrame()->DrawClone();

}

void axis1F(TH1F  *histo,
           TAxis *xaxis,
           TAxis *yaxis,
           char  *xtitle,
           char  *ytitle)
{
 histo->SetMarkerSize(0.5);
 histo->SetMarkerStyle(kFullCircle);

 histo->SetTitle("");

 xaxis = histo->GetXaxis();
 yaxis = histo->GetYaxis();

 xaxis->SetLabelFont(42);
 yaxis->SetLabelFont(42);
 xaxis->SetLabelOffset(0.005);
 yaxis->SetLabelOffset(0.005);
//xaxis->SetLabelSize(0.04);
 yaxis->SetLabelSize(0.04);

 xaxis->SetNdivisions(505);
 yaxis->SetNdivisions(505);

 xaxis->SetTitle(xtitle);
 yaxis->SetTitle(ytitle);
 xaxis->SetTitleColor(kBlack);
 yaxis->SetTitleColor(kBlack);
 xaxis->SetTitleFont(42);
 yaxis->SetTitleFont(42);
 xaxis->SetTitleOffset(1.0);
 yaxis->SetTitleOffset(1.3);
//xaxis->SetTitleSize(0.045);
//yaxis->SetTitleSize(0.045);
 yaxis->CenterTitle(kTRUE);
}

///////////////////////////
///////////////////////////
///////////////////////////
void fitDiMuon(TH1F* histoMC[], TH1F* histoDATA[], TString Extra){

   gStyle->SetOptFit(1111);
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1); 
  /// print Mass reco in data, MC and MC sim     
  Double_t yScale[Nhist_inv];
  Double_t xScale[Nhist_inv];
  Double_t eyScale[Nhist_inv];
  Double_t exScale[Nhist_inv];
  Double_t mass_Data[Nhist_inv];
  Double_t mass_Data_err[Nhist_inv];
  Double_t mass_MC[Nhist_inv];
  Double_t mass_MC_err[Nhist_inv];
  Double_t res_Data[Nhist_inv];
  Double_t res_Data_err[Nhist_inv];
  Double_t res_MC[Nhist_inv];
  Double_t res_MC_err[Nhist_inv];

  //TF1* fitBWnonrel = new TF1("fitBWnonrel", BWnonrel, -0.1, 0.1, 3);
  //TF1* fitGauss = new TF1("fitGauss", Gauss, -0.1, 0.1, 3);
  TF1* fitFuncVoigtian = new TF1("fitFuncVoigtian", FuncVoigtian, 80., 110., 4);
  //TF1* fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 60., 120., 7);
  TF1* fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 60., 120., 8);
  fitFuncVoigtianBG ->SetParName(0, "Norm");
  fitFuncVoigtianBG ->SetParName(1, "mass");
  fitFuncVoigtianBG ->SetParName(2, "width");
  fitFuncVoigtianBG ->SetParName(3, "resol");
  fitFuncVoigtianBG ->SetParName(4, "BG4");
  fitFuncVoigtianBG ->SetParName(5, "BG5");
  fitFuncVoigtianBG ->SetParName(6, "BG6");
  fitFuncVoigtianBG ->SetParName(7, "BG7");
  //Double_t par[6];
  //TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 6);

   TString gifname_inv[NPThist_inv][NETAhist_tag];
   TString gifname_invDATA[NPThist_inv][NETAhist_tag];
   TString gifmethodinv = "Voigtian";
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1);
   for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
   for(int iPT = 0; iPT < NPThist_inv; iPT++){
     //gifname_inv[iPT] = "test.png";
     gifname_inv[iPT][iETA_tag] = Form("plots/SingleMuMassInv%deta%dpt"+Extra+"MC_PtCorr%dSmearPt%d_"+gifmethodinv+"_pt%d_etatag%d", NETAhist_inv, NPThist_inv, MuCorr, Ismear, iPT, iETA_tag);
     gifname_invDATA[iPT][iETA_tag] = Form("plots/SingleMuMassInv%deta%dpt"+Extra+"DATA_PtCorr%dSmearPt%d_"+gifmethodinv+"_pt%d_etatag%d", NETAhist_inv, NPThist_inv, MuCorr, Ismear, iPT, iETA_tag);
    gifname_inv[iPT][iETA_tag] = gifname_inv[iPT][iETA_tag];
    gifname_invDATA[iPT][iETA_tag] = gifname_invDATA[iPT][iETA_tag];
     //cout << "test" <<endl;
     //cout << gifname_inv[iPT] <<endl;
     TCanvas *c21 = new TCanvas("c21","Resolution mass",3000,1700);
     TCanvas *c21_data = new TCanvas("c21_data","Resolution mass",3000,1700);
     //c21-> Divide(4,2);
     //c21_data-> Divide(4,2);
     //c21-> Divide(2,1);
     //c21_data-> Divide(2,1);
     c21-> Divide(4,2);
     c21_data-> Divide(4,2);
     int Neta = NETAhist_inv;
     //if(Neta > 7) {Neta = 7; cout << "more then 8 bins, change canvas size" << endl;}
     //if(Neta > 2) {Neta = 2; cout << "more then 2 bins, change canvas size" << endl;}
     if(Neta > 8) {Neta = 8; cout << "more then 8 bins, change canvas size" << endl;}
     for(int iETA = 0; iETA < Neta; iETA++){
       int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;

       //fitFuncVoigtian ->SetParameters(200., 90., 2.495, 1.);
       //fitFuncVoigtian->SetParLimits(1, 80., 100.);
       //fitFuncVoigtian->FixParameter(2, 2.495);
       //fitFuncVoigtian->SetParLimits(3, 0.2, 5.);

       fitFuncVoigtianBG->SetParameters(200., 90., 2.495, 1., 2., -0.01, 0.);
       if (PTbin_inv[iPT] < 20) fitFuncVoigtianBG->SetParameters(200., 90.9, 2.495, 1.9, -19., -0.6, -0.004);
       //fitFuncVoigtianBG->SetParameters(2000., 90., 2.495, 1., -20., -0.6, -0.0036);
       fitFuncVoigtianBG->SetParLimits(0, 0., 2000000.);
       fitFuncVoigtianBG->SetParLimits(1, 80., 100.);
       fitFuncVoigtianBG->FixParameter(2, 2.495);
       fitFuncVoigtianBG->SetParLimits(3, 0.2, 5.);
       //fitFuncVoigtianBG->SetParLimits(4, -200., 200.);
       fitFuncVoigtianBG->SetParLimits(5, -10., 10.);//BG5 
       fitFuncVoigtianBG->SetParLimits(6, -1., 1.);//BG6
       fitFuncVoigtianBG->SetParLimits(7, -1., 1.);//BG7

       ///////////////////////////// 
       c21 -> cd(iETA+1);
       histoMC[iK] -> Fit(fitFuncVoigtianBG,"RLE");//Fit invmass sim with Scale Facotor
       //histoMC[iK] -> Fit(fitFuncVoigtianBG,"RL");//Fit invmass sim with Scale Facotor
       Double_t par0 = fitFuncVoigtianBG ->GetParameter(0);
       Double_t par1 = fitFuncVoigtianBG ->GetParameter(1);// mass
       Double_t par1err = fitFuncVoigtianBG ->GetParError(1);// error mass
       Double_t par2 = fitFuncVoigtianBG ->GetParameter(2);
       Double_t par3 = fitFuncVoigtianBG ->GetParameter(3);// resolution
       Double_t par3err = fitFuncVoigtianBG ->GetParError(3);// error resolution
       Double_t par4 = fitFuncVoigtianBG ->GetParameter(4);
       Double_t par5 = fitFuncVoigtianBG ->GetParameter(5);
       Double_t par6 = fitFuncVoigtianBG ->GetParameter(6);
       Double_t par7 = fitFuncVoigtianBG ->GetParameter(7);
       //TF1* fitFuncVoigtianBGDraw = new TF1("fitFuncVoigtianBGDraw", FuncVoigtianBG, 60., 120., 7);
       TF1* fitFuncVoigtianBGDraw = new TF1("fitFuncVoigtianBGDraw", FuncVoigtianBG, 60., 120., 8);
       fitFuncVoigtianBGDraw -> SetParameter(0, par0);
       fitFuncVoigtianBGDraw -> SetParameter(1, par1);
       fitFuncVoigtianBGDraw -> SetParameter(2, par2);
       fitFuncVoigtianBGDraw -> SetParameter(3, par3);
       fitFuncVoigtianBGDraw -> SetParameter(4, par4);
       fitFuncVoigtianBGDraw -> SetParameter(5, par5);
       fitFuncVoigtianBGDraw -> SetParameter(6, par6);
       fitFuncVoigtianBGDraw -> SetParameter(7, par7);

       //fitFuncVoigtianBG -> SetLineColor(2);//red
       histoMC[iK] -> SetLineColor(kBlue);
       histoMC[iK] -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
       histoMC[iK] -> GetYaxis()->SetTitle("Entries");
       histoMC[iK] -> Draw("hist");
       fitFuncVoigtianBGDraw -> Draw("same");
       //histoMC[iK] -> Draw("esame");
       //TF1 *fpol1 = new TF1 ("fpol1", "[3]*exp([0]-[1]*x+[2]*x*x)", 60., 120.);
       TF1 *fpol1 = new TF1 ("fpol1", "exp([0]-[1]*x+[2]*x*x+[3]*x*x*x)", 60., 120.);
       fpol1 -> SetParameter(0,par4);
       fpol1 -> SetParameter(1,par5);
       fpol1 -> SetParameter(2,par6);
       fpol1 -> SetParameter(3,par7);
       //fpol1 -> SetParameter(3,par0);
       fpol1 ->SetLineColor(4);
       fpol1 -> Draw("same");
        TLegend* histMCinfo = SetLegend(.56,.6,1.,.7);
        histMCinfo->AddEntry(histoMC[iK],Form("MC reco + Fit"),"l");
        histMCinfo->AddEntry(fitFuncVoigtianBGDraw, "Fit ","l");
        //histMCinfo->AddEntry(histoMC[iK], "MC reco","lep");
        histMCinfo-> Draw("same");



       /////////////////
       //fitFuncVoigtianBG->FixParameter(4, par4); //fix from MC because not enough statistics
       fitFuncVoigtianBG->FixParameter(5, par5); //fix from MC because not enough statistics
       fitFuncVoigtianBG->FixParameter(6, par6); //fix from MC because not enough statistics
       fitFuncVoigtianBG->FixParameter(7, par7); //fix from MC because not enough statistics
       c21_data -> cd(iETA+1);
       histoDATA[iK] -> Fit(fitFuncVoigtianBG,"RLE");
       //histoDATA[iK] -> Fit(fitFuncVoigtianBG,"RL");
       Double_t par0_data = fitFuncVoigtianBG ->GetParameter(0);

       Double_t par1_data = fitFuncVoigtianBG ->GetParameter(1);// mass
       Double_t par1err_data = fitFuncVoigtianBG ->GetParError(1);// error mass
       Double_t par3_data = fitFuncVoigtianBG ->GetParameter(3); // resolution
       Double_t par3err_data = fitFuncVoigtianBG ->GetParError(3);// error resolution

       Double_t par4_data = fitFuncVoigtianBG ->GetParameter(4);
       Double_t par5_data = fitFuncVoigtianBG ->GetParameter(5);
       Double_t par6_data = fitFuncVoigtianBG ->GetParameter(6);
       Double_t par7_data = fitFuncVoigtianBG ->GetParameter(7);
       TF1 *fpol1_data = new TF1 ("fpol1_data", "exp([0]-[1]*x+[2]*x*x+[3]*x*x*x)", 60., 120.);
       fpol1_data -> SetParameter(0,par4_data);
       fpol1_data -> SetParameter(1,par5_data);
       fpol1_data -> SetParameter(2,par6_data);
       fpol1_data -> SetParameter(3,par7_data);
       //fpol1_data -> SetParameter(3,par0_data);
       fpol1_data ->SetLineColor(4);

       histoDATA[iK] -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
       histoDATA[iK] -> GetYaxis()->SetTitle("Entries");
       histoDATA[iK] -> Draw("e");
       histoMC[iK] -> SetLineColor(kBlue);
       histoMC[iK] -> DrawNormalized("samehist", histoDATA[iK]->Integral());
       fpol1_data -> Draw("same");
        TLegend* histDATAinfo = SetLegend(.6,.6, 1.,.7);
        histDATAinfo->AddEntry(histoMC[iK],Form("MC reco "),"l");
        histDATAinfo->AddEntry(histoDATA[iK], "DATA reco + Fit","lep");
        histDATAinfo->AddEntry(fitFuncVoigtianBG, "Fit","l");
        histDATAinfo-> Draw("same");

       /////////////////
       //begin: calculate and fill scale factor
       xScale[iK] = iK+1;
       exScale[iK] = 0.5;
       yScale[iK] = 0;
       eyScale[iK] = 0;
       if(par3_data !=0 && par3 != 0){
           yScale[iK] = par3_data/par3;
           eyScale[iK] = yScale[iK]* pow(( par3err*par3err/par3/par3 + par3err_data*par3err_data/par3_data/par3_data  ),0.5);
       }
       //end: calculate and fill scale factor
       mass_Data[iK] = par1_data;
       mass_Data_err[iK] = par1err_data;
       mass_MC[iK] = par1;
       mass_MC_err[iK] = par1err;
       res_Data[iK] = par3_data;
       res_Data_err[iK] = par3err_data;
       res_MC[iK] = par3;
       res_MC_err[iK] = par3err;
     }
     c21 ->Print(gifname_inv[iPT][iETA_tag]+".png");
     c21 ->Print(gifname_inv[iPT][iETA_tag]+".root");
     c21_data ->Print(gifname_invDATA[iPT][iETA_tag]+".png");
     c21_data ->Print(gifname_invDATA[iPT][iETA_tag]+".root");
   }} //end iETA_tag and iPT
   ////////////////////
   for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
     TF1 *fScale = new TF1 ("fScale", "[0]", 1., Nhist_inv1 );
     //TGraphErrors* grScale = new TGraphErrors(Nhist_inv1,xScale, yScale, exScale, eyScale);
     TGraphErrors* grScale = new TGraphErrors(Nhist_inv1,&xScale[0], &yScale[(iETA_tag*NETAhist_inv*NPThist_inv)], &exScale[0], &eyScale[(iETA_tag*NETAhist_inv*NPThist_inv)]);
     grScale->SetTitle(Form("Data/MC "));
     grScale->SetMarkerColor(4);
     grScale->SetMarkerStyle(21);

     TString gifname_ScaleFaclor[NETAhist_tag];
     TString gifname_Mass[NETAhist_tag];
     TString gifname_Res[NETAhist_tag];
     TCanvas *cScale = new TCanvas("cScale","Resolution mass",800,600);
     cScale -> cd();

     TH2F *h2 = new TH2F("h","Axes",Nhist_inv1+1,0,Nhist_inv1+1,100,0.,2.0);
     h2->GetXaxis()->SetNdivisions(Nhist_inv1);
     h2->SetTitle(Form("Data/MC, Tag #mu^{-} #eta: %4.1f#divide%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
     if(Extra == "TagMuPlus")
        h2->SetTitle(Form("Data/MC, Tag #mu^{+} #eta: %4.1f#divide%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
     gifname_ScaleFaclor[iETA_tag]= Form("plots/SingleMuScaleFactor%deta%dpt"+Extra+"DATA_PtCorr%dSmearPt%d_etatag%d_"+gifmethodinv, NETAhist_inv, NPThist_inv, MuCorr, Ismear, iETA_tag);
     gifname_Mass[iETA_tag] = Form("plots/SingleMuZMassPeak%deta%dpt"+Extra+"DATA_PtCorr%dSmearPt%d_etatag%d_"+gifmethodinv, NETAhist_inv, NPThist_inv, MuCorr, Ismear, iETA_tag);
     gifname_Res[iETA_tag] = Form("plots/SingleMuZResPeak%deta%dpt"+Extra+"DATA_PtCorr%dSmearPt%d_etatag%d_"+gifmethodinv, NETAhist_inv, NPThist_inv, MuCorr, Ismear, iETA_tag);
   for(int iPT = 0; iPT < NPThist_inv; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;
      int iKin = iPT + iETA*NPThist_inv; 
      if (iPT == 0) h2->GetXaxis()->SetBinLabel(iKin, Form("p_{T}: %4.0f-%4.0f GeV, #eta: %4.1f#divide%4.1f", PTbin_inv[0], PTbin_inv[NPThist_inv], ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
   }}

     h2->Draw();
     for(int iPT = 0; iPT < NPThist_inv; iPT++){
     for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        int iK = iPT + iETA*NPThist_inv+1;
        if (iPT == 0){
           float Xm = h2 -> GetBinLowEdge(iK+1);
           TLine *lineGrid = new TLine(Xm,0.,Xm,2.0);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");  
        }
     }}

     //grScale->Draw("ALP");
     grScale->Draw("LP");
     grScale->Fit("fScale", "R");
     cScale ->Print(gifname_ScaleFaclor[iETA_tag]+".png");
     cScale ->Print(gifname_ScaleFaclor[iETA_tag]+".root");
  /// end: print Mass reco in data, MC and MC sim     

//////////////////////
/// mass bihavior:

   //TGraphErrors* grMass = new TGraphErrors(Nhist_inv1,xScale, mass_Data, exScale, mass_Data_err);
   TGraphErrors* grMass = new TGraphErrors(Nhist_inv1,&xScale[0], &mass_Data[(iETA_tag*NETAhist_inv*NPThist_inv)], &exScale[0], &mass_Data_err[(iETA_tag*NETAhist_inv*NPThist_inv)]);
   grMass->SetTitle(Form("Z Mass peak"));
   grMass->SetMarkerColor(4);// blue
   grMass->SetMarkerStyle(21);
   //TGraphErrors* grMassMC = new TGraphErrors(Nhist_inv1,xScale, mass_MC, exScale, mass_MC_err);
   TGraphErrors* grMassMC = new TGraphErrors(Nhist_inv1,&xScale[0], &mass_MC[(iETA_tag*NETAhist_inv*NPThist_inv)], &exScale[0], &mass_MC_err[(iETA_tag*NETAhist_inv*NPThist_inv)]);
   grMassMC->SetTitle(Form("Z Mass peak"));
   grMassMC->SetMarkerColor(2);//red
   grMassMC->SetMarkerStyle(20);

     TCanvas *cScale1 = new TCanvas("cScale1","Resolution mass",800,600);
     cScale1 -> cd();

   //TH2F *h2Mass = new TH2F("h","Axes",Nhist_inv+1,0,Nhist_inv+1,100,89.,92.);
   //TH2F *h2Mass = new TH2F("h","Axes",Nhist_inv1+1,0,Nhist_inv1+1,100,90.,91.5);
   TH2F *h2Mass = new TH2F("h","Axes",Nhist_inv1+1,0,Nhist_inv1+1,100,89.,93.);
        //h2Mass ->  GetXaxis() ->SetNdivisions(505);// n = n1 + 100*n2 + 10000*n3   
   h2Mass->GetXaxis()->SetNdivisions(Nhist_inv1);
   h2Mass->SetTitle(Form("Z Mass peak, Tag #mu^{-} #eta: %4.1f#divide%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   if(Extra == "TagMuPlus")
      h2Mass->SetTitle(Form("Z Mass peak, Tag #mu^{+} #eta: %4.1f#divide%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   for(int iPT = 0; iPT < NPThist_inv; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;
      int iKin = iPT + iETA*NPThist_inv+1;
      //h2Mass->GetXaxis()->SetBinLabel(iK, Form("p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      if (iPT == 0){
          h2Mass->GetXaxis()->SetBinLabel(iKin, Form("p_{T}: %4.0f-%4.0f GeV, #eta: %4.1f#divide%4.1f", PTbin_inv[0], PTbin_inv[NPThist_inv],ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      }
   }}

     h2Mass->Draw();

     for(int iPT = 0; iPT < NPThist_inv; iPT++){
     for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        int iK = iPT + iETA*NPThist_inv+1;
        if (iPT == 0){
           float Xm = h2Mass -> GetBinLowEdge(iK+1);
           //TLine *lineGrid = new TLine(Xm,90.,Xm,91.5);
           TLine *lineGrid = new TLine(Xm,89.,Xm,93.);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");  
        }
     }}

     //grScale->Draw("ALP");
     grMass->Draw("LP");
     fScale->SetLineColor(4); // blue
     grMass->Fit("fScale", "R");
     grMassMC->Draw("LP");
     //grMassMC->Fit("fScale", "R");
     TLegend* legMass = SetLegend(.2,.2,0.6,.35);
     legMass->AddEntry(grMass, " Data","ep");
     legMass->AddEntry(fScale, " Fit Data","l");
     legMass->AddEntry(grMassMC, " MC","ep");
     legMass -> Draw("same");
     cScale1 ->Print(gifname_Mass[iETA_tag]+".png");
     cScale1 ->Print(gifname_Mass[iETA_tag]+".root");

///////////////////
/// res bihavior:

   //TGraphErrors* grRes = new TGraphErrors(Nhist_inv1,xScale, res_Data, exScale, res_Data_err);
   TGraphErrors* grRes = new TGraphErrors(Nhist_inv1,&xScale[0], &res_Data[(iETA_tag*NETAhist_inv*NPThist_inv)], &exScale[0], &res_Data_err[(iETA_tag*NETAhist_inv*NPThist_inv)]);
   grRes->SetTitle(Form(" Resolution near Z peak"));
   grRes->SetMarkerColor(4);
   grRes->SetMarkerStyle(21);
   //TGraphErrors* grResMC = new TGraphErrors(Nhist_inv1,xScale, res_MC, exScale, res_MC_err);
   TGraphErrors* grResMC = new TGraphErrors(Nhist_inv1,&xScale[0], &res_MC[(iETA_tag*NETAhist_inv*NPThist_inv)], &exScale[0], &res_MC_err[(iETA_tag*NETAhist_inv*NPThist_inv)]);
   grResMC->SetTitle(Form(" Resolution near Z peak"));
   grResMC->SetMarkerColor(2);//red
   grResMC->SetMarkerStyle(20);

     TCanvas *cScale2 = new TCanvas("cScale2","Resolution res",800,600);
     cScale2 -> cd();

   //TH2F *h2Res = new TH2F("h","Axes",Nhist_inv1+1,0,Nhist_inv1+1,100,0.,3.);
   TH2F *h2Res = new TH2F("h","Axes",Nhist_inv1+1,0,Nhist_inv1+1,100,0.,4.);
   h2Res->GetXaxis()->SetNdivisions(Nhist_inv1);
   h2Res->SetTitle(Form("Resolution near Z peak, Tag #mu^{-} #eta: %4.1f#divide%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   if(Extra == "TagMuPlus")
     h2Res->SetTitle(Form("Resolution near Z peak, Tag #mu^{+} #eta: %4.1f#divide%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   for(int iPT = 0; iPT < NPThist_inv; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;
      int iKin = iPT + iETA*NPThist_inv+1;      
      //h2Res->GetXaxis()->SetBinLabel(iK, Form("p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f#divide%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      if (iPT == 0) h2Res->GetXaxis()->SetBinLabel(iKin, Form("p_{T}: %4.0f-%4.0f GeV, #eta: %4.1f#divide%4.1f", PTbin_inv[0], PTbin_inv[NPThist_inv],ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
   }}

     h2Res->Draw();
     for(int iPT = 0; iPT < NPThist_inv; iPT++){
     for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        int iK = iPT + iETA*NPThist_inv+1;
        if (iPT == 0){
           float Xm = h2Res -> GetBinLowEdge(iK+1);
           //TLine *lineGrid = new TLine(Xm,0.,Xm,3.);
           TLine *lineGrid = new TLine(Xm,0.,Xm,4.);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");  
        }
     }}

     //grScale->Draw("ALP");
     grRes->Draw("LP");
     fScale->SetLineColor(4); // blue
     //grRes->Fit("fScale", "R");
     grResMC->Draw("LP");
     TLegend* legRes = SetLegend(.2,.2,0.6,.3);
     legRes->AddEntry(grRes, " Data","ep");
     //legRes->AddEntry(fScale, " Fit Data","l");
     legRes->AddEntry(grResMC, " MC","ep");
     legRes -> Draw("same");
     cScale2 ->Print(gifname_Res[iETA_tag]+".png");
     cScale2 ->Print(gifname_Res[iETA_tag]+".root");
     delete h2Res;
     delete h2Mass;
     delete h2; 
///////////////////

  } // end iETA_tag 
}
///////////////////////////
///////////////////////////
///////////////////////////
