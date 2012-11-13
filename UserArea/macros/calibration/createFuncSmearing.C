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
// to write to file
#include <iostream>
#include <fstream>

//include Rochechester correction
#include<rochcor.h>
#include<rochcor2012.h>
#include "MuScleFitCorrector.h"
// include smearing tool
#include <SmearingTool.h>

#include "./DataFormat.h"
using namespace std;

  int MuCorr = 1; // 0 - no muon correction, 
                  // 1 - Rochester Correction,
                  // 2 - MuscleFit correcton 
  int Ismear = 0; // 0 - take pt reco from MC
                  // 1 - smear pt with own tools using pt gen post FSR  
  
  TString RunYear = "2012"; // 2012; 2011A; 2011B;  

  //int RunYear = 2012; // 2011 or 2012
  const float PTbin[] = {25., 30., 35., 40., 45., 50., 70., 100., 150., 300.}; //default
  const float ETAbin[] = {-2.1, -1.6, -1.2, -0.8, 0., 0.8, 1.2, 1.6, 2.1};
  const int NPThist = (sizeof(PTbin)/sizeof(float)-1);
  const int NETAhist = (sizeof(ETAbin)/sizeof(float)-1);
  const int Nhist = NPThist*NETAhist;


  const float ScaleFactor = 1.;
  //const float ScaleFactor = 1.1;
  // variation Asig2 by 3 sigma;
  const float ScaleAsig2 = 0;// by default + 0 sigma of Asig2;
  //const float ScaleAsig2 = 3.;// + 3 sigma of Asig2;
  //const float ScaleAsig2 = -3.;// + 3 sigma of Asig2;
   TString Extra = "";
   //TString Extra = "_AsigPlus3sigma";
   //TString Extra = "_AsigMinus3sigma";
   //TString Extra = "_MATGRAPH";

  Double_t MASS_MUON = 0.105658367;    //GeV/c2
  Double_t mean[Nhist];
  Double_t sig1[Nhist];
  Double_t sig2[Nhist];
  Double_t Asig2[Nhist]; 
  Double_t ERRmean[Nhist];
  Double_t ERRsig1[Nhist];
  Double_t ERRsig2[Nhist];
  Double_t ERRAsig2[Nhist]; 
  Double_t ResRMS[Nhist]; 
  Double_t ErrResRMS[Nhist]; 
// ---- CP error ----

TH1F* hmuonRes[Nhist];
TH1F* hmuonResNonCorr[Nhist];

// unfolding matrix

TLegend *_leg2;



void bookhistos();     // to book histograms
void printhistos();
Double_t DoubleGauss(Double_t*, Double_t* );
Double_t BWnonrel(Double_t*, Double_t* );
Double_t Gauss(Double_t*, Double_t* );
Double_t FuncVoigtian(Double_t*, Double_t* );
Double_t FuncVoigtianBG(Double_t*, Double_t* );
bool isKinTight_2011(_MuonInfo& muon,float pTcorr);
bool isKinTight_2012(_MuonInfo& muon, float pTcorr);
float getRelIso(_MuonInfo& muon, float pTcorr);

void DrawWithRatio(TCanvas *canvas, char *cTitle,
                  TH1F *hNum, TH1F *hDen);

void axis1F(TH1F  *histo,
           TAxis *xaxis,
           TAxis *yaxis,
           char  *xtitle,
           char  *ytitle);



class createFuncSmearing
{
public :
   //constructor
   createFuncSmearing(){}
   void main();
};

void createFuncSmearing::main(){


	gROOT->Clear();
  	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	userStyle();
	//modifiedStyle();
  	// ---- open the MC files ----
  	TChain* treeMC = new TChain("tree");
           if(RunYear == "2012" || RunYear == "2012ABCsmall") treeMC -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root");
           if(RunYear == "2011A" || RunYear == "2011B" || RunYear == "2011")treeMC -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-01/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1/minimal/DYToMuMu_minimal.root");

        TFile *theFile    =new TFile(Form("NtupleCreateFuncSmearing.root", ScaleFactor), "RECREATE");

        bookhistos();// histo init 
        
        float rMass;
        float rPt;
        float rEta;
        _MuonInfo reco1, reco2;
        _TrackInfo true1mu, true2mu;
        _genPartInfo genZpostFSR;
        treeMC->SetBranchAddress("reco1",&reco1);
        treeMC->SetBranchAddress("reco2",&reco2);
        treeMC->SetBranchAddress("recoCandMass",&rMass);
        treeMC->SetBranchAddress("recoCandPt",&rPt);
        treeMC->SetBranchAddress("recoCandEta",&rEta);
        treeMC->SetBranchAddress("genZpostFSR",&genZpostFSR);
        treeMC->SetBranchAddress("genM1ZpostFSR",&true1mu);
        treeMC->SetBranchAddress("genM2ZpostFSR",&true2mu);

        //////////////////////////////////////////////////////////// 
        //Rochester correction 
        //To get the central value of the momentum correction 
        rochcor *rmcor = new rochcor(); // make the pointer of rochcor class
        rochcor2012 *rmcor2012 = new rochcor2012(); // make the pointer of rochcor class
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
        //////////////////////////////////////////////////////////// 

        cout << "Nevent to process = " << treeMC->GetEntries() << endl;
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

                ////////////////////////
                //check that generator part is filled 
                if(nbad_tpt < 10 && (true1mu.pt < -900. || true2mu.pt < -900)){//rejection for MATGRAPH
                   cout << "CHECK ntuple, possible problem with fill gen level while RECO level is fine: true1mu.pt = " << true1mu.pt << endl; 
                   nbad_tpt++;
                }
                if(true1mu.pt < -900. || true2mu.pt < -900) continue;//rejection for MATGRAPH
                if(nbad_tMass < 10 && genZpostFSR.mass < -900.){//rejection for MATGRAPH
                   cout << "CHECK ntuple, possible EXTRA problem with fill gen level while RECO level is fine: tMass = " << genZpostFSR.mass << endl; 
                   nbad_tMass++;
                }
                if(genZpostFSR.mass < -900.) continue;//rejection for MATGRAPH

                ////////////////////////
                if (reco1.charge == reco2.charge) continue;
 
                TLorentzVector MuReco1, MuReco2, MuTrue1, MuTrue2;
                float MuTrue1charge = true1mu.charge;
                float MuTrue2charge = true2mu.charge;
                MuReco1.SetPtEtaPhiM(reco1.pt, reco1.eta, reco1.phi, MASS_MUON);
                MuReco2.SetPtEtaPhiM(reco2.pt, reco2.eta, reco2.phi, MASS_MUON);
                MuTrue1.SetPtEtaPhiM(true1mu.pt, true1mu.eta, true1mu.phi, MASS_MUON);
                MuTrue2.SetPtEtaPhiM(true2mu.pt, true2mu.eta, true2mu.phi, MASS_MUON);

                ////////////////////////
                Float_t muonRes_t1r1 = -10.;
                Float_t muonRes_t1r2 = -10.;
                Float_t muonRes_t2r1 = -10.;
                Float_t muonRes_t2r2 = -10.;
                if(true1mu.pt > 0.) muonRes_t1r1 = (MuReco1.Pt()-true1mu.pt)/true1mu.pt;
                if(true2mu.pt > 0.) muonRes_t2r2 = (MuReco2.Pt()-true2mu.pt)/true2mu.pt;
                if(true1mu.pt > 0.) muonRes_t1r2 = (MuReco2.Pt()-true1mu.pt)/true1mu.pt;
                if(true2mu.pt > 0.) muonRes_t2r1 = (MuReco1.Pt()-true2mu.pt)/true2mu.pt;
                float deltaR_t1r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco1);
                float deltaR_t1r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco2);
                float deltaR_t2r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco1);
                float deltaR_t2r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco2);
                //make matching between MuTrue1,2 and MuReco1,2:
                if(deltaR_t1r1 > deltaR_t1r2){
                     MuTrue1.SetPtEtaPhiM(true2mu.pt, true2mu.eta, true2mu.phi, MASS_MUON);
                     MuTrue1charge = true2mu.charge;
                     muonRes_t1r1 = muonRes_t2r1;
                }
                if(deltaR_t2r2 > deltaR_t2r1){
                     MuTrue2.SetPtEtaPhiM(true1mu.pt, true1mu.eta, true1mu.phi, MASS_MUON);
                     MuTrue2charge = true1mu.charge;
                     muonRes_t2r2 = muonRes_t1r2;
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
                // for 2012 v4.1: rochcor2012::momcor_mc( TLorentzVector& mu, float charge, float sysdev, int runopt, float& qter)

                if (MuCorr == 1){
                  float err_corr = 1.;
                  if(RunYear == "2011B"){
                     rmcor->momcor_mc(MuReco1, float(reco1.charge), float(0), 1, err_corr);
                     rmcor->momcor_mc(MuReco2, float(reco2.charge), float(0), 1, err_corr);
                  }
                  if(RunYear == "2011A"){
                     rmcor->momcor_mc(MuReco1, float(reco1.charge), float(0), 0, err_corr);
                     rmcor->momcor_mc(MuReco2, float(reco2.charge), float(0), 0, err_corr);
                  }
                  if(RunYear == "2012" || RunYear == "2012ABCsmall" || RunYear == "2012ABC"){
                     rmcor2012->momcor_mc(MuReco1, float(reco1.charge), float(0), 0, err_corr);
                     rmcor2012->momcor_mc(MuReco2, float(reco2.charge), float(0), 0, err_corr);
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
                Float_t muonResCorr_t1r1 = -10.;
                Float_t muonResCorr_t2r2 = -10.;
                if(MuTrue1.Pt() > 0.) muonResCorr_t1r1 = (MuReco1.Pt()-MuTrue1.Pt())/MuTrue1.Pt();
                if(MuTrue2.Pt() > 0.) muonResCorr_t2r2 = (MuReco2.Pt()-MuTrue2.Pt())/MuTrue2.Pt();
                TLorentzVector MuRecoCand = MuReco1 + MuReco2;
                float rMassCorr = MuRecoCand.M();

                if (rMassCorr <  60) continue;
                if (rMassCorr > 120) continue;
                //cout << "Loose selection event = " << k << " reco1.pt = " << reco1.pt << " reco1.eta = " << reco1.eta << endl;
                //     << " reco1.numValidTrackerHits = " reco1.numValidTrackerHits << " reco1.trackIsoSumPt = " << reco1.trackIsoSumPt 
                //     << "reco1.d0 = " << reco1.d0 << endl;



                float pTcorr1 = MuReco1.Pt();
                float pTcorr2 = MuReco2.Pt();
                if( RunYear == "2012" && (!isKinTight_2012(reco1, pTcorr1) || !isKinTight_2012(reco2, pTcorr2)) ) continue;
                if( (RunYear == "2011A" || RunYear == "2011B")&& (!isKinTight_2011(reco1, pTcorr1) || !isKinTight_2011(reco2, pTcorr2)) ) continue;

                for(int iPT = 0; iPT < NPThist; iPT++){
                for(int iETA = 0; iETA < NETAhist; iETA++){
                   int iK = iPT + iETA*NPThist;
                   if(MuTrue1.Pt() >= PTbin[iPT] && MuTrue1.Pt() < PTbin[iPT+1] 
                      && MuTrue1.Eta() >= ETAbin[iETA] && MuTrue1.Eta() < ETAbin[iETA+1]) {
                           hmuonRes[iK] -> Fill(muonResCorr_t1r1);
                           hmuonResNonCorr[iK] -> Fill(muonRes_t1r1);
                   }  
                   if(MuTrue2.Pt() >= PTbin[iPT] && MuTrue2.Pt() < PTbin[iPT+1] 
                      && MuTrue2.Eta() >= ETAbin[iETA] && MuTrue2.Eta() < ETAbin[iETA+1]){
                           hmuonRes[iK] -> Fill(muonResCorr_t2r2);  
                           hmuonResNonCorr[iK] -> Fill(muonRes_t2r2);
                   }
                }} 



        }

   //printhistos();
   theFile->cd();
   theFile->Write();
   theFile->Close();
}


///////////////////////////////////////////
void bookhistos(){

   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      hmuonRes[iK] = new TH1F(Form("hmuonRes%d",iK), Form("MC Drell-Yan, PT: %4.1f-%4.1f GeV, ETA: %4.1f#divide%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 80, -0.1, 0.1);
      hmuonResNonCorr[iK] = new TH1F(Form("hmuonResNonCorr%d",iK), Form("MC Drell-Yan, PT: %4.1f-%4.1f GeV, ETA: %4.1f#divide%4.1f", PTbin[iPT], PTbin[iPT+1], ETAbin[iETA], ETAbin[iETA+1]), 80, -0.1, 0.1);
   }}

}
///////////////////////////////////////////
void printhistos(){
   /// print muon resolution data, MC and MC sim
   gStyle->SetOptFit(1111);
   TString gifname[NPThist];
   //TString gifname_inv[NPThist];
   TString gifmethod = "DoubleGauss";
   //TString gifmethodinv = "DoubleGauss";
   //TString gifmethod = "BWnonrel";
   //TString gifmethodinv = "Voigtian";
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1); 

  TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 5);
  TF1* fitDoubleGauss2 = new TF1("fitDoubleGauss2", DoubleGauss, -0.1, 0.1, 5);
  for(int iPT = 0; iPT < NPThist; iPT++){
    gifname[iPT] = Form("plots/ResolutionPTScaleFactor%1.2f_"+gifmethod+"_%d", ScaleFactor, iPT);
    gifname[iPT] = gifname[iPT]+Extra;
    TCanvas *c20 = new TCanvas("c20","Resolution",3000,1700);
    TCanvas *c21 = new TCanvas("c21","Resolution mass",3000,1700);
    c20-> Divide(4,2);
    c21-> Divide(4,2);
    int Neta = NETAhist;
    if(Neta > 8) {Neta = 8; cout << "more then 8 bins, change canvas size" << endl;}
    for(int iETA = 0; iETA < Neta; iETA++){
      int iK = iPT + iETA*NPThist;
       c20 -> cd(iETA+1);
       if(ETAbin[iETA] < 0.) c20 -> cd(int(Neta/2) - iETA); 
       // 
       fitDoubleGauss->SetParameters(0., 0.85*hmuonRes[iK] -> GetRMS(), 0.07,0.85*hmuonRes[iK] -> GetRMS()+0.015, 1000.);
       //if(ETAbin[iETA] > 0.75 || ETAbin[iETA] < -0.85) fitDoubleGauss->SetParameters(0., hmuonRes[iK] -> GetRMS(), 0.035, 0.03, 1000.); 
       //if(PTbin[iPT] > 90.) fitDoubleGauss->SetParameters(0., 0.024, 0.08, 0.054, 160.); 
       //if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 0.75 || ETAbin[iETA] < -0.85) ) 
       if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 1.15 || ETAbin[iETA] < -1.25) ) 
           || (PTbin[iPT] > 145.) ){ 
          fitDoubleGauss->SetParameters(0., hmuonRes[iK] -> GetRMS(), 0.0, 0.03, 100.); 
          fitDoubleGauss->FixParameter(2,0.);//fix Asig2 to 0
          fitDoubleGauss->FixParameter(3,3.); // fix sig2 to any big value (sig1 < sig2)
          fitDoubleGauss->SetParLimits(1, 0.005, 0.1);//restrict sigma1
       }
       else{
          fitDoubleGauss->SetParLimits(2, 0.01, 0.4);//restrict Asig2
          fitDoubleGauss->SetParLimits(3,0.005, 0.1); //restrict sig2
          fitDoubleGauss->SetParLimits(1, 0.005, hmuonRes[iK] -> GetRMS());//restrict sigma1
       } 
       //fitDoubleGauss->FixParameter(0,0.); // fix mean of resolution
       fitDoubleGauss->SetParName(0,"mean");
       fitDoubleGauss->SetParName(1,"sig1");
       fitDoubleGauss->SetParName(2,"Asig2");
       fitDoubleGauss->SetParName(3,"sig2");
       fitDoubleGauss->SetParName(4,"Norm");

       hmuonRes[iK] -> Fit(fitDoubleGauss,"RLE");
       fitDoubleGauss2->SetParameters(fitDoubleGauss->GetParameter(0),fitDoubleGauss->GetParameter(1),fitDoubleGauss->GetParameter(2),fitDoubleGauss->GetParameter(3),fitDoubleGauss->GetParameter(4));
       //if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 0.75 || ETAbin[iETA] < -0.85) ) 
       if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 1.15 || ETAbin[iETA] < -1.25) ) 
           || (PTbin[iPT] > 145.) ){ 
          fitDoubleGauss2->FixParameter(2,0.);//fix Asig2 to 0
          fitDoubleGauss2->FixParameter(3,3.); // fix sig2 to any big value (sig1 < sig2)
          fitDoubleGauss2->SetParLimits(1, 0.005, 0.1);//restrict sigma1
       }
       else{
          fitDoubleGauss2->SetParLimits(2, 0.01, 0.4);//restrict Asig2
          fitDoubleGauss2->SetParLimits(3,0.005, 0.1); //restrict sig2
          fitDoubleGauss2->SetParLimits(1, 0.005, hmuonRes[iK] -> GetRMS());//restrict sigma1
       } 
       //fitDoubleGauss2->FixParameter(0,0.); // fix mean of resolution
       fitDoubleGauss2->SetParName(0,"mean");
       fitDoubleGauss2->SetParName(1,"sig1");
       fitDoubleGauss2->SetParName(2,"Asig2");
       fitDoubleGauss2->SetParName(3,"sig2");
       fitDoubleGauss2->SetParName(4,"Norm");

       hmuonRes[iK] -> Fit(fitDoubleGauss2,"RLE");

        TLegend* histinfo = SetLegend(.6,.57,1.,.73);
        histinfo->AddEntry(hmuonRes[iK], "MC reco + Fit","lep");
        histinfo->AddEntry(fitDoubleGauss2, Form("Fit, RMS = %4.3f#pm%4.3f", ResRMS[iK], ErrResRMS[iK]),"l");


       hmuonRes[iK] ->  GetXaxis() ->SetNdivisions(505);// n = n1 + 100*n2 + 10000*n3
       hmuonRes[iK] -> GetXaxis()->SetTitle("p_{T}^{reco}(#mu)/p_{T}^{gen}(#mu)-1");
       hmuonRes[iK] -> GetYaxis()->SetTitle("Entries");

       hmuonRes[iK] -> Draw("e");
       histinfo -> Draw("same");
       //c21 -> cd(iETA+1);
       //hDimuonRes[iK] -> Draw("e");
     }
     c20 ->Print(gifname[iPT]+".png");
     c20 ->Print(gifname[iPT]+".root");
     //c21 ->Print(gifname_inv[iPT]); 
   }
  /// end: print muon resolution data, MC and MC sim

}
///////////////////////////////////////////

bool isKinTight_2011(_MuonInfo& muon, float pTcorr) {
  bool isKinTight_2011=false;
  // minimum requirement to start investigating: to be loose
  if (!muon.isPFMuon) return isKinTight_2011;
  if (!muon.isGlobal) return isKinTight_2011;
  // acceptance cuts
  if (pTcorr < 25)              return isKinTight_2011; // pt cut
  if(getRelIso(muon, pTcorr) > 0.12)     return isKinTight_2011;

  if (fabs(muon.eta) > 2.1)      return isKinTight_2011; // eta cut
  if ( muon.normChiSquare >= 10) return isKinTight_2011;
  if (muon.numTrackerLayers < 6) return isKinTight_2011; // # hits in tracker
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2011;
  if (fabs(muon.d0_PV) >= 0.2) return isKinTight_2011;
  if (fabs(muon.dz_PV) > 5) return isKinTight_2011;
  if ( muon.numOfMatchedStations < 2  ) return isKinTight_2011;
  if ( muon.numValidMuonHits < 1 ) return isKinTight_2011;
  isKinTight_2011=true;
  return isKinTight_2011;
}

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

       //return (par[0] * step * sum * invsq2pi / par[3] + par[4] + par[5]*x[0]);
       return (par[0] * step * sum * invsq2pi / par[3] + exp(par[4]-par[5]*x[0]+par[6]*x[0]*x[0]));
       //return (par[0] * step * sum * invsq2pi / par[3] + par[0] * exp(par[4]-par[5]*x[0]+par[6]*x[0]*x[0]));
}

// Sum of background and peak function
Double_t DoubleGauss(Double_t *x, Double_t *par) 
{ 

    //return BWnonrel(x,par) + Gauss(x,&par[3]);
    Double_t dgauss = 0.;
    //dgauss = Gauss(x,par) + Gauss(x,&par[3]);

     // sig1 < sig2 always
     //if(par[1] < par[3])dgauss = par[4]* ( TMath::Gaus(x[0],par[0],par[1],"kTRUE") + par[2] * TMath::Gaus(x[0],par[0],par[3],"kTRUE") );
     //if(par[1] < par[3]){//normalized gauss par[2] could correlated with par[3] here 
     //          dgauss =  1/par[1]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]);
     //          dgauss = dgauss + par[2]/par[3]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[3]/par[3]);
     //          dgauss = par[4]*dgauss;
     //}
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

