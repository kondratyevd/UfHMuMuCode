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

#include "./DataFormat.h"
using namespace std;

  int RunYear = 2012; // 2011 or 2012
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



void createFuncSmearing(){


	gROOT->Clear();
  	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	userStyle();
	//modifiedStyle();
  	// ---- open the MC files ----
  	TChain* treeMC = new TChain("tree");
        treeMC -> AddFile("/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root");
        TFile *theFile    =new TFile(Form("NtupleCreateFuncSmearing.root", ScaleFactor), "RECREATE");

        //create txt file with fit output
        ofstream myfile ("FuncSmearing.txt");

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
                if (reco1.charge == reco2.charge) continue;
                if (rMass <  60) continue;
                if (rMass > 120) continue;
                //cout << "Loose selection event = " << k << " reco1.pt = " << reco1.pt << " reco1.eta = " << reco1.eta << endl;
                //     << " reco1.numValidTrackerHits = " reco1.numValidTrackerHits << " reco1.trackIsoSumPt = " << reco1.trackIsoSumPt 
                //     << "reco1.d0 = " << reco1.d0 << endl;

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

                TLorentzVector MuReco1, MuReco2, MuTrue1, MuTrue2;
                MuReco1.SetPtEtaPhiM(reco1.pt, reco1.eta, reco1.phi, MASS_MUON);
                MuReco2.SetPtEtaPhiM(reco2.pt, reco2.eta, reco2.phi, MASS_MUON);
                MuTrue1.SetPtEtaPhiM(true1mu.pt, true1mu.eta, true1mu.phi, MASS_MUON);
                MuTrue2.SetPtEtaPhiM(true2mu.pt, true2mu.eta, true2mu.phi, MASS_MUON);

                float pTcorr1 = MuReco1.Pt();
                float pTcorr2 = MuReco2.Pt();
                if( RunYear == 2012 && (!isKinTight_2012(reco1, pTcorr1) || !isKinTight_2012(reco2, pTcorr2)) ) continue;
                if( RunYear == 2011 && (!isKinTight_2011(reco1, pTcorr1) || !isKinTight_2011(reco2, pTcorr2)) ) continue;
                if (reco1.pt < 20 )      continue; // pt cut
                if (reco2.pt < 20 )      continue; // pt cut
                if ((reco1.trackIsoSumPt)/reco1.pt >=0.1) continue; // isolation
                if ((reco2.trackIsoSumPt)/reco2.pt >=0.1) continue; // isolation


                Float_t muonRes_t1r1 = -10.; 
                Float_t muonRes_t1r2 = -10.;
                Float_t muonRes_t2r1 = -10.; 
                Float_t muonRes_t2r2 = -10.;
                if(true1mu.pt > 0.) muonRes_t1r1 = (reco1.pt-true1mu.pt)/true1mu.pt;  
                if(true1mu.pt > 0.) muonRes_t1r2 = (reco2.pt-true1mu.pt)/true1mu.pt;  
                if(true2mu.pt > 0.) muonRes_t2r2 = (reco2.pt-true2mu.pt)/true2mu.pt;
                if(true2mu.pt > 0.) muonRes_t2r1 = (reco1.pt-true2mu.pt)/true2mu.pt;
                float deltaR_t1r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco1);  
                float deltaR_t1r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco2);  
                float deltaR_t2r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco1);  
                float deltaR_t2r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco2);  
                Float_t muonResCorr_t1r1 = muonRes_t1r1;
                Float_t muonResCorr_t2r2 = muonRes_t2r2;
                if(deltaR_t1r1 > deltaR_t1r2) muonResCorr_t1r1 = muonRes_t1r2;
                if(deltaR_t2r2 > deltaR_t2r1) muonResCorr_t2r2 = muonRes_t2r1;

                for(int iPT = 0; iPT < NPThist; iPT++){
                for(int iETA = 0; iETA < NETAhist; iETA++){
                   int iK = iPT + iETA*NPThist;
                   if(true1mu.pt >= PTbin[iPT] && true1mu.pt < PTbin[iPT+1] 
                      //&& fabs(true1mu.eta) >= ETAbin[iETA] && fabs(true1mu.eta) < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t1r1);  
                      && true1mu.eta >= ETAbin[iETA] && true1mu.eta < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t1r1);  
                   if(true2mu.pt >= PTbin[iPT] && true2mu.pt < PTbin[iPT+1] 
                      //&& fabs(true2mu.eta) >= ETAbin[iETA] && fabs(true2mu.eta) < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t2r2);  
                      && true2mu.eta >= ETAbin[iETA] && true2mu.eta < ETAbin[iETA+1]) hmuonRes[iK] -> Fill(muonResCorr_t2r2);  
                }} 



        }


  TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 5);
  TF1* fitDoubleGauss2 = new TF1("fitDoubleGauss2", DoubleGauss, -0.1, 0.1, 5);
  for(int iPT = 0; iPT < NPThist; iPT++){
  for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
       // 
       
       fitDoubleGauss->SetParameters(0., 0.85*hmuonRes[iK] -> GetRMS(), 0.07,0.85*hmuonRes[iK] -> GetRMS()+0.015, 1000.);
       //fitDoubleGauss->SetParameters(0., hmuonRes[iK] -> GetRMS(), 0.07,hmuonRes[iK] -> GetRMS()+0.015, 1000.);
       //if(ETAbin[iETA] > 0.75 || ETAbin[iETA] < -0.85) fitDoubleGauss->SetParameters(0., hmuonRes[iK] -> GetRMS(), 0.035, 0.03, 1000.); 
       //if(PTbin[iPT] > 90.) fitDoubleGauss->SetParameters(0., 0.024, 0.08, 0.054, 160.); 
       //if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 0.75 || ETAbin[iETA] < -0.85) ) 
       if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 1.15 || ETAbin[iETA] < -1.25) ) 
           || (PTbin[iPT] > 145) ){ 
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
       ResRMS[iK] = hmuonRes[iK] -> GetRMS();
       ErrResRMS[iK] = hmuonRes[iK] -> GetRMSError();
       mean[iK] = fitDoubleGauss2->GetParameter(0);  
       sig1[iK] = fitDoubleGauss2->GetParameter(1);  
       sig2[iK] = fitDoubleGauss2->GetParameter(3);  
       Asig2[iK] = fitDoubleGauss2->GetParameter(2);  
       ERRmean[iK] = fitDoubleGauss2->GetParError(0);  
       ERRsig1[iK] = fitDoubleGauss2->GetParError(1);  
       ERRsig2[iK] = fitDoubleGauss2->GetParError(3);  
       ERRAsig2[iK] = fitDoubleGauss2->GetParError(2);  
       //if(Asig2[iK] < 0.) Asig2[iK] = 0.;  
       //
  }}        

        if (myfile.is_open())
        {
          myfile << "   Float_t PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen);\n"; 
          myfile << "   Double_t DoubleGauss(Double_t*, Double_t* );\n";

          myfile << "   float PTbin[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
              if (iPT != (NPThist - 1) ) myfile << PTbin[iPT] << ", ";
              else myfile << PTbin[iPT]; 
          }
          myfile << "};\n";

          myfile << "   float ETATbin[]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
              if (iETA != (NETAhist - 1) ) myfile << ETAbin[iETA] << ", ";
              else myfile << ETAbin[iETA]; 
          }
          myfile << "};\n\n";
          myfile << "   const int NPThist = (sizeof(PTbin)/sizeof(float)-1);\n";
          myfile << "   const int NETAhist = (sizeof(ETAbin)/sizeof(float)-1);\n";
          myfile << "   const int Nhist = NPThist*NETAhist;\n";

          myfile << "   ////// Smearing parametrization for single muon:\n";  
          //////////////////////////////////////////
          myfile << "   float mean[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << mean[iK] << ", ";
              else myfile << mean[iK]; 
          }}
          myfile << "};\n";
          myfile << "   float sig1[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << sig1[iK] << ", ";
              else myfile << sig1[iK];
          }}
          myfile << "};\n";
          myfile << "   float sig2[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << sig2[iK] << ", ";
              else myfile << sig2[iK];
          }}
          myfile << "};\n";
          myfile << "   float Asig2[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << Asig2[iK] << ", ";
              else myfile << Asig2[iK];
          }}
          myfile << "};\n\n";
          //////////////////////////////////////////
          myfile << "   float ERRmean[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRmean[iK] << ", ";
              else myfile << ERRmean[iK];
          }}
          myfile << "};\n";
          myfile << "   float ERRsig1[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRsig1[iK] << ", ";
              else myfile << ERRsig1[iK];
          }}
          myfile << "};\n";
          myfile << "   float ERRsig2[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRsig2[iK] << ", ";
              else myfile << ERRsig2[iK];
          }}
          myfile << "};\n";
          myfile << "   float ERRAsig2[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRAsig2[iK] << ", ";
              else myfile << ERRAsig2[iK];
          }}
          myfile << "};\n\n";
          //////////////////////////////////////////
          myfile << "   float ResRMS[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ResRMS[iK] << ", ";
              else myfile << ResRMS[iK];
          }}
          myfile << "};\n";
          myfile << "   float ErrResRMS[]={";

          for(int iPT = 0; iPT < NPThist; iPT++){
          for(int iETA = 0; iETA < NETAhist; iETA++){
              int iK = iPT + iETA*NPThist;
              if(iETA == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ErrResRMS[iK] << ", ";
              else myfile << ErrResRMS[iK];
          }}
          myfile << "};\n\n";
          //////////////////////////////////////////
          myfile << "   ////// End smearing parametrization for single muon:\n\n";  

          myfile << "Float_t PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen){\n";
          myfile << "   ////// Make smearing for single muon from Gen level PostFSR:\n";  
          myfile << "   Float_t PTmuonSmear = -10.;\n";
          myfile << "   int iK_cand = -1;\n";
          myfile << "   for(int iPT = 0; iPT < NPThist; iPT++){\n";
          myfile << "   for(int iETA = 0; iETA < NETAhist; iETA++){\n";
          myfile << "      int iK = iPT + iETA*NPThist;\n";
          myfile << "      // smear muons which out of range PT and ETA\n";
          myfile << "      float PTcutDown = PTbin[iPT];\n";
          myfile << "      float PTcutUp = PTbin[iPT+1];\n";
          myfile << "      if(iPT == 0 && PTmuonGen < PTbin[0])PTcutDown = PTmuonGen-0.001;\n"; 
          myfile << "      if(iPT == (NPThist-1) && PTmuonGen > PTbin[NPThist])PTcutUp = PTmuonGen+0.001;\n"; 
          myfile << "      float ETAcutDown = ETAbin[iETA];\n";
          myfile << "      float ETAcutUp = ETAbin[iETA+1];\n";
          myfile << "      if(iETA == 0 && ETAmuonGen < ETAbin[0])ETAcutDown = ETAmuonGen-0.001;\n"; 
          myfile << "      if(iETA == (NETAhist-1) && ETAmuonGen >= ETAbin[NETAhist])ETAcutUp = ETAmuonGen+0.001;\n"; 
          myfile << "      if(PTmuonGen >= PTcutDown && PTmuonGen < PTcutUp\n";
          myfile << "      && ETAmuonGen >= ETAcutDown && ETAmuonGen < ETAcutUp) iK_cand = iK;\n";
          myfile << "   }}\n\n";

          myfile << "   ///// Set Scale Factor\n";
          myfile << "   float ScaleFactor = 1.;\n";
          myfile << "   if (fabs(ETAmuonGen) < 0.8) ScaleFactor = 1.2;\n";
          myfile << "   if (fabs(ETAmuonGen) >= 0.8 && fabs(ETAmuonGen) < 1.2) ScaleFactor = 1.15;\n";
          myfile << "   if (fabs(ETAmuonGen) >= 1.2) ScaleFactor = 1.12;\n";

          myfile << "   fitDoubleGauss->SetParameter(4,1.);\n";
          myfile << "   if(iK_cand > -1){\n";
          myfile << "      fitDoubleGauss->SetParameter(0,mean[iK_cand]);\n";
          myfile << "      fitDoubleGauss->SetParameter(1,ScaleFactor*sig1[iK_cand]);\n";
          myfile << "      fitDoubleGauss->SetParameter(2,Asig2[iK_cand]);\n";
          myfile << "      fitDoubleGauss->SetParameter(3,ScaleFactor*sig2[iK_cand]);\n";
          myfile << "      Double_t resSim = fitDoubleGauss->GetRandom();\n";
          myfile << "      PTmuonSmear = PTmuonGen*(1+resSim);\n";
          myfile << "   }\n";
          myfile << "   ////// End Smearing parametrization for single muon:\n\n";  
          myfile << "   return PTmuonSmear;\n";
          myfile << " }//end PTsmear function\n\n";

          myfile << "   Double_t DoubleGauss(Double_t *x, Double_t *par)\n";
          myfile << "   {\n";
          myfile << "       Double_t dgauss = 0.;\n";
          myfile << "        if(par[1] < par[3]){//not normalized gauss both gauss are  = 1 at x[0]=par[0]\n";  
          myfile << "                  dgauss =  exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]);\n";
          myfile << "                  dgauss = dgauss + par[2]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[3]/par[3]);\n";
          myfile << "                  dgauss = par[4]*dgauss;\n";
          myfile << "        }\n";
          myfile << "       return dgauss;\n";
          myfile << "   }\n";


          myfile.close();
        }
        else cout << "Unable to open file";

   printhistos();
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

