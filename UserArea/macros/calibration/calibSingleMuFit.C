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
#include "TAxis.h"
#include<modifiedStyle.C>
#include<userStyle.C>

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

  int FitReq = 1; // 0: fit from 60 to 120 GeV; 1: fit from 80to 100 GeV:
  int MuCorr = 1; // 0 - no mu correction
                  // 1 - Rochester correction
                  // 2 - MuscleFit correcton in data
  int Ismear = 0; // 1 - make smear by own 
  //TString ExtaInfo = ""; 
  //TString ExtaInfo = "2011AV00_01_01"; 
  //TString ExtaInfo = "2011BV00_01_01"; 
  //TString ExtaInfo = "2012AV00_01_01"; 
  TString ExtaInfo = "2012ABCsmallV00_01_01"; 
  //TString ExtaInfo = "2012ABCsmallTEST"; 
  int MCDraw = 1; // 0 - don't draw MC
                  // 1 - draw MC
  //const float PTbin[] = {20., 30., 40., 50., 60., 70., 100., 150., 300.};
  //const float PTbin[] = {20., 30., 40., 45., 50., 60., 70., 100.};
  const float PTbin[] = {20., 30., 35., 40., 45., 50., 60., 70., 100.}; //default
  const float ETAbin[] = {0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
  const int NPThist = (sizeof(PTbin)/sizeof(float)-1);
  const int NETAhist = (sizeof(ETAbin)/sizeof(float)-1);
  const int Nhist = NPThist*NETAhist;

  //const float PTbin_inv[] = {20., 30., 50., 100.};
  //const float ETAbin_inv[] = {0., 0.9, 2.1};
  //const float PTbin_inv[] = {0., 20., 30., 40., 50., 70., 100., 300.};
  const float PTbin_inv[] = {25., 30., 35., 40., 45., 50., 70., 150.};
  //const float ETAbin_inv[] = {0., 0.3, 0.8, 1.2, 1.6, 2.1};
  const float ETAbin_inv[] = {-2.1, -1.6, -1.2, -0.8, 0., 0.8, 1.2, 1.6, 2.1};
  const int NPHIbin_inv = 10; // 10 bins from 0 to 2Pi 
  //const float ETAbin_tag[] = {0., 0.8, 1.2, 2.1};
  //const float ETAbin_tag[] = {0., 1., 2.1};
  const float ETAbin_tag[] = {0., 0.8, 1.2, 2.1};

  const int NPThist_inv = (sizeof(PTbin_inv)/sizeof(float)-1);
  const int NETAhist_inv = (sizeof(ETAbin_inv)/sizeof(float)-1);
  const int NETAhist_tag = (sizeof(ETAbin_tag)/sizeof(float)-1);
  const int Nhist_inv = NPThist_inv*NETAhist_inv*NETAhist_tag;
  const int Nhist_inv1 = NPThist_inv*NETAhist_inv;


  //const float ScaleFactor = 1.;
  // variation Asig2 by 3 sigma;
  //const float ScaleAsig2 = 0;// by default + 0 sigma of Asig2;
  //const float ScaleAsig2 = 3.;// + 3 sigma of Asig2;
  //const float ScaleAsig2 = -3.;// + 3 sigma of Asig2;
   //TString Extra = "";
   //TString Extra = "_AsigPlus3sigma";
   //TString Extra = "_AsigMinus3sigma";
   //TString Extra = "_MATGRAPH";

  Double_t MASS_MUON = 0.105658367;    //GeV/c2

// ---- CP error ----


TH1F* hDimuonMassGEN[Nhist_inv];
TH1F* hDimuonMass[Nhist_inv];
TH1F* hDimuonMassDATA[Nhist_inv];
TH1F* hDimuonMassTagMuPlusGEN[Nhist_inv];
TH1F* hDimuonMassTagMuPlus[Nhist_inv];
TH1F* hDimuonMassTagMuPlusDATA[Nhist_inv];

TH1F* hPhiDimuonMass[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassGEN[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassDATA[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassTagMuPlus[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassTagMuPlusGEN[NPHIbin_inv*NETAhist_inv*NETAhist_tag];
TH1F* hPhiDimuonMassTagMuPlusDATA[NPHIbin_inv*NETAhist_inv*NETAhist_tag];

TLegend *_leg2;

//////////// functions:
void fitDiMuon(TH1F* histoMC[], TH1F* histoDATA[], TH1F* histoGEN[], TString AsFunc, TString Extra);
void printhistos();
Double_t BWnonrel(Double_t*, Double_t* );
Double_t Gauss(Double_t*, Double_t* );
Double_t FuncScale(Double_t*, Double_t* );
Double_t FuncVoigtian(Double_t*, Double_t* );
Double_t FuncVoigtianBG(Double_t*, Double_t* );
Double_t DoubleGauss(Double_t*, Double_t* );

void DrawWithRatio(TCanvas *canvas, char *cTitle,
                  TH1F *hNum, TH1F *hDen);

void axis1F(TH1F  *histo,
           TAxis *xaxis,
           TAxis *yaxis,
           char  *xtitle,
           char  *ytitle);


void calibSingleMuFit(){


	gROOT->Clear();
  	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	userStyle();
	//modifiedStyle();
  	// ---- open the MC files ----
        // open file
        //TFile *theFile    =new TFile(Form("calibSingleMu%1.0f.root", ScaleFactor), "RECREATE");
        //TFile *theFile    =new TFile("calibSingleMu1_deltaR.root", "READ");
        TFile *theFile= new TFile(Form("calibSingleMu"+ExtaInfo+"PtCorr%dSmearPt%dGood.root", MuCorr, Ismear), "READ");;
        theFile -> cd();

        // read histograms from theFile

        for(int iPT = 0; iPT < NPThist_inv; iPT++){
        for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
             int iK = iPT + iETA*NPThist_inv + iETA_tag*NETAhist_inv*NPThist_inv;
             hDimuonMassGEN[iK] = (TH1F*)theFile -> Get(Form("hDimuonMassGEN%d",iK)); 
             hDimuonMass[iK] = (TH1F*)theFile -> Get(Form("hDimuonMass%d",iK)); 
             hDimuonMassDATA[iK] = (TH1F*)theFile -> Get(Form("hDimuonMassDATA%d",iK)); 
             hDimuonMassTagMuPlusGEN[iK] = (TH1F*)theFile -> Get(Form("hDimuonMassTagMuPlusGEN%d",iK)); 
             hDimuonMassTagMuPlus[iK] = (TH1F*)theFile -> Get(Form("hDimuonMassTagMuPlus%d",iK)); 
             hDimuonMassTagMuPlusDATA[iK] = (TH1F*)theFile -> Get(Form("hDimuonMassTagMuPlusDATA%d",iK)); 
        }}}
        for(int iPT = 0; iPT < NPHIbin_inv; iPT++){
        for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
             int iK = iPT + iETA*NPHIbin_inv + iETA_tag*NETAhist_inv*NPHIbin_inv;
             ////
             hPhiDimuonMassGEN[iK] = (TH1F*)theFile -> Get(Form("hPhiDimuonMassGEN%d",iK)); 
             hPhiDimuonMass[iK] = (TH1F*)theFile -> Get(Form("hPhiDimuonMass%d",iK)); 
             hPhiDimuonMassDATA[iK] = (TH1F*)theFile -> Get(Form("hPhiDimuonMassDATA%d",iK)); 
             hPhiDimuonMassTagMuPlusGEN[iK] = (TH1F*)theFile -> Get(Form("hPhiDimuonMassTagMuPlusGEN%d",iK)); 
             hPhiDimuonMassTagMuPlus[iK] = (TH1F*)theFile -> Get(Form("hPhiDimuonMassTagMuPlus%d",iK)); 
             hPhiDimuonMassTagMuPlusDATA[iK] = (TH1F*)theFile -> Get(Form("hPhiDimuonMassTagMuPlusDATA%d",iK)); 
        }}}
             ////

        printhistos();
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

 // resolution as function of PT
 TString AsFunc = "PT";
 TString Extra = "TagMuMinus";
 fitDiMuon(hDimuonMass, hDimuonMassDATA, hDimuonMassGEN, AsFunc, Extra);
 Extra = "TagMuPlus";
 fitDiMuon(hDimuonMassTagMuPlus, hDimuonMassTagMuPlusDATA, hDimuonMassTagMuPlusGEN, AsFunc, Extra);

 // resolution as function of PHI 
 AsFunc = "PHI";
 Extra = "TagMuMinus";
 fitDiMuon(hPhiDimuonMass, hPhiDimuonMassDATA, hPhiDimuonMassGEN, AsFunc, Extra);
 Extra = "TagMuPlus";
 fitDiMuon(hPhiDimuonMassTagMuPlus, hPhiDimuonMassTagMuPlusDATA, hPhiDimuonMassTagMuPlusGEN, AsFunc, Extra);

 // resolution as fucntion of PHI, MC gen leve
 AsFunc = "PHIgen";
 Extra = "TagMuMinus";
 fitDiMuon(hPhiDimuonMassGEN, hPhiDimuonMassDATA, hPhiDimuonMassGEN, AsFunc, Extra);
 Extra = "TagMuPlus";
 fitDiMuon(hPhiDimuonMassTagMuPlusGEN, hPhiDimuonMassTagMuPlusDATA, hPhiDimuonMassTagMuPlusGEN, AsFunc, Extra);

 /////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////

}
///////////////////////////////////////////


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

Double_t FuncScale(Double_t* x, Double_t* par){


   //Double_t eff = par[0]*(TMath::Erf(par[1]*x[0]-par[2]) - 1.) + par[3];
   //Non rel. Breit-Wigner                               mass     width 
   //Double_t bw_nonrel = par[0]*(TMath::BreitWigner(x[0], par[1], par[2])); 
   Double_t Scale = 0;
   if(x[0] < par[3] || x[0] >= par[6] ) Scale = par[2]; // Endcap Minus and Plus
   if(x[0] >= par[4] && x[0] < par[5] ) Scale = par[0]; // Barrel Minus and Plus
   if((x[0] >= par[3] && x[0] < par[4]) || (x[0] >= par[5] && x[0] < par[6]) ) Scale = par[1]; // Overlap Minus and Plus

   return Scale;
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
       return (par[0] * step * sum * invsq2pi / par[3] + exp(par[4]-par[5]*x[0]+0.01*par[6]*x[0]*x[0]+0.0001*par[7]*x[0]*x[0]*x[0]));
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
void fitDiMuon(TH1F* histoMC[], TH1F* histoDATA[], TH1F* histoGEN[], TString AsFunc, TString Extra){

   gStyle->SetOptFit(1111);
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1); 

   int NFunc = NPThist_inv;
   if(AsFunc == "PHI" || AsFunc == "PHIgen") NFunc = NPHIbin_inv;
   cout << " ************  NFunc = " << NFunc << endl;

  /// print Mass reco in data, MC and MC sim     
  Double_t yScale[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t xScale[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t eyScale[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t exScale[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t mass_Data[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t mass_Data_err[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t mass_MC[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t mass_MC_err[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t res_Data[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t res_Data_err[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t res_MC[NFunc*NETAhist_inv*NETAhist_tag];
  Double_t res_MC_err[NFunc*NETAhist_inv*NETAhist_tag];

  //TF1* fitBWnonrel = new TF1("fitBWnonrel", BWnonrel, -0.1, 0.1, 3);
  //TF1* fitGauss = new TF1("fitGauss", Gauss, -0.1, 0.1, 3);
  TF1* fitFuncVoigtian = new TF1("fitFuncVoigtian", FuncVoigtian, 80., 110., 4);
  //TF1* fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 60., 120., 7);
  TF1* fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 60., 120., 8);
  if(FitReq == 1)fitFuncVoigtianBG = new TF1("fitFuncVoigtianBG", FuncVoigtianBG, 80., 100., 8);
  fitFuncVoigtianBG ->SetParName(0, "Norm");
  fitFuncVoigtianBG ->SetParName(1, "mass");
  fitFuncVoigtianBG ->SetParName(2, "width");
  fitFuncVoigtianBG ->SetParName(3, "resol");
  fitFuncVoigtianBG ->SetParName(4, "BG4");
  fitFuncVoigtianBG ->SetParName(5, "BG5");
  fitFuncVoigtianBG ->SetParName(6, "10E-2BG6");
  fitFuncVoigtianBG ->SetParName(7, "10E-4BG7");
  //Double_t par[6];
  //TF1* fitDoubleGauss = new TF1("fitDoubleGauss", DoubleGauss, -0.1, 0.1, 6);

   TString gifname_invGEN[NFunc][NETAhist_tag];
   TString gifname_inv[NFunc][NETAhist_tag];
   TString gifname_invDATA[NFunc][NETAhist_tag];
   //TString gifmethodinv = Form("CristalBall_FitReq%d", FitReq);
   TString gifmethodinv = Form("Voigtian_FitReq%d", FitReq);
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1);
   for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){
   //for(int iETA_tag = 0; iETA_tag < 1; iETA_tag++){//for TEST
   for(int iPT = 0; iPT < NFunc; iPT++){
   //for(int iPT = 0; iPT < 2; iPT++){ // for TEST
     //gifname_inv[iPT] = "test.png";
     gifname_invGEN[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MuEffMassInv%deta%dpt"+Extra+"GEN_"+gifmethodinv+"_pt%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
     gifname_inv[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MassInv%deta%dpt"+Extra+"MC_"+gifmethodinv+"_pt%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
     gifname_invDATA[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MassInv%deta%dpt"+Extra+"DATA_"+gifmethodinv+"_pt%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);

     if(AsFunc == "PHI"){
       gifname_invGEN[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_EffMassInv%deta%dphi"+Extra+"GEN_"+gifmethodinv+"_phi%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
       gifname_inv[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MassInv%deta%dphi"+Extra+"MC_"+gifmethodinv+"_phi%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
       gifname_invDATA[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MassInv%deta%dphi"+Extra+"DATA_"+gifmethodinv+"_phi%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
     }
     if(AsFunc == "PHIgen"){
       gifname_invGEN[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_EffMassInv%deta%dphi"+Extra+"GENgen_"+gifmethodinv+"_phi%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
       gifname_inv[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MassInv%deta%dphi"+Extra+"MCgen_"+gifmethodinv+"_phi%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
       gifname_invDATA[iPT][iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_MassInv%deta%dphi"+Extra+"DATAgen_"+gifmethodinv+"_phi%d_etatag%d", MuCorr, Ismear, NETAhist_inv, NFunc, iPT, iETA_tag);
     }
     //cout << "test" <<endl;
     //cout << gifname_inv[iPT] <<endl;
     TCanvas *c21_gen = new TCanvas("c21_gen","Resolution mass",3000,1700);
     TCanvas *c21 = new TCanvas("c21","Resolution mass",3000,1700);
     TCanvas *c21_data = new TCanvas("c21_data","Resolution mass",3000,1700);
     //c21-> Divide(4,2);
     //c21_data-> Divide(4,2);
     //c21-> Divide(2,1);
     //c21_data-> Divide(2,1);
     c21_gen-> Divide(4,2);
     c21-> Divide(4,2);
     c21_data-> Divide(4,2);
     int Neta = NETAhist_inv;
     //if(Neta > 7) {Neta = 7; cout << "more then 8 bins, change canvas size" << endl;}
     //if(Neta > 2) {Neta = 2; cout << "more then 2 bins, change canvas size" << endl;}
     if(Neta > 8) {Neta = 8; cout << "more then 8 bins, change canvas size" << endl;}
     for(int iETA = 0; iETA < Neta; iETA++){
     //for(int iETA = 0; iETA < Neta; iETA++){// for TEST
       int iK = iPT + iETA*NFunc + iETA_tag*NETAhist_inv*NFunc;
       cout << "iETA_tag = " << iETA_tag << " iETA = " <<iETA <<  " iPT = " << iPT << endl;

       // check efficiency for pT threshould
       histoGEN[iK]->Sumw2();// numinarator

       // make histo with pt > 25 GeV (sum bins in pt):
       TH1F *hMassSum = (TH1F*)histoGEN[iETA*NFunc + iETA_tag*NETAhist_inv*NFunc] ->Clone("hMassSum");
       hMassSum->Sumw2(); // denom.: histo with no pT cut (pT > 25 GeV)
       for(int iPTsum = 1; iPTsum < NFunc; iPTsum++){
         hMassSum->Add(histoGEN[iPTsum + iETA*NFunc + iETA_tag*NETAhist_inv*NFunc]);
       } 

       TH1F *hMassEff = (TH1F*)histoGEN[iK] ->Clone("hMassEff");
       hMassEff->Sumw2();
       hMassEff-> Divide(hMassEff,hMassSum, 1., 1., "B");
       c21_gen -> cd(iETA+1);
       hMassEff -> SetLineColor(kBlue);
       hMassEff -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
       if(AsFunc == "PT")hMassEff -> GetYaxis()->SetTitle(Form("Acc. of cut %4.1f #leq p_{T}^{#mu} < %4.1f GeV",PTbin_inv[iPT], PTbin_inv[iPT+1]));
       float phi = 2./NPHIbin_inv;
       if(AsFunc == "PHI" || AsFunc == "PHIgen")hMassEff -> GetYaxis()->SetTitle(Form("Acc. of cut %4.1f#pi #leq #phi^{#mu} < %4.1f#pi ",iPT*phi, (iPT+1)*phi));
       hMassEff -> Draw("e");

       ////////////////// fit inv mass 
       fitFuncVoigtianBG->                         SetParameters(3.7E4,  91.,  2.495, 0.93, -290., -8.3, -7.8, 2.4);
       if ( (ETAbin_inv[iETA] < 1. && ETAbin_tag[iETA_tag] >1)||(ETAbin_inv[iETA] > 1. && ETAbin_tag[iETA_tag] <1)  
          ) fitFuncVoigtianBG->SetParameters(8.6E3,  91.,  2.495, 1.3, -180., -5., -4.5, 1.3);
       if (PTbin_inv[iPT] < 40) fitFuncVoigtianBG->SetParameters(1.6E+4, 90.9, 2.495, 1.3,   30., 0.8, 0.95, -3.8E-1);
       if (PTbin_inv[iPT] < 40 && 
           ((ETAbin_inv[iETA] < 1. && ETAbin_tag[iETA_tag] >1)||(ETAbin_inv[iETA] > 1. && ETAbin_tag[iETA_tag] <1))  
          ) fitFuncVoigtianBG->SetParameters(3.4E+4, 90.9, 2.495, 1.24,   8.3, 0.53, 1.2, -0.7);

       fitFuncVoigtianBG-> SetParameters(1.9e5,  91.,  2.495, 1.25, -69., -1.7, -1., 0.); 
       //fitFuncVoigtianBG->SetParameters(2000., 90., 2.495, 1., -20., -0.6, -0.0036);
       fitFuncVoigtianBG->FixParameter(2, 2.495);
       if(FitReq == 1){
           fitFuncVoigtianBG-> SetParameters(1.8e5,  91.,  2.495, 1.2, 10., 0.05, 0., 0.); 
           fitFuncVoigtianBG->SetParLimits(0, 0., 5E6);
           fitFuncVoigtianBG->SetParLimits(1, 80., 100.);
           fitFuncVoigtianBG->SetParLimits(3, 0.2, 5.);
           //fitFuncVoigtianBG->SetParLimits(4, -200., 200.);
           fitFuncVoigtianBG->SetParLimits(5, -15., 15.);//BG5 
           //fitFuncVoigtianBG->SetParLimits(6, -15., 15.);//BG6
           //fitFuncVoigtianBG->SetParLimits(7, -10., 10.);//BG7
           fitFuncVoigtianBG->FixParameter(6, 0.);
           fitFuncVoigtianBG->FixParameter(7, 0.);
       }
       if(AsFunc == "PHIgen"){
           fitFuncVoigtianBG->SetParameter(3, 0.4);
           fitFuncVoigtianBG->SetParLimits(3, 0.001, 1.);
        }
       ///////////////////////////// 
       c21 -> cd(iETA+1);
       TH1F *histoMCcorr = (TH1F*)histoMC[iK] ->Clone("histoMC");
       histoMCcorr->Sumw2();
       if(AsFunc == "PT")histoMCcorr->Divide(hMassEff); // make acc. correction only for PT
       histoMCcorr -> Fit(fitFuncVoigtianBG,"RWLM");//Fit invmass sim with Scale Facotor
       histoMCcorr -> Fit(fitFuncVoigtianBG,"RWL");//Fit invmass sim with Scale Facotor
       cout << "XXXXXXXXXXX  End MC fit XXXXXXXXXXXXXX" << endl;
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
       if(FitReq == 1)fitFuncVoigtianBGDraw = new TF1("fitFuncVoigtianBGDraw", FuncVoigtianBG, 80., 100., 8);
       fitFuncVoigtianBGDraw -> SetParameter(0, par0);
       fitFuncVoigtianBGDraw -> SetParameter(1, par1);
       fitFuncVoigtianBGDraw -> SetParameter(2, par2);
       fitFuncVoigtianBGDraw -> SetParameter(3, par3);
       fitFuncVoigtianBGDraw -> SetParameter(4, par4);
       fitFuncVoigtianBGDraw -> SetParameter(5, par5);
       fitFuncVoigtianBGDraw -> SetParameter(6, par6);
       fitFuncVoigtianBGDraw -> SetParameter(7, par7);

       //fitFuncVoigtianBG -> SetLineColor(2);//red

       histoMCcorr -> SetLineColor(kBlue);
       histoMCcorr -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
       histoMCcorr -> GetYaxis()->SetTitle("Entries");
       histoMCcorr -> Draw("hist");
       fitFuncVoigtianBGDraw -> Draw("same");
       //histoMCcorr -> Draw("esame");
       //TF1 *fpol1 = new TF1 ("fpol1", "[3]*exp([0]-[1]*x+[2]*x*x)", 60., 120.);
       TF1 *fpol1 = new TF1 ("fpol1", "exp([0]-[1]*x+0.01*[2]*x*x+0.0001*[3]*x*x*x)", 60., 120.);
       if(FitReq == 1)fpol1 = new TF1 ("fpol1", "exp([0]-[1]*x+0.01*[2]*x*x+0.0001*[3]*x*x*x)", 80., 100.);
       fpol1 -> SetParameter(0,par4);
       fpol1 -> SetParameter(1,par5);
       fpol1 -> SetParameter(2,par6);
       fpol1 -> SetParameter(3,par7);
       //fpol1 -> SetParameter(3,par0);
       fpol1 ->SetLineColor(4);
       fpol1 -> Draw("same");
        TLegend* histMCinfo = SetLegend(.56,.6,1.,.7);
        histMCinfo->AddEntry(histoMCcorr,Form("MC reco + Fit"),"l");
        histMCinfo->AddEntry(fitFuncVoigtianBGDraw, "Fit ","l");
        //histMCinfo->AddEntry(histoMCcorr, "MC reco","lep");
        histMCinfo-> Draw("same");



       /////////////////
       //fitFuncVoigtianBG->FixParameter(4, par4); //fix from MC because not enough statistics
       fitFuncVoigtianBG-> SetParameter(0,3.e4); 
       if(AsFunc == "PHIgen"){
           fitFuncVoigtianBG->SetParameter(3, 1.3);
           fitFuncVoigtianBG->SetParLimits(3, 0.2, 5.);
        }
       fitFuncVoigtianBG->FixParameter(5, par5); //fix from MC because not enough statistics
       fitFuncVoigtianBG->FixParameter(6, par6); //fix from MC because not enough statistics
       fitFuncVoigtianBG->FixParameter(7, par7); //fix from MC because not enough statistics
       c21_data -> cd(iETA+1);
       TH1F *histoDATAcorr = (TH1F*)histoDATA[iK] ->Clone("histoDATA");
       histoDATAcorr->Sumw2();
       if(AsFunc == "PT")histoDATAcorr->Divide(hMassEff); // make acc. correction only for PT
       histoDATAcorr -> Fit(fitFuncVoigtianBG,"RWLM");
       histoDATAcorr -> Fit(fitFuncVoigtianBG,"RWL");
       cout << "XXXXXXXXXXX  End DATA fit XXXXXXXXXXXXXX" << endl;
       Double_t par0_data = fitFuncVoigtianBG ->GetParameter(0);

       Double_t par1_data = fitFuncVoigtianBG ->GetParameter(1);// mass
       Double_t par1err_data = fitFuncVoigtianBG ->GetParError(1);// error mass
       Double_t par3_data = fitFuncVoigtianBG ->GetParameter(3); // resolution
       Double_t par3err_data = fitFuncVoigtianBG ->GetParError(3);// error resolution

       Double_t par4_data = fitFuncVoigtianBG ->GetParameter(4);
       Double_t par5_data = fitFuncVoigtianBG ->GetParameter(5);
       Double_t par6_data = fitFuncVoigtianBG ->GetParameter(6);
       Double_t par7_data = fitFuncVoigtianBG ->GetParameter(7);
       TF1 *fpol1_data = new TF1 ("fpol1_data", "exp([0]-[1]*x+0.01*[2]*x*x+0.0001*[3]*x*x*x)", 60., 120.);
       if(FitReq == 1)fpol1_data = new TF1 ("fpol1_data", "exp([0]-[1]*x+0.01*[2]*x*x+0.0001*[3]*x*x*x)", 80., 100.);
       fpol1_data -> SetParameter(0,par4_data);
       fpol1_data -> SetParameter(1,par5_data);
       fpol1_data -> SetParameter(2,par6_data);
       fpol1_data -> SetParameter(3,par7_data);
       //fpol1_data -> SetParameter(3,par0_data);
       fpol1_data ->SetLineColor(4);

       histoDATAcorr -> GetXaxis()->SetTitle("M_{#mu#mu} [GeV]");
       histoDATAcorr -> GetYaxis()->SetTitle("Entries");
       histoDATAcorr -> Draw("e");
       histoMCcorr -> SetLineColor(kBlue);
       if(MCDraw == 1)histoMCcorr -> DrawNormalized("samehist", histoDATAcorr->Integral());
       fpol1_data -> Draw("same");
        TLegend* histDATAinfo = SetLegend(.6,.6, 1.,.7);
        histDATAinfo->AddEntry(histoMCcorr,Form("MC reco "),"l");
        histDATAinfo->AddEntry(histoDATAcorr, "DATA reco + Fit","lep");
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
     // print plot:s
/*
     c21_gen ->Print(gifname_invGEN[iPT][iETA_tag]+".png");
     c21_gen ->Print(gifname_invGEN[iPT][iETA_tag]+".root");
     c21 ->Print(gifname_inv[iPT][iETA_tag]+".png");
     c21 ->Print(gifname_inv[iPT][iETA_tag]+".root");
     c21_data ->Print(gifname_invDATA[iPT][iETA_tag]+".png");
     c21_data ->Print(gifname_invDATA[iPT][iETA_tag]+".root");
*/
   }
   } //end iETA_tag and iPT
   ////////////////////
   for(int iETA_tag = 0; iETA_tag < NETAhist_tag; iETA_tag++){

     TF1 *fScale = new TF1 ("fScale", "[0]", 1., NFunc*NETAhist_inv );

     //TGraphErrors* grScale = new TGraphErrors(NFunc*NETAhist_inv,xScale, yScale, exScale, eyScale);
     TGraphErrors* grScale = new TGraphErrors(NFunc*NETAhist_inv,&xScale[0], &yScale[(iETA_tag*NETAhist_inv*NFunc)], &exScale[0], &eyScale[(iETA_tag*NETAhist_inv*NFunc)]);
     grScale->SetTitle(Form("Data/MC "));
     grScale->SetMarkerColor(4);
     grScale->SetMarkerStyle(21);

     TString gifname_ScaleFaclor[NETAhist_tag];
     TString gifname_Mass[NETAhist_tag];
     TString gifname_Res[NETAhist_tag];
     TCanvas *cScale = new TCanvas("cScale","Resolution mass",800,600);
     cScale -> cd();

     TH2F *h2 = new TH2F("h2","Axes2",NFunc*NETAhist_inv+1,0,NFunc*NETAhist_inv+1,100,0.,2.0);
     h2->GetXaxis()->SetNdivisions(NFunc*NETAhist_inv);
     //h2->GetXaxis()->LabelsOption("d");
     h2->SetTitle(Form("Scale Factor Data/MC, Tag #mu^{-} #eta: %4.1f-%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
     if(Extra == "TagMuPlus")
        h2->SetTitle(Form("Scale Factor Data/MC, Tag #mu^{+} #eta: %4.1f-%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
     if(AsFunc == "PT"){
        gifname_ScaleFaclor[iETA_tag]= Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ScaleFactor%deta%dpt"+Extra+"DATA_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
        gifname_Mass[iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ZMassPeak%deta%dpt"+Extra+"DATA_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
        gifname_Res[iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ZResPeak%deta%dpt"+Extra+"DATA_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
     }
     if(AsFunc == "PHI"){
        gifname_ScaleFaclor[iETA_tag]= Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ScaleFactor%deta%dphi"+Extra+"DATA_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
        gifname_Mass[iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ZMassPeak%deta%dphi"+Extra+"DATA_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
        gifname_Res[iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ZResPeak%deta%dphi"+Extra+"DATA_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
     }
     if(AsFunc == "PHIgen"){
        gifname_ScaleFaclor[iETA_tag]= Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ScaleFactor%deta%dphi"+Extra+"DATAMCgen_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
        gifname_Mass[iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ZMassPeak%deta%dphi"+Extra+"DATAMCgen_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
        gifname_Res[iETA_tag] = Form("plots/SingleMu_"+ExtaInfo+"_MuCorr%dSmearPt%d_ZResPeak%deta%dphi"+Extra+"DATAMCgen_etatag%d_"+gifmethodinv, MuCorr, Ismear, NETAhist_inv, NFunc, iETA_tag);
     }
   for(int iPT = 0; iPT < NFunc; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NFunc + iETA_tag*NETAhist_inv*NFunc;
      int iKin = iPT + iETA*NFunc + 1; 
      if (iPT == 0){
          if(AsFunc == "PT")
              h2->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, p_{T} vary: %4.0f-%4.0f GeV", ETAbin_inv[iETA], ETAbin_inv[iETA+1], PTbin_inv[0], PTbin_inv[NFunc]));
          if(AsFunc == "PHI" || AsFunc == "PHIgen")
              h2->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, #phi vary: 0-2#pi", ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      }   
   }}

     h2->LabelsOption("d", "X");
     h2->Draw();
     float Xmin = h2->GetBinLowEdge(1);
     float Xmax = h2->GetBinLowEdge(1+NFunc*NETAhist_inv);
     float XEndPl, XEndM, XBarPl, XBarM; 
     for(int iPT = 0; iPT < NFunc; iPT++){
     for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        int iK = iPT + iETA*NFunc+1;
        if (iPT == 0){
           if(ETAbin_inv[iETA] < -1.15 && ETAbin_inv[iETA] > - 1.25) XEndM = h2 -> GetBinLowEdge(iK+1);
           if(ETAbin_inv[iETA] > 1.15 && ETAbin_inv[iETA] <  1.25) XEndPl = h2 -> GetBinLowEdge(iK+1);
           if(ETAbin_inv[iETA] < -0.75 && ETAbin_inv[iETA] > - 0.85) XBarM = h2 -> GetBinLowEdge(iK+1);
           if(ETAbin_inv[iETA] > 0.75 && ETAbin_inv[iETA] < 0.85) XBarPl = h2 -> GetBinLowEdge(iK+1);
           float Xm = h2 -> GetBinLowEdge(iK+1);
           TLine *lineGrid = new TLine(Xm,0.,Xm,2.0);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");  
        }
     }}

     TF1* fitScale = new TF1("fitScale", FuncScale, Xmin, Xmax, 7);
     fitScale ->SetParName(0, "Barrel");
     fitScale ->SetParName(1, "Overlap");
     fitScale ->SetParName(2, "Endcap");
     fitScale->SetParameter(0, 1.);
     fitScale->SetParameter(1, 1.);
     fitScale->SetParameter(2, 1.);
     fitScale->FixParameter(3, XEndM);
     fitScale->FixParameter(4, XBarM);
     fitScale->FixParameter(5, XBarPl);
     fitScale->FixParameter(6, XEndPl);

     //grScale->Draw("ALP");
     grScale->Draw("LP");
     fitScale->SetLineColor(2); // red
     grScale->Fit("fitScale", "LR");
     cScale ->Print(gifname_ScaleFaclor[iETA_tag]+".png");
     cScale ->Print(gifname_ScaleFaclor[iETA_tag]+".root");
  /// end: print Mass reco in data, MC and MC sim     

//////////////////////
/// mass bihavior:

   //TGraphErrors* grMass = new TGraphErrors(NFunc*NETAhist_inv,xScale, mass_Data, exScale, mass_Data_err);
   TGraphErrors* grMass = new TGraphErrors(NFunc*NETAhist_inv,&xScale[0], &mass_Data[(iETA_tag*NETAhist_inv*NFunc)], &exScale[0], &mass_Data_err[(iETA_tag*NETAhist_inv*NFunc)]);
   grMass->SetTitle(Form("Z Mass peak"));
   grMass->SetMarkerColor(4);// blue
   grMass->SetMarkerStyle(21);
   //TGraphErrors* grMassMC = new TGraphErrors(NFunc*NETAhist_inv,xScale, mass_MC, exScale, mass_MC_err);
   TGraphErrors* grMassMC = new TGraphErrors(NFunc*NETAhist_inv,&xScale[0], &mass_MC[(iETA_tag*NETAhist_inv*NFunc)], &exScale[0], &mass_MC_err[(iETA_tag*NETAhist_inv*NFunc)]);
   grMassMC->SetTitle(Form("Z Mass peak"));
   grMassMC->SetMarkerColor(2);//red
   grMassMC->SetMarkerStyle(20);

     TCanvas *cScale1 = new TCanvas("cScale1","Resolution mass",800,600);
     cScale1 -> cd();

   //TH2F *h2Mass = new TH2F("h","Axes",Nhist_inv+1,0,Nhist_inv+1,100,89.,92.);
   //TH2F *h2Mass = new TH2F("h","Axes",NFunc*NETAhist_inv+1,0,NFunc*NETAhist_inv+1,100,90.,91.5);
   TH2F *h2Mass = new TH2F("h2Mass","Axes2Mass",NFunc*NETAhist_inv+1,0,NFunc*NETAhist_inv+1,100,89.5,92.);
        //h2Mass ->  GetXaxis() ->SetNdivisions(505);// n = n1 + 100*n2 + 10000*n3   
   h2Mass->GetXaxis()->SetNdivisions(NFunc*NETAhist_inv);
   //h2Mass->GetXaxis()->LabelsOption("d");
   h2Mass->SetTitle(Form("Z Mass peak, Tag #mu^{-} #eta: %4.1f-%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   if(Extra == "TagMuPlus")
      h2Mass->SetTitle(Form("Z Mass peak, Tag #mu^{+} #eta: %4.1f-%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   for(int iPT = 0; iPT < NFunc; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NFunc + iETA_tag*NETAhist_inv*NFunc;
      int iKin = iPT + iETA*NFunc + 1;
      //h2Mass->GetXaxis()->SetBinLabel(iK, Form("p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f-%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      if (iPT == 0){
          if(AsFunc == "PT")
              h2Mass->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, p_{T} vary: %4.0f-%4.0f GeV", ETAbin_inv[iETA], ETAbin_inv[iETA+1], PTbin_inv[0], PTbin_inv[NFunc]));
          if(AsFunc == "PHI" || AsFunc == "PHIgen")
              h2Mass->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, #phi vary: 0-2#pi", ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      }
   }}

     h2Mass->LabelsOption("d", "X");
     h2Mass->Draw();

     for(int iPT = 0; iPT < NFunc; iPT++){
     for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        int iK = iPT + iETA*NFunc+1;
        if (iPT == 0){
           float Xm = h2Mass -> GetBinLowEdge(iK+1);
           //TLine *lineGrid = new TLine(Xm,90.,Xm,91.5);
           TLine *lineGrid = new TLine(Xm,89.5,Xm,92.);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");  
        }
     }}

     //grScale->Draw("ALP");
     grMass->Draw("LP");
     fScale->SetLineColor(4); // blue
     grMass->Fit("fScale", "R");
     if (MCDraw == 1)grMassMC->Draw("LP");
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

   //TGraphErrors* grRes = new TGraphErrors(NFunc*NETAhist_inv,xScale, res_Data, exScale, res_Data_err);
   TGraphErrors* grRes = new TGraphErrors(NFunc*NETAhist_inv,&xScale[0], &res_Data[(iETA_tag*NETAhist_inv*NFunc)], &exScale[0], &res_Data_err[(iETA_tag*NETAhist_inv*NFunc)]);
   grRes->SetTitle(Form(" Resolution near Z peak"));
   grRes->SetMarkerColor(4);
   grRes->SetMarkerStyle(21);
   //TGraphErrors* grResMC = new TGraphErrors(NFunc*NETAhist_inv,xScale, res_MC, exScale, res_MC_err);
   TGraphErrors* grResMC = new TGraphErrors(NFunc*NETAhist_inv,&xScale[0], &res_MC[(iETA_tag*NETAhist_inv*NFunc)], &exScale[0], &res_MC_err[(iETA_tag*NETAhist_inv*NFunc)]);
   grResMC->SetTitle(Form(" Resolution near Z peak"));
   grResMC->SetMarkerColor(2);//red
   grResMC->SetMarkerStyle(20);

     TCanvas *cScale2 = new TCanvas("cScale2","Resolution res",800,600);
     cScale2 -> cd();

   //TH2F *h2Res = new TH2F("h","Axes",NFunc*NETAhist_inv+1,0,NFunc*NETAhist_inv+1,100,0.,3.);
   TH2F *h2Res = new TH2F("h2Res","Axes2Res",NFunc*NETAhist_inv+1,0,NFunc*NETAhist_inv+1,100,0.,3.);
   h2Res->GetXaxis()->SetNdivisions(NFunc*NETAhist_inv);
   //h2Res->GetXaxis()->LabelsOption("d");
   h2Res->SetTitle(Form("Resolution near Z peak, Tag #mu^{-} #eta: %4.1f-%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   if(Extra == "TagMuPlus")
     h2Res->SetTitle(Form("Resolution near Z peak, Tag #mu^{+} #eta: %4.1f-%4.1f", ETAbin_tag[iETA_tag], ETAbin_tag[iETA_tag+1]));
   for(int iPT = 0; iPT < NFunc; iPT++){
   for(int iETA = 0; iETA < NETAhist_inv; iETA++){
      int iK = iPT + iETA*NFunc + iETA_tag*NETAhist_inv*NFunc;
      int iKin = iPT + iETA*NFunc + 1;      
      //h2Res->GetXaxis()->SetBinLabel(iK, Form("p_{T}: %4.1f-%4.1f GeV, #eta: %4.1f-%4.1f", PTbin_inv[iPT], PTbin_inv[iPT+1], ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      if (iPT == 0){
          if(AsFunc == "PT")
              h2Res->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, p_{T} vary: %4.0f-%4.0f GeV", ETAbin_inv[iETA], ETAbin_inv[iETA+1], PTbin_inv[0], PTbin_inv[NFunc]));
          if(AsFunc == "PHI" || AsFunc == "PHIgen")
              h2Res->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, #phi vary: 0-2#pi", ETAbin_inv[iETA], ETAbin_inv[iETA+1]));
      }     
   }}

     h2Res->LabelsOption("d", "X");
     h2Res->Draw();
     for(int iPT = 0; iPT < NFunc; iPT++){
     for(int iETA = 0; iETA < NETAhist_inv; iETA++){
        int iK = iPT + iETA*NFunc+1;
        if (iPT == 0){
           float Xm = h2Res -> GetBinLowEdge(iK+1);
           //TLine *lineGrid = new TLine(Xm,0.,Xm,3.);
           TLine *lineGrid = new TLine(Xm,0.,Xm,3.);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");  
        }
     }}

     //grScale->Draw("ALP");
     grRes->Draw("LP");
     fScale->SetLineColor(4); // blue
     //grRes->Fit("fScale", "R");
     if (MCDraw == 1)grResMC->Draw("LP");
     TLegend* legRes = SetLegend(.2,.2,0.6,.3);
     legRes->AddEntry(grRes, " Data","ep");
     //legRes->AddEntry(fScale, " Fit Data","l");
     legRes->AddEntry(grResMC, " MC","ep");
     legRes -> Draw("same");
     cScale2 ->Print(gifname_Res[iETA_tag]+".png");
     cScale2 ->Print(gifname_Res[iETA_tag]+".root");
///////////////////
     //delete histo
     delete h2;
     delete h2Mass;
     delete h2Res;

  } // end iETA_tag 
}
///////////////////////////
///////////////////////////
///////////////////////////
