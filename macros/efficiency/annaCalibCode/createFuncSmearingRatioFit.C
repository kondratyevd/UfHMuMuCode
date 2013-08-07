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

  int MuCorr = 1; // 0 - no mu correction
                  // 1 - Rochester correction
                  // 2 - MuscleFit correcton in data

  int Ismear = 2; // 0 - take pt reco from MC
                  // 1 - smear pt with own tools using pt gen post FSR  
                  // 2 - smear pt with PT_smear = PT_gen + (PT_reco-PT_gen)*SF

  TString RunYear = "2012"; // 2011A, 2011B, 2012ABCsmall, 2012
  TString ExtraInfo = "Zmumu";
  //TString ExtraInfo = "ZmumuTEST";

  //const float PTbin[] = {25., 30., 35., 40., 45., 50., 70., 100., 150., 300.}; //default
  const float PTbin[] = {25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 100., 150., 300.}; //default
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
  Double_t xScale[Nhist];
  Double_t exScale[Nhist];
  
// ---- CP error ----

TH1F* hmuonRes[Nhist];
TH1F* hmuonResNonCorr[Nhist];

// unfolding matrix

TLegend *_leg2;

Double_t DoubleGauss(Double_t*, Double_t* );
Double_t BWnonrel(Double_t*, Double_t* );
Double_t Gauss(Double_t*, Double_t* );
Double_t FuncVoigtian(Double_t*, Double_t* );
Double_t FuncVoigtianBG(Double_t*, Double_t* );

void DrawWithRatio(TCanvas *canvas, char *cTitle,
                  TH1F *hNum, TH1F *hDen);

void axis1F(TH1F  *histo,
           TAxis *xaxis,
           TAxis *yaxis,
           char  *xtitle,
           char  *ytitle);



void createFuncSmearingRatioFit(){


Double_t ResRMS_data[104]={
                   0.0210252, 0.0222618, 0.0229496, 0.0236158, 0.0241682, 0.0250542, 0.0260522, 0.0266811, 0.0276247, 0.0287379, 0.0311999, 0.0358001, 0.0460679,
                   0.0193413, 0.0203062, 0.0206666, 0.0211374, 0.0215897, 0.0220207, 0.0225298, 0.0228837, 0.0233832, 0.0241872, 0.0255785, 0.0291778, 0.0354947,
                   0.0162945, 0.017076, 0.017658, 0.0182125, 0.0187553, 0.0193422, 0.0198572, 0.0204406, 0.0209423, 0.0221594, 0.0243818, 0.0283692, 0.0344815,
                   0.0117885, 0.0124387, 0.0130438, 0.0136022, 0.0141865, 0.0147862, 0.0154244, 0.0158522, 0.0165687, 0.0174949, 0.0193388, 0.0227088, 0.0305626,
                   0.01174, 0.0124049, 0.0129916, 0.0135508, 0.0141062, 0.0147206, 0.0153064, 0.0159591, 0.0164582, 0.0175827, 0.0192714, 0.0230384, 0.0307916,
                   0.016459, 0.0172271, 0.017821, 0.0184107, 0.0189264, 0.0194864, 0.0201298, 0.0210477, 0.0214746, 0.0226302, 0.0240678, 0.0281555, 0.0349817,
                   0.0192183, 0.0201268, 0.0205685, 0.0210414, 0.0214438, 0.0219575, 0.0225097, 0.022977, 0.0233997, 0.0243226, 0.025159, 0.0285203, 0.035189,
                   0.0210412, 0.0223047, 0.0230774, 0.0237544, 0.0243305, 0.0251402, 0.0260615, 0.0272217, 0.0280872, 0.0286924, 0.0314587, 0.0355779, 0.0440838};
Double_t ErrResRMS_data[104]={
                   4.27261e-05, 3.62719e-05, 3.2461e-05, 2.9801e-05, 3.69279e-05, 6.34258e-05, 9.844e-05, 0.000138802, 0.00018875, 0.000194414, 0.00026799, 0.000462486, 0.00136084,
                   4.14054e-05, 3.46053e-05, 2.96425e-05, 2.54916e-05, 3.11152e-05, 5.19195e-05, 7.97885e-05, 0.00010987, 0.000146662, 0.000147962, 0.000198098, 0.00033101, 0.000838016,
                   3.36998e-05, 2.76506e-05, 2.30592e-05, 2.00811e-05, 2.53868e-05, 4.27269e-05, 6.46542e-05, 9.07566e-05, 0.000120641, 0.000126241, 0.000173309, 0.000292481, 0.000721503,
                   1.58402e-05, 1.23596e-05, 1.0647e-05, 1.00887e-05, 1.3131e-05, 2.21326e-05, 3.37858e-05, 4.7261e-05, 6.40956e-05, 6.6357e-05, 9.12449e-05, 0.000152549, 0.000406239,
                   1.58681e-05, 1.23332e-05, 1.06145e-05, 1.00759e-05, 1.3096e-05, 2.20771e-05, 3.36714e-05, 4.76875e-05, 6.38154e-05, 6.68607e-05, 9.1731e-05, 0.000156785, 0.000409646,
                   3.41162e-05, 2.7897e-05, 2.32854e-05, 2.03573e-05, 2.57352e-05, 4.31126e-05, 6.54879e-05, 9.3467e-05, 0.000123634, 0.000127649, 0.00017238, 0.000289204, 0.000724398,
                   4.13254e-05, 3.45499e-05, 2.97061e-05, 2.55815e-05, 3.10461e-05, 5.2393e-05, 7.95321e-05, 0.000111694, 0.000146679, 0.000150889, 0.000196315, 0.00032422, 0.00083594,
                   4.20859e-05, 3.6035e-05, 3.22826e-05, 2.96491e-05, 3.67896e-05, 6.28124e-05, 9.73557e-05, 0.000140655, 0.000189072, 0.000192285, 0.000268417, 0.000459999, 0.00133649};
//
Double_t ResRMS_Roch[104]={
                   0.0213705, 0.0226962, 0.0235394, 0.0243347, 0.025085, 0.0260635, 0.0271216, 0.0280823, 0.0289924, 0.0305302, 0.0332952, 0.0376531, 0.0455497,
                   0.0198919, 0.0210258, 0.0216035, 0.0223093, 0.0230206, 0.0236592, 0.0244445, 0.0251006, 0.0259285, 0.0272687, 0.0292082, 0.0340249, 0.0420899,
                   0.0169424, 0.0179766, 0.0187979, 0.0195693, 0.0203604, 0.0212509, 0.0220375, 0.0229798, 0.0237745, 0.025364, 0.0281083, 0.0338164, 0.0425713,
                   0.0126815, 0.0136332, 0.0145164, 0.0154004, 0.016304, 0.0172494, 0.0182192, 0.019038, 0.0201048, 0.0214842, 0.0243478, 0.0301223, 0.0398399,
                   0.0126466, 0.0136083, 0.0145067, 0.0153544, 0.0162298, 0.0171611, 0.0181081, 0.019107, 0.0200403, 0.0216637, 0.0241318, 0.0301855, 0.0410237,
                   0.0171139, 0.018116, 0.0189223, 0.0197665, 0.0205349, 0.0214153, 0.0224372, 0.0235047, 0.0244026, 0.0257633, 0.0279222, 0.0339378, 0.0431269,
                   0.019768, 0.0209026, 0.0215499, 0.0222376, 0.0228696, 0.0236667, 0.0245564, 0.025267, 0.0257635, 0.0273873, 0.0292322, 0.0342779, 0.0423192,
                   0.0214313, 0.0227803, 0.0237307, 0.0245372, 0.0252486, 0.0262597, 0.0273071, 0.0285252, 0.0295464, 0.0305564, 0.0335718, 0.0389545, 0.045871};
Double_t ErrResRMS_Roch[104]={
                   4.33694e-05, 3.69818e-05, 3.32975e-05, 3.07092e-05, 3.83271e-05, 6.59871e-05, 0.000102529, 0.00014608, 0.000198123, 0.000206625, 0.000285651, 0.000488056, 0.00136594,
                   4.25912e-05, 3.58366e-05, 3.09853e-05, 2.69059e-05, 3.31793e-05, 5.57895e-05, 8.65803e-05, 0.000120618, 0.000162658, 0.000166888, 0.000226318, 0.000386547, 0.00101252,
                   3.50425e-05, 2.91148e-05, 2.45475e-05, 2.15768e-05, 2.75596e-05, 4.6943e-05, 7.17708e-05, 0.000102021, 0.000136929, 0.00014444, 0.00020002, 0.000349759, 0.000902713,
                   1.705e-05, 1.35449e-05, 1.18481e-05, 1.1422e-05, 1.50904e-05, 2.58183e-05, 3.99065e-05, 5.67613e-05, 7.77655e-05, 8.14753e-05, 0.000114866, 0.000202688, 0.000536129,
                   1.71024e-05, 1.35281e-05, 1.18509e-05, 1.14171e-05, 1.50684e-05, 2.57382e-05, 3.98293e-05, 5.71002e-05, 7.76848e-05, 8.23841e-05, 0.000114903, 0.000205825, 0.000553667,
                   3.54823e-05, 2.93393e-05, 2.47258e-05, 2.18575e-05, 2.79232e-05, 4.73878e-05, 7.30128e-05, 0.000104409, 0.000140556, 0.000145298, 0.000199924, 0.000349596, 0.000902798,
                   4.25012e-05, 3.58837e-05, 3.11232e-05, 2.70364e-05, 3.3111e-05, 5.64773e-05, 8.67984e-05, 0.000122896, 0.000161489, 0.000169692, 0.000227918, 0.000390278, 0.00103003,
                   4.2815e-05, 3.68025e-05, 3.31953e-05, 3.06285e-05, 3.81863e-05, 6.56152e-05, 0.000102063, 0.000147547, 0.000199202, 0.000204869, 0.000286677, 0.000503825, 0.00141025};

//
Double_t ResRMS_MuScle[104]={
                   0.0230572, 0.0245172, 0.0251848, 0.025924, 0.026475, 0.0273933, 0.0284075, 0.0291166, 0.0300034, 0.0311525, 0.0338589, 0.0383194, 0.0464229,
                   0.0208961, 0.0219267, 0.0224306, 0.0229746, 0.0235315, 0.0240054, 0.0246707, 0.0252126, 0.0257935, 0.0264966, 0.028256, 0.0316068, 0.0397723,
                   0.0176111, 0.0185377, 0.0192324, 0.019857, 0.0204764, 0.0212094, 0.0217644, 0.0225642, 0.0230722, 0.0242172, 0.0265699, 0.030612, 0.0378821,
                   0.0125516, 0.0133186, 0.0140027, 0.0146559, 0.0153251, 0.0160269, 0.0167627, 0.0173163, 0.0181373, 0.0190065, 0.0210873, 0.0251277, 0.0330297,
                   0.0123303, 0.0130798, 0.0137527, 0.0143899, 0.0150248, 0.0156928, 0.0164156, 0.0171306, 0.0176599, 0.0188801, 0.0208324, 0.0249675, 0.0329407,
                   0.0176529, 0.0185427, 0.0192135, 0.019867, 0.0204566, 0.0211177, 0.0218933, 0.0228067, 0.0232545, 0.0244594, 0.0262433, 0.0308419, 0.0379268,
                   0.0212114, 0.0223181, 0.0229104, 0.023485, 0.02396, 0.0246268, 0.0253192, 0.025682, 0.0264111, 0.0273868, 0.0286062, 0.0325942, 0.0391716,
                   0.0232621, 0.0246446, 0.025523, 0.0261587, 0.0267714, 0.0275971, 0.0285945, 0.0296856, 0.0303119, 0.031332, 0.0338073, 0.0380099, 0.0463303};
Double_t ErrResRMS_MuScle[104]={
                   4.70751e-05, 3.99814e-05, 3.56414e-05, 3.27281e-05, 4.04689e-05, 6.93819e-05, 0.000107395, 0.000151604, 0.000205339, 0.000211108, 0.000291433, 0.000497529, 0.00139213,
                   4.49193e-05, 3.73884e-05, 3.21818e-05, 2.77114e-05, 3.39195e-05, 5.66102e-05, 8.73718e-05, 0.000121013, 0.000161773, 0.000162351, 0.000218783, 0.000359029, 0.000941636,
                   3.65429e-05, 3.00369e-05, 2.51224e-05, 2.18964e-05, 2.77181e-05, 4.68655e-05, 7.08829e-05, 0.000100166, 0.000133008, 0.000138017, 0.000189092, 0.000316109, 0.00079405,
                   1.69097e-05, 1.32414e-05, 1.14329e-05, 1.08718e-05, 1.4187e-05, 2.39929e-05, 3.6724e-05, 5.16264e-05, 7.01891e-05, 7.21323e-05, 9.95725e-05, 0.000168951, 0.000440985,
                   1.67074e-05, 1.3011e-05, 1.12398e-05, 1.0702e-05, 1.39515e-05, 2.35386e-05, 3.61169e-05, 5.12016e-05, 6.84861e-05, 7.18535e-05, 9.93056e-05, 0.000170095, 0.000439326,
                   3.67032e-05, 3.0042e-05, 2.51121e-05, 2.19722e-05, 2.78214e-05, 4.67277e-05, 7.12333e-05, 0.00010134, 0.000133983, 0.000138107, 0.000188309, 0.0003172, 0.000788774,
                   4.58077e-05, 3.83409e-05, 3.30966e-05, 2.85608e-05, 3.46958e-05, 5.87768e-05, 8.94779e-05, 0.000124926, 0.000165581, 0.000170095, 0.000223786, 0.000370915, 0.000933187,
                   4.67372e-05, 3.98445e-05, 3.57186e-05, 3.26662e-05, 4.05029e-05, 6.89909e-05, 0.000106817, 0.000153566, 0.0002044, 0.000210306, 0.0002893, 0.000492929, 0.00141636};

Double_t Ratio_MuScle_Roch[104];
Double_t ErrRatio_MuScle_Roch[104];
Double_t Ratio_MuScle_data[104];
Double_t ErrRatio_MuScle_data[104];
Double_t Ratio_Roch_data[104];
Double_t ErrRatio_Roch_data[104];
Double_t SF_data = 1.;
Double_t SF_Roch = 1.;
Double_t SF_MuScle = 1.;


   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;

      if(ETAbin[iETA] > -0.85 && ETAbin[iETA+1] < 0.85) SF_data = 1.21;//Barrel
      if((ETAbin[iETA] > -1.25 && ETAbin[iETA+1] < -0.75) || (ETAbin[iETA] > 0.75 && ETAbin[iETA+1] < 1.25) ) SF_data = 1.16; //Overlap 
      if(ETAbin[iETA+1] < -1.15 || ETAbin[iETA] > 1.15) SF_data = 1.11;//Endcap
      ResRMS_data[iK] = SF_data*ResRMS_data[iK];
      ErrResRMS_data[iK] = SF_data*ErrResRMS_data[iK];

      if(ETAbin[iETA] > -0.85 && ETAbin[iETA+1] < 0.85) SF_Roch = 0.92;//Barrel
      if((ETAbin[iETA] > -1.25 && ETAbin[iETA+1] < -0.75) || (ETAbin[iETA] > 0.75 && ETAbin[iETA+1] < 1.25) ) SF_Roch = 0.99; //Overlap 
      if(ETAbin[iETA+1] < -1.15 || ETAbin[iETA] > 1.15) SF_Roch = 1.045;//Endcap
      ResRMS_Roch[iK] = SF_Roch*ResRMS_Roch[iK];
      ErrResRMS_Roch[iK] = SF_Roch*ErrResRMS_Roch[iK];

      if(ETAbin[iETA] > -0.85 && ETAbin[iETA+1] < 0.85) SF_MuScle = 0.975;//Barrel
      if((ETAbin[iETA] > -1.25 && ETAbin[iETA+1] < -0.75) || (ETAbin[iETA] > 0.75 && ETAbin[iETA+1] < 1.25) ) SF_MuScle = 0.99; //Overlap 
      if(ETAbin[iETA+1] < -1.15 || ETAbin[iETA] > 1.15) SF_MuScle = 1.0;//Endcap
      ResRMS_MuScle[iK] = SF_MuScle*ResRMS_MuScle[iK];
      ErrResRMS_MuScle[iK] = SF_MuScle*ErrResRMS_MuScle[iK];

      // calculate ration MuScle over Roch. Corr. 
      Ratio_MuScle_Roch[iK] = ResRMS_MuScle[iK]/ResRMS_Roch[iK];
      ErrRatio_MuScle_Roch[iK] = Ratio_MuScle_Roch[iK]*sqrt(ErrResRMS_Roch[iK]*ErrResRMS_Roch[iK]/ResRMS_Roch[iK]/ResRMS_Roch[iK]+ ErrResRMS_MuScle[iK]*ErrResRMS_MuScle[iK]/ResRMS_MuScle[iK]/ResRMS_MuScle[iK]);
      // calculate ration MuScle over no Corr. 
      Ratio_MuScle_data[iK] = ResRMS_MuScle[iK]/ResRMS_data[iK];
      ErrRatio_MuScle_data[iK] = Ratio_MuScle_data[iK]*sqrt(ErrResRMS_data[iK]*ErrResRMS_data[iK]/ResRMS_data[iK]/ResRMS_data[iK]+ ErrResRMS_MuScle[iK]*ErrResRMS_MuScle[iK]/ResRMS_MuScle[iK]/ResRMS_MuScle[iK]);
      // calculate ration Roch. over no Corr. 
      Ratio_Roch_data[iK] = ResRMS_Roch[iK]/ResRMS_data[iK];
      ErrRatio_Roch_data[iK] = Ratio_Roch_data[iK]*sqrt(ErrResRMS_data[iK]*ErrResRMS_data[iK]/ResRMS_data[iK]/ResRMS_data[iK]+ ErrResRMS_Roch[iK]*ErrResRMS_Roch[iK]/ResRMS_Roch[iK]/ResRMS_Roch[iK]);

   }}

   /// print muon resolution data, MC and MC sim
        gROOT->Clear();
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
        gStyle->SetPalette(1);
        userStyle();
   gStyle->SetOptFit(1111);
   gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
   gStyle->SetOptTitle(1); 

/// print RMS ration information:
     TCanvas *cScale1 = new TCanvas("cScale1","Resolution mass",800,600);
     cScale1 -> cd();

   TH2F *h2 = new TH2F("h2","Axes2",NPThist*NETAhist+1,0,NPThist*NETAhist+1,100,0.8,1.1);
   h2->GetXaxis()->SetNdivisions(NPThist*NETAhist);
   h2->SetTitle(Form("Ratio RMS like Data: MuScleFit/Rochester"));
   h2 -> GetYaxis()->SetTitle("Ratio RMS");
   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      int iKin = iPT + iETA*NPThist + 1;
      if (iPT == 0){
              h2->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, p_{T} vary: %4.0f-%4.0f GeV", ETAbin[iETA], ETAbin[iETA+1], PTbin[0], PTbin[NPThist]));
      }
   }}

     h2->LabelsOption("d", "X");
     h2->Draw();
   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      xScale[iK] = iK+1;
      exScale[iK] = 0.5;
      if (iPT == 0){
           float Xm = h2 -> GetBinLowEdge(iK+1);
           TLine *lineGrid = new TLine(Xm,0.8,Xm,1.1);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");
      }
   }}
   TLine *line1 = new TLine(h2->GetBinLowEdge(1),1.,h2->GetBinLowEdge(NPThist*NETAhist),1.);
          //lineGrid -> SetLineStyle(kDotted);
          line1 -> SetLineColor(kRed);
          line1 -> Draw("same");

    TGraphErrors* grScale = new TGraphErrors(NPThist*NETAhist,&xScale[0], &Ratio_MuScle_Roch[0], &exScale[0], &ErrRatio_MuScle_Roch[0]);
    grScale->SetMarkerColor(4);
    grScale->SetMarkerStyle(21);

    grScale->Draw("LP");

     cScale1 ->Print(Form("plots/ResolutionRatioRMS2012_MuScle_Roch.png"));
     cScale1 ->Print(Form("plots/ResolutionRatioRMS2012_MuScle_Roch.root"));

    delete h2;
/// end: print RMS ration information
   TH2F *h2MuScle = new TH2F("h2MuScle","Axes2",NPThist*NETAhist+1,0,NPThist*NETAhist+1,100,0.8,1.1);
   h2MuScle->GetXaxis()->SetNdivisions(NPThist*NETAhist);
   h2MuScle->SetTitle(Form("Ratio RMS like Data: MuScleFit/no Correction"));
   h2MuScle -> GetYaxis()->SetTitle("Ratio RMS");
   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      int iKin = iPT + iETA*NPThist + 1;
      if (iPT == 0){
              h2MuScle->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, p_{T} vary: %4.0f-%4.0f GeV", ETAbin[iETA], ETAbin[iETA+1], PTbin[0], PTbin[NPThist]));
      }
   }}

     h2MuScle->LabelsOption("d", "X");
     h2MuScle->Draw();
   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      xScale[iK] = iK+1;
      exScale[iK] = 0.5;
      if (iPT == 0){
           float Xm = h2MuScle -> GetBinLowEdge(iK+1);
           TLine *lineGrid = new TLine(Xm,0.8,Xm,1.1);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");
      }
   }}
   TLine *line1MuScle = new TLine(h2MuScle->GetBinLowEdge(1),1.,h2MuScle->GetBinLowEdge(NPThist*NETAhist),1.);
          //lineGrid -> SetLineStyle(kDotted);
          line1MuScle -> SetLineColor(kRed);
          line1MuScle -> Draw("same");

    TGraphErrors* grMuScle = new TGraphErrors(NPThist*NETAhist,&xScale[0], &Ratio_MuScle_data[0], &exScale[0], &ErrRatio_MuScle_data[0]);
    grMuScle->SetMarkerColor(4);
    grMuScle->SetMarkerStyle(21);

    grMuScle->Draw("LP");

     cScale1 ->Print(Form("plots/ResolutionRatioRMS2012_MuScle_noCorr.png"));
     cScale1 ->Print(Form("plots/ResolutionRatioRMS2012_MuScle_noCorr.root"));

    delete h2MuScle;
///
/// end: print RMS ration information
   TH2F *h2Roch = new TH2F("h2Roch","Axes2",NPThist*NETAhist+1,0,NPThist*NETAhist+1,100,0.8,1.1);
   h2Roch->GetXaxis()->SetNdivisions(NPThist*NETAhist);
   h2Roch->SetTitle(Form("Ratio RMS like Data: RochFit/no Correction"));
   h2Roch -> GetYaxis()->SetTitle("Ratio RMS");
   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      int iKin = iPT + iETA*NPThist + 1;
      if (iPT == 0){
              h2Roch->GetXaxis()->SetBinLabel((iKin+1), Form("#eta: %4.1f-%4.1f, p_{T} vary: %4.0f-%4.0f GeV", ETAbin[iETA], ETAbin[iETA+1], PTbin[0], PTbin[NPThist])
);
      }
   }}

     h2Roch->LabelsOption("d", "X");
     h2Roch->Draw();
   for(int iPT = 0; iPT < NPThist; iPT++){
   for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      xScale[iK] = iK+1;
      exScale[iK] = 0.5;
      if (iPT == 0){
           float Xm = h2Roch -> GetBinLowEdge(iK+1);
           TLine *lineGrid = new TLine(Xm,0.8,Xm,1.1);
           lineGrid -> SetLineStyle(kDotted);
           lineGrid -> Draw("same");
      }
   }}
   TLine *line1Roch = new TLine(h2Roch->GetBinLowEdge(1),1.,h2Roch->GetBinLowEdge(NPThist*NETAhist),1.);
          //lineGrid -> SetLineStyle(kDotted);
          line1Roch -> SetLineColor(kRed);
          line1Roch -> Draw("same");

    TGraphErrors* grRoch = new TGraphErrors(NPThist*NETAhist,&xScale[0], &Ratio_Roch_data[0], &exScale[0], &ErrRatio_Roch_data[0]);
    grRoch->SetMarkerColor(4);
    grRoch->SetMarkerStyle(21);

    grRoch->Draw("LP");

     cScale1 ->Print(Form("plots/ResolutionRatioRMS2012_Roch_noCorr.png"));
     cScale1 ->Print(Form("plots/ResolutionRatioRMS2012_Roch_noCorr.root"));

    delete h2Roch;
///

}
///////////////////////////////////////////

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

