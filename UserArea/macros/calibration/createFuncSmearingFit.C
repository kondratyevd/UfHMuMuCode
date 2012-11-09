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

  int MuCorr = 0; // 0 - no mu correction
                  // 1 - Rochester correction
                  // 2 - MuscleFit correcton in data
  TString RunYear = "2012"; // 2011A, 2011B, 2012ABCsmall, 2012
  //TString ExtraInfo = "Zmumu";
  TString ExtraInfo = "ZmumuTESTSF1";

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

void printhistos();
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



void createFuncSmearingFit(){


	gROOT->Clear();
  	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetPalette(1);
	userStyle();
	//modifiedStyle();
  	// ---- open the MC files ----
        TFile *theFile= new TFile(Form("NtupleCreateFuncSmearing_"+ExtraInfo+RunYear+"PtCorr%dGood.root", MuCorr ), "READ");
        theFile -> cd();

  for(int iPT = 0; iPT < NPThist; iPT++){
  for(int iETA = 0; iETA < NETAhist; iETA++){
      int iK = iPT + iETA*NPThist;
      hmuonRes[iK] = (TH1F*)theFile -> Get(Form("hmuonRes%d",iK));
      hmuonResNonCorr[iK] = (TH1F*)theFile -> Get(Form("hmuonResNonCorr%d",iK));
    }}


        printhistos();

        //create txt file with fit output
        ofstream myfileH ("SmearingTool.txt");
        ofstream myfile ("FuncSmearing.txt");
        if (myfileH.is_open())
        {

          myfileH << "#include <iostream>\n";
          myfileH << "#include <TChain.h>\n";
          myfileH << "#include <TClonesArray.h>\n";
          myfileH << "#include <TString.h>\n";
          myfileH << "#include <map>\n";

          myfileH << "#include <TSystem.h>\n";
          myfileH << "#include <TROOT.h>\n";
          myfileH << "#include <TMath.h>\n";
          myfileH << "#include <TLorentzVector.h>\n";
          myfileH << "#include <TRandom3.h>\n";
          myfileH << "#include \"TF1.h\"\n";
          myfileH << "#include <stdio.h>\n";
          myfileH << "#include <math.h>\n";

          myfileH << "// class init\n";
          myfileH << "class SmearingTool {\n";
          myfileH << " public:\n";
          
          myfileH << "  //SmearingTool();\n";
          myfileH << "  //SmearingTool(int seed);\n";
          myfileH << "  SmearingTool();\n";
          myfileH << "  float PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen, TString ParVar = \"null\",float ParSig = 0);\n\n"; 

          myfileH << " private:\n";
          myfileH << "  static const float mean[" << Nhist << "];\n";
          myfileH << "  static const float sig1[" << Nhist << "];\n";
          myfileH << "  static const float sig2[" << Nhist << "];\n";
          myfileH << "  static const float Asig2[" << Nhist << "];\n";

          myfileH << "  static const float ERRmean[" << Nhist << "];\n";
          myfileH << "  static const float ERRsig1[" << Nhist << "];\n";
          myfileH << "  static const float ERRsig2[" << Nhist << "];\n";
          myfileH << "  static const float ERRAsig2[" << Nhist << "];\n";

          myfileH << "  static const float ResRMS[" << Nhist << "];\n";
          myfileH << "  static const float ErrResRMS[" << Nhist << "];\n\n";

          myfileH << "  static const float PTbin[" << NPThist+1 << "];\n";
          myfileH << "  static const float ETAbin[" << NETAhist+1 << "];\n\n";
          myfileH << "};\n";
          myfileH << "// end: class init\n";

          //////////////////////////////////////////
          //////////////////////////////////////////
          myfileH.close();
        }
        else cout << "Unable to open file";


        if (myfile.is_open())
        {
          myfile << "#include <SmearingTool.h>\n\n";

          myfile << "// couldn't be inside class due to root feacher during TF1 creation:\n";
          myfile << "Double_t DoubleGauss(Double_t*, Double_t* );\n\n";

          myfile << "   const float SmearingTool::PTbin[" << NPThist+1 << "]={";

          for(int iPT = 0; iPT < NPThist+1; iPT++){
              if (iPT != (NPThist) ) myfile << PTbin[iPT] << ", ";
              else myfile << PTbin[iPT]; 
          }
          myfile << "};\n";

          myfile << "   const float SmearingTool::ETAbin[" << NETAhist+1 << "]={";

          for(int iETA = 0; iETA < NETAhist+1; iETA++){
              if (iETA != (NETAhist) ) myfile << ETAbin[iETA] << ", ";
              else myfile << ETAbin[iETA]; 
          }
          myfile << "};\n\n";

          myfile << "   ////// Smearing parametrization for single muon:\n";  
          //////////////////////////////////////////
          //////////////////////////////////////////
          myfile << "   const float SmearingTool::mean[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << mean[iK] << ", ";
              else myfile << mean[iK]; 
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::sig1[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << sig1[iK] << ", ";
              else myfile << sig1[iK];
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::sig2[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << sig2[iK] << ", ";
              else myfile << sig2[iK];
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::Asig2[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << Asig2[iK] << ", ";
              else myfile << Asig2[iK];
          }}
          myfile << "};\n\n";
          //////////////////////////////////////////
          myfile << "   const float SmearingTool::ERRmean[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRmean[iK] << ", ";
              else myfile << ERRmean[iK];
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::ERRsig1[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRsig1[iK] << ", ";
              else myfile << ERRsig1[iK];
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::ERRsig2[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRsig2[iK] << ", ";
              else myfile << ERRsig2[iK];
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::ERRAsig2[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ERRAsig2[iK] << ", ";
              else myfile << ERRAsig2[iK];
          }}
          myfile << "};\n\n";
          //////////////////////////////////////////
          myfile << "   const float SmearingTool::ResRMS[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ResRMS[iK] << ", ";
              else myfile << ResRMS[iK];
          }}
          myfile << "};\n";
          myfile << "   const float SmearingTool::ErrResRMS[" << Nhist << "]={";

          for(int iETA = 0; iETA < NETAhist; iETA++){
          for(int iPT = 0; iPT < NPThist; iPT++){
              int iK = iPT + iETA*NPThist;
              if(iPT == 0) myfile << "\n                   ";
              if (iK != (NPThist - 1 + (NETAhist-1)*NPThist) ) myfile << ErrResRMS[iK] << ", ";
              else myfile << ErrResRMS[iK];
          }}
          myfile << "};\n\n";
          //////////////////////////////////////////
          myfile << "   ////// End smearing parametrization for single muon:\n\n";  

          myfile << "SmearingTool::SmearingTool(){\n";
          myfile << "}\n";


          myfile << "float SmearingTool::PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen, TString ParVar,float ParSig){\n";
          myfile << "   const int NPThist = (sizeof(PTbin)/sizeof(float)-1);\n";
          myfile << "   const int NETAhist = (sizeof(ETAbin)/sizeof(float)-1);\n";
          myfile << "   //const int Nhist = NPThist*NETAhist;\n";
          myfile << "   ////// Make smearing for single muon from Gen level PostFSR:\n";  
          myfile << "   float PTmuonSmear = -10.;\n";
          myfile << "   float meanVar = 0;\n";
          myfile << "   float sigVar = 0;//for sig1 and sig2\n";
          myfile << "   float Asig2Var = 0;\n";
          myfile << "   if(ParVar == \"mean\")meanVar = ParSig;\n";
          myfile << "   if(ParVar == \"sig1\" || ParVar == \"sig2\")sigVar = ParSig;\n";
          myfile << "   if(ParVar == \"Asig2Var\")Asig2Var = ParSig;\n";
        
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
          myfile << "   if (fabs(ETAmuonGen) < 0.8) ScaleFactor = sqrt(1.2);\n";
          myfile << "   if (fabs(ETAmuonGen) >= 0.8 && fabs(ETAmuonGen) < 1.2) ScaleFactor = sqrt(1.15);\n";
          myfile << "   if (fabs(ETAmuonGen) >= 1.2) ScaleFactor = sqrt(1.12);\n";

          myfile << "   TF1* fitDoubleGauss = new TF1(\"fitDoubleGauss\", DoubleGauss, -0.1, 0.1, 5);\n";
          myfile << "   fitDoubleGauss->SetParameter(4,1.);\n";
          myfile << "   if(iK_cand > -1){\n";
          myfile << "      float meanCorr = mean[iK_cand] + meanVar*ERRmean[iK_cand];\n"; 
          myfile << "      float sig1Corr = sig1[iK_cand] + sigVar*ERRsig1[iK_cand];\n";
          myfile << "      if(sig1Corr < 0.005) sig1Corr = 0.005;\n"; 
          myfile << "      float sig2Corr = sig2[iK_cand] + sigVar*ERRsig2[iK_cand];\n"; 
          myfile << "      if(sig2Corr < 0.005) sig2Corr = 0.005;\n"; 
          myfile << "      float Asig2Corr = Asig2[iK_cand] + Asig2Var*ERRAsig2[iK_cand];\n"; 
          myfile << "      if(Asig2Corr < 0.) Asig2Corr = 0.;\n"; 
          myfile << "      fitDoubleGauss->SetParameter(0,meanCorr);\n";
          myfile << "      fitDoubleGauss->SetParameter(1,ScaleFactor*sig1Corr);\n";
          myfile << "      fitDoubleGauss->SetParameter(2,Asig2Corr);\n";
          myfile << "      fitDoubleGauss->SetParameter(3,ScaleFactor*sig2Corr);\n";
          myfile << "      Double_t resSim = fitDoubleGauss->GetRandom();\n";
          myfile << "      PTmuonSmear = PTmuonGen*(1+resSim);\n";
          myfile << "   }\n";
          myfile << "   delete fitDoubleGauss;\n";
          myfile << "   ////// End Smearing parametrization for single muon:\n\n";  
          myfile << "   return PTmuonSmear;\n";
          myfile << " }//end PTsmear function\n\n";

          myfile << "   Double_t DoubleGauss(Double_t *x, Double_t *par)\n";
          myfile << "   {\n";
          myfile << "       double dgauss = 0.;\n";
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
    gifname[iPT] = Form("plots/ResolutionPTScaleFactor%1.2f_"+ExtraInfo+RunYear+"PtCorr%d_"+gifmethod+"_%d", ScaleFactor, MuCorr,iPT);
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
       if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 1.15 || ETAbin[iETA] < -1.25) ) 
           || (PTbin[iPT] > 145.) ){ 
       //if (PTbin[iPT] > 145.){
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
       if( (PTbin[iPT] > 65. && ( ETAbin[iETA] > 1.15 || ETAbin[iETA] < -1.25) ) 
           || (PTbin[iPT] > 145.) ){ 
       //if (PTbin[iPT] > 145.){
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

       hmuonRes[iK] ->  GetXaxis() ->SetNdivisions(505);// n = n1 + 100*n2 + 10000*n3
       hmuonRes[iK] -> GetXaxis()->SetTitle("p_{T}^{reco}(#mu)/p_{T}^{gen}(#mu)-1");
       hmuonRes[iK] -> GetYaxis()->SetTitle("Entries");

       hmuonRes[iK] -> Draw("e");
       hmuonResNonCorr[iK] -> SetLineColor(kBlue);
       hmuonResNonCorr[iK] -> Draw("histsame");

        TLegend* histinfo = SetLegend(.6,.57,1.,.73);
        histinfo->AddEntry(hmuonRes[iK], "MC reco + Fit","lep");
        histinfo->AddEntry(fitDoubleGauss2, Form("Fit, RMS = %4.3f#pm%4.3f", ResRMS[iK], ErrResRMS[iK]),"l");
        histinfo->AddEntry(hmuonResNonCorr[iK], "MC no Corr.","lep");

       histinfo -> Draw("same");
       //c21 -> cd(iETA+1);
       //hDimuonRes[iK] -> Draw("e");
     }
     c20 ->Print(gifname[iPT]+".png");
     c20 ->Print(gifname[iPT]+".root");
   }
  /// end: print muon resolution data, MC and MC sim

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

