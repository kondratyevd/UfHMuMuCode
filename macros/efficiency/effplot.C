
#include <set>
#include <time.h>

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

void effplot(bool isSave      = !true) {


  gROOT->Clear();
  //userStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptTitle(1);

TH1F* hEff = new TH1F("hEff","VBF BDT > 0 cut for 40k VBF Higgs events", 20, 0.1,0.3);
TH1F* hEffBDT = new TH1F("hEffBDT","VBF BDT > 0 for non-VBF Higgs", 20, 0.1,0.3);

TH1F* hEffGG = new TH1F("hEffGG","VBF BDT > 0 for 100K GF Higgs events", 10, 0.00,0.01);

Int_t Nsample = 10; // number of sub samples
// mva BDT cut > -0.04
//Float_t Eff_train[] = {0.181, 0.179, 0.173, 0.175, 0.171, 0.191, 0.169, 0.171, 0.179, 0.184}; // 10 subsamples for 40k events
//Float_t Eff_check[] = {0.175, 0.19, 0.185, 0.178, 0.178, 0.174, 0.181, 0.171, 0.173, 0.169}; // 10 subsamples for 40k events
// mva BDT cut > 0.0
Float_t Eff_train[] = {0.151, 0.151, 0.146, 0.146, 0.141, 0.156, 0.139, 0.148, 0.149, 0.16}; // 10 subsamples for 40k events
Float_t Eff_check[] = {0.143, 0.157, 0.155, 0.14, 0.145, 0.143, 0.145, 0.14, 0.147, 0.135}; // 10 subsamples for 40k events

Int_t NsampleGG = 20; // number of sub samples
// mva BDT cut > -0.04
//Float_t Eff_GF[] = { 0.00558, 0.00652, 0.00641, 0.00528, 0.00734, 0.00375, 0.00534, 0.00511, 0.00499, 0.00463, // trained 
//                     0.00672, 0.00725, 0.00433, 0.00528, 0.00643, 0.00569, 0.00552, 0.00527, 0.00528, 0.00366};// checked 
// mva BDT cut > 0.0
Float_t Eff_GF[] = {0.00259, 0.00415, 0.00427, 0.0033, 0.00408, 0.00154, 0.00277, 0.00184, 0.00317, 0.00397, // trained 
                    0.00359, 0.00448, 0.00216, 0.00296, 0.00402, 0.00416, 0.00297, 0.00373, 0.00396, 0.00285};// checked 

for(unsigned iS = 0; iS<Nsample ;iS++){

   hEff -> Fill(Eff_check[iS]);
   hEffBDT -> Fill(Eff_train[iS]);
}

for(unsigned iS = 0; iS<NsampleGG ;iS++){

   hEffGG -> Fill(Eff_GF[iS]);
}


hEff -> GetXaxis()->SetTitle("#epsilon for non-VBF Higgs"); 
hEff -> GetYaxis()->SetTitle("# entries"); 

hEffBDT -> SetLineColor(kRed);
hEffBDT -> SetMarkerStyle(20); 
hEffBDT -> SetMarkerColor(kRed); 

        TLegend* histinfo = new TLegend(.58,.57,.88,.73);
        histinfo->AddEntry(hEff, "cross check sample","l");
        histinfo->AddEntry(hEffBDT, "trained sample","p");


    TCanvas *c1 = new TCanvas("c1","Mass",800,600);
    c1 -> Divide(1,1);
    c1-> cd(1);
    hEff -> Draw();
    hEffBDT -> Draw("psame");
    histinfo -> Draw("same");
c1 ->Print("BDTVBFtrain_VBF.png");

        TLegend* histGG = new TLegend(.58,.57,.88,.73);
        histGG->AddEntry(hEffGG, "cross check sample","l");
    TCanvas *c2 = new TCanvas("c2","Mass",800,600);
    c2 -> Divide(1,1);
    c2-> cd(1);
    hEffGG -> Draw();
    histGG -> Draw("same");
c2 ->Print("BDTVBFtrain_GF.png");

}
