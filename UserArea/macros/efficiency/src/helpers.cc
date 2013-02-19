#include "helpers.h"
#include <cmath>
#include <iostream>
#include "TStyle.h"

float getRelIso(_MuonInfo& muon)
{
  // tracker iso cut
  //if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_2012;  

  float result = muon.sumChargedHadronPtR04 + 
    std::max(0.0,muon.sumNeutralHadronEtR04 + muon.sumPhotonEtR04 - 0.5*muon.sumPUPtR04);

  return result/muon.pt;
}

bool isKinTight_2012(_MuonInfo& muon) 
{

  bool isKinTight_2012=false;

  if (!muon.isGlobal)  return isKinTight_2012;
  if (!muon.isPFMuon) return isKinTight_2012;

  // acceptance cuts
  if (muon.pt < 25)         return isKinTight_2012; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight_2012; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_2012; // # hits in tracker

  if(getRelIso(muon) > 0.12)
      return isKinTight_2012;

  if (fabs(muon.d0_PV) > 0.2) return isKinTight_2012;
  if (fabs(muon.dz_PV) > 0.5) return isKinTight_2012;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2012;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2012;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012;
  if ( muon.normChiSquare > 10)     return isKinTight_2012;

  isKinTight_2012=true;
  return isKinTight_2012;
}

bool isKinTight_2011(_MuonInfo& muon) 
{

  bool isKinTight_2011=false;

  if (!muon.isGlobal)  return isKinTight_2011;
  if (!muon.isTracker) return isKinTight_2011;

  // acceptance cuts
  if (muon.pt < 25)         return isKinTight_2011; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight_2011; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 9) return isKinTight_2011; // # hits in tracker

  //isolation
  if(getRelIso(muon) > 0.12) return isKinTight_2011;

  if (fabs(muon.d0_PV) > 0.2) return isKinTight_2011;
  //if (fabs(muon.dz_PV) > 0.5) return isKinTight_2011;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2011;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2011;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2011;
  if ( muon.normChiSquare > 10)     return isKinTight_2011;

  isKinTight_2011=true;
  return isKinTight_2011;
}
bool isKinTight_2012_noAcc_noIso(_MuonInfo& muon)
{

  bool isKinTight_2012=false;

  if (!muon.isGlobal) return isKinTight_2012;
  if (!muon.isPFMuon) return isKinTight_2012;

  // NO acceptance cuts

  // beam spot cut
  if (fabs(muon.d0_PV) > 0.2) return isKinTight_2012;
  if (fabs(muon.dz_PV) > 0.5) return isKinTight_2012;

  // kinematic cuts
  if (muon.numTrackerLayers < 6)       return isKinTight_2012;
  if ( muon.numValidMuonHits  < 1 )    return isKinTight_2012;
  if ( muon.numValidPixelHits < 1 )    return isKinTight_2012;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012;
  if ( muon.normChiSquare > 10)        return isKinTight_2012;



  isKinTight_2012=true;
  return isKinTight_2012;
}

bool isKinTight_2011_noAcc_noIso(_MuonInfo& muon)
{

  bool isKinTight_2011=false;

  if (!muon.isGlobal)  return isKinTight_2011;
  if (!muon.isTracker) return isKinTight_2011;

  // No acceptance cuts

  // kinematic cuts
  if (muon.numTrackerLayers < 9) return isKinTight_2011; // # hits in tracker


  if (fabs(muon.d0_PV) > 0.2) return isKinTight_2011;
  //if (fabs(muon.dz_PV) > 0.5) return isKinTight_2011;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2011;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2011;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2011;
  if ( muon.normChiSquare > 10)     return isKinTight_2011;

  isKinTight_2011=true;
  return isKinTight_2011;
}


bool passPUJetID(int flag, PUJetID desiredLevel)
{
 bool result = ( flag & (1 << desiredLevel) ) != 0;
 // std::cout << "puJetID: " << flag << " level: " << desiredLevel << " pass: "<< result << std::endl;
 return result;
}

float smearMC(float trueVal, float recoVal, float calib, float smearRatio,TRandom random, bool debug)
{
  if (trueVal > 0.0)
  {
    float uncalibratedReco = recoVal+calib;
    float dif = uncalibratedReco-trueVal;
    dif *= smearRatio;
    float result = trueVal+dif;
    if (debug)
      std::cout << "true: " << trueVal << " \t\tuncalReco: " << uncalibratedReco << " \t\trecoCand: " << recoVal << " \t\tfinal: " << result << std::endl;
    return result;
  }
  else
  {
    return recoVal+calib;
  }
}

bool isHltMatched(_MuonInfo& muon1, _MuonInfo& muon2, std::vector<int> allowedPaths)
{
  std::vector<int>::const_iterator path = allowedPaths.begin();
  std::vector<int>::const_iterator endPath = allowedPaths.end();
  for(;path != endPath;path++)
  {
    //std::cout << "muon1 index "<<*path<<": "<<muon1.isHltMatched[*path]<<std::endl;
    //std::cout << "muon2 index "<<*path<<": "<<muon2.isHltMatched[*path]<<std::endl;
    if(muon1.isHltMatched[*path]==1)
        return true;
    if(muon2.isHltMatched[*path]==1)
        return true;
  }
  return false;
}

bool isHltMatched(_MuonInfo& muon1, std::vector<int> allowedPaths)
{
  std::vector<int>::const_iterator path = allowedPaths.begin();
  std::vector<int>::const_iterator endPath = allowedPaths.end();
  for(;path != endPath;path++)
  {
    //std::cout << "muon1 index "<<*path<<": "<<muon1.isHltMatched[*path]<<std::endl;
    if(muon1.isHltMatched[*path]==1)
        return true;
  }
  return false;
}


// useful for Jet Energy Resolution
float resolutionBias(float eta)
{
  // return 0;//Nominal!
  if(eta< 1.1) return 0.05; 
  if(eta< 2.5) return 0.10; 
  if(eta< 5) return 0.30;
  return 0;
}

// from Michele for systematic uncertainty
float jerCorr(float ptold, float oldgenpt, float etaold){

  float corrpt=ptold;
  //std::cout << "ptold=" << ptold << std::endl;
  //std::cout << "oldgenpt=" << oldgenpt << std::endl;
  //std::cout << "(fabs(ptold - oldgenpt)/ oldgenpt)= " << (fabs(ptold - oldgenpt)/ oldgenpt) << std::endl;

  if (oldgenpt>15. && (fabs(ptold - oldgenpt)/ oldgenpt)<0.5) {
    if (fabs(etaold)<1.1){
      Float_t scale  = 0.05;
      Float_t deltapt = (ptold - oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      corrpt *= ptscale;
//      std::cout << " 1 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;                                    

    }
    else if (fabs(etaold)>1.1 && fabs(etaold)<2.5){
      Float_t scale  = 0.10;
      Float_t deltapt = (ptold - oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      corrpt *= ptscale;
//      std::cout << " 2 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;                                    
    } else  if (fabs(etaold)>2.5 && fabs(etaold)<5.0){
      Float_t scale  = 0.30;
      Float_t deltapt = (ptold - oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      //Float_t ptscale =  (ptold+deltapt)/ptold;
      corrpt *= ptscale;
//      std::cout << " 3 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;
    }
    //std::cout << " final ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;
  }

  return corrpt;
}



float corrPtUp(float ptold, float oldgenpt, float etaold){

  float corrpt=ptold;
  if (oldgenpt>15. && (fabs(ptold - oldgenpt)/ oldgenpt)<0.5) {
    if (fabs(etaold)<1.1){
      Float_t scale  = 0.05;
      Float_t deltapt = (ptold- oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      corrpt *= ptscale;
      //std::cout << " 1 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;
    }
    else if (fabs(etaold)>1.1 && fabs(etaold)<2.5){
      Float_t scale  = 0.10;
      Float_t deltapt = (ptold- oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      corrpt *= ptscale;
      //std::cout << " 2 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl; 
    } else  if (fabs(etaold)>2.5 && fabs(etaold)<5.0){
      Float_t scale  = 0.20;
      Float_t deltapt = (ptold- oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      //Float_t ptscale =  (ptold+deltapt)/ptold;
      corrpt *= ptscale;
      std::cout << " 3 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;
    }
    //std::cout << " final ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;
  }

  return corrpt;
}

  
// from Michele for systematic uncertainty
float corrPtDown(float ptold, float oldgenpt, float etaold){

  float corrpt=ptold;
  if (oldgenpt>15. && (fabs(ptold - oldgenpt)/ oldgenpt)<0.5) {
    if (fabs(etaold)<1.1){
      Float_t scale  = -0.05;
      Float_t deltapt = (ptold- oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      corrpt *= ptscale;
      //std::cout << " 1 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;                                    

    }
    else if (fabs(etaold)>1.1 && fabs(etaold)<2.5){
      Float_t scale  = -0.10;
      Float_t deltapt = (ptold- oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      corrpt *= ptscale;
      //std::cout << " 2 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;                                    
    } else  if (fabs(etaold)>2.5 && fabs(etaold)<5.0){
      Float_t scale  = -0.20;
      Float_t deltapt = (ptold- oldgenpt)*scale;
      Float_t ptscale = TMath::Max(float(0.0),(ptold+deltapt)/ptold);
      //Float_t ptscale =  (ptold+deltapt)/ptold;                                                                                                                           
      corrpt *= ptscale;
      std::cout << " 3 ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;                                    
    }
    //std::cout << " final ptold, deltapt, ptcor, etaold" << ptold << ", " <<  deltapt  << ", " << corrpt  << ", " << etaold << std::endl;                                  
  }

  return corrpt;
}


int whichSelection(_MuonInfo& mu1, _MuonInfo& mu2,   
                   std::vector<int> allowedHLTPaths,
                   std::string& runPeriod,
                   _PFJetInfo jets,
                   bool passIncBDTCut, // keeping the complication of the 
                   bool passVBFBDTCut, // BDT out of the function
                   double  sigmasJEC,
                   double  jerUncertainty // 0 = no corr (default); +/-1 pos/neg corr.
                   )
{
  
  ////////////////////////////////////////////////////////////////////////
  // dimuon selection
  ////////////////////////////////////////////////////////////////////////

  bool pass = true;

  // muon ID
  if    (runPeriod == "7TeV") {
    if ( !isKinTight_2011(mu1) || !isKinTight_2011(mu2) ) pass = false;
  }
  else if (runPeriod == "8TeV"){
    if ( !isKinTight_2012(mu1) || !isKinTight_2012(mu2) ) pass = false;
  }
  else{
    return notSelected;
  }

  // trigger matching
  if (!isHltMatched(mu1,mu2,allowedHLTPaths)) pass = false;

  // opposite charge
  if (mu1.charge*mu2.charge != -1) pass = false;
  
  if (!pass) return notSelected;
  ////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////
  // Jet Part
  ////////////////////////////////////////////////////////////////////////
  
  bool goodJets = false;

  // JER from corrPtUp
  double jetPt0 = jerCorr(jets.pt[0],jets.genPt[0],jets.eta[0]);
  double jetPt1 = jerCorr(jets.pt[1],jets.genPt[1],jets.eta[1]);

  //jets.pt[0];//
  //jets.pt[1];//

  if (jerUncertainty > 0.) {
    jetPt0 = corrPtUp(jetPt0,jets.genPt[0],jets.eta[0]);
    jetPt1 = corrPtUp(jetPt1,jets.genPt[1],jets.eta[1]);
  }
  if (jerUncertainty < 0.) {
     jetPt0 = corrPtDown(jetPt0,jets.genPt[0],jets.eta[0]);
     jetPt1 = corrPtDown(jetPt1,jets.genPt[1],jets.eta[1]);
  }
  
//   std::cout << "jetPt0 = " << jets.pt[0] 
//             << ", jerCorr=" << jerCorr(jets.pt[0],jets.genPt[0],jets.eta[0]) 
//             << ", jetGen=" << jets.genPt[0]
//             << ", Up = " << corrPtUp(jets.pt[0],jets.genPt[0],jets.eta[0]) 
//             << ", Down = " << corrPtDown(jets.pt[0],jets.genPt[0],jets.eta[0]) << std::endl;

//   std::cout << "jetPt1 = " << jets.pt[1] 
//             << ", jerCorr=" << jerCorr(jets.pt[1],jets.genPt[1],jets.eta[1]) 
//             << ", jetGen=" << jets.genPt[1]
//             << ", Up = " << corrPtUp(jets.pt[1],jets.genPt[1],jets.eta[1]) 
//             << ", Down = " << corrPtDown(jets.pt[1],jets.genPt[1],jets.eta[1]) << std::endl;


  if (sigmasJEC != 0 && jets.jecUnc[0] > 0. && jets.jecUnc[0] < 1.) jetPt0 += (sigmasJEC*jets.jecUnc[0]*jetPt0);
  if (sigmasJEC != 0 && jets.jecUnc[1] > 0. && jets.jecUnc[1] < 1.) jetPt1 += (sigmasJEC*jets.jecUnc[1]*jetPt1);
  
  //jetPt0 = jets.pt[0] + (sigmasJEC*jets.jecUnc[0]*jets.pt[0]);//
  //jetPt1 = jets.pt[1] + (sigmasJEC*jets.jecUnc[1]*jets.pt[1]);//
  
  //if (sigmasJEC > 0) {
  //  std::cout << "jetPt0 = " << jets.pt[0] << ", jecUnc=" << jets.jecUnc[0] << ", corrected = " << jetPt0 << std::endl;
  //  std::cout << "jetPt1 = " << jets.pt[1] << ", jecUnc=" << jets.jecUnc[1] << ", corrected = " << jetPt1 << std::endl;
  //}

  if(jets.nJets>=2 && jetPt0>30.0 && jetPt1>30.0 && jets.eta[0]<4.7 && jets.eta[1]<4.7) goodJets = true;

  
  double    deltaEtaJets = -999;
  double  productEtaJets = -999;
  bool  jetInRapidityGap = false;
  int nJetsInRapidityGap = 0;
  double       diJetMass = -999;
  
  if(goodJets)
    {
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      pJet1.SetPtEtaPhiM(jets.pt[0],jets.eta[0],jets.phi[0],jets.mass[0]);
      pJet2.SetPtEtaPhiM(jets.pt[1],jets.eta[1],jets.phi[1],jets.mass[1]);
      TLorentzVector diJet = pJet1+pJet2;

      deltaEtaJets   = fabs(jets.eta[0]-jets.eta[1]);
      productEtaJets =      jets.eta[0]*jets.eta[1];

      // Seeing if there are jets in the rapidity gap
      float etaMax = jets.eta[0];
      float etaMin = 9999999.0;
      if(etaMax < jets.eta[1])
      {
          etaMax = jets.eta[1];
          etaMin = jets.eta[0];
      }
      else
      {
          etaMin = jets.eta[1];
      }

      for(unsigned iJet=2; (iJet < jets.nJets && iJet < 10);iJet++)
      {

        double jetPt = jerCorr(jets.pt[iJet],jets.genPt[iJet],jets.eta[iJet]);
        //jets.pt[iJet];//

        if (jerUncertainty > 0.) {
          jetPt = corrPtUp(jetPt,jets.genPt[iJet],jets.eta[iJet]);
        }
        if (jerUncertainty < 0.) {
          jetPt = corrPtDown(jetPt,jets.genPt[iJet],jets.eta[iJet]);
        }

        if (sigmasJEC != 0        && 
            jets.jecUnc[iJet] > 0 &&
            jets.jecUnc[iJet] < 1  ) jetPt += (sigmasJEC*jets.jecUnc[iJet]*jetPt);

        //jetPt = jets.pt[iJet] + (sigmasJEC*jets.jecUnc[iJet]*jets.pt[iJet]);//

//         if (sigmasJEC > 0) {
//           std::cout << "jetPt+ = " << jets.pt[iJet] << ", jecUnc=" << jets.jecUnc[iJet] << ", corrected = " << jetPt << std::endl;
//         }
//         if (sigmasJEC < 0) {
//           std::cout << "jetPt- = " << jets.pt[iJet] << ", jecUnc=" << jets.jecUnc[iJet] << ", corrected = " << jetPt << std::endl;
//         }


        if(jetPt > 30.0 && jets.eta[iJet] < 4.7)
        {
          if(jets.eta[iJet] < etaMax && jets.eta[iJet] > etaMin)
          {
            jetInRapidityGap = true;
            nJetsInRapidityGap++;
            
          }
        }
      }
      
      
      diJetMass = diJet.M();
    }
  
  // muon-muon categories
  bool isBB = false;
  bool isBO = false;
  bool isBE = false;
  bool isOO = false;
  bool isOE = false;
  bool isEE = false;
  if(fabs(mu1.eta)<0.8 && fabs(mu2.eta)<0.8)
    {
      isBB=true;
    }
  else if(
          (fabs(mu1.eta)<0.8 && fabs(mu2.eta)<1.6)
          || (fabs(mu1.eta)<1.6 && fabs(mu2.eta)<0.8)
          )
    {
      isBO=true;
    }
  else if(
          fabs(mu1.eta)<0.8 || fabs(mu2.eta)<0.8
          )
    {
      isBE=true;
    }
  else if(
          fabs(mu1.eta)<1.6 && fabs(mu2.eta)<1.6
          )
    {
      isOO=true;
    }
  else if(
          fabs(mu1.eta)<1.6 || fabs(mu2.eta)<1.6
          )
    {
      isOE=true;
    }
  else
    {
      isEE=true;
    }
  bool isNotBB = !isBB;


  int code = 0;

  bool isvbfPreselection = false;
  if ( diJetMass>300.0  && 
       deltaEtaJets>3.0 && 
       productEtaJets<0.0 /*&& 
                            nJetsInRapidityGap == 0 */) isvbfPreselection = true;
  
  if (isvbfPreselection) code += vbfPresel;
  else                   code += incPresel;

  

  //VBF BDT Cut Plots
  if (isvbfPreselection             && passVBFBDTCut) code += vbfPresel_passVBFBDTCut;
  if (isvbfPreselection && isBB     && passVBFBDTCut) code += vbfPresel_isBB_passVBFBDTCut;
  if (isvbfPreselection && isNotBB  && passVBFBDTCut) code += vbfPresel_isNotBB_passVBFBDTCut;
                                              
  //Inc BDT Cut Plots                                       
  if (!isvbfPreselection            && passIncBDTCut) code += incPresel_passIncBDTCut;
  if (!isvbfPreselection && isBB    && passIncBDTCut) code += incPresel_isBB_passIncBDTCut;
  if (!isvbfPreselection && isBO    && passIncBDTCut) code += incPresel_isBO_passIncBDTCut;
  if (!isvbfPreselection && isBE    && passIncBDTCut) code += incPresel_isBE_passIncBDTCut;
  if (!isvbfPreselection && isOO    && passIncBDTCut) code += incPresel_isOO_passIncBDTCut;
  if (!isvbfPreselection && isOE    && passIncBDTCut) code += incPresel_isOE_passIncBDTCut;
  if (!isvbfPreselection && isEE    && passIncBDTCut) code += incPresel_isEE_passIncBDTCut;
  if (!isvbfPreselection && isNotBB && passIncBDTCut) code += incPresel_isNotBB_passIncBDTCut;


  return code;
}

void setStyle()
{
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderSize(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasDefH(700);
  gStyle->SetCanvasDefW(700);

  gStyle->SetPadColor       (0);
  gStyle->SetPadBorderSize  (10);
  gStyle->SetPadBorderMode  (0);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin   (0.08);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadGridX       (0);
  gStyle->SetPadGridY       (0);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  gStyle->SetFrameFillStyle ( 0);
  gStyle->SetFrameFillColor ( 0);
  gStyle->SetFrameLineColor ( 1);
  gStyle->SetFrameLineStyle ( 0);
  gStyle->SetFrameLineWidth ( 1);
  gStyle->SetFrameBorderSize(10);
  gStyle->SetFrameBorderMode( 0);

  gStyle->SetNdivisions(505);

  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetMarkerStyle(20);

  gStyle->SetLabelSize(0.040,"X");
  gStyle->SetLabelSize(0.040,"Y");

  gStyle->SetLabelOffset(0.010,"X");
  gStyle->SetLabelOffset(0.010,"Y");

  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");

  gStyle->SetTitleSize(0.045,"X");
  gStyle->SetTitleSize(0.045,"Y");

  gStyle->SetTitleOffset(1.4,"X");
  gStyle->SetTitleOffset(1.4,"Y");

  gStyle->SetTextSize(0.055);
  gStyle->SetTextFont(42);

  gStyle->SetOptStat(0);
}
