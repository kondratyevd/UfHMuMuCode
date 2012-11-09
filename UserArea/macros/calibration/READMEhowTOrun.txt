******* Describtion of programms ******************

***************************************************
************ USER:                       **********
************ HOW TO use Smearing Fuction **********
***************************************************

Include in your code YourCode.C class: 
  SmearingTool.h             
  FuncSmearingZmumu2012PtCorr0.C

In root:

.L FuncSmearingZmumu2012PtCorr0.C+
.L YourCode.C+
YourCode pf;// running the make function
pf.main();

In YourCode.C:

#include <SmearingTool.h>

..........................

//Befor event loop:

// pointer to Smearing Tool:
SmearingTool *smearPT = new SmearingTool();

..........................

// inside event loop:

TLorentzVector MuTrue1; MuReco1;
float MuTrue1Charge;


..........fill it ........

/// float PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen, TString ParVar = "null",float ParSig = 0);

float ptSmear = smearPT -> PTsmear(MuTrue1.Pt(),  MuTrue1.Eta(), MuTrue1Charge);

// ParVar and ParSig use only if you want to smear ParVar with not central
// value but by ParSig sigma shift for systematic study (look in the code for
// detailes)

MuReco1.SetPtEtaPhiM(pt1Smear, reco1.eta, reco1.phi, MASS_MUON);// for example 

..........................


***************************************************
************ END: HOW TO use Smearing Fuction *****
***************************************************
********************************************************************************* 
***************************************************
*********  ADVANCE(usualy used ready Func): ******* 
*********  To Produce New Function          *******
*********  of Smear Pt of Muon by Own       *******
***************************************************
1. 
===================================================
============  createFuncSmearing.C  ===============
===================================================
in root:
.x createFuncSmearing.C++

== creat histogramms of MC resolution for Signle Muon   
== write them in ntuple

2.
===================================================
============  createFuncSmearingFit.C  ============
===================================================
in root:
.x createFuncSmearingFit.C++

== fit histogramms from createFuncSmearing.C
== produce 2 files:  
== SmearingTool.txt and FuncSmearing.txt

cp SmearingTool.txt SmearingTool.h
cp FuncSmearing.txt FuncSmearingZmumu2012PtCorr0.C // or other name
***************************************************
************ END Smear Pt of muon *****************
***************************************************

