******* Describtion of programms ******************

***************************************************
************ USER:                       **********
************ HOW TO use Smearing Fuction **********
***************************************************

Include in your code YourCode.C class: 
  SmearingTool.h             
  FuncSmearingZmumu2012PtCorr0.C -- no muon correction in DATA
  FuncSmearingZmumu2012PtCorr1.C -- DATA after Rochester correction

In root:

.L FuncSmearingZmumu2012PtCorr0.C+
or
.L FuncSmearingZmumu2012PtCorr1.C+

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

//float PTsmear(float PTmuonGen, float ETAmuonGen, float CHARGEmuonGen, float PTmuonReco, int Ismear, TString ParVar = "null",float ParSig = 0);

float pt1Smear = smearPT -> PTsmear(MuTrue1.Pt(), MuTrue1.Eta(),MuTrue1charge, MuReco1.Pt(), Ismear);

// of cause MuTrue and MuReco should be matched before running this code
// Ismear = 1 <- Smear MC with SF and RND (no RECO use), Ismear = 2 <- Smear with SF for RECO and GEN 
// ParVar and ParSig use only if you want to smear ParVar with not central
// value but by ParSig sigma shift for systematic study (look in the code for
// detailes)

MuReco1.SetPtEtaPhiM(pt1Smear, reco1.eta, reco1.phi, MASS_MUON);// for example 

..........................

///ALSO you have to add DoubleGaussian Function:
///unfotunatly my class doesn't see it if it is not in the YourCode.C
/// it should be defined GLOBALLY:   

/// in the beginning:
Double_t DoubleGauss(Double_t*, Double_t* );

...................

// fuction 
Double_t DoubleGauss(Double_t *x, Double_t *par)
{

    //return BWnonrel(x,par) + Gauss(x,&par[3]);
    Double_t dgauss = 0.;
    //dgauss = Gauss(x,par) + Gauss(x,&par[3]);

     if(par[1] < par[3]){//not normalized gauss both gauss are  = 1 at
x[0]=par[0]  
               dgauss =  exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[1]/par[1]);
               dgauss = dgauss +
par[2]*exp(-0.5*(x[0]-par[0])*(x[0]-par[0])/par[3]/par[3]);
               dgauss = par[4]*dgauss;
     }
     //if (par[2]> 0.) dgauss = dgauss/(1+par[2]);  

    return dgauss;
}


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
.x run_createFuncSmearing.C

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

