{
//gROOT -> ProcessLine(".L rochcor.C++"); // compile "rochcor.C" 
gROOT -> ProcessLine(".L rochcor_wasym_v3.C++"); // compile "rochcor.C"
gROOT -> ProcessLine(".L rochcor2012.C++"); // compile "rochcor.C" 
gROOT -> ProcessLine(".L FuncSmearingZmumu2012ChargePtCorr0.C++"); // compile Smearing MC function 
//gROOT -> ProcessLine(".L FuncSmearingZmumu2012PtCorr0.C++"); // compile Smearing MC function 
//gROOT -> ProcessLine(".L FuncSmearingZmumu2012PtCorr1.C++"); // compile Smearing MC function Rochester
//gROOT -> ProcessLine(".L FuncSmearingZmumu2012PtCorr2.C++"); // compile Smearing MC function MuscleFit
gROOT -> ProcessLine(".L createFuncSmearing.C++");	// compile your analysis code, "Anal.C"
createFuncSmearing pf;	// running the main function 
pf.main();
}
