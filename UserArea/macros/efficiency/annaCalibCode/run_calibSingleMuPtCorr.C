{
//gROOT -> ProcessLine(".L rochcor.C+"); // compile "rochcor.C" 
gROOT -> ProcessLine(".L rochcor_wasym_v3.C++"); // compile "rochcor.C" 
gROOT -> ProcessLine(".L rochcor2012.C++"); // compile "rochcor.C" 
gROOT -> ProcessLine(".L FuncSmearingZmumu2012PtCorr0.C++"); // compile Smearing MC function 
gROOT -> ProcessLine(".L calibSingleMuPtCorr.C++");	// compile your analysis code, "Anal.C"
calibSingleMuPtCorr pf;	// running the main function 
pf.main();
}
