#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/exception.hpp"
#include <boost/lexical_cast.hpp>
//Defines method of std::string that appends any type :-)
#define appendAny(a) append(boost::lexical_cast<std::string>(a))

#include <cstdio>

#define JETPUID

using namespace std;
using namespace boost;

int main(int argc, char *argv[])
{

   time_t timeStart = time(NULL);

   if (argc < 2)
   {
    cout << "A config file is required!!" << endl;
    return 1;
   }
   std::string configFileName = argv[1];
   // Trees are name "tree"

   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( configFileName.c_str() );
   outfileName.ReplaceAll(".cfg",".root");
   outfileName = TString("TMVA_") + outfileName;
   std::cout << "Outfile name: " << outfileName << endl;
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is 
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:!Color:!DrawProgressBar:AnalysisType=Classification" );

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   // Muon Only Analysis
    using namespace boost;

    ifstream infile(configFileName.c_str(),ifstream::in);
    if(!infile)
    {
      std::cout << "Error: Can't open config file...exiting." << std::endl;
      throw;
    }

    program_options::options_description optionDesc;
    optionDesc.add_options()
        ("mDiMu",program_options::value<int>(),"")
        ("ptDiMu",program_options::value<int>(),"")
        ("yDiMu",program_options::value<int>(),"")
        ("mDiJet",program_options::value<int>(),"")
        ("ptDiJet",program_options::value<int>(),"")
        ("yDiJet",program_options::value<int>(),"")
        ("ptMu1",program_options::value<int>(),"")
        ("ptMu2",program_options::value<int>(),"")
        ("etaMu1",program_options::value<int>(),"")
        ("etaMu2",program_options::value<int>(),"")
        ("ptJet1",program_options::value<int>(),"")
        ("ptJet2",program_options::value<int>(),"")
        ("etaJet1",program_options::value<int>(),"")
        ("etaJet2",program_options::value<int>(),"")
        ("cosThetaStar",program_options::value<int>(),"")
        ("cosThetaStarCS",program_options::value<int>(),"")
        ("deltaEtaJets",program_options::value<int>(),"")
        ("productEtaJets",program_options::value<int>(),"")
        ("nJetsInRapidityGap",program_options::value<int>(),"")
        
        ("deltaPhiJets",program_options::value<int>(),"")
        ("deltaRJets",program_options::value<int>(),"")
        ("deltaPhiHJ1",program_options::value<int>(),"")
        ("deltaEtaMuons",program_options::value<int>(),"")
        ("deltaPhiMuons",program_options::value<int>(),"")
        ("deltaRMuons",program_options::value<int>(),"")
    
        ("relIsoMu1",program_options::value<int>(),"")
        ("relIsoMu2",program_options::value<int>(),"")
        ("ht",program_options::value<int>(),"")
        ("nJets",program_options::value<int>(),"")
        ("htInRapidityGap",program_options::value<int>(),"")
        ("weightsDirName",program_options::value<std::string>(),"")
        ("vbfSelection",program_options::value<int>(),"")

        ("nVtx",program_options::value<int>(),"")
  
        ("sigFile",program_options::value<std::vector<std::string> >(),"")
        ("bakFile",program_options::value<std::vector<std::string> >(),"")
        ("sigWeight",program_options::value<std::vector<float> >(),"")
        ("bakWeight",program_options::value<std::vector<float> >(),"")

        ("nTrees",program_options::value<unsigned>(),"")
        ("nEventsMin",program_options::value<unsigned>(),"")

        ("mvaSignalEff",program_options::value<float>(),"")
        ("mvaCutVal",program_options::value<float>(),"")

        ("met",program_options::value<int>(),"")
#ifdef PTMISSINMVA
        ("ptmiss",program_options::value<int>(),"")
#endif
        ("puJetIDSimpleDiscJet1",program_options::value<int>(),"")
        ("puJetIDSimpleDiscJet2",program_options::value<int>(),"")
        ("puJetIDSimpleDiscJet3",program_options::value<int>(),"")
    ;

    program_options::variables_map optionMap;
    program_options::store(program_options::parse_config_file(infile, optionDesc), optionMap);
    program_options::notify(optionMap);    

    if (optionMap.count("mDiMu") && optionMap["mDiMu"].as<int>() == 1)
      factory->AddVariable("mDiMu","","",'F');
    else
      factory->AddSpectator("mDiMu","","",'F');

    if (optionMap.count("ptDiMu") && optionMap["ptDiMu"].as<int>() == 1)
      factory->AddVariable("ptDiMu","","",'F');
    else
      factory->AddSpectator("ptDiMu","","",'F');

    if (optionMap.count("yDiMu") && optionMap["yDiMu"].as<int>() == 1)
      factory->AddVariable("yDiMu","","",'F');
    else
      factory->AddSpectator("yDiMu","","",'F');

    if (optionMap.count("mDiJet") && optionMap["mDiJet"].as<int>() == 1)
      factory->AddVariable("mDiJet","","",'F');
    else
      factory->AddSpectator("mDiJet","","",'F');

    if (optionMap.count("ptDiJet") && optionMap["ptDiJet"].as<int>() == 1)
      factory->AddVariable("ptDiJet","","",'F');
    else
      factory->AddSpectator("ptDiJet","","",'F');

    if (optionMap.count("yDiJet") && optionMap["yDiJet"].as<int>() == 1)
      factory->AddVariable("yDiJet","","",'F');
    else
      factory->AddSpectator("yDiJet","","",'F');

    if (optionMap.count("ptMu1") && optionMap["ptMu1"].as<int>() == 1)
      factory->AddVariable("ptMu1","","",'F');
    else
      factory->AddSpectator("ptMu1","","",'F');

    if (optionMap.count("ptMu2") && optionMap["ptMu2"].as<int>() == 1)
      factory->AddVariable("ptMu2","","",'F');
    else
      factory->AddSpectator("ptMu2","","",'F');

    if (optionMap.count("etaMu1") && optionMap["etaMu1"].as<int>() == 1)
      factory->AddVariable("etaMu1","","",'F');
    else
      factory->AddSpectator("etaMu1","","",'F');

    if (optionMap.count("etaMu2") && optionMap["etaMu2"].as<int>() == 1)
      factory->AddVariable("etaMu2","","",'F');
    else
      factory->AddSpectator("etaMu2","","",'F');

    if (optionMap.count("ptJet1") && optionMap["ptJet1"].as<int>() == 1)
      factory->AddVariable("ptJet1","","",'F');
    else
      factory->AddSpectator("ptJet1","","",'F');

    if (optionMap.count("ptJet2") && optionMap["ptJet2"].as<int>() == 1)
      factory->AddVariable("ptJet2","","",'F');
    else
      factory->AddSpectator("ptJet2","","",'F');

    if (optionMap.count("etaJet1") && optionMap["etaJet1"].as<int>() == 1)
      factory->AddVariable("etaJet1","","",'F');
    else
      factory->AddSpectator("etaJet1","","",'F');

    if (optionMap.count("etaJet2") && optionMap["etaJet2"].as<int>() == 1)
      factory->AddVariable("etaJet2","","",'F');
    else
      factory->AddSpectator("etaJet2","","",'F');

    if (optionMap.count("cosThetaStar") && optionMap["cosThetaStar"].as<int>() == 1)
      factory->AddVariable("cosThetaStar","","",'F');
    else
      factory->AddSpectator("cosThetaStar","","",'F');

    if (optionMap.count("cosThetaStarCS") && optionMap["cosThetaStarCS"].as<int>() == 1)
      factory->AddVariable("cosThetaStarCS","","",'F');
    else
      factory->AddSpectator("cosThetaStarCS","","",'F');

    if (optionMap.count("deltaEtaJets") && optionMap["deltaEtaJets"].as<int>() == 1)
      factory->AddVariable("deltaEtaJets","","",'F');
    else
      factory->AddSpectator("deltaEtaJets","","",'F');

    if (optionMap.count("productEtaJets") && optionMap["productEtaJets"].as<int>() == 1)
      factory->AddVariable("productEtaJets","","",'F');
    else
      factory->AddSpectator("productEtaJets","","",'F');

    if (optionMap.count("nJetsInRapidityGap") && optionMap["nJetsInRapidityGap"].as<int>() == 1)
      factory->AddVariable("nJetsInRapidityGap","","",'F');
    else
      factory->AddSpectator("nJetsInRapidityGap","","",'F');

    if (optionMap.count("deltaEtaMuons") && optionMap["deltaEtaMuons"].as<int>() == 1)
      factory->AddVariable("deltaEtaMuons","","",'F');
    else
      factory->AddSpectator("deltaEtaMuons","","",'F');

    if (optionMap.count("deltaPhiMuons") && optionMap["deltaPhiMuons"].as<int>() == 1)
      factory->AddVariable("deltaPhiMuons","","",'F');
    else
      factory->AddSpectator("deltaPhiMuons","","",'F');

    if (optionMap.count("deltaRMuons") && optionMap["deltaRMuons"].as<int>() == 1)
      factory->AddVariable("deltaRMuons","","",'F');
    else
      factory->AddSpectator("deltaRMuons","","",'F');

    if (optionMap.count("deltaPhiJets") && optionMap["deltaPhiJets"].as<int>() == 1)
      factory->AddVariable("deltaPhiJets","","",'F');
    else
      factory->AddSpectator("deltaPhiJets","","",'F');

    if (optionMap.count("deltaRJets") && optionMap["deltaRJets"].as<int>() == 1)
      factory->AddVariable("deltaRJets","","",'F');
    else
      factory->AddSpectator("deltaRJets","","",'F');

    if (optionMap.count("deltaPhiHJ1") && optionMap["deltaPhiHJ1"].as<int>() == 1)
      factory->AddVariable("deltaPhiHJ1","","",'F');
    else
      factory->AddSpectator("deltaPhiHJ1","","",'F');

    if (optionMap.count("relIsoMu1") && optionMap["relIsoMu1"].as<int>() == 1)
      factory->AddVariable("relIsoMu1","","",'F');
    else
      factory->AddSpectator("relIsoMu1","","",'F');

    if (optionMap.count("relIsoMu2") && optionMap["relIsoMu2"].as<int>() == 1)
      factory->AddVariable("relIsoMu2","","",'F');
    else
      factory->AddSpectator("relIsoMu2","","",'F');

    if (optionMap.count("ht") && optionMap["ht"].as<int>() == 1)
      factory->AddVariable("ht","","",'F');
    else
      factory->AddSpectator("ht","","",'F');

    if (optionMap.count("nJets") && optionMap["nJets"].as<int>() == 1)
      factory->AddVariable("nJets","","",'F');
    else
      factory->AddSpectator("nJets","","",'F');

    if (optionMap.count("nVtx") && optionMap["nVtx"].as<int>() == 1)
      factory->AddVariable("nVtx","","",'F');
    else
      factory->AddSpectator("nVtx","","",'F');

    if (optionMap.count("htInRapidityGap") && optionMap["htInRapidityGap"].as<int>() == 1)
      factory->AddVariable("htInRapidityGap","","",'F');
    else
      factory->AddSpectator("htInRapidityGap","","",'F');

    if (optionMap.count("met") && optionMap["met"].as<int>() == 1)
      factory->AddVariable("met","","",'F');
    else
      factory->AddSpectator("met","","",'F');

#ifdef PTMISSINMVA
    if (optionMap.count("ptmiss") && optionMap["ptmiss"].as<int>() == 1)
      factory->AddVariable("ptmiss","","",'F');
    else
      factory->AddSpectator("ptmiss","","",'F');
#endif

    if (optionMap.count("puJetIDSimpleDiscJet1") && optionMap["puJetIDSimpleDiscJet1"].as<int>() == 1)
      factory->AddVariable("puJetIDSimpleDiscJet1","","",'F');
    else
      factory->AddSpectator("puJetIDSimpleDiscJet1","","",'F');

    if (optionMap.count("puJetIDSimpleDiscJet2") && optionMap["puJetIDSimpleDiscJet2"].as<int>() == 1)
      factory->AddVariable("puJetIDSimpleDiscJet2","","",'F');
    else
      factory->AddSpectator("puJetIDSimpleDiscJet2","","",'F');

    if (optionMap.count("puJetIDSimpleDiscJet3") && optionMap["puJetIDSimpleDiscJet3"].as<int>() == 1)
      factory->AddVariable("puJetIDSimpleDiscJet3","","",'F');
    else
      factory->AddSpectator("puJetIDSimpleDiscJet3","","",'F');
  
    std::string weightsDirName;
    if (optionMap.count("weightsDirName"))
        weightsDirName = optionMap["weightsDirName"].as<std::string>();
    else
    {
      std::cout << "Error: Config File: " <<configFileName
            << " does not have option 'weightsDirName', it is required, exiting." << std::endl;
      throw;
    }

    std::string cuts;
    if (optionMap.count("vbfSelection") && optionMap["vbfSelection"].as<int>() == 1)
      cuts = "mDiJet>300.0 && deltaEtaJets>3.0 && productEtaJets<0.0";
    else
      cuts = "!(mDiJet>300.0 && deltaEtaJets>3.0 && productEtaJets<0.0)";

   if (optionMap.count("sigFile")<1)
   {
      std::cout << "Error: Config File: " <<configFileName
            << " does not have any sigFiles, it is required, exiting." << std::endl;
      return 1;
   }
   if (optionMap.count("sigWeight")<1)
   {
      std::cout << "Error: Config File: " <<configFileName
            << " does not have any sigWeight, it is required, exiting." << std::endl;
      return 1;
   }
   if (optionMap.count("bakFile")<1)
   {
      std::cout << "Error: Config File: " <<configFileName
            << " does not have any bakFiles, it is required, exiting." << std::endl;
      return 1;
   }
   if (optionMap.count("bakWeight")<1)
   {
      std::cout << "Error: Config File: " <<configFileName
            << " does not have any bakFiles, it is required, exiting." << std::endl;
      return 1;
   }
   std::vector<std::string> sigFileList = optionMap["sigFile"].as<std::vector<std::string> >();
   std::vector<std::string> bakFileList = optionMap["bakFile"].as<std::vector<std::string> >();
   std::vector<float> sigWeightList = optionMap["sigWeight"].as<std::vector<float> >();
   std::vector<float> bakWeightList = optionMap["bakWeight"].as<std::vector<float> >();

   cout <<"nSigFileNames: " <<sigFileList.size() << endl;
   cout <<"nbakFileNames: " <<bakFileList.size() << endl;
   cout <<"nSigWNames: " <<sigWeightList.size() << endl;
   cout <<"nbakWNames: " <<bakWeightList.size() << endl;

   if(sigFileList.size() != sigWeightList.size())
   {
     cout << "Error: different number of signal filenames and weights" << endl;
     return 1;
   }
   if(bakFileList.size() != bakWeightList.size())
   {
     cout << "Error: different number of background filenames and weights" << endl;
     return 1;
   }

   std::vector<TFile*> tfileList;
   std::vector<TTree*> ttreeList;
   for(unsigned iFile=0;iFile<sigFileList.size();iFile++)
   {
      TFile* tmpfile = TFile::Open(sigFileList[iFile].c_str());
      float weight = sigWeightList[iFile];
      TTree* tmptree = (TTree*) tmpfile->Get("tree");
      factory->AddSignalTree(tmptree, weight);
      tfileList.push_back(tmpfile);
      ttreeList.push_back(tmptree);
   }
   for(unsigned iFile=0;iFile<bakFileList.size();iFile++)
   {
      TFile* tmpfile = TFile::Open(bakFileList[iFile].c_str());
      float weight = bakWeightList[iFile];
      TTree* tmptree = (TTree*) tmpfile->Get("tree");
      factory->AddBackgroundTree(tmptree, weight);
      tfileList.push_back(tmpfile);
      ttreeList.push_back(tmptree);
   }

    //std::string bdtOptions = "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning"
    unsigned BDTnTrees = 850;
    unsigned BDTnEventsMin = 150;
    if (optionMap.count("nTrees"))
        BDTnTrees = optionMap["nTrees"].as<unsigned>();
    if (optionMap.count("nEventsMin"))
        BDTnEventsMin = optionMap["nEventsMin"].as<unsigned>();
    std::string bdtOptions = "!H:!V:NTrees=";
    bdtOptions.appendAny(BDTnTrees);
    bdtOptions.appendAny(":nEventsMin=");
    bdtOptions.appendAny(BDTnEventsMin);
    bdtOptions.appendAny(
      ":MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning"
      );

  std::cout << "BDT option string: " << bdtOptions << std::endl;

   factory->PrepareTrainingAndTestTree( cuts.c_str(),cuts.c_str(),
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   //Event Weights
   factory->SetSignalWeightExpression    ("weight");
   factory->SetBackgroundWeightExpression("weight");

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5:CreateMVAPdfs" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                            bdtOptions.c_str());


   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:nEventsMin=150:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","GA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   time_t timePreTrain = time(NULL);
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   time_t timePreTest = time(NULL);
   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   time_t timePreEval = time(NULL);
   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   time_t timePostEval = time(NULL);

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   ////////////////////////////////////////////////////////////

   for(unsigned iFile=0;iFile<ttreeList.size();iFile++)
   {
     delete ttreeList[iFile];
   }
   for(unsigned iFile=0;iFile<tfileList.size();iFile++)
   {
     tfileList[iFile]->Close();
     delete tfileList[iFile];
   }

   ///////////////////////////////////////////////////////////

   try
   {
     filesystem::rename("weights",weightsDirName);
   }
   catch (const boost::filesystem::filesystem_error& ex)
   {
     std::cout << "Warning, rename reports Error: " << ex.what() << std::endl;
   }

   std::cout << "Times: \n";
   std::cout << "  Setup: "<<difftime(timePreTrain,timeStart)<<"\n";
   std::cout << "  Train: "<<difftime(timePreTest,timePreTrain)<<"\n";
   std::cout << "  Test: "<<difftime(timePreEval,timePreTest)<<"\n";
   std::cout << "  Eval: "<<difftime(timePostEval,timePreEval)<<"\n";
   std::cout << "  ------------------------"<<"\n";
   std::cout << "  Total: "<<difftime(time(NULL),timeStart)<<"\n";

   std::cout << "done.";

   return 0;
}
