#include <algorithm>
#include <set>
#include <limits.h>
#include <ctime>
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include <boost/lexical_cast.hpp>
//Defines method of std::string that appends any type :-)
#define appendAny(a) append(boost::lexical_cast<std::string>(a))

#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "DataFormats.h"
#include "helpers.h"
#include "mva.h"
#include "LumiReweightingStandAlone.h"
// muon efficiencies scale factors

#include "ScaleFactors.h"

#include "annaCalibCode/SmearingTool.h"

#define JETPUID
#define PUREWEIGHT
//#define SMEARING
#define ISMEAR 1
//#define ROCHESTER
//#define MUSCLEFIT

#ifdef ROCHESTER
#include "rochester/rochcor2012.h"
#include "rochester/rochcor.h"
#endif
#ifdef MUSCLEFIT
#include "musclefit/MuScleFitCorrector.h"
#endif

using namespace std;
using namespace boost;
  
struct HistStruct
{
  HistStruct();
  ~HistStruct();
  void Write(TFile* outfile, std::string directory);
  void FillPU(const MVA& mva, 
              const double weight,
              const double weightPUErrUp,
              const double weightPUErrDown,
              const double weightPUSystShift);

  std::vector<TH1F*> histVec;

  // PU
  TH1F* mDiMu;
  TH1F* mDiMuPUErrUp;
  TH1F* mDiMuPUErrDown;
  TH1F* mDiMuPUSystShift;

  // Selection Efficiency

  // The filling of these histograms is done "manually"
  TH1F* nDiMuGen;          // Denominator of the acceptance
  TH1F* nDiMuGenInAcc;     // Denominator of the efficiency 
  TH1F* nDiMuFullSel;      // Acc x Eff

  // muon pog SF
  TH1F* nDiMuFullSelSF;    // Acc x Eff after Scale Factors from Muon POG
  TH1F* nDiMuFullSelSFUp;  // Acc x Eff after Scale Factors from Muon POG Up
  TH1F* nDiMuFullSelSFDown;// Acc x Eff after Scale Factors from Muon POG Down
    
  // JEC
  TH1F* nDiMuFullSelJECUp;  // Acc x Eff Jet Energy Corr. Up
  TH1F* nDiMuFullSelJECDown;// Acc x Eff Jet Energy Corr. Down

  // JER
  TH1F* nDiMuFullSelJERUp;  // Acc x Eff Jet Energy Corr. Up
  TH1F* nDiMuFullSelJERDown;// Acc x Eff Jet Energy Corr. Down

};

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  gErrorIgnoreLevel = kError;
  time_t timeStart = time(NULL);

  const char* optionIntro = "H->MuMu Analyzer For Systematic Uncertainties Evaluation\n\nUsage: ./systematics [--help] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names, put output name first followed by all input file names")
      ("maxEvents,m",program_options::value<int>(), "Maximum Number of Events to Process")
      ("runPeriod,r",program_options::value<string>(), "Running Periods e.g. 7TeV, 8TeV")
  ;
  
  program_options::positional_options_description optionPos;
  optionPos.add("filenames",-1);
  
  program_options::variables_map optionMap;
  program_options::store(program_options::command_line_parser(argc, argv).options(optionDesc).positional(optionPos).run(), optionMap);
  program_options::notify(optionMap);    
  
  if (optionMap.count("help")) 
  {
      cout << optionDesc << "\n";
      return 1;
  }

  std::vector<std::string> filenames;
  vector<string>::const_iterator filename;
  std::string outputFileName;
  if (optionMap.count("filenames")>0)
  {
     filenames = optionMap["filenames"].as<vector<string> >();
     if(filenames.size()<2)
     {
       cout << "Error: Need both input file and output file names, exiting." << endl;
       return 1;
     }
     outputFileName = filenames[0];
     filenames.erase(filenames.begin());
  }
  else
  {
     cout << "Error: Input file name  and ouput file name arguments required, exiting." << endl;
     return 1;
  }

  int maxEvents = std::numeric_limits<int>::max();
  if (optionMap.count("maxEvents")) 
  {
      int tmp = optionMap["maxEvents"].as<int>();
      if (tmp > 0)
      {
        maxEvents = tmp;
      }
  }
  cout << "maxEvents = "<< maxEvents << "\n";

  string runPeriod = "";
  if (optionMap.count("runPeriod"))
  {
    runPeriod = optionMap["runPeriod"].as<string>();
  }
  //else
  //{
  //  cout << "Run Period not specified, exiting"<< endl;
  //  return 1;
  //}
  cout << "Run Period: " << runPeriod << endl;

  /////////////////////////////
  //////////// Setup //////////
  /////////////////////////////

  float minMmm = 70.0;
  float maxMmm = 160.0;
  //float minMmm = 110.0;
  //float maxMmm = 150.0;

  float minBlind = 120;
  float maxBlind = 130;

  float calib = -0.1;
  float calibSysSmear = 0.2;
  float resSmear = 1.169; // should all be around 1; a ratio of resolutions
  float resSysSmear = 0.2; // Error on that ratio

  string cfgNameInc = "inclusive_";
  string cfgNameVBF = "vbf_";
  cfgNameInc.append(runPeriod);
  cfgNameVBF.append(runPeriod);
  cfgNameInc.append(".cfg");
  cfgNameVBF.append(".cfg");
  cout << "cfgNames: " << cfgNameInc <<" \t"<< cfgNameVBF << endl;

  std::vector<int> allowedHLTPaths;
  allowedHLTPaths.push_back(0); //IsoMu24

  // Check to see if it is data
  bool isData = false;
  std::vector<std::string> dataWords;
  dataWords.push_back("Run2012");
  dataWords.push_back("Run2011");
  dataWords.push_back("SingleMu");
  dataWords.push_back("DoubleMu");
  std::vector<std::string>::const_iterator dataWord;
  for(dataWord = dataWords.begin(); dataWord != dataWords.end();dataWord++)
  {
    regex re(*dataWord);
    if(regex_search(filenames[0],re))
    {
        isData = true;
        break;
    }
  }

  // Check to see if it is signal
  bool isSignal = false;
  std::vector<std::string> signalWords;
  signalWords.push_back("Hmumu");
  std::vector<std::string>::const_iterator signalWord;
  for(signalWord = signalWords.begin(); signalWord != signalWords.end();signalWord++)
  {
    regex re(*signalWord);
    if(regex_search(filenames[0],re))
    {
        isSignal = true;
        break;
    }
  }

  // Run periods
  bool is2011A = false;
  bool is2011B = false;
  bool is2012A = false;
  bool is2012B = false;
  bool is2012C = false;
  bool is2012D = false;
  static const regex re2011A("2011A");
  static const regex re2011B("2011B");
  static const regex re2012A("2012A");
  static const regex re2012B("2012B");
  static const regex re2012C("2012C");
  static const regex re2012D("2012D");
  if(regex_search(filenames[0],re2011A))
    is2011A = true;
  if(regex_search(filenames[0],re2011B))
    is2011B = true;
  if(regex_search(filenames[0],re2012A))
    is2012A = true;
  if(regex_search(filenames[0],re2012B))
    is2012B = true;
  if(regex_search(filenames[0],re2012C))
    is2012C = true;
  if(regex_search(filenames[0],re2012D))
    is2012D = true;
  if(is2011A)
    std::cout << "This is 2011A\n";
  if(is2011B)
    std::cout << "This is 2011B\n";
  if(is2012A)
    std::cout << "This is 2012A\n";
  if(is2012B)
    std::cout << "This is 2012B\n";
  if(is2012C)
    std::cout << "This is 2012C\n";
  if(is2012D)
    std::cout << "This is 2012D\n";

  if(isData)
    std::cout << "This is a Real Data Sample\n";
  else
    std::cout << "This is a MC Sample\n";
  if(isSignal)
    std::cout << "This is a Signal Sample\n";

  ///////////////////////////////
  // Which Muon Selection to Use

  bool (*muonIdFuncPtr)(_MuonInfo&);
  if (runPeriod == "7TeV")
  {
    cout << "Using 2011 Tight Muon Selection\n";
    muonIdFuncPtr = &isKinTight_2011;
  }
  else
  {
    cout << "Using 2012 Tight Muon Selection\n";
    muonIdFuncPtr = &isKinTight_2012;
  }

  ////////////
  
  TChain * tree = new TChain("tree");

  cout << "Input File Names: \n"; 
  for(filename = filenames.begin();filename != filenames.end();filename++)
  {
    cout<<"  "<< *filename << endl;
    tree->AddFile(filename->c_str());
  }

  cout << "Output File Name: " << outputFileName << endl;
  TFile * outFile = new TFile(outputFileName.c_str(),"RECREATE");

  std::string trainingTreeFileName = "";

  //////////////////////////
  // Tree Branches

  _MuonInfo reco1, reco2;

  tree->SetBranchAddress("reco1", &reco1);
  tree->SetBranchAddress("reco2", &reco2);

  float recoCandMass, recoCandPt, recoCandY, recoCandPhi;

  tree->SetBranchAddress("recoCandMass", &recoCandMass);
  tree->SetBranchAddress("recoCandPt"  , &recoCandPt );
  tree->SetBranchAddress("recoCandY"  , &recoCandY );
  tree->SetBranchAddress("recoCandPhi"  , &recoCandPhi );

  float trueMass=-99999.0;
  if(!isData && tree->GetBranchStatus("trueMass"))
    tree->SetBranchAddress("trueMass", &trueMass);

  _TrackInfo reco1GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM1HpostFSR"))
    tree->SetBranchAddress("genM1HpostFSR", &reco1GenPostFSR);

  _TrackInfo reco2GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM2HpostFSR"))
    tree->SetBranchAddress("genM2HpostFSR", &reco2GenPostFSR);

  _PFJetInfo jets;
  tree->SetBranchAddress("pfJets",&jets);

#ifdef JETPUID
  float puJetFullDisc[10];
  float puJetSimpleDisc[10];
  float puJetCutDisc[10];

  tree->SetBranchAddress("puJetFullDisc",&puJetFullDisc);
  tree->SetBranchAddress("puJetSimpleDisc",&puJetSimpleDisc);
  tree->SetBranchAddress("puJetCutDisc",&puJetCutDisc);

  int puJetFullId[10];
  int puJetSimpleId[10];
  int puJetCutId[10];

  tree->SetBranchAddress("puJetFullId",&puJetFullId);
  tree->SetBranchAddress("puJetSimpleId",&puJetSimpleId);
  tree->SetBranchAddress("puJetCutId",&puJetCutId);
#endif

  int nPU=0;
#ifdef PUREWEIGHT
  if (!isData)
  {
    tree->SetBranchAddress("nPU",&nPU);
  }
#endif
  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo",&vertexInfo);
  _EventInfo eventInfo;
  tree->SetBranchAddress("eventInfo",&eventInfo);
  _MetInfo met;
  tree->SetBranchAddress("met",&met);


  // generator level useful for acceptance calculation
  _genPartInfo genHpostFSR;
  tree->SetBranchAddress("genHpostFSR",&genHpostFSR);

  _TrackInfo   genM1HpostFSR;
  _TrackInfo   genM2HpostFSR;

  tree->SetBranchAddress("genM1HpostFSR",&genM1HpostFSR);
  tree->SetBranchAddress("genM2HpostFSR",&genM2HpostFSR);

  // for the final selection to avoid duplicates
  typedef std::pair<int,int> pairOfInt;
  typedef std::pair<double,double> pairOfDouble;
  set<pairOfDouble> uniqueGeneratedEvents;
  set<pairOfDouble> uniqueGeneratedEventsInAcc;


  //////////////////////////
  // Histograms
  std::vector< HistStruct > hsVec;
  hsVec.clear();

  HistStruct histsPUIncPresel;       hsVec.push_back (histsPUIncPresel     ); 
  HistStruct histsPUVBFPresel;       hsVec.push_back (histsPUVBFPresel     );
                                     
  HistStruct histsPUIncBDTCut;       hsVec.push_back (histsPUIncBDTCut     );   
  HistStruct histsPUIncBDTCutBB;     hsVec.push_back (histsPUIncBDTCutBB   );
  HistStruct histsPUIncBDTCutBO;     hsVec.push_back (histsPUIncBDTCutBO   );
  HistStruct histsPUIncBDTCutBE;     hsVec.push_back (histsPUIncBDTCutBE   );
  HistStruct histsPUIncBDTCutOO;     hsVec.push_back (histsPUIncBDTCutOO   );
  HistStruct histsPUIncBDTCutOE;     hsVec.push_back (histsPUIncBDTCutOE   );
  HistStruct histsPUIncBDTCutEE;     hsVec.push_back (histsPUIncBDTCutEE   );
  HistStruct histsPUIncBDTCutNotBB;  hsVec.push_back (histsPUIncBDTCutNotBB);
                                     
  HistStruct histsPUVBFBDTCut;       hsVec.push_back (histsPUVBFBDTCut     );
  HistStruct histsPUVBFBDTCutBB;     hsVec.push_back (histsPUVBFBDTCutBB   ); 
  HistStruct histsPUVBFBDTCutNotBB;  hsVec.push_back (histsPUVBFBDTCutNotBB);

  
  //////////////////////////
  //for MVA

  std::vector<std::string> mvaConfigNames;
  mvaConfigNames.push_back(cfgNameInc);
  mvaConfigNames.push_back(cfgNameVBF);
  MVA mva(mvaConfigNames,trainingTreeFileName);

  TRandom3 random(1457);
  TRandom3 randomForSF(702);

  //////////////////////////
  //for PU reweighting

#ifdef PUREWEIGHT
  //2012
  reweight::LumiReWeighting lumiWeights          ("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root",
                                                  "pileupDists/PileUpHist2012ABC.root","pileup","pileup");
  reweight::LumiReWeighting lumiWeights_ErrUp    ("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root",
                                                  "pileupDists/PileUpHist2012ABC_ErrUp.root","pileup","pileup");
  reweight::LumiReWeighting lumiWeights_ErrDown  ("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root",
                                                  "pileupDists/PileUpHist2012ABC_ErrDown.root","pileup","pileup");
  reweight::LumiReWeighting lumiWeights_SystShift("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root",
                                                  "pileupDists/PileUpHist2012ABC_SystematicShift.root","pileup","pileup");
//   if (runPeriod == "7TeV")
//   {
//     cout << "Using 2011AB PU reweighting\n";
//     lumiWeights = reweight::LumiReWeighting("pileupDists/PileUpHistMCFall11.root","pileupDists/PileUpHist2011AB.root","pileup","pileup");
//   }
//   else
//   {
//     cout << "Using 2012ABC PU reweighting\n";
//   }
#endif

  /////////////////////////////
  // Muon Momentum Corrections

#ifdef MUSCLEFIT
  std::string mfInFile;
  if (isData)
  {
    if(runPeriod == "8TeV")
      mfInFile = "musclefit/MuScleFit_2011_DATA_42X.txt";
    else
      mfInFile = "musclefit/MuScleFit_2012_DATA_53X.txt";
  }
  else
  {
    if(runPeriod == "7TeV")
      mfInFile = "musclefit/MuScleFit_2011_MC_42X.txt";
    else
      mfInFile = "musclefit/MuScleFit_2012_MC_52X.txt";
  }
  MuScleFitCorrector* mfCorr = new MuScleFitCorrector(mfInFile);
#endif
#ifdef ROCHESTER
  rochcor2012* rCorr12 = new rochcor2012();
  rochcor* rCorr11 = new rochcor();
  int rochesterRun=0;
  if(is2011B)
    rochesterRun=1;
#endif

  /////////////////////////
  // Smearing
  SmearingTool *smearPT = new SmearingTool();

  const double SQRT2 = sqrt(2);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;

  unsigned reportEach=1000;
  if (nEvents/1000>reportEach)
    reportEach = nEvents/1000;

  unsigned testCounter = 0;
  string testString;
  
  float timeReading = 0.;
  float timeReadingAll = 0.;
  float timeProcessing = 0.;
  float timeFilling = 0.;
  time_t timeStartEventLoop = time(NULL);
  for(unsigned i=0; i<nEvents;i++)
  {
    if(i >= maxEvents)
      continue;
    time_t timeStartReading = time(NULL);
    tree->GetEvent(i);
    time_t timeStopReading = time(NULL);
    timeReadingAll += difftime(timeStopReading,timeStartReading);
    if (i % reportEach == 0) cout << "Event: " << i << endl;

    TLorentzVector reco1Vec;
    TLorentzVector reco2Vec;
    reco1Vec.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,0.105);
    reco2Vec.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,0.105);
#ifdef MUSCLEFIT
    TLorentzVector reco1Cor;
    TLorentzVector reco2Cor;
    reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,0.105);
    reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,0.105);
    mfCorr->applyPtCorrection(reco1Cor,reco1.charge);
    mfCorr->applyPtCorrection(reco2Cor,reco2.charge);
    //if (!isData && runPeriod=="8TeV")
    //  mfCorr->applyPtSmearing(reco1Cor,reco1.charge);
    //  mfCorr->applyPtSmearing(reco2Cor,reco2.charge);
    TLorentzVector diMuonCor = reco1Cor + reco2Cor;
    reco1.pt = reco1Cor.Pt();
    reco2.pt = reco2Cor.Pt();
    recoCandMass = diMuonCor.M();
    recoCandPt = diMuonCor.Pt();
    recoCandY = diMuonCor.Rapidity();
    recoCandPhi = diMuonCor.Phi();
    reco1Vec = reco1Cor;
    reco2Vec = reco2Cor;
#endif
#ifdef ROCHESTER
    TLorentzVector reco1Cor;
    TLorentzVector reco2Cor;
    reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,0.105);
    reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,0.105);
    float rochesterError=1.0; //1.0 if you don't care
    if (runPeriod == "7TeV")
    {
      if (isData)
      {
        rCorr11->momcor_data(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
        rCorr11->momcor_data(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
      }
      else
      {
        rCorr11->momcor_mc(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
        rCorr11->momcor_mc(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
      }
    }
    else
    {
      if (isData)
      {
        rCorr12->momcor_data(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
        rCorr12->momcor_data(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
      }
      else
      {
        rCorr12->momcor_mc(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
        rCorr12->momcor_mc(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
      }
    }
    TLorentzVector diMuonCor = reco1Cor + reco2Cor;
    reco1.pt = reco1Cor.Pt();
    reco2.pt = reco2Cor.Pt();
    recoCandMass = diMuonCor.M();
    recoCandPt = diMuonCor.Pt();
    recoCandY = diMuonCor.Rapidity();
    recoCandPhi = diMuonCor.Phi();
    reco1Vec = reco1Cor;
    reco2Vec = reco2Cor;
#endif

    TLorentzVector recoCandVec = reco1Vec+reco2Vec;

    float mDiMuResSigUp = recoCandMass;
    float mDiMuResSigDown = recoCandMass;
    float mDiMuResASigUp = recoCandMass;
    float mDiMuResASigDown = recoCandMass;

#ifdef SMEARING
    if(isSignal)
    {
      if(reco1GenPostFSR.pt<0.)
        cout << "Muon 1 Post FSR not valid!\n";
      if(reco2GenPostFSR.pt<0.)
        cout << "Muon 2 Post FSR not valid!\n";
      float ptReco1 = smearPT -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR);
      float ptReco2 = smearPT -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR);
      TLorentzVector reco1Vec;
      TLorentzVector reco2Vec;
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      TLorentzVector diMuonVec = reco1Vec + reco2Vec;

      reco1.pt = ptReco1;
      reco2.pt = ptReco2;
      recoCandMass = diMuonVec.M();
      recoCandPt = diMuonVec.Pt();
      recoCandY = diMuonVec.Rapidity();
      recoCandPhi = diMuonVec.Phi();

      //Systematics Time
      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR, "sig1",1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"sig1",1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResSigUp = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR,"sig1",-1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"sig1",-1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResSigDown = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR,"Asig2Var",1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"Asig2Var",1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResASigUp = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR,"Asig2Var",-1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"Asig2Var",-1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResASigDown = diMuonVec.M();
      
    }
#endif

    mva.resetValues();
    mva.mDiMu = recoCandMass;
    bool inBlindWindow = mva.mDiMu < maxBlind && mva.mDiMu > minBlind;

    bool blind = false;
  
//#ifdef BLIND
//    if (inBlindWindow && isData)
//        blind = true;
//#endif

    double weight            = 1.0;
    double weightPUErrUp     = 1.0;
    double weightPUErrDown   = 1.0;
    double weightPUSystShift = 1.0;


#ifdef PUREWEIGHT
    if (!isData)
    {
      weight            = lumiWeights.weight(nPU);
      weightPUErrUp     = lumiWeights_ErrUp.weight(nPU);
      weightPUErrDown   = lumiWeights_ErrDown.weight(nPU);
      weightPUSystShift = lumiWeights_SystShift.weight(nPU);

    }
#endif

    // std::cout << "run, event = " 
    //           << eventInfo.run << ", " 
    //           << eventInfo.event << std::endl;
    
    // fill all the gen events
    // std::cout << "genHpostFSR.mass=" << genHpostFSR.mass << std::endl;
    // std::cout << "genM1HpostFSR.pt          =" << genM1HpostFSR.pt << std::endl;; 
    // std::cout << "genM2HpostFSR.pt          =" << genM2HpostFSR.pt << std::endl;;
    // std::cout << "fabs( genM1HpostFSR.eta ) =" << fabs( genM1HpostFSR.eta ) << std::endl; 
    // std::cout << "fabs( genM2HpostFSR.eta ) =" << fabs( genM2HpostFSR.eta ) << std::endl;

    pairOfInt    runEvent(eventInfo.run,eventInfo.event);
    pairOfDouble massPt(genHpostFSR.mass,genHpostFSR.pt);
    //std::cout << "uniqueGeneratedEvents=" 
    //          << uniqueGeneratedEvents.insert( runEvent ).second
    //          << std::endl;

    if (uniqueGeneratedEvents.insert( massPt ).second) {
      if (genHpostFSR.mass < 0) continue;
      //std::cout << "filling histograms\n";
      for (int iHist=0; iHist<hsVec.size(); iHist++) {
        hsVec[iHist].nDiMuGen -> Fill(0.);
      }
    }
    

    // generator level acceptance selections
    // mass gen cut 
    if (uniqueGeneratedEventsInAcc.insert( massPt ).second) {
      for (int iHist=0; iHist<hsVec.size(); iHist++) {
        
        if (genHpostFSR.mass < minMmm || genHpostFSR.mass > maxMmm)
          continue;
        
        if ( genM1HpostFSR.pt < 25 ) continue; 
        if ( genM2HpostFSR.pt < 25 ) continue;
        if ( fabs( genM1HpostFSR.eta ) > 2.1 ) continue; 
        if ( fabs( genM2HpostFSR.eta ) > 2.1 ) continue;
        
        hsVec[iHist].nDiMuGenInAcc -> Fill(0.);
      }
    }
    
    // dimuon selection
    if (!((*muonIdFuncPtr)(reco1)) || !((*muonIdFuncPtr)(reco2)))
      continue;
    
    // trigger matching
    if (!isHltMatched(reco1,reco2,allowedHLTPaths))
      continue;

    // opposite charge
    if (reco1.charge*reco2.charge != -1)
        continue;

    // mass cut
    if (mva.mDiMu < minMmm || mva.mDiMu > maxMmm)
        continue;

    _MuonInfo muon1;
    _MuonInfo muon2;
    if(reco1.pt>reco2.pt)
    {
        muon1 = reco1;
        muon2 = reco2;
    }
    else
    {
        muon1 = reco2;
        muon2 = reco1;
    }

    mva.weight = weight;
    mva.met = met.pt;
    mva.nPU = nPU;

    mva.mDiMuResSigUp = mDiMuResSigUp;
    mva.mDiMuResSigDown = mDiMuResSigDown;
    mva.mDiMuResASigUp = mDiMuResASigUp;
    mva.mDiMuResASigDown = mDiMuResASigDown;

    mva.ptMu1=muon1.pt;
    mva.ptMu2=muon2.pt;
    mva.etaMu1=muon1.eta;
    mva.etaMu2=muon2.eta;
    mva.deltaEtaMuons=fabs(muon1.eta-muon2.eta);
    mva.relIsoMu1 = getRelIso(muon1);
    mva.relIsoMu2 = getRelIso(muon2);

    mva.ptDiMu = recoCandPt;
    mva.yDiMu = recoCandY;

    float mDiMuCalibUp = mva.mDiMu+calibSysSmear;
    float mDiMuCalibDown = mva.mDiMu-calibSysSmear;

    float mDiMuResUp = smearMC(trueMass,recoCandMass,calib,resSmear+resSysSmear,random);
    float mDiMuResDown = smearMC(trueMass,recoCandMass,calib,resSmear-resSysSmear,random);

    bool isBB = false;
    bool isBO = false;
    bool isBE = false;
    bool isOO = false;
    bool isOE = false;
    bool isEE = false;
    if(fabs(muon1.eta)<0.8 && fabs(muon2.eta)<0.8)
    {
        isBB=true;
    }
    else if(
        (fabs(muon1.eta)<0.8 && fabs(muon2.eta)<1.6)
            || (fabs(muon1.eta)<1.6 && fabs(muon2.eta)<0.8)
        )
    {
        isBO=true;
    }
    else if(
        fabs(muon1.eta)<0.8 || fabs(muon2.eta)<0.8
        )
    {
        isBE=true;
    }
    else if(
        fabs(muon1.eta)<1.6 && fabs(muon2.eta)<1.6
        )
    {
        isOO=true;
    }
    else if(
        fabs(muon1.eta)<1.6 || fabs(muon2.eta)<1.6
        )
    {
        isOE=true;
    }
    else
    {
        isEE=true;
    }
    bool isNotBB = !isBB;
    

    //////////////////////////////////////////
    //Computing CosTheta*

    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    pMuon1.SetPtEtaPhiM(muon1.pt,muon1.eta,muon1.phi,0.105);
    pMuon2.SetPtEtaPhiM(muon2.pt,muon2.eta,muon2.phi,0.105);
    TLorentzVector diMuon = pMuon1+pMuon2;

    mva.deltaPhiMuons = pMuon1.DeltaPhi(pMuon2);
    mva.deltaRMuons = pMuon1.DeltaR(pMuon2);

    TLorentzVector starMuon1 = pMuon1;
    TLorentzVector starMuon2 = pMuon2;
    TVector3 boost = diMuon.BoostVector();
    starMuon1.Boost(-boost);
    starMuon2.Boost(-boost);

    //if (muon1.charge>0)
    //std::cout << "Run: " << eventInfo.run << " lumi: " << eventInfo.lumi << " event: " << eventInfo.event << std::endl;
    if ((int) (eventInfo.event) % 2 == 0)
    {
        TVector3 directionOfBoost = starMuon1.BoostVector();
        mva.cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }
    else
    {
        TVector3 directionOfBoost = starMuon2.BoostVector();
        mva.cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }

    //////////////////////////////////////////
    //Computing CosTheta* Collins-Soper

    //std::cout << "muon1 charge: " << muon1.charge << "muon2 charge: "<<muon2.charge << std::endl;
    if (muon1.charge != muon2.charge)
    {
      // p1 is lepton
      // p2 is anti-lepton
      float p1Plus=-1e15;
      float p2Plus=-1e15;
      float p1Minus=-1e15;
      float p2Minus=-1e15;
      if (muon1.charge < 0)
      {
        p1Plus  = (pMuon1.E()+pMuon1.Pz())/SQRT2;
        p1Minus = (pMuon1.E()-pMuon1.Pz())/SQRT2;
        p2Plus  = (pMuon2.E()+pMuon2.Pz())/SQRT2;
        p2Minus = (pMuon2.E()-pMuon2.Pz())/SQRT2;
      }
      else
      {
        p1Plus  = (pMuon2.E()+pMuon2.Pz())/SQRT2;
        p1Minus = (pMuon2.E()-pMuon2.Pz())/SQRT2;
        p2Plus  = (pMuon1.E()+pMuon1.Pz())/SQRT2;
        p2Minus = (pMuon1.E()-pMuon1.Pz())/SQRT2;
      }
      mva.cosThetaStarCS = diMuon.Pz()/fabs(diMuon.Pz()) * 
                               2*(p1Plus*p2Minus-p1Minus*p2Plus) / 
                               (diMuon.Mag()*sqrt(diMuon.Mag2()+diMuon.Pt()*diMuon.Pt()));
    }

    // Computing nVtx Valid
    //for(unsigned iVtx=0;iVtx<vertexInfo.nVertices;iVtx++)
    //{
    //  if(vertexInfo.isValid[iVtx])
    //  {
    //    mva.nVtx++;
    //  }
    //}
    mva.nVtx = vertexInfo.nVertices;


    // Jet Part
    for(unsigned iJet=0; (iJet < jets.nJets && iJet < 10);iJet++)
    {
        if(jets.pt[iJet] > 30.0)
        {
          mva.nJets++;
          mva.ht += jets.pt[iJet];
        }
    }

    bool goodJets = false;
    if(jets.nJets>=2 && jets.pt[0]>30.0 && jets.pt[1]>30.0)
        goodJets = true;


    if(goodJets)
    {
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      pJet1.SetXYZM(jets.px[0],jets.py[0],jets.pz[0],jets.mass[0]);
      pJet2.SetXYZM(jets.px[1],jets.py[1],jets.pz[1],jets.mass[1]);
      TLorentzVector diJet = pJet1+pJet2;

      double dEtaJets = fabs(jets.eta[0]-jets.eta[1]);
      double etaJetProduct = jets.eta[0]*jets.eta[1];
      mva.deltaPhiJets = pJet1.DeltaPhi(pJet2);
      mva.deltaRJets = pJet1.DeltaR(pJet2);

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
      bool jetInRapidityGap=false;
      for(unsigned iJet=2; (iJet < jets.nJets && iJet < 10);iJet++)
      {
        if(jets.pt[iJet] > 30.0)
        {
          if(jets.eta[iJet] < etaMax && jets.eta[iJet] > etaMin)
          {
            jetInRapidityGap = true;
            mva.nJetsInRapidityGap++;
            mva.htInRapidityGap += jets.pt[iJet];
          }
        }
      }

#ifdef JETPUID
      for(unsigned iJet=0; (iJet < jets.nJets && iJet < 10);iJet++)
      {
        if(jets.pt[iJet]>30.0)
        {
          if (iJet==0)
            mva.puJetIDSimpleDiscJet1 = puJetSimpleDisc[iJet];
          else if (iJet==1)
            mva.puJetIDSimpleDiscJet2 = puJetSimpleDisc[iJet];
          else if (iJet==2)
            mva.puJetIDSimpleDiscJet3 = puJetSimpleDisc[iJet];

          if (iJet==0)
            mva.puJetIDSimpleJet1 = passPUJetID(puJetSimpleId[iJet],puJetLoose);
          else if (iJet==1)
            mva.puJetIDSimpleJet2 = passPUJetID(puJetSimpleId[iJet],puJetLoose);
          else if (iJet==2)
            mva.puJetIDSimpleJet3 = passPUJetID(puJetSimpleId[iJet],puJetLoose);

          //if (iJet==0)
          // {
          //  std::cout << "Disc: " << puJetSimpleDisc[iJet] << " Raw Id: " << puJetSimpleId[iJet] << std::endl;
          //}
        }
      }
#endif

      mva.mDiJet = diJet.M();
      mva.yDiJet = diJet.Rapidity();
      mva.ptDiJet = diJet.Pt();
      mva.ptJet1 = pJet1.Pt();
      mva.ptJet2 = pJet2.Pt();
      mva.etaJet1 = pJet1.Eta();
      mva.etaJet2 = pJet2.Eta();
      mva.productEtaJets = etaJetProduct;
      mva.deltaEtaJets = dEtaJets;
      mva.ptmiss = (diJet+recoCandVec).Pt();

    }
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.

    if (!(inBlindWindow && isData))
      mva.writeEvent();


    bool vbfPreselection = mva.mDiJet>300.0 && mva.deltaEtaJets>3.0 && mva.productEtaJets<0.0 && mva.nJetsInRapidityGap == 0;
    //if(vbfPreselection)
    //  std::cout << "VBF Preselected!!";
    mva.vbfPreselection = vbfPreselection;


    bool pt0to30 = !vbfPreselection  && mva.ptDiMu <30.;
    bool pt30to50 = !vbfPreselection && mva.ptDiMu> 30. && mva.ptDiMu <50.;
    bool pt50to125 = !vbfPreselection && mva.ptDiMu> 50.  && mva.ptDiMu <125.;
    bool pt125to250 = !vbfPreselection && mva.ptDiMu> 125.  && mva.ptDiMu <250.;
    bool pt250 = !vbfPreselection && mva.ptDiMu >250.;

    float bdtValInc =  mva.getMVA(cfgNameInc,"BDT");
    float bdtValVBF =  mva.getMVA(cfgNameVBF,"BDT");
    float likeValInc =  mva.getMVA(cfgNameInc,"Likelihood");
    float likeValVBF =  mva.getMVA(cfgNameVBF,"Likelihood");
    bool passIncBDTCut = mva.getMVAPassBDTCut(cfgNameInc);
    bool passVBFBDTCut = mva.getMVAPassBDTCut(cfgNameVBF);
    mva.bdtValInc = bdtValInc;
    mva.bdtValVBF = bdtValVBF;

    time_t timeStartFilling = time(NULL);

    // F I L L I N G
    double randForSF = randomForSF.Rndm();
    double effWeight     = weightFromSF(randForSF,muon1,muon2, 0.,    0.,    0.   );
    double effWeightUp   = weightFromSF(randForSF,muon1,muon2, 0.005, 0.005, 0.005);
    double effWeightDown = weightFromSF(randForSF,muon1,muon2,-0.005,-0.005,-0.005);

    //std::cout << "weight From SF="      << effWeight     << std::endl;
    //std::cout << "weight From SF Up="   << effWeightUp   << std::endl;
    //std::cout << "weight From SF Down=" << effWeightDown << std::endl;

    int selectionMask = whichSelection(muon1,muon2,
                                       allowedHLTPaths,
                                       runPeriod,
                                       jets,
                                       passIncBDTCut,
                                       passVBFBDTCut);

    int selectionMaskJECUp = whichSelection(muon1,muon2,
                                            allowedHLTPaths,
                                            runPeriod,
                                            jets,
                                            passIncBDTCut,
                                            passVBFBDTCut,
                                            +1); // JEC +1 sigma


    int selectionMaskJECDown = whichSelection(muon1,muon2,
                                              allowedHLTPaths,
                                              runPeriod,
                                              jets,
                                              passIncBDTCut,
                                              passVBFBDTCut,
                                              -1); // JEC -1 sigma


    int selectionMaskJERUp = whichSelection(muon1,muon2,
                                            allowedHLTPaths,
                                            runPeriod,
                                            jets,
                                            passIncBDTCut,
                                            passVBFBDTCut,
                                            0,1); // JER +1 sigma
    
    int selectionMaskJERDown = whichSelection(muon1,muon2,
                                              allowedHLTPaths,
                                              runPeriod,
                                              jets,
                                              passIncBDTCut,
                                              passVBFBDTCut,
                                              0,-1); // JER -1 sigma


    if (selectionMask != 0) {
      std::cout << "\n\nprinting the mask\n";
     
      std::cout << "selectionMask=" 
                << selectionMask 
                << " = " << hex << selectionMask << std::endl;

      std::cout << "selectionMaskJECUp=" 
                << selectionMaskJECUp 
                << " = " << hex << selectionMaskJECUp << std::endl;

      std::cout << "selectionMaskJECDown=" 
                << selectionMaskJECDown 
                << " = " << hex << selectionMaskJECDown << std::endl;

      std::cout << "selectionMaskJERUp=" 
                << selectionMaskJERUp 
                << " = " << hex << selectionMaskJERUp << std::endl;

      std::cout << "selectionMaskJERDown=" 
                << selectionMaskJERDown 
                << " = " << hex << selectionMaskJERDown << std::endl;

//       std::cout << "selectionMask & vbfPresel = " << dec << (selectionMask & vbfPresel) << std::endl;
//       std::cout << "selectionMask & incPresel = " << dec << (selectionMask & incPresel) << std::endl;
//       std::cout << "selectionMask & vbfPresel_passVBFBDTCut = " << dec << (selectionMask & vbfPresel_passVBFBDTCut) << std::endl;
//       std::cout << "selectionMask & vbfPresel_isBB_passVBFBDTCut = " << dec << (selectionMask & vbfPresel_isBB_passVBFBDTCut) << std::endl;
//       std::cout << "selectionMask & vbfPresel_isNotBB_passVBFBDTCut = " << dec << (selectionMask & vbfPresel_isNotBB_passVBFBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_passIncBDTCut = " << dec << (selectionMask & incPresel_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isBB_passIncBDTCut = " << dec << (selectionMask & incPresel_isBB_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isBO_passIncBDTCut = " << dec << (selectionMask & incPresel_isBO_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isBE_passIncBDTCut = " << dec << (selectionMask & incPresel_isBE_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isOO_passIncBDTCut = " << dec << (selectionMask & incPresel_isOO_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isOE_passIncBDTCut = " << dec << (selectionMask & incPresel_isOE_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isEE_passIncBDTCut = " << dec << (selectionMask & incPresel_isEE_passIncBDTCut) << std::endl;
//       std::cout << "selectionMask & incPresel_isNotBB_passIncBDTCut = " << dec << (selectionMask & incPresel_isNotBB_passIncBDTCut) << std::endl;
    }

    //VBF Preselected Plots
    if ( (selectionMask & vbfPresel) == vbfPresel)
    {
      histsPUVBFPresel.FillPU(mva,
                              weight,weightPUErrUp,
                              weightPUErrDown,weightPUSystShift);
      
      // # of events after full selection
      histsPUVBFPresel.nDiMuFullSel      -> Fill(0.);
      histsPUVBFPresel.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUVBFPresel.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUVBFPresel.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & vbfPresel) == vbfPresel ) histsPUVBFPresel.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & vbfPresel) == vbfPresel ) histsPUVBFPresel.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & vbfPresel) == vbfPresel ) histsPUVBFPresel.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & vbfPresel) == vbfPresel ) histsPUVBFPresel.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================





    //Inc Preselected Plots
    if ( (selectionMask & incPresel) == incPresel)
    {
      histsPUIncPresel.FillPU(mva,
                              weight,weightPUErrUp,
                              weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncPresel.nDiMuFullSel      -> Fill(0.);
      histsPUIncPresel.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncPresel.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncPresel.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel) == incPresel ) histsPUIncPresel.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel) == incPresel ) histsPUIncPresel.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel) == incPresel ) histsPUIncPresel.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel) == incPresel ) histsPUIncPresel.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================





    //VBF BDT Cut Plots
    if ( (selectionMask & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut )
    {
      histsPUVBFBDTCut.FillPU(mva,
                              weight,weightPUErrUp,
                              weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUVBFBDTCut.nDiMuFullSel      -> Fill(0.);
      histsPUVBFBDTCut.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUVBFBDTCut.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUVBFBDTCut.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUVBFBDTCut.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUVBFBDTCut.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUVBFBDTCut.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUVBFBDTCut.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================


    if ( (selectionMask & vbfPresel_isBB_passVBFBDTCut) == vbfPresel_isBB_passVBFBDTCut )
    {
      histsPUVBFBDTCutBB.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUVBFBDTCutBB.nDiMuFullSel      -> Fill(0.);
      histsPUVBFBDTCutBB.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUVBFBDTCutBB.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUVBFBDTCutBB.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & vbfPresel_isBB_passVBFBDTCut) == vbfPresel_isBB_passVBFBDTCut ) histsPUVBFBDTCutBB.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & vbfPresel_isBB_passVBFBDTCut) == vbfPresel_isBB_passVBFBDTCut ) histsPUVBFBDTCutBB.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & vbfPresel_isBB_passVBFBDTCut) == vbfPresel_isBB_passVBFBDTCut ) histsPUVBFBDTCutBB.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & vbfPresel_isBB_passVBFBDTCut) == vbfPresel_isBB_passVBFBDTCut ) histsPUVBFBDTCutBB.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & vbfPresel_isNotBB_passVBFBDTCut) == vbfPresel_isNotBB_passVBFBDTCut )
    {
      histsPUVBFBDTCutNotBB.FillPU(mva,
                                   weight,weightPUErrUp,
                                   weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUVBFBDTCutNotBB.nDiMuFullSel      -> Fill(0.);
      histsPUVBFBDTCutNotBB.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUVBFBDTCutNotBB.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUVBFBDTCutNotBB.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & vbfPresel_isNotBB_passVBFBDTCut) == vbfPresel_isNotBB_passVBFBDTCut ) histsPUVBFBDTCutNotBB.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & vbfPresel_isNotBB_passVBFBDTCut) == vbfPresel_isNotBB_passVBFBDTCut ) histsPUVBFBDTCutNotBB.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & vbfPresel_isNotBB_passVBFBDTCut) == vbfPresel_isNotBB_passVBFBDTCut ) histsPUVBFBDTCutNotBB.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & vbfPresel_isNotBB_passVBFBDTCut) == vbfPresel_isNotBB_passVBFBDTCut ) histsPUVBFBDTCutNotBB.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================




    //Inc BDT Cut Plots
    if ( (selectionMask & incPresel_passIncBDTCut) == incPresel_passIncBDTCut )
    {
      histsPUIncBDTCut.FillPU(mva,
                              weight,weightPUErrUp,
                              weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCut.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCut.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCut.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCut.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUIncBDTCut.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUIncBDTCut.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUIncBDTCut.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & vbfPresel_passVBFBDTCut) == vbfPresel_passVBFBDTCut ) histsPUIncBDTCut.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================




    if ( (selectionMask & incPresel_isBB_passIncBDTCut) == incPresel_isBB_passIncBDTCut )
    {
      histsPUIncBDTCutBB.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutBB.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCutBB.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutBB.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutBB.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isBB_passIncBDTCut) == incPresel_isBB_passIncBDTCut ) histsPUIncBDTCutBB.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isBB_passIncBDTCut) == incPresel_isBB_passIncBDTCut ) histsPUIncBDTCutBB.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isBB_passIncBDTCut) == incPresel_isBB_passIncBDTCut ) histsPUIncBDTCutBB.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isBB_passIncBDTCut) == incPresel_isBB_passIncBDTCut ) histsPUIncBDTCutBB.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & incPresel_isBO_passIncBDTCut) == incPresel_isBO_passIncBDTCut )
    {
      histsPUIncBDTCutBO.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutBO.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCutBO.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutBO.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutBO.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isBO_passIncBDTCut) == incPresel_isBO_passIncBDTCut ) histsPUIncBDTCutBO.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isBO_passIncBDTCut) == incPresel_isBO_passIncBDTCut ) histsPUIncBDTCutBO.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isBO_passIncBDTCut) == incPresel_isBO_passIncBDTCut ) histsPUIncBDTCutBO.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isBO_passIncBDTCut) == incPresel_isBO_passIncBDTCut ) histsPUIncBDTCutBO.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & incPresel_isBE_passIncBDTCut) == incPresel_isBE_passIncBDTCut )
    {
      histsPUIncBDTCutBE.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutBE.nDiMuFullSel -> Fill(0.);
      histsPUIncBDTCutBE.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutBE.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutBE.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isBE_passIncBDTCut) == incPresel_isBE_passIncBDTCut ) histsPUIncBDTCutBE.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isBE_passIncBDTCut) == incPresel_isBE_passIncBDTCut ) histsPUIncBDTCutBE.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isBE_passIncBDTCut) == incPresel_isBE_passIncBDTCut ) histsPUIncBDTCutBE.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isBE_passIncBDTCut) == incPresel_isBE_passIncBDTCut ) histsPUIncBDTCutBE.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & incPresel_isOO_passIncBDTCut) == incPresel_isOO_passIncBDTCut ) 
    {
      histsPUIncBDTCutOO.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutOO.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCutOO.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutOO.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutOO.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isOO_passIncBDTCut) == incPresel_isOO_passIncBDTCut ) histsPUIncBDTCutOO.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isOO_passIncBDTCut) == incPresel_isOO_passIncBDTCut ) histsPUIncBDTCutOO.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isOO_passIncBDTCut) == incPresel_isOO_passIncBDTCut ) histsPUIncBDTCutOO.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isOO_passIncBDTCut) == incPresel_isOO_passIncBDTCut ) histsPUIncBDTCutOO.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & incPresel_isOE_passIncBDTCut) == incPresel_isOE_passIncBDTCut )
    {
      histsPUIncBDTCutOE.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutOE.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCutOE.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutOE.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutOE.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isOE_passIncBDTCut) == incPresel_isOE_passIncBDTCut ) histsPUIncBDTCutOE.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isOE_passIncBDTCut) == incPresel_isOE_passIncBDTCut ) histsPUIncBDTCutOE.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isOE_passIncBDTCut) == incPresel_isOE_passIncBDTCut ) histsPUIncBDTCutOE.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isOE_passIncBDTCut) == incPresel_isOE_passIncBDTCut ) histsPUIncBDTCutOE.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & incPresel_isEE_passIncBDTCut) == incPresel_isEE_passIncBDTCut )
    {
      histsPUIncBDTCutEE.FillPU(mva,
                                weight,weightPUErrUp,
                                weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutEE.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCutEE.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutEE.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutEE.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isEE_passIncBDTCut) == incPresel_isEE_passIncBDTCut )  histsPUIncBDTCutEE.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isEE_passIncBDTCut) == incPresel_isEE_passIncBDTCut )  histsPUIncBDTCutEE.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isEE_passIncBDTCut) == incPresel_isEE_passIncBDTCut )  histsPUIncBDTCutEE.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isEE_passIncBDTCut) == incPresel_isEE_passIncBDTCut )  histsPUIncBDTCutEE.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================



    if ( (selectionMask & incPresel_isNotBB_passIncBDTCut) == incPresel_isNotBB_passIncBDTCut )
    {
      histsPUIncBDTCutNotBB.FillPU(mva,
                                   weight,weightPUErrUp,
                                   weightPUErrDown,weightPUSystShift);

      // # of events after full selection
      histsPUIncBDTCutNotBB.nDiMuFullSel      -> Fill(0.);
      histsPUIncBDTCutNotBB.nDiMuFullSelSF    -> Fill(0.,effWeight);    
      histsPUIncBDTCutNotBB.nDiMuFullSelSFUp  -> Fill(0.,effWeightUp);  
      histsPUIncBDTCutNotBB.nDiMuFullSelSFDown-> Fill(0.,effWeightDown);
    }

    if ( (selectionMaskJECUp   & incPresel_isNotBB_passIncBDTCut) == incPresel_isNotBB_passIncBDTCut ) histsPUIncBDTCutNotBB.nDiMuFullSelJECUp  -> Fill(0.);  
    if ( (selectionMaskJECDown & incPresel_isNotBB_passIncBDTCut) == incPresel_isNotBB_passIncBDTCut ) histsPUIncBDTCutNotBB.nDiMuFullSelJECDown-> Fill(0.);

    if ( (selectionMaskJERUp   & incPresel_isNotBB_passIncBDTCut) == incPresel_isNotBB_passIncBDTCut ) histsPUIncBDTCutNotBB.nDiMuFullSelJERUp  -> Fill(0.);  
    if ( (selectionMaskJERDown & incPresel_isNotBB_passIncBDTCut) == incPresel_isNotBB_passIncBDTCut ) histsPUIncBDTCutNotBB.nDiMuFullSelJERDown-> Fill(0.);
    // =========================================================================





    timeReading += difftime(timeStopReading,timeStartReading);
    timeProcessing += difftime(timeStartFilling,timeStopReading);
    timeFilling += difftime(time(NULL),timeStartFilling);
  }// end event loop
  time_t timeEndEventLoop = time(NULL);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


  histsPUVBFPresel.Write(outFile,"PUVBFPresel");
  histsPUIncPresel.Write(outFile,"PUIncPresel");

  histsPUIncBDTCut  .Write(outFile,"PUIncBDTCut");
  histsPUIncBDTCutBB.Write(outFile,"PUIncBDTCutBB");
  histsPUIncBDTCutBO.Write(outFile,"PUIncBDTCutBO");
  histsPUIncBDTCutBE.Write(outFile,"PUIncBDTCutBE");
  histsPUIncBDTCutOO.Write(outFile,"PUIncBDTCutOO");
  histsPUIncBDTCutOE.Write(outFile,"PUIncBDTCutOE");
  histsPUIncBDTCutEE.Write(outFile,"PUIncBDTCutEE");
  histsPUIncBDTCutNotBB.Write(outFile,"PUIncBDTCutNotBB");

  histsPUVBFBDTCut     .Write(outFile,"PUVBFBDTCut");
  histsPUVBFBDTCutBB   .Write(outFile,"PUVBFBDTCutBB");
  histsPUVBFBDTCutNotBB.Write(outFile,"PUVBFBDTCutNotBB");


  cout <<"#######################\n" << endl;
  cout << "testCounter: "<< testCounter << endl;
  cout << "Total Time: "<<std::setprecision(3) << difftime(time(NULL),timeStart)<<"\n";
  cout << "Setup Time: "<<std::setprecision(3) <<difftime(timeStartEventLoop,timeStart)<<"\n";
  cout << "Event Loop Time: "<<std::setprecision(3) 
        <<difftime(timeEndEventLoop,timeStartEventLoop)<< ", "<<std::setprecision(3) 
        <<difftime(timeEndEventLoop,timeStartEventLoop)/(std::min(nEvents,(unsigned) maxEvents))*1000.
        <<" s / 1000 events or "
        <<(std::min(nEvents,(unsigned) maxEvents))/difftime(timeEndEventLoop,timeStartEventLoop)*3600./1.0e6
        <<"M events/hour \n";
  cout << "  Read Time: "<<std::setprecision(3) << timeReading << std::endl;
  cout << "  Proc Time: "<<std::setprecision(3) << timeProcessing << std::endl;
  cout << "  Fill Time: "<<std::setprecision(3) << timeFilling << std::endl;
  cout << "  All Read Time: "<<std::setprecision(3) << timeReadingAll << std::endl;
  cout << "Wrapup Time: "<<std::setprecision(3) <<difftime(time(NULL),timeEndEventLoop)<<"\n";
  cout << "analyzer done." << endl << endl;
  return 0;
}




HistStruct::HistStruct()
{

  // PU
  unsigned nMassBins = 800 ;
  float minMass      =   0.;
  float maxMass      = 400.;

  mDiMu = new TH1F("mDiMu","DiMuon Mass",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMu);

  mDiMuPUErrUp = new TH1F("mDiMuPUErrUp","DiMuon Mass PU ErrUp",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuPUErrUp);

  mDiMuPUErrDown = new TH1F("mDiMuPUErrDown","DiMuon Mass PU ErrDown",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuPUErrDown);

  mDiMuPUSystShift = new TH1F("mDiMuPUSystShift","DiMuon Mass PU Syst. Shift",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuPUSystShift);

  // Efficiency
  // We are simply counting  
  unsigned nEffBins = 1;
  float minEff      = 0.;
  float maxEff      = 1.;

  nDiMuGen = new TH1F("nDiMuGen","Counting Gen. Events",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuGen);

  nDiMuGenInAcc = new TH1F("nDiMuGenInAcc","Counting Gen. Events in Acc",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuGenInAcc);

  nDiMuFullSel = new TH1F("nDiMuFullSel","Counting Reco Events after Full Sel.",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSel);

  nDiMuFullSelSF = new TH1F("nDiMuFullSelSF","Counting Reco Events after Full Sel. and SF",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelSF);

  // varying all the systematic variations UP
  nDiMuFullSelSFUp = new TH1F("nDiMuFullSelSFUp","Counting Reco Events after Full Sel. and SF varied UP",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelSFUp);

  // varying all the systematic variations DOWN
  nDiMuFullSelSFDown = new TH1F("nDiMuFullSelSFDown","Counting Reco Events after Full Sel. and SF varied DOWN",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelSFDown);

  // JEC
  nDiMuFullSelJECUp   = new TH1F("nDiMuFullSelJECUp","Counting Reco Events after Full Sel. and JEC +1 sigma",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelJECUp);
  nDiMuFullSelJECDown = new TH1F("nDiMuFullSelJECDown","Counting Reco Events after Full Sel. and JEC -1 sigma",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelJECDown);

  // JER
  nDiMuFullSelJERUp   = new TH1F("nDiMuFullSelJERUp","Counting Reco Events after Full Sel. and JER +1 sigma",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelJERUp);
  nDiMuFullSelJERDown = new TH1F("nDiMuFullSelJERDown","Counting Reco Events after Full Sel. and JER -1 sigma",nEffBins,minEff,maxEff);
  histVec.push_back(nDiMuFullSelJERDown);

}

HistStruct::~HistStruct()
{
  //std::vector<TH1F*>::iterator hist;

  //for(hist = histVec.begin();hist != histVec.end(); hist++)
  //delete *hist;
}

void
HistStruct::Write(TFile* outfile, std::string directory)
{
  if(directory == "")
  {
    outfile->cd();
  }
  else
  {
    TDirectory* dir = outfile->mkdir(directory.c_str());
    dir->cd();
  }

  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    (*hist)->Write();

  outfile->cd();
}


void 
HistStruct::FillPU(const MVA& mva,
                   const double weight,
                   const double weightPUErrUp,
                   const double weightPUErrDown,
                   const double weightPUSystShift)
{
  
  mDiMu            -> Fill(mva.mDiMu, weight           );
  mDiMuPUErrUp     -> Fill(mva.mDiMu, weightPUErrUp    );
  mDiMuPUErrDown   -> Fill(mva.mDiMu, weightPUErrDown  );
  mDiMuPUSystShift -> Fill(mva.mDiMu, weightPUSystShift);

}

