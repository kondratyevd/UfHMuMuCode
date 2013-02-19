#include <algorithm>
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
#include <TRandom3.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "DataFormats.h"
#include "helpers.h"
#include "mva.h"
#include "LumiReweightingStandAlone.h"

#include "annaCalibCode/SmearingTool.h"

#include "src/ScaleFactors.h"
#include "src/ScaleFactors_2011.h"

//for Lorentz Vector tools
//#include<TMath.h>
//#include "Math/VectorUtil_Cint.h"
//#include "Math/GenVector/LorentzVector.h"
//#include <TLorentzVector.h>
///

#define JETPUID
#define PUREWEIGHT
#define SMEARING
#define ISMEAR 1 
//#define ROCHESTER
#define MUSCLEFIT

#ifdef ROCHESTER
#include "rochester/rochcor2012.h"
#include "rochester/rochcor.h"
#endif
#ifdef MUSCLEFIT
#include "musclefit/MuScleFitCorrector.h"
#endif

// for check doublets in Gen level:
#include <set>
typedef std::pair<float,float> pairOfDouble;


using namespace std;
using namespace boost;

  TString ExtraInfo = "Higgs";
  //TString ExtraInfo = "Zmumu";

  Double_t MASS_MUON = 0.105658367;    //GeV/c2

void fillMuonHist(TH1F* hist, _MuonInfo& mu1, _MuonInfo& mu2);
void printStationMiss(_MuonInfo& mu1, _MuonInfo& mu2, _EventInfo& eventInfo, std::string & testString, unsigned & testCounter);

struct HistStruct
{
  HistStruct();
  ~HistStruct();
  void Write(TFile* outfile, std::string directory);
  void Fill(const MVA& mva, bool blind);

  std::vector<TH1F*> histVec;
  std::vector<TH2F*> histVec2D;

  TH1F* mDiMu;
  TH1F* mDiMuResSigUp;
  TH1F* mDiMuResSigDown;
  TH1F* mDiMuResASigUp;
  TH1F* mDiMuResASigDown;

  TH1F* mDiJet;
  TH1F* ptDiMu;
  TH1F* ptDiJet;
  TH1F* yDiMu;
  TH1F* yDiJet;

  TH2F* yVptDiMu;
  TH2F* ptVmDiMu;
  TH2F* yVmDiMu;
  TH2F* phiVmDiMu;

  TH1F* ptMu1;
  TH1F* ptMu2;
  TH1F* ptJet1;
  TH1F* ptJet2;

  TH1F* etaMu1;
  TH1F* etaMu2;
  TH1F* etaJet1;
  TH1F* etaJet2;

  TH1F* deltaEtaJets;
  TH1F* deltaPhiJets;
  TH1F* deltaRJets;
  TH1F* deltaPhiHJ1;

  TH1F* deltaEtaMuons;
  TH1F* deltaPhiMuons;
  TH1F* deltaRMuons;

  TH1F* countsHist;
  TH1F* countsHist2;

  TH1F* cosThetaStar;
  TH1F* cosThetaStarCS;

  TH1F* puJetIDSimpleDiscJet1;
  TH1F* puJetIDSimpleDiscJet2;
  TH1F* puJetIDSimpleDiscJet3;

  TH1F* puJetIDSimpleJet1;
  TH1F* puJetIDSimpleJet2;
  TH1F* puJetIDSimpleJet3;

  TH1F* BDTHistMuonOnly;

  TH1F* BDTHistVBF;

  TH2F* BDTHistMuonOnlyVMass;
  TH2F* BDTHistVBFVMass;

  TH1F* relIsoMu1;
  TH1F* relIsoMu2;

  TH1F* nJets;
  TH1F* ht;
  TH1F* nJetsInRapidityGap;
  TH1F* htInRapidityGap;

  TH1F* nPU;
  TH1F* nVtx;
  TH1F* met;
  TH1F* ptmiss;
  TH1F* weight;
};

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  set<pairOfDouble> uniqueGeneratedEvents;

  gErrorIgnoreLevel = kError;
  time_t timeStart = time(NULL);

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("trainingTree,t",program_options::value<string>(), "Create Training Tree File with filename")
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

  //float minMmm = 70.0;
  //float maxMmm = 400.0;
  float minMmm = 110.0;
  float maxMmm = 150.0;

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
  //bool isSignal = false;
  bool isSignal = true;
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
    muonIdFuncPtr = &isKinTight_2011_noAcc_noIso;
  }
  else
  {
    cout << "Using 2012 Tight Muon Selection\n";
    muonIdFuncPtr = &isKinTight_2012_noAcc_noIso;
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
  bool trainingTreeRun=false;
  if (optionMap.count("trainingTree")) 
  {
      cout << "Training enabled" << "\n";
      trainingTreeFileName = optionMap["trainingTree"].as<string>();
      cout << "Training Tree File Name: " << trainingTreeFileName << "\n";
      trainingTreeRun=true;
  }

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

   /// Higgs Boson 
  _genPartInfo genHpostFSR;
  if(!isData && tree->GetBranchStatus("genHpostFSR"))
    tree->SetBranchAddress("genHpostFSR", &genHpostFSR);

   _TrackInfo reco1GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM1HpostFSR"))
    tree->SetBranchAddress("genM1HpostFSR", &reco1GenPostFSR);

  _TrackInfo reco2GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM2HpostFSR"))
    tree->SetBranchAddress("genM2HpostFSR", &reco2GenPostFSR);

   /// Z Boson
  _genPartInfo genZpostFSR;
  if(!isData && tree->GetBranchStatus("genZpostFSR"))
    tree->SetBranchAddress("genZpostFSR", &genZpostFSR);

  _TrackInfo recoZ1GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM1ZpostFSR"))
    tree->SetBranchAddress("genM1ZpostFSR", &recoZ1GenPostFSR);

  _TrackInfo recoZ2GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM2ZpostFSR"))
    tree->SetBranchAddress("genM2ZpostFSR", &recoZ2GenPostFSR);
   ///

  _PFJetInfo jets;
  tree->SetBranchAddress("pfJets",&jets);

#ifdef JETPUID
  float puJetFullDisc[10];
  float puJetSimpleDisc[10];
  float puJetCutDisc[10];

  tree->SetBranchAddress("puJetFullDisc",&puJetFullDisc);
  tree->SetBranchAddress("puJetSimpleDisc",&puJetSimpleDisc);
  tree->SetBranchAddress("puJetCutDisc",&puJetCutDisc);

  float puJetFullId[10];
  float puJetSimpleId[10];
  float puJetCutId[10];

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

  //////////////////////////
  // Histograms
  HistStruct hists;
  HistStruct histsBB;
  HistStruct histsBO;
  HistStruct histsBE;
  HistStruct histsOO;
  HistStruct histsOE;
  HistStruct histsEE;
  HistStruct histsNotBB;

  HistStruct hists4GeVWindow;

  HistStruct histsIncPresel; //
  HistStruct histsIncPreselBB;
  HistStruct histsIncPreselBO;
  HistStruct histsIncPreselBE;
  HistStruct histsIncPreselOO;
  HistStruct histsIncPreselOE;
  HistStruct histsIncPreselEE;
  HistStruct histsIncPreselNotBB;

  HistStruct histsVBFPresel; //
  HistStruct histsVBFPreselBB;
  HistStruct histsVBFPreselNotBB;

  HistStruct histsIncBDTCut;      //     
  HistStruct histsIncBDTCutBB;    //
  HistStruct histsIncBDTCutBO;    //
  HistStruct histsIncBDTCutBE;    //
  HistStruct histsIncBDTCutOO;    //
  HistStruct histsIncBDTCutOE;    //
  HistStruct histsIncBDTCutEE;    //
  HistStruct histsIncBDTCutNotBB; //
                                  
  HistStruct histsVBFBDTCut;      //
  HistStruct histsVBFBDTCutBB;    //   
  HistStruct histsVBFBDTCutNotBB; //

  HistStruct histsIncPreselDiMuPtL20;
  HistStruct histsVBFPreselDiMuPtL20;

  HistStruct histsIncPreselPUJETID;
  HistStruct histsVBFPreselPUJETID;

  HistStruct histsIncPreselPUJETIDForVeto;
  HistStruct histsVBFPreselPUJETIDForVeto;

  HistStruct histsVBFPreselPtMiss50Veto;

  HistStruct histsIncPreselPtG10;
  HistStruct histsIncPreselPtG10BB;
  HistStruct histsIncPreselPtG10BO;
  HistStruct histsIncPreselPtG10BE;
  HistStruct histsIncPreselPtG10OO;
  HistStruct histsIncPreselPtG10OE;
  HistStruct histsIncPreselPtG10EE;
  HistStruct histsIncPreselPtG10NotBB;

  //////////////////////////
  //for MVA

  std::vector<std::string> mvaConfigNames;
  mvaConfigNames.push_back(cfgNameInc);
  mvaConfigNames.push_back(cfgNameVBF);
  MVA mva(mvaConfigNames,trainingTreeFileName);

  TRandom3 random(1457);

  //////////////////////////
  //for PU reweighting

#ifdef PUREWEIGHT
  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012ABCD.root","pileup","pileup");
  if (runPeriod == "7TeV")
  {
    cout << "Using 2011AB PU reweighting\n";
    lumiWeights = reweight::LumiReWeighting("pileupDists/PileUpHistMCFall11.root","pileupDists/PileUpHist2011AB.root","pileup","pileup");
  }
  else
  {
    cout << "Using 2012ABCD PU reweighting\n";
  }
#endif

  /////////////////////////////
  // Muon Momentum Corrections

#ifdef MUSCLEFIT
  MuScleFitCorrector* mfCorr;
  TString mfInFile;
  if (isData)
  {
    if(runPeriod == "8TeV")
      mfInFile = "musclefit/MuScleFit_2012_DATA_53X.txt";
    else
      mfInFile = "musclefit/MuScleFit_2011_DATA_42X.txt";

    mfCorr = new MuScleFitCorrector(mfInFile);
  }
//  else
//  {
//    if(runPeriod == "7TeV")
//      mfInFile = "musclefit/MuScleFit_2011_MC_42X.txt";
//    else
//      mfInFile = "musclefit/MuScleFit_2012_MC_52X.txt";
//
//    mfCorr = new MuScleFitCorrector(mfInFile);
//  }
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

  TRandom3 randomForSF(123412845);

//////////////////////////////////////////////////////////////////////
  int counterGenBoson  = 0; // mass cut only
  int counterMinimCuts = 0; // mass, charge, vertex, cosmic
  int counterAccCuts   = 0; // pt, eta
  int counterIDCuts    = 0; // tight muon id
  int counterIsoCuts   = 0; // trk Iso
  int counterTrigCuts  = 0; // iso mu 24_eta2p1
  int counterJetSel    = 0; // jet selection
  int counterBDTvbfCut = 0; // VBF BDT cut
  int counterNonJetSel    = 0; // jet selection
  int counterDiPt10GeV = 0; // jet selection

  float EffMinimCuts = 0; // mass, charge, vertex, cosmic
  float EffAccCuts   = 0; // pt, eta
  float EffIDCuts    = 0; // tight muon id
  float EffIsoCuts   = 0; // trk Iso
  float EffTrigCuts  = 0; // iso mu 24_eta2p1
  float EffJetSel    = 0; // jet selection
  float EffBDTvbfCut = 0; // VBF BDT cut
  float EffNonJetSel    = 0; // jet selection
  float EffDiPt10GeV = 0; // pt(mumu) > = 10 GeV/c 

  float dEffMinimCuts = 0; // mass, charge, vertex, cosmic
  float dEffAccCuts   = 0; // pt, eta
  float dEffIDCuts    = 0; // tight muon id
  float dEffIsoCuts   = 0; // trk Iso
  float dEffTrigCuts  = 0; // iso mu 24_eta2p1
  float dEffJetSel    = 0; // jet selection
  float dEffBDTvbfCut = 0; // VBF BDT cut
  float dEffNonJetSel    = 0; // jet selection
  float dEffDiPt10GeV = 0; // pt(mumu) > = 10 GeV/c 
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

    if(ExtraInfo == "Zmumu"){
       if(recoZ1GenPostFSR.pt < -900. || recoZ2GenPostFSR.pt < -900) continue;//rejection fake in gen level 
       if(recoZ1GenPostFSR.eta < -900. || recoZ2GenPostFSR.eta < -900) continue;//rejection fake in gen level 
    }
    if(ExtraInfo == "Higgs"){
       if(reco1GenPostFSR.pt < -900. || reco2GenPostFSR.pt < -900) continue;//rejection fake in gen level 
       if(reco1GenPostFSR.eta < -900. || reco2GenPostFSR.eta < -900) continue;//rejection fake in gen level 
    }

    // calculate check for gen Mass
    float MassGenBoson;
          if(ExtraInfo == "Zmumu")MassGenBoson = genZpostFSR.mass;
          if(ExtraInfo == "Higgs")MassGenBoson = genHpostFSR.mass;

    // check for fhe doublets:
                pairOfDouble massPt(genHpostFSR.mass,genHpostFSR.pt);
                //if (uniqueGeneratedEvents.insert(massPt).second) continue;
                if(MassGenBoson >= minMmm && MassGenBoson < maxMmm && uniqueGeneratedEvents.insert(massPt).second){
                //if(MassGenBoson >= minMmm && MassGenBoson < maxMmm && uniqueGeneratedEvents.insert(massPt).second == false){
                    counterGenBoson++;
                    //hMassGenBoson ->Fill(MassGenBoson);
                }
    // additional selection cuts
    if (reco1.charge == reco2.charge) continue;
    if(reco1.pt < 0. || reco2.pt < 0) continue;//rejection fake in reco level 
    if(reco1.pt < -900. || reco2.pt < -900) continue;//rejection fake in reco level 
    if(reco1.eta < -900. || reco2.eta < -900) continue;//rejection fake in reco level

// make a matching between reco and gen muons
// very improtant for Smearing tool to not assign high pt to low pt muons and not get fakes
                TLorentzVector MuReco1, MuReco2, MuTrue1, MuTrue2;
                int MuTrue1charge, MuTrue2charge;
                float MuTrue1ptErr, MuTrue2ptErr;
                if(ExtraInfo == "Zmumu")MuTrue1charge = recoZ1GenPostFSR.charge;
                if(ExtraInfo == "Zmumu")MuTrue2charge = recoZ2GenPostFSR.charge;
                if(ExtraInfo == "Zmumu")MuTrue1ptErr = recoZ1GenPostFSR.ptErr;
                if(ExtraInfo == "Zmumu")MuTrue2ptErr = recoZ2GenPostFSR.ptErr;
                if(ExtraInfo == "Higgs")MuTrue1charge = reco1GenPostFSR.charge;
                if(ExtraInfo == "Higgs")MuTrue2charge = reco2GenPostFSR.charge;
                if(ExtraInfo == "Higgs")MuTrue1ptErr = reco1GenPostFSR.ptErr;
                if(ExtraInfo == "Higgs")MuTrue2ptErr = reco2GenPostFSR.ptErr;
                if(ExtraInfo == "Zmumu")MuTrue1.SetPtEtaPhiM(recoZ1GenPostFSR.pt, recoZ1GenPostFSR.eta, recoZ1GenPostFSR.phi, MASS_MUON);
                if(ExtraInfo == "Zmumu")MuTrue2.SetPtEtaPhiM(recoZ2GenPostFSR.pt, recoZ2GenPostFSR.eta, recoZ2GenPostFSR.phi, MASS_MUON);
                if(ExtraInfo == "Higgs")MuTrue1.SetPtEtaPhiM(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.phi, MASS_MUON);
                if(ExtraInfo == "Higgs")MuTrue2.SetPtEtaPhiM(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.phi, MASS_MUON);
                MuReco1.SetPtEtaPhiM(reco1.pt, reco1.eta, reco1.phi, MASS_MUON);
                MuReco2.SetPtEtaPhiM(reco2.pt, reco2.eta, reco2.phi, MASS_MUON);

////////////////////////
                float DiMuonPtNonCorr = (MuReco1+MuReco2).Pt();
                float muonRes_t1r1 = -10.;
                float muonRes_t1r2 = -10.;
                float muonRes_t2r1 = -10.;
                float muonRes_t2r2 = -10.;
                if(MuTrue1.Pt() > 0.) muonRes_t1r1 = (MuReco1.Pt()-MuTrue1.Pt())/MuTrue1.Pt();
                if(MuTrue2.Pt() > 0.) muonRes_t2r2 = (MuReco2.Pt()-MuTrue2.Pt())/MuTrue2.Pt();
                if(MuTrue1.Pt() > 0.) muonRes_t1r2 = (MuReco2.Pt()-MuTrue1.Pt())/MuTrue1.Pt();
                if(MuTrue2.Pt() > 0.) muonRes_t2r1 = (MuReco1.Pt()-MuTrue2.Pt())/MuTrue2.Pt();
                //float deltaR_t1r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco1);
                //float deltaR_t1r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue1, MuReco2);
                //float deltaR_t2r1 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco1);
                //float deltaR_t2r2 = ROOT::Math::VectorUtil::DeltaR(MuTrue2, MuReco2);
                float deltaR_t1r1 = MuTrue1.DeltaR(MuReco1);
                float deltaR_t1r2 = MuTrue1.DeltaR(MuReco2);
                float deltaR_t2r1 = MuTrue2.DeltaR(MuReco1);
                float deltaR_t2r2 = MuTrue2.DeltaR(MuReco2);
                float DR1 = deltaR_t1r1;   
                float DR2 = deltaR_t2r2;   
                //make matching between MuTrue1,2 and MuReco1,2:
                if(deltaR_t1r1 > deltaR_t1r2){
                     if(ExtraInfo == "Zmumu"){
                            MuTrue1.SetPtEtaPhiM(recoZ2GenPostFSR.pt, recoZ2GenPostFSR.eta, recoZ2GenPostFSR.phi, MASS_MUON);
                            MuTrue1charge = recoZ2GenPostFSR.charge;
                            MuTrue1ptErr  = recoZ2GenPostFSR.ptErr;
                     }
                     if(ExtraInfo == "Higgs"){
                            MuTrue1.SetPtEtaPhiM(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.phi, MASS_MUON);
                            MuTrue1charge = reco2GenPostFSR.charge;
                            MuTrue1ptErr  = reco2GenPostFSR.ptErr;
                     }
                     //cout << "switch gen muons: D11 = " << deltaR_t1r1 << " and D12 = " << deltaR_t1r2 << endl;
                     muonRes_t1r1 = muonRes_t2r1;
                     DR1          = deltaR_t2r1;
                }     
                if(deltaR_t2r2 > deltaR_t2r1){
                     if(ExtraInfo == "Zmumu"){
                            MuTrue2.SetPtEtaPhiM(recoZ1GenPostFSR.pt, recoZ1GenPostFSR.eta, recoZ1GenPostFSR.phi, MASS_MUON);
                            MuTrue2charge = recoZ1GenPostFSR.charge;
                            MuTrue2ptErr  = recoZ1GenPostFSR.ptErr;
                     }
                     if(ExtraInfo == "Higgs"){
                            MuTrue2.SetPtEtaPhiM(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.phi, MASS_MUON);
                            MuTrue2charge = reco1GenPostFSR.charge;
                            MuTrue2ptErr = reco1GenPostFSR.ptErr;
                     }
                     muonRes_t2r2 = muonRes_t1r2;
                     DR2          = deltaR_t1r2;
                }
                //change gen level to match gen1 -> reco1 && gen2->reco2 
                float PtGen1Orig = reco1GenPostFSR.pt;
                float PtGen2Orig = reco2GenPostFSR.pt;
                if(ExtraInfo == "Higgs"){
                  reco1GenPostFSR.pt = MuTrue1.Pt();  
                  reco1GenPostFSR.eta = MuTrue1.Eta();  
                  reco1GenPostFSR.phi = MuTrue1.Phi();  
                  reco1GenPostFSR.ptErr = MuTrue1ptErr;  
                  reco1GenPostFSR.charge = MuTrue1charge;  
                  reco2GenPostFSR.pt = MuTrue2.Pt();  
                  reco2GenPostFSR.eta = MuTrue2.Eta();  
                  reco2GenPostFSR.phi = MuTrue2.Phi();  
                  reco2GenPostFSR.ptErr = MuTrue2ptErr;  
                  reco2GenPostFSR.charge = MuTrue2charge;  
                }
// reject fake muon matches for smearing
#ifdef SMEARING
       if(DR1 > 0.5 || DR2 > 0.5) continue;
       if(fabs(muonRes_t1r1) > 0.4 || fabs(muonRes_t2r2) > 0.4) continue;
#endif
// end: make a matching between reco and gen muons   


    double weight = 1.0;
    if (isSignal)
    {
        double randForSF = randomForSF.Rndm();
        if (runPeriod == "7TeV")
            weight *= weightFromSF_2011(randForSF,reco1,reco2,0.,0.,0.);
        else
            weight *= weightFromSF(randForSF,reco1,reco2,0.,0.,0.);
    }
    float pTreco1Orig = reco1.pt;
    float pTreco2Orig = reco2.pt;
    float diMOrig     = recoCandMass; 

    TLorentzVector reco1Vec;
    TLorentzVector reco2Vec;
    reco1Vec.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
    reco2Vec.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
#ifdef MUSCLEFIT
    if (isData && runPeriod == "8TeV")
    {
      TLorentzVector reco1Cor;
      TLorentzVector reco2Cor;
      reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
      reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
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
    }
#endif
#ifdef ROCHESTER
    TLorentzVector reco1Cor;
    TLorentzVector reco2Cor;
    reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
    reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
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
    reco1Vec = recoCor1;
    reco2Vec = recoCor2;
#endif

    TLorentzVector recoCandVec = reco1Vec+reco2Vec;

    float mDiMuResSigUp = recoCandMass;
    float mDiMuResSigDown = recoCandMass;
    float mDiMuResASigUp = recoCandMass;
    float mDiMuResASigDown = recoCandMass;

#ifdef SMEARING
    if(isSignal && runPeriod == "8TeV")
    {
      if(reco1GenPostFSR.pt<0.)
        cout << "Muon 1 Post FSR not valid!\n";
      if(reco2GenPostFSR.pt<0.)
        cout << "Muon 2 Post FSR not valid!\n";
      float ptReco1 = smearPT -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR);
      float ptReco2 = smearPT -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR);
      TLorentzVector reco1Vec;
      TLorentzVector reco2Vec;
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
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
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResSigUp = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR,"sig1",-1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"sig1",-1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResSigDown = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR,"Asig2Var",1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"Asig2Var",1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResASigUp = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR,"Asig2Var",-1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR,"Asig2Var",-1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResASigDown = diMuonVec.M();
      
    }
#endif

    fillMuonHist(hists.countsHist2, reco1, reco2);
    //printStationMiss(reco1,reco2,eventInfo,testString,testCounter);

    mva.resetValues();
    mva.mDiMu = recoCandMass;
    bool inBlindWindow = mva.mDiMu < maxBlind && mva.mDiMu > minBlind;

    bool blind = false;
  
#ifdef BLIND
    if (inBlindWindow && isData)
        blind = true;
#endif

#ifdef PUREWEIGHT
    if (!isData)
    {
      weight *= lumiWeights.weight(nPU);
    }
#endif

    if (mva.mDiMu < minMmm || mva.mDiMu > maxMmm)
        continue;
    counterMinimCuts++;
    // acceptance cuts
    if (reco1.pt < 25)         continue;
    if (fabs(reco1.eta) > 2.1) continue;

    if (reco2.pt < 25)         continue;
    if (fabs(reco2.eta) > 2.1) continue;

    counterAccCuts++;

    if (!((*muonIdFuncPtr)(reco1)) || !((*muonIdFuncPtr)(reco2)))
          continue;
    counterIDCuts++;

    if (getRelIso(reco1) > 0.12) continue;
    if (getRelIso(reco2) > 0.12) continue;

    counterIsoCuts++;


    if (!isHltMatched(reco1,reco2,allowedHLTPaths))
        continue;
    counterTrigCuts++;
    if(fabs(muonRes_t1r1) > 0.3 || fabs(muonRes_t2r2) > 0.3 || DR1 > 0.4 || DR2 > 0.4){
       cout << "DR1 = " << DR1 << "  DR2 = " << DR2 << endl;
       cout << "Res1 = " << fabs(muonRes_t1r1) << " Res2 = " << fabs(muonRes_t2r2) <<endl;
       cout << "PtGen1Orig = " << PtGen1Orig << " PtGen1 = " << reco1GenPostFSR.pt << " MuTrue1.Pt = " << MuTrue1.Pt() << endl;   
       cout << "PtGen2Orig = " << PtGen2Orig << " PtGen2 = " << reco2GenPostFSR.pt << " MuTrue2.Pt = " << MuTrue2.Pt() << endl;   
       cout << "pTreco1Orig = " << pTreco1Orig << " pTreco1Cor = " << reco1.pt << endl;
       cout << "pTreco2Orig = " << pTreco2Orig << " pTreco2Cor = " << reco2.pt << endl;
       cout << "diMOrig = " << diMOrig << " mva.mDiMu = " << mva.mDiMu << endl;
       cout << "****************" << endl;
    } 


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

    bool inTrainingWindow = (mva.mDiMu < 160. && mva.mDiMu > 70.);
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
    pMuon1.SetPtEtaPhiM(muon1.pt,muon1.eta,muon1.phi,MASS_MUON);
    pMuon2.SetPtEtaPhiM(muon2.pt,muon2.eta,muon2.phi,MASS_MUON);
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

    //////////////////////////////////////////
    // Filling Hists

    if (!blind)
    {
      hists.mDiMu->Fill(mva.mDiMu, weight);
      hists.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      hists.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      hists.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      hists.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);
      hists.yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu), weight);
      hists.ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu, weight);
      hists.phiVmDiMu->Fill(mva.mDiMu,recoCandPhi, weight);
    }

    hists.yDiMu->Fill(mva.yDiMu, weight);
    hists.ptDiMu->Fill(mva.ptDiMu, weight);
    hists.yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu), weight);
    hists.ptMu1->Fill(mva.ptMu1, weight);
    hists.ptMu2->Fill(mva.ptMu2, weight);
    hists.etaMu1->Fill(mva.etaMu1, weight);
    hists.etaMu2->Fill(mva.etaMu2, weight);
    hists.cosThetaStar->Fill(mva.cosThetaStar, weight);
    hists.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);

    hists.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
    hists.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
    hists.deltaRMuons->Fill(mva.deltaRMuons, weight);

    hists.relIsoMu1->Fill(mva.relIsoMu1, weight);
    hists.relIsoMu2->Fill(mva.relIsoMu2, weight);

    hists.nPU->Fill(nPU, weight);
    hists.nVtx->Fill(mva.nVtx, weight);
    hists.met->Fill(met.pt, weight);
    hists.weight->Fill(weight);

    // Jet Part
    for(unsigned iJet=0; (iJet < jets.nJets && iJet < 10);iJet++)
    {
        if (jets.genPt[iJet]>0.0 && jets.pt[iJet]>15.)
          jets.pt[iJet] = jerCorr(jets.pt[iJet],jets.genPt[iJet],jets.eta[iJet]);
        if(jets.pt[iJet] > 30.0 && fabs(jets.eta[iJet])<4.7)
        {
          mva.nJets++;
          mva.ht += jets.pt[iJet];
        }
    }

    bool goodJets = false;
    if(jets.nJets>=2 && jets.pt[0]>30.0 && jets.pt[1]>30.0 && fabs(jets.eta[0])<4.7 && fabs(jets.eta[1])<4.7)
        goodJets = true;

//    if (mva.mDiMu > 140. && mva.mDiMu < 150. && goodJets)
//    {
//        testCounter++;
//        //std::cout <<eventInfo.run <<":"<<eventInfo.event <<"\n"<< std::endl;
//        testString.appendAny(eventInfo.run);
//        testString.append(":");
//        testString.appendAny(eventInfo.event);
//        testString.append("\n");
//    }


    //if(goodJets == false) continue;
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
      mva.deltaPhiHJ1 = pJet1.DeltaPhi(diMuon);

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
            mva.puJetIDSimpleJet1 = passPUJetID(int(puJetSimpleId[iJet]),puJetLoose);
          else if (iJet==1)
            mva.puJetIDSimpleJet2 = passPUJetID(int(puJetSimpleId[iJet]),puJetLoose);
          else if (iJet==2)
            mva.puJetIDSimpleJet3 = passPUJetID(int(puJetSimpleId[iJet]),puJetLoose);
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

      hists.mDiJet->Fill(mva.mDiJet, weight);
      hists.ptDiJet->Fill(mva.ptDiJet, weight);
      hists.yDiJet->Fill(mva.yDiJet, weight);
      hists.ptJet1->Fill(mva.ptJet1, weight);
      hists.ptJet2->Fill(mva.ptJet2, weight);
      hists.etaJet1->Fill(mva.etaJet1, weight);
      hists.etaJet2->Fill(mva.etaJet2, weight);
      hists.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      hists.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      hists.deltaRJets->Fill(mva.deltaRJets, weight);
      hists.deltaPhiHJ1->Fill(mva.deltaPhiHJ1, weight);
      hists.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      hists.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      hists.nJets->Fill(mva.nJets, weight);
      hists.ht->Fill(mva.ht, weight);
      hists.ptmiss->Fill(mva.ptmiss, weight);

      hists.puJetIDSimpleDiscJet1->Fill(mva.puJetIDSimpleDiscJet1,weight);
      hists.puJetIDSimpleDiscJet2->Fill(mva.puJetIDSimpleDiscJet2,weight);
      hists.puJetIDSimpleDiscJet3->Fill(mva.puJetIDSimpleDiscJet3,weight);

      hists.puJetIDSimpleJet1->Fill(mva.puJetIDSimpleJet1,weight);
      hists.puJetIDSimpleJet2->Fill(mva.puJetIDSimpleJet2,weight);
      hists.puJetIDSimpleJet3->Fill(mva.puJetIDSimpleJet3,weight);
    }
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.

    if (!(inBlindWindow && isData) && inTrainingWindow)
      mva.writeEvent();

    if (trainingTreeRun) //Skip Filling of histos when training Tree
        continue;

    bool vbfPreselection = mva.mDiJet>300.0 && mva.deltaEtaJets>3.0 && mva.productEtaJets<0.0;
    //if(vbfPreselection)
    //  std::cout << "VBF Preselected!!";
    mva.vbfPreselection = vbfPreselection;

    if(vbfPreselection) counterJetSel++;
    if(!vbfPreselection) counterNonJetSel++;
    if ((!vbfPreselection) && mva.ptDiMu >= 10.) counterDiPt10GeV++;

    if(vbfPreselection)
        hists.countsHist->Fill(6);
    else
        hists.countsHist->Fill(5);
    

    bool vbfVeryTight = false;
    bool vbfTight = false;
    bool vbfMedium = false;
    bool vbfLoose = false;
    if(vbfPreselection && mva.mDiJet>700.0 && mva.deltaEtaJets>5.)
    {
        vbfVeryTight=true;
    }
    else if(vbfPreselection && mva.mDiJet>400.0 && mva.deltaEtaJets>5.)
    {
        vbfTight=true;
    }
    else if(vbfPreselection && mva.mDiJet>400.0 && mva.deltaEtaJets>4.)
    {
        vbfMedium=true;
    }
    else if(vbfPreselection && mva.mDiJet>300.0 && mva.deltaEtaJets>3.)
    {
        vbfLoose=true;
    }

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
    if(vbfPreselection && passVBFBDTCut) counterBDTvbfCut++;
    mva.bdtValInc = bdtValInc;
    mva.bdtValVBF = bdtValVBF;

    time_t timeStartFilling = time(NULL);

    if (!blind)
    {
      if(!vbfPreselection)
      {
        hists.BDTHistMuonOnly->Fill(bdtValInc, weight);
        hists.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
  
      }
      else
      {
        hists.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
        hists.BDTHistVBF->Fill(bdtValVBF, weight);
      }
    }

    //4 GeV Window Plots
    if (mva.mDiMu < 127.0 && mva.mDiMu > 123.0)
    {
      if (!blind)
      {
        hists4GeVWindow.Fill(mva,blind);
      }
    }

    if (isBB)
    {
      histsBB.Fill(mva,blind);
    }

    if (isBO)
    {
      histsBO.Fill(mva,blind);
    }

    if (isBE)
    {
      histsBE.Fill(mva,blind);
    }

    if (isOO)
    {
      histsOO.Fill(mva,blind);
    }

    if (isOE)
    {
      histsOE.Fill(mva,blind);
    }

    if (isEE)
    {
      histsEE.Fill(mva,blind);
    }

    if (isNotBB)
    {
      histsNotBB.Fill(mva,blind);
    }

    //VBF Preselected Plots
    if (vbfPreselection)
    {
      histsVBFPresel.Fill(mva,blind);
    }

    if (vbfPreselection && isBB)
    {
      histsVBFPreselBB.Fill(mva,blind);
    }

    if (vbfPreselection && isNotBB)
    {
      histsVBFPreselNotBB.Fill(mva,blind);
    }

    //Inc Preselected Plots
    if (!vbfPreselection)
    {
      histsIncPresel.Fill(mva,blind);
    }

    if (!vbfPreselection && isBB)
    {
      histsIncPreselBB.Fill(mva,blind);
    }

    if (!vbfPreselection && isBO)
    {
      histsIncPreselBO.Fill(mva,blind);
    }

    if (!vbfPreselection && isBE)
    {
      histsIncPreselBE.Fill(mva,blind);
    }

    if (!vbfPreselection && isOO)
    {
      histsIncPreselOO.Fill(mva,blind);
    }

    if (!vbfPreselection && isOE)
    {
      histsIncPreselOE.Fill(mva,blind);
    }

    if (!vbfPreselection && isEE)
    {
      histsIncPreselEE.Fill(mva,blind);
    }

    if (!vbfPreselection && isNotBB)
    {
      histsIncPreselNotBB.Fill(mva,blind);
    }




    //VBF BDT Cut Plots
    if (vbfPreselection && passVBFBDTCut)
    {
      histsVBFBDTCut.Fill(mva,blind);
    }

    if (vbfPreselection && isBB && passVBFBDTCut)
    {
      histsVBFBDTCutBB.Fill(mva,blind);
    }

    if (vbfPreselection && isNotBB && passVBFBDTCut)
    {
      histsVBFBDTCutNotBB.Fill(mva,blind);
    }

    //Inc BDT Cut Plots
    if (!vbfPreselection && passIncBDTCut)
    {
      histsIncBDTCut.Fill(mva,blind);
    }

    if (!vbfPreselection && isBB && passIncBDTCut)
    {
      histsIncBDTCutBB.Fill(mva,blind);
    }

    if (!vbfPreselection && isBO && passIncBDTCut)
    {
      histsIncBDTCutBO.Fill(mva,blind);
    }

    if (!vbfPreselection && isBE && passIncBDTCut)
    {
      histsIncBDTCutBE.Fill(mva,blind);
    }

    if (!vbfPreselection && isOO && passIncBDTCut)
    {
      histsIncBDTCutOO.Fill(mva,blind);
    }

    if (!vbfPreselection && isOE && passIncBDTCut)
    {
      histsIncBDTCutOE.Fill(mva,blind);
    }

    if (!vbfPreselection && isEE && passIncBDTCut)
    {
      histsIncBDTCutEE.Fill(mva,blind);
    }

    if (!vbfPreselection && isNotBB && passIncBDTCut)
    {
      histsIncBDTCutNotBB.Fill(mva,blind);
    }


    if (!vbfPreselection && mva.ptDiMu < 20.0)
    {
      histsIncPreselDiMuPtL20.Fill(mva,blind);
    }

    if (vbfPreselection && mva.ptDiMu < 20.0)
    {
      histsVBFPreselDiMuPtL20.Fill(mva,blind);
    }

    if (vbfPreselection && mva.ptmiss < 50.0)
    {
      histsVBFPreselPtMiss50Veto.Fill(mva,blind);
    }

    if (!vbfPreselection && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10.Fill(mva,blind);
    }

    if (!vbfPreselection && isBB && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10BB.Fill(mva,blind);
    }

    if (!vbfPreselection && isBO && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10BO.Fill(mva,blind);
    }

    if (!vbfPreselection && isBE && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10BE.Fill(mva,blind);
    }

    if (!vbfPreselection && isOO && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10OO.Fill(mva,blind);
    }

    if (!vbfPreselection && isOE && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10OE.Fill(mva,blind);
    }

    if (!vbfPreselection && isEE && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10EE.Fill(mva,blind);
    }

    if (!vbfPreselection && isNotBB && mva.ptDiMu > 10.0)
    {
      histsIncPreselPtG10NotBB.Fill(mva,blind);
    }

///////////////////////////////////////////

    timeReading += difftime(timeStopReading,timeStartReading);
    timeProcessing += difftime(timeStartFilling,timeStopReading);
    timeFilling += difftime(time(NULL),timeStartFilling);
  }// end event loop
  time_t timeEndEventLoop = time(NULL);

//////////////////////////////////////////////////////////////////////
  std::cout << " ########################################## \n";
  std::cout << " Events after: \n\n";
  std::cout << " Gen Mass   Cuts = " << counterGenBoson << std::endl;
  std::cout << " Minimal    Cuts = " << counterMinimCuts << std::endl;
  std::cout << " pT/eta     Cuts = " << counterAccCuts   << std::endl;
  std::cout << " Muon ID    Cuts = " << counterIDCuts    << std::endl;
  std::cout << " Muon Iso   Cuts = " << counterIsoCuts   << std::endl;
  std::cout << " Trigger    Cuts = " << counterTrigCuts  << std::endl;
  std::cout << " VBF pres.  Cuts = " << counterJetSel  << std::endl;
  std::cout << " VBF BDT    Cuts = " << counterBDTvbfCut  << std::endl;
  std::cout << " non VBF pres.  Cuts = " << counterNonJetSel  << std::endl;
  std::cout << " pt(mumu)>10 GeV = " << counterDiPt10GeV  << std::endl;
  // calculate efficiency and binom. error:
  // calculate efficiency and binom. error:
  EffMinimCuts = float(counterMinimCuts)/float(counterGenBoson);
  EffAccCuts = float(counterAccCuts)/float(counterGenBoson);
  EffIDCuts = float(counterIDCuts)/float(counterGenBoson);
  EffIsoCuts = float(counterIsoCuts)/float(counterGenBoson);
  EffTrigCuts = float(counterTrigCuts)/float(counterGenBoson);
  EffJetSel = float(counterJetSel)/float(counterGenBoson);
  EffBDTvbfCut = float(counterBDTvbfCut)/float(counterGenBoson);
  EffNonJetSel = float(counterNonJetSel)/float(counterGenBoson);
  EffDiPt10GeV = float(counterDiPt10GeV)/float(counterGenBoson);

  dEffMinimCuts = sqrt( EffMinimCuts*(1-EffMinimCuts)/float(counterGenBoson) );
  dEffAccCuts = sqrt( EffAccCuts*(1-EffAccCuts)/float(counterGenBoson) );
  dEffIDCuts = sqrt( EffIDCuts*(1-EffIDCuts)/float(counterGenBoson) );
  dEffIsoCuts = sqrt( EffIsoCuts*(1-EffIsoCuts)/float(counterGenBoson) );
  dEffTrigCuts = sqrt( EffTrigCuts*(1-EffTrigCuts)/float(counterGenBoson) );
  dEffJetSel = sqrt( EffJetSel*(1-EffJetSel)/float(counterGenBoson) );
  dEffBDTvbfCut = sqrt( EffBDTvbfCut*(1-EffBDTvbfCut)/float(counterGenBoson) );
  dEffNonJetSel = sqrt( EffNonJetSel*(1-EffNonJetSel)/float(counterGenBoson) );
  dEffDiPt10GeV = sqrt( EffDiPt10GeV*(1-EffDiPt10GeV)/float(counterGenBoson) );
  std::cout << " ########################################## \n";
  std::cout << " Efficiency after selection: \n\n";
  cout.unsetf(ios::floatfield);            // floatfield not set
  cout.precision(3);
  std::cout << " Opposit charge, M_RECO = 110-150 GeV  = " << EffMinimCuts <<" +/- ";
  cout.precision(1);
  std::cout << dEffMinimCuts << std::endl;
  cout.precision(3);
  std::cout << " + pT > 25 GeV,|eta| < 2.1         Cut = " << EffAccCuts   <<" +/- ";
  cout.precision(1);
  std::cout << dEffAccCuts   << std::endl;
  cout.precision(3);
  std::cout << " + Tight Muon ID                   Cut = " << EffIDCuts    <<" +/- " ;
  cout.precision(1);
  std::cout <<  dEffIDCuts    << std::endl;
  cout.precision(3);
  std::cout << " + Muon Relative PF Isolation      Cut = " << EffIsoCuts   <<" +/- " ;
  cout.precision(1);
  std::cout << dEffIsoCuts   << std::endl;
  cout.precision(3);
  std::cout << " + Trigger HLT_Mu24Iso_eta2p1      Cut = " << EffTrigCuts  <<" +/- " ;
  cout.precision(1);
  std::cout << dEffTrigCuts  << std::endl;
  cout.precision(3);
  std::cout << " + VBF Jet preselection            Cut = " << EffJetSel    <<" +/- " ;
  cout.precision(1);
  std::cout << dEffJetSel << std::endl;
  cout.precision(3);
  std::cout << " + VBF BDT                         Cut = " << EffBDTvbfCut    <<" +/- " ;
  cout.precision(1);
  std::cout << dEffBDTvbfCut << std::endl;
  cout.precision(3);
  std::cout << " + non VBF Jet preselection         Cut = " << EffNonJetSel    <<" +/- " ;
  cout.precision(1);
  std::cout << dEffNonJetSel << std::endl;
  cout.precision(3);
  std::cout << " + pt(mumu)> 10 GeV/c, no VBF pres. and no BDT  Cut = " << EffDiPt10GeV    <<" +/- " ;
  cout.precision(1);
  std::cout << dEffDiPt10GeV << std::endl;
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  hists.Write(outFile,"");
  histsBB.Write(outFile,"BB");
  histsBO.Write(outFile,"BO");
  histsBE.Write(outFile,"BE");
  histsOO.Write(outFile,"OO");
  histsOE.Write(outFile,"OE");
  histsEE.Write(outFile,"EE");
  histsNotBB.Write(outFile,"NotBB");
  hists4GeVWindow.Write(outFile,"4GeVWindow");

  histsVBFPresel.Write(outFile,"VBFPresel");
  histsVBFPreselBB.Write(outFile,"VBFPreselBB");
  histsVBFPreselNotBB.Write(outFile,"VBFPreselNotBB");

  histsIncPresel.Write(outFile,"IncPresel");
  histsIncPreselBB.Write(outFile,"IncPreselBB");
  histsIncPreselBO.Write(outFile,"IncPreselBO");
  histsIncPreselBE.Write(outFile,"IncPreselBE");
  histsIncPreselOO.Write(outFile,"IncPreselOO");
  histsIncPreselOE.Write(outFile,"IncPreselOE");
  histsIncPreselEE.Write(outFile,"IncPreselEE");
  histsIncPreselNotBB.Write(outFile,"IncPreselNotBB");

  histsIncBDTCut.Write(outFile,"IncBDTCut");
  histsIncBDTCutBB.Write(outFile,"IncBDTCutBB");
  histsIncBDTCutBO.Write(outFile,"IncBDTCutBO");
  histsIncBDTCutBE.Write(outFile,"IncBDTCutBE");
  histsIncBDTCutOO.Write(outFile,"IncBDTCutOO");
  histsIncBDTCutOE.Write(outFile,"IncBDTCutOE");
  histsIncBDTCutEE.Write(outFile,"IncBDTCutEE");
  histsIncBDTCutNotBB.Write(outFile,"IncBDTCutNotBB");

  histsVBFBDTCut.Write(outFile,"VBFBDTCut");
  histsVBFBDTCutBB.Write(outFile,"VBFBDTCutBB");
  histsVBFBDTCut.Write(outFile,"VBFBDTCutNotBB");

  histsVBFPreselDiMuPtL20.Write(outFile,"VBFPreselDiMuPtL20");
  histsIncPreselDiMuPtL20.Write(outFile,"IncPreselDiMuPtL20");

  histsIncPreselPUJETID.Write(outFile,"IncPreselPUJETID");
  histsVBFPreselPUJETID.Write(outFile,"VBFPreselPUJETID");

  histsIncPreselPUJETIDForVeto.Write(outFile,"IncPreselPUJETIDForVeto");
  histsVBFPreselPUJETIDForVeto.Write(outFile,"VBFPreselPUJETIDForVeto");

  histsVBFPreselPtMiss50Veto.Write(outFile,"VBFPreselPtMiss50Veto");

  histsIncPreselPtG10.Write(outFile,"IncPreselPtG10");
  histsIncPreselPtG10BB.Write(outFile,"IncPreselPtG10BB");
  histsIncPreselPtG10BO.Write(outFile,"IncPreselPtG10BO");
  histsIncPreselPtG10BE.Write(outFile,"IncPreselPtG10BE");
  histsIncPreselPtG10OO.Write(outFile,"IncPreselPtG10OO");
  histsIncPreselPtG10OE.Write(outFile,"IncPreselPtG10OE");
  histsIncPreselPtG10EE.Write(outFile,"IncPreselPtG10EE");
  histsIncPreselPtG10NotBB.Write(outFile,"IncPreselPtG10NotBB");

  ofstream testOutFile;
  testOutFile.open("testEventNums.txt");
  testOutFile << testString;
  testOutFile.close();
  cout <<"#######################\n"<< testString <<"#######################\n" << endl;
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

void
fillMuonHist(TH1F* hist, _MuonInfo& mu1, _MuonInfo& mu2)
{

  hist->Fill(1.0);
  if (!mu1.isGlobal || !mu2.isGlobal)  return;
  hist->Fill(2.0);
  if (!mu1.isPFMuon || !mu2.isPFMuon) return;
  hist->Fill(3.0);

  // acceptance cuts
  if (mu1.pt < 25 || mu2.pt < 25)         return; // pt cut
  hist->Fill(4.0);
  if (fabs(mu1.eta) > 2.1 || fabs(mu2.eta) > 2.1) return; // eta cut
  hist->Fill(5.0);

  // kinematic cuts
  if (mu1.numTrackerLayers < 6 || mu2.numTrackerLayers < 6) return; // # hits in tracker
  hist->Fill(6.0);

  if(getRelIso(mu1) > 0.12 || getRelIso(mu2) > 0.12)
  {
      //cout << "Iso 1: "<< getRelIso(mu1) << "    Iso 2: " << getRelIso(mu2) << endl;
      return;
  }
  hist->Fill(7.0);

  if (fabs(mu1.d0_PV) > 0.2 || fabs(mu2.d0_PV) > 0.2) return;
  hist->Fill(8.0);
  if (fabs(mu1.dz_PV) > 0.5 || fabs(mu2.dz_PV) > 0.5) return;
  hist->Fill(9.0);

  if ( mu1.numValidMuonHits  < 1  || mu2.numValidMuonHits  < 1) return;
  hist->Fill(10.0);
  if ( mu1.numValidPixelHits < 1  || mu2.numValidPixelHits < 1) return;
  hist->Fill(11.0);
  if ( mu1.numOfMatchedStations < 2  || mu2.numOfMatchedStations < 2)
  {
      //cout << "Sta 1: "<<mu1.numOfMatchedStations << "    Sta 2: " << mu2.numOfMatchedStations << endl;
      return;
  }
  hist->Fill(12.0);
  if ( mu1.normChiSquare > 10 || mu2.normChiSquare > 10)     return;
  hist->Fill(13.0);

}

void
printStationMiss(_MuonInfo& mu1, _MuonInfo& mu2, _EventInfo& eventInfo, std::string & testString, unsigned & testCounter)
{

  if (!mu1.isGlobal || !mu2.isGlobal)  return;
  if (!mu1.isPFMuon || !mu2.isPFMuon) return;

  // acceptance cuts
  if (mu1.pt < 25 || mu2.pt < 25)         return; // pt cut
  if (fabs(mu1.eta) > 2.1 || fabs(mu2.eta) > 2.1) return; // eta cut

  // kinematic cuts
  if (mu1.numTrackerLayers < 6 || mu2.numTrackerLayers < 6) return; // # hits in tracker

  if(getRelIso(mu1) > 0.12 || getRelIso(mu2) > 0.12)
  {
      //cout << "Iso 1: "<< getRelIso(mu1) << "    Iso 2: " << getRelIso(mu2) << endl;
      return;
  }

  if (fabs(mu1.d0_PV) > 0.2 || fabs(mu2.d0_PV) > 0.2) return;
  if (fabs(mu1.dz_PV) > 0.5 || fabs(mu2.dz_PV) > 0.5) return;

  if ( mu1.numValidMuonHits  < 1  || mu2.numValidMuonHits  < 1) return;
  if ( mu1.numValidPixelHits < 1  || mu2.numValidPixelHits < 1) return;
  if ( mu1.numOfMatchedStations < 2  || mu2.numOfMatchedStations < 2)
  {
        testCounter++;
        //std::cout <<eventInfo.run <<":"<<eventInfo.event <<"\n"<< std::endl;
        testString.appendAny(eventInfo.run);
        testString.append(":");
        testString.appendAny(eventInfo.event);
        //testString.append("  #  ");
        //testString.appendAny(mu1.numOfMatchedStations);
        //testString.append("  ");
        //testString.appendAny(mu2.numOfMatchedStations);
        testString.append("\n");
      
      return;
  }
  if ( mu1.normChiSquare > 10 || mu2.normChiSquare > 10)     return;

}

HistStruct::HistStruct()
{

  unsigned nMassBins = 800;
  float minMass = 0.;
  float maxMass = 400.;
  unsigned nMVABins = 200;
  mDiMu = new TH1F("mDiMu","DiMuon Mass",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMu);

  mDiMuResSigUp = new TH1F("mDiMuResSigUp","DiMuon Mass Systematic Shift Up: Sigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResSigUp);
  mDiMuResSigDown = new TH1F("mDiMuResSigDown","DiMuon Mass Systematic Shift Down: Sigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResSigDown);
  mDiMuResASigUp = new TH1F("mDiMuResASigUp","DiMuon Mass Systematic Shift Up: ASigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResASigUp);
  mDiMuResASigDown = new TH1F("mDiMuResASigDown","DiMuon Mass Systematic Shift Down: ASigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResASigDown);

  mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);
  histVec.push_back(mDiJet);

  ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  histVec.push_back(ptDiMu);

  ptDiJet = new TH1F("ptDiJet","DiJet Pt",250,0,1000);
  histVec.push_back(ptDiJet);

  yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  histVec.push_back(yDiMu);
  yDiJet = new TH1F("yDiJet","DiJet Rapidity",100,-5,5);
  histVec.push_back(yDiJet);

  yVptDiMu = new TH2F("yVptDiMu","DiMuon Rapidity v. p_{T}",250,0,500,100,0,4);
  histVec2D.push_back(yVptDiMu);
  ptVmDiMu = new TH2F("ptVmDiMu","DiMuon p_{T} v. Mass",nMassBins,minMass,maxMass,250,0,250);
  histVec2D.push_back(ptVmDiMu);
  yVmDiMu = new TH2F("yVmDiMu","DiMuon |y| v. Mass",nMassBins,minMass,maxMass,100,0,4);
  histVec2D.push_back(yVmDiMu);
  phiVmDiMu = new TH2F("phiVmDiMu","DiMuon #phi v. Mass",nMassBins,minMass,maxMass,100,0,3.2);
  histVec2D.push_back(phiVmDiMu);

  ptMu1 = new TH1F("ptMu1","Leading Muon Pt",100,0,200);
  histVec.push_back(ptMu1);
  ptMu2 = new TH1F("ptMu2","Sub-Leading Muon Pt",100,0,200);
  histVec.push_back(ptMu2);
  ptJet1 = new TH1F("ptJet1","Leading Jet Pt",200,0,1000);
  histVec.push_back(ptJet1);
  ptJet2 = new TH1F("ptJet2","Sub-Leading Jet Pt",200,0,1000);
  histVec.push_back(ptJet2);

  etaMu1 = new TH1F("etaMu1","Leading Muon #eta",50,-2.5,2.5);
  histVec.push_back(etaMu1);
  etaMu2 = new TH1F("etaMu2","Sub-Leading Muon #eta",50,-2.5,2.5);
  histVec.push_back(etaMu2);
  etaJet1 = new TH1F("etaJet1","Leading Jet #eta",50,-5.0,5.0);
  histVec.push_back(etaJet1);
  etaJet2 = new TH1F("etaJet2","Sub-Leading Jet #eta",50,-5.0,5.0);
  histVec.push_back(etaJet2);

  deltaEtaJets = new TH1F("deltaEtaJets","#Delta#eta Jets",50,0.0,10.0);
  histVec.push_back(deltaEtaJets);
  deltaPhiJets = new TH1F("deltaPhiJets","#Delta#phi Jets",50,0.0,3.2);
  histVec.push_back(deltaPhiJets);
  deltaRJets = new TH1F("deltaRJets","#Delta R Jets",50,0.0,10.0);
  histVec.push_back(deltaRJets);
  deltaPhiHJ1 = new TH1F("deltaPhiHJ1","#Delta #phi Leading Jet Dimuon",50,0.0,3.2);
  histVec.push_back(deltaPhiHJ1);

  deltaEtaMuons = new TH1F("deltaEtaMuons","#Delta#eta Jets",50,0.0,10.0);
  histVec.push_back(deltaEtaMuons);
  deltaPhiMuons = new TH1F("deltaPhiMuons","#Delta#phi Jets",50,0.0,3.2);
  histVec.push_back(deltaPhiMuons);
  deltaRMuons = new TH1F("deltaRMuons","#Delta R Jets",50,0.0,10.0);
  histVec.push_back(deltaRMuons);

  countsHist = new TH1F("countsHist","Event Counts",10,0.0,10.0);
  countsHist->GetXaxis()->SetBinLabel(1,"total");
  countsHist->GetXaxis()->SetBinLabel(2,"2#mu ID");
  countsHist->GetXaxis()->SetBinLabel(3,"HLT");
  countsHist->GetXaxis()->SetBinLabel(4,"Charge");
  countsHist->GetXaxis()->SetBinLabel(5,"m_{#mu#mu}");
  countsHist->GetXaxis()->SetBinLabel(6,"Inc Pre");
  countsHist->GetXaxis()->SetBinLabel(7,"VBF Pre");
  histVec.push_back(countsHist);

  countsHist2 = new TH1F("countsHist2","Event Counts",14,0.0,14.0);
  countsHist2->GetXaxis()->SetBinLabel(1,"total");
  countsHist2->GetXaxis()->SetBinLabel(2,"total");
  countsHist2->GetXaxis()->SetBinLabel(3,"global");
  countsHist2->GetXaxis()->SetBinLabel(4,"PF");
  countsHist2->GetXaxis()->SetBinLabel(5,"pt");
  countsHist2->GetXaxis()->SetBinLabel(6,"eta");
  countsHist2->GetXaxis()->SetBinLabel(7,"tracker");
  countsHist2->GetXaxis()->SetBinLabel(8,"iso");
  countsHist2->GetXaxis()->SetBinLabel(9,"d0");
  countsHist2->GetXaxis()->SetBinLabel(10,"dz");
  countsHist2->GetXaxis()->SetBinLabel(11,"#mu hits");
  countsHist2->GetXaxis()->SetBinLabel(12,"pixel");
  countsHist2->GetXaxis()->SetBinLabel(13,"stations");
  countsHist2->GetXaxis()->SetBinLabel(14,"#chi^2");
  histVec.push_back(countsHist2);

  cosThetaStar = new TH1F("cosThetaStar","cos(#theta^{*})",50,-1.,1.);
  histVec.push_back(cosThetaStar);
  cosThetaStarCS = new TH1F("cosThetaStarCS","cos(#theta^{*}_{CS})",50,-1.,1.);
  histVec.push_back(cosThetaStarCS);

  puJetIDSimpleDiscJet1 = new TH1F("puJetIDSimpleDiscJet1","PU Jet ID--Simple Discriminator Leading Jet",50,-1.,1.);
  histVec.push_back(puJetIDSimpleDiscJet1);
  puJetIDSimpleDiscJet2 = new TH1F("puJetIDSimpleDiscJet2","PU Jet ID--Simple Discriminator Sub-Leading Jet",50,-1.,1.);
  histVec.push_back(puJetIDSimpleDiscJet2);
  puJetIDSimpleDiscJet3 = new TH1F("puJetIDSimpleDiscJet3","PU Jet ID--Simple Discriminator 3rd Leading Jet",50,-1.,1.);
  histVec.push_back(puJetIDSimpleDiscJet3);

  puJetIDSimpleJet1 = new TH1F("puJetIDSimpleJet1","PU Jet ID--Simple Loose Leading Jet",2,0,2);
  histVec.push_back(puJetIDSimpleJet1);
  puJetIDSimpleJet1->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet1->GetXaxis()->SetBinLabel(2,"Pass");
  puJetIDSimpleJet2 = new TH1F("puJetIDSimpleJet2","PU Jet ID--Simple Loose Sub-Leading Jet",2,-0,2);
  histVec.push_back(puJetIDSimpleJet2);
  puJetIDSimpleJet2->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet2->GetXaxis()->SetBinLabel(2,"Pass");
  puJetIDSimpleJet3 = new TH1F("puJetIDSimpleJet3","PU Jet ID--Simple Loose 3rd Leading Jet",2,-0,2);
  histVec.push_back(puJetIDSimpleJet3);
  puJetIDSimpleJet3->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet3->GetXaxis()->SetBinLabel(2,"Pass");

  BDTHistMuonOnly = new TH1F("BDTHistMuonOnly","BDT Discriminator",nMVABins,-1,1);
  histVec.push_back(BDTHistMuonOnly);

  BDTHistVBF = new TH1F("BDTHistVBF","BDT Discriminator",nMVABins,-1,1);
  histVec.push_back(BDTHistVBF);

  BDTHistMuonOnlyVMass = new TH2F("BDTHistMuonOnlyVMass","BDT Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(BDTHistMuonOnlyVMass);
  BDTHistVBFVMass = new TH2F("BDTHistVBFVMass","BDT Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(BDTHistVBFVMass);

  relIsoMu1 = new TH1F("relIsoMu1","",1000,0,10.0);
  histVec.push_back(relIsoMu1);
  relIsoMu2 = new TH1F("relIsoMu2","",1000,0,10.0);
  histVec.push_back(relIsoMu2);

  nJets = new TH1F("nJets","",11,0,11);
  histVec.push_back(nJets);
  ht = new TH1F("ht","",200,0,2000);
  histVec.push_back(ht);
  nJetsInRapidityGap = new TH1F("nJetsInRapidityGap","",11,0,11);
  histVec.push_back(nJetsInRapidityGap);
  htInRapidityGap = new TH1F("htInRapidityGap","",200,0,2000);
  histVec.push_back(htInRapidityGap);

  nPU = new TH1F("nPU","",100,0,100);
  histVec.push_back(nPU);
  nVtx = new TH1F("nVtx","",100,0,100);
  histVec.push_back(nVtx);
  met = new TH1F("met","",160,0,800);
  histVec.push_back(met);
  ptmiss = new TH1F("ptmiss","",160,0,800);
  histVec.push_back(ptmiss);
  weight = new TH1F("weight","",500,0,5.0);
  histVec.push_back(weight);

  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    (*hist)->Sumw2();
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    (*hist2D)->Sumw2();
}

HistStruct::~HistStruct()
{
  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    delete *hist;
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    delete *hist2D;
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
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    (*hist2D)->Write();

  outfile->cd();
}

void 
HistStruct::Fill(const MVA& mva, bool blind)
{
  yDiMu->Fill(mva.yDiMu, mva.weight);
  ptDiMu->Fill(mva.ptDiMu, mva.weight);
  ptMu1->Fill(mva.ptMu1, mva.weight);
  ptMu2->Fill(mva.ptMu2, mva.weight);
  etaMu1->Fill(mva.etaMu1, mva.weight);
  etaMu2->Fill(mva.etaMu2, mva.weight);
  cosThetaStar->Fill(mva.cosThetaStar, mva.weight);
  cosThetaStarCS->Fill(mva.cosThetaStarCS, mva.weight);
  deltaPhiMuons->Fill(mva.deltaPhiMuons, mva.weight);
  deltaEtaMuons->Fill(mva.deltaEtaMuons, mva.weight);
  deltaRMuons->Fill(mva.deltaRMuons, mva.weight);
  relIsoMu1->Fill(mva.relIsoMu1, mva.weight);
  relIsoMu2->Fill(mva.relIsoMu2, mva.weight);
  nPU->Fill(mva.nPU, mva.weight);
  nVtx->Fill(mva.nVtx, mva.weight);
  met->Fill(mva.met, mva.weight);
  ptmiss->Fill(mva.ptmiss, mva.weight);

  mDiJet->Fill(mva.mDiJet, mva.weight);
  ptDiJet->Fill(mva.ptDiJet, mva.weight);
  yDiJet->Fill(mva.yDiJet, mva.weight);
  ptJet1->Fill(mva.ptJet1, mva.weight);
  ptJet2->Fill(mva.ptJet2, mva.weight);
  etaJet1->Fill(mva.etaJet1, mva.weight);
  etaJet2->Fill(mva.etaJet2, mva.weight);
  deltaEtaJets->Fill(mva.deltaEtaJets, mva.weight);
  deltaPhiJets->Fill(mva.deltaPhiJets, mva.weight);
  deltaRJets->Fill(mva.deltaRJets, mva.weight);
  deltaPhiHJ1->Fill(mva.deltaPhiHJ1, mva.weight);
  nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, mva.weight);
  htInRapidityGap->Fill(mva.htInRapidityGap, mva.weight);
  nJets->Fill(mva.nJets, mva.weight);
  ht->Fill(mva.ht, mva.weight);

  puJetIDSimpleDiscJet1->Fill(mva.puJetIDSimpleDiscJet1,mva.weight);
  puJetIDSimpleDiscJet2->Fill(mva.puJetIDSimpleDiscJet2,mva.weight);
  puJetIDSimpleDiscJet3->Fill(mva.puJetIDSimpleDiscJet3,mva.weight);

  puJetIDSimpleJet1->Fill(mva.puJetIDSimpleJet1,mva.weight);
  puJetIDSimpleJet2->Fill(mva.puJetIDSimpleJet2,mva.weight);
  puJetIDSimpleJet3->Fill(mva.puJetIDSimpleJet3,mva.weight);

  yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu),mva.weight);

  if (!blind)
  {
    mDiMu->Fill(mva.mDiMu, mva.weight);
    mDiMuResSigUp->Fill(mva.mDiMuResSigUp, mva.weight);
    mDiMuResSigDown->Fill(mva.mDiMuResSigDown, mva.weight);
    mDiMuResASigUp->Fill(mva.mDiMuResASigUp, mva.weight);
    mDiMuResASigDown->Fill(mva.mDiMuResASigDown, mva.weight);

    yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu),mva.weight);
    ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu,mva.weight);

    if(!mva.vbfPreselection)
    {
      BDTHistMuonOnly->Fill(mva.bdtValInc, mva.weight);
      BDTHistMuonOnlyVMass->Fill(mva.mDiMu, mva.bdtValInc, mva.weight);
    }
    else
    {
      BDTHistVBF->Fill(mva.bdtValVBF, mva.weight);
      BDTHistVBFVMass->Fill(mva.mDiMu, mva.bdtValVBF, mva.weight);
    }
  }
}
