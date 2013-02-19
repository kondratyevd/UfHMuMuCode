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

using namespace std;
using namespace boost;

struct HistStruct
{
  HistStruct(float& pt, float& eta, float& phi, float& ptInv, float& q);
  ~HistStruct();
  void Write(TFile* outfile, std::string directory);

  std::vector<TH1F*> histVec;
  std::vector<TH2F*> histVec2D;

  TH1F* qOverPt;
  TH1F* qOverPtPlus;
  TH1F* qOverPtMinus;

  TH2F* qOverPtVEtaPlus;
  TH2F* qOverPtVPtPlus;
  TH2F* qOverPtVPhiPlus;

  TH2F* qOverPtVEtaMinus;
  TH2F* qOverPtVPtMinus;
  TH2F* qOverPtVPhiMinus;

  TTree* tree;
};

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  gErrorIgnoreLevel = kError;
  time_t timeStart = time(NULL);

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names, put output name first followed by all input file names")
      ("maxEvents,m",program_options::value<int>(), "Maximum Number of Events to Process")
      ("runPeriod,r",program_options::value<string>(), "Running Perios e.g. 7TeV, 8TeV")
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
  float maxMmm = 500.0;
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

  //////////////////////////
  // Histograms

  float pt;
  float eta;
  float phi;
  float ptInv;
  float q;
  HistStruct hists(pt,eta,phi,ptInv,q);


  TRandom3 random(1457);

  //////////////////////////
  //for PU reweighting

#ifdef PUREWEIGHT
  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012ABC.root","pileup","pileup");
  if (runPeriod == "7TeV")
  {
    cout << "Using 2011AB PU reweighting\n";
    lumiWeights = reweight::LumiReWeighting("pileupDists/PileUpHistMCFall11.root","pileupDists/PileUpHist2011AB.root","pileup","pileup");
  }
  else
  {
    cout << "Using 2012ABC PU reweighting\n";
  }
#endif

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
  
    double weight = 1.0;
#ifdef PUREWEIGHT
    if (!isData)
    {
      weight = lumiWeights.weight(nPU);
    }
#endif

    if (!((*muonIdFuncPtr)(reco1)) || !((*muonIdFuncPtr)(reco2)))
          continue;


    if (!isHltMatched(reco1,reco2,allowedHLTPaths) && isData)
        continue;


    if (reco1.charge*reco2.charge != -1)
        continue;


    if (recoCandMass < minMmm || recoCandMass > maxMmm)
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

    hists.qOverPt->Fill(muon1.charge/muon1.pt,weight);
    hists.qOverPt->Fill(muon2.charge/muon2.pt,weight);

    if (muon1.charge>0)
    {
    hists.qOverPtPlus->Fill(1./muon1.pt,weight);
    hists.qOverPtVEtaPlus->Fill(muon1.eta,1./muon1.pt,weight);
    hists.qOverPtVPtPlus->Fill(muon1.pt,1./muon1.pt,weight);
    hists.qOverPtVPhiPlus->Fill(muon1.phi,1./muon1.pt,weight);
    }

    if (muon2.charge>0)
    {
    hists.qOverPtPlus->Fill(1./muon2.pt,weight);
    hists.qOverPtVEtaPlus->Fill(muon2.eta,1./muon2.pt,weight);
    hists.qOverPtVPtPlus->Fill(muon2.pt,1./muon2.pt,weight);
    hists.qOverPtVPhiPlus->Fill(muon2.phi,1./muon2.pt,weight);
    }

    if (muon1.charge<0)
    {
    hists.qOverPtMinus->Fill(1./muon1.pt,weight);
    hists.qOverPtVEtaMinus->Fill(muon1.eta,1./muon1.pt,weight);
    hists.qOverPtVPtMinus->Fill(muon1.pt,1./muon1.pt,weight);
    hists.qOverPtVPhiMinus->Fill(muon1.phi,1./muon1.pt,weight);
    }

    if (muon2.charge<0)
    {
    hists.qOverPtMinus->Fill(1./muon2.pt,weight);
    hists.qOverPtVEtaMinus->Fill(muon2.eta,1./muon2.pt,weight);
    hists.qOverPtVPtMinus->Fill(muon2.pt,1./muon2.pt,weight);
    hists.qOverPtVPhiMinus->Fill(muon2.phi,1./muon2.pt,weight);
    }

    pt = muon1.pt;
    eta = muon1.eta;
    phi = muon1.phi;
    ptInv = 1./muon1.pt;
    q = muon1.charge;
    hists.tree->Fill();
    pt = muon2.pt;
    eta = muon2.eta;
    phi = muon2.phi;
    ptInv = 1./muon2.pt;
    q = muon2.charge;
    hists.tree->Fill();

  }// end event loop
  time_t timeEndEventLoop = time(NULL);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  hists.Write(outFile,"");
  outFile->Close();

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

HistStruct::HistStruct(float& pt, float& eta, float& phi, float& ptInv, float& q)
{

  qOverPt = new TH1F("qOverPt","",200,-0.05,0.05);
  histVec.push_back(qOverPt);

  qOverPtPlus = new TH1F("qOverPtPlus","",100,0.00,0.05);
  histVec.push_back(qOverPtPlus);
  qOverPtMinus = new TH1F("qOverPtMinus","",100,0.00,0.05);
  histVec.push_back(qOverPtMinus);

  qOverPtVEtaPlus = new TH2F("qOverPtVEtaPlus","",10,-2.4,2.4,100,0.0,0.05);
  histVec2D.push_back(qOverPtVEtaPlus);
  qOverPtVEtaMinus = new TH2F("qOverPtVEtaMinus","",10,-2.4,2.4,100,0.0,0.05);
  histVec2D.push_back(qOverPtVEtaMinus);

  qOverPtVPtPlus = new TH2F("qOverPtVPtPlus","",10,0,200,100,0.0,0.05);
  histVec2D.push_back(qOverPtVPtPlus);
  qOverPtVPtMinus = new TH2F("qOverPtVPtMinus","",10,0,200,100,0.0,0.05);
  histVec2D.push_back(qOverPtVPtMinus);

  qOverPtVPhiPlus = new TH2F("qOverPtVPhiPlus","",10,-3.14,3.14,100,0.0,0.05);
  histVec2D.push_back(qOverPtVPhiPlus);
  qOverPtVPhiMinus = new TH2F("qOverPtVPhiMinus","",10,-3.14,3.14,100,0.0,0.05);
  histVec2D.push_back(qOverPtVPhiMinus);

  tree = new TTree("tree","tree");
  tree->Branch("pt",&pt,"pt/F");
  tree->Branch("eta",&eta,"eta/F");
  tree->Branch("phi",&phi,"phi/F");
  tree->Branch("ptInv",&ptInv,"ptInv/F");
  tree->Branch("q",&q,"q/F");

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
  delete tree;
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
  tree->Write();

  outfile->cd();
}
