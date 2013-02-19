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

#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>

#include <limits.h>

#define JETPUID
#define PUREWEIGHT
#define VERBOSE

using namespace std;
using namespace boost;

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names, put output name first followed by all input file names")
      ("eventStr,e",program_options::value<vector<string> >(), "Run:Event Strings")
      ("eventFile,f",program_options::value<string>(), "File with one run:event per line")
      ("maxEvents,m",program_options::value<unsigned>(), "Max Events to Print")
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

  vector<std::string> filenames;
  vector<string>::const_iterator filename;
  if (optionMap.count("filenames")>0)
  {
     filenames = optionMap["filenames"].as<vector<string> >();
     if(filenames.size()<1)
     {
       cout << "Error: Need input file names, exiting." << endl;
       return 1;
     }
  }
  else
  {
     cout << "Error: Input file name arguments required, exiting." << endl;
     return 1;
  }

  unsigned nEventsDesired = 0;
  vector<std::string> eventStrs;
  vector<string>::const_iterator eventStr;
  if (optionMap.count("eventStr")>0)
  {
     eventStrs = optionMap["eventStr"].as<vector<string> >();
     if(eventStrs.size()<1)
     {
       cout << "Error: No Wanted events found, use -e <run>:<event> option." << endl;
       return 1;
     }
  }
  else if(optionMap.count("eventFile"))
  {
     string eventFn = optionMap["eventFile"].as<string>();
     ifstream eventF(eventFn.c_str());
     string tmpLine;
     while(eventF.good())
     {
      using namespace boost;
      getline(eventF,tmpLine);
      trim_left_if(tmpLine,is_any_of(" \t\r\n"));
      trim_right_if(tmpLine,is_any_of(" \t\r\n"));
      if (tmpLine.size() <1)
        continue;
      if (tmpLine[0] =='#')
        continue;
      eventStrs.push_back(tmpLine);
     }
     eventF.close();
  }
  else if(optionMap.count("maxEvents"))
  {
     nEventsDesired = optionMap["maxEvents"].as<unsigned>();
  }
  else
  {
     cout << "Error: No Wanted events found, use -e <run>:<event> option." << endl;
     return 1;
  }

  vector<unsigned> runs;
  vector<unsigned> events;
  if (nEventsDesired==0)
  {
    nEventsDesired = eventStrs.size();
    runs = vector<unsigned>(eventStrs.size());
    events = vector<unsigned>(eventStrs.size());
    unsigned tmpMissed = 0;
    for (unsigned i=0;i<nEventsDesired;i++)
    {
      std::vector<std::string> splitStrings;
      boost::split(splitStrings, eventStrs[i], boost::is_any_of(":"));
      if(splitStrings.size() != 2)
      {
        tmpMissed++;
        cerr << "Error: Didn't Correctly Parse Event String: "<< eventStrs[i] << endl;
      }
      runs[i] = boost::lexical_cast<unsigned>(splitStrings[0]);
      events[i] = boost::lexical_cast<unsigned>(splitStrings[1]);
    }
    nEventsDesired -= tmpMissed;
  }

  ofstream outTxt;
  outTxt.open("testEventPrints.txt");

  /////////////////////////////
  //////////// Setup //////////
  /////////////////////////////

  std::vector<int> allowedHLTPaths;
  allowedHLTPaths.push_back(0); //IsoMu24_v11
  allowedHLTPaths.push_back(2); //IsoMu24_v12

  ////////////
  
  TChain * tree = new TChain("tree");

  cout << "Input File Names: \n"; 
  for(filename = filenames.begin();filename != filenames.end();filename++)
  {
    cout<<"  "<< *filename << endl;
    tree->AddFile(filename->c_str());
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
  if(tree->GetBranchStatus("trueMass"))
    tree->SetBranchAddress("trueMass", &trueMass);

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
  if(tree->GetBranchStatus("nPU"))
    tree->SetBranchAddress("nPU",&nPU);
#endif
  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo",&vertexInfo);
  _EventInfo eventInfo;
  //TBranch* eventInfoBranch = tree->GetBranch("eventInfo");
  //eventInfoBranch->SetAddress(&eventInfo);
  tree->SetBranchAddress("eventInfo",&eventInfo);
  _MetInfo met;
  tree->SetBranchAddress("met",&met);

cout.precision(5);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;

  unsigned reportEach=1000;
  if (nEvents/1000>reportEach)
    reportEach = nEvents/1000;

  unsigned foundCounter = 0;
  
  for(unsigned i=0; i<nEvents;i++)
  {
    tree->GetEvent(i);
    if(events.size()!=0)
    {
      bool wantedEvent = false;
      for(unsigned j=0;j<nEventsDesired;j++)
      {
        if (eventInfo.run==runs[j] && eventInfo.event==events[j])
        {
          wantedEvent=true;
          break;
        }
      }
      if(!wantedEvent)
        continue;
    }

    _MuonInfo muon1=reco1;
    _MuonInfo muon2=reco2;
    if(reco1.pt<reco2.pt)
    {
        muon1 = reco2;
        muon2 = reco1;
    }

    unsigned nJets = 0;
    for(unsigned j=0;j<10;j++)
    {
      if (jets.pt[j]<30.0)
      {
        break;
      }
      nJets++;
    }

    cout << "Run:Event (Lumi) = "<<eventInfo.run<<":"<<eventInfo.event<<" ("<<eventInfo.lumi<<")\n"
         << "  Dimuon M = "<<recoCandMass<<" \t Pt = "<<recoCandPt << " \t Y = "<<recoCandY <<"\n"
         << "   Lead Mu Pt = "<< muon1.pt << " \t Eta = "<<muon1.eta << " \t Phi = "<<muon1.phi << " \t Q = "<<muon1.charge <<"\n"
         << "   2nd Mu Pt  = "<< muon2.pt << " \t Eta = "<<muon2.eta << " \t Phi = "<<muon2.phi << " \t Q = "<<muon2.charge <<"\n"
         << "     Nvtx = "<< vertexInfo.nVertices <<" \t nJets = "<<nJets<<"\n"
#ifdef VERBOSE
         << "   Lead Mu Trg = "<< isHltMatched(muon1,allowedHLTPaths) << " \t 2nd Mu Trg = "<< isHltMatched(muon2,allowedHLTPaths) <<"\n"
         << "   Lead Mu ID = "<< isKinTight_2012(muon1) << "\n"
         << "     isGlobal = "<< muon1.isGlobal << " \t isPFMuon = "<<muon1.isPFMuon <<"\n"
         << "     numTrackerLayers = "<< muon1.numTrackerLayers << "\n"
         << "     Iso = "<< getRelIso(muon1) << "\n"
         << "     |d0| = "<< fabs(muon1.d0_PV) << " \t |dz| = " << fabs(muon1.dz_PV) << "\n"
         << "     MuHits = "<< muon1.numValidMuonHits << " \t PixHits = "<<muon1.numValidPixelHits<<"\n"
         << "     Stations = "<< muon1.numOfMatchedStations << " \t chi^2/ndf = "<<muon1.normChiSquare<<"\n"
         << "   2nd Mu ID  = "<< isKinTight_2012(muon2) << "\n"
         << "     isGlobal = "<< muon2.isGlobal << " \t isPFMuon = "<<muon2.isPFMuon <<"\n"
         << "     numTrackerLayers = "<< muon2.numTrackerLayers << "\n"
         << "     Iso = "<< getRelIso(muon2) << "\n"
         << "     |d0| = "<< fabs(muon2.d0_PV) << " \t |dz| = " << fabs(muon2.dz_PV) << "\n"
         << "     MuHits = "<< muon2.numValidMuonHits << " \t PixHits = "<<muon2.numValidPixelHits<<"\n"
         << "     Stations = "<< muon2.numOfMatchedStations << " \t chi^2/ndf = "<<muon2.normChiSquare<<"\n"
         << "   Lead Mu Iso = "<< getRelIso(muon1) << " \t 2nd Mu Iso = "<< getRelIso(muon2) <<"\n"
#endif
         << "   Lead Jet Pt = "<< jets.pt[0] << " \t Eta = "<<jets.eta[0] << " \t Phi = "<<jets.phi[0] << " \t CSV = "<<-999 <<"\n"
         << "   2nd Jet Pt  = "<< jets.pt[1] << " \t Eta = "<<jets.eta[1] << " \t Phi = "<<jets.phi[1] << " \t CSV = "<<-999 <<"\n"
         << "   MET  = "<< met.pt << " \t MET Phi = "<<met.phi <<"\n"
         << endl;

  outTxt << "Run:Event (Lumi) = "<<eventInfo.run<<":"<<eventInfo.event<<" ("<<eventInfo.lumi<<")\n"
         << "  Dimuon M = "<<recoCandMass<<" \t Pt = "<<recoCandPt << " \t Y = "<<recoCandY <<"\n"
         << "   Lead Mu Pt = "<< muon1.pt << " \t Eta = "<<muon1.eta << " \t Phi = "<<muon1.phi << " \t Q = "<<muon1.charge <<"\n"
         << "   2nd Mu Pt  = "<< muon2.pt << " \t Eta = "<<muon2.eta << " \t Phi = "<<muon2.phi << " \t Q = "<<muon2.charge <<"\n"
         << "     Nvtx = "<< vertexInfo.nVertices <<" \t nJets = "<<nJets<<"\n"
#ifdef VERBOSE
         << "   Lead Mu Trg = "<< isHltMatched(muon1,allowedHLTPaths) << " \t 2nd Mu Trg = "<< isHltMatched(muon2,allowedHLTPaths) <<"\n"
         << "   Lead Mu ID = "<< isKinTight_2012(muon1) << "\n"
         << "     isGlobal = "<< muon1.isGlobal << " \t isPFMuon = "<<muon1.isPFMuon <<"\n"
         << "     numTrackerLayers = "<< muon1.numTrackerLayers << "\n"
         << "     Iso = "<< getRelIso(muon1) << "\n"
         << "     |d0| = "<< fabs(muon1.d0_PV) << " \t |dz| = " << fabs(muon1.dz_PV) << "\n"
         << "     MuHits = "<< muon1.numValidMuonHits << " \t PixHits = "<<muon1.numValidPixelHits<<"\n"
         << "     Stations = "<< muon1.numOfMatchedStations << " \t chi^2/ndf = "<<muon1.normChiSquare<<"\n"
         << "   2nd Mu ID  = "<< isKinTight_2012(muon2) << "\n"
         << "     isGlobal = "<< muon2.isGlobal << " \t isPFMuon = "<<muon2.isPFMuon <<"\n"
         << "     numTrackerLayers = "<< muon2.numTrackerLayers << "\n"
         << "     Iso = "<< getRelIso(muon2) << "\n"
         << "     |d0| = "<< fabs(muon2.d0_PV) << " \t |dz| = " << fabs(muon2.dz_PV) << "\n"
         << "     MuHits = "<< muon2.numValidMuonHits << " \t PixHits = "<<muon2.numValidPixelHits<<"\n"
         << "     Stations = "<< muon2.numOfMatchedStations << " \t chi^2/ndf = "<<muon2.normChiSquare<<"\n"
         << "   Lead Mu Iso = "<< getRelIso(muon1) << " \t 2nd Mu Iso = "<< getRelIso(muon2) <<"\n"
#endif
         << "   Lead Jet Pt = "<< jets.pt[0] << " \t Eta = "<<jets.eta[0] << " \t Phi = "<<jets.phi[0] << " \t CSV = "<<-999 <<"\n"
         << "   2nd Jet Pt  = "<< jets.pt[1] << " \t Eta = "<<jets.eta[1] << " \t Phi = "<<jets.phi[1] << " \t CSV = "<<-999 <<"\n"
         << "   MET  = "<< met.pt << " \t MET Phi = "<<met.phi<<"\n"
         << endl;

    foundCounter++;
    if(foundCounter >= nEventsDesired)
        break;
  }// end event loop

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  outTxt.close();
  cout << "Found "<<foundCounter<<"/"<<nEventsDesired << endl;
  cout << "done." << endl << endl;
  return 0;
}
