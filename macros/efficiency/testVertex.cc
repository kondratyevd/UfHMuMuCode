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

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "DataFormats.h"
#include "helpers.h"
#include "mva.h"
#include "LumiReweightingStandAlone.h"

#include "boost/program_options.hpp"
#include "boost/regex.hpp"

#include <limits.h>

#define JETPUID
//#define PUREWEIGHT

using namespace std;
using namespace boost;

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  int maxEvents = std::numeric_limits<int>::max();
  maxEvents = 100;
  cout << "maxEvents = "<< maxEvents << "\n";
  
  TChain * tree = new TChain("tree");

  tree->AddFile("/data/uftrig01b/jhugon/hmumu/privateSignalV3/ggHmumu125.root");

  //////////////////////////
  // Tree Branches
  _PFJetInfo jets;
  //tree->SetBranchAddress("pfJets",&jets);

  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo",&vertexInfo);


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;

  unsigned reportEach=1000;
  if (nEvents/1000>reportEach)
    reportEach = nEvents/1000;
  
  for(unsigned i=0; i<nEvents;i++)
  {
    if(i >= maxEvents)
      continue;
    tree->GetEvent(i);
    if (i % reportEach == 0) cout << "Event: " << i << endl;

    std::cout << "vertexInfo: "
        << "\n  nVertices: " << vertexInfo.nVertices << std::endl;

  }// end event loop

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
  cout << "analyzer done." << endl << endl;
  return 0;
}
