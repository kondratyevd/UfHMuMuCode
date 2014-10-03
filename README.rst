UFDiMuonsAnalyzer
=================

Creates "Stage 1" Ntuples for UF H->mumu analysis.

- Designed to run on 13 TeV miniAOD data and MC
- Tested on CMSSW_7_0_9 (slc6_amd64_gcc481)

Setup
-----

In your CMSSW "src/" directory, run:

::

  git clone -b CMSSW709_miniAOD https://github.com/jhugon/UfHMuMuCode.git UserArea
  cd UserArea
  nice scram b -j8

Notes
-----

- Not sure if trigger or gen-particle matching works
- Probably should validate that jet selection and cross-cleaning with muons works
- Remember to check that one of the muons passes the trigger.  That may be the only trigger check.
- JEC uncertainty commented out.  Couldn't get the appropriate data out of DB (remember delete statement when you put it back in)
- pre and post FSR muons don't actually mean pre and post FSR.  

  - For Muons: pre means in hard process (status 23).  post means final state, i.e. after Pythia 8 does all of its stuff: UE beam-beam remnant pt boost, ISR, FSR, and anything else it does (status 1).
  - For Bosons: pre means in the hard process (status 22). post means after UE beam-beam remnant pt boost, ISR, and FSR (status 62).

signalGeneration
=================

Cards and scripts to generate H->mumu signal samples
