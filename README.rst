UFDiMuonsAnalyzer
=================

Creates "Stage 1" Ntuples for UF H->mumu analysis.

- Designed to run on 13 TeV miniAOD data and MC
- Tested on CMSSW_7_0_9 (slc6_amd64_gcc481)

Setup
-----

You must use SLC6 with CMSSW_7_X.  The HPC, lxplus, and melrose.ihepa.ufl.edu all run SLC6.

In Florida, put the following in your ``.bashrc`` to make CMSSW work:

::

  export SCRAM_ARCH=slc6_amd64_gcc481
  export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
  . $VO_CMS_SW_DIR/cmsset_default.sh
  . $VO_CMS_SW_DIR/crab3/crab.sh


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

