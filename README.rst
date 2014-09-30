UfHMuMuCode
============

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

Remember to check that one of the muons passes the trigger.  That may be the only trigger check.
