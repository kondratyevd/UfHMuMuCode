UfHMuMuCode
============

Creates "Stage 1" Ntuples for UF H->mumu analysis.

- Designed to run on 13 TeV miniAOD data and MC
- Tested on CMSSW_7_0_9 (slc6_amd64_gcc481)

Setup
-----

In your CMSSW "src/" directory, run:

::

  git clone https://github.com/jhugon/UfHMuMuCode.git UserArea
  cd UserArea
  git checkout CMSSW709_miniAOD
  nice scram b -j8

