
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/MetHelper.h"

void FillMetInfo( MetInfo& _metInfo,
		  const edm::Handle < pat::METCollection >& mets,
		  const edm::Event& iEvent ) {

  _metInfo.init();
  if ( !mets.isValid() ) {
    std::cout << "No valid MET collection" << std::endl;
    return;
  }

  for (pat::METCollection::const_iterator met = mets->begin(),
         metsEnd = mets->end(); met !=metsEnd; ++met) {

    // // Filter list from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
    // // Incomplete! - AWB 14.11.16
    // if ( !iEvent->HBHENoiseFilter()                    ) break;
    // if ( !iEvent->HBHENoiseIsoFilter()                 ) break;
    // if ( !iEvent->EcalDeadCellTriggerPrimitiveFilter() ) break;
    // if ( !iEvent->goodVertices()                       ) break;
    // if ( !iEvent->eeBadScFilter()                      ) break;
    // if ( !iEvent->globalTightHalo2016Filter()          ) break;
    
    _metInfo.px    = met->px();
    _metInfo.py    = met->py();
    _metInfo.pt    = met->pt();
    _metInfo.phi   = met->phi();
    _metInfo.sumEt = met->sumEt();

    break;  // Only use the first "MET" in vector of "METs"
    
  }
  
}
