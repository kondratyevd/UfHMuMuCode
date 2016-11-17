
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _MetInfo::init() {
  px    = -999;
  py    = -999;
  pt    = -999;
  phi   = -999;
  sumEt = -999;
}

void _MetInfo::fill( const edm::Handle < pat::METCollection >& mets,
		     const edm::Event& iEvent ) {

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
    
    px    = met->px();
    py    = met->py();
    pt    = met->pt();
    phi   = met->phi();
    sumEt = met->sumEt();

    break;  // Only use the first "MET" in vector of "METs"
    
  }
  
}
