
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/FillGenHelper.h"

void fillBosonAndMuDaughters(const reco::Candidate& boson,
			     _GenPartInfo& _genGpreFSR,  _GenPartInfo& _genM1GpreFSR,  _GenPartInfo& _genM2GpreFSR,
			     _GenPartInfo& _genGpostFSR, _GenPartInfo& _genM1GpostFSR, _GenPartInfo& _genM2GpostFSR,
			     _GenPartInfo& _genZpreFSR,  _GenPartInfo& _genM1ZpreFSR,  _GenPartInfo& _genM2ZpreFSR,
			     _GenPartInfo& _genZpostFSR, _GenPartInfo& _genM1ZpostFSR, _GenPartInfo& _genM2ZpostFSR,
			     _GenPartInfo& _genHpreFSR,  _GenPartInfo& _genM1HpreFSR,  _GenPartInfo& _genM2HpreFSR,
			     _GenPartInfo& _genHpostFSR, _GenPartInfo& _genM1HpostFSR, _GenPartInfo& _genM2HpostFSR,
			     _GenPartInfo& _genWpreFSR,  _GenPartInfo& _genMWpreFSR,
			     _GenPartInfo& _genWpostFSR, _GenPartInfo& _genMWpostFSR ) {
  
  // photon, Z, W, H
  if (boson.status() == 62) {
    if (abs(boson.pdgId()) < 22 || abs(boson.pdgId()) > 25)
      return;
  }
  // Technically 21 is an incoming particle from the feynman diagram q, anti-q -> lept, lept 
  // with no intermediate gen particle. I think this is q, anti-q -> gamma* -> lept, lept.
  // Whenever there is no Z in the dyJetsToLL sample there is this q, anti-q -> lept, lept 
  // where the quark and the antiquark have the same two leptons as daughters.
  // We will have to reconstruct the values for the off shell gamma
  else if (boson.status() != 21)
    return;

  _GenPartInfo bosonInfo;
  _GenPartInfo mu1preFSR;
  _GenPartInfo mu1postFSR;
  _GenPartInfo mu2preFSR;
  _GenPartInfo mu2postFSR;

  bosonInfo.init();
  mu1preFSR.init();
  mu1postFSR.init();
  mu2preFSR.init();
  mu2postFSR.init();

  bosonInfo.fill(boson);

  TLorentzVector l1, l2, mother;
  bool moreThanOneLeptPair = false;

  // Get the daughter muons for the boson
  for (unsigned int i = 0; i < boson.numberOfDaughters(); i++) {
    const reco::Candidate* daughter = boson.daughter(i);
    
    // Get information about the lepton daughters to reconstruct the virtual photon later
    if (daughter->pdgId() == 11 || daughter->pdgId() == 13 || daughter->pdgId() == 15) {
      for (unsigned int j = 0; j < boson.numberOfDaughters(); j++) { // Loop from i? - AWB 08.11.16
	// We already know you can't make a lepton pair with yourself, so skip this one if it's the case
	if(i == j) continue;
	const reco::Candidate* daughter2 = boson.daughter(j);
	// Found a pair of opposite signed lepton daughters
	if(daughter->pdgId() == -1*daughter2->pdgId()) {
	  // We already found a pair, l1 AND l2 were initialized already
	  if (l1.M() != 0 && l2.M() != 0) moreThanOneLeptPair = true;
	  l1.SetPtEtaPhiM(daughter->pt(),  daughter->eta(),  daughter->phi(),  daughter->mass());
	  l2.SetPtEtaPhiM(daughter2->pt(), daughter2->eta(), daughter2->phi(), daughter2->mass());
	}
      }
    }

    // Status 23 muon, intermediate particle from a decay
    if (daughter->pdgId() == 13 && daughter->status() == 23) {
      // We have an intermediate status 23 muon, save the intermediate values as preFSR
      mu1preFSR.fill( *daughter );
      
      // If it did not radiate then the post and pre are the same
      mu1postFSR = mu1preFSR;

      // If it did radiate, get the postFSR final state, status 1 version 
      // of this muon and overwrite the postFSR quantities
      for (unsigned int i = 0; i < daughter->numberOfDaughters(); i++) {
	const reco::Candidate* postFSRcand = daughter->daughter(i);
	if (postFSRcand->pdgId() == 13 && daughter->status() == 1)
	  mu1postFSR.fill( *postFSRcand );
      }
    }

    // Status 23 antimuon, intermediate particle from a decay
    else if (daughter->pdgId() == -13 && daughter->status() == 23) {
      // We have an intermediate status 23 muon, save the intermediate values as preFSR                                                        
      mu2preFSR.fill( *daughter );
      
      mu2postFSR = mu2preFSR;
      
      // If it did radiate, get the postFSR final state, status 1 version 
      // of this muon and overwrite the postFSR quantities
      for (unsigned int i = 0; i<daughter->numberOfDaughters(); i++) {
	const reco::Candidate* postFSRcand = daughter->daughter(i);
	if (postFSRcand->pdgId() == -13 && daughter->status() == 1)
	  mu2postFSR.fill( *postFSRcand );
      }
    }

    // Final state muon
    else if (daughter->pdgId() == 13 && daughter->status() == 1) {
      // No intermediate status 23 muon that radiated only final state status 1, so pre and post are the same
      mu1preFSR.fill( *daughter );

      // No radiation, post and pre are the same                                                                                               
      mu1postFSR = mu1preFSR;
    }
    // Final state antimuon
    else if (daughter->pdgId() == -13 && daughter->status() == 1) {
      // No intermediate status 23 muon that radiated only final state status 1, so pre and post are the same                                  
      mu2preFSR.fill( *daughter );

      // No radiation, post and pre are the same
      mu2postFSR = mu2preFSR;
    }

  } // End loop: for (unsigned int i = 0; i < boson.numberOfDaughters(); i++)

  // No muons found, Z/gamma* decayed to other leptons
  if ( mu1preFSR.pt < 0 && mu1postFSR.pt < 0 && 
       mu2preFSR.pt < 0 && mu2postFSR.pt < 0 && 
       l1.M() !=0 && l2.M() != 0 ) {
    // Use the mass to let us know what the Z/gamma* decayed to
    mu1preFSR.pt     = -l1.M();
    mu1preFSR.eta    = -l1.M();
    mu1preFSR.phi    = -l1.M();;
    mu1preFSR.charge = -l1.M();
    
    mu2preFSR.pt     = -l2.M();
    mu2preFSR.eta    = -l2.M();
    mu2preFSR.phi    = -l2.M();;
    mu2preFSR.charge = -l2.M();
    
    mu1postFSR = mu1preFSR;
    mu2postFSR = mu2preFSR;
  }

  // Fill the appropriate boson and daughters
  if (boson.status() == 21) {
    // In DY this is an incoming quark, anti-quark annihilation
    // 21 could be any incoming particle in other samples though
    
    // If the virtual photon went to two leptons then reconstruct it
    if (l1.M() != 0 && l2.M() != 0) {
      mother = l1 + l2;

      bosonInfo.mass   = mother.M();
      bosonInfo.pt     = mother.Pt();
      bosonInfo.eta    = -111;
      bosonInfo.y      = mother.Rapidity();
      bosonInfo.phi    = mother.Phi();
      bosonInfo.charge = 0;

      // Not sure what to do if the virtual photon decayed to a bunch of leptons
      if (moreThanOneLeptPair) {
	bosonInfo.mass   = -333;
	bosonInfo.pt     = -333;
	bosonInfo.eta    = -333;
	bosonInfo.y      = -333;
	bosonInfo.phi    = -333;
	bosonInfo.charge = -333;
      }
    }
    else {
      bosonInfo.init();
    }

    _genGpreFSR = bosonInfo;
    _genM1GpreFSR = mu1preFSR;
    _genM2GpreFSR = mu2preFSR;
    
    _genGpostFSR = bosonInfo;
    _genM1GpostFSR = mu1postFSR;
    _genM2GpostFSR = mu2postFSR;

  } // End conditional: if (boson.status() == 21)

  // Z
  if(abs(boson.pdgId()) == 23) {
    _genZpreFSR = bosonInfo;
    _genM1ZpreFSR = mu1preFSR;
    _genM2ZpreFSR = mu2preFSR;

    _genZpostFSR = bosonInfo;
    _genM1ZpostFSR = mu1postFSR;
    _genM2ZpostFSR = mu2postFSR;
  }
  // W
  if(abs(boson.pdgId()) == 24) {
    _genWpreFSR = bosonInfo;
    if (bosonInfo.charge == mu1preFSR.charge) _genMWpreFSR = mu1preFSR;
    if (bosonInfo.charge == mu2preFSR.charge) _genMWpreFSR = mu2preFSR;
    
    _genWpostFSR = bosonInfo;
    if (bosonInfo.charge == mu1postFSR.charge) _genMWpreFSR = mu1postFSR;
    if (bosonInfo.charge == mu2postFSR.charge) _genMWpreFSR = mu2postFSR;
  }

  // H
  if(abs(boson.pdgId()) == 25) {
    _genHpreFSR = bosonInfo;
    _genM1HpreFSR = mu1preFSR;
    _genM2HpreFSR = mu2preFSR;
    
    _genHpostFSR = bosonInfo;
    _genM1HpostFSR = mu1postFSR;
    _genM2HpostFSR = mu2postFSR;
  }

} // End void FillGenHelper::fillBosonAndMuDaughters
