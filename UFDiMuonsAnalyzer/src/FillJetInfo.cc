
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _JetInfo::init() {

  nJets     = 0;
  nJetsCent = 0;
  nJetsFwd  = 0;
  for (unsigned int i = 0; i < arraySize; i++) {
    px[i]            = -999;
    py[i]            = -999;
    pz[i]            = -999;
    pt[i]            = -999;
    eta[i]           = -999;
    phi[i]           = -999;
    mass[i]          = -999;
    charge[i]        = -999;
    isB[i]           = -999;
    partonFlavour[i] = -999;

    chf[i]  = -999;
    nhf[i]  = -999;
    cef[i]  = -999;
    nef[i]  = -999;
    muf[i]  = -999;
    hfhf[i] = -999;
    hfef[i] = -999;

    cm[i]   = -999;
    chm[i]  = -999;
    nhm[i]  = -999;
    cem[i]  = -999;
    nem[i]  = -999;
    mum[i]  = -999;
    hfhm[i] = -999;
    hfem[i] = -999;

    jecFactor[i] = -999;
    jecUnc[i]    = -999;

    genMatched[i] = -999;
    genPx[i]      = -999;
    genPy[i]      = -999;
    genPz[i]      = -999;
    genPt[i]      = -999;
    genEta[i]     = -999;
    genPhi[i]     = -999;
    genMass[i]    = -999;
    genEMF[i]     = -999;
    genHadF[i]    = -999;
    genInvF[i]    = -999;
    genAuxF[i]    = -999;

    puid[i] = -999;
  }

}

void _JetInfo::fill(const pat::JetCollection jetsSelected, const std::vector<std::string> _btagNames) {

  nJets     = jetsSelected.size();
  nJetsCent = 0;
  nJetsFwd  = 0;

  for (int i = 0; i < nJets; i++) {
    
    if ( i >= int(arraySize) ) {
      std::cout << "Found " << i+1 << "th jet; only " << arraySize << " allowed in array" << std::endl;
      return;
    }

    pat::Jet jet = jetsSelected.at(i);

    if ( fabs( jet.eta() ) < 2.4 ) nJetsCent += 1;
    else                           nJetsFwd  += 1;

    px[i]            = jet.px();
    py[i]            = jet.py();
    pz[i]            = jet.pz();
    pt[i]            = jet.pt();
    eta[i]           = jet.eta();
    phi[i]           = jet.phi();
    mass[i]          = jet.mass();
    charge[i]        = jet.charge();
    partonFlavour[i] = jet.partonFlavour();
    
    chf[i]  = jet.chargedHadronEnergyFraction();
    nhf[i]  = jet.neutralHadronEnergyFraction();
    cef[i]  = jet.chargedEmEnergyFraction();
    nef[i]  = jet.neutralEmEnergyFraction();
    muf[i]  = jet.muonEnergyFraction();
    hfhf[i] = jet.HFHadronEnergyFraction();
    hfef[i] = jet.HFEMEnergyFraction();
    
    cm[i]   = jet.chargedMultiplicity();
    chm[i]  = jet.chargedHadronMultiplicity();
    nhm[i]  = jet.neutralHadronMultiplicity();
    cem[i]  = jet.electronMultiplicity();
    nem[i]  = jet.photonMultiplicity();
    mum[i]  = jet.muonMultiplicity();
    hfhm[i] = jet.HFHadronMultiplicity();
    hfem[i] = jet.HFEMMultiplicity();
    
    jecUnc[i]    = -1.0;
    jecFactor[i] = jet.jecFactor("Uncorrected");
    
    isB[i]  = jet.bDiscriminator(_btagNames[0]);
    puid[i] = jet.userFloat("pileupJetId:fullDiscriminant");
    
    const reco::GenJet* genJet = jet.genJet();
    if (genJet != NULL) {
      genMatched[i] = true;
      genPx[i]      = genJet->px();
      genPy[i]      = genJet->py();
      genPz[i]      = genJet->pz();
      genPt[i]      = genJet->pt();
      genEta[i]     = genJet->eta();
      genPhi[i]     = genJet->phi();
      genMass[i]    = genJet->mass();
      double genJetEnergy = genJet->energy();
      genEMF[i]     = genJet->emEnergy() / genJetEnergy;
      genHadF[i]    = genJet->hadEnergy() / genJetEnergy;
      genInvF[i]    = genJet->invisibleEnergy() / genJetEnergy;
      genAuxF[i]    = genJet->auxiliaryEnergy() / genJetEnergy;
    }
    else {
      genMatched[i] = false;
      genPx[i]      = -1;
      genPy[i]      = -1;
      genPz[i]      = -1;
      genPt[i]      = -1;
      genEta[i]     = -1;
      genPhi[i]     = -1;
      genMass[i]    = -1;
      genEMF[i]     = -1;
      genHadF[i]    = -1;
      genInvF[i]    = -1;
      genAuxF[i]    = -1;
    }
  }  // End loop: for (int i = 0; i < nJets; i++)

}


pat::JetCollection _JetInfo::select( const edm::Handle<pat::JetCollection>& jets, 
				     const JetCorrectorParameters& JetCorPar, const std::string _JES_syst,
				     const std::string _jet_ID, const double _jet_pT_min, const double _jet_eta_max ) {

  // Following https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
  // Modeled after "isGoodJet" in https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/src/MiniAODHelper.cc
  // VBF H-->bb analysis: http://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-16-003
  // Is PU jet ID used at some point in the sequence? (http://cds.cern.ch/record/1581583) Or just charged hadron subtraction (CHS)?
  // B-tagging docmuentation (CSV_v2) here: http://cds.cern.ch/record/2138504

  pat::JetCollection jetsSelected;
  jetsSelected.clear();

  if ( !jets.isValid() ) {
    std::cout << "No valid jet collection" << std::endl;
    return jetsSelected;
  }

  if (_jet_ID.find("loose") == std::string::npos && _jet_ID.find("tight") == std::string::npos)
    std::cout << "Jet ID is neither tight nor loose: " << _jet_ID
              << "\nWill not be used, no jet ID cuts applied" << std::endl;
  
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
  
  for (pat::JetCollection::const_iterator jet = jets->begin(), jetsEnd = jets->end(); jet !=jetsEnd; ++jet) {

    pat::Jet corr_jet = (*jet);
    
    // Apply jet energy scale systematics
    // Modeled after https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
    //           and https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/src/MiniAODHelper.cc#L374
    jecUnc->setJetEta( corr_jet.eta() );
    jecUnc->setJetPt ( corr_jet.pt()  );
    double uncert = jecUnc->getUncertainty(true);
    if      ( _JES_syst.find("JES_up")   != std::string::npos )
      corr_jet.scaleEnergy( 1. + uncert );
    else if ( _JES_syst.find("JES_down") != std::string::npos )
      corr_jet.scaleEnergy( 1. - uncert );
    else if ( _JES_syst.find("none")     == std::string::npos )
      std::cout << "Invalid jet energy systematic " << _JES_syst << ". Leaving uncorrected." << std::endl;
    
    if ( corr_jet.pt()          < _jet_pT_min  ) continue;
    if ( fabs( corr_jet.eta() ) > _jet_eta_max ) continue;
    
    bool isLoose = false;
    bool isTight = false;

    if ( fabs( corr_jet.eta() ) < 2.7 ) {

      isLoose = ( corr_jet.neutralHadronEnergyFraction() < 0.99 &&
		  corr_jet.neutralEmEnergyFraction()     < 0.99 &&
		  corr_jet.numberOfDaughters()           > 1 );

      isTight = ( isLoose                                   && 
		  corr_jet.neutralHadronEnergyFraction() < 0.90 &&
		  corr_jet.neutralEmEnergyFraction()     < 0.90 );

      if ( fabs( corr_jet.eta() ) < 2.4 ) {
	
	isLoose &= ( corr_jet.chargedHadronEnergyFraction() > 0.   &&
		     corr_jet.chargedMultiplicity()         > 0    &&
		     corr_jet.chargedEmEnergyFraction()     < 0.99 );
	
	isTight &= isLoose;
      }
    } // End if ( fabs( corr_jet.eta() ) < 2.7 )
    else if ( fabs( corr_jet.eta() ) < 3.0 ) {

      isLoose = ( corr_jet.neutralEmEnergyFraction() < 0.90 &&
		  corr_jet.neutralMultiplicity()     > 2 );
      isTight = isLoose;
    }
    else {
      isLoose = ( corr_jet.neutralEmEnergyFraction() < 0.90 &&
		  corr_jet.neutralMultiplicity()     > 10 );
      isTight = isLoose;
    }

    if (_jet_ID.find("loose") != std::string::npos && !isLoose) continue;
    if (_jet_ID.find("tight") != std::string::npos && !isTight) continue;
    
    jetsSelected.push_back(corr_jet);
  }
  
  return jetsSelected;
}


