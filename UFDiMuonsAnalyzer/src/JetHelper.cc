
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/JetHelper.h"

void FillJetInfos( JetInfos& _jetInfos, int& _nJetsFwd,
		   const pat::JetCollection jetsSelected, 
		   const std::vector<std::string> _btagNames) {

  _jetInfos.clear();
  _nJetsFwd  = 0;
  int nJets     = jetsSelected.size();
  int nJetsCent = 0;

  for (int i = 0; i < nJets; i++) {
    
    pat::Jet jet = jetsSelected.at(i);
    JetInfo _jetInfo;
    _jetInfo.init();

    _jetInfo.px            = jet.px();
    _jetInfo.py            = jet.py();
    _jetInfo.pz            = jet.pz();
    _jetInfo.pt            = jet.pt();
    _jetInfo.eta           = jet.eta();
    _jetInfo.phi           = jet.phi();
    _jetInfo.mass          = jet.mass();
    _jetInfo.charge        = jet.charge();
    _jetInfo.partonFlavour = jet.partonFlavour();
    
    _jetInfo.chf  = jet.chargedHadronEnergyFraction();
    _jetInfo.nhf  = jet.neutralHadronEnergyFraction();
    _jetInfo.cef  = jet.chargedEmEnergyFraction();
    _jetInfo.nef  = jet.neutralEmEnergyFraction();
    _jetInfo.muf  = jet.muonEnergyFraction();
    _jetInfo.hfhf = jet.HFHadronEnergyFraction();
    _jetInfo.hfef = jet.HFEMEnergyFraction();
    
    _jetInfo.cm   = jet.chargedMultiplicity();
    _jetInfo.chm  = jet.chargedHadronMultiplicity();
    _jetInfo.nhm  = jet.neutralHadronMultiplicity();
    _jetInfo.cem  = jet.electronMultiplicity();
    _jetInfo.nem  = jet.photonMultiplicity();
    _jetInfo.mum  = jet.muonMultiplicity();
    _jetInfo.hfhm = jet.HFHadronMultiplicity();
    _jetInfo.hfem = jet.HFEMMultiplicity();
    
    _jetInfo.jecUnc    = -1.0;
    _jetInfo.jecFactor = jet.jecFactor("Uncorrected");
    
    _jetInfo.isB  = jet.bDiscriminator(_btagNames[0]);
    _jetInfo.puid = jet.userFloat("pileupJetId:fullDiscriminant");
    
    const reco::GenJet* genJet = jet.genJet();
    if (genJet != NULL) {
      _jetInfo.genMatched = true;
      _jetInfo.genPx      = genJet->px();
      _jetInfo.genPy      = genJet->py();
      _jetInfo.genPz      = genJet->pz();
      _jetInfo.genPt      = genJet->pt();
      _jetInfo.genEta     = genJet->eta();
      _jetInfo.genPhi     = genJet->phi();
      _jetInfo.genMass    = genJet->mass();
      double genJetEnergy = genJet->energy();
      _jetInfo.genEMF     = genJet->emEnergy() / genJetEnergy;
      _jetInfo.genHadF    = genJet->hadEnergy() / genJetEnergy;
      _jetInfo.genInvF    = genJet->invisibleEnergy() / genJetEnergy;
      _jetInfo.genAuxF    = genJet->auxiliaryEnergy() / genJetEnergy;
    }
    else {
      _jetInfo.genMatched = false;
      _jetInfo.genPx      = -1;
      _jetInfo.genPy      = -1;
      _jetInfo.genPz      = -1;
      _jetInfo.genPt      = -1;
      _jetInfo.genEta     = -1;
      _jetInfo.genPhi     = -1;
      _jetInfo.genMass    = -1;
      _jetInfo.genEMF     = -1;
      _jetInfo.genHadF    = -1;
      _jetInfo.genInvF    = -1;
      _jetInfo.genAuxF    = -1;
    }
    
    if ( fabs( jet.eta() ) < 2.4 ) nJetsCent += 1;
    else                           _nJetsFwd += 1;
    _jetInfos.push_back( _jetInfo );

  }  // End loop: for (int i = 0; i < nJets; i++)

  if ( nJetsCent + _nJetsFwd != int(_jetInfos.size()) )
    std::cout << "Bizzare error: _jetInfos.size() = " << _jetInfos.size()
	      << ", nJetsCent = " << nJetsCent
	      << ", _nJetsFwd = " << _nJetsFwd << std::endl;
  
}


pat::JetCollection SelectJets( const edm::Handle<pat::JetCollection>& jets, 
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
