
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/MuonHelper.h"

void FillMuonInfos( MuonInfos& _muonInfos, 
		    const pat::MuonCollection muonsSelected,
		    const reco::Vertex primaryVertex, const int _nPV,
		    const edm::Handle<reco::BeamSpot>& beamSpotHandle,  
		    const edm::Event& iEvent, const edm::EventSetup& iSetup,
		    const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
		    const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
		    const std::vector<std::string> _triggerNames, const double _muon_trig_dR, 
		    const bool _muon_use_pfIso, const double _muon_iso_dR, const bool _isData,
		    KalmanMuonCalibrator& _KaMu_calib, const bool _doSys_KaMu,
		    rochcor2016* _Roch_calib[201], const bool _doSys_Roch ) {
  
  double const MASS_MUON  = 0.105658367; // GeV/c^2

  _muonInfos.clear();
  int nMuons = muonsSelected.size();

  for (int i = 0; i < nMuons; i++) {

    pat::Muon muon = muonsSelected.at(i);
    MuonInfo _muonInfo;
    _muonInfo.init();

    reco::Track track;
    // Do we want to use the inner tracker track by default for the Kalman corrections? - AWB 12.11.16
    if      ( muon.isGlobalMuon()  ) track = *(muon.globalTrack() );
    else if ( muon.isTrackerMuon() ) track = *(muon.innerTrack()  );
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return;
    }
    
    // Basic kinematics
    _muonInfo.charge = muon.charge();
    _muonInfo.pt     = muon.pt();
    _muonInfo.ptErr  = track.ptError();
    _muonInfo.eta    = muon.eta();
    _muonInfo.phi    = muon.phi();
    
    // Basic quality
    _muonInfo.isTightID  =  MuonIsTight  ( muon, primaryVertex );
    _muonInfo.isMediumID =  MuonIsMedium ( muon );
    _muonInfo.isLooseID  =  MuonIsLoose  ( muon );
    
    // Basic isolation
    if ( _muon_use_pfIso )
      _muonInfo.relIso = MuonCalcRelIsoPF ( muon, _muon_iso_dR );
    else
      _muonInfo.relIso = MuonCalcRelIsoTrk( muon, _muon_iso_dR );
    
    // Basic vertexing? - AWB 12.11.16

    // // Basic trigger (which triggers are in the list? - AWB 12.11.16)
    // for (unsigned int iTrigger = 0; iTrigger < _triggerNames.size(); iTrigger++) {
    //   if (iTrigger >= triggerArraySize) {
    // 	std::cout << "Found " << i+1 << "th muon trigger; only " << triggerArraySize << " allowed in array" << std::endl;
    // 	continue;
    //   }
    //   _muonInfo.isHltMatched[iTrigger] = IsHltMatched( muon, iEvent, iSetup, _triggerObjsHandle, _triggerResultsHandle,  
    // 							  _triggerNames[iTrigger], _muon_trig_dR );
    //   _muonInfo.hltEff[iTrigger] = MuonCalcTrigEff( muon, _nPV, _triggerNames[iTrigger] );
    // }


    // // -------------------------
    // // ONLY SAVE BASIC VARIABLES - AWB 12.11.16
    // // -------------------------
    // if ( _outputLevel == "Slim") 
    //   continue;

    // Inner track kinematics
    if (muon.isTrackerMuon()) {
      _muonInfo.trkPt    = muon.innerTrack()->pt();
      _muonInfo.trkPtErr = muon.innerTrack()->ptError();
      _muonInfo.trkEta   = muon.innerTrack()->eta();
      _muonInfo.trkPhi   = muon.innerTrack()->phi();

      // Kalman-calibrated pT
      double pt_KaMu    = muon.innerTrack()->pt();
      double ptErr_KaMu = muon.innerTrack()->ptError();
      double eta_KaMu   = muon.innerTrack()->eta();
      double phi_KaMu   = muon.innerTrack()->phi();
      int charge_KaMu   = muon.innerTrack()->charge();
      double pt_sys_up_KaMu, pt_sys_down_KaMu, pt_clos_up_KaMu, pt_clos_down_KaMu;

      CorrectPtKaMu( _KaMu_calib, _doSys_KaMu, 
      		     pt_KaMu, ptErr_KaMu, pt_sys_up_KaMu, 
      		     pt_sys_down_KaMu, pt_clos_up_KaMu, pt_clos_down_KaMu,
      		     eta_KaMu, phi_KaMu, charge_KaMu, _isData );

      _muonInfo.pt_KaMu           = pt_KaMu;
      _muonInfo.ptErr_KaMu        = ptErr_KaMu;
      _muonInfo.pt_sys_up_KaMu    = pt_sys_up_KaMu;
      _muonInfo.pt_sys_down_KaMu  = pt_sys_down_KaMu;
      _muonInfo.pt_clos_up_KaMu   = pt_clos_up_KaMu;
      _muonInfo.pt_clos_down_KaMu = pt_clos_down_KaMu;

      // Rochester-calibrated pT
      TLorentzVector mu_vec_Roch;
      mu_vec_Roch.SetPtEtaPhiM( muon.innerTrack()->pt(), muon.innerTrack()->eta(),
      				muon.innerTrack()->phi(), MASS_MUON );
      float q_term_Roch, pt_sys_up_Roch, pt_sys_down_Roch;
      int charge_Roch     = muon.innerTrack()->charge();
      int trk_layers_Roch = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement(); 

      CorrectPtRoch( _Roch_calib, _doSys_Roch, 
		     mu_vec_Roch, q_term_Roch, 
		     pt_sys_up_Roch, pt_sys_down_Roch,
      		     charge_Roch, trk_layers_Roch, _isData );

      _muonInfo.pt_Roch           = mu_vec_Roch.Pt();
      _muonInfo.q_term_Roch       = q_term_Roch;
      _muonInfo.pt_sys_up_Roch    = pt_sys_up_Roch;
      _muonInfo.pt_sys_down_Roch  = pt_sys_down_Roch;

    } // End if (muon.isTrackerMuon())

    // Particle flow kinematics
    if ( muon.isPFMuon() ) {
      reco::Candidate::LorentzVector pfmuon = muon.pfP4();
      _muonInfo.pfPt  = pfmuon.Pt();
      _muonInfo.pfEta = pfmuon.Eta();
      _muonInfo.pfPhi = pfmuon.Phi();
    }

    // Standard quality
    _muonInfo.isGlobal     = muon.isGlobalMuon();
    _muonInfo.isTracker    = muon.isTrackerMuon();
    _muonInfo.isStandAlone = muon.isStandAloneMuon();
    _muonInfo.isPF         = muon.isPFMuon();
    
    // Standard vertexing
    _muonInfo.d0_BS= muon.innerTrack()->dxy( beamSpotHandle->position() );
    _muonInfo.dz_BS= muon.innerTrack()->dz ( beamSpotHandle->position() );
    
    _muonInfo.d0_PV = track.dxy( primaryVertex.position() );
    _muonInfo.dz_PV = track.dz ( primaryVertex.position() );

    // Standard isolation
    _muonInfo.trackIsoSumPt     = muon.isolationR03().sumPt;
    _muonInfo.trackIsoSumPtCorr = muon.isolationR03().sumPt; // no correction with only 1 muon (??? - AWB 08.11.16)
    
    // Correct Iso calculation? - AWB 08.11.16
    double isovar = muon.isolationR03().sumPt;
    isovar += muon.isolationR03().hadEt; // tracker + HCAL
    isovar /= muon.pt(); // relative combine isolation
    _muonInfo.relCombIso = isovar;

    // // Standard trigger variables
    // for (unsigned int iTrigger = 0; iTrigger < _triggerNames.size(); iTrigger++) {
    //   if (iTrigger >= triggerArraySize) {
    // 	std::cout << "Found " << i+1 << "th muon trigger; only " << triggerArraySize << " allowed in array" << std::endl;
    // 	continue;
    //   }
    //   // _muonInfo.hltPt[iTrigger]  = ???; - AWB 12.11.16
    //   // _muonInfo.hltEta[iTrigger] = ???; - AWB 12.11.16
    //   // _muonInfo.hltPhi[iTrigger] = ???; - AWB 12.11.16
    // }

    // // -----------------------------
    // // DON'T SAVE ADVANCED VARIABLES - AWB 12.11.16
    // // -----------------------------
    // if ( _outputLevel != "Fat") 
    //   continue;
    
    // Advanced isolation
    _muonInfo.ecalIso = muon.isolationR03().emEt;
    _muonInfo.hcalIso = muon.isolationR03().hadEt;

    if ( muon.isPFMuon() ) {

      _muonInfo.sumChargedHadronPtR03   = muon.pfIsolationR03().sumChargedHadronPt  ;
      _muonInfo.sumChargedParticlePtR03 = muon.pfIsolationR03().sumChargedParticlePt;
      _muonInfo.sumNeutralHadronEtR03   = muon.pfIsolationR03().sumNeutralHadronEt  ;
      _muonInfo.sumPhotonEtR03          = muon.pfIsolationR03().sumPhotonEt         ;
      _muonInfo.sumPUPtR03              = muon.pfIsolationR03().sumPUPt             ;
      
      _muonInfo.sumChargedHadronPtR04   = muon.pfIsolationR04().sumChargedHadronPt  ;
      _muonInfo.sumChargedParticlePtR04 = muon.pfIsolationR04().sumChargedParticlePt;
      _muonInfo.sumNeutralHadronEtR04   = muon.pfIsolationR04().sumNeutralHadronEt  ;
      _muonInfo.sumPhotonEtR04          = muon.pfIsolationR04().sumPhotonEt         ;
      _muonInfo.sumPUPtR04              = muon.pfIsolationR04().sumPUPt             ;
    }

    _muonInfos.push_back( _muonInfo );
  } // End loop: for (unsigned int i = 0; i < nMuons; i++)

} // End function: void FillMuonInfos
    

pat::MuonCollection SelectMuons( const edm::Handle<pat::MuonCollection>& muons, 
				 const reco::Vertex primaryVertex, const std::string _muon_ID,
				 const double _muon_pT_min, const double _muon_eta_max, const double _muon_trig_dR, 
				 const bool _muon_use_pfIso, const double _muon_iso_dR, const double _muon_iso_max ) {
  
  // Following https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Identification
  
  pat::MuonCollection muonsSelected;
  muonsSelected.clear();
  
  if ( !muons.isValid() ) {
    std::cout << "No valid muon collection" << std::endl;
    return muonsSelected;
  }
  
  for (pat::MuonCollection::const_iterator muon = muons->begin(), muonsEnd = muons->end(); muon !=muonsEnd; ++muon) {

    if ( muon->pt()        < _muon_pT_min  ) continue;
    if ( fabs(muon->eta()) > _muon_eta_max ) continue;
    
    if ( _muon_ID.find("loose")  != std::string::npos && !MuonIsLoose ( (*muon) )                ) continue;
    if ( _muon_ID.find("medium") != std::string::npos && !MuonIsMedium( (*muon) )                ) continue;
    if ( _muon_ID.find("tight")  != std::string::npos && !MuonIsTight ( (*muon), primaryVertex ) ) continue;

    if ( _muon_use_pfIso ) {
      if ( MuonCalcRelIsoPF ( (*muon), _muon_iso_dR ) > _muon_iso_max ) continue;
    } else {
      if ( MuonCalcRelIsoTrk( (*muon), _muon_iso_dR ) > _muon_iso_max ) continue;
    }
    
    muonsSelected.push_back(*muon);
  }
  
  return muonsSelected;
}

bool MuonIsLoose ( const pat::Muon muon ) {
  bool _isLoose = ( muon.isPFMuon() && ( muon.isGlobalMuon() || muon.isTrackerMuon() ) );
  if ( _isLoose  != muon::isLooseMuon(muon)  )
    std::cout << "Manual muon isLoose = " << _isLoose << ", muon::isLooseMuon = " 
	      << muon::isLooseMuon(muon) << ". Using manual version ..." << std::endl;
  return _isLoose;
}

bool MuonIsMedium( const pat::Muon muon ) {
  bool _goodGlob = ( muon.isGlobalMuon()                           && 
		     muon.globalTrack()->normalizedChi2() < 3      && 
		     muon.combinedQuality().chi2LocalPosition < 12 && 
		     muon.combinedQuality().trkKink < 20           ); 
  bool _isMedium = ( MuonIsLoose( muon )                                            && 
		     muon.innerTrack()->validFraction() > 0.8                       && 
		     muon::segmentCompatibility(muon) > (_goodGlob ? 0.303 : 0.451) ); 
  if ( _isMedium  != muon::isMediumMuon(muon)  )
    std::cout << "Manual muon isMedium = " << _isMedium << ", muon::isMediumMuon = " 
	      << muon::isMediumMuon(muon) << ". Using manual version ..." << std::endl;
  return _isMedium;
}

bool MuonIsTight( const pat::Muon muon, const reco::Vertex primaryVertex ) {
  bool _isTight  = ( muon.isGlobalMuon() && muon.isPFMuon()                            && 
		     muon.globalTrack()->normalizedChi2() < 10.                        &&
		     muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0      && 
		     muon.numberOfMatchedStations() > 1                                &&
		     // fabs(muon.muonBestTrack()->dB()) < 0.2                            &&  // How do we implement this? - AWB 14.11.16 
		     fabs(muon.muonBestTrack()->dz( primaryVertex.position() )) < 0.5  &&
		     muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0      &&
		     muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 );
  
  if ( _isTight  != muon::isTightMuon(muon, primaryVertex)  )
    std::cout << "Manual muon isTight = " << _isTight << ", muon::isTightMuon = " 
	      << muon::isTightMuon(muon, primaryVertex) << ". Using manual version ..." << std::endl;
  return _isTight;
}

double MuonCalcRelIsoPF( const pat::Muon muon, const double _muon_iso_dR ) {
  // Following https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
  double iso = -999.;
  if ( _muon_iso_dR == 0.4 ) {
    iso  = muon.pfIsolationR04().sumChargedHadronPt;
    iso += std::max( 0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt );
  } 
  else if ( _muon_iso_dR == 0.3 ) {
    iso  = muon.pfIsolationR03().sumChargedHadronPt;
    iso += std::max( 0., muon.pfIsolationR03().sumNeutralHadronEt + muon.pfIsolationR03().sumPhotonEt - 0.5*muon.pfIsolationR03().sumPUPt );
  }
  else {
    std::cout << "Muon isolation dR = " << _muon_iso_dR << ", not 0.3 or 0.4.  Cannot compute, setting to -999." << std::endl;
    return iso;
  }
  return (iso / muon.pt() );
}

double MuonCalcRelIsoTrk( const pat::Muon muon, const double _muon_iso_dR ) {
  // Following https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
  double iso = -999.;
  // if ( _muon_iso_dR == 0.4 )
  //   iso = muon.isolationR04().sumPt;  // Apparently not an option
  if ( _muon_iso_dR == 0.3 )
    iso = muon.isolationR03().sumPt;
  else {
    std::cout << "Muon isolation dR = " << _muon_iso_dR << ", not 0.3.  Cannot compute, setting to -999." << std::endl;
    return iso;
  }
  return (iso / muon.pt() );
}

double MuonCalcTrigEff( const pat::Muon muon, const int _nPV, const std::string _triggerName ) {

  // Guestimates based on https://indico.cern.ch/event/570616/contributions/2359949/attachments/1365003/2067370/Muon_HLT_Budapest_WS.pdf
  
  double _HLT_eff_vs_nPV = 0.99 - 0.04 * pow( (_nPV/40.), 2 );
  double _HLT_eff_vs_pT  = 1.;  // Depends on muon.pt() near the threshold? - AWB 14.11.16

  double _L1T_eff_vs_eta = -999.;
  if      ( fabs(muon.eta()) < 1.2 ) _L1T_eff_vs_eta = 0.96;
  else if ( fabs(muon.eta()) < 1.6 ) _L1T_eff_vs_eta = 0.94;
  else if ( fabs(muon.eta()) < 2.1 ) _L1T_eff_vs_eta = 0.85;
  else                               _L1T_eff_vs_eta = 0.82;
  
  if ( _triggerName.find("Iso") != std::string::npos ) {
    return _L1T_eff_vs_eta * _HLT_eff_vs_nPV * _HLT_eff_vs_pT; }
  else
    return _L1T_eff_vs_eta * _HLT_eff_vs_pT;
}

bool IsHltMatched( const pat::Muon& muon, const edm::Event& iEvent, const edm::EventSetup& iSetup,
		   const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
		   const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
		   const std::string _desiredTriggerName, const double _muon_trig_dR ) {
  
  // Check if the muon is the one firing the HLT path
  
  using namespace std;
  using namespace edm;
  using namespace pat;
  using namespace reco;
  using namespace trigger;
  
  const boost::regex rgx("_v[0-9]+");
  
  const TriggerNames &triggerNames = iEvent.triggerNames(*_triggerResultsHandle);
  
  const unsigned nTriggers = _triggerResultsHandle->size();
  for (unsigned iTrigger = 0; iTrigger < nTriggers; ++iTrigger) {
    const string triggerName = triggerNames.triggerName(iTrigger);
    string triggerNameStripped = boost::regex_replace(triggerName, rgx, "", boost::match_default | boost::format_sed);
    if ( _desiredTriggerName == triggerNameStripped && _triggerResultsHandle->accept(iTrigger) ) {
      stringstream debugString;
      debugString << "IsHltMatched: ";
      debugString << "'" << _desiredTriggerName << "'\n";
      debugString << "  Trigger " << iTrigger << ": " << triggerName << "(" << triggerNameStripped
		  << ") passed: " << _triggerResultsHandle->accept(iTrigger) << endl;
      
      for ( TriggerObjectStandAloneCollection::const_iterator trigObj = (*_triggerObjsHandle).begin();
	    trigObj != (*_triggerObjsHandle).end(); trigObj++) {
	
	TriggerObjectStandAlone tmpTrigObj(*trigObj); // Get rid of const which messes up unpackPathNames
	tmpTrigObj.unpackPathNames(triggerNames);
	bool isRightObj = tmpTrigObj.hasPathName(triggerName, true, true); // Name, is L3 filter accepted, is last filter
	if (isRightObj) {
	  debugString << "    TriggerObject:  " << tmpTrigObj.collection() << endl;
	  bool isMatched = ( deltaR(tmpTrigObj, muon) < _muon_trig_dR );
	  if (isMatched) {
	    debugString << "      is Matched*****" << endl;
	    LogVerbatim("UFHLTTests") << debugString.str();
	    return true;
	  }
	} // if isRightObj
      } // trigObj lookp
      LogVerbatim("UFHLTTests") << debugString.str();
    } // if desiredTrigger
  } // iTrigger loop
  
  return false;
}

