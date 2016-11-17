
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _MuonInfo::init() {
  nMuons     = 0;

  for (unsigned int i = 0; i < arraySize; i++) {
    isTracker[i]    = -999;
    isStandAlone[i] = -999;
    isGlobal[i]     = -999;
    
    isTightID[i]    = -999;
    isMediumID[i]   = -999;
    isLooseID[i]    = -999;
    
    charge[i] = -999;
    pt[i]     = -999;
    eta[i]    = -999;
    phi[i]    = -999;
    
    d0_BS[i] = -999;
    dz_BS[i] = -999;
    
    d0_PV[i] = -999;
    dz_PV[i] = -999;

    relIso[i] = -999;
    
    trackIsoSumPt[i]     = -999;
    trackIsoSumPtCorr[i] = -999;
    hcalIso[i]           = -999;
    ecalIso[i]           = -999;
    relCombIso[i]        = -999;
    
    isPF[i] = -999;
    
    pfPt[i]  = -999;
    pfEta[i] = -999;
    pfPhi[i] = -999;
    
    sumChargedHadronPtR03[i]   = -999;
    sumChargedParticlePtR03[i] = -999;
    sumNeutralHadronEtR03[i]   = -999;
    sumPhotonEtR03[i]          = -999;
    sumPUPtR03[i]              = -999;
    
    sumChargedHadronPtR04[i]   = -999;
    sumChargedParticlePtR04[i] = -999;
    sumNeutralHadronEtR04[i]   = -999;
    sumPhotonEtR04[i]          = -999;
    sumPUPtR04[i]              = -999;
    
    for (unsigned int iTrigger = 0; iTrigger < triggerArraySize; iTrigger++) {
      isHltMatched[i][iTrigger] = -999;
      hltPt[i][iTrigger]        = -999;
      hltEta[i][iTrigger]       = -999;
      hltPhi[i][iTrigger]       = -999;
    }
  } // End loop: for (unsigned int i = 0; i < arraySize; i++)

} // End void _MuonInfo::init()

void _MuonInfo::fill( const pat::MuonCollection muonsSelected,
		      const reco::Vertex primaryVertex, const int _nPV,
		      const edm::Handle<reco::BeamSpot>& beamSpotHandle,  
		      const edm::Event& iEvent, const edm::EventSetup& iSetup,
		      const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
		      const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
		      const std::vector<std::string> _triggerNames, const double _muon_trig_dR, 
		      const bool _muon_use_pfIso, const double _muon_iso_dR ) {
  
  nMuons = muonsSelected.size();

  for (int i = 0; i < nMuons; i++) {

    if ( i >= int(arraySize) ) {
      std::cout << "Found " << i+1 << "th muon; only " << arraySize << " allowed in array" << std::endl;
      return;
    }

    pat::Muon muon = muonsSelected.at(i);

    reco::Track track;
    // Do we want to use the inner tracker track by default for the Kalman corrections? - AWB 12.11.16
    if      ( muon.isGlobalMuon()  ) track = *(muon.globalTrack() );
    else if ( muon.isTrackerMuon() ) track = *(muon.innerTrack()  );
    else {
      std::cout << "ERROR: The muon is NOT global NOR tracker ?!?\n";
      return;
    }

    // Basic kinematics
    charge[i] = muon.charge();
    pt[i]     = muon.pt();
    ptErr[i]  = track.ptError();
    eta[i]    = muon.eta();
    phi[i]    = muon.phi();

    // Basic quality
    isTightID[i]  =  IsTight  ( muon, primaryVertex );
    isMediumID[i] =  IsMedium ( muon );
    isLooseID[i]  =  IsLoose  ( muon );

    // Basic isolation
    if ( _muon_use_pfIso )
      relIso[i] = CalcRelIsoPF ( muon, _muon_iso_dR );
    else
      relIso[i] = CalcRelIsoTrk( muon, _muon_iso_dR );

    // Basic vertexing? - AWB 12.11.16

    // Basic trigger (which triggers are in the list? - AWB 12.11.16)
    for (unsigned int iTrigger = 0; iTrigger < _triggerNames.size(); iTrigger++) {
      if (iTrigger >= triggerArraySize) {
    	std::cout << "Found " << i+1 << "th muon trigger; only " << triggerArraySize << " allowed in array" << std::endl;
    	continue;
      }
      isHltMatched[i][iTrigger] = IsHltMatched( muon, iEvent, iSetup, _triggerObjsHandle, _triggerResultsHandle,  
						_triggerNames[iTrigger], _muon_trig_dR );
      hltEff[i][iTrigger] = CalcTrigEff( muon, _nPV, _triggerNames[iTrigger] );
    }


    // // -------------------------
    // // ONLY SAVE BASIC VARIABLES - AWB 12.11.16
    // // -------------------------
    // if ( _outputLevel == "Slim") 
    //   continue;

    // Inner track kinematics
    if (muon.isTrackerMuon()) {
      trkPt[i]    = muon.innerTrack()->pt();
      trkPtErr[i] = muon.innerTrack()->ptError();
      trketa[i]   = muon.innerTrack()->eta();
      trkPhi[i]   = muon.innerTrack()->phi();
    }

    // Particle flow kinematics
    if ( muon.isPFMuon() ) {
      reco::Candidate::LorentzVector pfmuon = muon.pfP4();
      pfPt[i]  = pfmuon.Pt();
      pfEta[i] = pfmuon.Eta();
      pfPhi[i] = pfmuon.Phi();
    }

    // Standard quality
    isGlobal[i]     = muon.isGlobalMuon();
    isTracker[i]    = muon.isTrackerMuon();
    isStandAlone[i] = muon.isStandAloneMuon();
    isPF[i]         = muon.isPFMuon();
    
    // Standard vertexing
    d0_BS[i]= muon.innerTrack()->dxy( beamSpotHandle->position() );
    dz_BS[i]= muon.innerTrack()->dz ( beamSpotHandle->position() );
    
    d0_PV[i] = track.dxy( primaryVertex.position() );
    dz_PV[i] = track.dz ( primaryVertex.position() );

    // Standard isolation
    trackIsoSumPt[i]    = muon.isolationR03().sumPt;
    trackIsoSumPtCorr[i]= muon.isolationR03().sumPt; // no correction with only 1 muon (??? - AWB 08.11.16)
    
    // Correct Iso calculation? - AWB 08.11.16
    double isovar = muon.isolationR03().sumPt;
    isovar += muon.isolationR03().hadEt; // tracker + HCAL
    isovar /= muon.pt(); // relative combine isolation
    relCombIso[i] = isovar;

    // Standard trigger variables
    for (unsigned int iTrigger = 0; iTrigger < _triggerNames.size(); iTrigger++) {
      if (iTrigger >= triggerArraySize) {
    	std::cout << "Found " << i+1 << "th muon trigger; only " << triggerArraySize << " allowed in array" << std::endl;
    	continue;
      }
      // hltPt[i][iTrigger]  = ???; - AWB 12.11.16
      // hltEta[i][iTrigger] = ???; - AWB 12.11.16
      // hltPhi[i][iTrigger] = ???; - AWB 12.11.16
    }

    // // -----------------------------
    // // DON'T SAVE ADVANCED VARIABLES - AWB 12.11.16
    // // -----------------------------
    // if ( _outputLevel != "Fat") 
    //   continue;
    
    // Advanced isolation
    ecalIso[i]          = muon.isolationR03().emEt;
    hcalIso[i]          = muon.isolationR03().hadEt;

    if ( muon.isPFMuon() ) {

      sumChargedHadronPtR03[i]   = muon.pfIsolationR03().sumChargedHadronPt  ;
      sumChargedParticlePtR03[i] = muon.pfIsolationR03().sumChargedParticlePt;
      sumNeutralHadronEtR03[i]   = muon.pfIsolationR03().sumNeutralHadronEt  ;
      sumPhotonEtR03[i]          = muon.pfIsolationR03().sumPhotonEt         ;
      sumPUPtR03[i]              = muon.pfIsolationR03().sumPUPt             ;
      
      sumChargedHadronPtR04[i]   = muon.pfIsolationR04().sumChargedHadronPt  ;
      sumChargedParticlePtR04[i] = muon.pfIsolationR04().sumChargedParticlePt;
      sumNeutralHadronEtR04[i]   = muon.pfIsolationR04().sumNeutralHadronEt  ;
      sumPhotonEtR04[i]          = muon.pfIsolationR04().sumPhotonEt         ;
      sumPUPtR04[i]              = muon.pfIsolationR04().sumPUPt             ;
    }

  } // End loop: for (unsigned int i = 0; i < nMuons; i++)
    
} // End void _MuonInfo::fill()
    

pat::MuonCollection _MuonInfo::select( const edm::Handle<pat::MuonCollection>& muons, 
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
    
    if ( _muon_ID.find("loose")  != std::string::npos && !IsLoose ( (*muon) )                ) continue;
    if ( _muon_ID.find("medium") != std::string::npos && !IsMedium( (*muon) )                ) continue;
    if ( _muon_ID.find("tight")  != std::string::npos && !IsTight ( (*muon), primaryVertex ) ) continue;

    if ( _muon_use_pfIso ) {
      if ( CalcRelIsoPF ( (*muon), _muon_iso_dR ) > _muon_iso_max ) continue;
    } else {
      if ( CalcRelIsoTrk( (*muon), _muon_iso_dR ) > _muon_iso_max ) continue;
    }
    
    muonsSelected.push_back(*muon);
  }
  
  return muonsSelected;
}

bool _MuonInfo::IsLoose ( const pat::Muon muon ) {
  bool _isLoose = ( muon.isPFMuon() && ( muon.isGlobalMuon() || muon.isTrackerMuon() ) );
  if ( _isLoose  != muon::isLooseMuon(muon)  )
    std::cout << "Manual muon isLoose = " << _isLoose << ", muon::isLooseMuon = " 
	      << muon::isLooseMuon(muon) << ". Using manual version ..." << std::endl;
  return _isLoose;
}

bool _MuonInfo::IsMedium( const pat::Muon muon ) {
  bool _goodGlob = ( muon.isGlobalMuon()                           && 
		     muon.globalTrack()->normalizedChi2() < 3      && 
		     muon.combinedQuality().chi2LocalPosition < 12 && 
		     muon.combinedQuality().trkKink < 20           ); 
  bool _isMedium = ( IsLoose( muon )                                                && 
		     muon.innerTrack()->validFraction() > 0.8                       && 
		     muon::segmentCompatibility(muon) > (_goodGlob ? 0.303 : 0.451) ); 
  if ( _isMedium  != muon::isMediumMuon(muon)  )
    std::cout << "Manual muon isMedium = " << _isMedium << ", muon::isMediumMuon = " 
	      << muon::isMediumMuon(muon) << ". Using manual version ..." << std::endl;
  return _isMedium;
}

bool _MuonInfo::IsTight( const pat::Muon muon, const reco::Vertex primaryVertex ) {
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

double _MuonInfo::CalcRelIsoPF( const pat::Muon muon, const double _muon_iso_dR ) {
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

double _MuonInfo::CalcRelIsoTrk( const pat::Muon muon, const double _muon_iso_dR ) {
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

double _MuonInfo::CalcTrigEff( const pat::Muon muon, const int _nPV, const std::string _triggerName ) {

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

bool _MuonInfo::IsHltMatched( const pat::Muon& muon, const edm::Event& iEvent, const edm::EventSetup& iSetup,
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

