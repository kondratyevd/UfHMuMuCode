
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/DataFormats.h"

void _ElectronInfo::init() {
  nElectrons = 0;
  for (unsigned int i = 0; i < arraySize; i++) {
    isPF[i]               = -999;
    isTightID[i]          = -999;
    isMediumID[i]         = -999;
    isLooseID[i]          = -999;
    isVetoID[i]           = -999;
    passConversionVeto[i] = -999;
    
    charge[i] = -999;
    pt[i]     = -999;
    eta[i]    = -999;
    phi[i]    = -999;
    
    d0_PV[i] = -999;
    dz_PV[i] = -999;
    
    missingInnerHits[i] = -999;
    
    sumChargedHadronPtR03[i] = -999;
    sumNeutralHadronEtR03[i] = -999;
    sumPhotonEtR03[i]        = -999;
    sumPUPtR03[i]            = -999;
  }
}

void _ElectronInfo::fill( const pat::ElectronCollection electronsSelected,
			  const reco::Vertex primaryVertex, const edm::Event& iEvent,
			  const edm::Handle< edm::ValueMap<bool> >& ele_id_veto, const edm::Handle< edm::ValueMap<bool> >& ele_id_loose,
			  const edm::Handle< edm::ValueMap<bool> >& ele_id_medium, const edm::Handle< edm::ValueMap<bool> >& ele_id_tight ) {

  nElectrons = electronsSelected.size();
  
  for (int i = 0; i < nElectrons; i++) {
    
    if ( i >= int(arraySize) ) {
      std::cout << "Found " << i+1 << "th electron; only " << arraySize << " allowed in array" << std::endl;
      return;
    }

    pat::Electron ele  = electronsSelected.at(i);

    // Basic kinematics
    charge[i] = ele.charge();
    pt[i]     = ele.pt();
    eta[i]    = ele.eta();
    phi[i]    = ele.phi();

    // Basic quality
    isPF[i]       = ele.isPF();
    // // Not sure how to access IDs ... code below does not work.  Appears to require edm::View. - AWB 15.11.16
    // const auto ele_ptr = electronsSelected.ptrAt(i);
    // isTightID[i]  = (*ele_id_tight )[ele_ptr] && PassKinematics(ele, primaryVertex);
    // isMediumID[i] = (*ele_id_medium)[ele_ptr] && PassKinematics(ele, primaryVertex);
    // isLooseID[i]  = (*ele_id_loose )[ele_ptr] && PassKinematics(ele, primaryVertex);
    // isVetoID[i]   = (*ele_id_veto  )[ele_ptr] && PassKinematics(ele, primaryVertex);

    // Basic isolation
    relIso[i] = CalcRelIsoPF_DeltaBeta( ele );

    // Basic vertexing?

    // Efficiencies and scale factors?
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale

    // // -------------------------
    // // ONLY SAVE BASIC VARIABLES - AWB 12.11.16
    // // -------------------------
    // if ( _outputLevel == "Slim")
    //   continue;

    // Standard quality
    passConversionVeto[i] = ele.passConversionVeto();
    missingInnerHits[i]   = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    // Standard vertexing
    d0_PV[i] = ele.gsfTrack()->dxy( primaryVertex.position() );
    dz_PV[i] = ele.gsfTrack()->dz ( primaryVertex.position() );

    // Standard isolation
    sumChargedHadronPtR03[i] = ele.pfIsolationVariables().sumChargedHadronPt;
    sumNeutralHadronEtR03[i] = ele.pfIsolationVariables().sumNeutralHadronEt;
    sumPhotonEtR03[i]        = ele.pfIsolationVariables().sumPhotonEt;
    sumPUPtR03[i]            = ele.pfIsolationVariables().sumPUPt;
    
  } // End loop: for (int i = 0; i < nElectrons; i++)
  
}


pat::ElectronCollection _ElectronInfo::select( const edm::Handle<edm::View<pat::Electron>>& electrons, const reco::Vertex primaryVertex,
					       const edm::Handle< edm::ValueMap<bool> >& ele_id_veto, const edm::Handle< edm::ValueMap<bool> >& ele_id_loose,
					       const edm::Handle< edm::ValueMap<bool> >& ele_id_medium, const edm::Handle< edm::ValueMap<bool> >& ele_id_tight,
					       const std::string _electron_ID, const double _electron_pT_min, const double _electron_eta_max ) {
  
  // Main Egamma POG page: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPOG
  // Following https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
  //       and https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
  // Modeled after https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc
  
  pat::ElectronCollection electronsSelected;
  electronsSelected.clear();

  if ( !electrons.isValid() ) {
    std::cout << "No valid electron collection" << std::endl;
    return electronsSelected;
  }
  if ( !ele_id_veto.isValid() || !ele_id_loose.isValid() || !ele_id_medium.isValid() || !ele_id_tight.isValid() ) {
    std::cout << "No valid electron ID" << std::endl;
    return electronsSelected;
  }

  if ( _electron_ID.find("veto")   == std::string::npos && _electron_ID.find("loose") == std::string::npos && 
       _electron_ID.find("medium") == std::string::npos && _electron_ID.find("tight") == std::string::npos )
    std::cout << "Electron ID is neither tight, medium, loose, nor tight: " << _electron_ID
              << "\nWill not be used, no electron ID cuts applied" << std::endl;

  for (size_t i = 0; i < electrons->size(); ++i) {
    
    const auto ele = electrons->ptrAt(i);

    if ( ele->pt()          < _electron_pT_min  ) continue;
    if ( fabs( ele->eta() ) > _electron_eta_max ) continue;

    bool _isVeto   = (*ele_id_veto  )[ele] && PassKinematics(*ele, primaryVertex);
    bool _isLoose  = (*ele_id_loose )[ele] && PassKinematics(*ele, primaryVertex);
    bool _isMedium = (*ele_id_medium)[ele] && PassKinematics(*ele, primaryVertex);
    bool _isTight  = (*ele_id_tight )[ele] && PassKinematics(*ele, primaryVertex);

    if (_electron_ID.find("veto")   != std::string::npos && !_isVeto)   continue;
    if (_electron_ID.find("loose")  != std::string::npos && !_isLoose)  continue;
    if (_electron_ID.find("medium") != std::string::npos && !_isMedium) continue;
    if (_electron_ID.find("tight")  != std::string::npos && !_isTight)  continue;

    electronsSelected.push_back(*ele);

  }
  
  return electronsSelected;
}

bool _ElectronInfo::PassKinematics( const pat::Electron ele, const reco::Vertex primaryVertex ) {

  double SC_eta = -999.;
  if ( ele.superCluster().isAvailable() )
    SC_eta = fabs( ele.superCluster()->position().eta() );
  bool isBarrel = ( SC_eta != -999 && SC_eta < 1.479 );
  bool inCrack  = ( SC_eta != -999 && SC_eta > 1.4442 && SC_eta < 1.5660 );
  
  double dXY = -999.;
  double dZ  = -999.;
  if ( ele.gsfTrack().isAvailable() ) {
    dXY = fabs( ele.gsfTrack()->dxy( primaryVertex.position() ) );
    dZ  = fabs( ele.gsfTrack()->dz ( primaryVertex.position() ) );
  }
  
  // From https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria
  // Different than https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/src/MiniAODHelper.cc#L890
  bool passDXY = dXY != -999 && ( isBarrel ? (dXY < 0.05) : (dXY < 0.10 ) );
  bool passDZ  = dZ  != -999 && ( isBarrel ? (dZ  < 0.10) : (dZ  < 0.20 ) );

  // What about passMVAId53x? no_exp_inner_trkr_hits? myTrigPresel? - AWB 15.11.16
  return (passDXY && passDZ && ele.passConversionVeto() && !inCrack);
  // Does ele->passConversionVeto() return the same thing as ConversionTools::hasMatchedConversion? - AWB 15.11.16
  // https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc#L358
}  
  
double _ElectronInfo::CalcRelIsoPF_DeltaBeta( const pat::Electron ele ) {
  // Following https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
  // Using Delta Beta corrections - simpler, only slightly worse performance
  
  double iso_charged    = ele.pfIsolationVariables().sumChargedHadronPt;  // Is this correct? should be chargedHadronIso?   - AWB 14.11.16
  double iso_neutral    = ele.pfIsolationVariables().sumNeutralHadronEt;    // Is this correct? should be sumNeutralHadronPt? - AWB 14.11.16
  double iso_photon     = ele.pfIsolationVariables().sumPhotonEt;
  double iso_PU_charged = ele.pfIsolationVariables().sumPUPt;
  
  double rel_iso = iso_charged;
  rel_iso += std::max( 0.0, iso_neutral + iso_photon - (0.5 * iso_PU_charged) );
  rel_iso /= ele.pt();
  
  return rel_iso;
}
