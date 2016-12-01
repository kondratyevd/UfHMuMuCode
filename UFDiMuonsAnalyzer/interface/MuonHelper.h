
#ifndef MUON_HELPER
#define MUON_HELPER

#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/MuonInfo.h"

void FillMuonInfos( MuonInfos& _muonInfos, 
		    const pat::MuonCollection muonsSelected,
		    const reco::Vertex primaryVertex, const int _nPV,
		    const edm::Handle<reco::BeamSpot>& beamSpotHandle,
		    const edm::Event& iEvent, const edm::EventSetup& iSetup,
		    const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
		    const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
		    const std::vector<std::string> _triggerNames, const double _muon_trig_dR,
		    const bool _muon_use_pfIso, const double _muon_iso_dR );

pat::MuonCollection SelectMuons( const edm::Handle<pat::MuonCollection>& muons,
				 const reco::Vertex primaryVertex, const std::string _muon_ID,
				 const double _muon_pT_min, const double _muon_eta_max, const double _muon_trig_dR,
				 const bool _muon_use_pfIso, const double _muon_iso_dR, const double _muon_iso_max );

bool MuonIsLoose ( const pat::Muon muon );
bool MuonIsMedium( const pat::Muon muon );
bool MuonIsTight ( const pat::Muon muon, const reco::Vertex primaryVertex );

double MuonCalcRelIsoPF ( const pat::Muon muon, const double _muon_iso_dR );
double MuonCalcRelIsoTrk( const pat::Muon muon, const double _muon_iso_dR );
double MuonCalcTrigEff  ( const pat::Muon muon, const int _nPV, const std::string _triggerName );

bool IsHltMatched( const pat::Muon& mu, const edm::Event& iEvent, const edm::EventSetup& iSetup,
		   const edm::Handle<pat::TriggerObjectStandAloneCollection>& _triggerObjsHandle,
		   const edm::Handle<edm::TriggerResults>& _triggerResultsHandle,
		   const std::string _desiredTriggerName, const double _muon_trig_dR );

#endif  // #ifndef MUON_HELPER                                                                                                                  
