
#ifndef JET_HELPER
#define JET_HELPER

#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/CommonIncludes.h"
#include "UfHMuMuCode/UFDiMuonsAnalyzer/interface/JetInfo.h"

void FillJetInfos( JetInfos& _jetInfos, int& _nJetsFwd,
		   const pat::JetCollection jetsSelected, 
		   const std::vector<std::string> _btagNames );

pat::JetCollection SelectJets( const edm::Handle<pat::JetCollection>& jets,
			       const JetCorrectorParameters& JetCorPar, const std::string _JES_syst,
			       const std::string _jet_ID, const double _jet_pT_min, const double _jet_eta_max );

#endif  // #ifndef JET_HELPER
