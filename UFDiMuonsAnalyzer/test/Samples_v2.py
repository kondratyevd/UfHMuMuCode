# root://xrootd-cms.infn.it/

# This is my current collection of datasets for H->MuMu.
# It currently does not contain all of the necessary datasets but has enough to play around with
# and get started.

class sample:
    def __init__(self, name="", dir="", files=[], numevents=0, globaltag="", jsonfiles=[], isData=False):
        self.name = name
        self.dir = dir     # DAS directory
        self.numevents = numevents
        self.files = files
        self.globaltag = globaltag
        self.jsonfiles = jsonfiles
        self.isData = isData


sample_array = [];

# =======================================================================================================
# ------------------------------- DATA ------------------------------------------------------------------
# =======================================================================================================

# old json files
jsonlist_v0 = ['sample_file_lists/data/json/Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON_MuonPhys.txt',
            'sample_file_lists/data/json/Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON.txt'] 

jsonlist = ['sample_file_lists/data/json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_MuonPhys_v2.txt',
            'sample_file_lists/data/json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'] 
# CERN
# CMSSW 7_4_5
doubleMuon_RunA_AOD = sample(name="doubleMuon_RunA_AOD", 
                             dir="/DoubleMuon/Run2015A-PromptReco-v1/AOD", 
                             files = open('sample_file_lists/data/doubleMuon_RunA_AOD.files').read().splitlines(), 
                             numevents=348, 
                             isData = True)

# CERN, FL
# CMSSW 7_4_7
doubleMuon_RunB_AOD = sample(name="doubleMuon_RunB_AOD", 
                             dir="/DoubleMuon/Run2015B-PromptReco-v1/AOD", 
                             files = open('sample_file_lists/data/doubleMuon_RunB_AOD.files').read().splitlines(),
                             numevents=1631653,
                             globaltag = '74X_dataRun2_Prompt_v1',
                             jsonfiles = jsonlist[:],
                             isData = True)

# CERN, FL
# CMSSW 7_4_7
doubleMuon_RunB_MINIAOD = sample(name="doubleMuon_RunB_MINIAOD", 
                                 dir="/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD", 
                                 files = open('sample_file_lists/data/doubleMuon_RunB_MINIAOD.files').read().splitlines(),
                                 numevents=1631653,
                                 globaltag = '74X_dataRun2_Prompt_v1',
                                 jsonfiles = jsonlist[:],
                                 isData = True)


sample_array.append(doubleMuon_RunA_AOD)
sample_array.append(doubleMuon_RunB_AOD)
sample_array.append(doubleMuon_RunB_MINIAOD)


# =======================================================================================================
# ------------------------------- SIGNAL ----------------------------------------------------------------
# =======================================================================================================

# ///////////////////////////////////////////////////////////////////////////////////////////////////////
#---- Gluon Gluon Fusion --------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////////////////

ggToHToMuMu_PU40bx50 = sample( name="ggToHToMuMu_PU40bx50", 
                               dir="/GluGluToHToMuMu_M-125_13TeV-powheg-pythia6/Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/MINIAODSIM", 
                               files = open('sample_file_lists/signal/ggToHToMuMu_PU40bx50.files').read().splitlines(),
                               numevents=99994,
                               globaltag = 'PLS170_V6AN2')


ggToHToMuMu_PU20bx25 = sample( name="ggToHToMuMu_PU20bx25", 
                               dir="/GluGluToHToMuMu_M-125_13TeV-powheg-pythia6/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM",
                               files = open('sample_file_lists/signal/ggToHToMuMu_PU20bx25.files').read().splitlines(),
                               numevents=99994,
                               globaltag = 'POSTLS170_V5')

sample_array.append(ggToHToMuMu_PU40bx50)
sample_array.append(ggToHToMuMu_PU20bx25)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////
#---- Vector Boson Fusion --------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////////////////

vbf_HToMuMu_PU40bx50 = sample(name="vbf_HToMuMu_PU40bx50", 
                              dir="/VBF_HToMuMu_M-125_13TeV-powheg-pythia6/Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/MINIAODSIM",
                              files = open('sample_file_lists/signal/vbf_HToMuMu_PU40bx50.files').read().splitlines(),
                              numevents=95968,
                              globaltag = 'PLS170_V6AN2')

vbf_HToMuMu_PU20bx25 = sample(name="vbf_HToMuMu_PU20bx25", 
                              dir="/VBF_HToMuMu_M-125_13TeV-powheg-pythia6/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM",
                              files = open('sample_file_lists/signal/vbf_HToMuMu_PU20bx25.files').read().splitlines(),
                              numevents=95968,
                              globaltag = 'POSTLS170_V5')

sample_array.append(vbf_HToMuMu_PU40bx50)
sample_array.append(vbf_HToMuMu_PU20bx25)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////
#---- W/Z to Higgs --------------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////////////////

wh_zh_HToMuMu_PU40bx50 = sample(name="wh_zh_HToMuMu_PU40bx50", 
                                dir="/WH_ZH_HToMuMu_M-125_13TeV_pythia6/Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/MINIAODSIM",
                                files = open('sample_file_lists/signal/vbf_HToMuMu_PU20bx25.files').read().splitlines(),
                                numevents=100000,
                                globaltag = 'PLS170_V6AN2')

wh_zh_HToMuMu_PU20bx25 = sample(name="wh_zh_HToMuMu_PU20bx25", 
                                dir="/WH_ZH_HToMuMu_M-125_13TeV_pythia6/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM",
                                files = open('sample_file_lists/signal/vbf_HToMuMu_PU20bx25.files').read().splitlines(),
                                numevents=100000,
                                globaltag = 'POSTLS170_V5')

sample_array.append(wh_zh_HToMuMu_PU40bx50)
sample_array.append(wh_zh_HToMuMu_PU20bx25)

# =======================================================================================================
# ------------------------------- BACKGROUND ------------------------------------------------------------
# =======================================================================================================

# ///////////////////////////////////////////////////////////////////////////////////////////////////////
# ---- DRELL YANN ---------------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////////////////
dy_jetsToLL_PU40bx50 = sample(name="dy_jetsToLL_PU40bx50", 
                              dir="/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/MINIAODSIM",
                              files = open('sample_file_lists/bg/dy_jetsToLL_PU40bx50.files').read().splitlines(),
                              numevents=2651053,
                              globaltag = 'PLS170_V6AN2')

dy_jetsToLL_PU20bx25 = sample(name="dy_jetsToLL_PU20bx25", 
                              dir="/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM",
                              files = open('sample_file_lists/bg/dy_jetsToLL_PU20bx25.files').read().splitlines(),
                              numevents=2725222,
                              globaltag = 'POSTLS170_V5');

sample_array.append(dy_jetsToLL_PU40bx50)
sample_array.append(dy_jetsToLL_PU20bx25)


# ///////////////////////////////////////////////////////////////////////////////////////////////////////
#---- T TBAR --------------------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////////////////

#----------------------------------------------------------------------------------------------------
#/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU40bx25_POSTLS170_V7-v2/MINIAODSIM
#----------------------------------------------------------------------------------------------------


ttJets_PU40bx25 = sample(name="ttJets_PU40bx25", 
                         dir="/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU40bx25_POSTLS170_V7-v2/MINIAODSIM",
                         files = open('sample_file_lists/bg/ttJets_PU40bx25.files').read().splitlines(),
                         numevents=25484264,
                         globaltag = 'POSTLS170_V7');

#----------------------------------------------------------------------------------------------------
#/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v2/MINIAODSIM
#----------------------------------------------------------------------------------------------------


ttJets_PU20bx25 = sample(name="ttJets_PU20bx25", 
                         dir="/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v2/MINIAODSIM",
                         files = open('sample_file_lists/bg/ttJets_PU20bx25.files').read().splitlines(),
                         numevents=25474122,
                         globaltag = 'POSTLS170_V5')

sample_array.append(ttJets_PU40bx25)
sample_array.append(ttJets_PU20bx25)
