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
jsonlist = ['sample_file_lists/data/json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_MuonPhys_v4.txt',
            'sample_file_lists/data/json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'] 

# CERN, FL
# CMSSW 7_4_7
singleMuon_RunBPrompt_MINIAOD = sample(name="singleMuon_RunBPrompt_AOD", 
                             dir="/SingleMu/Run2015B-PromptReco-v1/MINIAOD", 
                             files = open('sample_file_lists/data/singleMuon_RunBPrompt_MINIAOD.files').read().splitlines(),
                             numevents=504556,
                             globaltag = '74X_dataRun2_Prompt_v0',
                             jsonfiles = jsonlist[:],
                             isData = True)

doubleMuon_RunBPrompt_MINIAOD = sample(name="doubleMuon_RunBPrompt_MINIAOD", 
                                 dir="/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD", 
                                 files = open('sample_file_lists/data/doubleMuon_RunBPrompt_MINIAOD.files').read().splitlines(),
                                 numevents=1631653,
                                 globaltag = '74X_dataRun2_Prompt_v1',
                                 jsonfiles = jsonlist[:],
                                 isData = True)

sample_array.append(singleMuon_RunBPrompt_MINIAOD)
sample_array.append(doubleMuon_RunBPrompt_MINIAOD)


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

dy_ZToMuMu_asympt50 = sample(name="dy_ZToMuMu_asympt50", 
                              dir="/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM",
                              files = open('sample_file_lists/bg/dy_ZToMuMu_asympt50.files').read().splitlines(),
                              numevents=2895638,
                              globaltag = 'MCRUN2_74_V9A')

dy_jetsToLL_asympt50 = sample(name="dy_jetsToLL_asympt50", 
                              dir="/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM",
                              files = open('sample_file_lists/bg/dy_jetsToLL_asympt50.files').read().splitlines(),
                              numevents=19925500,
                              globaltag = 'MCRUN2_74_V9A')

dy_jetsToLL_asympt25 = sample(name="dy_jetsToLL_asympt25", 
                              dir="/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM",
                              files = open('sample_file_lists/bg/dy_jetsToLL_asympt25.files').read().splitlines(),
                              numevents=28825132,
                              globaltag = 'MCRUN2_74_V9');

sample_array.append(dy_ZToMuMu_asympt50)
sample_array.append(dy_jetsToLL_asympt50)
sample_array.append(dy_jetsToLL_asympt25)


# ///////////////////////////////////////////////////////////////////////////////////////////////////////
#---- T TBAR --------------------------------------------------------------------------------------------
# ///////////////////////////////////////////////////////////////////////////////////////////////////////

ttJets_asympt50 = sample(name="ttJets_asympt50", 
                         dir="TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM",
                         files = open('sample_file_lists/bg/ttJets_asympt50.files').read().splitlines(),
                         numevents=4994250,
                         globaltag = 'MCRUN2_74_V9A');


ttJets_asympt25 = sample(name="ttJets_asympt25", 
                         dir="TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM",
                         files = open('sample_file_lists/bg/ttJets_asympt25.files').read().splitlines(),
                         numevents=42730273,
                         globaltag = 'MCRUN2_74_V9')

sample_array.append(ttJets_asympt50)
sample_array.append(ttJets_asympt25)
