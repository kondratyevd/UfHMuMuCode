#!/usr/bin/env python

import os

class mcsample(object):
    def __init__(self, name, dataset, ntuple, caf_folder_extd, global_tag, nevents, cross_section, energy):
        self.name = name
        self.dataset = dataset
        self.ntuple = ntuple
        self.caf_folder_extd = caf_folder_extd
        self.global_tag = global_tag
        self.nevents = nevents
        self.cross_section = cross_section
        self.energy = energy

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV for xsecs (all below in pb)
# or from PREP
MCSamples = [
    mcsample('dyjetstoll',
           '/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYJetsToLL', 'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 30459503, 3503.71, 8),
    mcsample('ttjets',
           '/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
           'TTJets', 'TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3',
           'START53_V15', 6923750, 225.197, 8),
    mcsample('ww',
           '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'WW', 'WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 9900431, 54.838, 8),
    mcsample('wz',
           '/WZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'WZ', 'WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 10000283, 33.21, 8),
    mcsample('zz',
           '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'ZZ', 'ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 9799908, 17.654, 8),
    mcsample('dytotautau',
           '/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYToTauTau', 'DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 3295238, 5745.25/3., 8),
#    mcsample('dytomumu',
#           '/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
#           'DYToMuMu', 'DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1',
#           'START53_V15', 3293740, 5745.25/3., 8),
    mcsample('dytomumu',
           '/DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYToMuMu', 'DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 3293740, 5745.25/3., 8),
    mcsample('wjetstolnu',
           '/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
           'WJetsToLNu', 'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2',
           'START53_V15', 57330800, 30400.0, 8),
    mcsample('qcd',
           '/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
           'QCD_Pt_20_MuEnrichedPt_15', 'QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3',
           'START53_V15',20284602, 3.64E8*3.7E-4, 8),

    # special MC datasets
    mcsample('dy1jetstoll',
           '/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DY1JetsToLL', 'DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 24045248, 9999, 8),
    mcsample('dy2jetstoll',
           '/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
           'DY2JetsToLL', 'DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1',
           'START53_V15', 21852156, 9999, 8),
    mcsample('dy3jetstoll',
           '/DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DY3JetsToLL', 'DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 11015445, 9999, 8),
    mcsample('dy4jetstoll',
           '/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DY4JetsToLL', 'DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 6402827, 9999, 8),


    mcsample('dyjetstoll_ptz-100',
           '/DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
           'DYJetsToLL_PtZ-100', 'DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2',
           'START53_V15', 4930773, 93.8, 8),
    mcsample('dyjetstoll_ptz-50to70',
           '/DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYJetsToLL_PtZ-50to70', 'DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 1313395, 52.31, 8),
    mcsample('dyjetstoll_ptz-70to100',
           '/DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
           'DYJetsToLL_PtZ-70to100', 'DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2',
           'START53_V15', 2662137, 34.1, 8),

    mcsample('dy2jetstoll',
           '/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DY2JetsToLL', 'DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 2352304, 181.0, 8),
    mcsample('dy3jetstoll',
           '/DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DY3JetsToLL', 'DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 11015445, 51.1, 8),
    mcsample('dy4jetstoll',
           '/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DY4JetsToLL', 'DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V15', 6402827, 23.04, 8),

     # 2011 Samples
    mcsample('dyjetstoll2011',
           '/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'DYJetsToLL', 'DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 36264432, 3048., 7),
    mcsample('ttjets2011',
           '/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'TTJets', 'TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 59444088, 157.5, 7),
    mcsample('ww2011',
           '/WW_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'WW', 'WW_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 4225916,  43.0, 7),
    mcsample('wz2011',
           '/WZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'WZ', 'WZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 4265243, 18.2, 7),
    mcsample('zz2011',
           '/ZZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'ZZ', 'ZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 4191045, 5.9, 7),
    mcsample('dytotautau2011',
           '/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'DYToTauTau', 'DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 19937479, 1666., 7),
    mcsample('dytomumu2011',
           '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'DYToMuMu', 'DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 29743564, 1666., 7),
    mcsample('wjetstolnu2011',
           '/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'WJetsToLNu', 'WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 81345381,  27770.0, 7 ),
    mcsample('qcd2011',
           '/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'QCD_Pt_20_MuEnrichedPt_15', 'QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13',25080241, 0.0002855 * 296600000., 7 )
]


class datasample(object):
    def __init__(self, name, dataset, ntuple, caf_folder_extd, global_tag, json_file, energy):
        self.name = name
        self.dataset = dataset
        self.ntuple = ntuple
        self.caf_folder_extd = caf_folder_extd
        self.global_tag = global_tag
        self.json_file = json_file
        self.energy = energy


DATASamples = [
    # single mu datasets
    # see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2012Analysis#Analysis_using_the_Golden_JSON_f
    datasample('SingleMuRun2012A',
               '/SingleMu/Run2012A-13Jul2012-v1/AOD', 'SingleMuRun2012A-13Jul2012-v1', 'SingleMuRun2012A-13Jul2012-v1',
               'FT_53_V6_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt',
               8),
    datasample('SingleMuRun2012A-recover-06Aug2012-v1',
               '/SingleMu/Run2012A-recover-06Aug2012-v1/AOD', 'SingleMuRun2012A-recover-06Aug2012-v1', 'SingleMuRun2012A-recover-06Aug2012-v1',
               'FT_53_V6C_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt',
               8),
    datasample('SingleMuRun2012B',
               '/SingleMu/Run2012B-13Jul2012-v1/AOD', 'SingleMuRun2012B-13Jul2012-v1', 'SingleMuRun2012B-13Jul2012-v1',
               'FT_53_V6_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt',
               8),
    datasample('SingleMuRun2012C-24Aug2012-v1',
               '/SingleMu/Run2012C-24Aug2012-v1/AOD', 'SingleMuRun2012C-24Aug2012-v1', 'SingleMuRun2012C-24Aug2012-v1',
               'FT_53_V10_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt',
               8),
    datasample('SingleMuRun2012C',
               '/SingleMu/Run2012C-PromptReco-v2/AOD', 'SingleMuRun2012C-PromptReco-v2', 'SingleMuRun2012C-PromptReco-v2',
               'GR_P_V42_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt',
               8),
    datasample('SingleMuRun2012D',
               '/SingleMu/Run2012D-PromptReco-v1/AOD', 'SingleMuRun2012D-PromptReco-v1', 'SingleMuRun2012D-PromptReco-v1',
               'GR_P_V42_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt',
               8),

    # single mu 2011 datasets
    datasample('SingleMuRun2011A',
               '/SingleMu/Run2011A-08Nov2011-v1/AOD', 'SingleMuRun2011A-08Nov2011-v1', 'SingleMuRun2011A-08Nov2011-v1',
               'FT44_V9_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt',
               7),
    datasample('SingleMuRun2011B',
               '/SingleMu/Run2011B-19Nov2011-v1/AOD', 'SingleMuRun2011B-19Nov2011-v1', 'SingleMuRun2011B-19Nov2011-v1',
               ' FT44_V11_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt',
               7),

    #double mu datasets
    datasample('DoubleMuRun2012A',
               '/DoubleMu/Run2012A-13Jul2012-v1/AOD', 'DoubleMuRun2012A-13Jul2012-v1', 'DoubleMuRun2012A-13Jul2012-v1',
               'FT_53_V6_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt',
               8),
    datasample('DoubleMuRun2012A-recover-06Aug2012-v1',
               '/DoubleMu/Run2012A-recover-06Aug2012-v1/AOD', 'DoubleMuRun2012A-recover-06Aug2012-v1', 'DoubleMuRun2012A-recover-06Aug2012-v1',
               'FT_53_V6C_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt',
               8),
    datasample('DoubleMuRun2012B',
               '/DoubleMu/Run2012B-13Jul2012-v4/AOD', 'DoubleMuRun2012B-13Jul2012-v4', 'DoubleMuRun2012B-13Jul2012-v4',
               'FT_53_V6_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt',
               8),
    datasample('DoubleMuRun2012C-24Aug2012-v1',
               '/DoubleMu/Run2012C-24Aug2012-v1/AOD', 'DoubleMuRun2012C-24Aug2012-v1', 'DoubleMuRun2012C-24Aug2012-v1',
               'FT_53_V10_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt',
               8),
    datasample('DoubleMuRun2012C',
               '/DoubleMu/Run2012C-PromptReco-v2/AOD', 'DoubleMuRun2012C-PromptReco-v2', 'DoubleMuRun2012C-PromptReco-v2',
               'GR_P_V42_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt',
               8),
    datasample('DoubleMuRun2012D',
               '/DoubleMu/Run2012D-PromptReco-v1/AOD', 'DoubleMuRun2012D-PromptReco-v1', 'DoubleMuRun2012D-PromptReco-v1',
               'GR_P_V42_AN3',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt',
               8)
]

# these are the triggers' strings
#  data
HLT_Single2011A = 'HLT_IsoMu24,HLT_Mu30'
HLT_Single2011B = 'HLT_IsoMu24_eta2p1,HLT_Mu40_eta2p1'
HLT_Single2012A = 'HLT_IsoMu24_eta2p1,HLT_Mu40_eta2p1'
HLT_Single = 'HLT_IsoMu24,HLT_Mu40_eta2p1'
HLT_Double = 'HLT_Mu17_Mu8,HLT_Mu17_TkMu8'
#  mc
HLT_MC_8TeV     = 'HLT_IsoMu24,HLT_Mu17_Mu8,HLT_Mu17_TkMu8' #For 2012 Samples
HLT_MC_7TeV     = 'HLT_IsoMu24,HLT_Mu17_Mu8' #For 2011 Samples
HLT_MC_Ele_7TeV = 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL'
HLT_MC_Ele_8TeV = 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL'
# ele
HLT_Ele2011A = 'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL'
HLT_Ele = 'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL'

# cmssw version
cmssw_8TeV = 'CMSSW_5_3_5' #2012
cmssw_7TeV = 'CMSSW_4_4_5' # 2011
tag   = 'V00-01-09'

# these are the general location paths:
# higgs/<CMSSW Version>/<Ntupler Tag>/Ntuples/<Data or MC>/
