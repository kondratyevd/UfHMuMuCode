#!/usr/bin/env python

import os

class mcsample(object):
    def __init__(self, name, dataset, ntuple, caf_folder_extd, global_tag, nevents, cross_section):
        self.name = name
        self.dataset = dataset
        self.ntuple = ntuple
        self.caf_folder_extd = caf_folder_extd
        self.global_tag = global_tag
        self.nevents = nevents
        self.cross_section = cross_section


# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV for xsecs (all below in pb)
# or from PREP
MCSamples = [
    mcsample('dyjetstoll',
           '/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYJetsToLL', 'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 30459503, 3503.71),
    mcsample('ttjets',
           '/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
           'TTJets', 'TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3',
           'START53_V7E', 6923750, 225.197),
    mcsample('ww',
           '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'WW', 'WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 9900431, 54.838),
    mcsample('wz',
           '/WZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'WZ', 'WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 10000283, 33.21),
    mcsample('zz',
           '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'ZZ', 'ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 9799908, 17.654),
    mcsample('dytotautau',
           '/DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYToTauTau', 'DYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 3295238, 5745.25/3.),
#    mcsample('dytomumu',
#           '/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
#           'DYToMuMu', 'DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1',
#           'START53_V7E', 3293740, 5745.25/3.),
    mcsample('dytomumu',
           '/DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYToMuMu', 'DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 3293740, 5745.25/3.),
    mcsample('wjetstolnu',
           '/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
           'WJetsToLNu', 'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2',
           'START53_V7E', 57330800, 30400.0),
    mcsample('qcd',
           '/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
           'QCD_Pt_20_MuEnrichedPt_15', 'QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3',
           'START53_V7E',20284602, 3.64E8*3.7E-4),

    # special MC datasets
    mcsample('dyjetstoll_ptz-100',
           '/DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
           'DYJetsToLL_PtZ-100', 'DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2',
           'START53_V7E', 4930773, 93.8),
    mcsample('dyjetstoll_ptz-50to70',
           '/DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
           'DYJetsToLL_PtZ-50to70', 'DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1',
           'START53_V7E', 1313395, 52.31),
    mcsample('dyjetstoll_ptz-70to100',
           '/DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
           'DYJetsToLL_PtZ-70to100', 'DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2',
           'START53_V7E', 2662137, 34.1),

     # 2011 Samples
    mcsample('dyjetstoll2011',
           '/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'DYJetsToLL', 'DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 36264432, 3048.),
    mcsample('ttjets2011',
           '/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'TTJets', 'TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 59444088, 157.5),
    mcsample('ww2011',
           '/WW_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'WW', 'WW_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 4225916,  43.0),
    mcsample('wz2011',
           '/WZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'WZ', 'WZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 4265243, 18.2),
    mcsample('zz2011',
           '/ZZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'ZZ', 'ZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 4191045, 5.9),
    mcsample('dytotautau2011',
           '/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'DYToTauTau', 'DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 19937479, 1666.),
    mcsample('dytomumu2011',
           '/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'DYToMuMu', 'DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 29743564, 1666.),
    mcsample('wjetstolnu2011',
           '/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'WJetsToLNu', 'WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13', 81345381,  27770.0 ),
    mcsample('qcd2011',
           '/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START44_V9B-v1/AODSIM',
           'QCD_Pt_20_MuEnrichedPt_15', 'QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_Fall11-PU_S6_START44_V9B-v1',
           'START44_V13',25080241, 0.0002855 * 296600000. )
]


class datasample(object):
    def __init__(self, name, dataset, ntuple, caf_folder_extd, global_tag, json_file):
        self.name = name
        self.dataset = dataset
        self.ntuple = ntuple
        self.caf_folder_extd = caf_folder_extd
        self.global_tag = global_tag
        self.json_file = json_file


DATASamples = [
    # single mu datasets
    datasample('SingleMuRun2012A',
               '/SingleMu/Run2012A-13Jul2012-v1/AOD', 'SingleMuRun2012A-13Jul2012-v1', 'SingleMuRun2012A-13Jul2012-v1',
               'FT_53_V6_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'),
    datasample('SingleMuRun2012B',
               '/SingleMu/Run2012B-13Jul2012-v1/AOD', 'SingleMuRun2012B-13Jul2012-v1', 'SingleMuRun2012B-13Jul2012-v1',
               'FT_53_V6_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'),
    datasample('SingleMuRun2012C-24Aug2012-v1',
               '/SingleMu/Run2012C-24Aug2012-v1/AOD', 'SingleMuRun2012C-24Aug2012-v1', 'SingleMuRun2012C-24Aug2012-v1',
               'FT_53_V10_AN2',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'),
    datasample('SingleMuRun2012C',
               '/SingleMu/Run2012C-PromptReco-v2/AOD', 'SingleMuRun2012C-PromptReco-v2', 'SingleMuRun2012C-PromptReco-v2',
               'GR_P_V41_AN2',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt'),

    # single mu 2011 datasets
    datasample('SingleMuRun2011A',
               '/SingleMu/Run2011A-08Nov2011-v1/AOD', 'SingleMuRun2011A-08Nov2011-v1', 'SingleMuRun2011A-08Nov2011-v1',
               'FT44_V9_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt'),
    datasample('SingleMuRun2011B',
               '/SingleMu/Run2011B-19Nov2011-v1/AOD', 'SingleMuRun2011B-19Nov2011-v1', 'SingleMuRun2011B-19Nov2011-v1',
               ' FT44_V11_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt'),

    #double mu datasets
    datasample('DoubleMuRun2012A',
               '/DoubleMu/Run2012A-13Jul2012-v1/AOD', 'DoubleMuRun2012A-13Jul2012-v1', 'DoubleMuRun2012A-13Jul2012-v1',
               'FT_53_V6_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'),
    datasample('DoubleMuRun2012B',
               '/DoubleMu/Run2012B-13Jul2012-v4/AOD', 'DoubleMuRun2012B-13Jul2012-v4', 'DoubleMuRun2012B-13Jul2012-v4',
               'FT_53_V6_AN1',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'),
    datasample('DoubleMuRun2012C-24Aug2012-v1',
               '/DoubleMu/Run2012C-24Aug2012-v1/AOD', 'DoubleMuRun2012C-24Aug2012-v1', 'DoubleMuRun2012C-24Aug2012-v1',
               'FT_53_V10_AN2',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'),
    datasample('DoubleMuRun2012C',
               '/DoubleMu/Run2012C-PromptReco-v2/AOD', 'DoubleMuRun2012C-PromptReco-v2', 'DoubleMuRun2012C-PromptReco-v2',
               'GR_P_V41_AN2',
               '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt')
]

# these are the triggers' strings
#  data
HLT_Single2011A = 'HLT_IsoMu24,HLT_Mu40'
HLT_Single2011B = 'HLT_IsoMu24_eta2p1,HLT_Mu40_eta2p1'
HLT_Single2012A = 'HLT_IsoMu24_eta2p1,HLT_Mu40_eta2p1'
HLT_Single = 'HLT_IsoMu24,HLT_Mu40_eta2p1'
HLT_Double = 'HLT_Mu17_Mu8,HLT_Mu17_TkMu8'
#  mc
HLT_MC     = 'HLT_IsoMu24,HLT_Mu17_Mu8,HLT_Mu17_TkMu8' #For 2012 Samples
#HLT_MC     = 'HLT_IsoMu24,HLT_Mu17_Mu8,HLT_Mu17_TkMu8' #For 2011 Samples

# cmssw version
cmssw = 'CMSSW_5_3_3_patch3' #2012
#cmssw = 'CMSSW_4_4_5' # 2011
tag   = 'V00-01-01'

# these are the general location paths
caf_pathMC   = 'higgs/%s/%s/Ntuples/MC/'   % (cmssw,tag)
caf_pathData = 'higgs/%s/%s/Ntuples/Data/' % (cmssw,tag) 
