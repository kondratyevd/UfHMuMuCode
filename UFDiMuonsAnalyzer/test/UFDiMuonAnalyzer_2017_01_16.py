# =============================================================#
# UFDiMuonsAnalyzer                                             #
# =============================================================#
# Makes stage1 trees.                                          #
# Adds a cleaner vector of jets to each event.                 #
# Originally Made by Justin Hugon. Edited by Andrew Carnes.    #
#                                                              #
################################################################ 

# /////////////////////////////////////////////////////////////
# Load some things
# /////////////////////////////////////////////////////////////

#Add new comment to test push

import FWCore.ParameterSet.Config as cms

process = cms.Process("UFDiMuonsAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

## Geometry and Detector Conditions (needed for a few patTuple production steps)

## Correct geometry to use?  What about GeometryExtended2016_cff or GeometryExtended2016Reco_cff? - AWB 16.01.17
## https://github.com/cms-sw/cmssw/tree/CMSSW_8_0_X/Configuration/Geometry/python
process.load("Configuration.Geometry.GeometryIdeal_cff")  

# ## Geometry according to Tim Cox, used by Jia Fu Low
# ##   https://indico.cern.ch/event/588469/contributions/2372672/subcontributions/211968/attachments/1371248/2079893/
# ##   2016-11-14_coordinate_conversion_v1.pdf
# from Configuration.AlCa.autoCond import autoCond
# process.load('Configuration.StandardSequences.GeometryDB_cff')
# process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")
# process.preferFakeAlign = cms.ESPrefer("FakeAlignmentSource")


# /////////////////////////////////////////////////////////////
# Get a sample from our collection of samples
# /////////////////////////////////////////////////////////////

from python.Samples_Moriond17 import SingleMu_2016C as samp
# from python.Samples_Moriond17 import H2Mu_gg as samp
# from python.Samples_Moriond17 import ZJets_AMC as samp
# from python.Samples_Moriond17 import ZJets_MG_HT_2500_inf as samp

if samp.isData:
    print '\nRunning over data sample %s' % samp.name
else:
    print '\nRunning over MC sample %s' % samp.name
print '  * From DAS: %s' % samp.DAS

# /////////////////////////////////////////////////////////////
# global tag, automatically retrieved from the imported sample
# /////////////////////////////////////////////////////////////

print '\nLoading Global Tag: ' + samp.GT
process.GlobalTag.globaltag = samp.GT

# # /////////////////////////////////
# # Additional jet energy corrections
# # /////////////////////////////////

# ## See https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
# ## and https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections

# from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
# process.jec = cms.ESSource('PoolDBESSource',
#     CondDBSetup,
#     connect = cms.string('sqlite:data/JEC/Spring16_23Sep2016V2_MC.db'),
#     toGet = cms.VPSet(
#         # cms.PSet(
#         #     record = cms.string('JetCorrectionsRecord'),
#         #     tag    = cms.string('JetCorrectorParametersCollection_Fall15_V2_DATA_AK4PFchs'),
#         #     label  = cms.untracked.string('AK4PFchs')
#         # ),
#         cms.PSet(
#             record = cms.string('JetCorrectionsRecord'),
#             # record = cms.string('JetCorrectorParametersCollection'),  ## Produces run-time error
#             tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016V2_MC_AK4PFchs'),
#             # label  = cms.untracked.string('AK4PFchs')
#             label  = cms.untracked.string('slimmedJets')
#         ),
#         # ...and so on for all jet types you need
#     )
# )

# # Add an ESPrefer to override JEC that might be available from the global tag
# process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')


# /////////////////////////////////////////////////////////////
# ------------ PoolSource -------------
# /////////////////////////////////////////////////////////////
readFiles = cms.untracked.vstring();
# Get list of files from the sample we loaded
readFiles.extend(samp.files);

# readFiles.extend(['file:dy_jetsToLL.root']);

# readFiles.extend(['root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root']);

# ## From /GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM
# readFiles.extend(['/store/user/abrinke1/HiggsToMuMu/samples/GluGlu_HToMuMu_M125_13TeV_powheg_pythia8/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/12B931FE-CD3A-E611-9844-0025905C3D98.root'])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source("PoolSource",fileNames = readFiles)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()

# use a JSON file when locally executing cmsRun
if samp.isData:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = samp.JSON).getVLuminosityBlockRange()
    # process.source.lumisToProcess = LumiList.LumiList(filename = 'data/JSON/bad_evt.txt').getVLuminosityBlockRange()


# /////////////////////////////////////////////////////////////
# Save output with TFileService
# /////////////////////////////////////////////////////////////

process.TFileService = cms.Service("TFileService", fileName = cms.string("SingleMu_2016C_test.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("GluGlu_HToMuMu_M125_GEN_test.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("ZJets_AMC_GEN_Roch_test.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("ZJets_MG_HT_2500_inf_test_100.root") )

# /////////////////////////////////////////////////////////////
# Load UFDiMuonsAnalyzer
# /////////////////////////////////////////////////////////////

if samp.isData:
  process.load("UfHMuMuCode.UFDiMuonsAnalyzer.UFDiMuonsAnalyzer_cff")
else:
  process.load("UfHMuMuCode.UFDiMuonsAnalyzer.UFDiMuonsAnalyzer_MC_cff")

process.dimuons = process.DiMuons.clone()
# process.dimuons.jetsTag    = cms.InputTag("cleanJets")
process.dimuons.isVerbose  = cms.untracked.bool(False)
process.dimuons.doSys      = cms.bool(True)
process.dimuons.doSys_KaMu = cms.bool(False)
process.dimuons.doSys_Roch = cms.bool(True)
process.dimuons.slimOut    = cms.bool(True)

# # /////////////////////////////////////////////////////////////
# # Bad event flags
# # /////////////////////////////////////////////////////////////

# ## Following https://github.com/MiT-HEP/NeroProducer/blob/master/Nero/test/testNero.py
# ## See also https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2786/2.html
# process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
# process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
# process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
# process.BadPFMuonFilter.taggingMode = cms.bool(True) ## Accept all events, will just flag

# /////////////////////////////////////////////////////////////
# Electron Cut Based IDs
# /////////////////////////////////////////////////////////////

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

## First need to run: git cms-merge-topic ikrav:egm_id_80X_v1
## https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_8_0
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# /////////////////////////////////////////////////////////////
# Updated Jet Energy Scale corrections
# /////////////////////////////////////////////////////////////

## Following https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
##   - Last check that procedure was up-to-date: March 10, 2017 (AWB)
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

print 'samp.isData = %d' % samp.isData

if samp.isData:
    JEC_to_apply = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
else:
    JEC_to_apply = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', JEC_to_apply, 'None')
    )

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)


# /////////////////////////////////////////////////////////////
# Corrected MET (and jets?)
# /////////////////////////////////////////////////////////////

## Following https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription - AWB 01.03.17
##   - Last check that procedure was up-to-date: March 10, 2017 (AWB)

## First need to run: git cms-merge-topic cms-met:METRecipe_8020 -u
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

## If you only want to re-correct and get the proper uncertainties
runMetCorAndUncFromMiniAOD(process, isData=samp.isData)


# # /////////////////////////////////////////////////////////////
# # Save output tree
# # /////////////////////////////////////////////////////////////

# outCommands = cms.untracked.vstring('keep *')

# process.treeOut = cms.OutputModule("PoolOutputModule",
#                                    fileName = cms.untracked.string("GluGlu_HToMuMu_M125_JEC_10k_tree.root"),
#                                    outputCommands = outCommands
#                                    )

# process.treeOut_step = cms.EndPath(process.treeOut) ## Keep output tree
    

# /////////////////////////////////////////////////////////////
# Set the order of operations
# /////////////////////////////////////////////////////////////
    
print 'About to run the process path'

process.p = cms.Path( # process.BadPFMuonFilter *
                      process.egmGsfElectronIDSequence * 
                      process.jecSequence *
                      process.fullPatMetSequence *
                      process.dimuons )
# process.schedule = cms.Schedule(process.p, process.treeOut_step)

# ## Kill electrons for now - AWB 10.11.16
# process.p = cms.Path(process.dimuons)
