# Standard import for CMSSW files
import FWCore.ParameterSet.Config as cms

# Import VarParsing
from FWCore.ParameterSet.VarParsing import VarParsing

# implement input from the shell
options = VarParsing ('analysis')


options.register ('whichSample',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "Choose which sample to process from the samples array.")

options.parseArguments()

# Must define a CMSSW process named process
process = cms.Process("ReadDatasets");
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# The files to read in
readFiles = cms.untracked.vstring()


from Samples_v3 import background as s

print ""
print "Copying from " + s[options.whichSample].name + ", " + s[options.whichSample].dir
print ""

readFiles.extend(s[options.whichSample].files);

print "=== Files in the sample ==="
print readFiles
print ""

# Which files to get
process.source = cms.Source ("PoolSource", fileNames = readFiles)

# Save the events to a file
process.Out = cms.OutputModule("PoolOutputModule",
#         outputCommands = cms.untracked.vstring("drop *", "keep recoTracks_*_*_*"),
         fileName = cms.untracked.string (s[options.whichSample].name+".root")
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)
