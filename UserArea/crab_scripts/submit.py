#!/usr/bin/env python

import sys
import os
import commands
import Samples

# for the options
import optparse
parser = optparse.OptionParser()

# define the options
parser.add_option('-s', '--sample',
                     help='specify the sample to run over (mandatory)',
                     dest='sample',
                     action='store')

# get the options
(opts, args) = parser.parse_args()


# Sanity check: Making sure all mandatory options appeared
mandatories = ['sample']
for m in mandatories:
    if not opts.__dict__[m]:
        print "at least a mandatory option is missing\n"
        parser.print_help()
        exit(-1)

isMC = True
if ( 'SingleMu' in opts.sample):
    isMC = False
if ( 'DoubleMu' in opts.sample):
    isMC = False


# define the command to send
command = ''          
if isMC:          
    for s in Samples.MCSamples:
        if (s.name == opts.sample):
            hltString = ""
            hltEle = ""
            cmsswString = ""
            if s.energy==7:
              hltString = Samples.HLT_MC_7TeV
              hltEle = Samples.HLT_MC_Ele_7TeV
              cmsswString = Samples.cmssw_7TeV
            elif s.energy==8:
              hltString = Samples.HLT_MC_8TeV
              hltEle = Samples.HLT_MC_Ele_8TeV
              cmsswString = Samples.cmssw_8TeV
            else:
              print("Error: No MC HLT Trigger or CMSSW version for Energy: '%s' TeV" % s.energy)
              sys.exit(1)
            cafFolder = "higgs/%s/%s/Ntuples/%s/" % (cmsswString,Samples.tag,"MC")
            command += './createCrab.py -c %s -d %s -n %s -g %s -e %s -f %s --hlt_muon %s --hlt_ele %s ' \
                       % (cmsswString,
                          s.dataset,
                          s.ntuple,
                          s.global_tag,
                          s.energy,
                          cafFolder+s.caf_folder_extd,
                          hltString,
                          hltEle)
            
else:
    for s in Samples.DATASamples:
        if (s.name == opts.sample):
            cmsswString = ""
            if s.energy==7:
              cmsswString = Samples.cmssw_7TeV
            elif s.energy==8:
              cmsswString = Samples.cmssw_8TeV
            else:
              print("Error: No Data CMSSW version for Energy: '%s' TeV" % s.energy)
              sys.exit(1)
            cafFolder = "higgs/%s/%s/Ntuples/%s/" % (cmsswString,Samples.tag,"Data")
            command += './createCrab.py -c %s -d %s -n %s -g %s -e %s -f %s -j %s ' \
                       % (cmsswString,
                          s.dataset,
                          s.ntuple,
                          s.global_tag,
                          s.energy,
                          cafFolder+s.caf_folder_extd,
                          s.json_file) 

    if ('DoubleMu' in command):
        command += '--hlt_muon %s' % Samples.HLT_Double
    elif ('SingleMu' and '2011A' in command):
        command += '--hlt_muon %s --hlt_ele %s' % (Samples.HLT_Single2011A, Samples.HLT_Ele2011A)
    elif ('SingleMu' and '2011B' in command):
        command += '--hlt_muon %s --hlt_ele %s' % (Samples.HLT_Single2011B, Samples.HLT_Ele)
    elif ('SingleMu' and '2012A' in command):
        command += '--hlt_muon %s --hlt_ele %s' % (Samples.HLT_Single2012A, Samples.HLT_Ele)
    else:
        command += '--hlt_muon %s --hlt_ele %s' % (Samples.HLT_Single, Samples.HLT_Ele)

    
print '\n %s \n' % command
