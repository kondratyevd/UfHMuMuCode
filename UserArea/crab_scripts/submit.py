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
            command += './createCrab.py -c %s -d %s -n %s -g %s -f %s ' \
                       % (Samples.cmssw,
                          s.dataset,
                          s.ntuple,
                          s.global_tag,
                          Samples.caf_pathMC+s.caf_folder_extd)
                         
            
else:
    for s in Samples.DATASamples:
        if (s.name == opts.sample):
            command += './createCrab.py -c %s -d %s -n %s -g %s -f %s -j %s ' \
                       % (Samples.cmssw,
                          s.dataset,
                          s.ntuple,
                          s.global_tag,
                          Samples.caf_pathData+s.caf_folder_extd,
                          s.json_file) 


if (isMC):
    command += '--hlt %s' % Samples.HLT_MC
else:
    if ('DoubleMu' in command):
        command += '--hlt %s' % Samples.HLT_Double
    elif ('SingleMu' and '2012A' in command):
        command += '--hlt %s' % Samples.HLT_Single2012A
    else:
        command += '--hlt %s' % Samples.HLT_Single
    
print '\n %s \n' % command
