#!/usr/bin/env python

import sys
import os
import commands

# for the options
import optparse
parser = optparse.OptionParser()

# define the options
parser.add_option('-c', '--cmssw',
                  help='specify the CMSSW version (mandatory)',
                  dest='CMSSW',
                  action='store')

parser.add_option('-d', '--dataset',
                  help='specify the dataset to run over (mandatory)',
                  dest='dataset',
                  action='store')

parser.add_option('-n', '--ntuple',
                  help='specify the ntuple name (mandatory)',
                  dest='ntupleName',
                  action='store')

parser.add_option('-f', '--caffolder',
                  help='specify the dataset to run over (mandatory)',
                  dest='cafFolder',
                  action='store')

parser.add_option('-j', '--json',
                  help='specify the json file: isData flag will be set to true',
                  dest='jsonFile',
                  action='store')

# get the options
(opts, args) = parser.parse_args()


# Sanity check: Making sure all mandatory options appeared
mandatories = ['CMSSW', 'dataset', 'ntupleName', 'cafFolder']
for m in mandatories:
    if not opts.__dict__[m]:
        print "at least mandatory option is missing\n"
        parser.print_help()
        exit(-1)

isData=False
crab_file = 'crabMC.cfg'
cmssw_py_file = 'UFDiMuonAnalyzer.py'

if opts.jsonFile is not None:
   print '\nRunning over a DATA dataset'
   isData=True
   crab_file = 'crabData.cfg'
else:
   print '\nRunning over a Monte Carlo dataset'


## print opts.cafFolder.split("/")

# create the cafPath: user dependent
username = commands.getoutput('whoami')
cafPath  = '/store/user/'
cafPath += username

# Assuming the naming of the caffolder is cmssw/blabla/Ntuples/MC(Data)/cafFolder
localFolder=''
localFolder += opts.cafFolder.split("/")[-3]
localFolder += opts.cafFolder.split("/")[-2]
localFolder += opts.cafFolder.split("/")[-1]
#localFolder.join(seq)


print '1) creating caf Folder %s\n' % opts.cafFolder
os.system('cmsMkdir %s/%s' % (cafPath,opts.cafFolder) )

print '2) creating local Folder %s\n' % localFolder
os.system('mkdir %s' % localFolder)


print 'Prepare the configuration files:'
#print opts.cafFolder.replace("/","\/")
print ' a) crab.cfg'

if isData:
   os.system("cat crabTemplate/%s " \
          "| sed -e \'s/yourDataset/%s/g\' " \
          "| sed -e \'s/yourNtuple/%s/g\' " \
          "| sed -e \'s/yourCafFolder/%s/g\' " \
          "| sed -e \'s/yourJsonFile/%s/g\' " \
          " > %s/crab.cfg" % (crab_file,
                              opts.dataset.replace("/","\/"),
                              opts.ntupleName,
                              opts.cafFolder.replace("/","\/"), 
                              opts.jsonFile.split("/")[-1],
                              localFolder) )

   print '  --> copy json file'
   os.system( 'cp %s %s/' % (opts.jsonFile,localFolder) )

else:
    os.system("cat crabTemplate/%s " \
              "| sed -e \'s/yourDataset/%s/g\' " \
              "| sed -e \'s/yourNtuple/%s/g\' " \
              "| sed -e \'s/yourCafFolder/%s/g\' " \
              " > %s/crab.cfg" % (crab_file,
                                  opts.dataset.replace("/","\/"),
                                  opts.ntupleName,
                                  opts.cafFolder.replace("/","\/"), 
                                  localFolder) )

    
print ' b) UFDiMuonAnalyzer.py'

if isData:
    os.system("cat crabTemplate/%s " \
              "| sed -e \'s/yourNtuple/%s/g\' " \
              "| sed -e \'s/thisIsData = False/thisIsData = True/g\' " \
              "> %s/UFDiMuonAnalyzer.py"
              % (cmssw_py_file, opts.ntupleName, localFolder) )
else:
    os.system('cat crabTemplate/%s | sed -e \'s/yourNtuple/%s/g\' > %s/UFDiMuonAnalyzer.py'
              % (cmssw_py_file, opts.ntupleName, localFolder) )
    

print ' c) howToRunCrab'
os.system("cat crabTemplate/howToRunCrab " \
          "| sed -e 's/CMSSW/cms%s/g' > %s/howToRunCrab" % ( opts.CMSSW.replace("CMSSW_",""),
                                                             localFolder)
          )

print ' d) README'
os.system("cp crabTemplate/README %s/" % (localFolder) )


print '\n\nDone: you are ready to submit your jobs'
print     '---------------------------------------'
print 'cd', localFolder
print 'crab -create -submit\n'

