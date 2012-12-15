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
                  help='specify the caf folder to write to (mandatory)',
                  dest='cafFolder',
                  action='store')

parser.add_option('-g', '--globaltag',
                  help='specify the globaltag (mandatory)',
                  dest='globalTag',
                  action='store')

parser.add_option('--hlt_muon', 
                  help='specify the list of muon triggers (mandatory)',
                  dest='triggerstring',
                  action='store')

parser.add_option('--hlt_ele', 
                  help='specify the ele trigger (mandatory)',
                  dest='triggerele',
                  action='store')

parser.add_option('-j', '--json',
                  help='specify the json file: isData flag will be set to true',
                  dest='jsonFile',
                  action='store')

parser.add_option("-e",'--energy',
                  help='Flag to Specify LHC Energy of Dataset',
                  dest='energy',
                  action='store')

# get the options
(opts, args) = parser.parse_args()


# Sanity check: Making sure all mandatory options appeared
mandatories = ['CMSSW', 'dataset', 'ntupleName', 'cafFolder','globalTag', 'triggerstring', 'triggerele']
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

is2011=False
if opts.energy=="7" or opts.energy=="7.0":
  is2011=True

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
          "| sed -e \'s/digiovan/%s/g\' " \
          "| sed -e \'s/yourDataset/%s/g\' " \
          "| sed -e \'s/yourNtuple/%s/g\' " \
          "| sed -e \'s/yourCafFolder/%s/g\' " \
          "| sed -e \'s/yourJsonFile/%s/g\' " \
          " > %s/crab.cfg" % (crab_file,
                              os.environ["USER"],
                              opts.dataset.replace("/","\/"),
                              opts.ntupleName,
                              opts.cafFolder.replace("/","\/"), 
                              opts.jsonFile.split("/")[-1],
                              localFolder) )

   print '  --> copy json file'
   os.system( 'cp %s %s/' % (opts.jsonFile,localFolder) )

else:
    os.system("cat crabTemplate/%s " \
          "| sed -e \'s/digiovan/%s/g\' " \
              "| sed -e \'s/yourDataset/%s/g\' " \
              "| sed -e \'s/yourNtuple/%s/g\' " \
              "| sed -e \'s/yourCafFolder/%s/g\' " \
              " > %s/crab.cfg" % (crab_file,
                                  os.environ["USER"],
                                  opts.dataset.replace("/","\/"),
                                  opts.ntupleName,
                                  opts.cafFolder.replace("/","\/"), 
                                  localFolder) )

    
print ' b) UFDiMuonAnalyzer.py'

if isData:
    os.system("cat crabTemplate/%s " \
              "| sed -e \'s/yourNtuple/%s/g\' " \
              "| sed -e \'s/GLOBALTAG/%s/g\' " \
              "| sed -e \'s/thisIsData = False/thisIsData = True/g\' " \
              "| sed -e \'s/TRIGGERLIST/%s/g\' " \
              "| sed -e \'s/TRIGGERELE/%s/g\' " \
              "| sed -e \'s/thisIs2011 = False/thisIs2011 = %s/g\' " \
              "> %s/UFDiMuonAnalyzer.py"
              % (cmssw_py_file,
                 opts.ntupleName,
                 opts.globalTag,
                 opts.triggerstring.replace(",","\",\""),
                 opts.triggerele.replace(",","\",\""),
                 is2011,
                 localFolder)
              )
else:
    os.system("cat crabTemplate/%s " \
              "| sed -e \'s/yourNtuple/%s/g\' " \
              "| sed -e \'s/GLOBALTAG/%s/g\' " \
              "| sed -e \'s/TRIGGERLIST/%s/g\' " \
              "| sed -e \'s/TRIGGERELE/%s/g\' " \
              "| sed -e \'s/thisIs2011 = False/thisIs2011 = %s/g\' " \
              "> %s/UFDiMuonAnalyzer.py"
              % (cmssw_py_file,
                 opts.ntupleName,
                 opts.globalTag,
                 opts.triggerstring.replace(",","\",\""),
                 opts.triggerele.replace(",","\",\""),
                 is2011,
                 localFolder)
              )
    

print ' c) howToRunCrab'
os.system("cat crabTemplate/howToRunCrab " \
          "| sed -e 's/CMSSW/cms%s/g' " \
          "> %s/howToRunCrab"
          % ( opts.CMSSW.replace("CMSSW_",""),
              localFolder)
          )

print ' d) README'
os.system("cp crabTemplate/README %s/" % (localFolder) )


print '\n\nDone: you are ready to submit your jobs'
print     '---------------------------------------'
print 'cd', localFolder
print 'crab -create -submit\n'

