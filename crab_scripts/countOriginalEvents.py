#!/usr/bin/env python

import sys
import os
import commands

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./countOriginalEvents.py folder crab folder\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

if (len(sys.argv) != 3):
   Usage()

filesFolder = sys.argv[1]

files=commands.getoutput("ls %s | grep root | sort" % filesFolder)
files=files.split("\n")

crabFolder  = sys.argv[2]

nEvents = 0
for file in files:
   #print file
   file = file.split('_')
   #print file[-1], file[-2], file[-3], file[-4], file[-5]
   #print "grep \"TrigReport Events\" %s/res/CMSSW_%s.stdout" % (crabFolder, file[-3]) 
   os.system("grep \"TrigReport Events\" %s/res/CMSSW_%s.stdout" % (crabFolder, file[-3]) )

   ncands = commands.getoutput("grep \"TrigReport Events\" %s/res/CMSSW_%s.stdout | awk '{print $5}'" % (crabFolder, file[-3]) )
   nEvents += int(ncands)

print "number of events: %s" % int(nEvents)   
