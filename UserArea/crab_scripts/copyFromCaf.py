#!/usr/bin/env python

import sys
import os
import commands

# ------------------------------------------------------------------------------
def Usage():
   print 'Wrong syntax: \n'
   print '   ./copyFromCaf.py CafFolder LocalFolder\n'
   sys.exit()

# ------------------------------------------------------------------------------

if (len(sys.argv) != 3):
   print 'you have passed %d arguments' % int(len(sys.argv)-1) 
   Usage()

CafFolder = '%s' % sys.argv[1]
LocalFolder  = sys.argv[2]

print 'Copying from folder: '
print  ' %s' % CafFolder 

if CafFolder[-1]=="/":
   CafFolder = CafFolder[0:-1]

castorFiles=commands.getoutput("cmsLs %s | grep root | awk '{print $5}' | sort" % CafFolder)
castorFiles=castorFiles.split('\n')

localFiles=commands.getoutput("ls %s | grep root | sort | awk '{print \"%s/\"$1}'" % (LocalFolder, CafFolder) )
localFiles = localFiles.split('\n')

## Compute The Difference
diffList=list(set(castorFiles)-set(localFiles))
 
files = []
for file in diffList:
   files.append(file.strip())

# ------------------------------------------------------------------------------
# Copy
# ------------------------------------------------------------------------------
for fileName in files:
    print 'cmsStage %s %s' % (fileName, LocalFolder)
    os.system('cmsStage %s %s' % (fileName, LocalFolder))
