#!/usr/bin/env python

import sys
import os
import commands

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '  findMissingFile.py folder \"number of files\"\n'
   sys.exit()
    
# ------------------------------------------------------------------------------

if (len(sys.argv) != 3):
   Usage()

filesFolder = sys.argv[1]

files=commands.getoutput("ls %s | grep root | sort" % filesFolder)
files=files.split("\n")

nFiles  = sys.argv[2]

fileNames = files[0].split('_')
fileBase = fileNames[0]

#print fileBase

hypoList = []
for i in range(1,int(nFiles)+1):
   hypoList.append( '%s_%s' % (fileBase,i) ) 

foundList = []
for file in files:
##    #print file
   file = file.split('_')
   foundList.append( '%s_%s' % (file[0],file[1]) )

print "Hypo List has %s files" % len(hypoList)
print "Found List has %s files" % len(foundList)

diffList=list(set(hypoList)-set(foundList))

print "Missing Files are"
for file in diffList:
   print file
