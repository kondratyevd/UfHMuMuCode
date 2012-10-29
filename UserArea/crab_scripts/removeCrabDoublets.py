#!/usr/bin/env python

import sys
import os
import commands

# ------------------------------------------------------------------------------

def Usage():
   print 'Wrong syntax: \n'
   print '   ./removeCrabDoublets.py folder [action(0)]\n'
   print 'if action!=0 then remove the files otherwise just display doublets\n'
   sys.exit()
    
# ------------------------------------------------------------------------------


filesFolder  = ''
action = 0

if (len(sys.argv) == 2):
   filesFolder =  sys.argv[1]
elif (len(sys.argv) == 3):
   filesFolder =  sys.argv[1]
   action      =  sys.argv[2]
else:
   Usage()


files=commands.getoutput("ls %s | grep root | sort" % filesFolder)
files=files.split("\n")

## ------------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------------
previousFile = 0
index = 0

nFiles = len(files)

for index in range(nFiles):
   if (index < 1):
      continue

   previousFile = files[index-1].split('_')
   currentFile = files[index].split('_')

   #print currentFile
   #if ( (previousFile[0] == currentFile[0]) and
   #     (previousFile[1] == currentFile[1]) ):

   #print previousFile[-3], currentFile[-3]
   
   if (previousFile[-3] == currentFile[-3]):
   #     (previousFile[1] == currentFile[1]) ):
      #and
      #  (previousFile[2] == currentFile[2])   ):
   #if ( (previousFile[2] == currentFile[2]) and
   #     (previousFile[3] == currentFile[3])   ):

      previousFileSize = commands.getoutput("ls -l %s/%s | awk '{print $5}'" % (filesFolder, files[index-1]) )
      currentFileSize  = commands.getoutput("ls -l %s/%s | awk '{print $5}'" % (filesFolder, files[index])   )
  #    previousFileSize = commands.getoutput("ls -l %s/%s | awk '{print $6}'" % (filesFolder, files[index-1]) )
  #    currentFileSize  = commands.getoutput("ls -l %s/%s | awk '{print $6}'" % (filesFolder, files[index])   )
      
      print previousFileSize, files[index-1]
      print currentFileSize,  files[index]

      if (currentFileSize >=  previousFileSize):
         print "rm %s/%s" % (filesFolder, files[index-1])
         if (action !=0):
            os.system("rm %s/%s" % (filesFolder, files[index-1]))
      else:
         print "rm %s/%s" % (filesFolder, files[index])
         if (action !=0):
            os.system("rm %s/%s" % (filesFolder, files[index]))
                      
