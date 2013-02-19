#!/bin/tcsh

setenv ROOTSYS /data/uftrig01b/jhugon/libraries/root
setenv PATH $PATH:$ROOTSYS/bin
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH:$ROOTSYS/lib
setenv PYTHONPATH $PYTHONPATH:$ROOTSYS/lib:$ROOTSYS

#export PATH=/data/uftrig01b/jhugon/libraries/tskim08-02-02/bin:$PATH
#export TS_META_DATA=/data/uftrig01b/jhugon/libraries/tskim08-02-02/share/UFMetaData.txt
