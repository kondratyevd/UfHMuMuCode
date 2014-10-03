#! /bin/bash
#
#PBS -r n

##Job settings
#PBS -N Powheg_GluGludefault_125
#PBS -o Powheg_GluGludefault_125
#PBS -e Powheg_GluGludefault_125
#PBS -m a
#PBS -M jhugon@cern.ch

#Multiple Job Submission:
#Jobs will have a variable called $PBS_ARRAYID
#that will be one of the following numbers

##Job Resources
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=3000mb

# initialize environment for worker
cd /home/jhugon/scratchhpc/powheg/powheg20120804/gg_H/my8TeV125GeVHiggs/
source /home/jhugon/setupGen.sh


#####################################
#####Begin Real Job##################
#####################################

echo `date`
echo "Running Powheg:"

/home/jhugon/scratchhpc/powheg/powheg20120804/gg_H/my8TeV125GeVHiggs//../pwhg_main >& logfile 

echo `date`
echo "Done!"

