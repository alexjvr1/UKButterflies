#!/bin/bash
#PBS -N D2.mod.concat.idx  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)

#define variables

OUTNAME="D2.modExp.regionsALL2.unfolded"

#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd


#Merge files across all regions within a sample set
time ~/bin/angsd/misc/realSFS cat -outnames $OUTNAME *idx
