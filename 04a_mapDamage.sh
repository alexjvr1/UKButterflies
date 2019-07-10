#!/bin/bash
# -----------------------------------------

# (c) Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified: 03/10/2018 15:39:48

# Description:
# MapDamage to assess bias in transitions/transversions in the dataset. And to recalibrate .bam files


#PBS -N D2.MapDmg  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-87

#run job in working directory
cd $PBS_O_WORKDIR 

# Load modules
module load languages/python-2.7.6
module load languages/R-3.0.2

# Define variables
mapDamage="$HOME/.local/bin/mapDamage"
RefSeq="../RefGenome/*fasta"
NAME=$(sed "${PBS_ARRAYID}q;d" bamfiles.names)

##Script
echo "mapDamage2 started" >> map.log
echo "---------------" >> map.log

sample_name=`echo ${NAME} | awk -F "_L00" '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time $mapDamage --merge-reference-sequences -i ${NAME} -r $RefSeq >> mapdamage.log"
time $mapDamage --merge-reference-sequences -i ${NAME} -r $RefSeq
