#!/bin/bash
#PBS -N fastQC.AH  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=1:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-38


##Runs fastqc on paired end data listed in R1.names and R2.names

#run job in working directory
cd $PBS_O_WORKDIR 


#Load modules
module load apps/fastqc-0.11.5

#Define variables

NAME1=$(sed "${PBS_ARRAYID}q;d" R1.names)
NAME2=$(sed "${PBS_ARRAYID}q;d" R2.names)


##fastQC using array

sample_name=`echo ${NAME1} | awk -F "_R" '{print $1}'`
echo "[fastQC running for] $sample_name"
printf "\n"
echo "time fastqc --noextract ${NAME1} ${NAME2}" >> fastQC.log
time fastqc --noextract ${NAME1} ${NAME2}
