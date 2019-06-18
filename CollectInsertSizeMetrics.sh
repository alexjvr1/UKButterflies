#!/bin/bash
#PBS -N TA.InsertSize_Museum  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.
#PBS -j oe  #concatenates error and output files (with prefix job1)


#run job in working directory
cd $PBS_O_WORKDIR
pwd
#cd WPA #uncomment when running locally
#pwd

#Load modules
module load apps/picard-2.20.0
module load languages/java-jdk-1.8.1-44

#Define variables

for i in $(ls *bam); do java -jar /cm/shared/apps/Picard-2.20.0/picard.jar CollectInsertSizeMetrics \
	I=$i \
      O=$i.insert_size_metrics.txt \
      H=$i.insert_size_histogram.pdf \
      M=0.5
