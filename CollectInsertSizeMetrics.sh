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



# Check path to jar file after loading package in home directory
# eg:
# module load apps/picard-2.20.0
# echo $PATH
# /cm/shared/languages/Java-JDK-8u144/jdk1.8.0_144/bin:/cm/shared/apps/moab/7.2.9/sbin:/cm/shared/apps/moab/7.2.9/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/sbin:/usr/sbin:/opt/dell/srvadmin/bin:/cm/shared/apps/MplusDemo:.:/cm/shared/apps/torque/4.2.4.1/bin:/cm/shared/apps/torque/4.2.4.1/sbin:/cm/shared/apps/Picard-2.20.0
