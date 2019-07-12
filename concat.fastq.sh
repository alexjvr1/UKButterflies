#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 12/07/2019 05:49 
###########################################

## Concatenates fastq reads from a list of file names
## R1 and R2 are processed separately


#PBS -N G2.OS.concat.cutadapt  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-66

#run job in working directory
cd $PBS_O_WORKDIR

#Define variables

#create files with sample names listed for the 33 samples
#ls ../01a_museum2_cutadapt_reads/*R1*gz > samplenames.museum2.R1
#ls ../01a_museum2_cutadapt_reads/*R2*gz > samplenames.museum2.R2
#ls ../01a_museum_cutadapt_reads/*R1*gz >> samplenames.museum1.R1
#ls ../01a_museum_cutadapt_reads/*R2*gz >> samplenames.museum1.R2
#remove all the extra sample names from the last two files. Make sure the sample names are in the same order in all files.
#remove path before sample names

NAME.MUS1.R1=$(sed "${PBS_ARRAYID}q;d" samplenames.museum1.R1)
NAME.MUS1.R2=$(sed "${PBS_ARRAYID}q;d" samplenames.museum1.R2)
NAME.MUS2.R1=$(sed "${PBS_ARRAYID}q;d" samplenames.museum1.R1)
NAME.MUS2.R2=$(sed "${PBS_ARRAYID}q;d" samplenames.museum1.R2)


##Concat R1 fastq files

sample_name.mus1.R1=`echo ${NAME.MUS1.R1} | awk -F "_190115" '{print $1}'`
sample_name.mus2.R1=`echo ${NAME.MUS2.R1} | awk -F "_190627" '{print $1}'`
echo "[concatenating] $sample_name.mus1.R1 and $sample_name.mus2.R1"
printf "\n"
echo "time cat ../01a_museum_cutadapt_reads/${NAME.MUS1.R1} ../01a_museum2_cutadapt_reads/${NAME.MUS2.R1} > ${NAME.MUS2.R1}.concat.fastq.gz" >> concat.mus.R2.log
time cat ../01a_museum_cutadapt_reads/${NAME.MUS1.R2} ../01a_museum2_cutadapt_reads/${NAME.MUS2.R2} > ${NAME.MUS2.R2}.concat.fastq.gz
