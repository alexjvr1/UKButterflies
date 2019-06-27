#!/bin/bash
#PBS -N BWA_mod1  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-52

#run job in working directory
cd $PBS_O_WORKDIR 

#Load modules
module load apps/bwa-0.7.15
module load apps/samtools-1.8

#Define variables

RefSeq=Lymantria_monacha_v1.0.fasta
total_files=`find 01a_modern_cutadapt_reads/ -name '*.fastq.gz' | wc -l`
#It is more efficient to run this hashed code in local directory before submitting to queue
#ls 01a_modern_cutadapt_reads/*R1*fastq.gz >> R1.modern.names
#sed -i s:01a_modern_cutadapt_reads/::g R1.modern.names
#ls 01a_modern_cutadapt_reads/*R2*fastq.gz >> R2.modern.names
#sed -i s:01a_modern_cutadapt_reads/::g R2.modern.names
#mkdir 02a_modern_mapped


NAME1=$(sed "${PBS_ARRAYID}q;d" R1.modern.names)
NAME2=$(sed "${PBS_ARRAYID}q;d" R2.modern.names)

echo "mapping started" >> map.log
echo "---------------" >> map.log

##Check if Ref Genome is indexed by bwa
if [[ ! RefGenome/$RefSeq.fai ]]
then 
	echo $RefSeq" not indexed. Indexing now"
	bwa index RefGenome/$RefSeq
else
	echo $RefSeq" indexed"
fi


##Map with BWA MEM and output sorted bam file

sample_name=`echo ${NAME1} | awk -F "_190115" '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time bwa mem RefGenome/$RefSeq 01a_modern_cutadapt_reads/${NAME1} 01a_modern_cutadapt_reads/${NAME2} | samtools sort -o 02a_modern_mapped/${NAME1}.bam" >> map.log
time bwa mem RefGenome/$RefSeq 01a_modern_cutadapt_reads/${NAME1} 01a_modern_cutadapt_reads/${NAME2} | samtools sort -o 02a_modern_mapped/${NAME1}.bam
