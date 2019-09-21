#!/bin/bash
#PBS -N D2_GL1_mus2.ARRAY  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=20:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-100 #array job

#Set filters
MININD="18"
MINMAF="0.05"
MINQ="20"
minMAPQ="20"


#run job in working directory
cd $PBS_O_WORKDIR 

#load modules
module load languages/gcc-6.1
angsd=~/bin/angsd/angsd

#Define variables
REGIONFILE=$(sed "${PBS_ARRAYID}q;d" regionsALL.names)


#estimate SFS for modern expanding population using ANGSD

time $angsd -b mapped.museum.names -minQ $MINQ -minMapQ $minMAPQ -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rf ${REGIONFILE} -GL 1 -minInd $MININD -out D2.GL1.mus.${REGIONFILE} -doSaf 1 -anc ../RefGenome/*fasta

#angsd -bam ~/D3_Pararge_aegeria/mapped.modern.exp.names -GL 2 -doMaf 2 -minInd $MININD -minMaf $MINMAF -minQ $MINQ -minMapQ $minMAPQ -doGlf 2 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -out RES1 -doMajorMinor 1
