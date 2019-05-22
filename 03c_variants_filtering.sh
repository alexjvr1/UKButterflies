#!/bin/bash
##Specify resources
#PBS -N 03c_G3.VariantFiltering  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=3:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)


#load modules
module load apps/bcftools-1.8

#run job in working directory
cd $PBS_O_WORKDIR


#software and raw bcf file
BCFTOOLS='bcftools';
BCFRAW="/newhome/aj18951/G3_Hesperia_comma/03_variants/intermed.bcf/OUTF.concat.bcf";
OUTDIR="/newhome/aj18951/G3_Hesperia_comma/03_variants/filtered_variant_files_March2019";
NTOT="104";
NMOD="";
NMUS=""

#NTOT= total number of samples

mkdir -p $OUTDIR >& /dev/null;

### We first apply overall filters

# remove any other variants than SNPs
# remove private variants (here just different SNPs than reference)
# remove multiallelic SNPs
# remove any SNPs for with we do not have data for at least 50% of individuals
# remove SNPs with too less coverage (loose already taken care)
# remove SNPs with too much coverage - here more than 10 time the expected number of reads for a coverage of 5X
# remove SNPs with MAF < 0.01
# remove variant with quality call below 20

NSAMPLES=$NTOT; #here the combined numbers for both contemporary and museum populations
EXPCOVIND=5;
MAXDPMULT=10; # multiplier for upper limit on coverage

COV='0.5';
MINDP=10;
MAXDP=$(($NSAMPLES * $EXPCOVIND * $MAXDPMULT));
MAF='0.01';
QS=20;


OUTALL=$OUTDIR/variants.bial.noindel.qs"$QS".cov"$(perl -e 'print ('$COV'*100);')".mdp"$MINDP"Mdp"$MAXDP".maf"$MAF".bcf;
echo $OUTALL

EXCLUDE='MAF<'$MAF' | QUAL<'$QS' | TYPE!="snp" | N_ALT>1 | (AN/2)<('$NSAMPLES'*'$COV') | SUM(FMT/DP)<'$MINDP' | SUM(FMT/DP)>'$MAXDP' | AC/AN=1'

$BCFTOOLS filter -S . -O u -e 'FMT/DP=0' $BCFRAW | \
$BCFTOOLS filter -e "$EXCLUDE" -O u | \
$BCFTOOLS view -O b -o $OUTALL

$BCFTOOLS index $OUTALL;

# Then we filter every population independantly

BCFIN=$OUTALL;

### Contemporary

# remove any SNPs for with we do not have data for at least 40% of individuals
# remove SNPs with too less coverage (loose already taken care)
# remove SNPs with too much coverage - here more than 10 time the expected number of reads for a coverage of 5X
# remove variant with quality call below 20

IDLIST="/newhome/aj18951/G3_Hesperia_comma/03_variants/modern_samples_list.dsv";

NSAMPLES=$NMOD; #here the number of modern samples
EXPCOVIND=5;
MAXDPMULT=10; # multiplier for coverage uplimit
 
COV='0.4';
MINDP=10;
MAXDP=$(($NSAMPLES * $EXPCOVIND * $MAXDPMULT));

OUTCONTEMP=$OUTDIR/modern_variants.bial.noindel.qs"$QS".cov"$(perl -e 'print ('$COV'*100);')".mdp"$MINDP"Mdp"$MAXDP".bcf;
echo $OUTCONTEMP;

 EXCLUDE='(AN/2)<('$NSAMPLES'*'$COV') | SUM(FMT/DP)<'$MINDP' | SUM(FMT/DP)>'$MAXDP'';

$BCFTOOLS filter -S . -O u -e 'FMT/DP=0' $BCFIN | \
$BCFTOOLS view -S $IDLIST -O u | \
$BCFTOOLS filter -e "$EXCLUDE" -O u | \
$BCFTOOLS view -O b -o $OUTCONTEMP

$BCFTOOLS index $OUTCONTEMP;

### Museum

# remove any SNPs for with we do not have data for at least 40% of individuals
# remove SNPs with too less coverage (loose already taken care)
# remove SNPs with too much coverage - here more than 10 time the expected number of reads for a coverage of 5X
# remove variant with quality call bellow 20
 
IDLIST="/newhome/aj18951/G3_Hesperia_comma/03_variants/museum_samples_list.dsv";

NSAMPLES=$NMUS; #here the combined numbers for both contemporary and museum populations
EXPCOVIND=5; 
MAXDPMULT=10; # multiplier for coverage uplimit
  
COV='0.4';
MINDP=10;
MAXDP=$(($NSAMPLES * $EXPCOVIND * $MAXDPMULT));
 
OUTMUSEUM=$OUTDIR/museum_variants.bial.noindel.qs"$QS".cov"$(perl -e 'print ('$COV'*100);')".mdp"$MINDP"Mdp"$MAXDP".bcf;
echo $OUTMUSEUM;
 
EXCLUDE='(AN/2)<('$NSAMPLES'*'$COV') | SUM(FMT/DP)<'$MINDP' | SUM(FMT/DP)>'$MAXDP'';
  
$BCFTOOLS filter -S . -O u -e 'FMT/DP=0' $BCFIN | \
$BCFTOOLS view -S $IDLIST -O u | \
$BCFTOOLS filter -e "$EXCLUDE" -O u | \
$BCFTOOLS view -O b -o $OUTMUSEUM

$BCFTOOLS index $OUTMUSEUM;

echo ""
echo "Done!"
