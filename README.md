# UKButterflies
NERC UK Butterfly project


## Pipeline for Velocity project from raw data to filtered variants: 

The scripts for the final pipeline are outlined below. The scripts are organised in 3 folders within a main folder. 

Use this pipeline by running scripts in the pipeline/ folder in the numbered order. These scripts generate submission scripts for 
each step that can be submitted to the BlueCrystal p3 queue (i.e. qsub script.sh). 

Inputs for each step should be submitted via the command line. 

    
    |
    -----> pipeline  
    
            This contains all the scripts that generate submission scripts for BlueCrystal. Options can be specified in the command line. 
    
    |
    -----> wrapper
        
            Scripts called by pipeline scripts. Wrapper scripts for generating Queue request and specifying inputs from command line for                the tools to be called.     
            
    |
    -----> tools
            
            Scripts of tools or functions used in each step. These are called by the wrapper script
   
            
            

### 00. Quality control

We're using [fastQC](https://wiki.gacrc.uga.edu/wiki/FastQC) and [multiQC](https://multiqc.info) for an initial assessment of data quality. 

fastQC is installed on BlueCrystal, but multiQC needs to be installed locally: 

```
#run in home directory on bluecrystal

module load languages/python-2.7.10

pip install multiqc
```

1. Navigate to the folder with all the raw data. 

2. Make sure the raw R1 and R2 files are all in one folder (or set the paths in the script/names file). 

3. create an R1.names and R2.names file listing all the forward and reverse reads. 

4. edit the fastQC script [00_fastQC_ARRAY.sh](https://github.com/alexjvr1/UKButterflies/blob/master/00_fastQC_ARRAY.sh): change the job name and the number of array jobs (to the number of indivs in the folder). 

5. Run the fastQC script (qsub 00_fastQC_ARRAY.sh). You can check on the job progress by looking at the output files or using "qstat -u username -t"

6. Once this is complete use multiQC to concatenate all the output files: 

```
multiQC .
```

This will create a combined html to compare outputs across all samples. It also creates a folder (multiqc_data) with all the parsable info in txt files. 




### 01. Trimming adapter and poor quality sequence



We're removing all sequences that are shorter than 20bp and 3' quality trimmed to remove bases with PHRED quality score of < 20 with Cutadapt.  

If you're running this on BlueCrystal, you'll have to install cutadapt locally first
```
module load languages/python-anaconda3-5.2.0

#install cutadapt in your home directory using the web instructions
pip3 install --user --upgrade cutadapt

#Check that this cutadapt works
~/.local/bin/cutadapt --help


##Check if this directory is in your PATH:
echo $PATH

##And add to PATH if it isn't yet
PATH="$PATH:~/.local/bin/"

##Now you can run cutadapt directly
cutadapt --help
```

#### *pipeline*

01a_museum_cutadapt_filtering_trimming.sh

01a_modern_cutadapt_filtering_trimming.sh

01b_museum_trimmomatic_filtering_trimming.sh

01b_modern_trimmomatic_filtering_trimming.sh


#### *wrapper*


01a_parallel_cutadapt_bluecp3.sh


Edit the generated script above to submit from your home directory: 

```
1. Set all paths to your home directory if necessary. 

2. Adjust the number of threads (PBS -t 1-xx) to equal the number of individuals to be analysed. 

3. Check that any empty arguments have been removed from the cutadapt command

4. You might have to set the path to cutadapt to find your local version
```

To incorporate new data (e.g. resequencing of some individuals to increase mean depth), new fastq files need to be adapter trimmed. Fastq files are concatenated after this using the script [concat.fastq.sh](https://github.com/alexjvr1/UKButterflies/blob/master/concat.fastq.sh)

Reseq data are kept in the following folders: 
```
00_raw_data_museum2
01a_museum2_cutadapt_reads
01a_mus.concat_cutadapt_reads  ## concatednated museum1 and museum2
02a_museum2_mapped  ##see below
```




### 02. Map to reference genome with BWA mem

It is more efficient to run this code in local directory before submitting to queue
```
#Index the reference genome if needed. Check whether the *fasta.fai* file exists in the SpeciesName/RefGenome/ folder in your local directory. If not, run the indexing code. 

#index reference genome
module load apps/bwa-0.7.15
bwa index RefGenome/*fasta


#Create files with input names
ls 01a_museum_cutadapt_reads/*R1*fastq.gz >> R1.museum.names
sed -i s:01a_museum_cutadapt_reads/::g R1.museum.names

ls 01a_museum_cutadapt_reads/*R2*fastq.gz >> R2.museum.names
sed -i s:01a_museum_cutadapt_reads/::g R2.museum.names


#make output directories
mkdir 02a_museum_mapped
mkdir 02a_modern_mapped

#Check that you're pointing to the correct reference genome

#Check that the file separator makes sense. This will be different for the musuem and modern samples because the samples are named differently. On the line: 
##sample_name=`echo ${NAME1} | awk -F "_L007" '{print $1}'`
#Change the -F "xxx" according to the file names. 
#e.g the above works for files named as follows: 
#HS-01-2016-26_L007_cutadapt_filtered_R2.fastq.gz
#we want only the first part of this name to carry through. 
```


#### *pipeline*

02_MapwithBWAmem.ARRAY_museum.sh

02_MapwithBWAmem.ARRAY_modern.sh

* Check that everything has mapped correctly by checking the file sizes. If the mapping is cut short (e.g. by exceeding the requested walltime) the partial bam file will look complete and can be indexed. But the bam file size will be small (~500kb) and empty when you look at it. 


```
#To determine file size

du -sh *bam   

#To see bam file
module load apps/bcftools-1.8
bcftools view file.bam | head
```

Check the output with samtools flagstat
```
module load apps/samtools-1.8
samtools flagstat file.bam

#make a flagstat log file for all of the samples
for i in $(ls *bam); do ls $i >>flagstat.log && samtools flagstat $i >> flagstat.log; done
```

Index the bam files with the script [index.bamfiles.sh](https://github.com/alexjvr1/UKButterflies/blob/master/index.bamfiles.sh)


###### Reseq Data: Mapping

The concatenated mus1 and mus2 fastq files are mapped with bwa mem and stored in the following folder: 
```
02a_museum2_mapped 
```

Continue on with variant calling as below and with genotype likelihood estimation as described in the ANGSD pipeline. 

## Variant calling

After mapping the pipeline will split into 1) variant calling for final vcf with genotypes called for each individual, and 2) using likelihoods associated with variants for downstream analyses. 


### 03. Variant calling with samtools and bcftools

We use samtools mpileup to stack all reads for a locus from all individuals in the dataset. The likelihood of a variant at a particular site is determined based on the dataset-wide depth of two alleles. Next an individual genotype likelihood is calculated based on the depth of coverage within that individual. The final output is a vcf file with locus likelihood and depths for each individual's genotype. 

*pipeline*

The variant calling step benefits greatly from the multi-threading and array capabilites on BlueCrystal. We break up the job by calling variants for each genomic region separately. This script will generate a single submission script for the entire dataset which will refer to a newly generated file called "regions" which contains the names of all the genomic regions (contig-xx, scaffold-xx). We have to split this script up because the maximum number of array jobs that we can run for a single script on BlueCrystal is 100.


##### 1. First generate the submission script and the regions file. Remember to change the variables in this script for the species you're working on: 

[03a_variant_calling_bluecp3.sh](https://github.com/alexjvr1/UKButterflies/blob/master/03a_variant_calling_bluecp3.sh)

This generates a regions file where each "region" contains a max of 50 contigs/scaffolds and is max 4634570bp long. 

The script needs to be modified depending on what files we're calling variants on. If calling independently for museum or modern, change the input file (-i) to point to the appropriate folder (e.g. 02a_modern_mapped)

When calling on museum and modern together, create a folder called mapped and move all the mapped reads and indexes to this folder. Change the input file to mapped: 
```
mkdir mapped
mv 02a_modern_mapped/* mapped/
mv 02a_museum_mapped/* mapped/
```

Run the script above. 


##### 2. Split the regions file into several files with 100 regions in each

```
split -l 100 regions regions

```

This will generate several files called regionsaa, regionsab.. etc. The number will depend on how fragmented the reference genome is. 


##### 3. Create a submission script for each of these regions files. 

You'll need to replace 1) the script name 2) the job name in the PBS script, 3) the regions file referred to, 4) the job prefix. 

Count the number of subset regions files: 
```
ls regions* >> regions.names

#remove the top line "regions", as we only want to count the subset files
#you can use nano for this.

#then count the number of subset files
cat regions.names |wc -l

##remember to keep only the subset regions names in the regions.names file.
```



*1. Create this many copies of the submission script*
```
for i in {1..31}; do cp var_calling.20190211-161117.smsjob.sh "var_calling.20190211-161117.smsjob$i.sh";done
```
Where 31 should be replaced by number of regions files. 


*2-4 Replace the variables in each script*

This script doesn't work.. 
```
for i in $(ls var_calling.20*21.*sh); do while read -r a; do sed -i "s/regions/${a}/" $i; done; done < regions.names
```


*5 Correct number of threads requested*

Check how many threads need to be run for the total script by looking at one of the scripts
```
grep "PBS -t" var_calling..sh
```
This should return something like
```
var_calling.20190603-155815.smsjob10.sh:#PBS -t 1-1800
```
The total number of threads is 1800, but we can only run 100 at a time. Hence the multiple submission scripts.. 

Modify PBS -t for all of the subset scripts. 
```
for i in $(ls var_calling*sh); do sed -i 's:1-xxx:1-100:g' $i; done

```


*Check this has worked correctly*
```
grep "REGIONS=" var_calling*sh
grep regions var_calling*sh
grep "PBS -N"
grep "PBS -t"
```



##### 4. Submit all the scripts

Due to computational limits you can submit only ~6 jobs at a time. 

```
qsub var_calling.smsjob1.sh

```

Check on the status of the run
```
qstat -u bluecrystal.username

```
Keep an eye on the queue and keep submitting more jobs until all the bcf files have been created. 


##### 5. Concatenate all the bcf files

We can pool about 200 bcffiles at a time, so we're splitting this job up into batches again and then have a final concatenation where we combine sets of 200 bcffiles. 

First move all of the files into a tmp folder
```
mkdir tmp
mv job* tmp/
```

Make a file listing all of the bcf files. And split this into several files 500 lines long each. 

```
ls tmp/*bcf >> bcflist

split -l 200 bcflist bcflist.batch
```
Rename these from bcflist.batchaa, bcflist.batchab to bcflist.batch1, bcflist.batch2, etc. 

Make sure there is a submission script for each of these batches using [03b_summary_variant_calling_batch2.sh](https://github.com/alexjvr1/UKButterflies/blob/master/03b_summary_variant_calling_batch2.sh) as a template. Remeber to check all the file paths. 

Submit these to queue. 

Check that they're all the expected size, and that all the bcf files were indexed. Once all the OUTF.b.batchxx.bcf files have been created properly, we can concatenate them all together into a final xx.raw.bcf file.

#### NB
```
I've had a problem with these scripts timing out with no reported error. e.g. for G3 all longer scaffolds (>15k) were excluded from the concatenated bcfs. 
Hence I've decreased the bcflist.batch length from 500 to 200 lines. 

I used batches of 500 for triplet G and D. 
```





Make a list of all the interim bcf files
```
ls OUTF.*bcf >> bcflist.ALL

```

Modify the previous script to point to the bcflist.ALL file instead of one of the batch files. 


*##### Problems Encountered*

A1
```
empty bcf files were not deleted, thus the bcftools concat script terminated. 
Empty bcfs (i.e. where no variants are called) should be deleted by the script and have been for all other species. Perhaps the run timed out? But I received no error. 

Soln: 
Delete all empty bcf files before running concat

Find empty files:
find . -size 0
```

C2
```
Several 03_variants/tmp/*report.txt files with: 

"0 variants called for 0 samples/individuals"

These were all from job20 and job21. I will rerun these two jobs and replace the file

This seems like the variant calling didn't complete. 

In addition, concatenation of bcflist.batch10 failed with the following error: 
Failed to parse header: 03_variants/tmp/job200084.variants.raw.bcf

As I'm rerunning job20 and 21, I'll see if this error is resolved. 
```

C3
```
I checked whether I get the same size bcf file from two independent concatenations (from raw bcfs in tmp file to final C3.raw.bcf). 

This worked. So I am satisfied that we're not running into trouble during the concatenation steps. 
```

General
```
Nodes seem to get stuck on BlueCrystal. If a job looks like it's taking much longer than anticipated this could be a problem. This was particularly problematic during variant calling and could be the cause of some of the errors encountered above. 
This also happened when concatenating bcf files. 
```

Scratch space
```
The standard scratch space is 1.8Gb per node. This gets filled up when running the full dataset (which includes the reseq data). There's alternative scratch space available. 
Add 

export TMPDIR=/local

to .bashrc file in home directory
```


### 03c. SNP filtering

#### Part1: Filtering by population

Once we have the raw bcf file, we can look at the data and apply our filters. 

This script independently filters the museum and modern samples by splitting the bcf file. We have to create a list of museum and modern sample names for each of the datasets. For the expanding triplets we could also create a third population file to represent the expanding populations. For now we're keeping both modern populations together. 

Use the [03c_variants_filtering.sh](https://github.com/alexjvr1/UKButterflies/blob/master/03c_variants_filtering.sh) script: Remember to change the variables for the species of interest: 

1. Point to the raw bcf file

2. Create the 03_variants/museum_samples_list.dsv and 03_variants/modern_samples_list.dsv files

```
module load apps/bcftools-1.8
bcftools query -l xx.raw.bcf >> xx.samples_list.dsv
cat xx.samples_list.dsv
```

copy and paste the museum sample names into museum_samples_list.dsv using nano. Same with modern samples. 

3. Change the sample numbers in the file in three places: total, modern, museum

4. Change the output file according to the month

e.g 

OUTDIR="/newhome/bzzjrb/D2_Maniola_jurtina/03_variants/filtered_variant_files_June2019"


Once this has run, determine the number of variants discoverd for each of the populations separately: 
```
module load apps/bcftools-1.8
module load apps/vcftools-0.1.12b

bcftools view -Ov modern..bcf >> modern...vcf
bcftools view -Ov museum..bcf >> museum...vcf

#Number of variants in modern only file after filtering
vcftools --vcf modern..vcf

#Variants in museum only file after filtering
vcftools --vcf museum.vcf
```



#### Part2: Intersect the two files

Find the intersection between these datasets. I.e. Keep only variants present in both files. 

```
cd filtered_variant_files_xxx
module load apps/bcftools-1.8
bcftools isec -p dir museum_variants.bial.noindel.xxx.bcf modern_variants.bial.noindel.xxx.bcf 
```

This command creates a four vcf files and a README.txt file in dir/
```
#The README file tells you what each of the vcf files is
cd dir
cat README.txt

module load apps/vcftools-0.1.12b
vcftools --vcf 0002.vcf
vcftools --vcf 0003.vcf
```


To see how many variants overlap with the annotated genome: 

Modify the gff file if necessary. Bedtools doesn't like 0-based coordinates. Sam's code to turn 0 start coords into 1
```
awk -F $'\t' 'BEGIN {OFS = FS};{ if($4==0) {print $1,$2,$3,1,$5,$6,$7,$8,$9} else {print}}' ../RefGenome/RefGenome.gff > ../RefGenome/RefGenome.mod.gff 
```

And mask: 
```
module load apps/bedtools2  
bedtools subtract -header -a input.vcf -b ../RefGenome/RefGenome.mod.gff > input.masked.vcf
```




### 04. ANGSD pipeline

a. MapDamage: re-calibrate bam files after correcting for Cytosine deamination anticipated in museum samples

b. Downsample modern data: To create comparative datasets we downsample the modern data to have ~the same number of sequences (i.e. same depth) as the museum data) 

c. ANGSD - call variants: pipeline for estimating variant likelihoods across all sites. 


#### 04a. MapDamage


[MapDamage2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694634/) is a package used to estimate and correct for Cytosine deamination (or any other transition/transversion bias in the data). This is a problem anticipated for ancient DNA, and possibly for museum data. 

This needs to be locally installed on BlueCrystal. Follow the instructions in the [tutorial](https://ginolhac.github.io/mapDamage/)

1. Create a folder in the species directory called mapped.reseq.mapdamage. 

2. Move all museum bam files to this folder. Use the concatenated files with reseq data included where available. Sample numbers should correspond to the total museum samples in the LibraryPrep spreadsheet. 

3. Create a file listing all the bamfiles

```
ls *bam > bamfiles.mus.names

```

4. Copy the script [04a_mapDamage_museum.sh](https://github.com/alexjvr1/UKButterflies/blob/master/04a_mapDamage_museum.sh) to the mapped.reseq.mapdamage folder. Change the job name, the number of threads, and check the path to the reference genome. 

5. Submit to queue. 

6. [Analyse output stats](https://github.com/alexjvr1/UKButterflies/blob/master/04a_MapDamageOutputs.md)


#### 04b. Downsample modern data

Due to the difference in sample quality between museum and modern samples, mean coverage is much higher for the modern data. This may bias the confidence in variant calls downstream. To avoid this problem I will downsample the modern data to the same mean depth as the museum data. 

First filter the bam files to include only reads with PHRED quality >20 and properly paired reads using the [04b1_Filter_museum_bam_pp.PHRED20.md](https://github.com/alexjvr1/UKButterflies/blob/master/04b1_Filter_museum_bam_pp.PHRED20.md) script. 


Use samtools flagstat to calculate the number of properly paired reads in the recalibrated and filtered museum files. 

```
module load apps/samtools-1.8
for i in $(ls results*/*flt.bam); do ls $i >> mus.flagstat.log && samtools flagstat $i >> mus.flagstat.log; done
```

Do the same for the modern samples. 

Enter these data in the "Rescaled.ProperlyPaired.Q20" column in the Velocity_MapingStatsPerSpecies_AJvR_20190604.xlsx sheet on Dropbox. Calculate the mean number of museum reads and the proportion of modern reads to downsample to. 

Use the [04b_Downsample.sh](https://github.com/alexjvr1/UKButterflies/blob/master/04b_Downsample.sh) script to downsample the modern bam files. Remember to change the job name and the PROP variables and create the input file listing all the modern bams. 



#### 04c. Estimate Genotype Likelihoods

When all the samples have been pre-processed we can estimate genotype likelihoods using ANGSD. 

First we'll estimate the site frequency spectrum for each population separately. In this first step we'll also filter sequences for quality. 
This step is run as an array by splitting the genome up into 100 regions. 

Next we'll intersect the datasets from museum and modern populations to include only loci found in both datasets. 

Then we'll re-estimate the SFS using only these sites. 

Finally we can calculate Tajima's D and Fst for these sites. 




### 05. Analyses

3.1. Estimate genetic diversity & population structure

3.2. Identify outlier loci 

3.3. Order outliers in genome

3.4. Identify potential genes under selection




## Multi-species comparison


