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

   #### *pipeline*

00_fastqc_raw_museum.sh

00_fastqc_raw_modern.sh

   #### *wrapper*

parallel_fastqc_bcp3.sh



### 01. Trimming adapter and poor quality sequence

#### *pipeline*

We're removing all sequences that are shorter than 20bp and 3' quality trimmed to remove bases with PHRED quality score of < 20 with Cutadapt.  


01a_museum_cutadapt_filtering_trimming.sh

01a_modern_cutadapt_filtering_trimming.sh

01b_museum_trimmomatic_filtering_trimming.sh

01b_modern_trimmomatic_filtering_trimming.sh


#### *wrapper*


01a_parallel_cutadapt_bluecp3.sh



### 02. Map to reference genome with BWA mem

```
#It is more efficient to run this code in local directory before submitting to queue

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


## Variant calling

After mapping the pipeline will split into 1) variant calling for final vcf with genotypes called for each individual, and 2) using likelihoods associated with variants for downstream analyses. 


### 03. Variant calling with samtools and bcftools

We use samtools mpileup to stack all reads for a locus from all individuals in the dataset. The likelihood of a variant at a particular site is determined based on the dataset-wide depth of two alleles. Next an individual genotype likelihood is calculated based on the depth of coverage within that individual. The final output is a vcf file with locus likelihood and depths for each individual's genotype. 

*pipeline*

The variant calling step benefits greatly from the multi-threading and array capabilites on BlueCrystal. We break up the job by calling variants for each genomic region separately. This script will generate a single submission script for the entire dataset which will refer to a newly generated file called "regions" which contains the names of all the genomic regions (contig-xx, scaffold-xx). We have to split this script up because the maximum number of array jobs that we can run for a single script on BlueCrystal is 100.

##### 1. First generate the submission script and the regions file. Remember to change the variables in this script for the species you're working on: 

03a_variant_calling_bluecp3.sh


##### 2. Split the regions file into several files with 100 regions in each

```
split -l 100 regions regions

```

This will generate several files called regionsaa, regionsab.. etc. The number will depend on how fragmented the reference genome is. 


###### 3. Create a submission script for each of these regions files. 

You'll need to replace 1) the job name in the PBS script, 2) the regions file referred to, 3) the job prefix, 4) script name. 








### 2. Raw data to SNPs

2.1. Demultiplex 






### 3. Analyses

3.1. Estimate genetic diversity & population structure

3.2. Identify outlier loci 

3.3. Order outliers in genome

3.4. Identify potential genes under selection




## Multi-species comparison


