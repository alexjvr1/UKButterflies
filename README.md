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

With cutadapt we're removing all sequences that are shorter than 20bp and 3' quality trimmed to remove bases with PHRED quality score of < 20. 


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



### 2. Raw data to SNPs

2.1. Demultiplex 






### 3. Analyses

3.1. Estimate genetic diversity & population structure

3.2. Identify outlier loci 

3.3. Order outliers in genome

3.4. Identify potential genes under selection




## Multi-species comparison


