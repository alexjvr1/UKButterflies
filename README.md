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

Edit the scripts below to submit from your home directory: 

1. Set all paths to your home directory if necessary. 

2. Adjust the number of threads (PBS -t 1-xx) to equal the number of individuals to be analysed. 

3. Check that any empty arguments have been removed from the cutadapt command

4. You might have to set the path to cutadapt to find your local version




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

Check the output with samtools flagstat
```
module load apps/samtools-1.8
samtools flagstat file.bam

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


##### 3. Create a submission script for each of these regions files. 

You'll need to replace 1) the script name 2) the job name in the PBS script, 3) the regions file referred to, 4) the job prefix. 

Count the number of subset regions files: 
```
ls regions* >> regions.names
cat regions.names |wc -l
##remember to keep only the subset regions names in the regions.names file.
```



*1. Create this many copies of the submission script*
```
for i in {1..31}; do cp var_calling.20190211-161117.smsjob.sh "var_calling.20190211-161117.smsjob$i.sh";done
```
Where 31 should be replaced by number of regions files. 


*2-4 Replace the variables in each script*

```
for i in $(ls var_calling.20*sh); do while read -r a; do sed -i "s/regions/${a}/" $i; done; done < regions.names
```





*Check this has worked correctly*
```
grep "REGIONS=" var_calling*sh
grep regions var_calling*sh
grep 
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

We can pool about 500 bcffiles at a time, so we're splitting this job up into batches again and then have a final concatenation where we combine sets of 500 bcffiles. 

First move all of the files into a tmp folder
```
mkdir tmp
mv job* tmp/
```

Make a file listing all of the bcf files. And split this into several files 500 lines long each. 

```
ls tmp/*bcf >> bcflist

split -l 500 bcflist bcflist.batch
```
Rename these from bcflist.batchaa, bcflist.batchab to bcflist.batch1, bcflist.batch2, etc. 

Make sure there is a submission script for each of these batches using [03b_summary_variant_calling_batch2.sh](https://github.com/alexjvr1/UKButterflies/blob/master/03b_summary_variant_calling_batch2.sh) as a template. Remeber to check all the file paths. 

Submit these to queue. 

Check that they're all the expected size, and that all the bcf files were indexed. Once all the OUTF.b.batchxx.bcf files have been created properly, we can concatenate them all together into a final xx.raw.bcf file. 

Make a list of all the interim bcf files
```
ls OUTF.*bcf >> bcflist.ALL

```

Modify the previous script to point to the bcflist.ALL file instead of one of the batch files. 




### 03c. SNP filtering

Once we have the raw bcf file, we can look at the data and apply our filters. 

Use the [03c_variants_filtering.sh](https://github.com/alexjvr1/UKButterflies/blob/master/03c_variants_filtering.sh) script: Remember to change the variables for the species of interest. 

This script independently filters the museum and modern samples by splitting the bcf file. We have to create a list of museum and modern sample names for each of the datasets. For the expanding triplets we could also create a third population file to represent the expanding populations. For now we're keeping both modern populations together. 















### 04. Analyses

3.1. Estimate genetic diversity & population structure

3.2. Identify outlier loci 

3.3. Order outliers in genome

3.4. Identify potential genes under selection




## Multi-species comparison


