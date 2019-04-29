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
   
            
            

### 0. Quality control

*pipeline*

00_fastqc_raw_museum.sh

00_fastqc_raw_modern.sh

*wrapper*

parallel_fastqc_bcp3.sh


### 1. Reference genome


### 2. Raw data to SNPs

2.1. Demultiplex 


### 3. Analyses

3.1. Estimate genetic diversity & population structure

3.2. Identify outlier loci 

3.3. Order outliers in genome

3.4. Identify potential genes under selection


## Multi-species comparison


