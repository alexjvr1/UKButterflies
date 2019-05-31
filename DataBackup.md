# Where all the data are stored at UoB 

Data are backed up in long term storage on the RDSF in two folders: 

/projects/Butterfly_genome_analysis

/projects/Adaptation_in_rainforest_flies


## Data locations

#### */projects/Butterfly_genome_analysis/alex*

A1_Lymantria_monacha

B1_Cupido_minimus

C1_Aricia_artaxerxes

C3_Aricia_agestis (Velocity data distinct from JB & MdJ's Pop gen data - see below)

D3_ararge_aegeria_data

AriciaAgestsis (population genomics and linkage map data - JB & MdJ)


#### */projects/Butterfly_genome_analysis/butterfly_data*

rawseqModern2

rawseqMuseum1_Jan2019

rawseq_Pararge_aegeria

##### */projects/Butterfly_genome_analysis/butterfly_data/RefGenomes*

Hesperia_comma

Hipparchia_semele

Maniola_jurtina

Ochlodes_sylvanus

Plebejus_argus

Thymalicus_acteon


##### */projects/Adaptation_in_rainforest_flies/velocity_project/*

rawseqModern1_2018

rawseqMuseum2_2019

A1_Lymantria_monacha

G1_Thymelicus_acteon

G2_Ochlodes_sylvanus



##### */projects/Adaptation_in_rainforest_flies/velocity_project/Genome_fastqs/*

Raw data for each species

Aphantopus_hyperantus  

Eilema_griseola

Ochlodes_sylvanus

Xanthorhoe_fluctuata

Aricia_agestis

Hesperia_comma

Polyommatus_coridon

Aricia_artaxerxes

Maniola_jurtina

temp.sh

Cupido_minimus

Miltochrista_miniata

Thymelicus_acteon





See below for a description of folders. 

## File organisation

#### 1. Raw data

All raw output from the sequencing facility will be kept as is in a folder and are named according to the library batch. 

e.g. Modern1, Museum1. 

#### 2. Reference genomes

Raw data for reference genomes have been kept as received - raw data is organised by species. 

Reference genomes for each species is stored in their respective species folders. 

#### 3. 20 species

Data for each species will be kept in a species folder ordered by triplet name (see shared google doc). 

Data within each folder will be organised as follows


#### *00_raw_reads_modern* and *00_raw_reads_museum*

Raw fastq files for modern and museum samples. 

These folders will be emptied after the cutadapt reads have been produced to save space. Raw data can be recovered from the raw data folders. See shared doc to identify relevant libraries. 


#### *01a_modern_cutadapt_reads* and *01a_museum_cutadapt_reads*

Adapter trimmed fastq samples for modern and musuem. 

#### *02a_modern_mapped* and *02a_museum_mapped* and *mapped*

Mapped and indexed (.bam and .bai) files for modern and musuem. Where samples have been further processed, all these data have been moved into a single folder called mapped leaving the 02a.. folders empty. 

#### *03_variants*

Folder where samtools mpileup and bcftools call variant calling takes place. 

The */tmp* folder contains the raw bcf files for each region. 

The *xaa*.. or *regionsaa* files list the independently processed regions. 

The intermediate and final raw combined bcf files are found directly in the 03_variants folder

The *filtered_variant_files_xxx* folder contains the filtered museum and modern bcf files

The */filtered_variant_files_xxx/dir/* folder contains the modern, museum, and combined bcf files that contain only the intersecting datasets. 

#### *RefGenome*

Reference genome, index, and gff file for each species. 


#### *ANGSD*

All ANGSD analyses. 

