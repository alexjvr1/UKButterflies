# Genetic diversity of all triplets 

Plots of raw and filtered variant call files. 

All data in Dropbox file Velocity_MappingStatsPerSpecies_AJvR_20190604.xlsx

Plots drawn here: 

/Users/alexjvr/2018.postdoc/NercButterflies/G1.2.3.stats_03052019/Velocity_MappingStatsPerSpecies_Data


## 1. Individual Depth

Depth stats for raw data were gathered from raw vcf files (e.g. A1.raw.vcf) generated from variant call scripts. 

Depth stats for filtered data is based on the dataset intersecting between museum and modern samples (isec). These are located in: SpeciesName/03_variants/filtered_variant_files_date/dir/

```
VCFtools (v0.1.12b)
© Adam Auton and Anthony Marcketta 2009

Process Variant Call Format files

For a list of options, please go to:
	http://vcftools.sourceforge.net/docs.html

Questions, comments, and suggestions should be emailed to:
	vcftools-help@lists.sourceforge.net


vcftools --vcf file.vcf --depth
```


## 2. Individual missingness

Missingness stats for raw data from raw vcf files (e.g. A1.raw.vcf). These files need to be recoded to change missing variants to ".". 
```
bcftools filter -S . -O u -e 'FMT/DP=0' Species.raw.bcf |bcftools view -O b -o Species.raw.withmissing.bcf

vcftools --bcf Species.raw.withmissing.bcf --missing-indv
```


## 3. Expected heterozygosity




