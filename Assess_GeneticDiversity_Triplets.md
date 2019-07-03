# Genetic diversity of all triplets 

Plots of raw and filtered variant call files. 

All data in Dropbox file Velocity_MappingStatsPerSpecies_AJvR_20190604.xlsx


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

