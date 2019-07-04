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

Figure: 
```
ggplot(data, aes(x=Pop2, y=Missingness.isec)) + geom_boxplot(aes(colour=Pop))+ggtitle("Missingness in isec dataset")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![alt_txt][missingness.flt]

[missingness.flt]:https://user-images.githubusercontent.com/12142475/60621036-97fbcd00-9dd4-11e9-9c0e-009d1e327260.png


## 3. Expected heterozygosity

Calculated using vcftools 
```
vcftools --bcf isec.bcf --het

```

Museum samples consistently have lower heterozygosity. I need to check why this is. I'll compare results with ANGSD outputs. 


![alt_txt][ExpHet]

[ExpHet]:https://user-images.githubusercontent.com/12142475/60665418-4e56c500-9e5c-11e9-9a7f-92a57cac8a0c.png


## 4. Variants vs Depth

Is raw depth a predictor of the number of variants in the final dataset? 

```
ggplot(data, aes(x=Depth.Raw, y=Nr.variants.isec, group=Species)) + geom_point(aes(shape=Pop, colour=Species))+ggtitle("Raw depth vs number of variants")

ggplot(data[which(data$Pop=="Museum"&data$Depth.Raw<3),], aes(x=Depth.Raw, y=Nr.variants.isec, group=Species)) + geom_point(aes(colour=Species))+ggtitle("Museum: Raw depth vs number of variants")

ggplot(data, aes(x=Depth.Raw, y=Nr.variants.before.isec, group=Species)) + geom_point(aes(shape=Pop, colour=Species))+ggtitle("Raw depth vs number of variants (non intersecting)")
```

![alt_txt][variants.vs.depth]

[variants.vs.depth]:https://user-images.githubusercontent.com/12142475/60619574-e313e100-9dd0-11e9-9479-0baf80c281f0.png


![alt_txt][museum.depth.vs.variants]

[museum.depth.vs.variants]:https://user-images.githubusercontent.com/12142475/60689702-29903b00-9eb8-11e9-9a75-8608e4f35f13.png


![alt_txt][non.isec.variants.vs.DP]

[non.isec.variants.vs.DP]:https://user-images.githubusercontent.com/12142475/60663910-73e1cf80-9e58-11e9-8ede-45f854b27024.png
