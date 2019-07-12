# PopGen Analyses

Two datasets: 

1. Variants called with samtools mpileup and bcftools call

2. Estimates of genotype likelihood with ANGSD


### Planned analyses: 

1. Heterozygosity in sliding window

2. Expected Het museum vs modern

3. Fis museum vs modern

4. SFS comparison
 
5. Pop structure: PCA

6. Pop structure: Fst


## 5. Pop Structure: PCA

#### Based on variants

Use PCAdapt on bluecrystal

1. Filter and convert the isec vcf file

```
module load apps/bcftools-1.8
module load apps/vcftools-0.1.12b

###merge the two bcf files in species/03_variants/filteredxxx/dir/

bcftools view 0002.vcf -O b > 0002.bcf
bcftools view 0003.vcf -O b > 0003.bcf
bcftools index 0002.bcf
bcftools index 0003.bcf

bcftools merge -m id 0002.bcf 0003.bcf -O b > C3.isec.mus.mod.bcf

##check that you have the correct number of variants and that none of that none of these are now multi-allelic
##The following two files should have the same number of loci. The second file should have all the individuals (mod +mus)

vcftools --bcf 0002.bcf
vcftools --bcf C3.isec.mus.mod.bcf 

vcftools --bcf C3.isec.mus.mod.bcf --max-alleles 2

#### convert bcf to vcf

bcftools view -O v C3.isec.mus.mod.bcf > C3.isec.mus.mod.vcf

### filter for missingness

vcftools --vcf C3.isec.mus.mod.vcf --missing-indv

#This outputs out.imiss with missingness frequencies for all individuals. Use awk and sort to find all individuals with >50% missingness

awk -F "\t" '{print $1"\t"$5}' out.imiss | sort -k 2

##paste the names of these indivs into a file called indivs2remove
##Then proceed with the filtering

vcftools --vcf C3.isec.mus.mod.vcf --max-missing 0.8 --remove indivs2remove --recode --recode-INFO-all --out C3.isec.mus.mod.flt

##create a poplist with indiv name, sampling location and sampling grid square number. Call this file C3.poplist.forpcadapt
##Check that the order of this info is the same as for the filtered vcf file you'll use in pcadapt
##You can get sample order with: 

bcftools query -l C3.isec.mus.mod.vcf
```


2. Install pcadapt locally. The package was built in R.3.5


```
module load languages/R-3.5-ATLAS-gcc-7.1.0
```

The rest of this code is run in R:
```
####R script for pca drawn with PCAdapt

#install.packages("pcadapt")
library(pcadapt)

### read in popfile
poplist <- read.table("C3.poplist.forpcadapt", sep="\t", header=F)
colnames(poplist) <- c("sample", "pop", "gridsquare", "pop2")
tail(poplist)
           sample                             pop gridsquare  pop2
93 AAg-19-2016-35 Core.Bedfont Lakes Country Park   TQ078723 Mod.C
94 AAg-19-2016-36 Core.Bedfont Lakes Country Park   TQ078723 Mod.C
95 AAg-19-2016-37 Core.Bedfont Lakes Country Park   TQ078723 Mod.C
96 AAg-19-2016-38 Core.Bedfont Lakes Country Park   TQ078723 Mod.C
97 AAg-19-2016-39 Core.Bedfont Lakes Country Park   TQ078723 Mod.C
98 AAg-19-2016-40 Core.Bedfont Lakes Country Park   TQ078723 Mod.C

### read in vcf file
## For very large files you'll have to convert the vcf to plink and read in the ped file
## Check that the number of indivs and loci make sense

#C3 <- read.pcadapt("C3.isec.mus.mod.vcf", type="vcf") #The file was too big for this
C3 <- read.pcadapt("C3.forpcadapt.plink.ped", type=ped)

Summary:

	- input file:				C3.forpcadapt.plink.ped
	- output file:				/tmp/Rtmpp7eOXQ/file865c7745d8.pcadapt

	- number of individuals detected:	98
	- number of loci detected:		2614262


## Run PCA and plot variance explained by the first 20 principal components
x <- pcadapt(input=C3, K=20)
pdf("C3.pca.screeplot.pdf")
plot(x, option="screeplot")
dev.off()


## Plot pca coloured by population
pdf("C3.pca.pdf")
plot(x, option = "scores", pop = poplist$pop)
dev.off()

```

#### C1 results

![alt_txt][C1.screeplot]

[C1.screeplot]:https://user-images.githubusercontent.com/12142475/61144641-4b4b7c80-a4cd-11e9-92d6-1d700eed02fd.png

![alt_txt][C1.pca]

[C1.pca]:https://user-images.githubusercontent.com/12142475/61144657-5acac580-a4cd-11e9-8208-2c39b9caecf2.png

coordinates = museum samples. Grid reference points are the two modern populations. 
![alt_txt][C1.sampling]

[C1.sampling]:https://user-images.githubusercontent.com/12142475/61144679-72a24980-a4cd-11e9-9916-c0c60df47f07.png


![alt_txt][C2.screeplot]

[C2.screeplot]

![alt_txt][C2.pca]

[C2.pca]


![alt_txt][C3.screeplot]

[C3.screeplot]

![alt_txt][C3.pca]

[C3.pca]
