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
## For very large files you'll have to convert the vcf to plink and read in the ped file. See above for keeping chromosome names. 
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

### A1 results (declining species)

Lymantria monacha screeplot 

![alt_txt][LM.scree]

[LM.scree]:https://user-images.githubusercontent.com/12142475/61230838-ca85be00-a722-11e9-95a4-f425785eca5c.png



![alt_txt][LM.pca]

[LM.pca]:https://user-images.githubusercontent.com/12142475/61230837-c9ed2780-a722-11e9-930f-3169022b1e2a.png

Map: Green flags= Museum, Red = Modern

![alt_txt][map]

[map]:https://user-images.githubusercontent.com/12142475/61228054-6a404d80-a71d-11e9-8f0b-a535fa2d08f3.png



### C1 results (declining species)

![alt_txt][C1.screeplot]

[C1.screeplot]:https://user-images.githubusercontent.com/12142475/61144641-4b4b7c80-a4cd-11e9-92d6-1d700eed02fd.png

![alt_txt][C1.pca]

[C1.pca]:https://user-images.githubusercontent.com/12142475/61144657-5acac580-a4cd-11e9-8208-2c39b9caecf2.png

coordinates = museum samples. Grid reference points are the two modern populations. 
![alt_txt][C1.sampling]

[C1.sampling]:https://user-images.githubusercontent.com/12142475/61144679-72a24980-a4cd-11e9-9916-c0c60df47f07.png


### C2 results (stable) 

![alt_txt][C2.screeplot]

[C2.screeplot]:https://user-images.githubusercontent.com/12142475/61146705-90be7880-a4d2-11e9-9fa7-db2c611c990d.png

![alt_txt][C2.pca]

[C2.pca]:https://user-images.githubusercontent.com/12142475/61146733-a0d65800-a4d2-11e9-9ad0-c218235e8490.png


### C3 results (expanding)

![alt_txt][C3.screeplot]

[C3.screeplot]:https://user-images.githubusercontent.com/12142475/61146787-c2374400-a4d2-11e9-97cc-d9c82157aa02.png

![alt_txt][C3.pca]

[C3.pca]:https://user-images.githubusercontent.com/12142475/61146809-ce230600-a4d2-11e9-86d8-077f32def424.png


### D1 results (contracting)

![alt txt][D1.screeplot]

[D1.screeplot]:https://user-images.githubusercontent.com/52965134/61298722-a1bf0080-a7d6-11e9-936a-54f775dc42ec.png

![alt txt][D1.pca]

[D1.pca]:https://user-images.githubusercontent.com/52965134/61298855-ddf26100-a7d6-11e9-9584-e03b06f165b7.png

Map: Green= Museum, Red= Beaulieu Heath, Blue= Telegraph Hill, Pink= Turf Hill

Baulieu Heathm Telegraph Hill, Turf Hill are modern populations. 

![alt txt][D1.sampling]

[D1.sampling]:https://user-images.githubusercontent.com/52965134/61299182-856f9380-a7d7-11e9-91ee-9266c5583fad.png


### D2 results (stable)

![alt txt][D2.screeplot]

[D1.screeplot]:https://user-images.githubusercontent.com/52965134/61299545-32e2a700-a7d8-11e9-90a6-9d640a95f8d5.png

![alt txt][D2.pca]

[D2.pca]:https://user-images.githubusercontent.com/52965134/61299461-121a5180-a7d8-11e9-8661-11dc304a9f6f.png

Map: Green= Museum, Red= Hemsbury, Blue= Stoke Woods

Hemsbury and Stoke Woods are modern populations.

![alt txt][D2.sampling]

[D2.sampling]:https://user-images.githubusercontent.com/52965134/61299758-9bca1f00-a7d8-11e9-9074-0a881288c910.png


### D3 results (expanding)

![alt txt][D3.screeplot]

[D3.screeplot]:https://user-images.githubusercontent.com/52965134/61299931-f2375d80-a7d8-11e9-8b92-9eaf1875bb15.png

![alt txt][D3.pca]

[D3.pca]:https://user-images.githubusercontent.com/52965134/61300059-33c80880-a7d9-11e9-97eb-514cc4d7715f.png

Map: Turquoise= Museum, Red= Carverel copse, Brown= Hadleigh Railway Walk, Green= Martin Down NNR, Blue= Snapes Wood, Purple= Stour Wood, Pink= Wolves Wood

Carverel copse, Hadleigh Railway Walk, Martin Down NNR, Snapes Wood, Stour Wood and Wolves Wood are modern populations. 

![alt txt][D3.sampling1]

[D3.sampling1]:https://user-images.githubusercontent.com/52965134/61300370-bd77d600-a7d9-11e9-9638-f6ba855998c6.png

![alt txt][D3.sampling2]

[D3.sampling2]:https://user-images.githubusercontent.com/52965134/61300471-eac48400-a7d9-11e9-8322-5a8d201426a4.png

![alt txt][D3.sampling3]

[D3.sampling3]:https://user-images.githubusercontent.com/52965134/61300522-0596f880-a7da-11e9-9242-00e344b2990e.png


### G1 results (declining)

![alt txt][G1.screeplot]

[G1.screeplot]:https://user-images.githubusercontent.com/52965134/61300896-c2895500-a7da-11e9-86cd-104766b0f60d.png

![alt txt][G1.pca]

[G1.pca]:https://user-images.githubusercontent.com/52965134/61301140-3b88ac80-a7db-11e9-9226-e32c7e4e4f8b.png

Map: Blue= Museum, Red= Durlston Head Country Park, Green= Lulworth Cove

![alt txt][G1.sampling]

[G1.sampling]:https://user-images.githubusercontent.com/52965134/61301337-98846280-a7db-11e9-86a0-4da2282ae8fd.png


### G2 results (Stable)

![alt txt][G2.screeplot]

[G2.screeplot]:https://user-images.githubusercontent.com/52965134/61302441-adfa8c00-a7dd-11e9-8e26-b05bb96e64ee.png

![alt txt][G2.pca]

[G2.pca]:https://user-images.githubusercontent.com/52965134/61303568-a76d1400-a7df-11e9-9db6-badebdbe2c5c.png

Map: Blue= Museum, Red=Havering Country Park, Green= Kemsing Access Land

Havering Country Park and Kemsing Access Land are modern Populations

![alt txt][G2.sampling]

[G2.sampling]:https://user-images.githubusercontent.com/52965134/61303745-ef8c3680-a7df-11e9-8e0c-342493a2e9d7.png


### G3 results (Expanding)

![alt txt][G3.screeplot]

[G3.screeplot]:https://user-images.githubusercontent.com/52965134/61303868-24988900-a7e0-11e9-9ea7-59a9462ac3d8.png

![alt txt][G3.pca]

[G3.pca]:https://user-images.githubusercontent.com/52965134/61304137-94a70f00-a7e0-11e9-8929-24230e598d37.png

Map: Purple= Museum, Red=Brockham Quarry, Brown= Buckland Hills, Green= Chantry Hill, Turquiose- Cissbury Ring, Blue= Denbies Hillside, Pink= Newtimber Hill

Brockham Quarry, Buckland Hills, Chantry Hill, Cissbury Ring, Denbies Hillside and Newtimber Hill are all modern populations.

![alt txt][G3.sampling]

[G3.sampling]:https://user-images.githubusercontent.com/52965134/61305083-19465d00-a7e2-11e9-9cef-b8f88c6cc170.png

![alt txt][G3.sampling1]

[G3.sampling1]:https://user-images.githubusercontent.com/52965134/61305176-45fa7480-a7e2-11e9-9383-e0623d9d121c.png
