# README

### Species 

[Pararge aegeria](https://butterfly-conservation.org/butterflies/speckled-wood) (Speckled wood)

[Expanding species from Family: Browns](https://docs.google.com/spreadsheets/d/1G9r50W0VV_ANZ19rIvqZpXWFemy2MW76_iXuyBuCQGA/edit?ts=5bb3681e#gid=0)

Catterpillar foodplants: grasses

![alt_txt][Parargeaegeria.fig]

[Parargeaegeria.fig]:https://user-images.githubusercontent.com/12142475/46355235-55943c80-c658-11e8-8fa3-ab38f2e0433c.png


### Genome

Full annotated genome was provided by [Chris Wheat](https://www.su.se/english/profiles/cwhea-1.191315) currently based at Stockhold University. 


##### Cumulative size of the genome

Estimated genome size ~450Mb (?)
```
/panfs/panasas01/bisc/aj18951/D3_Pararge_aegeria/03_variants
/panfs/panasas01/bisc/aj18951/D3_Pararge_aegeria/RefGenome

##check the number of contigs in the Reference genome

cat Pararge_aegeria_v2.fa.fai  |wc -l
93568

##Check that all these contigs have been used during variant calling
grep "#contig" variants.raw.vcf |wc -l
26567

##Read the contig lengths into R and plot the lengths as a proportion of the total estimated genome length. This is done on my personal mac so dataframes were copied to local computer. 

/Users/alexjvr/2018.postdoc/NercButterflies/D3_Pararge_aegeria/RawDataStats/

R version 3.5.0 (2018-04-23) -- "Joy in Playing"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin15.6.0 (64-bit)

library(dplyr)
library(ggplot2)

genomesize <- 450000000
contig.length <- read.table("Pararge_aegeria_v2.fa.fai", header=F)
contig.length$prop.TotGenome <- (contig.length$V2/genomesize)*100
colnames(contig.length) <- c("scaffold", "length", "bytesindex.fasta", "length.fasta", "bytes", "prop.GenomeTot")
contig.length <- contig.length[order(-contig.length$prop.GenomeTot),]
head(contig.length)
contig.length$index <- (1:nrow(contig.length))
contig.length2 <- contig.length %>% mutate(cumsum=cumsum(prop.GenomeTot))

pdf("D3_Pararge_aegeria_CumulativeGenome_scaffolds.pdf")
ggplot(contig.length2, aes(x=index, y=cumsum)) + geom_line() + ggtitle("Pararge aegeria Reference Genome (predicted 450Mb)") + ylab("Cumulative fraction of genome (%)") + xlab("Number of scaffolds")
dev.off()


#And number of contigs >1Mb in size
pdf("D3_Pararge_aegeria_1MbScaffolds.pdf")
ggplot(contig.length2[which(contig.length$length>1000000),], aes(length/1000000)) + geom_histogram() + ggtitle("Pararge aegeria: frequency and size of contigs longer than 1Mb") + ylab("count")  +xlab("Contig size (Mbp)")
dev.off()

```


![alt_txt][D3_Cumulativegenomesize]

[D3_Cumulativegenomesize]:https://user-images.githubusercontent.com/12142475/53562560-903bcb80-3b49-11e9-9391-463d1d9f28ee.png


![alt_txt][D3_1Mbcontigs]

[D3_1Mbcontigs]:https://user-images.githubusercontent.com/12142475/53562581-9cc02400-3b49-11e9-9be7-67b590a73848.png


### Raw data

RDSF server:
/projects/Butterfly_genome_analysis/butterfly_data/rawseq/Pararge_aegeria

