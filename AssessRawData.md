# Assessment of Raw data


## MuseumI

Sequenced 8 lanes (1 species per lane) at the CGR in Liverpool. 

Data on RDSF server in Bristol: 
/projects/Butterfly_genome_analysis/butterfly_data/rawseq_Museum1_Jan2019/


Assess the sequencing based on FastQC output from CGR: 

/Users/alexjvr/2018.postdoc/NercButterflies/MuseumData/MuseumI_FastQC

### 1. Number of raw reads per lane

```
R version 3.5.0 (2018-04-23) -- "Joy in Playing"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin15.6.0 (64-bit)

library(ggplot2)

MusI.fastqc <- read.table("MuseumIFastQC_CGR", header=T)
MusI.fastqc$INDV <- gsub("/Project_ISaccheri_18649_[1-9]/", "", MusI.fastqc$Project.Sample)
MusI.fastqc$INDV <- gsub("Sample_[0-9][0-9]-", "", MusI.fastqc$INDV)
MusI.fastqc$Lane <- gsub("/Project_ISaccheri_18649_", "", MusI.fastqc$Project.Sample)
MusI.fastqc$Lane <- gsub("/Sample_[1:9].*", "", MusI.fastqc$Lane)

ggplot(MusI.fastqc, aes(x=Lane, y=raw.R1/1000000)) + geom_boxplot() + ylab("x Million reads") + xlab("Illumina Lane") + ggtitle("MuseumI: Number of raw reads (R1) per lane")
```

![alt_txt][MusIRawseq]

[MusIRawseq]:https://user-images.githubusercontent.com/12142475/53577049-dfdfbe80-3b6c-11e9-853e-eb6a246a3459.png


### 2. Does the number of reads correlate with DNA concentration in pool or with date sample was collected? 

```
DNA.ng <- read.table("MusI.DNA.ng", header=T)
MusI.tot <- merge(MusI.fastqc, DNA.ng)
head(MusI.tot)

             INDV  raw.R1  raw.R2 trim.R0 trim.R1 trim.R2  R0.
1 Cmin-20-1900-01 7031652 7031652  152297 6517743 6517743 1.15
2 Cmin-20-1900-02 7409096 7409096  168304 6488834 6488834 1.28
3 Cmin-20-1900-03 9880213 9880213  191121 9633719 9633719 0.98
4 Cmin-20-1900-04 4564900 4564900  175266 4234942 4234942 2.03
5 Cmin-20-1900-05 5304583 5304583  154609 5078740 5078740 1.50
6 Cmin-20-1900-06 4603192 4603192  155113 4279541 4279541 1.78
                                       Project.Sample Lane DNA.ng   DNA.conc
1 /Project_ISaccheri_18649_7/Sample_1-Cmin-20-1900-01    7   9.79 0.02128261
2 /Project_ISaccheri_18649_7/Sample_2-Cmin-20-1900-02    7  10.19 0.02215217
3 /Project_ISaccheri_18649_7/Sample_3-Cmin-20-1900-03    7  10.36 0.02252174
4 /Project_ISaccheri_18649_7/Sample_4-Cmin-20-1900-04    7  10.21 0.02219565
5 /Project_ISaccheri_18649_7/Sample_5-Cmin-20-1900-05    7   9.96 0.02165217
6 /Project_ISaccheri_18649_7/Sample_6-Cmin-20-1900-06    7  10.14 0.02204348
  SampleDate
1       1892
2       1898
3       1910
4       1917
5       1921
6       1921


ggplot(MusI.tot, aes(x=DNA.conc*100, y=raw.R1/1000000)) + geom_point() + ylab("Millions of reads") + xlab("% of pooled library") + ggtitle("read number is not correlated with starting DNA concentration")

ggplot(MusI.tot, aes(x=SampleDate, y=raw.R1/1000000)) + geom_point() + ylab("Millions of reads") + xlab("Sampling date") + ggtitle("read number is not correlated with age of sample") + theme(axis.text.x = element_text(angle = 270, hjust = 1))
```

![alt_txt][sampleConc.raw]

[sampleConc.raw]:https://user-images.githubusercontent.com/12142475/53579113-e5d79e80-3b70-11e9-9455-ccdfdd7e187b.png


![alt_txt][sampleAge.raw]

[sampleAge.raw]:https://user-images.githubusercontent.com/12142475/53578659-f76c7680-3b6f-11e9-95fc-129028073e7e.png


### 3. Distribution of sampling dates

```
ggplot(MusI.tot, aes(SampleDate, y=Species)) + geom_point()

```


![alt_txt][samplingDates.musi]

[samplingDates.musi]:https://user-images.githubusercontent.com/12142475/53578296-4c5bbd00-3b6f-11e9-8459-51447d87016c.png



## Choosing samples for resequencing

### May 2019

As the museum sample sequence depth limits the number of variants we can identify with our pipeline, we will resequence a subset of Museum samples for each species. This was trialled with Lymantria and increased the number of variants recovered almost ten fold. 

The choice of individuals will be based on the number of reads mapped per individual. We will target those samples that have mid-level reads as they will more readily be comparable to the well sequenced samples. The poorly sequenced samples will be ignored as we need only 18 individuals per population for confident population genetic analyses. 


