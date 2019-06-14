# Triplet G: some initial checks

<br>


<br>

## Contig length histogram (just for contigs with variants)

*Modern data*

![G_modern2](./images/G.modern_contigsWithVariants.png)


*Museum data*


![G_museum2](./images/G.museum_contigsWithVariants.png)

**Note:** Dashed line corresponds to contig length = 50 kb and *n* is the total number of contigs that have variants after the variant filtering. There are no variants for contigs larger than 50kb in G3.

<br>

## Contig length histogram (zoom in for contigs larger than 50kb)

*Modern data*

![G_modern3](./images/G.modern_50kbUPcontigsWithVariants.png)


*Museum data*


![G_museum3](./images/G.museum_50kbUPcontigsWithVariants.png)

**Note:** *n* is the total number of contigs larger than 50kb that have variants after the variant filtering. Similar histograms for both museum and and modern samples, but as we can see below the number of snps per contig is going to be different.  
 
 <br>

## Comparison between genome assembly and population data for the number of scaffolds >50 kb 

Species | Genome Assembly | Population data (after variant filtering)
------------ | ------------- |---------
G1 | 374 | 320
G2 | 421 | 364
G3 | 135 | 0



<br>

## Plotting number of SNPs per contig length 

*Modern data*

![G_modern1](./images/G.modern_nbSNPs_contigLength.png)


*Museum data*


![G_museum1](./images/G.museum_nbSNPs_contigLength.png)

<br>

**Note:** Red dashed line corresponds to the length of the largest contig.



<br>

>
----------------------------------------------------------
__*FINAL NOTE:*__ 

Both modern and museum G3 don't show any variants in contigs larger than ~ 12 kb. 

Problem with varaint calling/filtering?





