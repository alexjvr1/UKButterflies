# Output from MapDamage2

Cross-species comparison of MapDamage outputs, used to estimate Cytosine deamination in the museum samples. 


**Note

MapDamage only assessess properly paired, non-overlapping, and inward-facing reads. The samtools stats function outputs the number of inward facing pairs and a lot of other useful statistics for a bam file. 

In our case, MapDamage assessed only ~5% of the reads in the dataset. This corresponds to the number of inward facing read pairs from the samtools stats output.
For the museum data we have orientation information for only ~10% of the properly paired reads compared to ~50% for the modern data. I expect that all the overlapping read pairs are excluded from the analysis. I will need to incorporate these in the pipeline after the York meeting (Sept 2019). 


Once MapDamage has run all the outputs are written to results_samplename/*

We're interested in the misincorporations.txt files which count all the mutations relative to the reference sequence. 

1. Combine all the files together

2. Read into R

3. Sort by 

