# Temporal scans:

The present document describe the methodology used for generating data and analysis population data in
order to generate genome scans.

Most of the exploratory method is done with Lymantria monacha (black arches) at the moment. Other
species will be analysed soon hereafter.


## I. Samples and library preparation and sequencing:

Carl is still to hand me over a complete version of the protocol he used. But briefly, we used a tweaked
version of NEBnext protocol with custom barcode adaptors. The treatment of modern and museum
samples differed in two points.

Modern samples were usually extracted with Quiagen DNA easy kit and the DNA was sheared and 3
'PCR' cycles were realized to incorporate barcodes. It is in fact, the lowest number of cycles one can do to
generate working libraries with this protocol.
The museum samples were extracted using Quiagen Micro kit, DNA was NOT sheared, and we
performed 10 PCR cycles here.
Lymantria monacha's populations were sequenced in CGR at Liverpool on a Hiseq 4000 with 150 pe
cycles.


## II. Data quality, possible adjustment for museum sample sequencing:

Lymantria monacha modern samples look extremely good. No changes to operate there. Coverage may be
around 4X which is good and unexpected... It does not look like we have massive adaptor contaminations.

Lymantria monacha museum samples look less good. Still some analysis to be performed here, but
quality of data degrade greatly after 100bp. Adaptor contamination is significant roof with ~80% of
sequences contaminated. Sequences will be quite short here.


## III. Trimming:

Trimming method changed from previous work here. Indeed, trimmomatic was not handling the massive
adaptor contamination well for museum samples. I therefore used cutadapt to do it. I developed a function
in Sheffield to do it in parallel on iceberg. I intend to adapt it on Bristol cluster. But cutadapt does not
have the nice features of trimmomatic when it comes to trimming and quality filtering. I therefore used
cutadapt first to remove adaptor contamination and then trimmomatic for further filtering


### Cutadapt:

I trimmed 6 adaptors (3 forward, 3 reverse):

```
forward
AGATCGGAAGAGCACACGTCTGAACTCCAGTC
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
reverse:

AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

The first ones (forward and reverse) are the sequence of the NEBnext adaptors, but I noticed we had a lot
of polyA and polyT in the sequences, so decided to trim them as well.


Here is an example command I used:

```
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC -a
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -a
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -A
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -A
CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -A
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA -m 21 -q 20 -o
./01tris_museum_cutadapt_trial/LM-01-1900-01_R1.fastq.gz -p ./01tris_museum_cutadapt_trial/LM-01-1900-01_R2.fastq.gz 00_raw_reads_museum/LM-01-1900-01_R1.fastq.gz 00_raw_reads_museum/LM-
01-1900-01_R2.fastq.gz
```


Cutadapt was used on paired files and an option was added to discard bases with phred score < 20 starting
from the 5' end and read pairs with a minimum length bellow 20 bp.
This provided us read pairs files that I used a input from trimmomatic.


### Trimmomatic:

Trimmomatic was run to 'finish' the trimming job and filter further for quality:

```
ILLUMINACLIP:/pub46/romain/Highlight/useful_files/adaptor_files/NEBNext_adaptors.fa:2:30:8:1:true
with ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip
threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
LEADING: 20, remove bases from start of the reads under this phredscore
TRALING: 20, remove bases from end of the reads under this phredscore
SLIDINGWINDOW:4:20, do sliding windows of 4 bp from end of reads and removes bases if under 20
phredscore threshold
MINLEN:20 discard any reads bellow 20 bp in length
AVGQUAL:20 discard any reads with average quality bellow 20
```

So here there is adaptor trimming (but I am not expecting to catch anything there), then base quality
trimming, length check and mean quality check.
I used trimmomatic version 0.36.


## IV. Reads mapping:

I will inspire myself greatly from Victor's rfseq pipeline.
