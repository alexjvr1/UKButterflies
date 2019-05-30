
#1.GREP CONTIG INFO FROM VCF FILE
#NOTE: Change G if using for another triplet

#unzip all vcf files
gunzip -d G*.vcf.gz

#Get the contig ids and length from vcf files
for i in {1..3}
do
  grep '##contig' G$i.*.vcf| cut -f3 -d=|cut -f1 -d, > temp.id
  grep '##contig' G$i.*.vcf|cut -f4 -d=|cut -f1 -d'>' > temp.length
  paste temp.id temp.length > G$i.id_length.txt
done

rm temp*
#2.CALCULATE number of snps per contig and per species  
#(USING BCFTOOLS: more efficient than my previous bash script. It runs in a few minutes): 

for i in {1..3}
do
tabix -p vcf G$i.modern_variants.bial.noindel.qs20.cov40.mdp10Mdp2200.vcf.gz #needs indexing
#extract contigs with snps
bcftools index -s G$i.modern_variants.bial.noindel.qs20.cov40.mdp10Mdp2200.vcf.gz > G$i.seqinfo_variants

#compare files:create a file for contigs without variants
awk 'FNR==NR{a[$1];next};!($1 in a)' G$i.seqinfo_variants G$i.id_length.txt > G$i.seqinfo_novariants
done

 
