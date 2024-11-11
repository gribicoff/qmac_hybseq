#!/bin/bash

grep -f regions-to-include.txt vcflist.txt > vcfs-to-include.txt

parallel 'bgzip {}' :::: vcfs-to-include.txt

while read NAME
do
  echo $NAME.gz >> vcfgzs-to-include.txt
done < vcfs-to-include.txt

parallel 'bcftools index -f {} > {}.csi' :::: vcfgzs-to-include.txt

split -400 vcfgzs-to-include.txt subset

for i in subset*
do
  echo $i >> subsetlist.txt
done

while read NAME
do
  bcftools concat --threads 100 -f $NAME -o "merged-$NAME".vcf
done < subsetlist.txt

for i in merged-subset*
do
  echo $i >> vcfpartlist.txt
done

while read NAME
do
  echo $NAME.gz >> vcfgzpartlist.txt
done < vcfpartlist.txt

parallel 'bgzip {}' :::: vcfpartlist.txt
parallel 'bcftools index -f {} > {}.csi' :::: vcfgzpartlist.txt

bcftools concat --threads 100 -f vcfgzpartlist.txt -o merged-genome.vcf | bgzip -
bcftools index -f merged-genome.vcf.gz > merged-genome.vcf.gz.csi

bcftools stats --threads 100 merged-genome.vcf.gz > merged-stats.txt

# manually copy samples (IN ORDER) from header to samples.txt

while read NAME
do
  grep $NAME abbrev-sample-species.txt >> sample-species-reheader.txt
done < samples.txt
bcftools reheader --threads 100 -s sample-species-reheader.txt merged-genome.vcf.gz > merged-genome-reheaded.vcf.gz
rm merged-genome.vcf.gz
mv merged-genome-reheaded.vcf.gz merged-genome.vcf.gz
bcftools index -f merged-genome.vcf.gz > merged-genome.vcf.gz.csi
sed 's/.*_//g' sample-species-reheader.txt | sort | uniq -c > species-count.txt
