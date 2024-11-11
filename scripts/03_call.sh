#!/bin/bash

# Index reference genome (required for GATK)

samtools faidx /home/gabe/qmac/analyses/align/index/QlobataReference.fna

bedtools makewindows -g ~/qmac/analyses/align/index/QlobataReference.fna.fai -w 100000 > parfb.bed

while read a b c
do
  echo $a":"$b"-"$c >> parfb.txt
done < parfb.bed

while read NAME
do
  echo "markdup-$NAME".bam >> finalbamlist.txt
done < ~/qmac/inds.txt

samtools merge -@ 100 -b finalbamlist.txt merged-files.bam
samtools index -b merged-files.bam merged-files.bam.bai

parallel 'freebayes -f ~/qmac/analyses/align/index/QlobataReference.fna \
-r {} --min-mapping-quality 1 merged-files.bam > fb-{}.vcf' :::: parfb.txt

while read NAME
do
  echo "fb-$NAME".vcf >> vcflist.txt
done < parfb.txt
