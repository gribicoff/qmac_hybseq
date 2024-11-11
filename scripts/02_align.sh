#!/bin/bash

# Align reads against reference genome
for t1 in *R1_001.fastq.gz
do
  o1=`echo ${t1%_R*}`
  t2=${t1%%R1_001.fastq.gz}"R2_001.fastq.gz"
  bwa mem -M -t 50 /home/gabe/qmac/analyses/align/index/QlobataReference.fna \
  $t1 $t2 > "aligned-$o1".sam | samtools view -b > $o1.bam
  rm "aligned-$o1".sam
done

while read NAME
do
  samtools fixmate -m -@ 104 -O bam "rg-$NAME".bam "fixed-$NAME".bam
done < ~/qmac/inds.txt

while read NAME
do
  samtools sort -m 4G -@ 104 -O bam -o "sorted-$NAME".bam "fixed-$NAME".bam
done < ~/qmac/inds.txt

while read NAME
do
  samtools markdup -@ 104 -O bam "sorted-$NAME".bam "markdup-$NAME".bam
done < ~/qmac/inds.txt

while read NAME
do
  samtools index -b -@ 104 "markdup-$NAME".bam "markdup-$NAME".bai
done < ~/qmac/inds.txt

mkdir align_qc
while read NAME
do
  samtools stats -@104 "markdup-$NAME".bam > ./align_qc/"markdup-$NAME".bam_stats
done < ~/qmac/inds.txt
multiqc ./align_qc

# Calculate alignment stats (stored in sorted_sum_*un*paired.txt);
for i in *.bam
do
  j=`echo ${i%.*}`
  echo $j >> sorted_sum_paired.txt
  samtools flagstat $i >> sorted_sum_paired.txt
done

#Calculate mean read depth per library (stored in sorted_depth_sum_*un*paired.txt)
for i in *.bam
do
  j=`echo ${i%.*}`
  echo $j >> sorted_sum_depth_paired.txt
  samtools depth $i | awk '{sum+=$3} END { print "Average = ",sum/NR}' >> sorted_sum_depth_paired.txt
done
