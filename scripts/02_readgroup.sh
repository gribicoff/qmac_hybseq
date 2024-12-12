#!/bin/bash

#Pull flowcell ID/lane info from read header of each R1 fastq file

while read NAME
do
  zcat "trimmed-paired-$NAME"_R1_001.fastq.gz | head -1 >> prelimrg.txt
done < ~/qmac/inds.txt

#Get SM read group field

while read NAME
do
  i=`echo ${NAME%_S*}`
  echo $i >> sm.txt
done < ~/qmac/inds.txt

#Get ID read group field

while read NAME
do
  i=`echo ${NAME% *}`
  j=`echo ${i#*:}`
  k=`echo ${j#*:}`
  l=`echo ${k%:*}`
  m=`echo ${l%:*}`
  n=`echo ${m%:*}`
  o=`echo ${n/:/.}`
  echo $o >> id.txt
done < prelimrg.txt

#Get PU read group field

paste -d "." id.txt sm.txt > pu.txt

#Get LB read group field

while read NAME
do
  echo "$NAME"_01 >> lb.txt
done < ~/qmac/inds.txt

#Assign readgroups
while IFS="," read a b c d e
do
  samtools addreplacerg -@ 100 -r ID:$b -r LB:$c -r PL:ILLUMINA \
  -r PU:$d -r SM:$e -o "rg-$a".bam "$a".bam
done < rg.csv

samtools addreplacerg -@ 80 -r ID:000000000-BY462.1 -r LB:JCB-21_S42_L001_01 -r PL:ILLUMINA \
-r PU:000000000-BY462.1.QUE002637 -r SM:QUE002637 -o rg-JCB-21_S42_L001.bam JCB-21_S42_L001.bam

samtools addreplacerg -@ 80 -r ID:000000000-BPKMG.1 -r LB:WT_6980_EGX_P1_S31_L001_01 -r PL:ILLUMINA \
-r PU:000000000-BPKMG.1.QUE002937 -r SM:QUE002937 -o rg-WT_6980_EGX_P1_S31_L001.bam WT_6980_EGX_P1_S31_L001.bam
WT_6980_EGX_P1_S31_L001,000000000-BPKMG.1,WT_6980_EGX_P1_S31_L001_01,000000000-BPKMG.1.QUE002937,QUE002937
