#!/bin/bash

## PRE-TRIM QC

# run fastqc on raw fastq files to assess read quantity/quality before trimming

cd ~/qmac/data
mkdir ~/qmac/data/fastqc
mkdir ~/qmac/data/fastqc/fastqc_multi
fastqc *fastq.gz -o ./fastqc

# summarize with multiqc

multiqc ./fastqc/*.zip -o ./fastqc/fastqc_multi

## TRIMMING WITH TRIMMOMATIC

# make list of fastq read files

if [ ! -d "~/qmac/analyses" ]; then
  mkdir ~/qmac/analyses
fi
mkdir ~/qmac/analyses/trimming
cd ~/qmac/analyses/trimming

ls ~/qmac/data/*R1_001.fastq.gz | sed 's/_R1.*//g' > readfiles.txt

# get Illumina adapter sequences for trimming reads

wget https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa

# trim reads

while read NAME
do
  trimmomatic PE -threads 30 "${NAME}_R1_001.fastq.gz"  "${NAME}_R2_001.fastq.gz"  "trimmed_paired_${NAME}_R1_001.fastq.gz" "trimmed_unpaired-${NAME}_R1_001.fastq.gz" "trimmed-paired-${NAME}_R2_001.fastq.gz" "trimmed-unpaired-${NAME}_R2_001.fastq.gz" ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40
done < readfiles.txt

## POST-TRIM QC

# run fastqc to assess read quantity/quality per read file after trimming

mkdir postqc
mkdir postqc/postqc_multi
for i in paired/trimmed*
do
  fastqc $i -o ./postqc
done

# summarize with multiqc

multiqc ./postqc/*.zip -o ./postqc_multi

