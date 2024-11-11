#!/bin/bash

# preliminary handling, filtering: break up by chromosome, filter records

parallel 'bcftools view -O v -r {} --threads 7 ~/qmac/analyses/calling/merged-genome.vcf.gz > {}-trunc.vcf' :::: included_regions.txt
parallel 'vcfallelicprimitives -kg {}-trunc.vcf > dec-{}.vcf' :::: included_regions.txt
parallel 'vt normalize -r ~/qmac/analyses/align/index/QlobataReference.fna dec-{}.vcf -o norm-{}.vcf' :::: included_regions.txt
parallel 'vcftools --vcf norm-{}.vcf --remove-indels --min-alleles 2 --max-alleles 2 --max-missing .2 --minQ 30 --recode --recode-INFO-all --stdout > prefiltered-{}.vcf' :::: included_regions.txt
parallel 'bcftools sort -O z -o sorted-prefiltered-{}.vcf.gz prefiltered-{}.vcf' :::: included_regions.txt
parallel 'bcftools index -f sorted-prefiltered-{}.vcf.gz > sorted-prefiltered-{}.vcf.gz.csi' :::: included_regions.txt

while read NAME
do
  echo "sorted-prefiltered-$NAME.vcf.gz" >> sortedchrlist.txt
done < included_regions.txt

bcftools concat --threads 100 -f sortedchrlist.txt -O z -o sorted-prefiltered-genome.vcf.gz

# depth, mac filters

vcftools --gzvcf sorted-prefiltered-genome.vcf.gz --mac 2 --minDP 5 --min-meanDP 8 --max-meanDP 60 --maxDP 60 --recode --recode-INFO-all --stdout | gzip -c > dpmacfilt.vcf.gz
# depth, filter
vcftools --gzvcf dpmacfilt.vcf.gz --max-missing .5 --recode --recode-INFO-all --stdout | gzip -c > geno50pc.vcf.gz
vcftools --gzvcf geno50pc.vcf.gz --missing-indv --out geno50pc
awk 'NR > 1 && $5 > 0.9 {print $1}' geno50pc.imiss > geno50pc.Remindv
vcftools --gzvcf geno50pc.vcf.gz --remove geno50pc.Remindv --recode --recode-INFO-all --stdout | gzip -c > ind90pc.vcf.gz
vcftools --gzvcf ind90pc.vcf.gz --max-missing .6 --recode --recode-INFO-all --stdout | gzip -c > geno60pc.vcf.gz

vcftools --gzvcf geno60pc.vcf.gz --missing-indv --out geno60pc
awk 'NR > 1 && $5 > 0.7 {print $1}' geno60pc.imiss > geno60pc.Remindv
vcftools --gzvcf geno60pc.vcf.gz --remove geno60pc.Remindv --recode --recode-INFO-all --stdout | gzip -c > ind70pc.vcf.gz
vcftools --gzvcf ind70pc.vcf.gz --max-missing .7 --recode --recode-INFO-all --stdout | gzip -c > geno70pc.vcf.gz
vcftools --gzvcf geno70pc.vcf.gz --missing-indv --out geno70pc
awk 'NR > 1 && $5 > 0.5 {print $1}' geno70pc.imiss > geno70pc.Remindv
vcftools --gzvcf geno70pc.vcf.gz --remove geno70pc.Remindv --recode --recode-INFO-all --stdout | gzip -c > ind50pc.vcf.gz

zcat ind50pc.vcf.gz | vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" | vcffilter -f "SAF > 1 & SAR > 1 & RPR > 1 & RPL > 1" | bcftools view --threads 10 -e 'QUAL/FMT/DP<20 & FMT/QA/FMT/AO<30' -O z -o itfilt_infofields.vcf.gz -

vcftools --gzvcf itfilt_infofields.vcf.gz --missing-indv --out itfilt_infofields
awk 'NR > 1 && $5 > 0.25 {print $1}' itfilt_infofields.imiss > itfilt_infofields.Remindv
vcftools --gzvcf itfilt_infofields.vcf.gz --remove itfilt_infofields.Remindv --recode --recode-INFO-all --stdout | gzip -c > ind25pc.vcf.gz
vcftools --gzvcf ind25pc.vcf.gz --max-missing .8 --recode --recode-INFO-all --stdout > geno80pc.vcf

tabix geno80pc.vcf.gz
bcftools stats -r NW_022155175.1 geno80pc.vcf.gz

# change QUE002099 from bic to ??? (WIP)

bcftools query -l geno80pc.vcf.gz > vcf_inds.txt
sed -i 's/QUE002099_bic/QUE002099_???/g' vcf_inds.txt
bcftools reheader -s vcf_inds.txt
