#!/bin/bash

mkdir ~/qmac/analyses/ldprune_qc
cd ~/qmac/analyses/ldprune_qc

printf 'it_filt\nprune_it_filt_50snps\nprune_it_filt_50kb\nprune_it_filt_500kb\nprune_it_filt_longestgene\n' > filterlist.txt

# generate PLINK bed file from filtered VCF, calculate MAF stats

cp ~/qmac/analyses/it_filt/geno80pc.vcf filtered_snps_unpruned.vcf
plink --vcf ~/qmac/analyses/it_filt/filtered_snps_unpruned.vcf --double-id --allow-extra-chr --set-missing-var-ids "@:#" --make-bed --out it_filt
plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --freq --out it_filt
bgzip filtered_snps_unpruned.vcf
tabix filtered_snps_unpruned.vcf.gz

cat sra_inds.txt | while read LINE
do
  printf "$(grep $LINE coords_noSHF.csv)\n"
done | cut -f2-3 -d','

# prune using 50 SNP sliding window, calculate MAF stats

plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --indep-pairwise 50 10 0.1 --out prune_it_filt_50snps
plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --extract prune_it_filt_50snps.prune.in --make-bed --out prune_it_filt_50snps
plink --bfile prune_it_filt_50snps --double-id --allow-extra-chr --set-missing-var-ids @:"#" --freq --out prune_it_filt_50snps

# prune using 50 kb sliding window, calculate MAF stats

plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --indep-pairwise 50 kb 10 0.1 --out prune_it_filt_50kb
plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --extract prune_it_filt_50kb.prune.in --make-bed --out prune_it_filt_50kb
plink --bfile prune_it_filt_50kb --double-id --allow-extra-chr --set-missing-var-ids @:"#" --freq --out prune_it_filt_50kb

# prune using 0.5 kb sliding window, calculate MAF stats

plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --indep-pairwise 1 kb 10 0.1 --out prune_it_filt_1kb
plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --extract prune_it_filt_1kb.prune.in --make-bed --out prune_it_filt_1kb
plink --bfile prune_it_filt_05kb --double-id --allow-extra-chr --set-missing-var-ids @:"#" --freq --out prune_it_filt_1kb

# thin to 1 SNP per 500 bp

plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --bp-space 500 --make-bed --out thin_it_filt_500bp

# prune using 500kb sliding window, calculate MAF stats

plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --indep-pairwise 500 kb 10 0.1 --out prune_it_filt_500kb
plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --extract prune_it_filt_500kb.prune.in --make-bed --out prune_it_filt_500kb
plink --bfile prune_it_filt_500kb --double-id --allow-extra-chr --set-missing-var-ids @:"#" --freq --out prune_it_filt_500kb

# prune using length of longest gene (42,601 bp), calculate MAF stats

plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --indep-pairwise 42.6 kb 10 0.1 --out prune_it_filt_longestgene
plink --bfile it_filt --double-id --allow-extra-chr --set-missing-var-ids @:"#" --extract prune_it_filt_longestgene.prune.in --make-bed --out prune_it_filt_longestgene
plink --bfile prune_it_filt_longestgene --double-id --allow-extra-chr --set-missing-var-ids @:"#" --freq --out prune_it_filt_longestgene

# get summary across filtered datasets

for NAME in it_filt prune_it_filt_50snps prune_it_filt_50kb prune_it_filt_500kb prune_it_filt_longestgene
do
    printf "$NAME,$(wc -l $NAME.bim | cut -f1 -d' '),$(awk '{sum+=$5} END {print sum/NR}' $NAME.frq)\n"
done | sed '1i dataset,num_SNPs,mean_MAF' > prunestats.txt

## ADMIXTURE analysis on filtered SNP datasets

# prep thinned ADMIXTURE analysis

sed -i 's/N._//g; s/\.1//g' thin_it_filt_500bp.bim

mkdir ~/qmac/analyses/ldprune_qc/thin_admix
cd ~/qmac/analyses/ldprune_qc/thin_admix

parallel 'for j in {1..50}; do admixture -s time ../thin_it_filt_500bp.bed {} --cv -j10 > thin_it_filt_500bp_K{}_run${j}.out; mv thin_it_filt_500bp.{}.Q thin_it_filt_500bp_K{}_run${j}.Q; mv thin_it_filt_500bp.{}.P thin_it_filt_500bp_K{}_run${j}.P; done' ::: {2..8}

for j in {2..8}
do
  for i in {1..50}
  do
    printf "thin_it_filt_500bp_K${j}_run${i}.out\t$(grep '^Loglikelihood' thin_it_filt_500bp_K${j}_run${i}.out | sed 's/.*\s//g')\n"
  done | awk 'NR==1{max=$2;line=$1} $2>max{max=$2;line=$1} END{print line}'
done | sed 's/out$/Q/g' > thin_it_filt_500bp_bestadmixrun_perK.txt


cat thin_it_filt_500bp_bestadmixrun_perK.txt | while read LINE
do
  printf "$(grep 'CV' ${LINE/Q/out} | sed 's/).*//g;s/.*=//g')\t$(grep 'CV' ${LINE/Q/out} | sed 's/.*: //g')\n"
done | sed '1i K\tcv_err' > thin_it_filt_500bp_cverror.txt

sed 's/QUE002099_bic/QUE002099_mue/g' ../onesnppergene/admix/indslist.txt > indslist.txt

conda activate R
Rscript ~/qmac/scripts/cverror_plot.R "thin_it_filt_500bp_cverror"
Rscript ~/qmac/scripts/plot_admix_acrossK.R "indslist.txt" "thin_it_filt_500bp_bestadmixrun_perK.txt" "~/qmac/analyses/new_analyses/figures/" --vanilla
conda deactivate

# prep PLINK files for ADMIXTURE - edit chromosome names in bim files

while read NAME
do
  sed -i 's/N._//g; s/\.1//g' $NAME.bim
done < filterlist.txt

mkdir ~/qmac/analyses/ldprune_qc/admix
cd ~/qmac/analyses/ldprune_qc/admix

# run admixture with CV, 50 runs per value of K

parallel 'for i in {2..8}; do for j in {1..50}; do admixture -s time ../{}.bed $i --cv -j10 > {}_K${i}_run${j}.out; mv {}.${i}.Q {}_K${i}_run${j}.Q; mv {}.${i}.P {}_K${i}_run${j}.P; done; done' :::: ../filterlist.txt

# find best run per K

while read NAME
do
  for j in {2..8}
  do
    for i in {1..50}
    do
      printf "${NAME}_K${j}_run${i}.out\t$(grep '^Loglikelihood' ${NAME}_K${j}_run${i}.out | sed 's/.*\s//g')\n"
    done | awk 'NR==1{max=$2;line=$1} $2>max{max=$2;line=$1} END{print line}'
  done | sed 's/out$/Q/g' > ${NAME}_bestadmixrun_perK.txt
done < ../filterlist.txt

# find CV error per best run per value of K

while read NAME
do
  cat ${NAME}_bestadmixrun_perK.txt | while read LINE
  do
    grep 'CV' ${LINE/Q/out}
  done | sed '1i K\tcv_err' > ${NAME}_cverror.txt
done < ../filterlist.txt

# graph CV error as a function of K

conda activate R
while read NAME
do
  Rscript ~/qmac/scripts/cverror_plot.R "${NAME}_cverror"
done < ../filterlist.txt
conda deactivate

## one SNP per gene tests

mkdir ~/qmac/analyses/ldprune_qc/onesnppergene
cd ~/qmac/analyses/ldprune_qc/onesnppergene

# copy file with gene positions

# get genes with chr positions for all SNPs in pruned dataset ('--indep-pairwise 50 10 0.1')

grep -f <(cut -f1 ../prune_it_filt_50snps.bim | sort -u) <(cut -f1-4 ~/qmac/analyses/blast/sorted-Qlobatagenes-PCG.bed) | sed 's/N._//g; s/\.1//g' > allgenes_chrbp.txt

# get set file with lists of SNPs within each gene region

plink --bfile ../prune_it_filt_50snps --double-id --allow-extra-chr --set-missing-var-ids "@:#" --make-set allgenes_chrbp.txt --write-set --out allgenes

# get number of SNPs per gene

grep '^T' allgenes.set | while read LINE
do
  printf "$(sed "0,/$LINE/d;/END/Q" allgenes.set | wc -l)\t$LINE\n"
done | sort -nr > numsnpspergene.txt

# randomly subsample one SNP per gene (10 subsets)

for i in {1..10}
do
  grep '^T' allgenes.set | while read LINE
  do
    sed "0,/$LINE/d;/END/Q" allgenes.set | sort -R | head -1
  done > onesnppergene_sample_${i}_variantlist.txt
done

grep '^T' ../../ldprune_qc/onesnppergene/allgenes.set | while read LINE
do
  sed "0,/$LINE/d;/END/Q" ../../ldprune_qc/onesnppergene/allgenes.set | wc -l
done

ls onesnppergene*variantlist.txt | sort -n -t "_" -k3 | cut -f1-3 -d "_" > samplelist.txt

# make new PLINK bed files for each subset

while read NAME
do
  plink --bfile ../prune_it_filt_50snps --double-id --allow-extra-chr --set-missing-var-ids "@:#" --extract "${NAME}_variantlist.txt" --make-bed --out "$NAME"
done < samplelist.txt

# run ADMIXTURE from K = 2 to K = 8 with 50 replicates per value of K

mkdir admix
cd admix

parallel 'for i in {2..8}; do for j in {1..50}; do admixture -s time ../{}.bed $i --cv -j4 > {}_K${i}_run${j}.out; mv {}.${i}.Q {}_K${i}_run${j}.Q; mv {}.${i}.P {}_K${i}_run${j}.P; done; done' :::: ../samplelist.txt

# find best run per K

while read NAME
do
  for j in {2..8}
  do
    for i in {1..50}
    do
      printf "${NAME}_K${j}_run${i}.out\t$(grep '^Loglikelihood' ${NAME}_K${j}_run${i}.out | sed 's/.*\s//g')\n"
    done | awk 'NR==1{max=$2;line=$1} $2>max{max=$2;line=$1} END{print line}'
  done | sed 's/out$/Q/g' > ${NAME}_bestadmixrun_perK.txt
done < ../samplelist.txt

# find CV error per best run per value of K

while read NAME
do
  cat ${NAME}_bestadmixrun_perK.txt | while read LINE
  do
    grep 'CV' ${LINE/Q/out} | sed 's/^.*=//g; s/).* /\t/g'
  done | sed '1i K\tcv_err' > ${NAME}_cverror.txt
done < ../samplelist.txt

# graph CV error as a function of K

conda activate R
while read NAME
do
  Rscript ~/qmac/scripts/cverror_plot.R "${NAME}_cverror"
done < ../samplelist.txt
conda deactivate

# create masterlist of all CV errors across subsets

cat ../samplelist.txt | while read NAME
do
  sed '1d' ${NAME}_cverror.txt | sed "s/^/${NAME##*_}\t/g"
done | sed '1i subset\tK\tcv_err' > allsubsets_cverror.txt

# graph CV error as a function of K across subsets

conda activate R
Rscript ~/qmac/scripts/cverror_stack.R allsubsets_cverror
conda deactivate

# generate admixture plots

cut -f1 -d' ' ~/qmac/analyses/ldprune_qc/prune_it_filt_50snps.fam > indslist.txt

conda activate R
while read NAME
do
  Rscript ~/qmac/scripts/plot_admix_acrossK.R "indslist.txt" "${NAME}_bestadmixrun_perK.txt" "~/qmac/analyses/ldprune_qc/onesnppergene/" --vanilla
done < ../samplelist.txt
conda deactivate

conda activate R
while read NAME
do
  Rscript ~/qmac/scripts/plot_admix_acrossK.R "indslist.txt" "${NAME}_bestadmixrun_perK.txt" "~/qmac/analyses/ldprune_qc/onesnppergene/" --vanilla
done < ../samplelist.txt
conda deactivate