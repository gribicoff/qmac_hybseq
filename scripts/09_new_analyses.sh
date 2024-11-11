#!/bin/bash

# copy files to new

mkdir ~/qmac/analyses/new_analyses
cd ~/qmac/analyses/new_analyses
mkdir figures
# cp 

# get genes containing SNPs

grep -f <(cut -f1 thin_500bp.bim | sort -u) <(cut -f1-4 ~/qmac/analyses/blast/sorted-Qlobatagenes-PCG.bed) | sed 's/N._//g; s/\.1//g' > allgenes_chrbp_thin_500bp.txt

plink --bfile thin_500bp --double-id --allow-extra-chr --set-missing-var-ids "@:#" --make-set allgenes_chrbp_thin_500bp.txt --write-set --out allgenes
grep -E '^[0-9]' allgenes.set | wc -l # number of genes falling 

# make VCF file of SNPs

plink --bfile thin_500bp --double-id --allow-extra-chr --set-missing-var-ids @:"#" --recode vcf-iid --keep-allele-order --out thin_500bp
cut -f1,4 thin_500bp.bim | sed 's/^044/NC_044/g;s/^022/NW_022/g;s/\t/.1\t/g' > snpsforvcf.txt
vcftools --vcf ~/qmac/analyses/it_filt/geno80pc.vcf --positions snpsforvcf.txt --recode --recode-INFO-all --stdout > thin_500bp.vcf
vcftools --vcf thin_500bp.vcf --site-mean-depth --stdout > thin_500bp_sitedepth.txt
conda activate R
Rscript ~/qmac/scripts/depthstats.R "thin_500bp_sitedepth.txt"
conda deactivate

# # get summaries on genotyping, depth, etc

# conda activate R
# Rscript ~/qmac/scripts/

# format admixtable

conda activate R
Rscript ~/qmac/scripts/admixtable_format.R K5.Q K6.Q indslist.txt
conda deactivate

# generate admixture plots

sed -n '/divergences/,/CV/{//!p}' ~/qmac/analyses/ldprune_qc/thin_admix/$(grep 'K6' ~/qmac/analyses/ldprune_qc/thin_admix/thin_it_filt_500bp_bestadmixrun_perK.txt | sed 's/Q$/out/') | sed 's/\t$//g;s/\t/,/g;1b;s/$/,0/g' | sed '1s/$/,Pop5/;2s/$/,,,,,/;3s/$/,,,,/;4s/$/,,,/;5s/$/,,/;6s/$/,/' > admix_K6_fstmatrix.csv
cut -f1-2 K6_clusterlabels.txt | sed '1d;s/\t/,/g' | while read NAME
do
  sed -i "s/${NAME#*,}/${NAME%,*}/g" admix_K6_fstmatrix.csv
done

conda activate R
Rscript ~/qmac/scripts/plot_admix_K5.R "K5_admixtable_allplinkinds.txt" "./figures/"
Rscript ~/qmac/scripts/plot_admix_K6.R "K6_admixtable_allplinkinds.txt" "./figures/"
Rscript ~/qmac/scripts/plot_admix_K6_fst.R "admix_K6_fstmatrix.csv" "./figures/"
conda deactivate

# NMDS

plink --bfile thin_500bp --genome --allow-extra-chr --double-id --set-missing-var-ids "@:#" --out thin_500bp
plink --bfile thin_500bp --read-genome thin_500bp.genome --cluster --mds-plot 3 --allow-extra-chr --double-id --set-missing-var-ids "@:#" --out thin_500bp
plink --bfile thin_500bp --read-genome thin_500bp.genome --cluster --mds-plot 2 --allow-extra-chr --double-id --set-missing-var-ids "@:#" --out thin_500bp_2ax

conda activate R_envdata
Rscript ~/qmac/scripts/mds_snps.R "thin_500bp.mds" "K6_admixtable_allplinkinds.txt" "./figures/thin_500bp_mdsplot.pdf"
conda deactivate

# PCA

plink --bfile thin_500bp --double-id --allow-extra-chr --set-missing-var-ids "@:#" --pca --out thin_500bp

conda activate R_envdata
Rscript ~/qmac/scripts/pca_snps.R "thin_500bp.eigenval" "thin_500bp.eigenvec" "./figures/thin_500bp_pcaplot.pdf"
conda deactivate

# join coordinates with updated admixture proportions

conda activate R
Rscript ~/qmac/scripts/join_coords.R "~/qmac/coords_noSHF.csv" "K6_admixtable_allplinkinds.txt"
conda deactivate

# join environmental data with admixture/coordinate data

conda activate R
Rscript ~/qmac/scripts/merge_envdata_new.R
conda deactivate

# conduct global Moran's I test

conda activate R_envdata
Rscript ~/qmac/scripts/moran_i_mac.R "qmac_formodeling_alldata.csv"
conda deactivate
mv moranI_qmac.txt figures
mv moranI_qmac_permutationdist.pdf figures

# run GLMMs

conda activate R_envdata
Rscript "admix_formodeling_alldata.csv"
Rscript ~/qmac/scripts/model_glmm.R "admix_formodeling_alldata.csv" "qmac_formodeling_alldata.csv"

## NewHybrids

mkdir newhybrids
cd newhybrids
conda activate R_newhybrids

# subset individuals for determining pairwise SNP panels - pairwise list of individuals from 2 spp with >0.95 ADMIXTURE cluster assignment to determined species (essentially unadmixed)

Rscript ~/qmac/scripts/inds_for_hybriddetective.R ../K6_admixtable_allplinkinds.txt

# recode pairwise PLINK ped files from unpruned bedfile master list of SNPs, convert from ped to genepop format

cp ~/qmac/analyses/newhybrids/ped2genepop.spid .
while read NAME
do
    plink --bfile ~/qmac/analyses/ldprune_qc/it_filt --double-id --allow-extra-chr --set-missing-var-ids "@:#" --keep-fam inds_for_hybriddetective_$NAME.txt --recode 12 --out unpruned_$NAME
    cut -d' ' -f1 --complement unpruned_$NAME.ped | paste -d'_' <(cut -d' ' -f1 unpruned_$NAME.ped | cut -d'_' -f2) - | paste -d' ' <(cut -d' ' -f1 unpruned_$NAME.ped | cut -d'_' -f2) - > tmp_unpruned_$NAME.ped && mv -f tmp_unpruned_$NAME.ped unpruned_$NAME.ped
    sed "s/PED_PARSER_MAP_FILE_QUESTION=/&unpruned_$NAME.map/" ped2genepop.spid > ped2genepop_$NAME.spid
    PGDSpider2-cli -inputfile unpruned_$NAME.ped -inputformat PED -outputfile genepop_unpruned_$NAME.txt -outputformat GENEPOP -spid ped2genepop_$NAME.spid
done < speciespairs.txt

# get the 200 most differentiated unlinked loci for each species pair and move to separate files

while read NAME
do
  Rscript ~/qmac/scripts/snps_for_newhybrids.R genepop_unpruned_${NAME}.txt --vanilla 2>&1 | tee snps_for_newhybrids_$NAME.log
  sed '1d' genepop_unpruned_${NAME}_200_Loci_Panel.txt | head -200 > snps_for_newhybrids_${NAME}.txt
done < speciespairs.txt

# get all individuals with Q values assigned >.95 in species pairs (admixed and unadmixed)

Rscript ~/qmac/scripts/inds_for_newhybrids.R ../K6_admixtable_allplinkinds.txt

# subset master PLINK bed file for individuals, loci for each species pair, convert to NewHybrids file format

cp ~/qmac/analyses/newhybrids/ped2newhybrids.spid .
sed -i 's/QUE002099_bic/QUE002099_mue/g' ~/qmac/analyses/ldprune_qc/it_filt.fam
mkdir run_nh
# while read NAME
# do
#   plink --bfile ~/qmac/analyses/ldprune_qc/it_filt --double-id --allow-extra-chr --set-missing-var-ids "@:#" --keep-fam inds_for_newhybrids_${NAME}.txt --extract snps_for_newhybrids_${NAME}.txt --recode 12 --out ${NAME}_nhinput
#   sed "s/PED_PARSER_MAP_FILE_QUESTION=/&${NAME}_nhinput.map/" ped2newhybrids.spid > ped2newhybrids_$NAME.spid
#   mkdir run_nh/$NAME
#   for i in {1..3}
#   do
#     mkdir run_nh/$NAME/run$i
#   done
#   PGDSpider2-cli -inputfile ${NAME}_nhinput.ped -inputformat PED -outputfile run_nh/$NAME/${NAME}.newhybrids -outputformat NEWHYBRIDS -spid ped2newhybrids_$NAME.spid
# done < speciespairs.txt

while read NAME
do
  plink --bfile ~/qmac/analyses/ldprune_qc/it_filt --double-id --allow-extra-chr --set-missing-var-ids "@:#" --keep-fam inds_for_newhybrids_${NAME}.txt --extract snps_for_newhybrids_${NAME}.txt --recode 12 --out ${NAME}_nhinput
  sed "s/PED_PARSER_MAP_FILE_QUESTION=/&${NAME}_nhinput.map/" ped2newhybrids.spid > ped2newhybrids_$NAME.spid
  mkdir run_nh_new/$NAME
  for i in {1..4}
  do
    mkdir run_nh_new/$NAME/run$i
  done
  PGDSpider2-cli -inputfile ${NAME}_nhinput.ped -inputformat PED -outputfile run_nh_new/$NAME/${NAME}.newhybrids -outputformat NEWHYBRIDS -spid ped2newhybrids_$NAME.spid
done < speciespairs.txt

# make sure genotype frequency matrix is in newhybrids directory

cp ~/qmac/analyses/newhybrids/GtypeFreqExtBxs.txt .
cut -f1 GtypeFreqExtBxs.txt | sed '1d' | paste -s > GtypeFreqExtBxs_headers.txt

# run newhybrids in parallel

# parallel 'cd run_nh/{1}/run{2} && ~/newhybrids/newhybrids-no-gui-linux.exe -d ../{1}.newhybrids -c ../../../GtypeFreqExtBxs.txt --burn-in 300000 --num-sweeps 600000 --seeds {2}133 93{2}7 --no-gui --print-traces Pi 5 2>&1 | tee {1}_run{2}.log' :::: speciespairs.txt ::: {1..3}
parallel 'cd run_nh_new/{1}/run{2} && ~/newhybrids/newhybrids-no-gui-linux.exe -d ../{1}.newhybrids -c ../../../GtypeFreqExtBxs.txt --burn-in 300000 --num-sweeps 600000 --seeds {2}133 93{2}7 --no-gui --print-traces Pi 5 2>&1 | tee {1}_run{2}.log' :::: speciespairs.txt ::: {1..4}

# parallel 'cd run_nh_new/alb_mue/run{} && ~/newhybrids/newhybrids-no-gui-linux.exe -d ../alb_mue.newhybrids -c ../../../GtypeFreqExtBxs.txt --burn-in 300000 --num-sweeps 600000 --seeds {}133 93{}7 --no-gui --print-traces Pi 5 2>&1 | tee alb_mue_run{}.log' ::: {1..4}

# rename output files (WIP)

# while read NAME
# do
#   for i in run1 run2 run3
#   do
#     for j in $(find run_nh/$NAME/${i} -name "aa-*")
#     do
#       mv $j $(echo $j | sed "s/aa-/${NAME}_${i}_/g")
#     done
#   done
# done < speciespairs.txt

while read NAME
do
  for i in run1 run2 run3 run4
  do
    for j in $(find run_nh_new/$NAME/${i} -name "aa-*")
    do
      mv $j $(echo $j | sed "s/aa-/${NAME}_${i}_/g")
    done
  done
done < speciespairs.txt

# create input files for traceplots

# while read NAME
# do
#   for i in run1 run2 run3
#   do
#     grep 'PI_TRACE' run_nh/$NAME/${i}/${NAME}_${i}.log | sed '1d' | cut -d':' -f2- | awk -v run=$i '{print $0"\t"run}'
#   done | sed '1i itnum\tPure_0\tPure_1\tF1\tF2\tP0_F1\tP1_F1\tP0_P0F1\tP1_P1F1\tP0_P0P0F1\tP1_P1P1F1\trunnum' > run_nh/${NAME}_pitraces.txt
# done < speciespairs.txt

while read NAME
do
  for i in run1 run2 run3 run4
  do
    grep 'PI_TRACE' run_nh_new/$NAME/${i}/${NAME}_${i}.log | sed '1d' | cut -d':' -f2- | awk -v run=$i '{print $0"\t"run}'
  done | sed '1i itnum\tPure_0\tPure_1\tF1\tF2\tP0_F1\tP1_F1\tP0_P0F1\tP1_P1F1\tP0_P0P0F1\tP1_P1P1F1\trunnum' > run_nh_new/${NAME}_pitraces.txt
done < speciespairs.txt

# plot traces for pi (estimated distribution of genotype frequency classes)

# conda activate R_newhybrids
# while read NAME
# do
#   Rscript ~/qmac/scripts/plot_nhpitrace.R run_nh/"${NAME}_pitraces.txt" ../figures/"pitraceplot_${NAME}.pdf"
# done < speciespairs.txt
# conda deactivate

conda activate R_newhybrids
while read NAME
do
  Rscript ~/qmac/scripts/newhybrids_plotpitrace.R run_nh_new/"${NAME}_pitraces.txt" "${NAME}_pitraceplot.pdf"
done < speciespairs.txt
conda deactivate

# merge traceplots into single PDF

pdfunite $(cat speciespairs.txt | tr '\n' ' \n' | sed 's/ /_pitraceplot.pdf /g;s/ $/\n/g') ../figures/allpairs_pitraceplots.pdf
rm *pitraceplot.pdf

# Rscript ~/qmac/scripts/plot_nhpitrace.R run_nh/"mac_alb_pitraces.txt" ../figures/"pitraceplot_mac_alb.pdf"

# format hybrid class PP files

while read NAME
do
  paste <(sed '1i ind_ID' inds_for_newhybrids_$NAME.txt) <(cut --complement -f1-2 run_nh_new/$NAME/run1/${NAME}_run1_PofZ.txt | sed '1d' | cat GtypeFreqExtBxs_headers.txt -) > ${NAME}_nh_rawlabels.txt
done < speciespairs.txt

# relabel classes with species

conda activate R_newhybrids
while read NAME
do
  Rscript ~/qmac/scripts/newhybrids_fixlabels.R ${NAME}_nh_rawlabels.txt ${NAME}_nh.txt
  rm ${NAME}_nh_rawlabels.txt
done < speciespairs.txt
conda deactivate

# plot newhybrids barplots

conda activate R_newhybrids
while read NAME
do
  Rscript ~/qmac/scripts/newhybrids_plot.R $NAME ${NAME}_nh.txt ${NAME}_newhybrids.pdf
done < speciespairs.txt
conda deactivate

# merge pdfs

pdfunite $(cat speciespairs.txt | tr '\n' ' \n' | sed 's/ /_newhybrids.pdf /g;s/ $/\n/g') ../figures/allpairs_newhybrids.pdf
rm *newhybrids.pdf

# get NewHybrids text summaries

conda activate R_newhybrids
Rscript ~/qmac/scripts/newhybrids_summary.R speciespairs.txt ../figures/
conda deactivate

## SRA

diff <(cut -f1 -d '_' indslist_allsampsfb.txt | sort) <(cut -f2 -d ',' OakSpmTable-Extract-20240904-HybSeqMetadata.csv | sed '1d' | sort)


## jackknife distributions

mkdir subsample_dists
cd subsample_dists

# get genes with chr positions for all SNPs

grep -f <(cut -f1 ../thin_500bp.bim | sort -u) <(cut -f1-4 ~/qmac/analyses/blast/sorted-Qlobatagenes-PCG.bed) | sed 's/N._//g; s/\.1//g' > allgenes_chrbp.txt

# get set file with lists of SNPs within each gene region

plink --bfile ../thin_500bp --double-id --allow-extra-chr --set-missing-var-ids "@:#" --make-set allgenes_chrbp.txt --write-set --out allgenes

# get number of SNPs per gene

grep '^T' allgenes.set | while read LINE
do
  printf "$(sed "0,/$LINE/d;/END/Q" allgenes.set | wc -l)\t$LINE\n"
done | sed '1i n_snps\tgene' > numsnpspergene.txt

# get corresponding functional class for all genes

for i in other bud drought icehos
do
  sed "s/$/\t$i/g" ~/qmac/analyses/subsample_dists/${i}genes_ids.txt
done | sed '1i gene\tclass' > allgenes_ids.txt

# # get genes per class with at least one SNP for subsampling

# conda activate R_envdata
# Rscript ~/qmac/scripts/get_subsamplegenes.R
# conda deactivate R_envdata

# get genes per class with at least 3 SNPs for subsampling

conda activate R_envdata
Rscript ~/qmac/scripts/get_subsamplegenes_3snps.R
conda deactivate R_envdata

# annotate SNPs in target genes with gene name

csplit -kfs allgenes.set '/^T/' '{*}'
rm s00
for i in s*
do
  mv $i $(sed -n '1p' $i)_snps.txt
done
for i in *snps.txt
do
  sed "1d;s/$/,${i/_snps.txt/}/g" $i | head -n -2
done > allsnps_annotated.txt
rm *snps.txt

# mkdir subset_admix
mkdir subset_fst

# generate 10000 jackknife replicates (20 genes subsampled per class, 1 variant subsampled per gene) for pairwise Fst and admixture distributions

# for j in subset_fst subset_admix
# do
#   for NAME in othergenes budgenes droughtgenes
#   do
#     for i in {1..10000}
#     do
#       sort -R ${NAME}_idsforsub.txt | head -20 | while read LINE
#       do
#         grep $LINE allsnps_annotated.txt | sort -R | head -1
#       done | cut -f1 -d',' > $j/${NAME}_subset${i}.txt
#     done
#   done
# done

# # parallelized version

# parallel -j96 'for j in subset_fst subset_admix; do for NAME in othergenes budgenes droughtgenes; do sort -R "$NAME"_idsforsub.txt | head -20 | while read LINE; do grep "$LINE" allsnps_annotated.txt | sort -R | head -1; done | cut -f1 -d',' > "$j"/"$NAME"_subset{}.txt; done; done' ::: {1..10000}

# parallelized version

parallel -j96 'for NAME in othergenes budgenes droughtgenes; do sort -R "$NAME"_idsforsub.txt | head -15 | while read LINE; do grep "$LINE" allsnps_annotated.txt | sort -R | head -3; done | cut -f1 -d',' > subset_fst/"$NAME"_subset{}.txt; done' ::: {1..10000}

# create subset lists

for i in subset_fst/*subset*.txt
do
  echo ${i%.txt} | sed 's/.*\///g'
done > samplinglist.txt

# create file of admixed individuals to remove from distributions

# conda activate R
# Rscript ~/qmac/scripts/indsbelowQ95K6.R --vanilla
# conda deactivate

# create plink bed files for each subset, remove admixed individuals (conspecific Q < 0.95) identified

# for j in subset_fst subset_admix
# do
#   while read NAME
#   do
#     plink --bfile ../thin_500bp --double-id --allow-extra-chr --set-missing-var-ids '@:#' --keep-fam ../indsaboveQ95K6.txt --make-bed --extract $j/${NAME}.txt --out $j/$NAME
#   done < samplinglist.txt
# done

# parallelized version

parallel -j96 'plink --bfile ../thin_500bp --double-id --allow-extra-chr --set-missing-var-ids @:# --keep-fam ../indsaboveQ95K6.txt --make-bed --extract subset_fst/{}.txt --out subset_fst/{}' :::: samplinglist.txt

# generate pairwise Fst jackknife distribution

cut -f1 -d' ' subset_fst/othergenes_subset1.fam | cut -f2 -d'_' > poplist.txt

conda activate R_Fstat
Rscript ~/qmac/scripts/fstdist.R poplist.txt subset_fst/ ../figures/pairwise_fst_jackknife_10000.pdf ../figures/pairwise_fst_pvals.txt
conda deactivate

# admix/Fst jackknife distribution

mkdir subset_admix/admix
cd subset_admix/admix

# run admixture with 10 replicates (K=6)

parallel -j96 'for i in {1..10}; do admixture -s time ../subset_admix/plink_{}.bed 6 > {}.Run$i.out; mv plink_{}.6.Q {}.Run$i.Q; mv plink_{}.6.P {}.Run$i.P; done' :::: ../samplinglist.txt

# find run with highest loglikelihood per pseudoreplicate and place in file

while read NAME
do
  for i in {1..10}
  do
    printf "${NAME}.Run$i.out\t$(grep '^Loglikelihood' ${NAME}.Run${i}.out | sed 's/.*\s//g')\n"
  done | awk 'NR==1{max=$2;line=$0} $2>max{max=$2;line=$0} END{print line}' >> all_bestadmixruns.txt
done < ../samplinglist.txt

sed 's/\.out.*//g' all_bestadmixruns.txt > all_bestadmixruns_fnameonly.txt

# separate best runs by gene class into namefiles

for i in othergenes budgenes droughtgenes
do
  grep "$i" all_bestadmixruns_fnameonly.txt > "${i}_bestadmixruns.txt"
done

# align all best runs

conda activate R
Rscript alignk_within_nomerge.R all_bestadmixruns_fnameonly.txt
conda deactivate

# manually assign clusters and print as list in textfile (ditto for gene classes)

printf "mac1_frac\nbic_frac\nste_frac\nmac2_frac\nmue_frac\nalb_frac\n" > clusterlabels.txt
printf "othergenes\nbudgenes\ndroughtgenes\n" > geneclasslist.txt

# generate plot of admix distribution

Rscript ~/qmac/scripts/admixdist.R all_bestadmixruns_fnameonly.txt geneclasslist.txt clusterlabels.txt ~/qmac/poplist.txt admix_dist_5000reps

# # for X-ORIGIN

## FEEMS reanalysis

# filter PLINK bed file to only Q. macrocarpa individuals with conspecific Q > 0.95, remove invariant sites

mkdir feems
cd feems
cat ../qmac_aboveK95.txt | while read LINE
do
  grep "$LINE" ../admixcoords.csv | cut -f1,9-10 -d',' | sed 's/,/ /g'
done > qmac_feems_indscoords.txt
plink --bfile ../thin_500bp --double-id --allow-extra-chr --set-missing-var-ids '@:#' --keep-fam <(cut -f1 -d' ' qmac_feems_indscoords.txt) --maf --make-bed --out qmac_feems
awk '{print $3,$2}' qmac_feems_indscoords.txt > qmac.coord
python ~/qmac/scripts/qmac_feems.py 2>&1 | tee qmac_feems.log
mv qmac_feems*pdf ../figures
cd ..

## admixture maps

# add site to admixcoords.csv

sed '1d' admixcoords.csv | cut -f1 -d',' | while read LINE
do
  printf "$(grep $LINE admixcoords.csv),$(grep $LINE ~/qmac/analyses/envdata/sitedata/sitedata_allfbindividuals.csv | cut -f2 -d',')\n"
done | sed "1i $(head -1 admixcoords.csv),site" > admixsitecoords.csv

# generate map of species ranges and sampling sites

conda activate R_envdata
Rscript ~/qmac/scripts/map_sppranges_sites.R "admixsitecoords.csv" "./figures/"
conda deactivate

# generate table of sampling sites

conda activate R_envdata
Rscript ~/qmac/scripts/make_site_table.R "admixsitecoords.csv" "./figures/"
conda deactivate

# generate table of individuals (for supplements)

conda activate R_envdata
Rscript ~/qmac/scripts/make_ind_table.R "K6_admixtable_allplinkinds.txt" "admixsitecoords.csv" "quercus_hybseq_individuals_metadata.csv" "./figures/"
conda deactivate

# generate all species admixture, Q. macrocarpa population structure maps

conda activate R_envdata
Rscript ~/qmac/scripts/map_admix_allspp.R "admixsitecoords.csv" "./figures/"
Rscript ~/qmac/scripts/map_admix_onlymac.R "admixsitecoords.csv" "./figures/"
conda deactivate

# generate chromosome map of SNPs/target loci

conda activate R_envdata
Rscript ~/qmac/scripts/plot_snpstargets.R "thin_500bp.bim" "./figures/"
conda deactivate