#!/bin/bash

# WIP!!!

mkdir ~/qmac/analyses/envdata
cd ~/qmac/analyses/envdata
conda activate R_envdata

# make sure that file with individual IDs and coordinates is in directory - append species endings to IDs

# cut -d',' -f1,1 ~/qmac/coords_noSHF.csv | sed '1d' | while read LINE
# do
#   grep "$LINE" ~/qmac/indslist_allsampsfb.txt
# done | sed '1i ind_ID' | paste - <(awk -F',' '{OFS=",";print $2,$3}' ~/qmac/coords_noSHF.csv) > ~/qmac/tmp_coords_noSHF.csv
# rm ~/qmac/coords_noSHF.csv
# mv ~/qmac/tmp_coords_noSHF.csv ~/qmac/coords_noSHF.csv
# cp ~/qmac/coords_noSHF.csv .

# download bioclim/elevation data, extract data for coords

mkdir climdata
cd climdata

wget "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip"
unzip wc2.1_30s_bio.zip
wget "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip"
unzip wc2.1_30s_elev.zip
# find . -name "*.tif" -printf '%f\n' > rasterlist.txt

Rscript ~/qmac/scripts/get_bioclimdata.R ../coords_noSHF.csv
cd ..

# obtain soil data from SDA

mkdir soildata
cd soildata

Rscript ~/qmac/scripts/get_SDA_from_GPS_Paolucci.R "~/qmac/analyses/envdata/coords_noSHF.csv" "~/qmac/analyses/envdata/soildata"

cd ..

# subset to only individuals in admixture analysis

# sed '1d;s/_.*//g' ~/qmac/.txt | while read LINE
# do
#   grep "$LINE" FinalResults_SDA_from_GPS.csv
# done | cat <(sed -n '1p' FinalResults_SDA_from_GPS.csv) - > soildata_noSHF.csv
#
# manually inspect for NAs, erroneous coordinates (eg, water)

# change coordinates for soil data, upload new dataset

# obtain number of sympatric species per sample from GBIF

mkdir sympsppdata
cd sympsppdata

Rscript ~/qmac/scripts/get_sympsppnum.R ../coords_noSHF.csv ../coords_noSHF_sympsppnum.csv

cd ..

# obtain distance from range edge/centroid (make sure shapefiles are in directory)

mkdir rangedata
cd rangedata

find . -name "*little.shp" -printf '%f\n' > rangelist.txt
# Rscript ~/qmac/scripts/get_distfromcentroid.R ../coords_noSHF
Rscript ~/qmac/scripts/get_distfromrangeedge.R ../coords_noSHF.csv ../coords_noSHF_distfromrangeedge.csv
Rscript ~/qmac/scripts/get_sympsppbin.R ../coords_noSHF.csv ../coords_noSHF_sympsppbin.csv

cd ..

# download cleaned FIA data for building SDM - see https://gitlab.com/meireles/compile_fia_data/-/tree/master/data/clean?ref_type=heads for data, which are assumed to have been loaded into R and then exported as CSV files

Rscript merge_occdata.R

# run ENMeval R scripts - memory usage has been a consistent issue, with the kernel killing the R session

# for i in mac alb bic ste mue
# do
#   Rscript ~/qmac/scripts/maxent_run.R $i --vanilla --verbose 2>&1 | tee ENMeval_$i.log
# done

# get occurrence data, background data from GBIF and apply preliminary filters

Rscript ~/qmac/scripts/get_rawgbif_targetbg.R

# manually clean in QGIS...

# import manually curated occurrence shapefiles and export as CSVs

Rscript ~/qmac/scripts/get_cleanedgbif_targetbg.R --vanilla

# create convex hulls around occurrence points, remove lakes, crop rasters, downsample background points, and export input files for Maxent

Rscript ~/qmac/scripts/maxent_prep_targetbg.R --vanilla

# run Maxent across combinations of tuning parameters via ENMeval

for i in mac alb bic ste mue
do
  Rscript ~/qmac/scripts/maxent_run_final.R $i --vanilla --verbose 2>&1 | tee ENMeval_$i.log
done

for i in mac alb bic ste mue
do
  Rscript ~/qmac/scripts/maxent_process.R $i --vanilla --verbose 2>&1 | tee modres_$i.log
done

# generate maps of habitat suitability across ENA

Rscript ~/qmac/scripts/maxent_plot.R --vanilla --verbose

conda deactivate

# obtain site data - requires "sitedata" subdirectory with textfile containing site codes for each individual in dataset

cd sitedata
grep -f <(cut -d'_' -f1 ~/qmac/coords_noSHF.csv) sitedata_allfbindividuals.csv > ../coords_noSHF_sitedata.csv
cd ..

# copy modified Q matrix with admixture proportions to directory, run R script to 

sed 's/\t/,/g' ~/qmac/.txt > .csv
Rscript ~/qmac/scripts/clean_admixenvdata.R .csv ~/qmac/analyses/model_admix/mactot95_inds.txt

# run R script to merge envdata into single csv

# find . -name "coords_noSHF_*" -printf '%f\n' > filelist.txt
ls coords_noSHF_* > filelist.txt
Rscript ~/qmac/scripts/merge_envdata.R filelist.txt admixcoords_formodeling.csv maccoords_formodeling.csv

# raw_predictormatrix.csv

cut -f1 ~/qmac/.txt | sed '1d' | while read LINE
do
  printf "$(grep $LINE coords_noSHF_bioclim.csv),$(grep $LINE coords_noSHF_SDA_onlyadmix.csv | cut --complement -d',' -f1-4,6-8,11-12,27-28),$(grep $LINE coords_noSHF_sympsppnum.csv | cut -d',' -f4),$(grep $LINE coords_noSHF_distfromrangeedge.csv | cut -d',' -f4),$(echo $LINE | cut -d'_' -f2),$(grep $LINE ~/qmac/admixtableK6_allplinkinds.txt | cut -f3-8 --output-delimiter=','),$(grep $LINE ~/qmac/admixtableK6_allplinkinds.txt | awk '{print $7+$8}')\n"
done | sed "1i $(head -1 coords_noSHF_bioclim.csv),$(head -1 coords_noSHF_SDA_onlyadmix.csv | cut --complement -d',' -f1-4,6-8,11-12,27-28),num_sympspp,dist_m,sppcode,$(head -1 ~/qmac/admixtableK6_allplinkinds.txt | cut -f3-8 --output-delimiter=','),mactot_frac" > _allenvdata.csv
