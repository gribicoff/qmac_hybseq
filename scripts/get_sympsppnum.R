args <- commandArgs(trailingOnly=T)

library(sf)
library(rgbif)
library(CoordinateCleaner)
library(dplyr)

# args = c("../coords_noSHF.csv",fiel path to output csv,)

# taxonKey: 2878213 (mac), 2879737 (alb), 2878470 (bic), 8879961 (ste), 2879947 (mue)

sppkeys <- c("2878213","2879737","2878470","8879961","2879947")

# download occurrence data from GBIF

gbif.dl <- occ_download(
  pred_in("taxonKey",sppkeys),
  pred("hasCoordinate",T),
  pred("hasGeospatialIssue",F),
  pred("occurrenceStatus","present"),
  pred_in("country",c("US","CA")),
  pred_not(pred_in("basisOfRecord",c("fossilSpecimen","livingSpecimen"))),
  format="SIMPLE_CSV"
)

occ_download_wait(gbif.dl)

gbif.df <- occ_download_import(occ_download_get(gbif.dl))

# clean coordinates

gbif.df <- setNames(gbif.df,tolower(names(gbif.df)))
gbif.df <- filter(gbif.df,year >= 1900)
gbif.df <- filter(gbif.df,coordinateprecision < 0.01 | is.na(coordinateprecision))
gbif.df <- filter(gbif.df,coordinateuncertaintyinmeters < 1000 | is.na(coordinateuncertaintyinmeters))
gbif.df <- filter(gbif.df,!coordinateuncertaintyinmeters %in% c(301,3036,999,9999))
gbif.df <- cc_zero(gbif.df,buffer=1,value="clean",verbose=T)
gbif.df <- cc_cen(gbif.df,buffer=2000,value="clean",test="both",verbose=T)
gbif.df <- cc_cap(gbif.df,buffer=2000,value="clean",verbose=T)
gbif.df <- cc_inst(gbif.df,buffer=2000,value="clean",verbose=T)
gbif.df <- cc_sea(gbif.df,value="clean")
gbif.df_nourb <- cc_urb(gbif.df,value="clean")
gbif.df <- gbif.df[,c("decimallongitude","decimallatitude","species")]
names(gbif.df) <- c("long","lat","spp")
gbif.df$spp <- substr(gsub("^.* ","",gbif.df$spp),1,3)
gbif.df_nourb <- gbif.df_nourb[,c("decimallongitude","decimallatitude","species")]
names(gbif.df_nourb) <- c("long","lat","spp")
gbif.df_nourb$spp <- substr(gsub("^.* ","",gbif.df_nourb$spp),1,3)

# write CSVs of cleaned GBIF occurrence data with and without urban areas included

write.csv(gbif.df,"gbif_allspp.csv",quote=F,row.names=F)
write.csv(gbif.df_nourb,"gbif_allspp_nourb.csv",quote=F,row.names=F)

# read in coordinates with 20km buffered radius - convert to Albers Equal Area Conic proj for eastern US (same as Little map shapefiles)

proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"

all.df <- read.csv(args[1],header=T)
crds <- st_transform(st_as_sf(all.df,coords=c("long","lat"),crs=4326),crs=proj)
samp.sf <- st_buffer(crds,dist=20000)

# find points (as row numbers) within each buffered 20km radius, find unique spp, remove conspecific, count number of sympatric spp

gbif.sf <- st_transform(st_as_sf(gbif.df,coords=c("long","lat"),crs=4326),crs=proj)
sympspp <- st_intersects(samp.sf,gbif.sf) # WIP
names(sympspp) <- all.df$ind_ID
sympspp <- lapply(sympspp,function(x) {unique(gbif.df[c(x),]$spp)})
# extracts last n characters from string
substr_end <- function(x,n) {
  substr(x,nchar(x)-n+1,nchar(x))
}
sympspp <- sapply(names(sympspp),function(x) {sympspp[[x]] <- subset(sympspp[[x]],sympspp[[x]] != substr_end(x,3))},USE.NAMES=T)
sympspp <- lapply(sympspp,length)

# convert to dataframe and write to table

sympspp.df <- data.frame(ind_ID=names(sympspp),lat=all.df$lat,long=all.df$long,num_sympspp=do.call(rbind,c(sympspp,use.names=F)))
write.csv(sympspp.df,args[2],quote=F,row.names=F)
