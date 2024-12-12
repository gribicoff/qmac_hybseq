library(sf)
library(rgbif)
library(CoordinateCleaner)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)

# GBIF taxonKeys: 2878213 (mac), 2879737 (alb), 2878470 (bic), 8879961 (ste), 2879947 (mue)

sppkeys <- c("2878213","2879737","2878470","8879961","2879947")
spp <- c("mac","alb","bic","ste","mue")

# download occurrence data from GBIF - make sure GBIF API key is set

gbif_dl <- occ_download(
  pred_in("taxonKey",sppkeys),
  pred_in("country",c("US","CA")),
  pred("HAS_GEOSPATIAL_ISSUE",FALSE),
  pred("HAS_COORDINATE",TRUE),
  pred("OCCURRENCE_STATUS","PRESENT"),
  pred_not(pred_in("BASIS_OF_RECORD",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
  format="SIMPLE_CSV"
)
occ_download_wait(gbif_dl)

# import occurrence data from GBIF query into R and clean records with CoordinateCleaner

gbif.df <- gbif_dl %>%
  occ_download_get() %>%
  occ_download_import() %>%
  setNames(tolower(names(.))) %>%
  dplyr::filter(year >= 1900) %>%
  dplyr::filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>%
  dplyr::filter(coordinateuncertaintyinmeters < 1000 | is.na(coordinateuncertaintyinmeters)) %>%
  dplyr::filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>%
  cc_zero(buffer = 1,value = "clean",verbose = TRUE) %>%
  cc_cen(buffer = 2000,value = "clean",test = "both",verbose = TRUE) %>%
  cc_cap(buffer = 2000,value = "clean",verbose = TRUE) %>%
  cc_inst(buffer = 2000,value = "clean",verbose = TRUE) %>%
  cc_sea(value = "clean") %>%
  cc_urb(value = "clean") %>%
  dplyr::select(decimallongitude,decimallatitude,species) %>%
  dplyr::rename(long = decimallongitude, lat = decimallatitude) %>%
  dplyr::mutate(species = substr(gsub("^.* ","",species),1,3))

# create set of dataframes filtered for each species, write each to separate CSV

occ.ls <- sapply(spp,function(x) {
  gbif.df %>% dplyr::filter(species == x)
},USE.NAMES = TRUE,simplify = FALSE)
sapply(spp,function(x) {
  write.csv(occ.ls[[x]],file = paste0("rawgbifoccs_",x,".csv"),row.names = FALSE,quote = FALSE)
})

# plot raw occurrence data, one map each per species

uscan <- ne_states(c("United States of America","Canada"),returnclass="sf")
uscan <- st_crop(uscan,c(xmin=-140,xmax=-55,ymin=20,ymax=80))
uscan <- st_transform(uscan,crs=4326)

sapply(spp,function(x) {
  p <- ggplot() +
    geom_sf(data=uscan) +
    geom_point(data = occ.ls[[x]],aes(x=long,y=lat),size=0.1) +
    ggtitle(paste0("Q",toupper(x)))
  p %>% ggsave(filename=paste0("rawgbifoccs_",x,".png"))
})

# save raw occurrence datapoints as shapefiles to be manually curated in QGIS

sapply(spp,function(x) {
  occ.ls[[x]] %>% dplyr::select(long,lat) %>%
  st_as_sf(crs=4326,coords=c("long","lat")) %>%
  st_write(paste0("rawgbifoccs_",x,".shp"))
})


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

# convert to dataframe and write to CSVs

sympspp.df <- data.frame(ind_ID=names(sympspp),lat=all.df$lat,long=all.df$long,num_sympspp=do.call(rbind,c(sympspp,use.names=F)))
write.csv(sympspp.df,args[2],quote=F,row.names=F)
