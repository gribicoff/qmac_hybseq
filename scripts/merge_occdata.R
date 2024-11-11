args <- commandArgs(trailingOnly=T)

library(sf)

# REWRITE!!!

# args <- c("~/qmac/analyses/envdata/sympsppdata/fiadata.rds","~/qmac/analyses/envdata/sympsppdata/gbif_allspp.csv","~/qmac/analyses/envdata/rangedata/")

# read in FIA, GBIF occurrence data

fia <- readRDS(args[1])
gbif <- read.csv(args[2],header=T)

occ_fia <- fia[["trees"]][,c("LON","LAT","SPCD")]
colnames(occ_fia) <- c("long","lat","spp")
rownames(occ_fia) <- NULL

# store FIA spp codes (in order of spp vector)

spp <- c("mac","alb","bic","ste","mue")
fiacodes <- c("823","802","804","835","826")

# create list of dataframes with deduplicated species-specific FIA occurrences, replaced spp codes

occ_fia <- lapply(fiacodes,function(x) {
  unique(occ_fia[occ_fia$spp == x,-3])
})
names(occ_fia) <- spp
occ_fia <- sapply(spp,function(x) {
  occ_fia[[x]]$spp <- rep(x,nrow(occ_fia[[x]]));
  return(occ_fia[[x]])
},USE.NAMES=T,simplify=F)

# deduplify GBIF occurrences, convert GBIF df to list of species-specific dfs

gbif_occ <- sapply(spp,function(x) {
  unique(gbif[gbif$spp == x,])
},USE.NAMES=T,simplify=F)

# read in range shapefiles

proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
shp <- sapply(spp,function(x) {
  st_cast(st_set_crs(st_read(paste0(args[3],x,"_little.shp")),proj),"MULTIPOLYGON")
},USE.NAMES=T,simplify=F)
shp[["bic"]] <- shp[["bic"]][-12,] # removes Lake Winnebago
shp[["mac"]] <- shp[["mac"]][-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively
shp[["alb"]] <- shp[["alb"]][-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively

# add 50km buffer around ranges

shp <- lapply(shp,function(x) {
  st_union(st_buffer(x,50000))
})

# get binary of whether GBIF occurrence points fall inside buffered ranges

gbif_bin <- sapply(spp,function(x) {
  as.vector(st_intersects(st_transform(st_as_sf(gbif_occ[[x]],coords=c("long","lat"),crs=4326),crs=proj),shp[[x]],sparse=F))
},USE.NAMES=T,simplify=F)

# extract points outside and inside buffered range, respectively, and write to CSV file

gbif_ex <- sapply(spp,function(x) {
  gbif_occ[[x]][-which(gbif_bin[[x]]),]
},USE.NAMES=T,simplify=F)
gbif_in <- sapply(spp,function(x) {
  gbif_occ[[x]][which(gbif_bin[[x]]),]
},USE.NAMES=T,simplify=F)

gbif_ex.df <- do.call(rbind,gbif_ex)
rownames(gbif_ex.df) <- NULL
gbif_in.df <- do.call(rbind,gbif_in)
rownames(gbif_in.df) <- NULL

write.csv(gbif_ex.df,"gbif_exclude.csv",row.names=F,quote=F)
write.csv(gbif_in.df,"gbif_include.csv",row.names=F,quote=F)

# merge GBIF in buffered range with FIA occurrences, convert to df, write to CSV

occ_data <- sapply(spp,function(x) {
  rbind(occ_fia[[x]],gbif_in[[x]])
},USE.NAMES=T,simplify=F)

# remove conspicuously out-of-range FIA points

occ_data[["ste"]] <- subset(occ_data[["ste"]],!(lat > 42 | long < -85 & lat > 40.75))
occ_data[["bic"]] <- subset(occ_data[["bic"]],!(lat < 35 | lat > 48 | long > -78 & lat < 38))
occ_data[["alb"]] <- subset(occ_data[["alb"]],!(lat > 48 | long < -97))
occ_data[["mue"]] <- subset(occ_data[["mue"]],!(long > -84 & lat < 31))

occ.df <- do.call(rbind,occ_data)
rownames(occ.df) <- NULL
write.csv(occ.df,"occ_include.csv",row.names=F,quote=F)

# METHOD 2 (WIP): apply buffered range filter to all data

# # merge deduplicated GBIF occurrences with FIA occurrences, replace species code
#
# occ_data <- sapply(spp,function(x) {
#   rbind(occ_data[[x]],unique(gbif[gbif$spp == x,]))[,-3]
# },USE.NAMES=T,simplify=F)
# occ_data <- sapply(spp,function(x) {
#   occ_data[[x]]$sppcode <- rep(x,nrow(occ_data[[x]]));
#   return(occ_data[[x]])
# },USE.NAMES=T,simplify=F)
#
# # read in range shapefiles
#
# proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
# shp <- lapply(spp,function(x) {
#   st_cast(st_set_crs(st_read(paste0(args[3],"/",spp,"_little.shp")),crs=proj),"MULTIPOLYGON")
# })
#
# shp <- sapply(spp,function(x) {
#   st_cast(st_set_crs(st_read(paste0("~/qmac/analyses/envdata/rangedata","/",x,"_little.shp")),proj),"MULTIPOLYGON")
# },USE.NAMES=T,simplify=F)
# shp[["bic"]] <- shp[["bic"]][-12,] # removes Lake Winnebago
# shp[["mac"]] <- shp[["mac"]][-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively
# shp[["alb"]] <- shp[["alb"]][-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively
#
# # add 50km buffer around ranges
#
# shp <- lapply(shp,function(x) {
#   st_union(st_buffer(x,50000))
# })
#
# # get binary of whether points fall inside buffered ranges
#
# occ_data_bin <- sapply(spp,function(x) {
#   as.vector(st_intersects(st_transform(st_as_sf(occ_data[[x]],coords=c("long","lat"),crs=4326),crs=proj),shp[[x]],sparse=F))
# },USE.NAMES=T,simplify=F)
#
# # get points outside and inside of buffered range, respectively
#
# occ_ex <- sapply(spp,function(x) {
#   occ_data[[x]][-which(occ_data_bin[[x]]),]
# },USE.NAMES=T,simplify=F)
# occ_in <- sapply(spp,function(x) {
#   occ_data[[x]][which(occ_data_bin[[x]]),]
# },USE.NAMES=T,simplify=F)
#
# # convert to dataframes with species as factor variable
#
# occ_in.df <- do.call(rbind,occ_in)
# rownames(occ_in.df) <- NULL
# occ_ex.df <- do.call(rbind,occ_ex)
# rownames(occ_ex.df) <- NULL
#
# # write included/excluded points to respective csv files
#
# write.csv(occ_in.df,"occ_include.csv",row.names=F,quote=F)
# write.csv(occ_ex.df,"occ_exclude.csv",row.names=F,quote=F)
