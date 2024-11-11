args <- commandArgs(trailingOnly=T)

# args <- c("occ_formaxent.csv","../coords_noSHF.csv")

set.seed(1234)

library(spThin)
library(sp)
library(raster)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(maxnet)
library(tidverse)
library(dismo)
# library(dismo)

# read in occurrence data and create list of species-specific occurrence dataframes

spp <- c("mac","alb","bic","ste","mue")
occ_all <- read.csv(args[1],header=T)
occ_byspp <- sapply(spp,function(x) {
  occ_all[occ_all$spp == x,]
},USE.NAMES=T,simplify=F)

# pre-thin occurrence records to 20000 for Q. alba - without thinning, the number of records requires too much memory for spThin to run

occ_byspp[["alb"]] <- occ_byspp[["alb"]][sample(nrow(occ_byspp[["alb"]]),20000),]

# thin occurrence points to 20km apart using spThin algorithm - run 200 replicates

thin.ls <- lapply(occ_byspp,function(x) {
  thin(loc.data=x,
    lat.col="lat",long.col="long",
    spec.col="spp",thin.par=20,reps=200,
    locs.thinned.list.return=T,
    write.files=F,write.log.file=T,
    log.file="occ_thinning_log.txt",
    verbose=T)
})

# get thinned dataframe with highest number of records for each species

occ_thinned <- lapply(thin.ls,function(x) {
  x[[1]]
})
occ_thinned <- sapply(spp,function(x) {
  occ_thinned[[x]]$spp <- x
  names(occ_thinned[[x]]) <- c("long","lat","spp")
  return(occ_thinned[[x]])
},USE.NAMES=T,simplify=F)

# convert list of species-specific dataframes into a single dataframe and write to CSV file

thinned.df <- do.call("rbind",occ_thinned)
rownames(thinned.df) <- NULL
write.csv(thinned.df,"occ_thinned_allspp.csv",row.names=F,quote=F)

# read in environmental rasters and clip to study extent (200 km buffer around Little range edges)

env_vars <- c(paste0("bio_",1:19),"elev")
rast.ls <- sapply(env_vars,function(x) {
  rast(paste0("~/qmac/analyses/envdata/climdata/wc2.1_30s_",x,".tif"))
},USE.NAMES=T,simplify=F)

# stack climate/elevation rasters, transform CRS to aea for NA (NAME)

stackrast <- rast(rast.ls)
stackrast <- crop(stackrast,ext(-110,-60,25,55)) # pre-crop to make reprojection faster
proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
stackrast <- project(stackrast,proj)

# read in map for US/Canada

namap <- ne_states(c("United States of America","Canada","Mexico"),returnclass="sf")
namap <- st_crop(namap,c(xmin=-110,xmax=-60,ymin=20,ymax=55))
namap <- project(vect(namap),proj)

# get polygons of lakes to remove from rasters

lakes <- ne_download(scale=10,type="lakes",category="physical",returnclass="sf")
lake_in <- c(grep("Lake Winnipeg|Lake Winnebago|Lake Manitoba|Lake of the Woods",lakes$name),grep("Great Lakes",lakes$name_alt))
lakes <- lakes[lake_in,] # retain only Lake Winnebago, Great Lakes, Lake Manitoba, Lake Winnipeg, Lake of the Woods
lakes <- project(vect(lakes),proj)

# mask freshwater features from rasters

stackrast <- mask(stackrast,lakes,inverse=T)

# get environmental data for each presence point - first convert dataframe to terra spatvector object

prspoints <- lapply(occ_thinned,function(x) {
  project(vect(x,geom=c("long","lat"),crs="epsg:4326"),proj)
})

# extract data from rasters at points

prsdata <- lapply(prspoints,function(x) {
  na.omit(extract(stackrast,x))
})

# define study regions by 200km buffered Little ranges

# read in Little ranges

range.ls <- sapply(spp,function(x) {
  st_cast(st_read(paste0("~/qmac/analyses/envdata/rangedata/",x,"_little.shp"),crs=proj),"MULTIPOLYGON")
},USE.NAMES=T,simplify=F)
range.ls[["bic"]] <- range.ls[["bic"]][-12,] # removes Lake Winnebago
range.ls[["mac"]] <- range.ls[["mac"]][-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively
range.ls[["alb"]] <- range.ls[["alb"]][-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively
range.ls <- lapply(range.ls,function(x) {
  st_union(st_buffer(x,200000))
})
range.ls <- lapply(range.ls,vect)

# crop stacked raster to species-specific study regions

stack_byspp <- lapply(range.ls,function(x) {
  mask(stackrast,x)
})

# randomly sample 10000 background points from cropped rasters

# sample points

bgpoints <- lapply(stack_byspp,function(x) {
  setNames(as.data.frame(randomPoints(raster(x),10000)),c("long","lat"))
})

# convert into terra spatialvector object and set CRS

bgpoints <- lapply(bgpoints,function(x) {
  vect(x,geom=c("long","lat"),crs=proj)
})

# convert background points to longlat and export CSV as single dataframe with species codes

bgpts_save <- lapply(bgpoints,function(x) {
  setNames(as.data.frame(project(x,"epsg:4326"),geom="XY"),c("long","lat"))
})
bgpts_save <- sapply(spp,function(x) {
  bgpts_save[[x]]$spp <- x
  return(bgpts_save[[x]])
},USE.NAMES=T,simplify=F)
bgpoints.df <- do.call("rbind",bgpts_save)
rownames(bgpoints.df) <- NULL
write.csv(bgpoints.df,"bgpoints.csv",row.names=F,quote=F)

# extract raster values for points

bgdata <- lapply(bgpoints,function(x) {
  extract(stackrast,x)
})

# add presence/background code and combine species-specific data for maxnet call

bgdata <- lapply(bgdata,function(x) {
  x$ID <- 0
  colnames(x)[1] <- "presence"
  return(x)
})
prsdata <- lapply(prsdata,function(x) {
  x$ID <- 1
  colnames(x)[1] <- "presence"
  return(x)
})
SDMdata <- sapply(spp,function(x) {
  rbind(prsdata[[x]],bgdata[[x]])
},USE.NAMES=T,simplify=F)

# run Maxent for each species (implemented in maxnet package)

trainmodels <- lapply(SDMdata,function(x) {
  maxnet(x$presence,x[,-1])
})

# predict occurrence (cloglog link) across ENA

predmodels <- lapply(trainmodels,function(x) {
  predict(stackrast,x,type="cloglog",na.rm=T)
})

# write raster files of habitat suitability for each species across ENA

sapply(spp,function(x) {
  writeRaster(predmodels[[x]],paste0(x,"_maxent_cloglog.tif"),overwrite=T,gdal=c("COMPRESS_LZW"))
})

crds <- read.csv(args[2],header=T)
crds <- project(vect(crds,geom=c("long","lat"),crs="epsg:4326"),proj)


# bg.ls <- lapply(range.ls,function(x) {
#   spatSample(x,10000)
# })
# bgdata <- lapply(bg.ls,function(x) {
#   na.omit(extract(stackrast,x))
# })
# bgdata <- lapply(bg.ls,function(x) {
#   extract(stackrast,x)
# })
#
# sapply(spp,function(x) {
#   writeVector(range.ls[[x]],paste0(x,"_studyregion.shp"))
# })

# randomly sample 10000 ba

# # use all occurrence records (across species) for 100 km buffered convex hull defining study extent
#
# occ_hull <- buffer(project(convHull(vect(occ_thinned,geom=c("long","lat"),crs="epsg:4326")),proj),100000)
# stackrast <- mask(stackrast,stack_byspp)







# read in Little range shapefiles

# proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
# range.ls <- sapply(spp,function(x) {
#   st_read(paste0("~/qmac/analyses/envdata/rangedata/",x,"_little.shp"))
# })

# generate background points for Maxent - first read in CSV files
