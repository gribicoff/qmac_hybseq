args <- commandArgs(trailingOnly=T)

# args <- "../coords_noSHF.csv"

set.seed(1234)

library(spThin)
library(sp)
library(raster)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidyverse)

spp <- c("mac","alb","bic","ste","mue")

# import cleaned occurrence data from GBIF - CSVs and shapefiles

occ.shp <- sapply(spp,function(x) {
    st_read(paste0("cleanedgbifoccs_",x,".shp"))
},simplify = FALSE)
occ.ls <- sapply(spp,function(x) {
    read.csv(file = paste0("cleanedgbifoccs_",x,".csv"))
},simplify = FALSE)

occ.ls <- sapply(spp,function(x) {
    occ.ls[[x]] %>% dplyr::mutate(species = x)
},simplify = FALSE)

# thin records using spThin - 20km and 50km thinned datasets

# 20km

thin20.ls <- lapply(occ.ls,function(x) {
    x %>%
    thin(lat.col = "lat",
        long.col = "long",
        spec.col = "species",
        thin.par=20,
        reps=200,
        locs.thinned.list.return = TRUE,
        write.files = FALSE,
        write.log.file = TRUE,
        log.file = "occ_thinning_20km.log",
        verbose = TRUE)
})

# 50 km

# thin50.ls <- lapply(occ.ls,function(x) {
#     x %>%
#     thin(lat.col = "lat",
#         long.col = "long",
#         spec.col = "species",
#         thin.par=50,
#         reps=200,
#         locs.thinned.list.return = TRUE,
#         write.files = FALSE,
#         write.log.file = TRUE,
#         log.file = "occ_thinning_50km.log",
#         verbose = TRUE)
# })

# select thinned dataset with greatest number of retained points per species

thin20.ls <- lapply(thin20.ls,function(x) {
   x[[1]] %>% dplyr::rename(long = Longitude, lat = Latitude)
})
# thin50.ls <- lapply(thin50.ls,function(x) {
#    x[[1]] %>% dplyr::rename(long = Longitude, lat = Latitude)
# })

# write thinned datasets to CSV files

sapply(spp,function(x) {
  thin20.ls[[x]] %>%
  write.csv(file = paste0("thinnedgbifoccs_",x,"_20km.csv"),row.names = FALSE, quote = FALSE)
  # thin50.ls[[x]] %>%
  # dplyr::rename(long = Longitude, lat = Latitude) %>%
  # write.csv(file = paste0("thinnedgbifoccs_",x,"_50km.csv"),row.names = FALSE, quote = FALSE)
})

# plot thinned occurrence data

uscan <- ne_states(c("United States of America","Canada"),returnclass="sf")
uscan <- st_crop(uscan,c(xmin=-110,xmax=-60,ymin=20,ymax=55))
uscan <- st_transform(uscan,crs=4326)

sapply(spp,function(x) {
  p <- ggplot() +
    geom_sf(data=uscan) +
    geom_point(data = thin20.ls[[x]],aes(x=long,y=lat),size=0.1) +
    ggtitle(paste0("Q",toupper(x)," 20km Thinning"))
  p %>% ggsave(filename=paste0("thinnedgbifoccs_",x,"_20km.png"))
  # p <- ggplot() +
  #   geom_sf(data=uscan) +
  #   geom_point(data = thin50.ls[[x]],aes(x=long,y=lat),size=0.1) +
  #   ggtitle(paste0("Q",toupper(x)," 50km Thinning"))
  # p %>% ggsave(filename=paste0("thinnedgbifoccs_",x,"_50km.png"))
})

# read in individual bioclim rasters and stack into single raster, transform CRS to NA aea (???)

env_vars <- paste0("bio_",1:19)
rast.ls <- sapply(env_vars,function(x) {
  rast(paste0("~/qmac/analyses/envdata/climdata/wc2.1_30s_",x,".tif"))
},simplify = FALSE)
stackrast <- rast(rast.ls)
stackrast <- crop(stackrast,ext(-110,-60,25,55)) # pre-crop to make reprojection faster
proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
stackrast <- project(stackrast,proj)

# get polygons of lakes to remove from rasters

lakes <- ne_download(scale=10,type="lakes",category="physical",returnclass="sf")
lake_in <- c(grep("Lake Winnipeg|Lake Winnebago|Lake Manitoba|Lake of the Woods",lakes$name),grep("Great Lakes",lakes$name_alt))
lakes <- lakes[lake_in,] # retain only Lake Winnebago, Great Lakes, Lake Manitoba, Lake Winnipeg, Lake of the Woods
lakes <- project(vect(lakes),proj)

# mask freshwater features from rasters

stackrast <- mask(stackrast,lakes,inverse = TRUE)

# define study regions by 200km convex hull around occurrence points for each species

hulls <- lapply(occ.shp,function(x) {
  x %>% st_cast("MULTIPOINT") %>%
  st_transform(crs=proj) %>%
  st_convex_hull %>%
  st_buffer(200000) %>%
  st_union %>%
  vect
})

# crop stacked rasters to hull extents

stack.ls <- lapply(hulls,function(x) mask(stackrast,x))

# transform CRS of presence points

prspts20 <- lapply(thin20.ls,function(x) {
  vect(x,geom=c("long","lat"),crs="epsg:4326") %>%
  project(proj) %>%
  terra::crds(df=TRUE) %>%
  dplyr::rename(X = x, Y = y)
})

# prspts50 <- lapply(thin50.ls,function(x) {
#   vect(x,geom=c("long","lat"),crs="epsg:4326") %>%
#   project(proj) %>%
#   terra::crds(df=TRUE) %>%
#   dplyr::rename(X = x, Y = y)
# })

# get 25000 background points from hulls for each species

bgpts <- lapply(stack.ls,function(x) {
  raster(x) %>%
  randomPoints(25000) %>%
  as.data.frame %>%
  dplyr::rename(X = 1, Y = 2)
})

# save ENMevaluate inputs: terra SpatRasters as tif files, presence and background points as CSV files (for each species)

lapply(spp,function(x) {
  writeRaster(stack.ls[[x]],
              paste0("stackraster_",x,".tif"),
              overwrite=TRUE)
  write.csv(prspts20[[x]],
            paste0("maxentprspts20_",x,".csv"),
            quote = FALSE,
            row.names = FALSE)
  write.csv(bgpts[[x]],
            paste0("maxentbgpts_",x,".csv"),
            quote = FALSE,
            row.names = FALSE)
})

# partition presence/background points into testing and training groups using checkerboard2 method in ENMeval

# checks <- sapply(spp,function(x) {
#   get.block(prspts20[[x]],bgpts[[x]],orientation = "lat_lon")
# },simplify = FALSE)

# # plot partitions

# lapply(spp,function(x) {
#   png(paste0("partitions_",x,".png"))
#   print(evalplot.grps(envs = raster(stack.ls[[x]]),
#                       pts = bgpts[[x]],
#                       pts.grp = checks[[x]]$bg.grp,
#                       pts.size = 0.25))
#   dev.off()
# })

# specify Maxent feature classes to be tested

ftclasses <- c("L","Q","H","LQ","QH","LH","LQH","LQHP","LQHPT")
# fts <- c("L","Q","H","P","T")
# ftcombn <- function(n) {
#   combn(fts,n,function(x) paste0(x,collapse=""))}
# ftclasses <- sapply(1:5,ftcombn) %>% unlist

# run ENMevaluate for each 

# library(doSNOW)
# max20 <- sapply(spp,function(x) {
#   ENMevaluate(occs = prspts20[[x]],
#               envs = stack.ls[[x]],
#               bg = bgpts[[x]],
#               partitions = "block",
#               algorithm = "maxent.jar",
#               tune.args = list(
#                 fc = ftclasses,
#                 rm = seq(0.5,5,0.5)),
#               quiet = FALSE)
# },simplify = FALSE)

max20_mac <- ENMevaluate(occs = prspts20[["mac"]],
              envs = stack.ls[["mac"]],
              bg = bgpts[["mac"]],
              partitions = "block",
              algorithm = "maxent.jar",
              tune.args = list(
                fc = ftclasses,
                rm = seq(0.5,5,0.5)),
              parallel = TRUE,
              numCores = 32,
              quiet = FALSE)

max20_alb <- ENMevaluate(occs = prspts20[["alb"]],
              envs = stack.ls[["alb"]],
              bg = bgpts[["alb"]],
              partitions = "block",
              algorithm = "maxent.jar",
              tune.args = list(
                fc = ftclasses,
                rm = seq(0.5,5,0.5)),
              parallel = TRUE,
              numCores = 32,
              quiet = FALSE)

max20_bic <- ENMevaluate(occs = prspts20[["bic"]],
              envs = stack.ls[["bic"]],
              bg = bgpts[["bic"]],
              partitions = "block",
              algorithm = "maxent.jar",
              tune.args = list(
                fc = ftclasses,
                rm = seq(0.5,5,0.5)),
              parallel = TRUE,
              numCores = 32,
              quiet = FALSE)

max20_ste <- ENMevaluate(occs = prspts20[["ste"]],
              envs = stack.ls[["ste"]],
              bg = bgpts[["ste"]],
              partitions = "block",
              algorithm = "maxent.jar",
              tune.args = list(
                fc = ftclasses,
                rm = seq(0.5,5,0.5)),
              parallel = TRUE,
              numCores = 32,
              quiet = FALSE)

max20_mue <- ENMevaluate(occs = prspts20[["mue"]],
              envs = stack.ls[["mue"]],
              bg = bgpts[["mue"]],
              partitions = "block",
              algorithm = "maxent.jar",
              tune.args = list(
                fc = ftclasses,
                rm = seq(0.5,5,0.5)),
              parallel = TRUE,
              numCores = 32,
              quiet = FALSE)

max20 <- sapply(spp,function(x) {
  ENMevaluate(occs = prspts20[[x]],
              envs = stack.ls[[x]],
              bg = bgpts[[x]],
              partitions = "block",
              algorithm = "maxent.jar",
              tune.args = list(
                fc = ftclasses,
                rm = seq(0.5,5,0.5)),
              parallel = TRUE,
              numCores = 32,
              quiet = FALSE)
},simplify = FALSE)











# extract data from rasters at points

# prsdata20 <- lapply(prspts20,function(x) {
#   na.omit(extract(stackrast,x))
# })
# prsdata50 <- lapply(prspts50,function(x) {
#   na.omit(extract(stackrast,x))
# })

# sample points

bgpts <- lapply(stack.ls,function(x) {
  raster(x) %>%
  randomPoints(25000) %>%
  as.data.frame %>%
  dplyr::rename(X = 1, Y = 2) %>%
  vect(geom = c("long","lat"),crs=proj)
})

# convert spatvectors of background points to dataframes

bgpts.df <- lapply(bgpts,function(x) {
  terra::crds(x,df = TRUE) %>%
  dplyr::rename(X = x, Y = y)
})

# convert background points to longlat and export CSV files/maps for each species

bgpts.df <- lapply(bgpts,function(x) {
  x %>% terra::project("epsg:4326") %>%
  terra::as.data.frame(geom = "XY") %>%
  dplyr::rename(long = 1, lat = 2)
})

sapply(spp,function(x) {
  bgpts.df[[x]] %>%
  write.csv(file = paste0("maxentbgpts_",x,".csv"),quote = FALSE,row.names = FALSE)
  p <- ggplot() +
    geom_sf(data=uscan) +
    geom_point(data = bgpts.df[[x]],aes(x=long,y=lat),size=0.1) +
    ggtitle(paste0("Q",toupper(x)," Background Points"))
  p %>% ggsave(filename=paste0("maxentbgpts_",x,".png"))
})

# extract environmental raster values for background points

bgdata <- lapply(bgpts,function(x) terra::extract(stackrast,x))

# add presence/background binary codes for presence and background points

bgdata <- lapply(bgdata,function(x) {
  x %>% mutate(ID = 0) %>% dplyr::rename(presence = ID)
})
prsdata20 <- lapply(prsdata20,function(x) {
  x %>% mutate(ID = 1) %>% dplyr::rename(presence = ID)
})
# prsdata50 <- lapply(prsdata50,function(x) {
#   x %>% mutate(ID = 1) %>% dplyr::rename(presence = ID)
# })

# combine presence/background dataframes

SDMdata20 <- sapply(spp,function(x) {
  rbind(prsdata20[[x]],bgdata[[x]])
},simplify = FALSE)
# SDMdata50 <- sapply(spp,function(x) {
#   rbind(prsdata50[[x]],bgdata[[x]])
# },simplify = FALSE)

# partition data into groups for testing and training with 10-fold CV - see ENMeval vignette

prspts20.df <- lapply(prspts20,function(x) {
  x %>% terra::project("epsg:4326") %>%
  terra::as.data.frame(geom = "XY") %>%
  dplyr::rename(long = 1, lat = 2)
})

stacksamp <- lapply(stack.ls,function(x) {
  terra::project(x,"epsg:4326")
})

CVgrp <- sapply(spp,function(x) {
  get.checkerboard2(occs = prspts20.df[[x]],bg = bgpts.df[[x]],kfolds = 10)
},simplify = FALSE)

# run Maxent for each species (maxnet implementation)

trainmod20 <- lapply(SDMdata20,function(x) maxnet(x$presence,x[,-1]))
trainmod50 <- lapply(SDMdata50,function(x) maxnet(x$presence,x[,-1]))

# predict occurrence (cloglog link) across eastern North America

predmod20 <- lapply(trainmod20,function(x) {
  predict(stackrast,x,type="cloglog",na.rm = TRUE)
})
predmod50 <- lapply(trainmod50,function(x) {
  predict(stackrast,x,type="cloglog",na.rm = TRUE)
})

# write raster files of habitat suitability for each species across ENA

sapply(spp,function(x) {
  writeRaster(predmod20[[x]],paste0(x,"_maxent20km_cloglog.tif"),overwrite=T,gdal=c("COMPRESS_LZW"))
  writeRaster(predmod50[[x]],paste0(x,"_maxent50km_cloglog.tif"),overwrite=T,gdal=c("COMPRESS_LZW"))
})

# plot suitability maps

pdf("maxent_suitability_cloglog.pdf")
sapply(spp,function(x) {
  plot(predmod20[[x]],main = paste0("Q",toupper(x)," 20km"))
  plot(predmod50[[x]],main = paste0("Q",toupper(x)," 50km"))
},simplify = FALSE)
dev.off()

# stack SDMs of all species into single raster

predstack <- sapply(spp,function(x) {
  rast(paste0(x,"_maxent_cloglog.tif"))
},simplify = FALSE)
predstack <- rast(predstack)

# import coordinate data, convert to terra spatvector, convert CRS

# crds_proj <- read.csv()

crds <- read.csv(args) %>% mutate(species = gsub(".*_","",ind_ID))
crds_proj <- crds %>% vect(geom = c("long","lat"),crs = "epsg:4326") %>% project(proj)

# extract suitability values (cloglog) for each sample coordinate

crdsdata <- extract(predstack,crds_proj,layer = crds_proj$species)
extract(predstack,crds_proj)
# combine habitat suitability with original coordinates and write dataframe to CSV

