set.seed(1234)

library(sp)
library(raster)
library(terra)
library(sf)
library(dismo)
library(spThin)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidyverse)

spp <- c("mac", "alb", "bic", "ste", "mue")
proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"

# import manually cleaned occurrence data

occ.ls <- sapply(spp, function(x) {
    read.csv(paste0("cleanedoccs_",x,".csv"), header = TRUE) %>%
    dplyr::mutate(species = x)
}, simplify = FALSE)

# thin occurrence records to 20km on per-species basis

thin.ls <- sapply(spp, function(x) {
    occ.ls[[x]] %>%
        thin(lat.col = "lat",
            long.col = "long",
            spec.col = "species",
            thin.par = 20,
            reps = 200,
            locs.thinned.list.return = TRUE,
            write.files = FALSE,
            write.log.file = TRUE,
            log.file = paste0("occ_thinning_",x,".log"),
            verbose = TRUE)
}, simplify = FALSE)

# select thinned dataset with greatest number of retained points per species

thin.ls <- sapply(spp, function(x) {
    thin.ls[[x]][[1]] %>%
        dplyr::rename(long = Longitude, lat = Latitude)
}, simplify = FALSE)

# write thinned occurrence datasets to csv files

lapply(spp, function(x) {
    write.csv(thin.ls[[x]], file = paste0("thinnedoccs_", x, ".csv"), row.names = FALSE, quote = FALSE)
})

# read in individual bioclim rasters and stack into single raster, transform CRS to NAAEA

env_vars <- paste0("bio_", 1:19)
rast.ls <- sapply(env_vars, function(x) {
    rast(paste0("~/qmac/analyses/envdata/climdata/wc2.1_30s_", x, ".tif"))
},simplify = FALSE)
stackrast <- rast(rast.ls)
stackrast <- crop(stackrast, ext(-110,-60,25,55)) # pre-crop to make reprojection faster
stackrast <- project(stackrast, proj) # reproject

# get polygons of lakes to remove from rasters

lakes <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf")
lake_in <- c(grep("Lake Winnipeg|Lake Winnebago|Lake Manitoba|Lake of the Woods", lakes$name), grep("Great Lakes", lakes$name_alt))
lakes <- lakes[lake_in,] %>% # retain only Lake Winnebago, Great Lakes, Lake Manitoba, Lake Winnipeg, Lake of the Woods
    vect %>%
    project(proj)

# mask freshwater features from rasters

stackrast <- mask(stackrast, lakes, inverse = TRUE)

# export stacked raster

writeRaster(stackrast, "stackraster_allspp.tif", overwrite = TRUE)

# define study regions by 250km convex hull around occurrence points for each species

hulls <- sapply(spp, function(x) {
    occ.ls[[x]] %>% 
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_cast("MULTIPOINT") %>%
    st_transform(crs = proj) %>%
    st_convex_hull %>%
    st_buffer(250000) %>%
    st_union %>%
    vect
}, simplify = FALSE)

# mask freshwater features from hull extents

hulls <- lapply(hulls, function(x) terra::erase(x, lakes))

# crop stacked rasters to hull extents

stack.ls <- lapply(hulls, function(x) terra::mask(stackrast, x))

# generate background points - randomly sample 20,000 points across study area rasters (hulls)

bg.ls <- sapply(spp, function(x) {
    stack.ls[[x]] %>%
        spatSample(size = 20000, method = "random", replace = FALSE, na.rm = TRUE, xy = TRUE) %>%
        dplyr::mutate(long = x, lat = y, .keep = "none")
}, simplify = FALSE)

# write background points as csv files of coordinates

lapply(spp, function(x) {
    bg.ls[[x]] %>%
        vect(geom = c("long", "lat"), crs = proj) %>%
        project("epsg:4326") %>%
        as.data.frame(geom = "XY") %>%
        dplyr::rename(long = x, lat = y) %>%
        write.csv(file = paste0("bgpoints_",x,".csv"), row.names = FALSE, quote = FALSE)
})

## all inputs for Maxent/ENMeval are now available in directory - 
## per-species cleaned/thinned presence points (CRS must be transformed to NAAEA),
## per-species background points (CRS must be transformed to NAAEA),
## single stacked raster of bioclim variables across ALL of ENA (already transformed to NAAEA) 