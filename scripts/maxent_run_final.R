options(java.parameters = "-Xmx32g")
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)

# args <- "sppname"

library(rJava)
library(terra)
library(ENMeval)

# read in presence points, background points, stacked raster - convert CRS of coordinates

proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
stackrast <- rast("stackraster_allspp.tif")
bgpts <- read.csv(paste0("bgpoints_",args[1],".csv"), header = TRUE) %>%
    vect(geom = c("long", "lat"), crs = "epsg:4326") %>%
    project(proj) %>%
    as.data.frame(geom = "XY") %>%
    dplyr::rename(long = x, lat = y)
prspts <- read.csv(paste0("thinnedoccs_",args[1],".csv"), header = TRUE) %>%
    vect(geom = c("long", "lat"), crs = "epsg:4326") %>%
    project(proj) %>%
    as.data.frame(geom = "XY") %>%
    dplyr::rename(long = x, lat = y)

# specify Maxent feature classes and RMs to be tested

ftclasses <- c("H", "LQH", "LQHP")
lambdas <- c(1, 2, 5, 10, 20)

# run ENMeval

max <- ENMevaluate(
    occs = prspts,
    envs = stackrast,
    bg = bgpts,
    partitions = "checkerboard1",
    algorithm = "maxent.jar",
    tune.args = list(
        fc = ftclasses,
        rm = lambdas),
    parallel = TRUE,
    numCores = 32,
    quiet = FALSE)

# save ENMeval object

saveRDS(max, paste0("ENMeval_",args[1],".rds"))