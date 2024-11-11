args <- commandArgs(trailingOnly = TRUE)

# args <- "sppcode"

library(raster)
library(terra)
library(ENMeval)
library(tidyverse)

# read in ENMeval object

max <- readRDS(paste0("ENMeval_",args[1],".rds"))

# subset evaluation results across models

max.res <- eval.results(max)

# select optimal model based on test omission rate/AUC calculations

opt.AUC <- max.res %>% 
    dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
    dplyr::filter(auc.val.avg == max(auc.val.avg))

# select optimal model based on (delta) AICc

opt.AICc <- max.res %>% 
    dplyr::filter(delta.AICc == 0)

# get model details

mod.AUC <- eval.models(max)[[opt.AUC$tune.args]]
mod.AICc <- eval.models(max)[[opt.AICc$tune.args]]

# get prediction rasters for selected models, transform CRS, export as tif

pred.AUC <- eval.predictions(max)[[opt.AUC$tune.args]] %>%
    raster::projectRaster(crs = 4326) %>%
    rast %>%
    terra::project("epsg:4326") %>%
    terra::writeRaster(paste0("habsuitAUC_",args[1],".tif"), overwrite = TRUE)
pred.AICc <- eval.predictions(max)[[opt.AICc$tune.args]] %>%
    raster::projectRaster(crs = 4326) %>%
    rast %>%
    terra::project("epsg:4326") %>%
    terra::writeRaster(paste0("habsuitAICc_",args[1],".tif"), overwrite = TRUE)