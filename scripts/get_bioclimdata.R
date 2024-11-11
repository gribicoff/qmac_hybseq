args <- commandArgs(trailingOnly=TRUE)

# args <- "../coords_noSHF.csv"

library(terra)
library(tidyverse)

# read in coordinates

crds <- read.csv(args[1],header=T)
crds <- vect(crds,geom=c("long","lat"),crs="epsg:4326")

# read in, crop bioclim/elevation/PET rasters - PET from https://www.nature.com/articles/s41597-022-01493-1

biovars <- c(paste0("bio_",1:19),"elev")
biostack <- rast(sapply(biovars,function(x) {
  rast(paste0("wc2.1_30s_",x,".tif"))
},simplify = FALSE))
biostack <- c(biostack,rast("et0_v3_yr.tif"))
names(biostack)[[length(names(biostack))]] <- "PET"
names(biostack) <- gsub("_","",names(biostack))

# extract raster values for coordinates

biodata <- terra::extract(biostack,crds)
biodata[,1] <- crds$ind_ID

# calculate moisture index I_m as 100*(MAP - PET)/PET

biodata <- biodata %>%
  mutate(moist_index = 100*((bio12-PET)/PET)) %>%
  select(-PET)
names(biodata)[1] <- "ind_ID"

# write csv file

write.csv(biodata,file="../coords_noSHF_bioclim.csv",quote=F,row.names=F)