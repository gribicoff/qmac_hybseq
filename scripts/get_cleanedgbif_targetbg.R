library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)

spp <- c("mac","alb","bic","ste","mue")

# read in occurrence data cleaned manually in QGIS, convert shapefiles into dataframes of coordinates

occ.shp <- sapply(spp,function(x) {
    st_read(paste0("cleanedgbifoccs_",x,".shp"))
},simplify = FALSE)
occ.ls <- lapply(occ.shp,function(x) {
    x %>% st_coordinates %>% as.data.frame %>% dplyr::rename(long = 1, lat = 2)
})

# plot points for each species

uscan <- ne_states(c("United States of America","Canada"),returnclass="sf")
uscan <- st_crop(uscan,c(xmin=-140,xmax=-55,ymin=20,ymax=80))
uscan <- st_transform(uscan,crs=4326)

sapply(spp,function(x) {
  p <- ggplot() +
    geom_sf(data=uscan) +
    geom_point(data = occ.ls[[x]],aes(x=long,y=lat),size=0.1) +
    ggtitle(paste0("Q",toupper(x)))
  p %>% ggsave(filename=paste0("cleanedgbifoccs_",x,".png"))
})

# write CSV files of occurrence data

sapply(spp,function(x) {
    write.csv(occ.ls[[x]],file=paste0("cleanedgbifoccs_",x,".csv"),row.names = FALSE,quote = FALSE)
})
