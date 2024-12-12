library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidyverse)

spp <- c("mac", "alb", "bic", "ste", "mue")

# read in Maxent projections for each species

habsuitAUC.ls <- sapply(spp, function(x) {
    terra::rast(paste0("habsuitAUC_",x,".tif")) %>%
    dplyr::rename(habsuit = 1)
}, simplify = FALSE)
habsuitAICc.ls <- sapply(spp, function(x) {
    terra::rast(paste0("habsuitAICc_",x,".tif")) %>%
    dplyr::rename(habsuit = 1)
}, simplify = FALSE)

# read in map of US/CA/MX states

statemap <- ne_states(c("United States of America", "Canada", "Mexico"), returnclass = "sf") %>%
    st_crop(c(xmin = -130, xmax = -45, ymin = 20, ymax = 60)) %>%
    st_transform(crs = 4326)

# generate plots

plotAUC.ls <- sapply(spp, function(x) {
    ggplot() +
        geom_sf(data = statemap) +
        coord_sf(
            crs = 4326,
            xlim = c(-110, -60),
            ylim = c(25, 55),
            expand = FALSE) +
        geom_raster(
            data = habsuitAUC.ls[[x]], 
            aes(x = x, y = y, fill = habsuit), 
            alpha = 0.9) +
        xlab(NULL) +
        ylab(NULL) +
        scale_fill_gradientn(
            na.value = "transparent",
            colors = rev(grDevices::terrain.colors(50)), 
            name = "Habitat\nSuitability") +
        ggtitle(paste0(toupper(x)," - AUC")) +
        theme_minimal()
}, simplify = FALSE)
plotAICc.ls <- sapply(spp, function(x) {
    ggplot() +
        geom_sf(data = statemap) +
        coord_sf(
            crs = 4326,
            xlim = c(-110, -60), 
            ylim = c(25, 55),
            expand = FALSE) +
        geom_raster(
            data = habsuitAICc.ls[[x]], 
            aes(x = x, y = y, fill = habsuit), 
            alpha = 0.9) +
        xlab(NULL) +
        ylab(NULL) +
        scale_fill_gradientn(
            na.value = "transparent",
            colors = rev(grDevices::terrain.colors(50)),
            name = "Habitat\nSuitability") +
        ggtitle(paste0(toupper(x)," - AICc")) +
        theme_minimal()
}, simplify = FALSE)

# print maps to pdf

pdf("maxent_habsuitENAwhiteoaks.pdf")
lapply(spp, function(x) {
    print(plotAUC.ls[[x]])
    print(plotAICc.ls[[x]])
})
dev.off()