args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/qmac/analyses/new_analyses/admix_formodeling_alldata.csv", "/path/to/outputdir/"

library(tidyverse)
library(scatterpie)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)

# read in coord/admix df, shapefile as sf objects, assign proj data/convert from polygon to multipolygon

spp <- c("mac", "alb", "bic", "ste", "mue")
adcrds <- read.csv(args[1], header = TRUE) %>%
  mutate(mac_frac = mac1_frac + mac2_frac, .keep = "unused") %>%
  dplyr::select(
    ind_ID,
    sppcode,
    contains("frac"),
    lat,
    long,
    site)
shp.ls <- sapply(spp, function(x) {
  st_read(paste0("~/qmac/analyses/envdata/rangedata/", x, "_little.shp"), crs = "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu") %>%
  st_cast("MULTIPOLYGON")
}, simplify = FALSE)
shp.ls[["bic"]] <- shp.ls[["bic"]][-12,] # removes Lake Winnebago
shp.ls[["alb"]] <- shp.ls[["alb"]] [-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively
shp.ls[["mac"]] <- shp.ls[["mac"]][-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively

# subset by species, make new column with admixed fraction for mapping, group individuals by site

ad.ls <- sapply(spp, function(x) {
  adcrds %>%
    dplyr::filter(sppcode == x) %>%
    dplyr::mutate(admix_frac = 1 - .data[[paste0(x, "_frac")]]) %>%
    dplyr::select(., paste0(x, "_frac"), admix_frac, lat, long, site) %>%
    group_by(site) %>%
    summarise(across(everything(), mean), n_site = n())
}, simplify = FALSE)

# load/convert multipolygon of US+CAN, crop and set crs

uscan <- ne_states(c("United States of America", "Canada"), returnclass = "sf") %>%
  st_crop(c(xmin = -120, xmax = -60, ymin = 20, ymax = 60)) %>%
  st_transform(crs = 4326)

# mac range: "#1aa3ff", mac1_frac: "#0059b3", mac2_frac: "#4dffff"

# plot maps of individuals and admixture fractions

sppcol <- c("#1aa3ff", "#ee90f371", "#059c07", "#e2d219", "#d65900d2") # colors for species ranges
names(sppcol) <- spp
admixcol <- list("#046ed7", "#6f16f48f", "#06f10abb", "#fbff00", "#f99b03c9") # colors for spp-specific ancestry fractions
names(admixcol) <- spp

pdf(paste0(args[2], "admixmap_perspp.pdf"))
sapply(spp, function(x) {
  ggplot() +
    geom_sf(data = uscan) +
    geom_sf(data = shp.ls[[x]], alpha = 0.35, fill = sppcol[x]) +
    coord_sf(crs = 4236, xlim = c(-103, -66), ylim = c(28, 52)) +
    geom_scatterpie(data = ad.ls[[x]], aes(x = long, y = lat, r = 0.2*log2(n_site+1)), cols = c(paste0(x, "_frac"), "admix_frac"), linewidth = 0.1, alpha = 0.9) +
    scale_fill_manual(values = c(admixcol[[x]], "red"), labels = c(toupper(x), "ADMIX")) +
    labs(xlab = NULL, ylab = NULL, fill = "Cluster") +
    theme_minimal()
}, simplify = FALSE)
dev.off()