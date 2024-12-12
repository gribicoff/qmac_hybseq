args <- commandArgs(trailingOnly = TRUE)

# args <- c("admixsitecoords.csv", "/path/to/outputdir/")

library(tidyverse)
library(scatterpie)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)

# read in coord/admix df, subset only mac, make new column with non-mac fraction for mapping

adcrds <- read.csv(args[1], header = TRUE) %>%
  dplyr::mutate(
    sppcode = as.factor(gsub(".*_", "", ind_ID)),
    nonmac_frac = 1 - (mac1_frac + mac2_frac)) %>%
    dplyr::filter(sppcode == "mac") %>%
    dplyr::select(mac1_frac, mac2_frac, nonmac_frac, lat, long, site)
adcrds$nonmac_frac[adcrds$nonmac_frac <= 1e-3] <- 0 # convert very small admix fractions to 0 (artifact from pophelper cluster alignment output)
adcrds$mac1_frac[adcrds$mac1_frac <= 1e-3] <- 0 # convert very small admix fractions to 0 (artifact from pophelper cluster alignment output)
adcrds$mac2_frac[adcrds$mac2_frac <= 1e-3] <- 0 # convert very small admix fractions to 0 (artifact from pophelper cluster alignment output)

# group individuals by site

adcrds <- adcrds %>%
  group_by(site) %>%
  summarise(across(everything(), mean), n_site = n())

# read in shapefile as sf objects, assign proj data/convert from polygon to multipolygon

mac_shp <- st_read("~/qmac/analyses/envdata/rangedata/mac_little.shp", crs = "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu") %>%
  st_cast("MULTIPOLYGON") %>%
  dplyr::slice(-c(2,3,18)) # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively

# load/convert multipolygon of US+CAN, crop and convert crs

uscan <- ne_states(c("United States of America", "Canada"), returnclass = "sf") %>%
  st_crop(c(xmin = -120, xmax = -60, ymin = 20, ymax = 60)) %>%
  st_transform(crs = 4326)

# mac range: "#1aa3ff", mac1_frac: "#0059b3", mac2_frac: "#4dffff"

# plot map of Q. macrocarpa population structure

pdf(paste0(args[2], "qmacstructmap.pdf"))
ggplot() +
  geom_sf(data = uscan) +
  geom_sf(data = mac_shp, alpha = 0.35, fill = "#1aa3ff") +
  coord_sf(crs = 4236, xlim = c(-103, -73.5), ylim = c(34.75, 50)) +
  geom_scatterpie(data = adcrds, aes(x = long, y = lat, r = 0.2*log2(n_site+1)), cols = c("mac1_frac", "mac2_frac", "nonmac_frac"), linewidth = 0.1, alpha = 0.75) +
  scale_fill_manual(values = c("#0059b3", "#4dffff", "red"), labels = c("MAC1", "MAC2", "NONMAC")) +
  labs(xlab = NULL, ylab = NULL, fill = "Cluster") +
  theme_minimal() +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size=7))
dev.off()