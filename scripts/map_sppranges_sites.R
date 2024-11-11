args <- commandArgs(trailingOnly = TRUE)

# args <- c("admixsitecoords.csv", "./figures/")

library(ggpattern)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidyverse)

spp <- c("mac", "alb", "bic", "ste", "mue")
proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"

# read in basemap

namap <- ne_states(c("United States of America", "Canada", "Mexico"), returnclass = "sf")
namap <- namap %>%
  st_crop(c(xmin = -130, xmax = -50, ymin = 10, ymax = 70))

# read in shapefiles

rlist <- sapply(spp, function(x) {
  st_read(paste0("~/qmac/analyses/envdata/rangedata/",x,"_little.shp"), crs = proj) %>%
  st_cast("MULTIPOLYGON")
}, simplify = FALSE)

rlist[["bic"]] <- rlist[["bic"]][-12,] # removes Lake Winnebago
rlist[["mac"]] <- rlist[["mac"]][-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively
rlist[["alb"]] <- rlist[["alb"]][-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively

# read in csv with sites/coordinates

sites.df <- read.csv(args[1], header = TRUE) %>%
  dplyr::select(site, lat, long) %>%
  group_by(site) %>%
  dplyr::summarise(lat = mean(lat), long = mean(long), num_inds = n())

# generate plot and write to pdf

pdf(paste0(args[2],"sites_sppranges_little.pdf"))
ggplot() +
  geom_sf(data = namap) +
  geom_sf(data = rlist[["alb"]], alpha = 0.4, fill = "#3aa0bf", linetype = "blank") +
  geom_sf(data = rlist[["bic"]], alpha = 0.4, fill = "#5900e2", linetype = "blank") +
  geom_sf(data = rlist[["mue"]], alpha = 0.4, fill = "#ff4000", linetype = "blank") +
  geom_sf(data = rlist[["ste"]], alpha = 0.4,fill = "#ffbf00", linetype = "blank") +
  geom_sf(data = rlist[["mac"]], alpha = 0, fill = "white", linetype = "solid", linewidth = 0.5) +
  geom_point(data = sites.df, aes(x = long, y = lat, size = 0.2*log2(num_inds + 1)), shape = 21, color = "#ff2222", stroke = 0.75) +
  coord_sf(crs = 4326, xlim = c(-105, -65), ylim = c(25, 51.5)) +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()
