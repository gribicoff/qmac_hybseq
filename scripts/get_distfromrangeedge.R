library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)

# read in shapefiles as sf objects, assign proj data/convert from polygon to multipolygon

spp <- c("mac", "bic", "alb", "ste", "mue")
proj <- "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"

shp.ls <- sapply(spp, function(x) {
  st_read(paste0("~/qmac/analyses/envdata/rangedata/",x, "_little.shp"), crs = proj) %>%
  st_cast("MULTIPOLYGON")
}, simplify = FALSE)

shp.ls[["bic"]] <- shp.ls[["bic"]][-12,] # removes Lake Winnebago
shp.ls[["mac"]] <- shp.ls[["mac"]][-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively
shp.ls[["alb"]] <- shp.ls[["alb"]][-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively

# get coastline, NA lakes from Natural Earth, buffer by 10km

coast_shp <- ne_download(
  scale = 10,
  type = "coastline",
  category = "physical",
  returnclass = "sf")
coast_shp <- coast_shp %>%
  st_make_valid %>%
  st_crop(ymin = 25, ymax = 50, xmin = -100, xmax = -65)

lakes_shp <- ne_download(
  scale = 10,
  type = "lakes",
  category = "physical",
  returnclass = "sf")
lakes_shp <- lakes_shp[c(9:12,18,20,28,466:468),] # retain only Lake Winnebago, Great Lakes, Lake Manitoba, Lake Winnipeg, Lake of the Woods

coast_shp <- coast_shp %>%
  st_transform(crs = proj) %>%
  st_make_valid %>%
  st_buffer(10000) %>%
  st_union

lakes_shp <- lakes_shp %>%
  st_transform(crs = proj) %>%
  st_make_valid %>%
  st_buffer(10000) %>%
  st_union

# convert range multipolygons to multilinestrings, crop lines that overlap with buffered lakeshore/coastline

line.ls <- sapply(spp, function(x) {
  shp.ls[[x]] %>%
  st_cast("MULTILINESTRING") %>%
  st_union %>%
  st_difference(., coast_shp) %>%
  st_difference(., lakes_shp)
}, simplify = FALSE)

# read in sample coords, convert CRS

all.df <- read.csv("~/qmac/coords_noSHF.csv", header = TRUE) %>%
  dplyr::mutate(spp = gsub(".*_", "", ind_ID)) %>%
  dplyr::select(ind_ID, spp, long, lat)
crds <- all.df %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = proj)

# read in F1 coords, convert CRS

F1.df <- read.table("~/qmac/F1_inds.txt", header = TRUE) %>%
  left_join(all.df, by = join_by(ind_ID)) %>%
  drop_na %>%
  dplyr::select(ind_ID, P1, P2, long, lat)
F1crds <- F1.df %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = proj)

# find nearest distance from each point to corresponding cropped range edge, add to dataframe with coordinates and ind_IDs

crds <- crds %>%
  dplyr::filter(!(ind_ID %in% F1.df$ind_ID)) %>% # remove F1s
  rowwise %>%
  dplyr::mutate(dist_from_rangeedge = st_distance(geometry, line.ls[[spp]])) %>%
  ungroup %>%
  as.data.frame %>%
  dplyr::select(ind_ID, dist_from_rangeedge) %>%
  dplyr::mutate(dist_from_rangeedge = as.numeric(dist_from_rangeedge))

F1crds <- F1crds %>%
  rowwise %>%
  dplyr::mutate(dist_from_rangeedge = min(
    st_distance(geometry, line.ls[[P1]]),
    st_distance(geometry, line.ls[[P2]])
  )) %>% # retain smallest distance to range edge calculated for both parent species
  ungroup %>%
  as.data.frame %>%
  dplyr::select(ind_ID, dist_from_rangeedge) %>%
  dplyr::mutate(dist_from_rangeedge = as.numeric(dist_from_rangeedge))

# join distances with coordinate dataframe

out.df <- rbind(crds, F1crds) %>%
  full_join(all.df, ., by = join_by(ind_ID)) %>%
  dplyr::select(-spp)

# write to CSV file

write.csv(out.df, "~/qmac/analyses/envdata/coords_noSHF_distfromrangeedge_NEW.csv", quote = FALSE, row.names = FALSE)
