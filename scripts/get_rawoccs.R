library(rgbif)
library(CoordinateCleaner)
library(tidyverse)

# GBIF taxonKeys: 2878213 (mac), 2879737 (alb), 2878470 (bic), 8879961 (ste), 2879947 (mue)

sppkeys <- c("2878213", "2879737", "2878470", "8879961", "2879947")
spp <- c("mac", "alb", "bic", "ste", "mue")

# download occurrence data from GBIF - make sure GBIF API key is set

gbif_dl <- occ_download(
  pred_in("taxonKey",sppkeys),
  pred_in("country",c("US","CA")),
  pred("HAS_GEOSPATIAL_ISSUE",FALSE),
  pred("HAS_COORDINATE",TRUE),
  pred("OCCURRENCE_STATUS","PRESENT"),
  pred_not(pred_in("BASIS_OF_RECORD",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN"))),
  format="SIMPLE_CSV"
)
occ_download_wait(gbif_dl)

# import occurrence data from GBIF query into R and clean records with CoordinateCleaner

gbif.df <- gbif_dl %>%
  occ_download_get() %>%
  occ_download_import() %>%
  setNames(tolower(names(.))) %>%
  dplyr::filter(year >= 1900) %>%
  dplyr::filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>%
  dplyr::filter(coordinateuncertaintyinmeters < 1000 | is.na(coordinateuncertaintyinmeters)) %>%
  dplyr::filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>%
  cc_zero(buffer = 1, value = "clean", verbose = TRUE) %>%
  cc_cen(buffer = 2000, value = "clean", test = "both", verbose = TRUE) %>%
  cc_cap(buffer = 2000, value = "clean", verbose = TRUE) %>%
  cc_inst(buffer = 2000, value = "clean", verbose = TRUE) %>%
  cc_sea(value = "clean") %>%
  cc_urb(value = "clean") %>%
  dplyr::select(decimallongitude, decimallatitude, species) %>%
  dplyr::rename(long = decimallongitude, lat = decimallatitude) %>%
  dplyr::mutate(species = substr(gsub("^.* ", "", species), 1, 3))

# create set of dataframes filtered for each species, write each to separate CSV

gbif.ls <- sapply(spp, function(x) {
  gbif.df %>% 
    dplyr::filter(species == x) %>%
    dplyr::select(-species) %>%
    distinct(long, lat)
}, simplify = FALSE)

# import all FIA data

fiacodes <- list(mac = "823",
    alb = "802",
    bic = "804",
    ste = "835",
    mue = "826")
fia_raw <- readRDS("~/qmac/analyses/envdata/fia/trees.rds")

# filter FIA plots by species codes and add as dataframes to list

fia.ls <- sapply(spp, function(x) {
  fia_raw %>%
    dplyr::select(LAT, LON, SPCD) %>%
    dplyr::filter(SPCD == fiacodes[[x]]) %>%
    dplyr::mutate(long = LON, lat = LAT, .keep = "none") %>%
    remove_rownames %>%
    distinct(long, lat)
}, simplify = FALSE)

# merge GBIF and FIA datasets into single dataframe per species

allocc.ls <- sapply(spp, function(x) {
  rbind(gbif.ls[[x]], fia.ls[[x]]) %>%
  distinct(long, lat)
}, simplify = FALSE)

# write datasets to csv files for manual curation

lapply(spp, function(x) {
  write.csv(allocc.ls[[x]], file = paste0("rawoccs_", x, ".csv"), row.names = FALSE, quote = FALSE)
})