args <- commandArgs(trailingOnly = TRUE)

# args <- c("K6_admixtable_allplinkinds.txt", "admixsitecoords.csv", "quercus_hybseq_individuals_metadata.csv", "quercus_collectiondates.csv", "./figures/")

library(tidyverse)

# read in tables with georeferenced and non-georeferenced genotyped individuals

samps.df <- read.table(args[1], header = TRUE, na.strings = c("NA", "")) %>%
    dplyr::mutate(
        ind_ID = gsub("_.*", "", ind_ID),
        species = case_when(
            sppcode == "mac" ~ "Q. macrocarpa",
            sppcode == "alb" ~ "Q. alba",
            sppcode == "bic" ~ "Q. bicolor",
            sppcode == "ste" ~ "Q. stellata",
            sppcode == "mue" ~ "Q. muehlenbergii"),
        .keep = "none")
sites.df <- read.csv(args[2], header = TRUE, na.strings = c("NA", "")) %>%
    dplyr::mutate(ind_ID = gsub("_.*", "", ind_ID)) %>%
    dplyr::select(ind_ID, lat, long, site)

# join dataframes

samps.df <- left_join(samps.df, sites.df, by = "ind_ID")

# import and filter individual data

inds.df <- read.csv(args[3], header = TRUE, na.strings = c("NA", "")) %>%
    dplyr::mutate(
        ind_ID = Specimen.CODE,
        extraction_code = HYB.SEQ.EXTRACT,
        georef_source = Georeference_source,
        country = country.georef,
        state_province = state_province_georef,
        collector = collectors,
        collection_number = CollNum,
        .keep = "none")

# import dates, join dataframes

dates.df <- read.csv(args[4], header = TRUE, na.strings = c("NA", "")) %>%
    dplyr::rename(collection_date = collectionDate)

# join dataframes

samps.df <- left_join(samps.df, inds.df, by = "ind_ID") %>%
    dplyr::distinct() %>%
    dplyr::left_join(., dates.df, by = "ind_ID") %>%
    dplyr::mutate(
        georeferenced = ifelse(is.na(lat), "N", "Y"),
        collection_date = case_when(
            grepl("IL-PPWS-NS", collection_number) ~ "7/23/2015",
            collection_number == "IA-MG-269" ~ "6/14/2017",
            grepl("Salina", collection_number) ~ "2017",
            .default = collection_date)) %>%
    dplyr::relocate(georeferenced, .after = "species") %>%
    dplyr::mutate(across(everything(), ~ str_replace_all(.x, ",", ";")))

# write to table

write.csv(samps.df, paste0(args[5], "ind_table.csv"), row.names = FALSE, quote = FALSE)