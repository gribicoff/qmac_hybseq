args <- commandArgs(trailingOnly = TRUE)

# args <- c("admixsitecoords.csv", "./figures/")

library(tidyverse)

sites.df <- read.csv(args[1], header = TRUE) %>%
    dplyr::mutate(species = case_when(
        sppcode == "mac" ~ "Q. macrocarpa",
        sppcode == "alb" ~ "Q. alba",
        sppcode == "bic" ~ "Q. bicolor",
        sppcode == "ste" ~ "Q. stellata",
        sppcode == "mue" ~ "Q. muehlenbergii"
    ), .keep = "unused") %>%
    dplyr::select(species, site, lat, long) %>%
    group_by(site, species) %>%
    # dplyr::summarise(lat = mean(lat), long = mean(long), num_inds = n()) %>%
    dplyr::summarise(lat = mean(lat), long = mean(long), n = n(), n = paste(n, unique(species), collapse = " ")) %>%
    group_by(site) %>%
    dplyr::summarise(lat = mean(lat), long = mean(long), num_inds = paste(n, collapse = "; ")) %>%
    ungroup %>%
    dplyr::mutate(across(everything(), ~ str_replace_all(.x, ",", ";")))

# write to textfile

write.csv(sites.df, paste0(args[2], "site_table.csv"), row.names = FALSE, quote = FALSE)