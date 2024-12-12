library(tidyverse)

# read in expanded Q-matrix

admix.df <- read.table("K6_admixtable_allplinkinds.txt", header = TRUE) %>%
    dplyr::mutate(
        current_determination = gsub(".*_", "", ind_ID),
        specimen_code = gsub("_.*", "", ind_ID),
        .keep = "unused") %>%
    dplyr::select(-sppcode)

# read in newhybrids table

nh.df <- read.table("figures/newhybrids_summary_allinds.txt", header = TRUE) %>%
    dplyr::select(ind_ID, ind_class) %>%
    dplyr::rename(specimen_code = ind_ID, hybrid_class = ind_class)

# read in collection numbers

colnum.df <- read.csv("quercus_hybseq_individuals_metadata.csv", header = TRUE) %>%
    dplyr::select(Specimen.CODE, CollNum) %>% 
    dplyr::rename(specimen_code = Specimen.CODE, collection_number = CollNum)

# join dataframes

all.df <- left_join(admix.df, nh.df, by = "specimen_code") %>%
    left_join(., colnum.df, by = "specimen_code") %>%
    dplyr::select(
        collection_number,
        specimen_code,
        current_determination,
        hybrid_class,
        contains("frac")) %>%
    dplyr::mutate(
        current_determination = case_when(
            current_determination == "mac" ~ "Q. macrocarpa",
            current_determination == "alb" ~ "Q. alba",
            current_determination == "bic" ~ "Q. bicolor",
            current_determination == "ste" ~ "Q. stellata",
            current_determination == "mue" ~ "Q. muehlenbergii",
            .default = NA),
        hybrid_class = case_when(
            hybrid_class == "Pure" ~ "Unadmixed",
            hybrid_class == "Bx_undetermined" ~ "Bx; generation undetermined",
            hybrid_class == "Pure_Bx_undetermined" ~ "Unadmixed or Bx; undetermined",
            .default = hybrid_class))

# write to csv

write.csv(all.df, "hybseq_oak_genotypetable_GRNov2024.csv", row.names = FALSE, quote = FALSE)