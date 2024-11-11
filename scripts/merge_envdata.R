# args <- c("filelist.txt","admixcoords_formodeling.csv","maccoords_formodeling.csv")

library(tidyverse)

# read in environmental data, Q value files

envfiles <- scan("~/qmac/analyses/envdata/filelist.txt", what = "character")
env.ls <- lapply(envfiles, function(x) read.csv(x, header = TRUE))
admix.df <- read.csv("~/qmac/analyses/envdata/admixcoords_formodeling.csv", header = TRUE)
mac.df <- read.csv("~/qmac/analyses/envdata/maccoords_formodeling.csv", header = TRUE)
admix.ls <- c(list(admix.df),env.ls)
mac.ls <- c(list(mac.df),env.ls)

# combine envdata and admixdata into single dataframes for full dataset, Q. macrocarpa popstruct

admix.all <- Reduce(function(...) left_join(..., by = "ind_ID"), admix.ls)
mac.all <- Reduce(function(...) left_join(..., by = "ind_ID"), mac.ls)


admix.df <- admix.df %>%
    dplyr::mutate(
        con_frac = case_when(
            sppcode == "mac" ~ mac_frac,
            sppcode == "alb" ~ alb_frac,
            sppcode == "bic" ~ bic_frac,
            sppcode == "ste" ~ ste_frac,
            sppcode == "mue" ~ mue_frac),
        admix_frac = 1 - con_frac,
        sppcode = sppcode,
        .keep = "unused")

admix.df %>%
    group_by(sppcode) %>%
    summarise(
        n_tot = n(),
        n_ad = sum(admix_frac > 0.05),
        prop = n_ad/n_tot 
    )

# remove variables with excessive missing data (> 25%)

admix.all <- admix.all %>% 
    select(c(
        ind_ID,
        sppcode,
        species_frac,
        admix_frac,
        lat,
        long,
        site,
        dist_from_rangeedge,
        sympspp_bin,
        sympspp_num,
        hab_suit))

mac.all <- mac.all %>% 
    select(-c(
        min_surf_texture,
        mukey,
        musym,
        muname,
        areasymbol,
        ecec_r,
        dist_from_rangeedge,
        hab_suit,
        surf_texture,
        avg_soildepth,
        sympspp_num,
        sympspp_bin))

# replace sppcode for F1s

F1s <- read.table("~/qmac/F1_inds.txt", header = TRUE) %>%
    dplyr::select(ind_ID) %>%
    unlist %>%
    unname

admix.all <- admix.all %>%
    dplyr::mutate(sppcode = case_when(
        ind_ID %in% F1s ~ "F1",
        .default = sppcode
    ))

# write to csv

write.csv(mac.all, "~/qmac/analyses/envdata/qmac_formodeling_alldata.csv", quote = FALSE, row.names = FALSE)
write.csv(admix.all, "~/qmac/analyses/envdata/admix_formodeling_alldata.csv", quote = FALSE, row.names = FALSE)