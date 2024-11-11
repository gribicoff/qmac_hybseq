# args <- c("filelist.txt","admixcoords_formodeling.csv","maccoords_formodeling.csv")

library(tidyverse)

# read in environmental data, Q value files

envfiles <- paste0("~/qmac/analyses/envdata/", scan("~/qmac/analyses/envdata/filelist.txt", what = "character"))
env.ls <- lapply(envfiles, function(x) read.csv(x, header = TRUE))
admix.df <- read.csv("~/qmac/analyses/new_analyses/admixcoords.csv", header = TRUE)

# combine envdata and admixdata into single dataframes

all.df <- Reduce(function(...) left_join(..., by = "ind_ID"), c(list(admix.df), env.ls))

# filter variables

all.admix <- all.df %>% 
    dplyr::mutate(
        mac_frac = mac1_frac + mac2_frac,
        species_frac = case_when(
            sppcode == "mac" ~ mac_frac,
            sppcode == "alb" ~ alb_frac,
            sppcode == "bic" ~ bic_frac,
            sppcode == "ste" ~ ste_frac,
            .default = mue_frac
        ),
        admix_frac = 1 - species_frac
    ) %>%
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

all.mac <- all.df %>% 
    dplyr::filter(sppcode == "mac" & mac1_frac + mac2_frac > 0.95) %>%
    select(-c(
        alb_frac,
        bic_frac,
        ste_frac,
        mue_frac,
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

all.admix <- all.admix %>%
    dplyr::mutate(sppcode = case_when(
        ind_ID %in% F1s ~ "F1",
        .default = sppcode
    ))

# write to csv

write.csv(all.mac, "~/qmac/analyses/new_analyses/qmac_formodeling_alldata.csv", quote = FALSE, row.names = FALSE)
write.csv(all.admix, "~/qmac/analyses/new_analyses/admix_formodeling_alldata.csv", quote = FALSE, row.names = FALSE)