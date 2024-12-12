args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/qmac/analyses/new_analyses/figures/")

library(tidyverse)
library(data.table)
library(coda)

spp <- c("mac", "alb", "bic", "ste", "mue")
sppcombn <- combn(spp, 2, paste, collapse = "_")

# read in MCMC traceplots for mixing proportion parameters pi_g

traces.ls <- sapply(sppcombn, function(x) {
    fread(paste0("~/qmac/analyses/new_analyses/newhybrids/run_nh/", x, "_pitraces.txt"), header = TRUE) %>%
    dplyr::rename(
        Pure_P0 = Pure_0,
        Pure_P1 = Pure_1,
        Bx1_P0 = P0_F1,
        Bx2_P0 = P0_P0F1,
        Bx3_P0 = P0_P0P0F1,
        Bx1_P1 = P1_F1,
        Bx2_P1 = P1_P1F1,
        Bx3_P1 = P1_P1P1F1)
}, simplify = FALSE)

traces.ls <- sapply(sppcombn, function(x) {
    sapply(paste0("run", 1:4), function(z) {
        traces.ls[[x]] %>%
            dplyr::filter(runnum == z) %>%
            dplyr::slice_tail(n = 120000)
    }, simplify = FALSE)
}, simplify = FALSE)

# align hybrid class parameters across runs

# read in dataframes of hybrid class posterior probabilities for pairwise NewHybrids analyses (PofZ files)

PofZ.ls <- sapply(sppcombn, function(x) {
    sapply(paste0("run", 1:4), function(z) {
        read.table(paste0("~/qmac/analyses/new_analyses/newhybrids/run_nh/ess/", x, "_nh_rawlabels_", z, ".txt"), header = TRUE) %>%
            dplyr::mutate(spp = gsub(".*_", "", ind_ID)) %>%
            relocate(spp, .after = 1) %>%
            dplyr::rename(
                Pure_P0 = Pure_0,
                Pure_P1 = Pure_1,
                Bx1_P0 = P0_F1,
                Bx2_P0 = P0_P0F1,
                Bx3_P0 = P0_P0P0F1,
                Bx1_P1 = P1_F1,
                Bx2_P1 = P1_P1F1,
                Bx3_P1 = P1_P1P1F1)
    }, simplify = FALSE)
}, simplify = FALSE)

# find species that correspond to P0/P1 labels for each run for each species pair

spplab.ls <- sapply(sppcombn, function(x) {
    sapply(paste0("run", 1:4), function(i) {
        sapply(c("Pure_P0", "Pure_P1"), function(z) {
            PofZ.ls[[x]][[i]] %>%
                dplyr::filter(!!rlang::sym(z) == max(!!rlang::sym(z))) %>%
                slice_head(n = 1) %>%
                pull(spp) %>%
                data.frame(label = gsub(".*_", "", z), species = ., sppcombn = x, run = i)
        }, simplify = FALSE) %>% do.call(rbind, .)
    }, simplify = FALSE)
}, simplify = FALSE)

# replace P0/P1 with species labels, reorder columns

traces.ls <- sapply(sppcombn, function(x) {
    sapply(paste0("run", 1:4), function(i) {
        colnames(traces.ls[[x]][[i]]) <- colnames(traces.ls[[x]][[i]]) %>%
            gsub("P0", unname(unlist(subset(spplab.ls[[x]][[i]], label == "P0", select = species))), .) %>%
            gsub("P1", unname(unlist(subset(spplab.ls[[x]][[i]], label == "P1", select = species))), .)
        traces.ls[[x]][[i]] %>%
            dplyr::relocate(contains(gsub("_.*", "", x)), .before = contains(gsub(".*_", "", x)))
    }, simplify = FALSE)
}, simplify = FALSE)

# convert dataframes to coda 'mcmc' objects

mcmc.ls <- sapply(sppcombn, function(x) {
    sapply(paste0("run", 1:4), function(i) {
        traces.ls[[x]][[i]] %>%
            dplyr::select(-c(runnum, itnum)) %>%
            mcmc(thin = 5)
    }, simplify = FALSE) %>%
        as.mcmc.list
}, simplify = FALSE)

# calculate ESS across MCMC chains

ess.ls <- sapply(sppcombn, function(x) {
    mcmc.ls[[x]] %>%
        effectiveSize
}, simplify = FALSE)

# summarize ESS by pairwise comparison by hybrid class parameter

ess.df <- do.call(rbind, ess.ls) %>%
    as.data.frame %>%
    rownames_to_column(var = "sppcombn") %>%
    pivot_longer(cols = -sppcombn, values_to = "ESS", names_to = "pi_g") %>%
    dplyr::mutate(ESS = round(ESS))

# write out csv of parameter ESSs

write.csv(ess.df, paste0(args[1], "newhybrids_pi_ess.csv"), row.names = FALSE, quote = FALSE)