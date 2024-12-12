args <- commandArgs(trailingOnly = TRUE)

# args <- c("speciespairs.txt", "~/qmac/analyses/new_analyses/figures/")

library(tidyverse)

spp <- c("mac", "alb", "bic", "ste", "mue")
pairs <- scan(args[1], what = "character")

# read in hybrid class posterior probability tables for each pairwise species comparison

nh.ls <- sapply(pairs, function(x) {
    read.table(paste0(x, "_nh.txt"), header = TRUE) %>%
        dplyr::mutate(
            Pure = !!rlang::sym(paste0("Pure_", str_split_1(x, "_")[1])) + !!rlang::sym(paste0("Pure_", str_split_1(x, "_")[2])),
            Bx1 = !!rlang::sym(paste0("Bx1_", str_split_1(x, "_")[1])) + !!rlang::sym(paste0("Bx1_", str_split_1(x, "_")[2])),
            Bx2 = !!rlang::sym(paste0("Bx2_", str_split_1(x, "_")[1])) + !!rlang::sym(paste0("Bx2_", str_split_1(x, "_")[2])),
            Bx3 = !!rlang::sym(paste0("Bx3_", str_split_1(x, "_")[1])) + !!rlang::sym(paste0("Bx3_", str_split_1(x, "_")[2])),
            ind_class = factor(case_when(
                Pure >= 0.9 ~ "Pure",
                F1 >= 0.9 ~ "F1",
                F2 >= 0.9 ~ "F2",
                Bx1 >= 0.9 ~ "Bx1",
                Bx2 >= 0.9 ~ "Bx2",
                Bx3 >= 0.9 ~ "Bx3",
                Bx1 + Bx2 + Bx3 >= 0.9 & Bx1 < 0.9 & Bx2 < 0.9 & Bx3 < 0.9 ~ "Bx_undetermined",
                Bx1 + Bx2 + Bx3 + Pure >= 0.9 & Bx1 < 0.9 & Bx2 < 0.9 & Bx3 < 0.9 & Pure < 0.9 ~ "Pure_Bx_undetermined",
                .default = "Other"), levels = c("Pure", "F1", "F2", "Bx1", "Bx2", "Bx3", "Bx_undetermined", "Pure_Bx_undetermined", "Other")),
            parent_spp = ifelse(ind_class == "Pure", spp, x),
            spp = factor(spp, levels = str_split_1(x, "_")),
            ind_ID = gsub("_.*", "", ind_ID)) %>%
        dplyr::select(ind_ID, spp, parent_spp, Pure, F1, F2, Bx1, Bx2, Bx3, ind_class)
}, simplify = FALSE)

# get summary of classes for each species pair (not broken down by species)

sum.ls <- sapply(pairs, function(x) {
    nh.ls[[x]] %>%
    group_by(ind_class) %>%
    dplyr::summarise(count = n()) %>%
    complete(ind_class, fill = list(count = 0))
}, simplify = FALSE)

# get summary of classes for each species pair (broken down by species)

sumspp.ls <- sapply(pairs, function(x) {
    nh.ls[[x]] %>%
    group_by(ind_class, spp) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    tidyr::complete(ind_class, spp, fill = list(count = 0))
}, simplify = FALSE)

# collapse all pairwise dataframes into single dataframe - for duplicate individuals, retain the observation with the lowest purebred PP (all late-gen Bxs)

all.df <- do.call(rbind, nh.ls) %>%
    remove_rownames %>%
    group_by(ind_ID) %>%
    arrange(Pure, .by_group = TRUE) %>%
    distinct(ind_ID, .keep_all = TRUE) %>%
    ungroup

# get list of F1s for excluding from admixture GLMMs, write to table

F1s.df <- all.df %>%
    dplyr::filter(ind_class == "F1") %>%
    rowwise %>%
    dplyr::mutate(
        ind_ID = paste0(ind_ID, "_", spp),
        .keep = "none") %>%
    ungroup

write.table(F1s.df, "~/qmac/analyses/new_analyses/F1_inds.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# summarize distribution of hybrid classes across dataset

sum.df <- all.df %>%
    mutate(ind_class = factor(case_when(
        Pure >= 0.9 ~ "Pure",
        F1 >= 0.9 ~ "F1",
        F2 >= 0.9 ~ "F2",
        Bx1 >= 0.9 ~ "Bx1",
        Bx2 >= 0.9 ~ "Bx2",
        Bx3 >= 0.9 ~ "Bx3",
        Bx1 + Bx2 + Bx3 >= 0.9 & Bx1 < 0.9 & Bx2 < 0.9 & Bx3 < 0.9 ~ "Bx_undetermined",
        Bx1 + Bx2 + Bx3 + Pure >= 0.9 & Bx1 < 0.9 & Bx2 < 0.9 & Bx3 < 0.9 & Pure < 0.9 ~ "Pure_Bx_undetermined",
        .default = "Unclassified"), levels = c("Pure", "F1", "F2", "Bx1", "Bx2", "Bx3", "Bx_undetermined", "Pure_Bx_undetermined", "Unclassified")),
        .keep = "unused") %>%
    group_by(ind_class, .drop = FALSE) %>%
    dplyr::summarise(count = n())

# write tables of individuals and summaries, pairwise summaries

write.table(all.df, paste0(args[2], "newhybrids_summary_allinds.txt"), row.names = FALSE, quote = FALSE)
write.table(sum.df, paste0(args[2], "newhybrids_summary_allspp.txt"), row.names = FALSE, quote = FALSE)
lapply(pairs, function(x) write.table(sum.ls[[x]], paste0(args[2],"newhybrids_summary_",x,"_nospp.txt"), row.names = FALSE, quote = FALSE))
lapply(pairs, function(x) write.table(sumspp.ls[[x]], paste0(args[2],"newhybrids_summary_",x,"_byspp.txt"), row.names = FALSE, quote = FALSE))
