args <- commandArgs(trailingOnly = TRUE)

# args <- c()

library(tidyverse)

spp <- c("mac", "alb", "bic", "ste", "mue")
spp_pairs <- combn(spp, 2) %>%
    as.data.frame %>%
    summarise(across(everything(), ~ paste0(.x, collapse = "_"))) %>%
    as.vector %>%
    unlist %>%
    unname

# read in NewHybrids outputs with posterior probabilities for each class

nh.ls <- sapply(spp_pairs, function(x) {
    read.table(paste0(x,"_admix_nh.txt"), header = TRUE) %>%
        dplyr::select(!(contains("frac"))) %>%
        dplyr::rename_with(~ gsub("nh\\.", "", .x), contains("nh")) %>%
        dplyr::mutate(
            purebred = !!rlang::sym(paste0("pure_",gsub("_.*", "", x))) + !!rlang::sym(paste0("pure_",gsub(".*_", "", x))),
            Bx1 = !!rlang::sym(paste0(gsub("_.*", "", x),"_F1")) + !!rlang::sym(paste0(gsub(".*_", "", x),"_F1")),
            Bx2 = !!rlang::sym(paste0(gsub("_.*", "", x),"_",gsub("_.*", "", x),"F1")) + !!rlang::sym(paste0(gsub(".*_", "", x),"_",gsub(".*_", "", x),"F1")),
            Bx3 = !!rlang::sym(paste0(gsub("_.*", "", x),"_",gsub("_.*", "", x),gsub("_.*", "", x),"F1")) + !!rlang::sym(paste0(gsub(".*_", "", x),"_",gsub(".*_", "", x),gsub(".*_", "", x),"F1")),
            .keep = "unused"
        ) %>%
        dplyr::select(ind_ID, spp, purebred, F1, F2, Bx1, Bx2, Bx3)
}, simplify = FALSE)

# collapse all pairwise dataframes into single dataframe - for duplicate individuals, retain the observation with the lowest purebred PP (all late-gen Bxs)

nh.df <- do.call(rbind, nh.ls)
rownames(nh.df) <- NULL
nh.df <- nh.df %>%
    group_by(ind_ID) %>%
    arrange(purebred, .by_group = TRUE) %>%
    distinct(ind_ID, .keep_all = TRUE)

# summarize distribution of hybrid classes across dataset

sum.df <- nh.df %>%
    ungroup %>%
    mutate(ind_class = case_when(
        purebred >= 0.9 ~ "pure",
        F1 >= 0.9 ~ "F1",
        F2 >= 0.9 ~ "F2",
        Bx1 >= 0.9 ~ "Bx1",
        Bx2 >= 0.9 ~ "Bx2",
        Bx3 >= 0.9 ~ "Bx3",
        Bx1 + Bx2 + Bx3 >= 0.9 & Bx1 < 0.9 & Bx2 < 0.9 & Bx3 < 0.9 ~ "Bx_noclass",
        Bx1 + Bx2 + Bx3 + purebred >= 0.9 & Bx1 < 0.9 & Bx2 < 0.9 & Bx3 < 0.9 & purebred < 0.9 ~ "purebred_Bx1",
        .default = "unknown"),
        .keep = "unused") %>%
    group_by(ind_class) %>%
    dplyr::summarise(count = n())

# write tables of individuals and summaries

write.csv(nh.df, "~/qmac/newhybrids_allinds.csv", row.names = FALSE, quote = FALSE)
write.csv(sum.df, "~/qmac/newhybrids_summary.csv", row.names = FALSE, quote = FALSE)