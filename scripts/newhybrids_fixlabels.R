args <- commandArgs(trailingOnly = TRUE)

# args <- c("/path/to/name_nh_rawlabels.txt", "/path/to/name_nh.txt")

library(tidyverse)

# read in dataframe of hybrid class posterior probabilities for pairwise NewHybrids analyses

nh.df <- read.table(args[1], header = TRUE) %>%
    dplyr::mutate(spp = as.factor(gsub(".*_", "", ind_ID))) %>%
    relocate(spp, .after = 1)

# find species that correspond to P0/P1 labels

spplabs <- sapply(c("Pure_0", "Pure_1"), function(x) {
    nh.df %>%
        dplyr::filter(!!rlang::sym(x) == max(!!rlang::sym(x))) %>%
        slice_head(n = 1) %>%
        pull(spp) %>%
        as.character
})

# replace 0/1 with species labels

colnames(nh.df) <- colnames(nh.df) %>%
    gsub("Pure_0", "Pure_P0", .) %>%
    gsub("Pure_1", "Pure_P1", .) %>%
    gsub("P0_F1", "Bx1_P0", .) %>%
    gsub("P0_P0F1", "Bx2_P0", .) %>%
    gsub("P0_P0P0F1", "Bx3_P0", .) %>%
    gsub("P1_F1", "Bx1_P1", .) %>%
    gsub("P1_P1F1", "Bx2_P1", .) %>%
    gsub("P1_P1P1F1", "Bx3_P1", .) %>%
    gsub("P0", spplabs[1], .) %>%
    gsub("P1", spplabs[2], .)

# write out new table

write.table(nh.df, args[2], quote = FALSE, row.names = FALSE)