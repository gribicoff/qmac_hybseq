args <- commandArgs(trailingOnly = TRUE)

# args <- c("K5.Q", "K6.Q", "indslist.txt")

library(tidyverse)

# read in list of individuals

inds <- scan(args[3], what = "character")

## K = 5

# read in list of individuals, K = 5 Q matrix and assign species codes

K5.df <- read.table(args[1], header = FALSE) %>%
    dplyr::mutate(
        ind_ID = inds,
        sppcode = as.factor(gsub(".*_", "", ind_ID))
    ) %>%
    relocate(-where(is.numeric))

# assign ancestral cluster names

colsK5 <- sapply(c("QUE002637_mac", "QUE002124_alb", "QUE001366_bic", "QUE001839_ste", "QUE001930_mue"), function(x) {
    K5.df %>%
        dplyr::filter(ind_ID == x) %>%
        mutate(col = paste0("V", max.col(.[,-(1:2)])[1]),
        .keep = "none") %>%
        unlist %>% unname
}, USE.NAMES = FALSE)

K5.df <- K5.df %>%
    dplyr::rename(
        mac_frac = colsK5[1],
        alb_frac = colsK5[2],
        bic_frac = colsK5[3],
        ste_frac = colsK5[4],
        mue_frac = colsK5[5]) %>%
    dplyr::select(
        ind_ID,
        sppcode,
        mac_frac,
        alb_frac,
        bic_frac,
        ste_frac,
        mue_frac)

# write out admixtable

write.table(K5.df, "K5_admixtable_allplinkinds.txt", quote = FALSE, row.names = FALSE)

## K = 6

# read in list of individuals, K = 6 Q matrix and assign species codes

K6.df <- read.table(args[2], header = FALSE) %>%
    dplyr::mutate(
        ind_ID = inds,
        sppcode = as.factor(gsub(".*_", "", ind_ID))
    ) %>%
    relocate(-where(is.numeric))

# assign ancestral cluster names

colsK6 <- sapply(c("QUE002637_mac", "QUE001937_mac", "QUE002124_alb", "QUE001366_bic", "QUE001839_ste", "QUE001930_mue"), function(x) {
    K6.df %>%
        dplyr::filter(ind_ID == x) %>%
        mutate(col = paste0("V", max.col(.[,-(1:2)])[1]),
        .keep = "none") %>%
        unlist %>% unname
}, USE.NAMES = FALSE)

K6.df <- K6.df %>%
    dplyr::rename(
        mac1_frac = colsK6[1],
        mac2_frac = colsK6[2],
        alb_frac = colsK6[3],
        bic_frac = colsK6[4],
        ste_frac = colsK6[5],
        mue_frac = colsK6[6]) %>%
    dplyr::select(
        ind_ID,
        sppcode,
        mac1_frac,
        mac2_frac,
        alb_frac,
        bic_frac,
        ste_frac,
        mue_frac)

# write out admixtable

write.table(K6.df, "K6_admixtable_allplinkinds.txt", quote = FALSE, row.names = FALSE)

# create and write dataframe of column labels

clusts.df <- data.frame(
    cluster = c("mac_n", "mac_s", "alb", "bic", "ste", "mue"),
    qmatrix_column = colsK6) %>%
    dplyr::mutate(admix_pop = paste0("Pop", as.numeric(gsub("V", "", qmatrix_column)) - 1)) %>%
    dplyr::select(cluster, admix_pop, qmatrix_column)
write.table(clusts.df, "K6_clusterlabels.txt", quote = FALSE, row.names = FALSE, sep = "\t")