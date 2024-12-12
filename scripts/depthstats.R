args <- commandArgs(trailingOnly = TRUE)

# args <- depthfile

library(tidyverse)

read.table(args[1], header = TRUE) %>%
    dplyr::mutate(mean_depth = MEAN_DEPTH, .keep = "none") %>%
    summarise(
        n_snps = n(),
        depth_mean = mean(mean_depth),
        depth_sd = sd(mean_depth)
    ) %>%
    print(n = "all")