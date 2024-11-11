args <- commandArgs(trailingOnly = TRUE)

# args <- c("admix_K6_fstmatrix.csv", "/path/to/output/dir/")

library(tidyverse)
library(Matrix)

# read in fst matrix and convert to long dataframe

clusters <- c("MAC_N", "MAC_S", "ALB", "BIC", "STE", "MUE")
fst.df <- read.csv(args[1], header = TRUE) %>%
    dplyr::rename(var1 = 1) %>%
    column_to_rownames("var1") %>%
    as.matrix %>%
    forceSymmetric(uplo = "L") %>%
    as.matrix %>%
    as.data.frame %>%
    rownames_to_column("var1") %>%
    dplyr::select(var1, tolower(clusters)) %>%
    arrange(factor(var1, levels = tolower(clusters))) %>%
    column_to_rownames("var1") %>%
    as.matrix
fst.df[lower.tri(fst.df, diag = TRUE)] <- NA
fst.df <- fst.df %>%
    as.data.frame %>%
    rownames_to_column("var1") %>%
    pivot_longer(names_to = "var2", values_to = "fst", cols = -1) %>%
    dplyr::mutate(
        across(where(is.character), toupper),
        var1 = factor(var1, levels = clusters),
        var2 = factor(var2, levels = rev(clusters)),
        fst = case_when(
            fst == 0 ~ NA,
            .default = fst)) %>%
    dplyr::filter(!(var1 == "MUE" | var2 == "MAC_N"))
    # dplyr::filter(!is.na(fst))
    # arrange(var1) %>%
    # group_by(var1) %>%
    # arrange(var2, .by_group = TRUE) %>%
    # ungroup

# generate pairwise Fst heatmap

fst.plot <- ggplot(data = fst.df, aes(x = factor(var1, levels = clusters[-length(clusters)]), y = factor(var2, levels = rev(clusters[-1])), fill = fst)) +
    geom_tile(color = NA) +
    geom_text(aes(label = fst), col = "white") +
    xlab(NULL) +
    ylab(NULL) +
    scale_fill_viridis_c(name = bquote(F[ST]), na.value = "white") +
    scale_x_discrete(position = "bottom") +
    theme_minimal() +
    theme(axis.ticks.length.x = unit(0, units = "cm"))

# write to pdf

pdf(paste0(args[2], "admixK6_fstheatmap.pdf"))
fst.plot
dev.off()