args <- commandArgs(trailingOnly = TRUE)

# args <- c("eigenvalues", "eigenvectors", "outputname.pdf")

library(tidyverse)

# make dataframe of variance explained by each principal component for screeplot

scree.df <- data.frame(
    pc = as.factor(1:length(scan(args[1], what = "numeric"))),
    ve = as.numeric(scan(args[1], what = "numeric"))/sum(as.numeric(scan(args[1], what = "numeric")))
)

# make screeplot

scree.plot <- ggplot() +
    geom_col(data = dplyr::slice_head(scree.df, n = 10), aes(x = pc, y = ve), fill = "lightblue") +
    xlab("Principal Component") +
    ylab("Proportion of Variance") +
    ggtitle("Scree Plot of Thinned SNP PCA") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_minimal()

# read in dataframe of eigenvectors

pca.df <- read.table(args[2], header = FALSE) %>%
    dplyr::select(-2) %>%
    dplyr::rename(ind_ID = 1) %>%
    dplyr::rename_with(~ paste0("PC", as.numeric(gsub("V", "", .x)) - 2), where(is.numeric)) %>%
    dplyr::mutate(sppcode = as.factor(gsub(".*_", "", ind_ID))) %>%
    dplyr::relocate(sppcode, .after = ind_ID)

pca12.plot <- pca.df %>%
    ggplot(aes(x = PC1, y = PC2, color = sppcode)) +
    geom_point(alpha = 0.8) +
    xlab(paste0("PC1 (", round(100*scree.df[1,2], 3), "%)")) +
    ylab(paste0("PC2 (", round(100*scree.df[2,2], 3), "%)")) +
    theme_minimal()

pca13.plot <- pca.df %>%
    ggplot(aes(x = PC1, y = PC3, color = sppcode)) +
    geom_point(alpha = 0.8) +
    xlab(paste0("PC1 (", round(100*scree.df[1,2], 3), "%)")) +
    ylab(paste0("PC3 (", round(100*scree.df[3,2], 3), "%)")) +
    theme_minimal()

pca23.plot <- pca.df %>%
    ggplot(aes(x = PC2, y = PC3, color = sppcode)) +
    geom_point(alpha = 0.8) +
    xlab(paste0("PC2 (", round(100*scree.df[2,2], 3), "%)")) +
    ylab(paste0("PC3 (", round(100*scree.df[3,2], 3), "%)")) +
    theme_minimal()


pca14.plot <- pca.df %>%
    ggplot(aes(x = PC1, y = PC4, color = sppcode)) +
    geom_point(alpha = 0.8) +
    xlab(paste0("PC1 (", round(100*scree.df[1,2], 3), "%)")) +
    ylab(paste0("PC4 (", round(100*scree.df[4,2], 3), "%)")) +
    theme_minimal()

pdf(args[3])
scree.plot
pca12.plot
pca13.plot
pca23.plot
pca14.plot
dev.off()