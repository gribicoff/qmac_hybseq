args <- commandArgs(trailingOnly = TRUE)

# args <- c("/path/to/plink.mds", "/path/to/admixtable.txt", "/path/to/plot.pdf")

library(tidyverse)
library(scatterpie)

spp <- c("mac", "alb", "bic", "ste", "mue")

# read in MDS data

mds.df <- read.table(args[1], header = TRUE) %>%
    mutate(
        ind_ID = FID,
        sppcode = as.factor(gsub(".*_", "", ind_ID)),
        axis1 = C1,
        axis2 = C2,
        axis3 = C3,
        .keep = "none")

# read in K = 6 admixture data

admix.df <- read.table(args[2], header = TRUE) %>%
    mutate(
        ind_ID = ind_ID,
        sppcode = as.factor(gsub(".*_", "", ind_ID)),
        mac1_frac = mac1_frac,
        mac2_frac = mac2_frac,
        admix_frac = 1 - (mac1_frac + mac2_frac),
        .keep = "none") %>%
    dplyr::filter(sppcode == "mac" & admix_frac < 0.05) %>%
    dplyr::select(ind_ID, mac1_frac, mac2_frac, admix_frac) %>%
    left_join(mds.df, by = "ind_ID")

# plot NMDS

pdf(args[3])
ggplot() + 
    geom_point(data = mds.df, aes(x = axis1, y = axis2, fill = factor(sppcode, levels = spp)), stroke = NA, pch = 21, size = 2.5) + 
    scale_fill_manual(
        values = c("#74b9ff", "#ee90f3", "#06f10abb", "#e7e420", "#f96503"),
        labels = spp,
        guide = NULL) +
    theme_light() +
    xlab("NMDS Axis 1") +
    ylab("NMDS Axis 2")
ggplot() + 
    geom_point(data = mds.df, aes(x = axis1, y = axis3, fill = factor(sppcode, levels = spp)), stroke = NA, pch = 21, size = 2.5) + 
    scale_fill_manual(
        values = c("#74b9ff", "#ee90f3", "#06f10abb", "#e7e420", "#f96503"),
        labels = spp,
        guide = NULL) +
    theme_light() +
    xlab("NMDS Axis 1") +
    ylab("NMDS Axis 3")
ggplot() + 
    geom_point(data = mds.df, aes(x = axis2, y = axis3, fill = factor(sppcode, levels = spp)), stroke = NA, pch = 21, size = 2.5) + 
    scale_fill_manual(
        values = c("#74b9ff", "#ee90f3", "#06f10abb", "#e7e420", "#f96503"),
        labels = spp,
        guide = NULL) +
    theme_light() +
    xlab("NMDS Axis 2") +
    ylab("NMDS Axis 3")
ggplot() +
    geom_scatterpie(data = admix.df, aes(x = axis1, y = axis2), cols = c("mac1_frac", "mac2_frac", "admix_frac"), linewidth = 0.1, alpha = 0.75, pie_scale = 0.5) +
    scale_fill_manual(values = c("#0059b3", "#4dffff", "red"), labels = c(expression(MAC[N]), expression(MAC[S]), "ADMIX")) +
    labs(fill = "ADMIXTURE\nCluster") +
    xlab("NMDS Axis 1") +
    ylab("NMDS Axis 2") +
    theme_light() + 
    coord_fixed(ratio = 1) +
    theme(legend.title.align = 0.5) 
ggplot() +
    geom_scatterpie(data = admix.df, aes(x = axis1, y = axis3), cols = c("mac1_frac", "mac2_frac", "admix_frac"), linewidth = 0.1, alpha = 0.75, pie_scale = 0.5) +
    scale_fill_manual(values = c("#0059b3", "#4dffff", "red"), labels = c("MAC_N", "MAC_S", "ADMIX")) +
    labs(fill = "ADMIXTURE\nCluster") +
    xlab("NMDS Axis 1") +
    ylab("NMDS Axis 3") +
    theme_light() +
    theme(legend.title.align = 0.5) 
ggplot() +
    geom_scatterpie(data = admix.df, aes(x = axis2, y = axis3), cols = c("mac1_frac", "mac2_frac", "admix_frac"), linewidth = 0.1, alpha = 0.75, pie_scale = 0.5) +
    scale_fill_manual(values = c("#0059b3", "#4dffff", "red"), labels = c("MAC_N", "MAC_S", "ADMIX")) +
    labs(fill = "ADMIXTURE\nCluster") +
    xlab("NMDS Axis 2") +
    ylab("NMDS Axis 3") +
    theme_light() +
    theme(legend.title.align = 0.5) 
# ggplot() + 
#     geom_point(data = mds.df, aes(x = axis1, y = axis2, fill = factor(sppcode, levels = spp)), stroke = NA, pch = 21, size = 2.5) + 
#     scale_fill_manual(
#         values = c("#74b9ff", "#6f16f48f", "#06f10abb", "#fbff00", "#f99b03c9"),
#         labels = spp,
#         guide = NULL) +
#     theme_light() +
#     xlab("NMDS Axis 1") +
#     ylab("NMDS Axis 3")
dev.off()


