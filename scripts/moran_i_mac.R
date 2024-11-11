args <- commandArgs(trailingOnly = TRUE)

set.seed(1234)

library(spdep)
library(adegenet)
library(tidyverse)

# read in Q. macrocarpa coordinates/Q values, group by site

mac.df <- read.csv(args[1], header = TRUE) %>%
    dplyr::select(mac1_frac, long, lat, site) %>%
    group_by(site) %>%
    summarise(mac1_frac = mean(mac1_frac),
        long = mean(long),
        lat = mean(lat))

# calculate inverse distance weights

mac.listw <- mac.df %>% 
    dplyr::select(long, lat) %>%
    chooseCN(ask = FALSE, type = 7, a = 1, dmin = 0, plot.nb = FALSE)

# calculate Moran's I

mac.moran.param <- moran.test(mac.df$mac1_frac, mac.listw, alternative = "greater")

# calculate Moran's I with permutation test

mac.moran.perm <- moran.mc(mac.df$mac1_frac, mac.listw, alternative = "greater", nsim = 99999)

# plot permutation distribution alongside test statistic calculated

pdf("moranI_qmac_permutationdist.pdf")
ggplot() +
    geom_density(aes(x = mac.moran.perm$res), fill = "yellow", alpha = 0.75) +
    geom_vline(xintercept = mac.moran.perm$statistic, col = "red", alpha = 1) +
    xlab("Moran's I") +
    ylab(NULL) +
    ggtitle("Permutation Test of Moran's I for\nQuercus macrocarpa population structure") +
    theme_light()
dev.off()

# print test summaries (parametric and permutation) as textfile

sink("moranI_qmac.txt")
print(mac.moran.param)
print(mac.moran.perm)
sink()