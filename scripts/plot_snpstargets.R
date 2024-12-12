args <- commandArgs(trailingOnly = TRUE)

# args <- c("/path/to/plink.bim", "/path/to/outputdir/")

library(tidyverse)

# read in SNP, target locus positions 

snps <- read.table(args[1], header = FALSE) %>%
    dplyr::select(1,4) %>%
    dplyr::rename(chr = 1, pos = 2) %>%
    dplyr::mutate(chr = case_when(chr == 44904 ~ "NC_044904.1",
                                chr == 44905 ~ "NC_044905.1",
                                chr == 44906 ~ "NC_044906.1",
                                chr == 44907 ~ "NC_044907.1",
                                chr == 44908 ~ "NC_044908.1",
                                chr == 44909 ~ "NC_044909.1",
                                chr == 44910 ~ "NC_044910.1",
                                chr == 44911 ~ "NC_044911.1",
                                chr == 44912 ~ "NC_044912.1",
                                chr == 44913 ~ "NC_044913.1",
                                chr == 44914 ~ "NC_044914.1",
                                chr == 44915 ~ "NC_044915.1",
                                chr == 22154709 ~ "NW_022154709.1",
                                chr == 22155009 ~ "NW_022155009.1"))
genes <- read.table("~/qmac/analyses/blast/sorted-Qlobatagenes-PCG.bed", header = FALSE) %>%
    dplyr::select(1:3) %>%
    dplyr::rename(chr = 1, start = 2, stop = 3)

# read in chromosome/scaffold names, lengths

chrs <- read.table("~/qmac/analyses/align/index/QlobataChrom.sizes", header = FALSE) %>%
    dplyr::rename(chr = 1, length = 2) %>%
    dplyr::filter(chr %in% unique(genes$chr))

# generate plots of target genes and SNPs on chromosomes

plotchr <- ggplot() +
    geom_segment(data = chrs, aes(x = chr, xend = chr, y = 0, yend = length), lineend = "round", col = "lightgrey", alpha = 0.75, size = 3) +
    geom_segment(data = genes, aes(x = chr, xend = chr, y = start - 50000, yend = stop + 50000), col = "red", size = 5) +
    geom_point(data = snps, aes(x = chr, y = pos), col = "black", size = 0.025) +
    scale_x_discrete(guide = guide_axis(angle = 90),
                    labels = c(paste0("Chr",1:12), chrs$chr[13:16])) +
    scale_y_continuous(labels = function(x)x/1000000) +
    xlab(NULL) +
    ylab("Mb") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme_classic()

pdf(paste0(args[2], "snplocusmap.pdf"))
print(plotchr)
dev.off()
