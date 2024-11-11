args <- commandArgs(trailingOnly = TRUE)

# args <- c(indlist,"onesnppergene_sample_X_bestadmixrun_perK.txt","destination/path/for/pdf")

library(tidyverse)
library(grid)
library(gridExtra)

spp <- c("mac", "alb", "bic", "ste", "mue")

# read in list of individuals, create dataframe of individuals and species

inds.df <- read.table(args[1], header = FALSE) %>%
  dplyr::rename(ind_ID = 1) %>%
  dplyr::mutate(sppcode = factor(gsub(".*_", "", ind_ID), level = spp))

# read in Q matrices across K, add characters

admix.ls <- scan(args[2], what = "character") %>%
  lapply(function(x) {
    read.table(x, header = FALSE) %>%
      cbind(inds.df, .) %>%
      pivot_longer(names_to = "cluster", values_to = "Q", cols = contains("V"))
  }) %>%  
  setNames(paste0("K", 2:8))

# read in ordering of individuals for admixture plot (consistent across all values of K and runs)

ind_order <- read.table("~/qmac/ind_order_admixplot_fromK5.txt", header = FALSE)[,1]

# plot admixture results

plot.ls <- sapply(admix.ls, function(df) {
  ggplot(data = df, aes(x = factor(ind_ID, level = ind_order), y = Q, fill = cluster)) +
    geom_col(width = 1) +
    facet_grid(~ sppcode, scales = "free", space = "free") +
    scale_fill_brewer(type = "qual", palette = "Set3", guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    ylab(NULL) + 
    xlab(NULL) +
    ggtitle(paste0("K = ", length(unique(df$cluster)))) +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      axis.text.x = element_text(angle = 90, size = 1, margin = margin(0.5, 0, 0, 0)),
      axis.ticks.length.x = unit(0, units = "cm"),
      axis.text.y = element_text(size = 4),
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.title = element_text(size = 9, margin = margin(0, 0, 0, 0)),
      plot.margin = unit(c(0.05, 0.1, 0.05, 0.25), "cm"))
}, simplify = FALSE)

# save to PDF

grid.arrange(grobs = plot.ls, ncol = 1, nrow = length(plot.ls), left = textGrob("Ancestry Coefficient (Q-Value)", rot = 90, vjust = 1, gp = gpar(col = "black", fontsize = 7))) %>%
  ggsave(paste0(args[3], str_extract(args[2], ".*(?<=[0-9])"), "_admixplot_acrossK.pdf"), plot = .)