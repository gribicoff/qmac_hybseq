args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

# args <- "cverror_filename_stem"

# read in CV error table

cv <- read.table(paste0(args[1], ".txt"), header = TRUE)

# plot CV error

pdf(paste0(args[1], "_plot.pdf"))
ggplot(data = cv, aes(x = K, y = cv_err)) +
    geom_line(col = "red", linewidth = 1) +
    geom_point(col = "#be0404") +
    xlab("K (Ancestral Clusters)") +
    ylab("CV Error") +
    scale_x_continuous(breaks = seq(min(cv$K), max(cv$K), 1)) +
    theme_light()
dev.off()
