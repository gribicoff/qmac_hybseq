args <- commandArgs(trailingOnly=TRUE)

library("tidyverse", "ggplot2", "dplyr")

var_freq <- read_delim(paste(args, ".frq", sep = ""), delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
pdf(paste(args, "-freq.pdf", sep = ""))
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "yellow", colour = "black", alpha = 0.25)
a + theme_light()
dev.off()
