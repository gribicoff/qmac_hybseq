args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

ind_miss <- read_delim(paste(args, ".imiss", sep = ""), delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
rm_ind_miss <- ind_miss[ind_miss$fmiss >= .35,]
pdf(paste(args, "-ind-missing-rm.pdf", sep = ""))
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "blue", colour = "black", alpha = 0.25)
a + theme_light()
dev.off()
ind_to_rm_list <- as.vector(rm_ind_miss$ind)
capture.output(cat(ind_to_rm_list, sep = "\n"), file = paste(args, "-inds-to-remove.txt", sep = ""))
perc_retained <- (1 - length(rm_ind_miss$ind) / length(ind_miss$ind)) * 100
capture.output(paste(perc_retained, "% of individuals retained.", sep = ""), file = paste(args, "-percent-inds-retained.txt", sep = ""))
