args <- commandArgs(trailingOnly = TRUE)

# args <- c("traceplot_input.txt", "/path/to/output.pdf")

library(tidyverse)

# read in textfile with traceplots for species pairs, convert to long format

trace.df <- read.table(args[1], header = TRUE) %>%
    pivot_longer(cols = -c(itnum, runnum), names_to = "class", values_to = "value") %>%
    dplyr::mutate(across(c(class, runnum), as.factor)) %>%
    dplyr::rename(sweep = itnum, run = runnum)

# plot traces

trace.plot <- ggplot(trace.df, aes(x = sweep, y = value, col = class)) +
    geom_line(alpha = 0.5, linewidth = 0.25) + 
    facet_wrap(~ run) +
    scale_color_discrete(guide = "none") +
    ggtitle(toupper(gsub("_", " - " , str_extract(args[1], ".*(?=_)")))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_minimal()

pdf(args[2])
print(trace.plot)
dev.off()


