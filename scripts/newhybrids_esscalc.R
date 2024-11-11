args <- commandArgs(trailingOnly = TRUE)

# args <- c("", "")

library(tidyverse)
library(data.table)
library(coda)

spp <- c("mac", "alb", "bic", "ste", "mue")
sppcombn <- combn(spp, 2, paste, collapse = "_")

# read in MCMC traceplots for population-averaged hybrid class frequencies (???)

traces.ls <- sapply(sppcombn, function(x) {
    fread(paste0("~/qmac/analyses/new_analyses/newhybrids/run_nh_new/", x, "_pitraces.txt"), header = TRUE)
}, simplify = FALSE)

