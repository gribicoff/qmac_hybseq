args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/qmac/analyses/new_analyses/K6_admixtable_allplinkinds.txt")

library(tidyverse)

# create vector of pairwise species combinations

spp <- c("mac", "alb", "bic", "ste", "mue")
sppcombn <- combn(spp, 2, paste0, collapse = "_")
pw.df <- 
pw.df <- setNames(as.data.frame(combn(spp,2)),sppcombn)

# read in admixture Q matrix

admix.df <- read.table(args[1], header = TRUE) %>%
  dplyr::mutate(
    sppcode = as.factor(gsub(".*_", "", ind_ID)),
    mac_frac = mac1_frac + mac2_frac) %>%
  dplyr::select(-c(mac1_frac, mac2_frac))

# subset individuals for pairwise hybridization detection - include only individuals with summed cluster assignment from species pair >0.95

pwinds.ls <- sapply(sppcombn, function(x) {
  admix.df %>%
    dplyr::filter(sppcode %in% str_split_1(x, "_") & get(paste0(str_split_1(x, "_")[1],"_frac")) + get(paste0(str_split_1(x, "_")[2],"_frac")) > 0.95) %>%
    dplyr::select(ind_ID)
}, simplify = FALSE)

# write sets of individuals to txt files

lapply(sppcombn, function(x) write.table(pwinds.ls[[x]], paste0("inds_for_newhybrids_",x,".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE))