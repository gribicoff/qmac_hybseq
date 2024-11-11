args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
spp <- c("mac", "alb", "bic", "ste", "mue")

# args <- "admixtableK6_allplinkinds.txt"

# subset unadmixed individuals (Q > .95) for each species for hybriddetective's SNP panel selection procedure

pwinds.df <- read.table(args[1], header = TRUE) %>%
  rowwise %>%
  dplyr::mutate(
    sppcode = as.factor(gsub(".*_", "", ind_ID)),
    mac_frac = mac1_frac + mac2_frac) %>%
    ungroup %>%
    dplyr::select(-c(mac1_frac, mac2_frac)) %>%
    dplyr::filter(
      case_when(
        sppcode == "mac" ~ mac_frac > 0.95,
        sppcode == "alb" ~ alb_frac > 0.95,
        sppcode == "bic" ~ bic_frac > 0.95,
        sppcode == "ste" ~ ste_frac > 0.95,
        sppcode == "mue" ~ mue_frac > 0.95)) %>%
    dplyr::select(ind_ID, sppcode)

# write all individuals to file

write.table(pwinds.df[,1], "inds_for_hybriddetective_all.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# subset, write species pairs

sppcombn <- combn(spp, 2, paste0, collapse = "_")
combn.ls <- sapply(sppcombn, function(x) {
  pwinds.df %>%
    dplyr::filter(sppcode %in% str_split_1(x, "_")) %>%
    group_by(sppcode) %>%
    dplyr::arrange(desc(ind_ID), .by_group = TRUE) %>%
    ungroup %>%
    dplyr::select(ind_ID) %>%
    write.table(paste0("inds_for_hybriddetective_",x,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
})

# write list of species pairs

write.table(sppcombn, "speciespairs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
