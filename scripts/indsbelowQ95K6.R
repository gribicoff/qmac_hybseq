library(tidyverse)

# read in admixture fractions for individuals in K = 6 run

spp <- c("mac", "alb", "bic", "ste", "mue")

admix.df <- read.table("~/qmac/analyses/new_analyses/K6_admixtable_allplinkinds.txt", header = TRUE) %>%
    mutate(
        sppcode = gsub(".*_", "", ind_ID),
        mac_frac = mac1_frac + mac2_frac) %>%
    dplyr::select(-c(mac1_frac, mac2_frac))

# calculate heterospecific Q-values and filter individuals below, above thresholds

indsbelow95 <- sapply(spp, function(x) {
    admix.df %>%
        dplyr::filter(sppcode == x) %>%
        dplyr::mutate(admix_frac = 1 - .data[[paste0(x, "_frac")]]) %>%
        dplyr::filter(admix_frac > 0.05) %>%
        dplyr::select(ind_ID) %>%
        as.vector %>%
        unname
}, simplify = FALSE) %>%
    unlist %>%
    unname

indsabove95 <- admix.df %>%
    dplyr::filter(!(ind_ID %in% indsbelow95)) %>%
    dplyr::select(ind_ID)

# write textfiles of lists of individuals

write.table(indsabove95, "~/qmac/analyses/new_analyses/indsaboveQ95K6.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(indsbelow95, "~/qmac/analyses/new_analyses/indsbelowQ95K6.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)