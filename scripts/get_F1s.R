library(tidyverse)

# read in NewHybrids hybrid class assignment file, select F1s and filter to individual ID

F1.df <- read.csv("newhybrids_allinds.csv", header = TRUE) %>%
    dplyr::filter(F1 > 0.9) %>%
    dplyr::select(ind_ID)

# read in admixture Q matrix, get two parental species based on two clusters with highest Q values

admix.df <- read.table("admixtableK6_allplinkinds.txt", header = TRUE) %>%
    dplyr::filter(ind_ID %in% F1.df$ind_ID) %>%
    dplyr::select(-sppcode) %>%
    mutate(mac_frac = mac1_frac + mac2_frac, .keep = "unused") %>%
    rowwise %>%
    mutate(
        P1 = names(.[,-1])[order(c_across(contains("frac")), decreasing = TRUE)[1]],
        P2 = names(.[,-1])[order(c_across(contains("frac")), decreasing = TRUE)[2]],
        across(P1:P2, ~ gsub("_.*", "", .x))) %>%
    dplyr::select(ind_ID, P1, P2)

# write out table of F1s and parent species

write.table(admix.df, "F1_inds.txt", row.names = FALSE, quote = FALSE)