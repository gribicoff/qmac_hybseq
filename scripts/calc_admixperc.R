library(tidyverse)

# args <- ("/path/to/admixtableK6_allplinkinds.txt")

admix.df <- read.table("admixtableK6_allplinkinds.txt", header = TRUE)
admix.df <- admix.df %>%
    mutate(sppcode = as.factor(sppcode)) %>%
    mutate(sppcode = case_when(sppcode == 1 ~ "mac",
                               sppcode == 2 ~ "alb",
                               sppcode == 3 ~ "bic",
                               sppcode == 4 ~ "ste",
                               sppcode == 5 ~ "mue"))

# number of admixed individuals: 0.9 Q-value threshold

admix.df %>% dplyr::filter(sppcode == "mac" & mac1_frac + mac2_frac < .9 |
                  sppcode == "alb" & alb_frac < .9 |
                  sppcode == "bic" & bic_frac < .9 |
                  sppcode == "ste" & ste_frac < .9 |
                  sppcode == "mue" & mue_frac < .9) %>%
    nrow %>%
    print(paste0("For the Q = 0.9 threshold, there are ",.,"/",nrow(admix.df)," admixed individuals"))

# number of admixed individuals: 0.95 Q-value threshold

admix.df %>% dplyr::filter(sppcode == "mac" & mac1_frac + mac2_frac < .95 |
                  sppcode == "alb" & alb_frac < .95 |
                  sppcode == "bic" & bic_frac < .95 |
                  sppcode == "ste" & ste_frac < .95 |
                  sppcode == "mue" & mue_frac < .95) %>%
    nrow

# number of admixed individuals: 0.99 Q-value threshold

admix.df %>% dplyr::filter(sppcode == "mac" & mac1_frac + mac2_frac < .99 |
                  sppcode == "alb" & alb_frac < .99 |
                  sppcode == "bic" & bic_frac < .99 |
                  sppcode == "ste" & ste_frac < .99 |
                  sppcode == "mue" & mue_frac < .99) %>%
    nrow
