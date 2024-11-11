args <- commandArgs(trailingOnly = TRUE)

# args <- c("admixcoords_noSHF.csv","~/qmac/analyses/model_admix/mactot95_inds.txt")

library(tidyverse)

# read in modified Q matrix with individuals, species codes, coordinates

admix.df <- read.csv(args[1], header = TRUE)

# read in list of unadmixed (< .05 nonmac Q) Q. macrocarpa individuals

macinds <- scan(args[2], what = "character")

# replace species code

admix.df <- admix.df %>% mutate(sppcode = gsub(".*_","",ind_ID))

# create new df with only unadmixed Q. macrocarpa individuals

mac.df <- admix.df %>% dplyr::filter(ind_ID %in% macinds)

# replace Q matrix columns with column of species-cluster Q value, admixture proportion (summed Q values for all other species clusters) 

admix.df <- admix.df %>% mutate(species_frac = case_when(sppcode == "mac" ~ mac1_frac + mac2_frac,
                                                    sppcode == "alb" ~ alb_frac,
                                                    sppcode == "bic" ~ bic_frac,
                                                    sppcode == "ste" ~ ste_frac,
                                                    sppcode == "mue" ~ mue_frac),
                                admix_frac = 1 - species_frac) %>% select(ind_ID, sppcode, species_frac, admix_frac, lat, long)

# rescale mac1, mac2 Q values

mac.df <- mac.df %>% mutate(mac1_frac = mac1_frac/(mac1_frac + mac2_frac),
                            mac2_frac = 1 - mac1_frac) %>% select(ind_ID, mac1_frac, mac2_frac, lat, long)

# write dataframes to csv files

write.csv(admix.df,"admixcoords_formodeling.csv",row.names=F,col.names=T,quote=F)
write.csv(mac.df,"maccoords_formodeling.csv",row.names=F,col.names=T,quote=F)