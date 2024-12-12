set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)

# args <- c("poplist.txt", "/path/to/subsamples/", "/path/to/plot.pdf", "/path/to/pval.csv")

# command line arguments: poplist (three letter sp codes in order of each individual)

library(tidyverse)
library(ggtext)
library(gaston)
library(hierfstat)

# read in poplist/get paths to bed files

spp <- c("mac", "alb", "bic", "ste", "mue")
poplist <- scan(args[1], what = "factor")
pw <- combn(spp, 2, paste0, collapse = "_")
combn(spp, 2)
files <- Sys.glob(paste0(args[2],"*bed"))

# generate pairwise Fst estimates per bed fileset, convert to df, add column for gene class

fst.ls <- sapply(files, function(x) {
  read.bed.matrix(x) %>%
    as.matrix %>%
    pairwise.fst.dosage(poplist) %>%
    .[t(combn(spp, 2))] %>%
    data.frame(fst = ., pw_pop = pw, class = gsub(".*/", "", gsub("genes.*", "", x)))
}, simplify = FALSE)

# convert list of dfs to single df

fst.df <- do.call("rbind", fst.ls)

# plot Fst distributions

fst.plot <- ggplot(data = fst.df, aes(x = fst, fill = class)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ pw_pop, scales = "free" , labeller = labeller(pw_pop = ~ gsub("_", " - ", toupper(.x)))) +
  xlab(bquote(F[ST])) +
  ylab(NULL) +
  scale_fill_manual(
    values = c("grey", "#01a6d4", "#e3be05"),
    name = "Gene\nClass",
    labels = c("Other", "BBP", "DT")) +
  theme_minimal()

pdf(args[3])
print(fst.plot)
dev.off()

# calculate P-values for each comparison

sum.df <- fst.df %>%
  group_by(pw_pop) %>%
  summarise(
    bud_pval1 = sum(mean(fst[class == "bud"]) > fst[class == "other"])/length(fst[class == "other"]),
    bud_pval2 = sum(mean(fst[class == "bud"]) < fst[class == "other"])/length(fst[class == "other"]),
    drought_pval1 = sum(mean(fst[class == "drought"]) > fst[class == "other"])/length(fst[class == "other"]),
    drought_pval2 = sum(mean(fst[class == "drought"]) < fst[class == "other"])/length(fst[class == "other"])) %>%
  rowwise %>%
  mutate(
    bud_other = min(bud_pval1, bud_pval2),
    drought_other = min(drought_pval1, drought_pval2),
    .keep = "unused") %>%
  ungroup %>%
  pivot_longer(names_to = "comparison", values_to = "pval", cols = -pw_pop) %>%
  mutate(pval = 2*pval)

# write table of P-values to csv

write.csv(sum.df, args[4], row.names = FALSE, quote = FALSE)
