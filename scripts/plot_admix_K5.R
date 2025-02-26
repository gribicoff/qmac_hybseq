args <- commandArgs(trailingOnly = TRUE)

# args <- c("~/qmac/analyses/new_analyses/K5_admixtable_allplinkinds.txt", "/path/to/output/dir/")

library(tidyverse)

spp <- c("mac", "alb", "bic", "ste", "mue")

# read in admixture data and order levels for species assignment according to file

admix.df <- read.table(args[1], header = TRUE) %>%
    mutate(sppcode = factor(gsub(".*_", "", ind_ID), level = spp)) %>%
    mutate(
      across(contains("_frac"), ~ case_when(
        . <= 0.00001 ~ 0,
        . >= 0.99999 ~ 1,
        .default = .)))

# create vector of individual IDs sorted by species-specific Q value (ascending?) (WIP) and sort dataframe

ind_order <- lapply(spp, function(x) {
    admix.df %>%
        dplyr::filter(sppcode == x) %>%
        dplyr::select(ind_ID, paste0(x,"_frac")) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(!!rlang::sym(paste0(x,"_frac")))) %>%
        dplyr::select(ind_ID) %>%
        as.vector %>%
        unname
}) %>%
    unlist

# write ordering of individuals for admixture plotting of K = 6 - needs to be consistent across runs

write.table(ind_order, paste0(args[2], "ind_order_admixplot_fromK5.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# convert from wide to long format

admix.long <- admix.df %>%
  pivot_longer(names_to = "cluster", values_to = "Q", cols = contains("frac"))

# plot admixture results

admix.plotK5 <- ggplot(data = admix.long, aes(x = factor(ind_ID, level = ind_order), y = Q, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~ sppcode, scales = "free", space = "free") +
  scale_fill_manual(
    values = c(
      "mac_frac" = "#0059b3",
      "alb_frac" = "#a060ff76",
      "bic_frac" = "#06f10abb",
      "ste_frac" = "#fbff00",
      "mue_frac" = "#ff6a00"), 
    guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Ancestry Coefficient (Q-Value)") + 
  xlab(NULL) +
  theme(
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.text.x = element_text(angle = 90, size = 1, margin = margin(0.5, 0, 0, 0)),
    axis.ticks.length.x = unit(0, units = "cm"),
    strip.background = element_blank(),
    strip.text = element_blank())

# write pdf

ggsave(paste0(args[2], "admixture_K5.pdf"), admix.plotK5, width = 14, height = 4)
