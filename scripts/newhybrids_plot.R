args <- commandArgs(trailingOnly = TRUE)

# args = c("sp1_sp2", "/path/to/nh.txt", "/path/to/plot.pdf")

library(tidyverse)

# read in newhybrids data and order levels for species assignment according to file

pwspp <- str_split_1(args[1], "_")
nh.df <- read.table(args[2], header = TRUE) %>%
    dplyr::mutate(spp = factor(spp, levels = pwspp))

# create vector of individual IDs sorted by species-specific Q value (ascending?) (WIP) and sort dataframe, convert to long format

indorder <- nh.df %>%
    dplyr::mutate(
        purepp = case_when(
            spp == pwspp[1] ~ !!rlang::sym(paste0("Pure_", pwspp[1])),
            spp == pwspp[2] ~ !!rlang::sym(paste0("Pure_", pwspp[2])))) %>%
        group_by(spp) %>%
        arrange(desc(purepp), .by_group = TRUE) %>%
        pull(ind_ID)
nh.df <- nh.df %>%
    dplyr::mutate(ind_ID = factor(ind_ID, levels = indorder)) %>%
    pivot_longer(cols = -c(ind_ID, spp), names_to = "class", values_to = "pp")

# define ordering of hybrid class categories

classord <- c(
    paste0("Pure_", pwspp[1]),
    paste0("Pure_", pwspp[2]),
    "F1",
    "F2",
    paste0("Bx1_", pwspp[1]),
    paste0("Bx1_", pwspp[2]),
    paste0("Bx2_", pwspp[1]),
    paste0("Bx2_", pwspp[2]),
    paste0("Bx3_", pwspp[1]),
    paste0("Bx3_", pwspp[2]))

# set colors for each hybrid class

classcols <- c(
    "#006d2c",
    "#993404",
    "#818589",
    "#0096FF",
    "#b2e2e2",
    "#fed98e",
    "#66c2a4",
    "#fe9929",
    "#2ca25f",
    "#d95f0e")
names(classcols) <- classord

# create labels for hybrid classes

classlabs <- classord %>%
    gsub("_", " ", .) %>%
    gsub(pwspp[1], toupper(pwspp[1]), .) %>%
    gsub(pwspp[2], toupper(pwspp[2]), .)
names(classlabs) <- classord

# plot newhybrids results

nh.plot <- ggplot(data = nh.df, aes(x = ind_ID, y = pp, fill = class)) +
    geom_col(width = 1) +
    facet_grid(~ spp, scales = "free", space = "free", labeller = as_labeller(setNames(toupper(pwspp), pwspp))) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual("Hybrid Classes", values = classcols, labels = classlabs, breaks = classord) +
    guides(fill = guide_legend(ncol = 5)) +
    labs(y = "Posterior Probability", x = NULL) +
    theme(
        plot.margin = unit(c(0.1, 0, 0, 0), units = "cm"),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        axis.text.x = element_text(angle = 90, size = 2),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(face = "bold"),
        strip.background = element_blank())

# print NewHybrids barplot

pdf(args[3])
print(nh.plot)
dev.off()
