set.seed(515)

args <- commandArgs(trailingOnly = TRUE)

# args <- c("admix_formodeling_alldata.csv", "qmac_formodeling_alldata.csv")

library(caret)
library(tidyverse)
library(brms)
library(bayesplot)

## ADMIXTURE MODELING

# read in admixture dataset (Q values + predictors)

admix_raw <- read.csv(args[1], header = TRUE)

# rename, set factor variables

admix_proc <- admix_raw %>%
    dplyr::select(-c(ind_ID, species_frac, sympspp_bin, long)) %>%
    mutate(across(c(sppcode, site), as.factor))

# identify and remove highly (|r| > 0.75) correlated predictor variables

cor_var <- admix_proc %>% 
    dplyr::select(-c(sppcode, site)) %>% 
    cor(use = "complete.obs") %>% 
    findCorrelation(cutoff = 0.75, names = TRUE)
admix_proc <- admix_proc %>% 
    dplyr::select(-cor_var)

# remove individuals with missing data

admix_proc <- admix_proc %>% drop_na

## admixture without F1s - testing species fixed effect for significance

admix_noF1 <- admix_proc %>%
    dplyr::filter(!(sppcode == "F1"))

# center and scale numerical predictors

admix_noF1 <- admix_noF1 %>%
    mutate(across(where(is.numeric) & !contains("admix"), scale))

# specify model

frmla <- as.formula("admix_frac ~ lat + dist_from_rangeedge + hab_suit + sympspp_num + sppcode + (1|site)")

# fit model

mod <- brm(frmla, data = admix_noF1, family = Beta(), iter = 50000, thin = 2, chains = 15, cores = 15, control = list(adapt_delta = 0.9), silent = 0)

# save plots of parameter posterior distributions as pdf

color_scheme_set("green")
# parameter traceplots
trace_plot <- mcmc_trace(as.array(mod),
    pars = vars(contains("b_") & !contains("intercept"))) +
        xlab("Post-Burnin Iteration")
# parameter posterior distributions
ppd_plot <- mcmc_areas(mod,
    pars = vars(contains("b_") & !contains("intercept")),
    prob = 0.95,
    prob_outer = 0.99,
    point_est = "mean")
# parameter posterior intervals
ppi_plot <- mcmc_intervals(mod,
    pars = vars(contains("b_") & !contains("intercept")),
    prob = 0.95,
    prob_outer = 0.99,
    point_est = "mean")
# perform posterior predictive check
ppc <- posterior_predict(mod, cores = 10)
ppc_plot <- ppc_dens_overlay(admix_noF1$admix_frac, ppc[1:100,]) +
    xlim(0, 0.001)

pdf("brms_admixnoF1_summary.pdf")
trace_plot
ppd_plot
ppi_plot
ppc_plot
dev.off()

# save model fit summary as text file

sink("brms_admixnoF1_summary.txt")
print(summary(mod))
sink()

## admixture with F1s, no species fixed effect

# center and scale numerical predictors

admix_proc <- admix_proc %>%
    mutate(across(where(is.numeric) & !contains("admix"), scale))

# specify model

frmla <- as.formula("admix_frac ~ lat + dist_from_rangeedge + hab_suit + sympspp_num + (1|site)")

# fit model

mod <- brm(frmla, data = admix_proc, family = Beta(), iter = 50000, thin = 2, chains = 15, cores = 15, control = list(adapt_delta = 0.9), silent = 0)

# save plots of parameter posterior distributions as pdf

color_scheme_set("pink")
# parameter traceplots
trace_plot <- mcmc_trace(as.array(mod),
    pars = vars(contains("b_") & !contains("intercept"))) +
        xlab("Post-Burnin Iteration")
# parameter posterior distributions
ppd_plot <- mcmc_areas(mod,
    pars = vars(contains("b_") & !contains("intercept")),
    prob = 0.95,
    prob_outer = 0.99,
    point_est = "mean")
# parameter posterior intervals
ppi_plot <- mcmc_intervals(mod,
    pars = vars(contains("b_") & !contains("intercept")),
    prob = 0.95,
    prob_outer = 0.99,
    point_est = "mean")
# perform posterior predictive check
ppc <- posterior_predict(mod, cores = 10)
ppc_plot <- ppc_dens_overlay(admix_proc$admix_frac, ppc[1:100,]) +
    xlim(0, 0.001)

pdf("brms_admix_summary.pdf")
trace_plot
ppd_plot
ppi_plot
ppc_plot
dev.off()

# save model fit summary as text file

sink("brms_admix_summary.txt")
print(summary(mod))
sink()

## QMAC POPULATION STRUCTURE MODELING

# read in admixture/Q. macrocarpa population structure datasets (Q values + predictors)

mac_raw <- read.csv(args[2], header = TRUE)

# rename, reorder factor levels for flooding/ponding frequency, drainage class

levels(mac_raw$pondfreqcl_dc) <- list(
    none_pond="None",
    very_rare_pond="Very rare",
    rare_pond="Rare",
    occ_pond="Occasional",
    freq_pond="Frequent",
    very_freq_pond="Very Frequent")
levels(mac_raw$flodfreqcl_dc) <- list(
    none_flood="None",
    very_rare_flood="Very rare",
    rare_flood="Rare",
    occ_flood="Occasional",
    freq_flood="Frequent",
    very_freq_flood="Very Frequent")
levels(mac_raw$drclassdcd) <- list(
    exc_drain="Excessively drained",
    some_exc_drain="Somewhat excessively drained",
    well_drain="Well drained",
    mod_well_drain="Moderately well drained",
    some_poor_drain="Somewhat poorly drained",
    poor_drain="Poorly drained",
    very_poor_drain="Very poorly drained",
    subaq_drain="Subaqueous")

# rename, set factor variables, subset variables of ecological interest

mac_proc <- mac_raw %>%
    dplyr::rename(pond_freq = pondfreqcl_dc, flood_freq = flodfreqcl_dc, drain_class = drclassdcd, pH = ph1to1h2o_r) %>%
    mutate(across(c(pond_freq, flood_freq, drain_class, site), as.factor)) %>%
    dplyr::select(mac2_frac, bio4, bio5, bio12, moist_index, flood_freq, pond_freq, drain_class, site)

# identify and remove highly (|r| > 0.75) correlated predictor variables

cor_var <- mac_proc %>% 
    dplyr::select(-c(mac2_frac, where(is.factor))) %>%
    cor(use = "complete.obs") %>% 
    findCorrelation(cutoff = 0.75, names = TRUE)
mac_proc <- mac_proc %>% 
    dplyr::select(-cor_var)

# remove individuals with missing data

mac_proc <- mac_proc %>% drop_na

# center and scale numerical predictors

mac_proc <- mac_proc %>%
    mutate(across(where(is.numeric) & !contains("mac2"), scale))

# specify model

frmla <- as.formula("mac2_frac ~ bio4 + bio5 + bio12 + moist_index + flood_freq + pond_freq + drain_class + (1|site)")

# fit model

mod <- brm(frmla, data = mac_proc, family = Beta(), iter = 50000, thin = 2, chains = 15, cores = 15, control = list(adapt_delta = 0.9), silent = 0)


color_scheme_set("brightblue")
# parameter traceplots
trace_plot <- mcmc_trace(as.array(mod),
    pars = vars(contains("b_") & !contains("intercept"))) +
        xlab("Post-Burnin Iteration")
# parameter posterior distributions
ppd_plot <- mcmc_areas(mod,
    pars = vars(contains("b_") & !contains("intercept")),
    prob = 0.95,
    prob_outer = 0.99,
    point_est = "mean")
# parameter posterior intervals
ppi_plot <- mcmc_intervals(mod,
    pars = vars(contains("b_") & !contains("intercept")),
    prob = 0.95,
    prob_outer = 0.99,
    point_est = "mean")
# perform posterior predictive check
ppc <- posterior_predict(mod, cores = 10)
ppc_plot <- ppc_dens_overlay(mac_proc$mac2_frac, ppc[1:100,])

pdf("brms_qmacstruct_summary.pdf")
trace_plot
ppd_plot
ppi_plot
ppc_plot
dev.off()

# save model fit summary as text file

sink("brms_qmacstruct_summary.txt")
print(summary(mod))
sink()
