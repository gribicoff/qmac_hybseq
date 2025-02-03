# estimates overlap at the FIA plot level for 
# Q. macrocarpa vs other spp included in study
# data source: https://datadryad.org/stash/dataset/doi:10.5061/dryad.j4sf2
# original citation: Cavender-Bares, J., Kothari, S., Meireles, J. E., Kaproth, M. A., Manos, P. S., & Hipp, A. L. (2018). The role of diversification in community assembly of the oaks (Quercus L.) across the continental U.S. American Journal of Botany, 105(3), 565–586. https://doi.org/10.1002/ajb2.1049

if(!exists('dat')) 
    dat <- read.csv('analyses/new_analyses/fia/fia_subset_2025-01-24.csv')
if(!exists('dat.agg'))
    dat.agg <- aggregate(dat[grep('quercus', names(dat), value = T)], 
                        by = list(plotNum = dat$PLT_CN), FUN = sum)
albCol <- 2
bicCol <- 3
mueCol <- 4
steCol <- 5
macCol <- 6

overlaps <- c(
    mac = apply(dat.agg, 1, function(x) x[macCol] > 0) |> sum(),
    alb = apply(dat.agg, 1, function(x) x[albCol] > 0) |> sum(),
    mac_alb = 
        apply(dat.agg, 1, function(x) {x[macCol] > 0 & x[albCol] > 0}) |> sum(),
    mac_bic = 
        apply(dat.agg, 1, function(x) {x[macCol] > 0 & x[bicCol] > 0}) |> sum(),
    mac_mue = 
        apply(dat.agg, 1, function(x) {x[macCol] > 0 & x[mueCol] > 0}) |> sum(),
    mac_ste = 
        apply(dat.agg, 1, function(x) {x[macCol] > 0 & x[steCol] > 0}) |> sum(),
    alb_bic = 
        apply(dat.agg, 1, function(x) {x[albCol] > 0 & x[bicCol] > 0}) |> sum(),
    alb_mue = 
        apply(dat.agg, 1, function(x) {x[albCol] > 0 & x[mueCol] > 0}) |> sum(),
    alb_ste = 
        apply(dat.agg, 1, function(x) {x[albCol] > 0 & x[steCol] > 0}) |> sum()
    )

proportions <- c(
    by_mac = overlaps['mac_alb'] / overlaps['mac'],
    by_alb = overlaps['mac_alb'] / overlaps['alb'],
    overlaps['mac_mue'] / overlaps['mac'],
    overlaps['mac_ste'] / overlaps['mac'],
    overlaps['mac_bic'] / overlaps['mac'],
    overlaps['alb_mue'] / overlaps['alb'],
    overlaps['alb_bic'] / overlaps['alb'],
    overlaps['alb_ste'] / overlaps['alb']
)

print(overlaps)
print(round(proportions, 3))
writeLines(paste(names(proportions), as.character(round(proportions, 3))), 
            'analyses/new_analyses/fia/proportionalOverlaps-v1.txt')
writeLines(paste(names(overlaps), as.character(overlaps)),
            'analyses/new_analyses/fia/totalPlots-v1.txt')