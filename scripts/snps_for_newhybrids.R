args <- commandArgs(trailingOnly = TRUE)

library(hybriddetective)
library(genepopedit)

getTopLoc(args[1], LDpop = "Both", r2.threshold = 0.2,  ld.window = 250, panel.size = 200, where.PLINK = "~/miniconda3/pkgs/plink-1.90b6.21-h031d066_5/bin/", where.PGDspider = "~/miniconda3/envs/R_newhybrids/share/pgdspider-2.1.1.5-1/")


# gpvec <- Sys.glob(args[1])
# gp.ls <- lapply(gpvec,function(x) {
#   getTopLoc(x, LDpop = "Both", r2.threshold = 0.2,  ld.window = 250, panel.size = 200, where.PLINK = "~/miniconda3/pkgs/plink-1.90b6.21-h031d066_5/bin/", where.PGDspider = "~/miniconda3/envs/R_newhybrids/share/pgdspider-2.1.1.5-1/")
# })