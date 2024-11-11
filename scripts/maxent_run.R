options(java.parameters = "-Xmx32g")
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)

# args <- "sppname"

library(rJava)
library(terra)
library(ENMeval)

# read in presence points, background points, cropped environmental rasters

envrast <- rast(paste0("stackraster_",args[1],".tif"))
bgpts <- read.csv(paste0("maxentbgpts_",args[1],".csv"))
prspts20 <- read.csv(paste0("maxentprspts20_",args[1],".csv"),header = TRUE)

# specify Maxent feature classes to be tested

ftclasses <- c("L","Q","H","LQ","QH","LH","LQH","LQHP","LQHPT")

# run ENMeval

max20 <- ENMevaluate(occs = prspts20,
                        envs = envrast,
                        bg = bgpts,
                        partitions = "block",
                        algorithm = "maxent.jar",
                        tune.args = list(
                            fc = ftclasses,
                            rm = seq(0.5,5,0.5)),
                        parallel = TRUE,
                        numCores = 32,
                        quiet = FALSE)

# save ENMeval object

saveRDS(max20,paste0("ENMeval_",args[1],".rds"))