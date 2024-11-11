# MODIFIED FROM SCRIPT PROVIDED BY USDA-NRCS SOIL SCIENTIST ANDY PAOLUCCI

# install packages
install.packages('remotes', repos="https://cloud.R-project.org")
install.packages('aqp', repos="https://cloud.R-project.org")
install.packages('sp', repos="https://cloud.R-project.org")
install.packages('terra', repos='https://rspatial.R-universe.dev')
install.packages('rgdal', repos="https://cloud.R-project.org")
install.packages('raster', repos="https://cloud.R-project.org")
install.packages('rgeos', repos="https://cloud.R-project.org")
install.packages('dplyr', repos="https://cloud.R-project.org")

# IMPORTANT: MAKE SURE YOU ARE USING CORRECT VERSION OF soilDB
# might need to use development version of soilDB for .get_comonth_SDA function that summarizes flooding/ponding data
remotes::install_github("ncss-tech/soilDB", dependencies = FALSE, force = TRUE)
# if not use this
#install.packages("soilDB")

# load packages
library(aqp)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(dplyr)
library(soilDB)

# command line arguments: coord csv file path, working directory
# args = c("~/qmac/analyses/envdata/soildata/newcoords/new_soilcoords_noSHF.csv","~/qmac/analyses/envdata/soildata/newcoords/")


args <- commandArgs(trailingOnly=TRUE)

#### Load data and prepare data ####
#set working directory
setwd(args[2])

# Load csv with points
x <- read.csv(file=args[1],header=T)

#initalize spatial coordinates. lat and long should match the columns in the spreadsheet where the coordinates are stored
coordinates(x) <- ~ long + lat

#get/set spatial reference system
proj4string(x) <- '+proj=longlat +datum=WGS84 +no_defs'

# plot in map
# library(maps)
# map('state', xlim=c(-110, -70), ylim=c(30, 50))
# map.axes()
# plot(x, add=TRUE, pch=17, col='red')

# #convert to spatial points data frame for extracting frost free days later
x.spdf <- as(x, 'SpatialPointsDataFrame')

#Optional: write as shapefile so we can view data in ArcMap/ArcPro
#writeOGR(x.spdf, dsn='C:/geodata/projectdata', layer='SDA_Points', driver='ESRI Shapefile', overwrite_layer=TRUE)

# perform SDA query on collection of points
mu.data <- SDA_spatialQuery(x, what = 'geom')

# use local spatial intersection to link source + SDA data
x$mukey <- over(x, mu.data)$mukey

# list unique MUKEYS and remove NA value for querying SDA data
y <- unlist(list(na.omit(unique(x$mukey))))

# create dataframe for joining later
joining.df <- as.data.frame(x)

#### Horizon Properties ####
# Weighted average for 0-100 cm
# CaCO3 = caco3_r
# CEC-7 = cec7_r
# EC = ec_r
# ECEC = ecec_r
# gypsum = gypsum_r
# pH = ph1to1h2o_r
# % clay = claytotal_r
# % sand = sandtotal_r
# % silt = silttotal_r
# Bulk density - 1/3 bar = dbthirdbar_r
# Linear extensibility = lep_r
# Organic matter (%) = om_r
# Ksat = ksat_r
# AWC = awc_r
final.horizon.output <- get_SDA_property(c("caco3_r","cec7_r","ec_r", "ecec_r", "gypsum_r", "ph1to1h2o_r", "claytotal_r", "sandtotal_r", "silttotal_r", "dbthirdbar_r", "lep_r", "om_r", "ksat_r", "awc_r"),
                       method = "WEIGHTED AVERAGE",
                       top_depth = 0,
                       bottom_depth = 100,
                       include_minors = FALSE,
                       miscellaneous_areas = FALSE,
                       mukeys = y)

write.csv(final.horizon.output, "./WT_AVG_HorizonProperties0100_1_24_2023.csv")







#### Flooding/Ponding and Mapunit Data from MUAGGATT Table (drainage class, AWS 0-100cm) ####
# Drainage class dominant condition =  drclassdcd
# Available water storage  0-100cm weighted average = aws0100wta

# get component data from list of mukeys
comp.data <- get_component_from_SDA(WHERE=paste("mukey in (", paste(y, collapse =","),")"), duplicates = TRUE)

# remove desired columns for joining to output later
comp.data.sub <- comp.data[c("mukey", "cokey", "compname","comppct_r", "majcompflag")]

# get list of cokeys for running flooding/ponding frequency
comps <- as.character(comp.data.sub$cokey)

# run function using for-loop
comp.pond.data<-NULL
for (x in comps){
  tmp<-soilDB:::.get_comonth_from_SDA(cokey=x)
  comp.pond.data<-rbind(comp.pond.data, tmp)}

#inspect
head(comp.pond.data)

# merge flooding/ponding data in with component data
comp.pond.data <- merge(comp.data.sub, comp.pond.data, by.x='cokey', by.y='cokey')

# remove NA values
comp.pond.data<-comp.pond.data[!is.na(comp.pond.data$pondfreqcl),]

# set NA Component Percentages to zero
comp.pond.data["comppct_r"][is.na(comp.pond.data["comppct_r"])] <- 0

# dominant condition function.
dominant.condition <- function(dmu, attr) {
  # get needed vars
  mukey <- unique(dmu$mukey)
  d.attr <- dmu[[attr]]
  d.comppct <-  dmu[["comppct_r"]]

  # sum up comppct for each level of attr
  res <- aggregate(d.comppct, list(d.attr), sum)
  if(!nrow(res))
    res <- data.frame(a=character(0), b=character(0))

  # select the first occurence of maximum value
  #  might not be stable if two attr have same percentage
  suppressWarnings(res <- res[which(res$x == max(res$x, na.Rm=TRUE))[1],])
  names(res) <- c("dominant_condition", "dominant_condition_pct")

  # create data.frame output
  return(data.frame(mukey=mukey, res))}

# Calculate dominant condition for all mukeys
# # use mukey, a site-level attribute, to split components into mapunits
pond.mukey <- split(comp.pond.data, comp.pond.data$mukey)

# run dominant condition function: ponding frequency
pond.dom.condition <- do.call('rbind', lapply(pond.mukey, dominant.condition, "pondfreqcl"))

# set up col/row names
head(pond.dom.condition)
colnames(pond.dom.condition) <- c("mukey","pondfreqcl_dc","pondfreqcl_dc_pct")
rownames(pond.dom.condition) <- NULL

# run dominant condition function: flooding frequency
flood.dom.condition <- do.call('rbind', lapply(pond.mukey, dominant.condition, "flodfreqcl"))

# set up col/row names
head(flood.dom.condition)
colnames(flood.dom.condition) <- c("mukey","flodfreqcl_dc","flodfreqcl_dc_pct")
rownames(flood.dom.condition) <- NULL

# merge two dataframes
flood.pond.freq.dc <- merge(pond.dom.condition, flood.dom.condition, by.x='mukey', by.y='mukey')

# get muaggatt data to compare flooding frequency
muaggatt <- get_SDA_muaggatt(mukeys=y)

# filter desired columns
muaggatt.output <- muaggatt[c("mukey", "drclassdcd", "aws0100wta","flodfreqdcd", "pondfreqprs")]

# merge muaggatt data with flooding and ponding data
mapunit.data.output <- merge(flood.pond.freq.dc, muaggatt.output, by.x='mukey', by.y='mukey')
head(mapunit.data.output)

# check to make sure flodfreqdcd from muaggatt table matches flodfreqcl_dc in output
identical(mapunit.data.output$flodfreqdcd, mapunit.data.output$flodfreqcl_dc)
table(mapunit.data.output$flodfreqdcd, mapunit.data.output$flodfreqcl_dc)

# figure out which map units dont have matching values
# replace NA values with 'missing'
mapunit.data.output[is.na(mapunit.data.output)] <- "missing.data"

# create new column and add if values match Yes/No
mapunit.data.output$check <- ifelse(mapunit.data.output$flodfreqdcd==mapunit.data.output$flodfreqcl_dc,"Yes","No")

# subset records that dont match. Manually check to see why they don't match. From 1/24/2023 data set it looks like missmatches are due to erros with muaggatt data and not SDA flooding/ponding calculation
dont.match <-subset(mapunit.data.output, check =='No')

# subset mukey and flooding and ponding frequency dominant condition columns
mapunit.data.output <- mapunit.data.output[c("mukey", "pondfreqcl_dc", "flodfreqcl_dc", "drclassdcd", "aws0100wta")]

write.csv(mapunit.data.output, "./Flood.Pond.Drainage.AWS100cm_from_GPS_1_24_2023.csv")





#### Frost Free Days ####
library(raster)
# load rast file with frost free dats
ffd <- raster("./ffd_90_pct_800m.tif")

# plot raster
# plot(ffd, main="Frost Free Days")

#extract values
final.ffd.output <- raster::extract(ffd,     # raster layer
                                    x.spdf,   # Spatial points dataframe
                                    df=TRUE)  # return a dataframe?

# grab the point IDs for merging later
final.ffd.output$ind_ID <- x.spdf$ind_ID

final.ffd.output <- final.ffd.output[c('ffd_90_pct_800m', 'ind_ID')]

# rename columns
colnames(final.ffd.output) <- c('ffd', 'ind_ID')

# write to csv
write.csv(final.ffd.output, './PRISM_FFD_Extract_1_24_2023.csv')







#### Surface Texture ####
# get component data for mukeys in list
spc <- fetchSDA(WHERE=paste("mukey in (", paste(y, collapse =","),")"), duplicates = TRUE)

# create empty dataframe for storing results
output <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(output) <- c("mukey", "compname", "comppct_r", "cokey", "texture")

# create loop for pulling surface texture from each component in mukey
for (k in unique(spc$cokey)) {
  tmp <-subset(spc, spc@site$cokey==k)
  output[nrow(output)+1,]= c(tmp@site$mukey, tmp@site$compname, tmp@site$comppct_r, tmp@site$cokey, tmp@horizons$texture[1])}

# add surface texture to site level attribute
spc@site$surf_texture <- output$texture

# Calculate dominant condition for all mukeys
# # use mukey, a site-level attribute, to split components into mapunits
f.mukey <- aqp::split(spc, 'mukey')
# run dominant condition function
final.surface.texture <- do.call('rbind', lapply(f.mukey, dominant.condition, "surf_texture"))

colnames(final.surface.texture) <-c('mukey', 'surf_texture', 'pct')
rownames(final.surface.texture) <- NULL
surface.texture.df <- final.surface.texture[c('mukey', 'surf_texture')]





#### First Mineral Horizon Texture ####
# make a new column to store generalized reskind, and fill with NA
spc@horizons$min_text <- rep(NA, times=length(spc@horizons$texture))

# populate new column with texture
spc@horizons$min_text <- spc@horizons$texture
# populate organic textures as NA
spc@horizons$min_text[spc@horizons$texture %in% c('muck','mpt','peat','hpm','mpm','spm')] <- NA

#remove NA rows and put back into SPC
spc@horizons <- spc@horizons[!is.na(spc@horizons$min_text),]

# create empty dataframe for storing results
output2 <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(output2) <- c("mukey", "compname", "comppct_r", "cokey", "texture")

# create loop for pulling surface texture from each component in mukey
for (k in unique(spc$cokey)) {
  tmp <-subset(spc, spc@site$cokey==k)
  output2[nrow(output2)+1,]= c(tmp@site$mukey, tmp@site$compname, tmp@site$comppct_r, tmp@site$cokey, tmp@horizons$texture[1])}

# remove unwanted columns
output2 <- output2[c('cokey','texture')]

# add surface texture to site level attribute
spc@site <- merge(spc@site, output2, by.x='cokey', by.y='cokey')

# Calculate dominant condition for all mukeys
# # use mukey, a site-level attribute, to split components into mapunits
f.mukey <- aqp::split(spc, 'mukey')
# run dominant condition function
final.mineral.surface.texture <- do.call('rbind', lapply(f.mukey, dominant.condition, "texture"))
colnames(final.mineral.surface.texture) <- c('mukey', 'min_surf_texture', 'pct')
rownames(final.mineral.surface.texture) <- NULL

mineral.surface.texture.df <- final.mineral.surface.texture[c('mukey', 'min_surf_texture')]
both.surface.textures <- merge(surface.texture.df, mineral.surface.texture.df, by.x='mukey', by.y='mukey')

write.csv(both.surface.textures, './DOM_CD_SurfaceTexture_1_24_2023.csv')






####Total Available Water supply####
# sum(AWC X horizon thickness) for each horizon = total available water supply
# get component data for mukeys in list
spc.aws <- fetchSDA(WHERE=paste("mukey in (", paste(y, collapse =","),")"), duplicates = TRUE)

# calculate horizon thickness
spc.aws@horizons$hzdepb_r - spc.aws@horizons$hzdept_r -> spc.aws@horizons$hzthick

# calculate horizon Available water supply by multiplying thickness by aw
spc.aws@horizons$aws <- spc.aws@horizons$hzthick*spc.aws@horizons$awc_r

# set NA values to 0 <- might try removing NA values
spc.aws@horizons$aws <- spc.aws@horizons$aws %>% replace(is.na(.), 0)

# create empty dataframe for storing results
output.aws <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(output.aws) <- c("cokey", "aws")

# for each cokey sum AWS and add to new dataframe
for (k in unique(spc.aws$cokey)) {
  tmp <-subset(spc.aws, spc.aws@site$cokey==k)
  output.aws[nrow(output.aws)+1,]= c(tmp@site$cokey, sum(tmp@horizons$aws))}

# Replace 0 AWC values to NA since those are either bedrock components or components with missing data
output.aws["aws"][output.aws["aws"] == 0] <- NA

#inspect
print(output.aws)

#rename columns so aws isnt both horizon and site property
colnames(output.aws) <- c("cokey", "total_aws")

# join AWS to site level of soil profile collection
aqp::site(spc.aws) <- output.aws

# extract the site data
aws.df <- aqp::site(spc.aws)

# filter desired columns
aws.df <- aws.df[c("mukey", "cokey", "comppct_r","total_aws")]

# omit components with 'NA' values
aws.df <- na.omit(aws.df)

# convert total aws to numeric
aws.df$total_aws <- as.numeric(aws.df$total_aws)

# calculate weighted mean for each map unit
aws.summary <- as.data.frame(sapply(split(aws.df, aws.df$mukey), function(x) weighted.mean(x$total_aws, w = x$comppct_r)))

# add row ids to new column
mukey <- row.names(aws.summary)
rownames(aws.summary) <- NULL
final.aws.output <- cbind(mukey, aws.summary)

#rename columns
colnames(final.aws.output) <- c("mukey", "total_aws")

# round value to two decimal places
final.aws.output$total_aws<- round(final.aws.output$total_aws, 2)

# write to csv
write.csv(final.aws.output, './WT_AVG_AWS_WholeProfile_1_24_2023.csv')




#### Depth to Restrictive Layers ####
# get component data from list of mukeys
spc.Res <- fetchSDA(WHERE=paste("mukey in (", paste(y, collapse =","),")"), duplicates = TRUE)


# subset major components since minors likely are missing data
spc.Res <- subset(spc.Res, majcompflag=='Yes')
str(spc.Res)

# list unique diagnostic features
unique(spc.Res@diagnostic$featkind)
# list unique restriction kinds
unique(spc.Res@restrictions$reskind)
# list unique horizon designations
unique(spc.Res@horizons$hzname)
# list unique horizon textures
unique(spc.Res@horizons$texture)

# remove MISC Components
spc.Res<- subset(spc.Res, spc.Res@site$compkind != 'Miscellaneous area')

# Get depth to any component restrictions
# Generalize component restriction kinds
#Lithic = 'Lithic bedrock'
#Paralithic = 'Paralithic bedrock'
#Densic = 'Densic material' & 'Densic bedrock'
#Texture = 'Strongly contrasting textural stratification' & 'Abrupt textural change'
#Fragipan = 'Fragipan'
#Cemented  ='Manufactured layer'&'Petrocalcic'

   ####..Restrictive Features ####
# get restrictions
restrict.df<- spc.Res@restrictions
head(restrict.df)

# filter out columns: kind, top depth, cokey
restrict.df <- restrict.df[c("reskind", "resdept_r", "cokey")]

# make a new column to store generalized reskind, and fill with NA
restrict.df$gen_reskind <- rep(NA, times=length(restrict.df$reskind))

# generalize component restriction kinds
restrict.df$gen_reskind[restrict.df$reskind %in% c('Lithic bedrock')] <- 'Lithic'
restrict.df$gen_reskind[restrict.df$reskind %in% c('Paralithic bedrock')] <- 'Paralithic'
restrict.df$gen_reskind[restrict.df$reskind %in% c('Densic material', 'Densic bedrock')] <- 'Densic'
restrict.df$gen_reskind[restrict.df$reskind %in% c('Strongly contrasting textural stratification','Abrupt textural change')] <- 'Texture'
restrict.df$gen_reskind[restrict.df$reskind %in% c('Fragipan')] <- 'Fragipan'
restrict.df$gen_reskind[restrict.df$reskind %in% c('Manufactured layer','Petrocalcic')] <- 'Cemented'

# Get depth to restriction closest to surface for each unique cokey

# use cokey to split restrictions by component
restrict.cokey <- split(restrict.df, restrict.df$cokey)

# Function to get diagnostic feature with lower depth from surface horizon
min.Restriction <- function(x) {
  res <- x[which.min(x$resdept_r),]
  return(data.frame(res))}

# run function on all components
min.Restrict.df <- do.call('rbind', lapply(restrict.cokey, min.Restriction))
rownames(min.Restrict.df)<-NULL

colnames(min.Restrict.df) <-c('reskind', 'res_dept', 'cokey', 'gen_reskind')

# order and subset wanted columns
restriction.output <- min.Restrict.df[c("cokey","gen_reskind","res_dept")]

   ####..Diagnostic Features ####
# Get depth to any diagnostic features that are restrictive
# get diagnostics
# get restrictions
diagnostic.df<- spc.Res@diagnostic
head(diagnostic.df)

# filter out columns: kind, top depth, cokey
diagnostic.df <- diagnostic.df[c("featkind", "featdept_r", "cokey")]

# make a new column to store generalized featkind, and fill with NA
diagnostic.df$gen_featkind <- rep(NA, times=length(diagnostic.df))

# Generalize component diagnostic feature kinds
#Lithic = 'Lithic contact'
#Paralithic = 'Paralithic contact
#Densic = Densic contact
#Texture = 'Abrupt textural change','Strongly contrasting particle size class' <- tried this and got some soils with abrupt clay increases that didnt appear root restrictive
#Fragipan = 'Frangipan'
#Cemented = = 'Petrocalcic horizon'
#diagnostic.df$gen_featkind[diagnostic.df$featkind %in% c('Abrupt textural change', 'Strongly contrasting particle size class')] <- 'Texture'
diagnostic.df$gen_featkind[diagnostic.df$featkind %in% c('Paralithic contact')] <- 'Paralithic'
diagnostic.df$gen_featkind[diagnostic.df$featkind %in% c('Densic contact')] <- 'Densic'
diagnostic.df$gen_featkind[diagnostic.df$featkind %in% c('Lithic contact')] <- 'Lithic'
diagnostic.df$gen_featkind[restrict.df$featkind %in% c('Fragipan')] <- 'Fragipan'
diagnostic.df$gen_featkind[restrict.df$featkind %in% c('Petrocalcic horizon')] <- 'Cemented'

# subset rows that have generalized diagnostic populated
diagnostic.df  <- diagnostic.df[!is.na(diagnostic.df$gen_featkind), ]
nrow(diagnostic.df)

# # # use cokey to split diagnostics by component
diagnostic.cokey <- split(diagnostic.df, diagnostic.df$cokey)

# Function to get diagnostic feature with lower depth from surface horizon
min.diagnostic <- function(x) {
  res <- x[which.min(x$featdept_r),]
  return(data.frame(res))}

# run function on all components
min.diagnostic.df <- do.call('rbind', lapply(diagnostic.cokey, min.diagnostic))
rownames(min.diagnostic.df)<-NULL

# rename columns
colnames(min.diagnostic.df) <- c('featkind', 'diag_dept', 'cokey', 'gen_diag')
# reorder and subset wanted columns
diagnostic.output <- min.diagnostic.df[c("cokey","gen_diag","diag_dept")]

   ####..Horizon Designation ####
# Get depth to horizons that are root restrictive
# Cd: densic
# R: lithic
# Cr: paralithic
# Btx: fragipan
# Bkm/Bkkm: petrocalcic
# Generalize to match restrictions/diagnostic features
# make a new column to store generalized featkind, and fill with NA
spc.Res@horizons$gen_hzname<- rep(NA, times=length(spc.Res@horizons$hzname))
#Lithic = 'R','2R','3R','3Rt'
#Paralithic = 'Cr','2Cr'
#Densic = 'Cd','2Cd'
#Fragipan = 'Btx','2Btx'
#Cemented = = 'Bkm','Bkkm', Bkkm1','Bkkm2'

spc.Res@horizons$gen_hzname[spc.Res@horizons$hzname %in% c('R','2R','3R', '3Rt', 'R4')] <- 'Lithic'
spc.Res@horizons$gen_hzname[spc.Res@horizons$hzname %in% c('Cr','2Cr')] <- 'Paralithic'
spc.Res@horizons$gen_hzname[spc.Res@horizons$hzname %in% c('Cd','2Cd')] <- 'Densic'
spc.Res@horizons$gen_hzname[spc.Res@horizons$hzname %in% c('Btx','2Btx')] <- 'Fragipan'
spc.Res@horizons$gen_hzname[spc.Res@horizons$hzname %in% c('Bkm','Bkkm','Bkkm1','Bkkm2')] <- 'Cemented'

# isolate component horizon table
horizons.df <- spc.Res@horizons

# subset horizons that have generalized label for restriction kind
horizons.df<- horizons.df[!is.na(horizons.df$gen_hzname), ]

# filter wanted columns
horizons.df <- horizons.df[c("cokey", "gen_hzname", "hzdept_r")]

# use cokey to split restrictive horizons by component
horizons.cokey <- split(horizons.df, horizons.df$cokey)

# Function to get restrictive horizon with lower depth from surface horizon
min.horizon <- function(x) {
  res <- x[which.min(x$hzdept_r),]
  return(data.frame(res))}

# run function on all components
min.horizons.df <- do.call('rbind', lapply(horizons.cokey, min.horizon))
rownames(min.horizons.df)<-NULL

# rename columns
colnames(min.horizons.df) <- c('cokey', 'gen_desg', 'desg_dept')
# reorder and subset wanted columns
designation.output <- min.horizons.df[c('cokey', 'gen_desg', 'desg_dept')]

   ####..Texture ####
# Get depth to textures that are root restrictive
# wb: weathered bedrock
# uwb: unweathered bedrock
# cem-mat: cemented material
# br: bedrock
#
# make a new column to store generalized texture, and fill with NA
spc.Res@horizons$gen_texture<- rep(NA, times=length(spc.Res@horizons$texture))

# Generalize to match restrictions/diagnostic features
# Lithic = 'uwb' & 'br'
# Paralithic = 'wb'
# Cemented  = 'cem-mat'
spc.Res@horizons$gen_texture[spc.Res@horizons$texture %in% c('uwb')] <- 'Lithic'
spc.Res@horizons$gen_texture[spc.Res@horizons$texture %in% c('wb')] <- 'Paralithic'
spc.Res@horizons$gen_texture[spc.Res@horizons$texture %in% c('br')] <- 'Paralithic/Lithic'
spc.Res@horizons$gen_texture[spc.Res@horizons$texture %in% c('cem-mat')] <- 'Cemented'

# isolate component horizon table
texture.df <- spc.Res@horizons

# subset horizons that have generalized label for restriction kind
texture.df <- texture.df[!is.na(texture.df$gen_texture), ]

# filter wanted columns
texture.df <- texture.df[c("cokey", "gen_texture", "hzdept_r")]

# use cokey to split restrictive horizons by component
texture.cokey <- split(texture.df, texture.df$cokey)

# Function to get restrictive horizon with lower depth from surface horizon
min.texture <- function(x) {
  res <- x[which.min(x$hzdept_r),]
  return(data.frame(res))}

# run function on all components
min.texture.df <- do.call('rbind', lapply(texture.cokey, min.texture))
rownames(min.texture.df)<-NULL

# rename columns
colnames(min.texture.df) <- c('cokey', 'gen_texcl', 'texcl_dept')
# reorder and subset wanted columns
texture.output <- min.texture.df[c('cokey', 'gen_texcl', 'texcl_dept')]

   ####..Analyze Restriction Results####
#put all data frames into list
df_list <- list(restriction.output,
                diagnostic.output,
                designation.output,
                texture.output)

#merge all data frames in list
final.output <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

# Make new column to store soil depth cm
final.output$soildepth <- rep(NA, times=length(final.output$cokey))

# calculate minimum value between 4 methods of determining soil depth (restriction, diagnostic, horizon designation, texture). Ignore NA values.
final.output$soildepth<- pmin(final.output$res_dept, final.output$diag_dept, final.output$desg_dept, final.output$texcl_dept, na.Rm=TRUE)

# take cokey and soil depth columns
final.output<- final.output[c('cokey','soildepth')]

# add soil depth to component site data
aqp::site(spc.Res) <- final.output

# subset components that have NA for soil depth
depth.no.na <- subset(spc.Res, is.na(spc.Res@site$soildepth))

# isolate horizon data.
depth.horizons <- depth.no.na@horizons

#check: this shouldn't contain any restrictive horizons such as Cr, Cd, Btx, R, Bkkm
unique(depth.horizons$hzname)

# filter wanted columns
depth.horizons <- depth.horizons[c("cokey", "hzdepb_r")]

# use cokey to split horizons by component
depth.horizons.cokey <- split(depth.horizons, depth.horizons$cokey)

# create function to pick highest soil depth hzdepb_r value (deepest horizon lower boundary) and assign to each cokey
max.depth <- function(x) {
  res <- x[which.max(x$hzdepb_r),]
  return(data.frame(res))}

# apply function to components
max.depth.df <- do.call('rbind', lapply(depth.horizons.cokey, max.depth))
rownames(max.depth.df)<-NULL

# change column names
colnames(max.depth.df)<-c('cokey', 'soildepth')

# add to final.output
soildepth.df <- rbind(final.output, max.depth.df)

# get fresh set of component data  from original comp.data.sda file
# subset major components since minors likely are missing data
final.spc <- fetchSDA(WHERE=paste("mukey in (", paste(y, collapse =","),")"), duplicates = TRUE)

# subset major components
final.spc <- subset(final.spc, majcompflag=='Yes')

# remove MISC Components
final.spc <- subset(final.spc, final.spc@site$compkind != 'Miscellaneous area')

# add soil depth data as site level attribute
aqp::site(final.spc) <- soildepth.df

# calculate weighted average for each mukey
# extract the site data
final.spc.site <- aqp::site(final.spc)

# filter desired columns
final.spc.site <- final.spc.site [c("mukey", "cokey", "comppct_r","soildepth")]

# omit components with 'NA' values. Shouldnt have any but good check
final.spc.site <- na.omit(final.spc.site)

# convert soildepth to numeric
final.spc.site$soildepth <- as.numeric(final.spc.site$soildepth)

# calculate weighted mean for each map unit
depth.summary <- as.data.frame(sapply(split(final.spc.site, final.spc.site$mukey), function(x) weighted.mean(x$soildepth, w = x$comppct_r)))

# add row ids to new column
mukey <- row.names(depth.summary)
rownames(depth.summary) <- NULL
final.depth.output <- cbind(mukey, depth.summary)

#rename columns
colnames(final.depth.output) <- c("mukey", "avg_soildepth")
#change avg_soil depth column to integer
final.depth.output$avg_soildepth <- as.integer(final.depth.output$avg_soildepth)

# write csv with results
write.csv(final.depth.output, './WT_AVG_soildepth_1_24_2023.csv')

#### Export Results ####
# Inspect results from each section of the script
head(final.horizon.output)
head(mapunit.data.output)
head(both.surface.textures)
head(final.aws.output)
head(final.depth.output)

# join joining dataframe with coordinates and results from frost free day extraction
points.df <- cbind(final.ffd.output, joining.df)
# grab/order columns
points.df <- points.df[,c(2,4,5,6,1)]

# merge files together
z1 <- merge(final.horizon.output, mapunit.data.output, by.x='mukey', by.y='mukey', all=TRUE)
z2 <- merge(z1, both.surface.textures, by='mukey', all=TRUE)
z3 <- merge(z2, final.aws.output, by='mukey', all=TRUE)
complete.dataset <- merge(z3, final.depth.output, by='mukey', all=TRUE)

# merge point and result dataset
final.Results <- merge(points.df, complete.dataset, by.x='mukey', all=TRUE)

# replace missing.data with NA in final results
# replace missing.data in drainage class and available water storage with NA
final.Results$drclassdcd[final.Results$drclassdcd %in% c('missing.data')] <- NA
final.Results$aws0100wta[final.Results$aws0100wta %in% c('missing.data')] <- NA

# replace commas with semicolons to prevent issues with csv formatting

final.Results$muname <- gsub(",",";",final.Results$muname)

# write csv with final results
write.csv(final.Results, './FinalResults_SDA_from_GPS.csv',quote=F,row.names=F)
