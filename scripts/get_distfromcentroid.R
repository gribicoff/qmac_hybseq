args <- commandArgs(trailingOnly=TRUE)

library(sf)

# proj4 string: "+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"
# IMPORTANT: the Q. alba, mac, bic shapefiles contain some internal gaps as polygons - must remove from

# read in shapefiles as sf objects, assign proj data/convert from polygon to multipolygon

spp <- c("mac","bic","alb","ste","mue")
sppshp <- paste0(spp,"_shp")
proj <- c("+proj=aea","+lon_0=-82","+ellps=clrk66","+lat_1=38","+lat_2=42","+lat_0=40","+units=m","+x_0=0","+y_0=0","+datum=WGS84","+axis=enu")
for (i in spp) {
  assign(paste0(i,"_shp"),st_read(paste0(i,"_little.shp")))
  assign(paste0(i,"_shp"),st_set_crs(get(paste0(i,"_shp")),"+proj=aea ++lon_0=-82 +ellps=clrk66 +lat_1=38 +lat_2=42 +lat_0=40 +units=m +x_0=0 +y_0=0 +datum=WGS84 +axis=enu"))
  assign(paste0(i,"_shp"),st_cast(get(paste0(i,"_shp")),"MULTIPOLYGON"))
}

bic_shp <- bic_shp[-12,] # removes Lake Winnebago
mac_shp <- mac_shp[-c(2,3,18),] # removes Lake Winnipeg, Lake of the Woods, Lake Winnebago respectively
alb_shp <- alb_shp[-c(7,13,19,20),] # removes the Adirondacks, Lake Winnebago, Lake Erie/Saint Clair, the Catskills respectively

# find centroid of each range, convert to WGS84 geodetic coordinates

for (i in spp) {
  assign(paste0(i,"_centroid"),st_centroid(st_union(st_geometry(get(paste0(i,"_shp"))))))
  assign(paste0(i,"_centroid"),st_transform(get(paste0(i,"_centroid")),crs=4326))
}

ctrds <- data.frame()
for (i in spp) {
  ctrds <- rbind(ctrds,(st_coordinates(get(paste0(i,"_centroid")))))
}
ctrds$spp <- paste0(spp,"_centroid")
colnames(ctrds) <- c("long","lat","spp")

# write table of centroids for manual spot-checking of distance calculations

write.table(ctrds[,c(3,1,2)],file="./little_range_centroids.txt",row.names=F,quote=F,sep="\t")

# read sample coordinates, add ordering (for sorting later)

pts <- read.csv(paste0(args,".csv"))
pts$order <- rownames(pts)

# calculate distance between points and respective centroids
all_pts <- data.frame()
for (i in spp) {
  b <- pts[grep(i,pts$ind_ID),]
  b_sf <- st_as_sf(b,coords=c("long","lat"),crs=4326)
  b_dists <- st_distance(b_sf,get(paste0(i,"_centroid")))
  b$dist_from_ctrd_m <- b_dists
  assign(paste0(i,"_pts"),b)
  assign(paste0(i,"_pts_sf"),b_sf)
  assign(paste0(i,"_dists"),b_dists)
  all_pts <- rbind(all_pts,b)
}

all_pts$order <- as.numeric(all_pts$order)
data.table::setorder(all_pts,order)

write.csv(all_pts[,-4],file=paste0(args,"_distfromcentroid.csv"),row.names=F,quote=F)
