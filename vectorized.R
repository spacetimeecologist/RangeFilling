require(rgdal)
require(maptools)
require(rgeos)
require(plyr)
require(dismo)
require(foreach)
require(doParallel)

# Sets extent and resolution of analysis, creates bounding box polygon called bbox
ext <- extent(-130, -50, 25, 60)
res <- 5
pts <- coordinates
pts <- rbind(c(xmin(ext), ymin(ext)), c(xmin(ext), ymax(ext)),
             c(xmax(ext), ymax(ext)), c(xmax(ext), ymin(ext)))
pts <- rbind(pts, pts[1,])
bbox <- SpatialPolygons(list(Polygons(list(Polygon(pts)), 1)))
crs(bbox) <- crs('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 
                  +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')
bbox <- spTransform(bbox, CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 
                              +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m 
                              +no_defs'))

# Creates land shapefile from mexco, canada and united states
cntry <- readShapePoly('C:\\Users\\ben\\Documents\\Research\\RangeFilling\\GIS\\50m_cultural\\ne_50m_admin_0_countries.shp')
state <- readShapePoly('C:\\Users\\ben\\Documents\\Research\\RangeFilling\\GIS\\50m_cultural\\ne_50m_admin_1_states_provinces_lakes.shp')
mex <- cntry[which(cntry$sovereignt=='Mexico'),]
can <- state[which(state$sr_adm0_a3=='CAN'),]
usa <- state[which(state$sr_adm0_a3=='USA'),]
land <- mex+can+usa
crs(land) <- crs(bbox)
land <- crop(land, bbox)
remove(list=c('mex', 'can', 'usa'), envir=.GlobalEnv)

# gets environmental layers for SDMs using bbox, resolution, and land. West=T includes precip seasonality
envi <- getEnvi(bbox, res, land, west=T)
crs(envi) <- crs(bbox)

littlelist <- read.csv('C:\\Users\\ben\\Documents\\Research\\RangeFilling\\Data\\USGS_TreeAtlas\\USGS_littlelist.csv', stringsAsFactors=FALSE)

bientraits <- read.csv('C:\\Users\\ben\\Documents\\Research\\RangeFilling\\Data\\BIEN\\bien2_traits_complete_2012-11-22.txt\\bien2_traits_complete_2012-11-22.txt', sep='\t', stringsAsFactors=FALSE) 

trytraits <- read.csv(file='C:/Users/ben/Documents/Data/TRY/try473_07012015103839/473.txt',
                      sep='\t', stringsAsFactors=FALSE)


getlist <- function(list){
  
  colnames(list) <- c('Sciname', 'Height', 'GrowthForm', 'Tree')
  
  list$Tree <- NA
  list$Height <- lapply(list$Sciname, FUN=bien.height)
  list$GrowthForm <- lapply(list$Sciname, FUN=bien.form)
  list$Height <- as.numeric(list$Height)
  list$GrowthForm <- as.character(list$GrowthForm)
  list$Height[which(is.na(list$Height)==TRUE)] <- 0
  list$GrowthForm[which(list$GrowthForm=='character(0)')] <- 'null'
  
  list$Tree[which(list$Height>=10)] <- TRUE
  list$Tree[which(list$GrowthForm=='tree')] <- TRUE
  
  list <- list[which(list$Tree==TRUE),]
  return(list)
}

treelist <- getlist(littlelist)
treelist <- treelist[which(treelist$Tree==TRUE),]

# This block of code gets the path names for the UGGS tree ranges
shp_dirs <- list.dirs(path='C:/Users/ben/Documents/Research/RangeFilling/Data/USGS_TreeAtlas')
shp_dirs <- shp_dirs[-1]
shp_paths <- list.files(path=shp_dirs, pattern='.shp')
shp_paths.full <- paste(shp_dirs, shp_paths, sep='/')
spcodes <- list()
spfull <- list()
i <- 1
for(i in i:nrow(treelist)){
  Sciname <- treelist$Sciname[i]
  genus <- unlist(strsplit(Sciname, split=' '))[1]
  species <- unlist(strsplit(Sciname, split=' '))[2]
  spcodes[i] <- paste(substr(tolower(genus), 0, 4), substr(species, 0, 4), sep='')
  spfull[i] <- paste(genus, species, sep='')
}
matches <- match(spcodes, substr(shp_paths.full[], start=75, stop=82))
sp_paths.full <- shp_paths.full[matches]

# Creates rangeslist, a list of shapefiles of the USGS tree atlas ranges
rangeslist <- llply(sp_paths.full, readShapePoly, proj4string=CRS("+proj=aea 
+lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"),
                    delete_null_obj=TRUE, .progress="text")
names(rangeslist) <- spfull

# crops ranges list by bbox
ranges.crop <- llply(rangeslist, gIntersection,
                     spgeom2=bbox, byid=c(T,F), drop_not_poly=T, .progress="text")

# this block of code finds which ranges were outside of bbox, then removes them from ranges
j <- 1
nulls <- list()
for(j in j:length(ranges.crop)){
  nulls[j] <- is.null(ranges.crop[[j]])
}
ranges.crop <- ranges.crop[which(nulls==FALSE)]
rangeslist <- rangeslist[which(nulls==FALSE)]

# projects cropped ranges into equal area projection, and take centroids of total ranges
ranges.proj <- llply(ranges.crop, spTransform, 
                     CRSobj=CRS("+proj=aea +datum=WGS84 +units=m"), .progress="text")
r.centroids <- llply(rangeslist, gCentroid, .progress="text")

# this block of code removes species with ranges less than 10000km2 and those whose centroids
# are outside of extent
R.km2 <- ldply(ranges.proj, gArea)
R.coords <- ldply(r.centroids, coordinates)
R.lon <- data.frame(R.coords$x)
R.lat <- data.frame(R.coords$y)
remove(R.coords)
R.in <- ldply(.data=c(1:length(r.centroids)), 
              .fun=function(X){R.lat[X,]<ymax(ext)&&
                                 ymin(ext)<R.lat[X,]&&
                                 R.lon[X,]<xmax(ext)&&
                                 xmin(ext)<R.lon[X,]})
m <- 1 
R.pass <- data.frame(1)
for(m in m:nrow(R.in)){
  if(R.in[m,]==T&&(R.km2$V1[m]>10000)==T){
    R.pass[m,] <- TRUE
  }
  else{R.pass[m,]<-FALSE}
}

R.km2 <- data.frame(R.km2[which(R.pass==T),2], row.names=(R.km2[which(R.pass==T), 1]))
R.lon <- data.frame(R.lon[which(R.pass==T),], row.names=rownames(R.km2))
R.lat <- data.frame(R.lat[which(R.pass==T),], row.names=rownames(R.km2))

ranges.crop <- ranges.crop[which(R.pass==T)]



# pres is a list of presence points for every species
pres <- llply(.data=c(1:length(ranges.crop)), 
              .fun=function(X){X=names(ranges.crop[X])
                               set.seed(1991)
                               spsample(ranges.crop[[X]],
                                        n=sqrt(R.km2[row.names(R.km2)==X,1]),
                                        type='regular')},
              .progress='text')
names(pres) <- names(ranges.crop)

# ranges.buff takes a 100km buffer surrounding a species range which will be used for psuedoabsences
ranges.buff <- llply(.data=c(1:length(ranges.crop)), 
                     .fun=function(X){X=names(ranges.crop[X])
                                      aggregate(gBuffer(ranges.crop[[X]],
                                                        byid=T, width=1))},
                     .progress='text')
names(ranges.buff) <- names(ranges.crop)

# not ranges is the land area which is outside of ranges.buff
not.ranges <- llply(.data=c(1:length(ranges.crop)), 
                    .fun=function(X){X=names(ranges.buff[X])
                                     land-ranges.buff[[X]]},
                    .progress='text')
names(not.ranges) <- names(ranges.buff)

workers <- makeCluster(4)
registerDoParallel(workers)
clusterEvalQ(workers, c(library(sp), library(rgdal)))
clusterExport(workers, varlist=c('ranges.crop', 'R.km2', 'land'
), envir=.GlobalEnv)

p.abs <- llply(.data=c(1:length(ranges.crop)), 
               .fun=function(X){X=names(ranges.crop[X])
                                set.seed(1991)
                                spsample(not.ranges[[X]], 
                                         n=100+(100000/sqrt(R.km2[row.names(R.km2)==X,1])), 
                                         type='regular')},
               .progress='text')
names(p.abs) <- names(ranges.crop)
stopCluster(workers)


set.seed(1991)
group <- llply(pres, kfold, k=5)



k <- 1  
for (k in k:5) {
  print(k)
  train <- llply(.data=c(1:length(pres)), 
                 .fun=function(X){X=names(pres[X])
                                  pres[[X]][which(group[[X]]!=k)]})
  test <- llply(.data=c(1:length(pres)),
                .fun=function(X){X=names(pres[X]) 
                                 pres[[X]][which(group[[X]]==k)]})
  names(train) <- names(pres)
  names(test) <- names(pres)
  assign(paste('bc.p.', k, sep=''), 
         value=climateEnvelope2(test, train, p.abs, pres, envi))
}

# k <- 1  
# for (k in k:5) {
#   print(k)
#   train <- pres.mac[which(group!=k)]
#   test <- pres.mac[which(group==k)]
#   names(train) <- names(pres.mac)
#   names(test) <- names(pres.mac)
#   assign(paste('bc.p.', k, sep=''), 
#          value=climateEnvelope2(test, train, mac.abs, pres.mac, envi))
# }
# 


bc.p.sum <- llply(.data=c(1:length(bc.p.1)), .fun=function(X){X=names(bc.p.1[X])
                                                              sum(bc.p.1[[X]], 
                                                                  bc.p.2[[X]], 
                                                                  bc.p.3[[X]], 
                                                                  bc.p.4[[X]], 
                                                                  bc.p.5[[X]])},
                  .progress='text')
names(bc.p.sum) <- names(bc.p.1)
bc.p.sum <- llply(.data=c(1:length(bc.p.sum)), .fun=function(X){X=names(bc.p.sum[X])
                                                                bc.p.sum[[X]] > 4},
                  .progress='text')
names(bc.p.sum) <- names(bc.p.1)

bc.p.area <- bc.p.sum
names(bc.p.area) <- names(bc.p.sum)
bc.p.area <- llply(.data=c(1:length(bc.p.area)), 
                   .fun=function(X){X=names(bc.p.area[X])
                                    bc.p.area[[X]][bc.p.area[[X]]==FALSE] <- NA
                                    return(bc.p.area[[X]])},
                   .progress='text')
names(bc.p.area) <- names(bc.p.sum)

P.km2 <- ldply(.data=c(1:length(bc.p.area)), 
               .fun=function(X){X=names(bc.p.area[X])
                                cellarea <- area(bc.p.area[[X]], na.rm=TRUE, weights=F)
                                cellarea <- sum(cellarea[!is.na(cellarea)])
                                return(cellarea)},
               .progress='text')
row.names(P.km2) <- names(bc.p.area)

p.centroids <- llply(.data=c(1:length(bc.p.area)),
                     .fun=function(X){X=names(bc.p.area[X])
                                      cells <- xyFromCell(bc.p.area[[X]],
                                                          cell=which(values(bc.p.area[[X]])==T),
                                                          spatial=T)                
                                      centroid <- gCentroid(cells)
                                      return(centroid)},
                     .progress='text')
names(p.centroids) <- names(bc.p.area)

P.coords <- ldply(p.centroids, coordinates)
P.lon <- data.frame(P.coords$x, row.names=names(p.centroids))
P.lat <- data.frame(P.coords$y, row.names=names(p.centroids))
remove(P.coords)

P.shape <- ldply(.data=c(1:length(bc.p.area)),
                 .fun=function(X){X=names(bc.p.area[X])
                                  cells <- xyFromCell(bc.p.area[[X]],
                                                      cell=which(values(bc.p.area[[X]])==T),
                                                      spatial=T)                
                                  xlength <- xmax(extent(cells))-xmin(extent(cells))
                                  ylength <- ymax(extent(cells))-ymin(extent(cells))
                                  shape <- ylength/xlength
                                  return(shape)},
                 .progress='text')
row.names(P.shape) <- names(bc.p.area)

R.shape <- ldply(.data=c(1:length(ranges.crop)),
                 .fun=function(X){X=names(ranges.crop[X])              
                                  xlength <- xmax(extent(ranges.crop[[X]]))
                                  -xmin(extent(ranges.crop[[X]]))
                                  ylength <- ymax(extent(ranges.crop[[X]]))
                                  -ymin(extent(ranges.crop[[X]]))
                                  shape <- ylength/xlength
                                  return(shape)},
                 .progress='text')
row.names(R.shape) <- names(ranges.crop)


df <- data.frame(row.names=row.names(R.km2), R.km2, P.km2, (R.km2/P.km2)*100,
                 R.lat, R.lon, P.lat, P.lon, R.shape, P.shape, R.lon<(-100), (R.lon>(-84)&R.lat<(28)),
                 check.rows=T, stringsAsFactors=F)
colnames(df) <- c('R.km2', 'P.km2', 'RP', 'R.lat', 'R.lon', 'P.lat', 'P.lon',
                  'R.shape', 'P.shape', 'West', 'Everglades')

k <- 1
for(k in k:length(unique(trytraits$TraitName))){
  traitname <- unique(trytraits$TraitName)[k]
  if(nchar(traitname)>0){
    vectorname <- gsub(' ', '.', traitname)
    print(vectorname)
    assign(vectorname, 
           value=ldply(.data=c(1:length(ranges.crop)),
                       .fun=function(X){X=names(ranges.crop)[X]
                                        sptrait <- traits.try(X,
                                                              traitname,
                                                              trytraits)
                                        return(sptrait)},
                       .progress='text')
    )
  }
  k <- k + 1
}

traitdf <- data.frame(row.names=names(ranges.crop), Dispersal.syndrome, 
                      Seed.dry.mass, Seed.length, Seed.width, Seed.thickness,
                      Dispersal.unit.type, Dispersal.unit.dry.mass, 
                      Dispersal.unit.length, Dispersal.unit.width, 
                      Dispersal.unit.thickness., check.rows=T, 
                      stringsAsFactors=F)

