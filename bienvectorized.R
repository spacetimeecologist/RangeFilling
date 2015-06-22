require(rgdal)
require(maptools)
require(rgeos)
require(plyr)
require(dismo)
require(foreach)
require(doParallel)


bien.all <- read.csv('C:\\Users\\ben\\Documents\\Research\\RangeFilling\\Data\\BIEN\\All data (v2)\\All data (v2).csv', 
                     stringsAsFactors=F)

bien.gymno <- read.csv("C:\\Users\\ben\\Documents\\Research\\RangeFilling\\Data\\BIEN\\Gymnocladus dioica.csv",
                       stringsAsFactors=F)
colnames(bien.gymno) <- colnames(bien.all)
bien.gymno$Species <- gsub('Gymnocladus dioica', 'Gymnocladus dioicus', 
                           bien.gymno$Species)


bien.all <- rbind(bien.all, bien.gymno) 


bien.sp <- unique(bien.all$Species)
bien.sp <- gsub(' ', '', bien.sp)
bien.sp <- bien.sp[na.exclude(match(bien.sp, row.names(df)))]
bien.notsp <- bien.sp[which(is.na(match(bien.sp, row.names(df)))==T)]

bien.clean <- bien.all[-which((gsub(' ', '', bien.all$Species[]) %in% bien.notsp)==T),]



workers <- makeCluster(4)
registerDoParallel(workers)
clusterEvalQ(workers, c(library(raster), library(sp)))
clusterExport(workers, varlist=c('bien.clean', 'bbox'), envir=.GlobalEnv)

# This creates a list of spatial points (bien.pts) for every species within
# bien.clean. Duplicated points are removed for every species, and species
# with less than 100 points are not included.
bien.pts <- llply(.data=c(1:length(unique(bien.clean$Species))),
                  .fun=function(X){X <- unique(bien.clean$Species)[X]
                                   pts <- bien.clean[which(bien.clean$Species==X), ]
                                   pts <- pts[which(duplicated(pts[,3:4])==F),]
                                   pts <- SpatialPoints(pts[c(4,3)],
                                                        proj4string=CRS('+proj=aea 
                                                                        +lat_1=29.5 
                                                                        +lat_2=45.5 
                                                                        +lat_0=37.5 
                                                                        +lon_0=-96 
                                                                        +x_0=0 
                                                                        +y_0=0 
                                                                        +datum=NAD83 
                                                                        +units=m 
                                                                        +no_defs'))
                                   pts <- crop(pts, extent(bbox))
                                   if(length(pts)>=100){
                                     return(pts)}
                                   else{
                                     return(NULL)
                                   }},
                  .parallel=T
                  #.progress='text'
)

names(bien.pts) <- unique(bien.clean$Species)
stopCluster(workers)

bien.null <- list()
g <- 1
for(g in g:length(bien.pts)){
  bien.null[g] <- is.null(bien.pts[[g]])
}

names(bien.null) <- names(bien.pts)

bien.pts <- bien.pts[which(bien.null==F)]
names(bien.pts) <- names(bien.null[which(bien.null==F)])

bien.centroids <- llply(.data=c(1:length(bien.pts)),
                        .fun=function(X){X <- names(bien.pts[X])
                                         sp.centroid <- gCentroid(bien.pts[[X]])
                        return(sp.centroid)})

bien.coords <- ldply(bien.centroids, coordinates)
bien.lon <- data.frame(bien.coords$x)
bien.lat <- data.frame(bien.coords$y)

workers <- makeCluster(4)
registerDoParallel(workers)
clusterEvalQ(workers, c(library(raster), library(rgeos)))
clusterExport(workers, varlist=c('bien.pts'), envir=.GlobalEnv)

# this creates a buffer around every point and then turns all buffers into a 
# single polygon. This takes ~20 minutes when running parallel.
bien.buff <- llply(.data=c(1:length(bien.pts)),
                   .fun=function(X){X <- names(bien.pts)[X]
                                    buffer <- aggregate(gBuffer(bien.pts[[X]],
                                                                byid=T, width=1))
                                    return(buffer)},
                   .parallel=T)
names(bien.buff) <- names(bien.pts)

stopCluster(workers)

bien.not.ranges <- llply(.data=c(1:length(bien.buff)), 
                         .fun=function(X){X=names(bien.buff)[X]
                                          land-bien.buff[[X]]},
                         .progress='text')
names(bien.not.ranges) <- names(bien.buff)

bien.abs <- llply(.data=c(1:length(bien.not.ranges)),
                  .fun=function(X){X=names(bien.not.ranges)[X]
                                   bien.abs <- spsample(bien.not.ranges[[X]], 
                                                        n=100 + 
                                                          (100000/sqrt(gArea(bien.buff[[X]]))),
                                                        type='regular')
                                   return(bien.abs)},
                  .progress='text')
names(bien.abs) <- names(bien.not.ranges)
set.seed(1991)

group <- llply(bien.pts, kfold, k=5)
names(group) <- names(bien.pts)


k <- 1  
for (k in k:5) {
  print(k)
  train <- llply(.data=c(1:length(bien.pts)), 
                 .fun=function(X){X=names(bien.pts[X])
                                  bien.pts[[X]][which(group[[X]]!=k)]})
  test <- llply(.data=c(1:length(bien.pts)),
                .fun=function(X){X=names(bien.pts[X]) 
                                 bien.pts[[X]][which(group[[X]]==k)]})
  names(train) <- names(bien.pts)
  names(test) <- names(bien.pts)
  assign(paste('bc.p.', k, sep=''), 
         value=climateEnvelope2(test, train, p.abs=bien.abs, pres=bien.pts, envi))
}

bc.sum.bien <- llply(.data=c(1:length(bc.p.1)), .fun=function(X){X=names(bc.p.1[X])
                                                                 sum(bc.p.1[[X]], 
                                                                     bc.p.2[[X]], 
                                                                     bc.p.3[[X]], 
                                                                     bc.p.4[[X]], 
                                                                     bc.p.5[[X]])},
                     .progress='text')
names(bc.sum.bien) <- names(bc.p.1)

bc.sum.bien <- llply(.data=c(1:length(bc.sum.bien)), 
                     .fun=function(X){X=names(bc.sum.bien[X])
                                      bc.sum.bien[[X]] > 4},
                     .progress='text')
names(bc.sum.bien) <- names(bc.p.1)

bc.area.bien <- bc.sum.bien
names(bc.area.bien) <- names(bc.sum.bien)

bc.area.bien <- llply(.data=c(1:length(bc.area.bien)), 
                   .fun=function(X){X=names(bc.area.bien[X])
                                    bc.area.bien[[X]][bc.area.bien[[X]]==FALSE] <- NA
                                    return(bc.area.bien[[X]])},
                   .progress='text')
names(bc.area.bien) <- names(bc.sum.bien)

bc.km.bien <- ldply(.data=c(1:length(bc.area.bien)), 
               .fun=function(X){X=names(bc.area.bien[X])
                                cellarea <- area(bc.area.bien[[X]], na.rm=TRUE, weights=F)
                                cellarea <- sum(cellarea[!is.na(cellarea)])
                                return(cellarea)},
               .progress='text')
row.names(bc.km.bien) <- names(bc.area.bien)

bc.centroids.bien <- llply(.data=c(1:length(bc.area.bien)),
                     .fun=function(X){X=names(bc.area.bien[X])
                                      cells <- xyFromCell(bc.area.bien[[X]],
                                                          cell=which(values(bc.area.bien[[X]])==T),
                                                          spatial=T)                
                                      centroid <- gCentroid(cells)
                                      return(centroid)},
                     .progress='text')
names(bc.centroids.bien) <- names(bc.area.bien)

bc.centroid.coords <- ldply(bc.centroids.bien, coordinates)
bc.lon.bien <- data.frame(bc.centroid.coords$x, row.names=names(bc.centroids.bien))
bc.lat.bien <- data.frame(bc.centroid.coords$y, row.names=names(bc.centroids.bien))
remove(bc.centroid.coords)

bc.shape.bien <- ldply(.data=c(1:length(bc.area.bien)),
                 .fun=function(X){X=names(bc.area.bien[X])
                                  cells <- xyFromCell(bc.area.bien[[X]],
                                                      cell=which(values(bc.area.bien[[X]])==T),
                                                      spatial=T)                
                                  xlength <- xmax(extent(cells))-xmin(extent(cells))
                                  ylength <- ymax(extent(cells))-ymin(extent(cells))
                                  shape <- ylength/xlength
                                  return(shape)},
                 .progress='text')
row.names(bc.shape.bien) <- names(bc.area.bien)

bien.shape <- ldply(.data=c(1:length(bien.pts)),
                    .fun=function(X){X=names(bien.pts[X])              
                                     xlength <- xmax(extent(bien.pts[[X]]))
                                     -xmin(extent(bien.pts[[X]]))
                                     ylength <- ymax(extent(bien.pts[[X]]))
                                     -ymin(extent(bien.pts[[X]]))
                                     shape <- ylength/xlength
                                     return(shape)},
                    .progress='text')
row.names(bien.shape) <- names(bien.pts)


bien.df <- data.frame(row.names=names(bien.pts), 
                      R.km2[row.names(R.km2)%in%gsub(' ', '', names(bien.pts))==T,], 
                      bc.km.bien, 
                      (R.km2[row.names(R.km2)%in%gsub(' ', '', names(bien.pts))==T,]/bc.km.bien)*100,
                      bien.lat, bien.lon, bc.lat.bien, bc.lon.bien, bien.shape, bc.shape.bien,
                      bien.lon<(-100), (bien.lon>(-84)&bien.lat<(28)),
                      check.rows=T, stringsAsFactors=F)
colnames(bien.df) <- c('R.km2', 'bc.km.bien', 'RP', 'bien.lat', 'bien.lon', 'P.lat', 'P.lon',
                       'bien.shape', 'P.shape', 'West', 'Everglades')
