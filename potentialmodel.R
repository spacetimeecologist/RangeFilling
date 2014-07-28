getRP <- function(genus, species, envi) {
  require(rgdal)
  require(raster)
  require(dismo)
  require(XML)
  require(rJava)
  require(rgeos)
  # Loads occurence data from GBIF and turns it into spatial points. 
  set.seed(1991)
  pres <- gbif(genus=genus, species=paste(species, '*', sep=""), ext=extent(envi))
  pres <- pres[, 7:8]
  pres <- na.omit(pres)
  pres <- pres[, c(2, 1)]
  pres <- unique(pres)
  pres <- SpatialPoints(pres, bbox=as.matrix(ext))
  print('Species occurence data retrieved')
  
  # makes random background points within extent. seed is 
  set.seed(1991)
  back <- randomPoints(envi, n=75, pres, ext=ext, excludep=TRUE)
  back <- SpatialPoints(back)
  
  
  # divides presence into 5 groups. 4 are used for model training (building) and
  # 1 is used for model testing (evaluating)
  k <- 5
  set.seed(1991)
  group <- kfold(pres, k)
  train <- pres[group != 1]
  test <- pres[group == 1]
  print("Occurence data split into training and testing sets")
  
  # Functions used for model building. 'xx' is the type of model used. 'xx.pred'
  # is the suitability predictions predicted over the environment. 'xx.e' is the
  # evaluation object for the model. 'xx.tr' is the threshold value that separates
  # presence from absence. 'xx.prab' is the cells with greater suitablity than the
  # threshold value, respresenting the potential range. 'xx.area' is the amount
  # of potential area that 'xx.prab' defines.
  
  #   bc.e <- list()
  #   bc.tr <- list()
  #   bc.pred <- list()
  #   bc.prab <- list()
  #   i <- 1  
  #   for (i in 1:k) {
  #     train <- pres[group != i,]
  #     test <- pres[group == i,]
  #     bc <- bioclim(envi, train)
  #     bc.e[[i]] <- evaluate(p=test, a=back, bc, envi)
  #     bc.tr[[i]] <- threshold(bc.e[[i]], stat='spec_sens')
  #     bc.pred[[i]] <- predict(bc, envi)
  #     bc.prab[[i]] <- (bc.pred[[i]] > bc.tr[[i]])
  #   }
  #   bc.p <- sum(bc.prab[[1]], bc.prab[[2]], bc.prab[[3]], bc.prab[[4]], bc.prab[[5]])
  #   names(bc.prab) <- 'Bioclim'
  #   bc.area <- (bc.p >= 4)
  #   bc.area <- cellStats(bc.area, 'sum', na.rm=TRUE)*prod(res(bc.p))  
  #   print("Bioclim Model completed")
  #   
  #   dm.e <- list()
  #   dm.tr <- list()
  #   dm.pred <- list()
  #   dm.prab <- list()
  #   i <- 1  
  #   for (i in 1:k) {
  #     train <- pres[group != i,]
  #     test <- pres[group == i,]
  #     dm <- domain(envi, train)
  #     dm.e[[i]] <- evaluate(p=test, a=back, dm, envi)
  #     dm.tr[[i]] <- threshold(dm.e[[i]], stat='spec_sens')
  #     dm.pred[[i]] <- predict(dm, envi)
  #     dm.prab[[i]] <- (dm.pred[[i]] > dm.tr[[i]])
  #   }
  #   dm.p <- sum(dm.prab[[1]], dm.prab[[2]], dm.prab[[3]], dm.prab[[4]], dm.prab[[5]])
  #   names(dm.prab) <- 'Domain' 
  #   dm.area <- (dm.p >= 4)
  #   dm.area <- cellStats(dm.area, 'sum', na.rm=TRUE)*prod(res(dm.p)) 
  #   print("Domain Model completed")
  #   
  #   mh.e <- list()
  #   mh.tr <- list()
  #   mh.pred <- list()
  #   mh.prab <- list()
  #   i <- 1  
  #   for (i in 1:k) {
  #     train <- pres[group != i,]
  #     test <- pres[group == i,]
  #     mh <- mahal(envi, train)
  #     mh.e[[i]] <- evaluate(p=test, a=back, mh, envi)
  #     mh.tr[[i]] <- threshold(mh.e[[i]], stat='spec_sens')
  #     mh.pred[[i]] <- predict(mh, envi)
  #     mh.prab[[i]] <- (mh.pred[[i]] > mh.tr[[i]])
  #   }
  #   mh.p <- sum(mh.prab[[1]], mh.prab[[2]], mh.prab[[3]], mh.prab[[4]], mh.prab[[5]])
  #   names(mh.prab) <- 'Mahalanobis'
  #   mh.area <- (mh.p >=4)
  #   mh.area <- cellStats(mh.area, 'sum', na.rm=TRUE)*prod(res(mh.p))
  #   print("Mahalanobis Model completed")
  
  mx.e <- list()
  mx.tr <- list()
  mx.pred <- list()
  mx.prab <- list()
  i <- 1  
  for (i in 1:k) {
    train <- pres[group != i,]
    test <- pres[group == i,]
    mx <- maxent(envi, train, a=back)
    mx.e[[i]] <- evaluate(p=test, a=back, mx, envi)
    mx.tr[[i]] <- threshold(mx.e[[i]], stat='spec_sens')
    mx.pred[[i]] <- predict(mx, envi)
    mx.prab[[i]] <- (mx.pred[[i]] > mx.tr[[i]])
  }
  mx.p <- sum(mx.prab[[1]], mx.prab[[2]], mx.prab[[3]], mx.prab[[4]], mx.prab[[5]])
  names(mx.prab) <- 'Maxent'
  mx.area <- (mx.p >= 4)
  mx.area <- cellStats(mx.area, 'sum', na.rm=TRUE)*prod(res(mx.p))
  
  print("Maxent Model completed")
  print(mx.area)
  
  #  all.p <- c(bc.p, dm.p, mh.p, mx.p)
  
  all.p <- mx.p
  
  # code for calculating range filling
  #realized <- circles(pres, d=25000, lonlat=TRUE)
  #realized.p <- as(realized, 'SpatialPolygons')
  realized <- crop(shapefile(paste('C:/Users/ben/Documents/Data/USGS_TreeAtlas/',
                                   tolower(substr(genus, start=0, stop=4)),
                                   substr(species, start=0, stop=4), '/',
                                   tolower(substr(genus, start=0, stop=4)),
                                   substr(species, start=0, stop=4),'.shp', sep='')), ext)
  print("Little range found and retrieved")
  
  # creates raster and images of model results and puts them into a folder
  # based on the species' name and the time (stamp) which the model was ran.
#   stamp <- Sys.time()
#   print(stamp)
#   dir.create(path=paste('C:/Users/ben/Documents/Data/Output/', genus, species, 
#                         sep=""))
#   
#   dir.create(path=paste('C:/Users/ben/Documents/Data/Output/', genus, species, 
#                         '/', substr(Sys.time(), start=6, stop=10),'.', 
#                         substr(stamp, start=12, stop=13), 
#                         substr(stamp, start=15, stop=16),
#                         sep=""))
#   print('Directory Made')
#   
#   write.csv(pres, file=paste('C:/Users/ben/Documents/Data/Output/', genus, 
#                              species, '/', substr(Sys.time(), start=6, stop=10),
#                              '.', substr(stamp, start=12, stop=13), 
#                              substr(stamp, start=15, stop=16),
#                              '/presence.csv', sep=""))
  #   writeRaster(bc.p, 
  #               filename=paste('C:/Users/ben/Documents/Data/Output/', genus, 
  #                              species, '/', substr(stamp, start=6, stop=10),'.', 
  #                              substr(stamp, start=12, stop=13), 
  #                              substr(stamp, start=15, stop=16),'/bioclim.',
  #                              substr(stamp, start=12, stop=13), 
  #                              substr(stamp, start=15, stop=16),'.envi',
  #                              sep=""))
  #   print("Bioclim Raster made")
  #   writeRaster(dm.p, 
  #               filename=paste('C:/Users/ben/Documents/Data/Output/', genus,
  #                              species, '/', substr(stamp, start=6, stop=10),'.', 
  #                              substr(stamp, start=12, stop=13), 
  #                              substr(stamp, start=15, stop=16),'/domain.',
  #                              substr(stamp, start=12, stop=13), 
  #                              substr(stamp, start=15, stop=16),'.envi',
  #                              sep=""))
  #   print("Domain Raster made")
  #   writeRaster(mh.p, 
  #               filename=paste('C:/Users/ben/Documents/Data/Output/', genus, 
  #                              species, '/', substr(stamp, start=6, stop=10),'.', 
  #                              substr(stamp, start=12, stop=13), 
  #                              substr(stamp, start=15, stop=16),'/mahal.',
  #                              substr(stamp, start=12, stop=13), 
  #                              substr(stamp, start=15, stop=16),'.envi',
  #                              sep=""))
  #   print("Mahal Raster made")
#   writeRaster(mx.p, 
#               filename=paste('C:/Users/ben/Documents/Data/Output/', genus, 
#                              species, '/', substr(stamp, start=6, stop=10),'.', 
#                              substr(stamp, start=12, stop=13), 
#                              substr(stamp, start=15, stop=16),'/maxent.',
#                              substr(stamp, start=12, stop=13), 
#                              substr(stamp, start=15, stop=16),'.envi',
#                              sep=""))
#   print("Maxent Raster made")
#   
  #   png(width=1280, height=1280, 
  #       filename=paste('C:/Users/ben/Documents/Data/Output/', genus, species, 
  #                      '/', substr(stamp, start=6, stop=10),'.', 
  #                      substr(stamp, start=12, stop=13), 
  #                      substr(stamp, start=15, stop=16),'/bioclim.',
  #                      substr(stamp, start=12, stop=13), 
  #                      substr(stamp, start=15, stop=16),'.png',
  #                      sep=""),
  #       bg='white')
  #   bc.gg <- getggmap(bc.p, pres, realized)
  #   bc.gg <- bc.gg + ggtitle(paste('Bioclim for ', genus, ' ', species, 
  #                                  '\n Range Filled: ', 
  #                                  round(sum(area(realized)/bc.area)*100, 4),
  #                                  '%', sep="")) +
  #     theme(plot.title = element_text(size=24, face="bold"),
  #           axis.title.x = element_text(colour='#403E3B',size=20),
  #           axis.title.y = element_text(colour='#403E3B',size=20, angle=90),
  #           axis.text.x = element_text(colour='#403E3B', size=18),
  #           axis.text.y = element_text(colour='#403E3B', size=18))
  #   print(bc.gg)
  #   dev.off()
  #   print("Bioclim plot made")
  
  #   png(width=1280, height=1280, 
  #       filename=paste('C:/Users/ben/Documents/Data/Output/', genus, species, 
  #                      '/', substr(stamp, start=6, stop=10),'.', 
  #                      substr(stamp, start=12, stop=13), 
  #                      substr(stamp, start=15, stop=16),'/domain.',
  #                      substr(stamp, start=12, stop=13), 
  #                      substr(stamp, start=15, stop=16),'.png',
  #                      sep=""),
  #       bg='white')
  #   dm.gg <- getggmap(dm.p, pres, realized)
  #   dm.gg <- dm.gg + ggtitle(paste('Domain for ', genus, ' ', species,
  #                                  '\n Range Filled: ', 
  #                                  round(sum(area(realized)/dm.area)*100, 4),
  #                                  '%', sep="")) +
  #     theme(plot.title = element_text(size=24, face="bold"),
  #           axis.title.x = element_text(colour='#403E3B',size=20),
  #           axis.title.y = element_text(colour='#403E3B',size=20, angle=90),
  #           axis.text.x = element_text(colour='#403E3B', size=18),
  #           axis.text.y = element_text(colour='#403E3B', size=18))
  #   print(dm.gg)
  #   dev.off()
  #   print("Domain plot made")
  #   
  #   png(width=1280, height=1280, 
  #       filename=paste('C:/Users/ben/Documents/Data/Output/', genus, species, 
  #                      '/', substr(stamp, start=6, stop=10),'.', 
  #                      substr(stamp, start=12, stop=13), 
  #                      substr(stamp, start=15, stop=16),'/mahal.',
  #                      substr(stamp, start=12, stop=13), 
  #                      substr(stamp, start=15, stop=16),'.png',
  #                      sep=""),
  #       bg='white')
  #   mh.gg <- getggmap(mh.p, pres, realized)
  #   mh.gg <- mh.gg + ggtitle(paste('Mahalanobis for ', genus, ' ', species, 
  #                                  '\n Range Filled: ', 
  #                                  round(sum(area(realized)/mh.area)*100, 4),
  #                                  '%', sep="")) +
  #     theme(plot.title = element_text(size=24, face="bold"),
  #           axis.title.x = element_text(colour='#403E3B',size=20),
  #           axis.title.y = element_text(colour='#403E3B',size=20, angle=90),
  #           axis.text.x = element_text(colour='#403E3B', size=18),
  #           axis.text.y = element_text(colour='#403E3B', size=18))
  #   print(mh.gg)
  #   dev.off()
  #   print("Mahal plot made")
  #   
#   png(width=1280, height=1280, 
#       filename=paste('C:/Users/ben/Documents/Data/Output/', genus, species, 
#                      '/', substr(stamp, start=6, stop=10),'.', 
#                      substr(stamp, start=12, stop=13), 
#                      substr(stamp, start=15, stop=16),'/maxent.',
#                      substr(stamp, start=12, stop=13), 
#                      substr(stamp, start=15, stop=16),'.png',
#                      sep=""),
#       bg='white')
  mx.gg <- getggmap(mx.p, pres, realized)
  mx.gg <- mx.gg + ggtitle(paste('Maxent for ', genus, ' ', species, 
                                 '\n Range Filled: ', 
                                 round(sum(area(realized)/mx.area)*100, 4),
                                 '%', sep="")) +
    theme(plot.title = element_text(size=24, face="bold"),
          axis.title.x = element_text(colour='#403E3B',size=20),
          axis.title.y = element_text(colour='#403E3B',size=20, angle=90),
          axis.text.x = element_text(colour='#403E3B', size=18),
          axis.text.y = element_text(colour='#403E3B', size=18))
#   print(mx.gg)
#   dev.off()
plot(mx.gg)  
print("Maxent plot made")
  
  return(all.p)
}


