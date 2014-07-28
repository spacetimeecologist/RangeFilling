require(ggplot2)
require(maptools)
require(plyr)
require(raster)

# raster <- raster('C:/Users/ben/Documents/Data/Output/Yuccabrevifolia/06-24.1705/maxent.1705.envi')
# points <- read.csv('C:/Users/ben/Documents/Data/Output/Yuccabrevifolia/06-24.1705/presence.csv')
# points <- points[, 2:3]
# poly <- shapefile('C:/Users/ben/Documents/Data/USGS_TreeAtlas/yuccbrev/yuccbrev.shp')
# rp <- 14.504343567

getggmap <- function(raster, points, poly) {
  require(ggplot2)
  require(maptools)
  require(plyr)
  require(raster)
  require(rgeos)
    
  # Import map data
  cntry <- shapefile('C:/Users/ben/Documents/Data/50m_cultural/ne_50m_admin_0_countries.shp')
  coast <- shapefile('C:/Users/ben/Documents/Data/50m_physical/ne_50m_coastline.shp')
  lakes <- shapefile('C:/Users/ben/Documents/Data/50m_physical/ne_50m_lakes.shp')
  river <- shapefile('C:/Users/ben/Documents/Data/50m_physical/ne_50m_rivers_lake_centerlines.shp')
  state <- shapefile('C:/Users/ben/Documents/Data/50m_cultural/ne_50m_admin_1_states_provinces_lakes.shp')
  
  #   raster.nat <- crop(raster('C:/Users/ben/Documents/Data/50m_physical/GRAY_50M_SR_OB/GRAY_50M_SR_OB_crop.asc'), extent(raster))
  #   p.nat <- rasterToPoints(raster.nat)
  #   df.nat <- data.frame(p.nat)
  #   colnames(df.nat) <- c('Longitude', 'Latitude', 'Elevation')
  
  quick.subset <- function(x, longlat){
    
    # longlat should be a vector of four values: c(xmin, xmax, ymin, ymax)
    x@data$id <- rownames(x@data)
    
    x.f = fortify(x, region="id")
    x.join = join(x.f, x@data, by="id")
    
    x.subset <- subset(x.join, x.join$long > longlat[1] & x.join$long < longlat[2] &
                         x.join$lat > longlat[3] & x.join$lat < longlat[4])
    
    x.subset
  }
  
  domain <- c((xmin(raster)*1.1), (xmax(raster)*0.9), (ymin(raster)*0.9), (ymax(raster)*1.1))
  coast.subset <- quick.subset(coast, domain)
  lakes.subset <- quick.subset(lakes, domain)
  river.subset <- quick.subset(river, domain)
  state.subset <- quick.subset(state, domain)
  
  # Remove raster pixels outside of coast. Convert to data.frame.
  cntry.crop <- crop(cntry, extent(raster))
  raster.mask <- mask(raster, cntry.crop)
  p <- rasterToPoints(raster.mask) 
  df <- data.frame(p)
  colnames(df) <- c('long', 'lat', 'predictions')
  
  # Convert points to data.frame
  n <- data.frame(points)
  colnames(n) <- c('long', 'lat')
  
  # Convert poly to data.frame
  poly@data$id  <- rownames(poly@data)
  r <- fortify(poly, region='id')
  
  ggmap <- ggplot() + geom_tile(data=df, aes(x=long, y=lat, fill=factor(predictions))) + 
    scale_fill_manual(values= c( "5"="#69BE13", "4"="#85BE4B", "3"="#A5B398",
                                 "2"="#B3AB98", "1"="#BEB4A2", "0"="#CCBFAD"), 
                      breaks= c("5", "4", "3", "2", "1", "0"),
                      labels= c("Five", "Four", "Three", "Two", "One", "Zero")) +  
    geom_path(data=state.subset, aes(x=long, y=lat, group=group), colour='#878378', size=.2) +
    geom_path(data=river.subset, aes(x=long, y=lat, group=group), colour='#4F5069', size=.3) +
    geom_polygon(data=lakes.subset, aes(x=long, y=lat, group=group), fill='#727497') +
    geom_path(data=coast.subset, aes(x=long, y=lat, group=group), colour='#4F5069', size=.7) +    
    geom_path(data=lakes.subset, aes(x=long, y=lat, group=group), colour='#4F5069', size=.3) +
    geom_path(data=r, aes(x=long, y=lat, group=group), size=.5) +
    geom_point(data=n, aes(x=long, y=lat), color='#EB6200', size=3, shape=16, alpha=.5) +
    coord_quickmap(xlim=c(xmin(raster), xmax(raster)),
                   ylim=c(ymin(raster), ymax(raster))) +
    xlab("Longitude") + ylab("Latitude") +
    guides(fill=guide_legend(title='# of Models\nPredicting Presence')) +
    theme(panel.background = element_rect(fill = '#727497'),
          panel.grid.major = element_line(colour='#BFBAAB', size=.2),
          panel.grid.minor = element_blank(),
          panel.border= element_rect(colour='#403E3B', fill=NA))
  
  
  return(ggmap)
}