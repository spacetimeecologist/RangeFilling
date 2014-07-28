library(rgdal)
library(raster)

ext <- extent(-130, -95, 20, 50)

getEnvi <- function(ext) {
  require(rgdal)
  require(raster)
  # extent of analysis.
  ext <- ext
  
  # base climate layers to calculate evironmental variables from. bio is a set of 19
  # bioclim vars, tmean is mean temperature (hi and lo) for 12 months, prec is total
  # precipitation for 12 months, and pet is potential evapotranspiration for every month
  bio <- getData('worldclim', var='bio', res=5)
  tmean <- getData('worldclim', var='tmean', res=5)
  prec <- getData('worldclim', var='prec', res=5)
  pet <- stack(
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet1.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet2.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet3.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet4.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet5.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet6.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet7.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet8.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet9.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet10.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet11.crop.tif'), ext),
    crop(raster('C:/Users/ben/Documents/Data/PET_Monthy/pet12.crop.tif'), ext))
  
  # crops environmental layers to extent for efficiency. pet is already cropped
  bio <- crop(bio, ext, unstack=TRUE)
  tmean <- crop(tmean, ext, unstack=TRUE)
  prec <- crop(prec, ext, unstack=TRUE)
  
  # tmin is mean low temperatures for coldest quarter
  tmin <- (raster(bio, layer=11))/10
  
  # psea is precipitation seasonality
  psea <- ((raster(prec, layer=6) + raster(prec, layer=7) + raster(prec, layer=8))
          - (raster(prec, layer=12) + raster(prec, layer=1) + raster(prec, layer=2)))
  
  # gets and returns a raster of growing degree days based on t.mean and a base temp (tbase)
  tbase <- 5
  
  getGDD <- function(tmean, tbase){
    tbase <- tbase*10
    days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    tmean.base <- (tmean-tbase)
    values.to.zero <- c((min(minValue(tmean.base))-1), tbase, 0)
    tmean.base <- reclassify(tmean.base, values.to.zero)
    gdd <- overlay((tmean.base*days)/10, fun=sum, unstack=TRUE)
    
    return(gdd)
  }
  
  gdd <- getGDD(tmean, tbase)
  
  # gets and returns waterbalance (precip-potential evap)
  getWaterBal <- function(prec, pet){
    agg.multi <- (res(subset(prec, 1)))/res(pet)
    pet.agg <- aggregate(pet, agg.multi)
    waterbal <- (prec-pet.agg)
    waterbal.total <- overlay(waterbal, fun=sum, unstack=TRUE)
    
    return(waterbal.total)
  }
  
  waterbal <- getWaterBal(prec, pet)
  
  # envi is rasterstack of environmental layers used for predictions
  envi <- stack(gdd, tmin, psea, waterbal)
  
  return(envi)
}

envi <- getEnvi(ext)
