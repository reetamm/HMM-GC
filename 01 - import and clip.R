library(raster)
library(ncdf4)
library(rgdal)
library(sp)
library(ggplot2)
library(lubridate)
source('00 - functionDefns.R')

fnames <- list.files(pattern = "3B-DAY.MS.MRG.3IMERG.*", full.names = F, recursive = T, include.dirs = T) #gets file names from working directory; can specify path
precipdata <- data.frame(layer = NA)
iter <- 1
for (fname in fnames){
  print(iter)
  imerg <- nc_open(fname)
  data <- ncvar_get(imerg, 'precipitationCal')
  nc_close(imerg)
  imerg_raster <- raster::raster(x = as.matrix(data), xmn = long[1], xmx = long[NROW(long)], ymn = lat[1], ymx = lat[NROW(lat)] , crs=sp::CRS('+proj=longlat +datum=WGS84'))
  if(iter == 1){
    fname2 <- 'CBAY/Chesapeake_Bay_Watershed.shp'
    myExtent <- readOGR(fname2)
    myExtent <- spTransform(myExtent, CRS(proj4string(imerg_raster)))
    }
  iter <- iter+1
  imerg_crop <- mask(x = imerg_raster, mask = myExtent)
  #plot(imerg_crop)
  imerg_matrix <- rasterToPoints(imerg_crop)
  #dim(imerg_matrix)
  #head(imerg_matrix)
  precipdata <- cbind(precipdata, imerg_matrix[, 3])
}
latlong <- imerg_matrix[, -3]
write.table(latlong, 'latlong', col.names = F, row.names = F)
precipdata <- precipdata[,-1]
precipdata <- as.data.frame(t(precipdata))
date <- seq(ymd('2019-12-30'), ymd('2019-12-31'), by = 'days')
precipdata <- cbind(date,precipdata)

write.table(precipdata, 'datafiles/cbaytotal', row.names = F, col.names = F, sep = ' ')
cbaywet <- extractmonth(precipdata, 7:9)
write.table(cbaywet, 'datafiles/cbaywet', row.names = F, col.names = F, sep = ' ')
