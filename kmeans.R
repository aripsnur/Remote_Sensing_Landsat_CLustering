#Spatial analysis: composite_crop_mask_mosoiac

#Load Libraries

library(terra) #mostly for raster calculuation
library(gdata)
library(corrplot)
library(sf) #mostly for vector calculuation #simple feature
library(pracma)
library(stats)
library(ggplot2)
library(cluster)
library(factoextra)
library(raster)

#set path
path <- 'C:/Users/Admin/Desktop/LAMAR UNIVERSITY/MES CE/Spring24/Remote Sensing for CEE/Assignment 1/kmeans'
setwd(path)

#composit files
comp_files <- list.files(path = "C:/Users/Admin/Desktop/LAMAR UNIVERSITY/MES CE/Spring24/Remote Sensing for CEE/Assignment 1/lsat_composits", full.names = TRUE)
#pca files
pca_files <- list.files(path = "C:/Users/Admin/Desktop/LAMAR UNIVERSITY/MES CE/Spring24/Remote Sensing for CEE/Assignment 1/pca_score", full.names = TRUE)


#bands
bands <- c("aerosolx", "bluex", "greenx", "redx", "NIRx", "SWIR1x")

#number of year
years <- as.character(c(1984:201))


#PCA loop

for (i in seq_along(years)){
  composite <- rast(comp_files[i])
  pcadata <- read.csv(pca_files[i])
  pcadata <- as.matrix(pcadata)
  #write the bands as a data frame for additional analysis
  banddf <- as.data.frame(composite, xy = TRUE, geom = 'WKT')
  banddf <- na.omit(banddf)
  colnames(banddf) <- c('X', 'Y', bands)
  #Clamp values for contrast enhancement #histogram equalization
  LL <- 0.05
  UL <- 0.95
  bandnamec <- c()
  for (k in seq(1, length(bands), 1)){
    bname <- paste(bands[k], 'c', sep = "")
    xmin <- quantile(banddf[, (k+2)], LL, na.rm = T)
    xmax <- quantile(banddf[, (k+2)], UL, na.rm = T)
    x <- clamp(composite[[k]], lower = xmin, upper = xmax, values = T)
    mv('x', bname)
    bandnamec[k] <- bname
  }
  clamped <- c(aerosolxc, bluexc, greenxc, redxc, NIRxc, SWIR1xc)
  #create a dataframe
  bandcdf <- as.data.frame(clamped, xy = TRUE, geom = 'WKT')
  bandcdf <- na.omit(bandcdf)
  colnames(bandcdf) <- c('X', 'Y', bands, 'WKT') 
  
  #Reproject data to lat-lon
  crsaea <- crs(clamped, proj = T) #Get crs
  crs84 <- 4326 #EPSG code #Define the new crs
  
  banddf.SP <- st_as_sf(bandcdf, coords = c('X', 'Y'), crs = crsaea)
  banddf.SP$XAEA = st_coordinates(banddf.SP)[,1]
  banddf.SP$YAEA = st_coordinates(banddf.SP)[,2]
  
  
  bandcdf.SP <- st_transform(x = banddf.SP, crs = crs84) #reproject
  bandcdf.SP$Lon = st_coordinates(bandcdf.SP)[,1]
  bandcdf.SP$Lat = st_coordinates(bandcdf.SP)[,2]
  
  #perform clustering
  
  bpcaclust <- kmeans(pcadata[,1:3], centers=4, nstart = 10)
  clusts <- bpcaclust$cluster
  pcadat.clust <- data.frame(bandcdf.SP, clusts)
  pcadat.clust <- subset(pcadat.clust, select= c('Lon', 'Lat', 'clusts'))
  write.csv(pcadat.clust, paste('kmeans_', years[i],'.csv', sep = ""), row.names = F)
  
  
  #loop update
  print(paste0("loop ",i,' done', sep = ''))
}


























