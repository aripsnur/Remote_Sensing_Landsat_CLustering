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
path <- 'C:/Users/Admin/Desktop/LAMAR UNIVERSITY/MES CE/Spring24/Remote Sensing for CEE/Assignment 1/pca_score'
setwd(path)

#composit files
comp_files <- list.files(path = "C:/Users/Admin/Desktop/LAMAR UNIVERSITY/MES CE/Spring24/Remote Sensing for CEE/Assignment 1/lsat_composits", full.names = TRUE)

#bands
bands <- c("aerosolx", "bluex", "greenx", "redx", "NIRx", "SWIR2x")

#number of year
years <- as.character(c(1984:1984))


#PCA loop

for (i in seq_along(years)){
  composite <- rast(comp_files[i])
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
  clamped <- c(aerosolxc, bluexc, greenxc, redxc, NIRxc, SWIR2xc)
  #create a dataframe
  bandcdf <- as.data.frame(clamped, xy = TRUE, geom = 'WKT')
  bandcdf <- na.omit(bandcdf)
  colnames(bandcdf) <- c('X', 'Y', bands, 'WKT') 
  
  # Perform PCA with bands only
  pca <- prcomp(bandcdf[,3:8], scale = TRUE)
  
  # Get the principal components(scores)
  pcadata <- pca$x
  write.csv(pcadata[,1:3], paste('pcdata_',years[i],'.csv', sep = ""), row.names = F)
  
  #loop update
  print(paste0("loop ",i,' done', sep = ''))
}


























