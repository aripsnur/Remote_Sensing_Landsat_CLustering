#Spatial analysis of landsat ARD data
#Exploratory data analysis

#step1: Load Libraries

library(terra) #mostly for raster calculuation
library(gdata)
library(corrplot)
library(sf) #mostly for vector calculuation #simple feature
library(pracma)
library(stats)
library(ggplot2)
library(cluster)
library(factoextra)

Lsat_files <- list.files(path = "F:/Semester 4 Spring 2024/Remote_sensing/Landsat_2014_2023", full.names = TRUE)

bound <- vect("F:/Semester 4 Spring 2024/Remote_sensing/Data/Chambers.shp")

bound<-  project(bound, "epsg:4326", partial = F)
a<-crs(bound)

path <- 'F:/Semester 4 Spring 2024/Remote_sensing/Landsat_comb'
setwd(path)
bands <- c("aerosolx", "bluex", "greenx", "redx", "NIRx", "SWIR2x")
year <- as.character(c(2014:2014))

#Read the Chambers vector file
Cha <- vect("F:/Semester 4 Spring 2024/Remote_sensing/Data/chambers.gpkg") #terra and sf library
plot(Cha)


yearwise_file <- list.files(path = Lsat_files[1], full.names = TRUE)

tif_files1 <- list.files(path = yearwise_file[1], pattern = "._B..TIF", full.names = TRUE)

for (i in 1:length(bands)) {
  x <- rast(tif_files1[i])
  mv('x', bands[i])
}

comp1 <- c(aerosolx, bluex, greenx, redx, NIRx, SWIR2x)
plotRGB(comp1, r= 5, g=4, b=2, stretch = 'hist')


tif_files2 <- list.files(path = yearwise_file[2], pattern = "._B..TIF", full.names = TRUE)

for (i in 1:length(bands)) {
  x <- rast(tif_files2[i])
  mv('x', bands[i])
}

comp2 <- c(aerosolx, bluex, greenx, redx, NIRx, SWIR2x)
plotRGB(comp2, r= 5, g=4, b=2, stretch = 'hist')


mosaic_raster <- mosaic(comp1, comp2, fun ="mean")
plotRGB(mosaic_raster, r= 5, g=4, b=2, stretch = 'hist')


masked <- mask(mosaic_raster, Cha)
croped <- crop(masked,Cha)
plotRGB(masked, r= 5, g=4, b=2, stretch = 'hist')
plotRGB(croped, r= 5, g=4, b=2, stretch = 'hist')

croped <- list.files(pattern = "\\.TIF$", full.names = TRUE)
croped <-rast(croped)

writeRaster(croped, paste(year[1],'.TIF', sep = ""),
            filetype = "GTiff", overwrite = TRUE)

#step 7: Create Boxplot

bandnames <- c("aerosolx", "bluex", "greenx", "redx", "NIRx", "SWIR2x")
f <- boxplot(croped, axes = FALSE, outline = F, ylab = "Value",  notch = F)
axis(1, at = 1:length(bandnames), labels = bandnames, las=1)
axis(2)
title('DN for Various bands in Chambers')
grid()
box()


#write the bands as a dataframe for aditional analysis
banddf <- as.data.frame(croped, xy = TRUE, geom = 'WKT')
head(banddf)
banddf <- na.omit(banddf)
colnames(banddf) <- c('X', 'Y', bandnames)
head(banddf)
summary(banddf)

#cor
corbands <- cor(banddf[,3:8], method = 'pearson') #spearman
corrplot(corbands, method = 'number', type = 'lower', diag = F)

#histogram
hist(croped)

#Clamp values for contrast enhancement #histogram equalization
hist(banddf[, 3])
box()
LL <- 0.05
UL <- 0.95
bandnamec <- c()
for (i in seq(1, length(bandnames), 1)){
  bname <- paste(bands[i], 'c', sep = "")
  xmin <- quantile(banddf[, (i+2)], LL, na.rm = T)
  xmax <- quantile(banddf[, (i+2)], UL, na.rm = T)
  x <- clamp(croped[[i]], lower = xmin, upper = xmax, values = T)
  mv('x', bname)
  bandnamec[i] <- bname
}

clamped <- c(aerosolxc, bluexc, greenxc, redxc, NIRxc, SWIR2xc)
plotRGB(clamped, r= 4, g=3, b=2, stretch = 'lin')

#histogram
hist(clamped)

#create a dataframe
bandcdf <- as.data.frame(clamped, xy = TRUE, geom = 'WKT')
summary(bandcdf)
bandcdf <- na.omit(bandcdf)
colnames(bandcdf) <- c('X', 'Y', bandnames, 'WKT') 


#compute correlation between bands
corbands <- cor(bandcdf[,3:8], method = 'spearman')
corrplot(corbands, method = 'number', type = 'lower', diag = FALSE)

#Reproject data to lat-lon
crsaea <- crs(clamped, proj = T) #Get crs
crs84 <- 4326 #EPSG code #Define the new crs

banddf.SP <- st_as_sf(bandcdf, coords = c('X', 'Y'), crs = crsaea)
banddf.SP$XAEA = st_coordinates(banddf.SP)[,1]
banddf.SP$YAEA = st_coordinates(banddf.SP)[,2]
banddf.SP

bandcdf.SP <- st_transform(x = banddf.SP, crs = crs84) #reproject
bandcdf.SP$Lon = st_coordinates(bandcdf.SP)[,1]
bandcdf.SP$Lat = st_coordinates(bandcdf.SP)[,2]
bandcdf.SP

#create a csv file
bandcdf.SP <- subset(bandcdf.SP, select = -c(WKT, geometry))
write.csv(bandcdf.SP, 'bandcdf.csv', row.names = F)

# Perform PCA with bands only
pca <- prcomp(bandcdf[,3:8], scale = TRUE)
summary(pca)

# Get the first two principal components
pcadata <- pca$x
#write.csv(pcadata[,1:3], 'pcdata.csv', row.names = F)

# Look at the first few rows of the rotation matrix
head(pca$rotation)
pc1 <- pca$x[ ,2]
pc2 <- pca$x[ ,3]
pc3 <- pca$x[ ,4]

# Visualize the data points in the first two principal components
ggplot(bandcdf[, 3:8], aes(x = pc1, y = pc2, color = "darkblue")) +
  geom_point() +
  labs(title = " first two principal components",
       x = "PC1", y = "PC2")

#Elbow method
# Define the k range
k.range <- 1:10

# Loop through k values and calculate WCSS
wss <- vector(length = length(k.range))
for (i in 1:length(k.range)) {
  k <- k.range[i]
  model <- kmeans(pcadata[,1:3], centers = k, nstart = 10)
  wss[i] <- sum(model$withinss)
  print(paste0("loop ",i,', done!!'))
}

# Plot the Elbow Graph
plot(k.range, wss, type = "b", main = "Elbow Method", xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Squares (WCSS)")

x#perform clustering
bpcaclust <- kmeans(pcadata[,1:3], centers=5, nstart = 10)
summary(bpcaclust)
clusts <- bpcaclust$cluster
pcadat.clust <- data.frame(bandcdf.SP, clusts)
write.csv(pcadat.clust, 'bandcdfclus.csv', row.names = F)


fviz_cluster(bpcaclust, data = pcadata[,1:3],
             geom = "point",
             ellipse.type='convex')

clusts <- bpcaclust$cluster
pcadat.clust <- data.frame(bandcdf.SP, clusts)
pcadat.clust.sp <- st_as_sf(pcadat.clust, coords = c('Lon', 'Lat'), crs = crs84)
pcadat.clust.sp$Lon = st_coordinates(pcadat.clust.sp)[,1]
pcadat.clust.sp$Lat = st_coordinates(pcadat.clust.sp)[,2]
# Plot the points colored by cluster number

# Plot the points with colors based on cluster
cluster_colors <- c("red", "white", "darkgreen", "orange", "blue")
plot(pcadat.clust.sp$Lon, pcadat.clust.sp$Lat, 
     col = cluster_colors[as.factor(pcadat.clust$clusts)],
     xlab = "Longitude", ylab = "Latitude")
plot(bound, add = TRUE, border = "black", lwd = 5)
legend("topright", legend = levels(as.factor(pcadat.clust$clusts)), 
       col = cluster_colors,pch=1)

# Perform PCA with bands only
pcaX <- princomp(bandcdf[,3:8])
summary(pcax)

#hc <- hclust(dist(pcadata[,1:3]), method = "complete")

# Define chunk size (adjust based on your data and memory constraints)
chunk_size <- 4  # Adjust as needed

# Function to perform HCA on a chunk
cluster_chunk <- function(data_chunk) {
  # Calculate pairwise distances within the chunk
  distances <- dist(data_chunk)
  
  # Perform HCA using complete linkage
  clusters <- hclust(distances, method = "complete")
  
  # Return the cluster labels for the chunk
  return(clusters)
}

# Split the principal component data into chunks
data_chunks <- split(pcadata[, 1:3], seq(1, nrow(pcadata), chunk_size))

# Apply the clustering function to each chunk and store 
all_cluster <- lapply(data_chunks, FUN = cluster_chunk)

# Combine cluster from all chunks (assuming all_cluster_labels is a list of vectors)
final_clusters <- unlist(all_cluster)

clusters <- cutree(final_clusters, k=5)

#plot the results
fviz_cluster(list(data=df, cluster=clusters),
             geom = "point",
             ellipse.type='convex')
