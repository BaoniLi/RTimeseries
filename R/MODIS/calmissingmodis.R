library(trend)
library(utils)
library(raster)


setwd("E:/MOD/MIDET")


nrow = 1384
ncol = 1711

name1 <- paste('extract_2001169','ET_500m','tif',sep = ".")
name2 <- paste('extract_2001185','ET_500m','tif',sep = ".")



p = matrix(0,nrow,ncol)
p[] <- raster(name1)
q = matrix(0,nrow,ncol)
q[] <- raster(name2)
t = matrix(0,nrow,ncol)
for (i in (1:nrow))
{
  for (j in (1:ncol))
  {
    t[i,j] = (p[i,j]+q[i,j])/2
  }
}


GridTopology1km <- GridTopology(cellcentre.offset = c(110.2526646,25.98056863), cellsize = c(0.004811252243,0.004811252243),cells.dim = c(ncol,nrow))


SatBase.grid <- read.asciigrid("E:/YANGTZE/MID/ascii.txt") 
Sat.base <- as.data.frame(SatBase.grid)
colnames(Sat.base) <- c("ascii","xlon","ylat")
Sat.base <- Sat.base[order(Sat.base$xlon,-Sat.base$ylat),] 


dim(t) <- c(nrow*ncol,1)
t1 <- as.data.frame(t)
t1 <- cbind(Sat.base,t1)


######
ET05_1km.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(t1$V1))
ET05_1km.raster <- raster(ET05_1km.grid,layer=1,values=T)
plot(ET05_1km.raster)
tname <- paste('extract_2001177','ET_500m','tif',sep = ".")
writeRaster(ET05_1km.raster,tname)