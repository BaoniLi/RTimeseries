library(trend)
library(utils)
library(raster)


setwd("E:/YANGTZE/UP/ET")


nrow = 968
ncol = 1518

p = matrix(0,nrow,ncol)
q = matrix(0,nrow,ncol)
t = matrix(0,nrow,ncol)

name1 <- paste('2001','ET_500m','tif',sep = ".")
p[] = raster(name1)
name2 <- paste('2007','ET_500m','tif',sep = ".")
q[] = raster(name2)

for (i in (1:nrow))
{
  for (j in (1:ncol))
  {
    if (p[i,j]==0 | q[i,j]==0){
      t[i,j] = 0
    } else {
      t[i,j] = q[i,j]-p[i,j]
    }
  }
}


GridTopology1km <- GridTopology(cellcentre.offset = c(101.9422626,27.66348163), cellsize = c(0.004811252243,0.004811252243),cells.dim = c(ncol,nrow))
  
  
SatBase.grid <- read.asciigrid("E:/YANGTZE/UP/ascii.txt") 
Sat.base <- as.data.frame(SatBase.grid)
colnames(Sat.base) <- c("ascii","xlon","ylat")
Sat.base <- Sat.base[order(Sat.base$xlon,-Sat.base$ylat),] 
  
  
dim(t) <- c(nrow*ncol,1)
t1 <- as.data.frame(t)
t1 <- cbind(Sat.base,t1)
  
  
######
ET05_1km.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(t1$V1))
ET05_1km.raster <- raster(ET05_1km.grid,layer=1,values=T)
# plot(ET05_1km.raster)
tname <- paste('0107','tif',sep = ".")
writeRaster(ET05_1km.raster,tname)
str1 = paste('¼ÆËãÍê³É')
print(str1)