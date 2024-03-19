library(trend)
library(utils)
library(raster)


setwd("E:/MOD/ET")


nrow = 2906
ncol = 5320
nyear = 17
nd = 45


for (k in (1:nyear))
{
  year = 2000+k
  p = matrix(0,nrow,ncol)
  q = matrix(0,nrow,ncol)
  t = matrix(0,nrow,ncol)
  counter = matrix(0,nrow,ncol)
  for (l in (0:nd))
  {
    day1 = 1+8*l
    day <- formatC(day1,flag = 0, width = 3)
    name1 <- paste(year,day,sep = "")
    name <- paste(name1,'ET_500m','tif',sep = ".")
    p[] <- raster(name)
    for (i in (1:nrow))
    {
      for (j in (1:ncol))
      {
        if (p[i,j] >= 10000){
           counter[i,j] = counter[i,j]+0
        } else{
          counter[i,j] = counter[i,j]+1
          q[i,j] = q[i,j]+p[i,j]
        }
      }
    }
    # print(q[1,1])
    # print(day1)
  }
  for (i in (1:nrow))
  {
    for (j in (1:ncol))
    {
      if (counter[i,j] == 0){
        t[i,j] = NA
      } else {
        t[i,j] = q[i,j]/counter[i,j]*46*0.1 ###summation of the year
      }
    }
  }
  # print(t[1,1])
  GridTopology1km <- GridTopology(cellcentre.offset = c(97.35861443,21.14707802), cellsize = c(0.004811252243,0.004811252243),cells.dim = c(ncol,nrow))
  
  SatBase.grid <- read.asciigrid("F:/MODIS/ET/ascii.txt") 
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
  tname <- paste(year,'ET_500m','tif',sep = ".")
  writeRaster(ET05_1km.raster,tname)
  str1 = paste('第',year,'年计算完成',sep='')
  print(str1)
}