library(utils)
library(raster)


######
setwd("E:/YANGTZE/YRD/TEM4")
dir = list.files(pattern = "tif")
fns = Sys.glob(dir)
YEAR <- length(fns)  #文件夹里tif文件个数


######
nrow = 1345
ncol = 1496
p1 = array(0,dim=c(nrow,ncol,YEAR))
p2 = array(0,dim=c(nrow,ncol,YEAR))
corvalue = array(0,dim=c(nrow,ncol))
pvalue = array(0,dim=c(nrow,ncol))


#读取栅格数据到三维数组
for (i in (1:YEAR))
{
  d = matrix(0,nrow,ncol)
  # raster1 <- paste('RH','idw',i,'tif',sep = ".")
  d[] <- raster(fns[i])
  p1[,,i] <- d[]                        
}


######
setwd("E:/YANGTZE/YRD/ET1")
dir = list.files(pattern = "tif")
fns = Sys.glob(dir)
YEAR <- length(fns)  #文件夹里tif文件个数


#读取栅格数据到三维数组
for (i in (1:YEAR))
{
  d = matrix(0,nrow,ncol)
  # raster1 <- paste('TEM','idw',i,'tif',sep = ".")
  d[] <- raster(fns[i])
  p2[,,i] <- d[]                        
}


#每个栅格做两种变量时间序列的相关分析
for (j in (1:nrow))
{
  for (k in (1:ncol))
  {
    timeseries1 = p1[j,k,]
    timeseries2 = p2[j,k,]
    if (sum(is.na(timeseries1))>0|sum(is.na(timeseries2))>0){
      corvalue[j,k] = NA 
      pvalue[j,k] = NA
    } else {
      aa <- cor.test(timeseries1,timeseries2)
      corvalue[j,k] <- aa$estimate
      pvalue[j,k] <- aa$p.value
    }
  }
}


######
# write.csv(corvalue,"E:/YANGTZE/MID/RHVSTEMc.csv")
# write.csv(pvalue,"E:/YANGTZE/MID/RHVSTEMp.csv")


pthresh = 0.1
cor_p = matrix(0,nrow,ncol)


#画相关图
for (j in (1:nrow))
{
  for (k in (1:ncol))
  {
    if (is.na(corvalue[j,k])|is.na(pvalue[j,k])){
      cor_p[j,k] = 2
    } else if (pvalue[j,k]<pthresh & corvalue[j,k]>0){
      cor_p[j,k] = 0
    } else if (pvalue[j,k]<pthresh & corvalue[j,k]<0){
      cor_p[j,k] = 1
    } else {
      cor_p[j,k] = 2
    }
  }
}

# for (j in (1:nrow))
# {
#   for (k in (1:ncol))
#   {
#     if (pvalue[j,k]<pthresh & corvalue[j,k]>0){
#       cor_p[j,k] = 0
#     } else if (pvalue[j,k]<pthresh & corvalue[j,k]<0){
#       cor_p[j,k] = 1
#     } else {
#       cor_p[j,k] = 2
#     }
#   }
# }

######
GridTopology1km <- GridTopology(cellcentre.offset = c(115.7574026,28.01630263), cellsize = c(0.004811252243,0.004811252243),cells.dim = c(ncol,nrow))


######
SatBase.grid <- read.asciigrid("E:/YANGTZE/YRD/ascii.txt") 
Sat.base <- as.data.frame(SatBase.grid)
colnames(Sat.base) <- c("ascii","xlon","ylat")
Sat.base <- Sat.base[order(Sat.base$xlon,-Sat.base$ylat),] 


dim(cor_p) <- c(nrow*ncol,1)
dim(corvalue) <- c(nrow*ncol,1)
cor_p1 <- as.data.frame(cor_p)
cor_p1 <- cbind(Sat.base,cor_p1)
corvalue1 <- as.data.frame(corvalue)
corvalue1 <- cbind(Sat.base,corvalue1)

setwd("E:/YANGTZE/CORRELATION/YRD1")
######
ET05_1km.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(cor_p1$V1))
ET05_1km.raster <- raster(ET05_1km.grid,layer=1,values=T)
#plot(ET05_1km.raster)
newname <- paste('ETVSTEMP','tif',sep = ".")
writeRaster(ET05_1km.raster,newname)
ET05_1km.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(corvalue1$V1))
ET05_1km.raster <- raster(ET05_1km.grid,layer=1,values=T)
#plot(ET05_1km.raster)
newname <- paste('ETVSTEMC','tif',sep = ".")
writeRaster(ET05_1km.raster,newname)