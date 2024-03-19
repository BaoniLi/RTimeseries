library(utils)
library(raster)
library(sp)

######
setwd("E:/YANGTZE/ET0andRH")
mkp <- read.csv("mkpEA5.csv",header=F)
mkz <- read.csv("mkzEA5.csv",header=F)
pvalue = 0.05


nro <- nrow(mkp)
nco <- ncol(mkp)
trendpic = matrix(0,nro,nco)


#画趋势图
for (j in (1:nro))
{
  for (k in (1:nco))
  {
    if (mkp[j,k]<pvalue & mkz[j,k]>0){
      trendpic[j,k] = 0
    } else if (mkp[j,k]<pvalue & mkz[j,k]<0){
      trendpic[j,k] = 1
    } else {
      trendpic[j,k] = 2
    }
  }
}


# #写入栅格，注意排列顺序
# template_raster <- raster("E:/YANGTZE/YZDEM/subtif1km/s5.tif")
# trendpicraster = raster(trendpic,template_raster)
# plot(trendpicraster)
# newname <- paste('ET05','trend','tif',sep = ".")
# writeRaster(trendpicraster,newname)

######
GridTopology1km <- GridTopology(cellcentre.offset = c(115.5953, 24.45498), cellsize = c(0.00833,0.00833),cells.dim = c(nco,nro))
#cellcentre.offset：左下角栅格中心坐标（经度、纬度）；cellsize：栅格尺寸；dim：列数*行数；ascii为左下角栅格的左下角顶点坐标（最小经度最小纬度）
# Sat.cor <- as.data.frame(SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(vector(mode = "numeric",length = nro*nco))))
# Sat.cor <- Sat.cor[,-1]
# colnames(Sat.cor) <- c("xlon","ylat")

######
SatBase.grid <- read.asciigrid("s5ascii.txt") #ascii顺序是左上角开始，一行一行
Sat.base <- as.data.frame(SatBase.grid)
colnames(Sat.base) <- c("ascii","xlon","ylat")
Sat.base <- Sat.base[order(Sat.base$xlon,-Sat.base$ylat),] #改为和R矩阵一样，左上角开始，一列一列
#SatBase.raster <- raster(Sat.base,layer=1,values=T)
#plot(SatBase.raster)
#Sat.base1 <- merge(Sat.cor,Sat.base,by.x=c("xlon","ylat"),by.y=c("xlon","ylat"),all.x=T)



dim(trendpic) <- c(nro*nco,1)
trendpic1 <- as.data.frame(trendpic)
trendpic1 <- cbind(Sat.base,trendpic1)



######
ET05_1km.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(trendpic1$V1))#data须是数据框数据,填充进grid里
ET05_1km.raster <- raster(ET05_1km.grid,layer=1,values=T)
#plot(ET05_1km.raster)
newname <- paste('EA5','trend','tif',sep = ".")
writeRaster(ET05_1km.raster,newname)