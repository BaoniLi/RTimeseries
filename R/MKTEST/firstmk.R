library(trend)
library(utils)
library(raster)
library(sp)
library(rgdal)
library(gstat)
library(maptools)

######
setwd("E:/YANGTZE/YRD")
STA <- read.table("E:/YANGTZE/YZstation/yrd.txt", head=T)
indic <- read.table("RH.txt", head=F)
nsta = nrow(STA)
year = 17
######
mkp = matrix(0,nsta,1)
mkz = matrix(0,nsta,1)


#趋势检验
for (i in (1:nsta))
{
  timeseries = indic[i,]
  mk <- mk.test(timeseries)
  mkp[i,1] <- mk$p.value
  mkz[i,1] <- mk$statistic
}


#Z值插值
gauge_data_all <- mkz[,1]
date_len = ncol(gauge_data_all) #获取总列数,用做循环条件
for (i in seq(1,date_len,1)) {
  gauge_data = gauge_data_all[,i]
  
  ##设定站点数据的投影为WGS84，先lon
  dsp <- SpatialPoints(STA[,2:3], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  dsp <- SpatialPointsDataFrame(dsp,STA)
  
  #此段话并未投影成平面坐标,但不建议修改,不影响结果
  WGS84<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#设置参考系WGS84
  dsp1<-spTransform(dsp,WGS84)#将经纬度转成平面坐标，使用WGS参考系
  
  
  #控制插值分辨率,插入栅格模板
  template_raster<-raster("E:/YANGTZE/YZDEM/yrd/template.tif")
  
  
  #反距离权重(idw)插值
  gs <- gstat(formula=gauge_data~1, locations=dsp1,set=list(idp = 2))
  idw <- interpolate(template_raster, gs)
  date = i #输出年份
  newname <- paste('EA','idw',date,'tif',sep = ".")
  writeRaster(idw,newname)
  str1 = paste('第',date,'张插值完毕',sep='')
  print(str1)
}


nrow = nrow(template_raster)
ncol = ncol(template_raster)
d = matrix(0,nrow,ncol)
d[] = raster(template_raster)


pvalue1 = 0.05
pvalue2 = 0.1



trendpic = matrix(0,nrow,ncol)


#画趋势图
# for (j in (1:nrow))
# {
#   for (k in (1:ncol))
#   {
#     trendpic[j,k] = mkz[j,k]
#   }
# }
for (j in (1:nrow))
{
  for (k in (1:ncol))
  {
    if (mkp[j,k]<pvalue1 & mkz[j,k]<0){
      trendpic[j,k] = 0
    } else if (mkp[j,k]>=pvalue1 & mkp[j,k]<pvalue2 & mkz[j,k]<0){
      trendpic[j,k] = 1
    } else {
      trendpic[j,k] = 2
    }
  }
}

setwd("E:/YANGTZE/YRD")
######
GridTopology1km <- GridTopology(cellcentre.offset = c(115.7574026,28.01630263), cellsize = c(0.004811252243,0.004811252243),cells.dim = c(ncol,nrow))
#cellcentre.offset：左下角栅格中心坐标（经度、纬度）；cellsize：栅格尺寸；dim：列数*行数；ascii为左下角栅格的左下角顶点坐标（最小经度最小纬度）
# Sat.cor <- as.data.frame(SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(vector(mode = "numeric",length = nro*nco))))
# Sat.cor <- Sat.cor[,-1]
# colnames(Sat.cor) <- c("xlon","ylat")

######
SatBase.grid <- read.asciigrid("ascii.txt") #ascii顺序是左上角开始，一行一行
Sat.base <- as.data.frame(SatBase.grid)
colnames(Sat.base) <- c("ascii","xlon","ylat")
Sat.base <- Sat.base[order(Sat.base$xlon,-Sat.base$ylat),] #改为和R矩阵一样，左上角开始，一列一列
#SatBase.raster <- raster(Sat.base,layer=1,values=T)
#plot(SatBase.raster)
#Sat.base1 <- merge(Sat.cor,Sat.base,by.x=c("xlon","ylat"),by.y=c("xlon","ylat"),all.x=T)



dim(trendpic) <- c(nrow*ncol,1)
trendpic1 <- as.data.frame(trendpic)
trendpic1 <- cbind(Sat.base,trendpic1)



######
ET05_1km.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(trendpic1$V1))#data须是数据框数据,填充进grid里
ET05_1km.raster <- raster(ET05_1km.grid,layer=1,values=T)
#plot(ET05_1km.raster)
newname <- paste('RH14Z','trend','tif',sep = ".")
writeRaster(ET05_1km.raster,newname)