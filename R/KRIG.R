library(raster)
library(sp)
library(rgdal)
library(gstat)
library(maptools)

#设置工作空间
setwd("D:/datasets1/krigtmean")


#读取站点数据及待插值数据
gauge_data_all<-na.omit(read.table("D:/datasets1/krigtmean.txt"))
date_len = ncol(gauge_data_all) #获取总列数,用做循环条件
yzsta <- read.table("D:/datasets1/meteostatmean.txt",header=T)


for (i in seq(1,date_len,1)) {
  gauge_data = gauge_data_all[,i] # each column is a tif
   
  ##设定站点数据的投影为WGS84，先lon
  dsp <- SpatialPoints(yzsta[,2:3], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  dsp <- SpatialPointsDataFrame(dsp,yzsta)
  
  #此段话并未投影成平面坐标,但不建议修改,不影响结果
  WGS84<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#设置参考系WGS84
  dsp1<-spTransform(dsp,WGS84)#将经纬度转成平面坐标，使用WGS参考系
  #bound1<-spTransform(bound,WGS84)
  
  #控制插值分辨率,插入栅格模板
  template_raster<-raster("D:/CHM_PRE/krig.tif")
  #bound_raster<-rasterize(bound,template_raster)
  
  non_missing_indices <- !is.na(gauge_data)
  min_non_missing <- min(gauge_data[non_missing_indices], na.rm = TRUE)
  gauge_data[!non_missing_indices] <- min_non_missing
  
  # Calculate variogram
  v <- variogram(gauge_data ~ 1, data = dsp1)
  # plot(v,plot.number=T)
  v.fit<-fit.variogram(v,model=vgm(1,"Lin",0))
  #plot(v,v.fit)
  Grid<-as(template_raster,"SpatialGridDataFrame")#首先现将边界栅格转成空间网格
  kri<-krige(formula=gauge_data~1,model=v.fit,locations=dsp1,newdata=Grid,nmax=12, nmin=10)#location为已知点的坐标；newdata为需要插值的点的位置；nmax和nmin分别代表最多和最少搜索点的个数
  # spplot(kri["var1.pred"])
  
  date = i #输出年份
  newname <- paste(date,'tif',sep = ".")
  writeRaster(raster(kri["var1.pred"]),newname)
  
  #str1 = paste('第',date,'插值完毕',sep='')
  print(date)
}