library(raster)
library(sp)
library(rgdal)
library(gstat)
library(maptools)

#设置工作空间
setwd("E:/YANGTZE/UP/EA")

#读取站点数据及待插值数据
gauge_data_all<-na.omit(read.table("E:/YANGTZE/UP/EA.txt"))
date_len = ncol(gauge_data_all) #获取总列数,用做循环条件
yzsta <- read.table("E:/YANGTZE/YZstation/up.txt",header=T)

# #读取流域矢量边界(插值到空白栅格时需要掩膜裁剪)
# bound<-readOGR("YZDEM/subshp/s5.shp")
# #plot(bound,col="grey")

for (i in seq(1,date_len,1)) {
  gauge_data = gauge_data_all[,i]
  
  ##设定站点数据的投影为WGS84，先lon
  dsp <- SpatialPoints(yzsta[,2:3], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  dsp <- SpatialPointsDataFrame(dsp,yzsta)
  
  #此段话并未投影成平面坐标,但不建议修改,不影响结果
  WGS84<- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#设置参考系WGS84
  dsp1<-spTransform(dsp,WGS84)#将经纬度转成平面坐标，使用WGS参考系
  #bound1<-spTransform(bound,WGS84)
  
  #控制插值分辨率,插入栅格模板
  template_raster<-raster("E:/YANGTZE/YZDEM/up/template.tif")
  #bound_raster<-rasterize(bound,template_raster)
  
  #反距离权重(idw)插值
  gs <- gstat(formula=gauge_data~1, locations=dsp1,set=list(idp = 2))
  idw <- interpolate(template_raster, gs)
  #idwmask<-mask(idw,bound)
  #plot(idwmask)
  date = i #输出年份
  newname <- paste('EA','idw',date,'tif',sep = ".")
  writeRaster(idw,newname)
  str1 = paste('第',date,'年插值完毕',sep='')
  print(str1)
}
