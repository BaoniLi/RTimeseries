library(raster)
library(sp)
library(rgdal)
library(gstat)
library(maptools)

#设置工作空间
setwd("E:/YANGTZE/UP/VPD1")

#读取站点数据及待插值数据
gauge_data_all<-na.omit(read.table("E:/YANGTZE/UP/VPD1.txt"))
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
  date = i+1964 #输出年份
  newname <- paste(date,'VPD','tif',sep = ".")
  writeRaster(idw,newname)
  str1 = paste('第',date,'年插值完毕',sep='')
  print(str1)
}

# #将结果读取到csv表格.
# library(utils)
# dir=list.files(pattern = "tif")
# fns = Sys.glob(dir)
# YEAR <- length(fns)  #年数(等于生成的.tif文件的个数)
# p = matrix(0,25000000,YEAR)    #总共12*10个栅格，每一列存放一年的数据(剔除掉NA值即可)
# 
# for (i in (1:YEAR)) 
# {
#   d = matrix(0,5000,5000)
#   d[]<- raster(fns[i])
#   dim(d) <- c(nrow(d)*ncol(d),1)      #将多行多列数据转换为一列
#   p[,i]<-d[,1]                        #把每一列添加到汇总表中
# }
# write.csv(p,"D:/长江流域/griddata/111.csv") #生成的表格中把NA行删除，求均值即可得到面平均降雨量