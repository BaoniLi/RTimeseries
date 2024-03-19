library(trend)
library(utils)
library(raster)
library(sp)

######
setwd("E:/YANGTZE/MID/RH4")

dir = list.files(pattern = "tif")
fns = Sys.glob(dir)
YEAR <- length(fns)  #年数(等于.tif文件的个数)

######
nrow = 1384
ncol = 1711
p = array(0,dim=c(nrow,ncol,YEAR))
mkp = array(0,dim=c(nrow,ncol))
mkz = array(0,dim=c(nrow,ncol))


#读取栅格数据到三维数组
for (i in (1:YEAR))
{
  d = matrix(0,nrow,ncol)
  # raster1 <- paste('RH','idw',i,'tif',sep = ".") #####################################
  d[] <- raster(fns[i])
  p[,,i] <- d[]                        
}

#趋势检验
for (j in (1:nrow))
{
  for (k in (1:ncol))
  {
    timeseries = p[j,k,]
    if (sum(is.na(timeseries))>0){
      mkp[j,k] = NA
      mkz[j,k] = NA
    } else {
      mk <- mk.test(timeseries)
      mkp[j,k] <- mk$p.value
      mkz[j,k] <- mk$statistic
    }
  }
}


######
# write.csv(mkz,"E:/YANGTZE/MID/mkzRH_0.1.csv")
# write.csv(mkp,"E:/YANGTZE/MID/mkpRH_0.1.csv")


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
    if (is.na(mkp[j,k])|is.na(mkz[j,k])){
      trendpic[j,k] = 2
    } else if (mkp[j,k]<pvalue2 & mkz[j,k]<0){
      trendpic[j,k] = 0
    } else if (mkp[j,k]>=pvalue2 & mkz[j,k]<0){
      trendpic[j,k] = 1
    } else {
      trendpic[j,k] = 2
    }
  }
}

setwd("E:/YANGTZE/MID")
######
GridTopology1km <- GridTopology(cellcentre.offset = c(110.2526646,25.98056863), cellsize = c(0.004811252243,0.004811252243),cells.dim = c(ncol,nrow))
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
newname <- paste('RH11111111','trend','tif',sep = ".")
writeRaster(ET05_1km.raster,newname)