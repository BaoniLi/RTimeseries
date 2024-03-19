library(utils)
library(raster)


######
setwd("E:/YANGTZE/YRD/Q4")
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

pthresh = 0.1
cor_p = matrix(0,nrow,ncol)
v1 <- c(1)
v2 <- c(1)

#每个栅格做两种变量时间序列的相关分析
for (j in (250:500))
{
  for (k in (250:500))
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
    if (is.na(corvalue[j,k])|is.na(pvalue[j,k])){
      cor_p[j,k] = 2
    } else if (pvalue[j,k]<pthresh & corvalue[j,k]>0){
      cor_p[j,k] = 0
      v1 <- c(v1,timeseries1)
      v2 <- c(v2,timeseries2)
    } else if (pvalue[j,k]<pthresh & corvalue[j,k]<0){
      cor_p[j,k] = 1
    } else {
      cor_p[j,k] = 2
    }
  }
}


######
write.csv(v1,"E:/YANGTZE/scatter/YRDQ.csv")
write.csv(v2,"E:/YANGTZE/scatter/YRDQET.csv")