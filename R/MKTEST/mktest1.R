library(trend)
library(utils)
library(raster)


setwd("E:/YANGTZE/ET0andRH/VPD5")


dir = list.files(pattern = "tif")
fns = Sys.glob(dir)
YEAR <- length(fns)  #年数(等于.tif文件个数)
nrow = 1357
ncol = 758
t = nrow*ncol


p = matrix(0,t,YEAR) #每个栅格的时间序列是一行
mkp1 = matrix(0,t,1)
mkz1 = matrix(0,t,1)
trendpic = matrix(0,nrow,ncol)


#读取插值数据
for (i in (1:YEAR))
{
  d = matrix(0,nrow,ncol)
  d[] <- raster(fns[i])
  dim(d) <- c(t,1)      #将多行多列数据转换为一列
  p[,i] <- d[,1]
}


#趋势分析
for (l in (1:t))
{
  timeseries <- p[l,]
  mk <- mk.test(timeseries)
  mkp1[l,1] = mk$p.value
  mkz1[l,1] = mk$statistic
}
dim(mkp1) <- c(nrow,ncol)
dim(mkz1) <- c(nrow,ncol)


write.csv(mkp1,"mkpVPD51.csv")
write.csv(mkz1,"mkzVPD51.csv")