library(trend)
library(utils)
library(raster)

setwd("E:/YANGTZE/annualtmean/s5")

dir = list.files(pattern = "tif")
fns = Sys.glob(dir)
YEAR <- length(fns)  #年数(等于.tif文件的个数)
nrow = 1357
ncol = 758


mkp = array(0,dim=c(nrow,ncol))
mkz = array(0,dim=c(nrow,ncol))
trendpic = array(0,dim=c(nrow,ncol))



for (j in (280:290))
{
  for (k in (750:758))
  {
    timeseries <- numeric(0)
    for (i in (1:YEAR))
    {
      d = matrix(0,nrow,ncol)
      d[] <- raster(fns[i])
      timeseries <- c(timeseries,d[j,k])
    }
    mk <- mk.test(timeseries)
    mkp[j,k] <- mk$p.value
    mkz[j,k] <- mk$statistic
  }
}
write.csv(mkz,"mkz2.csv")