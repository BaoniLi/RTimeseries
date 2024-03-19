library(trend)
library(utils)
library(raster)
library(sp)

######
setwd("E:/YANGTZE/YRD/EA5")

dir = list.files(pattern = "tif")
fns = Sys.glob(dir)
YEAR <- length(fns)  #年数(等于.tif文件的个数)

######
nrow = 1345
ncol = 1496
p = array(0,dim=c(nrow,ncol,YEAR))
coe = array(0,dim=c(nrow,ncol))
x = c(1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017)
nn = 53

#读取栅格数据到三维数组
for (i in (1:YEAR))
{
  d = matrix(0,nrow,ncol)
  # raster1 <- paste('RH','idw',i,'tif',sep = ".") #####################################
  d[] <- raster(fns[i])
  p[,,i] <- d[]                        
}

multi = matrix(0,nrow,ncol)
for (j in (1:nrow))
{
  for (k in (1:ncol))
  {
    for (l in (1:YEAR))
    {
      y = p[j,k,l]
      if (is.na(y)){
        multi[j,k] = NA
      } else {
        multi[j,k] = multi[j,k]+y
      }
    }
  }
}


trendpic = matrix(0,nrow,ncol)



for (j in (1:nrow))
{
  for (k in (1:ncol))
  {
    trendpic[j,k] = multi[j,k]/YEAR
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
newname <- paste('MULTIEA','tif',sep = ".")
writeRaster(ET05_1km.raster,newname)