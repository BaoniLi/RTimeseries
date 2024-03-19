library(sp)
library(lattice)
library(spgwr)
library(gstat)
library(ncdf4)
library(GWmodel)
###读取数据
rain <- read.delim("F:/R/P.txt")
rainfall<-read.delim("F:/R/P03.txt")
point <- read.delim("F:/R/point5.txt")

coordinates(rainfall)=~X+Y
coordinates(point)=~XLon+Ylat

####################################
#    读取TRMM数据及重采样          #
####################################
dir = "F:/TRMM1/*.nc4" ##批量读取TRMM数据
fns = Sys.glob(dir)

P=matrix(0,260,2557)   #存放最终的P数据，一列存放一个文件
P2=matrix(0,260,8)  #存放每天的数据

###读取3小时数据
for (i in (0:2556)) 
{
  for(j in (1:8))
  {
  data<-nc_open(fns[8*i+j]) 
  P1 <- ncvar_get(data,"precipitation") #读取降雨数据
  dim(P1) <- c(nrow(P1)*ncol(P1),1)      #将多行多列数据转换为一列
  for(n in (1:260)){if(is.na (P1[n])) {P1[n]<-0}}
  
  P2[,j]<-P1[,1]*3
  nc_close(data)
  }
P[,i+1]<-apply(P2,1,sum)
}
P1[2,1]=3
###读取天数据
for (i in (1:2557)) {
  data<-nc_open(fns[i]) 
  P1 <- ncvar_get(data,"precipitation") #读取降雨数据
  dim(P1) <- c(nrow(P1)*ncol(P1),1)      #将多行多列数据转换为一列
  P[,i]<-P1[,1]
  nc_close(data)
}

##获得坐标信息
data<-nc_open(fns[1])
v1 <- data$var[[1]] #获取变量值
varsize <- v1$varsize  #获取变量各维数的大小
nx <- varsize[1] #分别赋值，获得各维数大小
ny<- varsize[2]
Xval<-matrix(0,nx,1)  #经纬度
Yval<-matrix(0,ny,1)

for( i in 1:nx )       #获取经度坐标
  Xval[i]<- ncvar_get( data, v1$dim[[1]]$name, start=i, count=1 )

for( i in 1:ny )      #获取维度坐标
  Yval[i]<- ncvar_get( data, v1$dim[[2]]$name, start=i, count=1 )

#获取输出数据的经纬度,R语言按列读取数据
n=1
coord<-matrix(0,nx*ny,2)
for(j in 1:ny)
  for (i in 1:nx)
  {
    coord[n,1]=Xval[i]      #存放经度数据
    coord[n,2]=Yval[j]      #存放纬度数据
    n=n+1
  }

for(i in 1:ny)
  for (j in 1:nx)
  {
    coord[n,1]=Yval[i]      #存放经度数据
    coord[n,2]=Xval[j]      #存放纬度数据
    n=n+1
  }

write.csv(coord,"E:/coord.csv")
write.csv(P,"E:/P.csv")


###验证TRMM精度
TRP1 <- read.delim("F:/R/TRMM3B42RT_D.txt") #验证原始和5km数据应该采用不同的网格点
TR<-TRP1[,3:2559]
TR<-TRP2[,c(1:2557)]
ss.PTR<-matrix(0,nrow(rain),nsta)

eva.PTR<-matrix(0,nrow(rain),6)
s.PTR<-matrix(0,nrow(rain),nsta)


for (j in (1:nrow(rain)))#空间精度
{
  for (i in (1:nsta))
  {
  
    s.PTR[j,i]<-rain[j,i]
    ss.PTR[j,i]<-TR[dis.min[i,1],j]
  }
    eva.PTR[j,1]<-sum(abs(s.PTR[j,]-ss.PTR[j,]))/nsta  #平均绝对值误差
    
    if(mean(s.PTR[j,])>0.5) #只对流域平均降雨量大于0.5mm的进行统计
    {
      eva.PTR[j,2]<-eva.PTR[j,1]/mean(s.PTR[j,])       #相对绝对值偏差
      eva.PTR[j,3]<-cor(s.PTR[j,],ss.PTR[j,],method = "pearson") #计算线性相关系数
    }
}

for (i in (1:nsta))   #时序精度
{
  
  eva.PTR[i,4]<-sum(abs(TR[dis.min[i,1],]-rain[,i]))/nrow(rain)
  eva.PTR[i,5]<-eva.PTR[i,4]/mean(rain[,i])
  eva.PTR[i,6]<-cor(t(TR[dis.min[i,1],]),rain[,i],method = "pearson") 
}
write.csv(eva.PTR,"E:/evaPTR.csv")

TR_OUT<-TR[dis.min,]
write.csv(TR_OUT,"E:/TRMM5km_54.csv")  #输出网格点对应的卫星数据

write.csv(TRP2,"E:/TRMM5km.csv")

#插值到5km*5km
TRP1 <- read.delim("F:/R/TRMM3B42RT_D.txt")
point <- read.delim("F:/R/point5.txt")

dis<-matrix(0,nrow(TRP1),1)
TRP2<-matrix(0,nrow(point),ncol(TRP1)-2)

for (i in (1:nrow(point)))
{
  dis[,1]<-sqrt((point[i,2]-TRP1[,1])^2+(point[i,3]-TRP1[,2])^2) #计算距离 
  
  dis.ord<-order(dis[,1],decreasing=FALSE)[1:9]      #寻找最邻近的9个网格点
  dis.ord1<-(1/dis[dis.ord,1])/sum(1/dis[dis.ord,1]) #计算权重值
  
  TRP<-TRP1[dis.ord[],c(3:ncol(TRP1))]*dis.ord1
  TRP2[i,]<-apply(TRP,2,sum)
}

################################
#降雨数据融合
#############################################
rainfall<-read.delim("F:/R/P03.txt")
rain<-read.delim("F:/R/P.txt")
pointTR <- read.delim("F:/R/pointTR_D.txt") #卫星网格

pointTR[,2:5]<-pointTR[,2:5]/1000 
rainfall[,2:5]<-rainfall[,2:5]/1000 

nsta<-20     #定义雨量站点数据
#在54个雨量站中，随机选取nsta个站点
nsta1<-matrix(0,nsta,1)

kc <- kmeans(rainfall[,4:5], nsta)
for (i in (1:nsta))
{
  nsta1[i,1]<-which(kc$cluster==i)[1]
}
nsta1<-c(35,17,27,7,6,9,1,19,14,15,49,20,3,5,21,22,4,10,8,24)


#nsta1[,1]<-c(1:54)
rainfall<-rainfall[nsta1,]
rain<-rain[,nsta1]

dis<-matrix(0,nrow(pointTR),1)
dis.min<-matrix(0,nsta,1)  #存放最小值所在的距离

#寻找雨量站所在的网格点
for (i in (1:nsta)) 
{
dis[,1]<-sqrt((rainfall[i,2]-pointTR[,2])^2+(rainfall[i,3]-pointTR[,3])^2)  
dis.min[i,1]<-which.min(dis[,1])
}

#背景残差
P.res<-matrix(0,nrow(rain),nsta)
P.TR54<-matrix(0,nrow(rain),nsta)
for (j in (1:nrow(rain)))
{
  for (i in (1:nsta))
  {
    P.res[j,i]<-rain[j,i]-TRP2[dis.min[i,1],j]
    P.TR54[j,i]<-TRP2[dis.min[i,1],j]
  }
}

#############################################
#   交叉验证GWR，GWRK融合数据精度   #
#############################################
rain.res<-P.res
coordinates(rainfall)=~X+Y
rainfall.res<-rainfall

eva.Pgwr<-matrix(0,nrow(rain.res),6) #采用三个评价指标计算每天的值
eva.Pgwrk<-matrix(0,nrow(rain.res),6)
s.Pgwr<-matrix(0,nsta,2) #存放预测点和实测点降雨数据，1列为预测，2列为实测
s.Pgwrk<-matrix(0,nsta,2)

ss.Pgwr<-matrix(0,nrow(rain),nsta)
ss.Pgwrk<-matrix(0,nrow(rain),nsta)

for(j in (1:nrow(rain.res)))
{ 
 j=2400
  rainfall.res$P<-c(t(rain.res[j,]))
  rainfall$P<-c(t(rain[j,]))
 i=54
 if(length(which((rainfall$P==0)))/nsta<0.95) #对≠0雨量站>80%统计,不是对残差≠0进行统计
 {  
    for (i in 1:nsta)
    {
      rain.m<-rainfall.res[-i,]  #用于模拟模型的数据
      rain.P<-rainfall.res[i,] #用于交叉验证的数据点
      
      ##GWRK，GWR降雨融合 
      rain.dist<-gw.dist(dp.locat=coordinates(rain.m))
      rain.b<-bw.gwr(P ~ VX+VY , data=rain.m,approach="CV",adaptive=TRUE, 
                      kernel = "bisquare", dMat= rain.dist)
      rain.gwr<-gwr.basic(P ~ VX+VY , data=rain.m, bw=rain.b,adaptive=TRUE,
                          kernel = "bisquare", dMat=rain.dist)
      rain.pre<-gwr.predict(P ~ VX+VY , data=rain.m, predictdata=rain.P, bw=9, kernel="bisquare",
                            adaptive=T,  dMat2=rain.dist)
                       
      ##计算拟合残差
      res<-rain.m
      res$res<-rain.m$P-rain.gwr$SDF$yhat                #计算拟合残差
      #plot(variogram(res ~ 1, res), fit.vgm,xlab="半变差函数",ylab="滞后距离")
      #对残差进行OK插值
      null.vgm <- vgm(var(res$res), "Exp",max(variogram(res~1, res)$dist),mean(variogram(res~1, res)$gamma[1:3])) # 采用固定的模型参数
      
      #fit.vgm <- fit.variogram(variogram(res ~ 1, res),model=null.vgm)
      
      res_ok <- krige(res~1, locations=res, newdata=rain.P, model=null.vgm)
       
   
      s.Pgwr[i,1]<-rain.pre$SDF$prediction+P.TR54[j,i]
      s.Pgwrk[i,1]<-rain.pre$SDF$prediction+res_ok$var1.pred+P.TR54[j,i] #预测点数值
      s.Pgwrk[i,2]<-rain[j,i]
      
      if(s.Pgwr[i,1]<0) {s.Pgwr[i,1]<-0}
      if(is.na(s.Pgwrk[i,1])) {s.Pgwrk[i,1]<-0}
      if(s.Pgwrk[i,1]<0) {s.Pgwrk[i,1]<-0} 
      
      ss.Pgwr[j,i]<-s.Pgwr[i,1]
      ss.Pgwrk[j,i]<-s.Pgwrk[i,1] 
    }
    ###计算反映空间精度的评价指标
    eva.Pgwrk[j,1]<-sum(abs(s.Pgwrk[,2]-s.Pgwrk[,1]))/nsta  #平均绝对值误差
    eva.Pgwr[j,1]<-sum(abs(s.Pgwrk[,2]-s.Pgwr[,1]))/nsta  #平均绝对值误差
    
    if(mean(s.Pgwrk[,2])>0.5) #只对流域平均降雨量大于0.5mm的进行统计
    {
    eva.Pgwr[j,2]<-eva.Pgwr[j,1]/mean(s.Pgwrk[,2])       #相对绝对值偏差
    eva.Pgwr[j,3]<-cor(s.Pgwrk[,2],s.Pgwr[,1],method = "pearson") 
    
    eva.Pgwrk[j,2]<-eva.Pgwrk[j,1]/mean(s.Pgwrk[,2])       #相对绝对值偏差
    eva.Pgwrk[j,3]<-cor(s.Pgwrk[,2],s.Pgwrk[,1],method = "pearson") #计算线性相关系数
    }
    
    }
} 

for (i in (1:nsta))
{

  eva.Pgwrk[i,4]<-sum(abs(ss.Pgwrk[,i]-rain[,i]))/nrow(rain)
  eva.Pgwrk[i,5]<-eva.Pgwrk[i,4]/mean(rain[,i])
  eva.Pgwrk[i,6]<-cor(ss.Pgwrk[,i],rain[,i],method = "pearson") 
  
  eva.Pgwr[i,4]<-sum(abs(ss.Pgwr[,i]-rain[,i]))/nrow(rain)
  eva.Pgwr[i,5]<-eva.Pgwr[i,4]/mean(rain[,i])
  eva.Pgwr[i,6]<-cor(ss.Pgwr[,i],rain[,i],method = "pearson") 
}

write.csv(eva.Pgwrk,"E:/evaPgwrk.csv")
write.csv(eva.Pgwr,"E:/evaPgwr.csv")

###################################################
#              交叉验证GWRK，GWR插值数据精度      #
###################################################
rainfall<-read.delim("F:/R/P03.txt")
rain<-read.delim("F:/R/P.txt")

rainfall[,2:5]<-rainfall[,2:5]/1000 
coordinates(rainfall)=~X+Y

#在54个雨量站中，随机选取nsta个站点
rainfall<-rainfall[c(nsta1),]
rain<-rain[,c(nsta1)]

eva.gwr<-matrix(0,nrow(rain),6) #采用三个评价指标计算每天的值
eva.gwrk<-matrix(0,nrow(rain),6)
s.gwr<-matrix(0,nsta,2) #存放预测点和实测点降雨数据，1列为预测，2列为实测
s.gwrk<-matrix(0,nsta,2)

ss.gwr<-matrix(0,nrow(rain),nsta)
ss.gwrk<-matrix(0,nrow(rain),nsta)

for(j in (1:nrow(rain)))   #数字需要修改，nrow(rain)
{ 
  rainfall$P<-c(t(rain[j,]))
  
  if(length(which((rainfall$P==0)))/nsta<0.95) #只对降雨不为0的雨量站大于90%的天数进行统计
  {  
    for (i in 1:nsta)
    {
      rain.m<-rainfall[-i,]  #用于模拟模型的数据
      rain.P<-rainfall[i,] #用于交叉验证的数据点
      
      ##GWRK，GWR降雨融合 
      rain.dist<-gw.dist(dp.locat=coordinates(rain.m))
      rain.b<-bw.gwr(P ~ VX+VY , data=rain.m,approach="CV",adaptive=TRUE, 
                     kernel = "bisquare", dMat= rain.dist)
      rain.gwr<-gwr.basic(P ~ VX+VY , data=rain.m, bw=rain.b,adaptive=TRUE,
                          kernel = "bisquare", dMat=rain.dist)
      rain.pre<-gwr.predict(P ~ VX+VY , data=rain.m, predictdata=rain.P, bw=rain.b, kernel="bisquare",
                            adaptive=T,  dMat2=rain.dist)
      
      ##计算拟合残差
      res<-rain.m
      res$res<-rain.m$P-rain.gwr$SDF$yhat                #计算拟合残差
      
      #对残差进行OK插值
      null.vgm <- vgm(var(res$res), "Exp",max(variogram(res~1, res)$dist),mean(variogram(res~1, res)$gamma[1:3])) # 采用固定的模型参数
      res_ok <- krige(res~1, locations=res, newdata=rain.P, model=null.vgm)  
      
      s.gwr[i,1]<-rain.pre$SDF$prediction
      s.gwrk[i,1]<-rain.pre$SDF$prediction+res_ok$var1.pred #预测点数值
      s.gwrk[i,2]<-rain.P$P
      
      if(s.gwr[i,1]<0) {s.gwr[i,1]<-0}
      if(is.na(s.gwrk[i,1])) {s.gwrk[i,1]<-0}
      if(s.gwrk[i,1]<0) {s.gwrk[i,1]<-0} 
      
      ss.gwr[j,i]<-s.gwr[i,1]
      ss.gwrk[j,i]<-s.gwrk[i,1] 
  }
    ###计算评价指标
    eva.gwrk[j,1]<-sum(abs(s.gwrk[,2]-s.gwrk[,1]))/nsta  #平均绝对值误差
   
    eva.gwr[j,1]<-sum(abs(s.gwrk[,2]-s.gwr[,1]))/nsta  #平均绝对值误差
    if(mean(s.gwrk[,2])>0.5) #只对流域平均降雨量大于0.5mm的进行统计
    {
    eva.gwr[j,2]<-eva.gwr[j,1]/mean(s.gwrk[,2])       #相对绝对值偏差
    eva.gwr[j,3]<-cor(s.gwrk[,2],s.gwr[,1],method = "pearson") 
    
    eva.gwrk[j,2]<-eva.gwrk[j,1]/mean(s.gwrk[,2])       #相对绝对值偏差
    eva.gwrk[j,3]<-cor(s.gwrk[,2],s.gwrk[,1],method = "pearson") #计算线性相关系数
    }
  }
 
} 

for (i in (1:nsta))
{
  eva.gwrk[i,4]<-sum(abs(ss.gwrk[,i]-rain[,i]))/nrow(rain)
  eva.gwrk[i,5]<-eva.gwrk[i,4]/mean(rain[,i])
  eva.gwrk[i,6]<-cor(ss.gwrk[,i],rain[,i],method = "pearson") 
  
  eva.gwr[i,4]<-sum(abs(ss.gwr[,i]-rain[,i]))/nrow(rain)
  eva.gwr[i,5]<-eva.gwr[i,4]/mean(rain[,i])
  eva.gwr[i,6]<-cor(ss.gwr[,i],rain[,i],method = "pearson") 
}

write.csv(eva.gwrk,"E:/evagwrk.csv")
write.csv(eva.gwr,"E:/evagwr.csv")

######################################################################
#                      GWRK计算降水融合数据                          #
######################################################################
rainfall<-read.delim("F:/R/P03.txt")

rainfall[,2:5]<-rainfall[,2:5]/1000 
rainfall<-rainfall[c(nsta1),]
coordinates(rainfall)=~X+Y

mer.gwr<-matrix(0,3208,nrow(rain.res))
mer.gwrk<-matrix(0,3208,nrow(rain.res))

rain.res<-P.res
rainfall.res<-rainfall

rain.P<-point #网格点数据
rain.P[,2:5]<-rain.P[,2:5]/1000
coordinates(rain.P)=~XLon+Ylat

for (i in 1:nrow(rain.res))  #总共2588天降雨数据，03-04两年731天数据
{
 
  rainfall.res$P<-c(t(rain.res[i,]))
  rain.m<-rainfall.res  #实测站点数据
  
  ##GWRK，GWR降雨融合 
  rain.dist<-gw.dist(dp.locat=coordinates(rain.m))
  #rain.b<-bw.gwr(P ~ VX+VY , data=rain.m,approach="CV",adaptive=TRUE, 
             #kernel = "bisquare", dMat= rain.dist)
  rain.gwr<-gwr.basic(P ~ VX+VY , data=rain.m, bw=19,adaptive=TRUE,
                      kernel = "bisquare", dMat=rain.dist)
  rain.pre<-gwr.predict(P ~ VX+VY , data=rain.m, predictdata=rain.P, bw=19, kernel="bisquare",
                        adaptive=T,  dMat2=rain.dist)
  
  ##计算拟合残差
  res<-rain.m
  res$res<-rain.m$P-rain.gwr$SDF$yhat                #计算拟合残差
 
  #对残差进行OK插值
  null.vgm <- vgm(var(res$res), "Exp",max(variogram(res~1, res)$dist),mean(variogram(res~1, res)$gamma[1:3])) # 采用固定的模型参数
  res_ok <- krige(res~1, locations=res, newdata=rain.P, model=null.vgm)
  
 
  mer.gwr[,i]<-rain.pre$SDF$prediction+TRP2[,i]
  mer.gwrk[,i]<-rain.pre$SDF$prediction+res_ok$var1.pred+TRP2[,i] #预测点数值
  if(is.na(mer.gwrk[,i])) {mer.gwrk[,i]<-mer.gwr[,i]}
}

write.csv(mer.gwrk,"E:/mer.gwrk.csv")
write.csv(mer.gwr,"E:/mer.gwr.csv") 

###验证融合数据精度
s.rain<-matrix(0,nrow(rain),nsta)
s.mer<-matrix(0,nrow(rain),nsta)

eva.mer<-matrix(0,nrow(rain),3)

for (j in (1:nrow(rain)))
{
  for (i in (1:nsta))
  {
      s.rain[j,i]<-rain[j,i]
      s.mer[j,i]<-mer.gwr[dis.min[i,1],j]
  }
    eva.mer[j,1]<-sum(abs(s.mer[j,]-s.rain[j,]))/nsta  #平均绝对值误差
    
    if(mean(s.rain[j,])>0.5) #只对流域平均降雨量大于0.5mm的进行统计
    {
      eva.mer[j,2]<-eva.mer[j,1]/mean(s.rain[j,])       #相对绝对值偏差
      eva.mer[j,3]<-cor(s.rain[j,],s.mer[j,],method = "pearson") #计算线性相关系数
    }
}

write.csv(eva.mer,"E:/eva.mer.csv")
mean(eva.mer[,1])

mer.gwrk[which(mer.gwrk<0)]<-0

mer.gwrk1<-apply(mer.gwrk,2,mean)
write.csv(mer.gwrk1,"E:/mer.gwrk.csv")

mer.gwr[which(mer.gwr<0)]<-0

mer.gwr1<-apply(mer.gwr,2,mean)
write.csv(mer.gwr1,"E:/mer.gwr.csv")

##look for subcatchment
point.sub<-read.delim("F:/R/子流域网格点/zhangshu.txt")
num<-matrix(0,nrow(point.sub),1)
for(j in (1:nrow(point.sub)))
  num[j,1]<- which(point[,2]==point.sub[j,2]&point[,3]==point.sub[j,3])

mer.sub<-mer.gwr[num,]
mer.subk<-mer.gwrk[num,]

mer.sub1<-apply(mer.subk,2,mean)
write.csv(mer.sub1,"E:/mer.sub.csv")



#######################################################
#             GWRK插值获得赣江流域03-09面降雨量       #
#######################################################
rain<-read.delim("F:/R/P.txt")

rainfall<-read.delim("F:/R/P03.txt")
rainfall[,2:5]<-rainfall[,2:5]/1000 
coordinates(rainfall)=~X+Y

rain.P<-point #未知数据点
rain.P[,2:5]<-rain.P[,2:5]/1000
coordinates(rain.P)=~XLon+Ylat

inter.gwr<-matrix(0,3208,nrow(rain))
inter.gwrk<-matrix(0,3208,nrow(rain))

for (i in 1:nrow(rain))  #总共2557天降雨数据，03-04两年731天数据
{
  
  rainfall$P<-c(t(rain[i,]))
  rain.m<-rainfall  #实测站点数据

  ##GWRK，GWR降雨插值 
  rain.dist<-gw.dist(dp.locat=coordinates(rain.m))
  rain.b<-bw.gwr(P ~ VX+VY , data=rain.m,approach="CV",adaptive=TRUE, 
                kernel = "bisquare", dMat= rain.dist)
  rain.gwr<-gwr.basic(P ~ VX+VY , data=rain.m, bw=rain.b,adaptive=TRUE,
                      kernel = "bisquare", dMat=rain.dist)
  rain.pre<-gwr.predict(P ~ VX+VY , data=rain.m, predictdata=rain.P, bw=rain.b, kernel="bisquare",
                        adaptive=T,  dMat2=rain.dist)
  
  ##计算拟合残差
  res<-rain.m
  res$res<-rain.m$P-rain.gwr$SDF$yhat                #计算拟合残差
  
  #对残差进行OK插值
  null.vgm <- vgm(var(res$res), "Exp",max(variogram(res~1, res)$dist),mean(variogram(res~1, res)$gamma[1:3])) # 采用固定的模型参数
  res_ok <- krige(res~1, locations=res, newdata=rain.P, model=null.vgm)
  
  inter.gwr[,i]<-rain.pre$SDF$prediction
  inter.gwrk[,i]<-rain.pre$SDF$prediction+res_ok$var1.pred #预测点数值
  if(is.na(inter.gwrk[,i])) {inter.gwrk[,i]<-inter.gwr[,i]}
}

inter.gwrk[which(inter.gwrk<0)]<-0


inter.subk<-inter.gwrk[num,]

inter.sub1<-apply(inter.subk,2,mean)
write.csv(inter.sub1,"E:/inter.sub.csv")

inter.gwr<-apply(inter.gwr,2,mean)
write.csv(inter.gwrk,"E:/gwrk.csv")
write.csv(inter.gwr,"E:/gwr.csv")


