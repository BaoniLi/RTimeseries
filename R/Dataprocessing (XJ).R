library(sp)
library(maptools)
library(lattice)
library(spgwr)
library(rgdal)
library(gstat)
library(ncdf4)
library(GWmodel)
library(RNetCDF)
library(ggplot2)
library(dplyr)
library(gwrr)
library(fields)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_181')
library(rJava) 
library(xlsx)
library(raster)
library(rasterNA)
library(parallel)
library(iterators)
library(lubridate)  #日期提取年、月、日
library(Rmisc)

rm(list=ls())       #清空所有环境变量
WorkingDirectory<-"D:/Precipitation merging"  #默认的当前文件夹即可
setwd(WorkingDirectory)
Filepath<-choose.dir(paste(WorkingDirectory,"/Study area",sep=""))     #需选择子文件夹，即目标研究区数据文件夹（输入：雨量计数据及研究边界；输出：处理及融合数据），如“E:\Precipitation merging\Study area\西江流域"
ElapsedTime <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
time1<-Sys.time()

###################################
#           读取雨量站数据        #
###################################
pdata<- read.table(paste(Filepath,"/Input data/Raingauge data/TargetPre.txt",sep=""))    #读取之后为数据框数据
pdata[,3]<-ifelse(nchar(pdata[,3])==1,paste("0",pdata[,3],sep=""),pdata[,3])
pdata[,4]<-ifelse(nchar(pdata[,4])==1,paste("0",pdata[,4],sep=""),pdata[,4])     #将原始txt文件中字符数为1的前面加0，ifelse为向量化的函数
pdata$Date<-paste(pdata[,2],pdata[,3],pdata[,4],sep = "")
pdata$Date<-as.Date(pdata$Date,format='%Y%m%d')
pdata<-pdata[,-c(2:4)]        #删除第二列到第四列
names(pdata)[1:2]<-c("STCD","p")    #对data数据框的第一二列变量进行重命名
startdate<-as.Date("2014/09/30")
enddate<-as.Date("2014/09/30")

ndays<-enddate-startdate+2
print(ndays)
# tt<-ts(1:ndays,frequency = 1,start=as.Date("2014/07/05"))   #构造时间序列 
dates<-seq(from=startdate-1,by=1,length.out = ndays)    #产生按天的date数据
# tt<-data.frame(dates,tt)    #整合成时间序列格式
# head(tt)
Hour<-c("12","15","18","21","00","03","06","09")
nHours<-length(Hour)

###################################  
#       输入TRMM数据提取范围      #      
###################################  
xmn<-21.375
xmx<-26.875    #######最小最大纬度
ymn<-102.125
ymx<-111.875    #######最小最大经度
nx<-(xmx-xmn)/0.25+1     #纬度的维数
ny<-(ymx-ymn)/0.25+1     #经度的维数
nx0<-(xmn+49.875)/0.25+1  #起始栅格纬度的序数    注意：-49.875不同于日数据的-59.875（日数据在-60~-50，50~60范围内有数据，但均为NaN）
ny0<-(ymn+179.875)/0.25+1   #起始栅格经度的序数
##########TRMM RT#################
nx1<-(xmn+59.875)/0.25+1  #起始栅格纬度的序数    注意：RT数据和日数据的-59.875
ny1<-(ymn+179.875)/0.25+1   #起始栅格经度的序数

####################################
#           读取TRMM_RT数据        #
####################################
###读取、存储日数据
P_TRMM.RT=matrix(0,nx*ny,ndays-1)   #存放最终的P数据，一列存放一个文件
for (i in 2:ndays){
  P<-matrix(0,ny,nx)
  for (j in 1:nHours){
    if (Hour[j]>11){
      Date<-format(dates[i-1],'%Y%m%d')
    }
    else {
      Date<-format(dates[i],'%Y%m%d')
    }
    Time0<-as.POSIXct(paste(Date,Hour[j],"00",sep=""),format='%Y%m%d%H%M')   #此处Time为时间格式，因为paste之后的字符串不是时间的标准格式，故用format控制字符串（文本）转换成时间格式
    Time0<-format(Time0,"%Y-%m-%d %H:%M")    #利用format将时间输出为字符串，此处Time（左）为字符串格式
    path<-paste(WorkingDirectory,"/Raw TRMM3B42RT data/3B42RT.",Date,Hour[j],".7R2.nc4",sep="")
    if (file.exists(path)==FALSE){
      path<-paste(WorkingDirectory,"/Raw TRMM3B42RT data/3B42RT.",Date,Hour[j],".7.nc4",sep="")
    }
    data<-nc_open(path)
    P.TRR_3hour<- ncvar_get(data,"precipitation",c(ny1,nx1),c(ny,nx)) #读取降雨数据  P1 <- ncvar_get(data,"precipitation",c(325,1129),c(24,40))
    # write.csv(P.TRR_3hour,paste(Filepath,"/Output data/Processed TRMM RT data/3-hour data/TRMM_RT",Date," ",Hour[j],"时.csv",sep=""))
    P.TRR_3hour<-ifelse(is.na(P.TRR_3hour),0,P.TRR_3hour)
    P<-P+P.TRR_3hour*3
    nc_close(data)
  }
  P<- P[,nx:1]#将列顺序倒置，这与输出图的方向有关
  write.csv(P,paste(Filepath,"/Output data/Processed TRMM RT data/Daily  data converting from 3-hour data/TRMM_RT",Date,".csv",sep=""))   #P为40行24列
  #P<-t(P)   #因为有下步，故无需转置
  dim(P) <- c(nrow(P)*ncol(P),1)      #将多行多列数据转换为一列,原P为40行（经度）24列（纬度），转换的时候按列读取，相当于完成了转置。写入ascii文件时和用P转置（24行40列）写入结果一样。
  P[is.na(P)]=-9999 #将NA替换为-9999，否则ascii数据读入到GIS会出错。
  P_TRMM.RT[,i-1]<-P[,1]
  fileCon.TRMMRT <- paste(Filepath,"/Output data/Processed TRMM RT data/Daily  data converting from 3-hour data/ascii/TRMM_RT",Date,".asc",sep="")
  writeLines('ncols\t40\nnrows\t23\nxllcorner\t102\nyllcorner\t21.25\ncellsize\t0.25\nNODATA_value\t-9999', fileCon.TRMMRT)
  write.table(P_TRMM.RT[,i-1], fileCon.TRMMRT,append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)
}


###################################
#             basemap            #
###################################

###############################GSMaP格式转换后掩模数据(官网转换软件)###########################
# DATA1<-as.data.frame(read.asciigrid("I:/12.txt"))
# DATA2<-as.data.frame(read.asciigrid("I:/23.txt"))
# A1.grid<-read.asciigrid("I:/12.txt")
# A1.raster<-raster(A1.grid,layer=1,values=T)
# A11<-as.data.frame(A1.raster)
#A11和DATA1数据顺序相同，说明raster转为data.frame时也是从左上角栅格逐行读取

############   25km卫星基准面   #############
#######获取卫星基准面（25km）的raster数据以及相应网格中心点的经度、纬度、高程、坡向、坡度等信息#######
GridTopology25km<-GridTopology(cellcentre.offset = c(102.125, 21.375), cellsize = c(0.25,0.25),cells.dim = c(40,23))#cellcentre.offset：左下角栅格中心点坐标（经度、纬度）；cellsize：栅格尺寸；dim：列数*行数
# SatBase.grid <- SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(vector(mode = "numeric",length = nx*ny)))#data需是数据框数据,填充进grid里
# SatBase.raster<-raster(SatBase.grid,layer=1,values=T)
Sat.cor<-as.data.frame(SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(vector(mode = "numeric",length = nx*ny))))
Sat.cor<-Sat.cor[,-1]
colnames(Sat.cor)<-c("Xlon","Ylat")
SatBase.grid<-read.asciigrid(paste(Filepath,"/Input data/SatBase.txt",sep = ""))
SatBase.raster<-raster(SatBase.grid,layer=1,values=T)
plot(SatBase.raster)
Sat.base<-as.data.frame(SatBase.grid)
#rasterToPointsNA(SatBase.raster)
#Sat.base<-as.data.frame(TRMM.coord)
Sat.base<-Sat.base[,-1]
ngrids<-nrow(Sat.base)
colnames(Sat.base)<-c("Xlon","Ylat")


############    1km水文模拟基准面   #############
#######获取水文模拟基准面（1km）的raster数据以及相应网格中心点的经度、纬度、高程、坡向、坡度等信息#######
#此处设置的是cellcentre.offset即左下角栅格中心点，而ascii则是xllcorner和yllcorner即左下角像元的左下角
####不能这样写，因为和ascii的xllcorner和yllcorner有所不同
#GridTopology1km<-GridTopology(cellcentre.offset = c(xcent, ycent), cellsize = c(0.0102,0.0102),cells.dim = c(294,294))#cellcentre.offset：左下角栅格中心点坐标（经度、纬度）；cellsize：栅格尺寸；dim：列数*行数
HydBase.grid<-read.asciigrid(paste(Filepath,"/Input data/Hydbase.txt",sep=""))  #法二，运用sp包读取ASCII格式的栅格文件，读出为SpatialGridDataFrame格式  
HydBase.raster <- raster(HydBase.grid,layer=1,values=T)    #raster格式
GridTopology1km<-HydBase.grid@grid     #一个grid文件包括（1）data：“dataframe”；（2）grid：“GridTopology”；（3）bbox；（4）proj4string：Formal class“CRS" with 1 slot
plot(HydBase.raster)
#读入文件为ascii格式，其数据实际为一列数，base.raster中data值显示为854887行1列的数据框（由下面ascii写入数据为表格一列亦可知，看起来多行多列或是为了利用文本空间，实际可看成逐行读取后形成的一个数据向量），由ascii读入的数据经计算后可直接写入ascii，因为其栅格顺序与读入一致
#SpatialGridDataFrame格式文件由三部分子文件组成：data（数据框，用以存储值）、grid(GridTopology类的S4对象)和proj4string（CRS类的S4对象），后两个用以存储坐标信息
Hyd.base <-as.data.frame(HydBase.grid)    #SpatialGridDataFrame格式转为data.frame格式，将多个子文件（data信息和坐标信息）存储在一个数据框里（data和经纬度分成三列），像元由点表示，data.frame的第二列是经度，第三列是纬度
colnames(Hyd.base)<-c("Isin","Xlon","Ylat")   #对数据框列名重命名，Y经度，X纬度
Xlon_1km.grid<-SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(Hyd.base$Xlon))
Xlon_1km.raster<-raster(Xlon_1km.grid,layer=1,values=T)
Ylat_1km.grid<-SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(Hyd.base$Ylat))
Ylat_1km.raster<-raster(Ylat_1km.grid,layer=1,values=T)
dem_1km.grid<-read.asciigrid(paste(Filepath,"/Input data/dem1km.txt",sep=""))
dem_1km<-as.data.frame(dem_1km.grid)
colnames(dem_1km)<-c("dem","Xlon","Ylat") 
dem_1km.raster<-raster(dem_1km.grid,layer=1,values=T)
aspect_1km.grid<-read.asciigrid(paste(Filepath,"/Input data/aspect1km.txt",sep=""))
aspect_1km<-as.data.frame(aspect_1km.grid)
colnames(aspect_1km)<-c("aspect","Xlon","Ylat")  
aspect_1km.raster<-raster(aspect_1km.grid,layer=1,values=T)
slope_1km.grid<-read.asciigrid(paste(Filepath,"/Input data/slope1km.txt",sep=""))
slope_1km.raster<-raster(slope_1km.grid,layer=1,values=T)
slope_1km<-as.data.frame(rasterToPointsNA(slope_1km.raster))
colnames(slope_1km)<-c("Xlon","Ylat","slope")
slope_1km$slope<-ifelse(slope_1km$slope<1,1,slope_1km$slope)
slope_1km.grid<-SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(slope_1km$slope))
slope_1km<-as.data.frame(slope_1km.grid)
colnames(slope_1km)<-c("slope","Xlon","Ylat")
min(slope_1km$slope)
slope_1km.raster<-raster(slope_1km.grid,layer=1,values=T)

####DEM处理
r<-raster(nrow=366,ncol=600,xmn=101.77,xmx=112.75,ymn=21.07,ymx=27.59)
dem_2km.raster<-resample(dem_1km.raster,r,method="bilinear")
dem_1km.raster<-resample(dem_2km.raster,dem_1km.raster,method="bilinear")
#plot(dem_1km.raster)
slope_2km.raster<-resample(slope_1km.raster,r,method="bilinear")
slope_1km.raster<-resample(slope_2km.raster,slope_1km.raster,method="bilinear")
#plot(slope_1km.raster)

Hyd.top<-inner_join(Hyd.base,dem_1km,by=c("Xlon","Ylat"))
Hyd.top<-inner_join(Hyd.top,aspect_1km,by=c("Xlon","Ylat"))
Hyd.top<-inner_join(Hyd.top,slope_1km,by=c("Xlon","Ylat"))

dem_25km.raster<-resample(dem_1km.raster,SatBase.raster,method="bilinear")
dem_25km<-as.data.frame(rasterToPointsNA(dem_25km.raster) )
colnames(dem_25km)<-c("Xlon","Ylat","dem") 
slope_25km.raster<-resample(slope_1km.raster,SatBase.raster,method="bilinear")
slope_25km<-as.data.frame(rasterToPointsNA(slope_25km.raster) )
colnames(slope_25km)<-c("Xlon","Ylat","slope") 

####################################
#           读取CMORPH数据         #
####################################
###读取、存储日数据，最终得到的栅格在数据中和TRMM一一对应，坐标一致。
nx0_CMO<-(59.875-xmx)/0.25+1  #起始栅格纬度的序数 59.875（北纬）（上，第一个）~ -59.875（南纬）（下，最后一个）——————先读高纬度
####y只有西经时
ny0_CMO<-(ymn-0.125)/0.25+1    #起始栅格经度的序数   0.125（东经0.125）（左，第一个）~179.875（东经179.875，中间）~ -0.125（西经-179.875，中间）~ 358.875（西径-0.125）（右，最后一个）----先读低经度
P_CMO=matrix(0,nx*ny,ndays-1)   #存放最终的P数据，一列存放一个文件
for (i in 2:ndays){
  P<-matrix(0,ny,nx)
  for (j in 1:nHours){
    if (Hour[j]>11){
      Date<-format(dates[i-1],'%Y%m%d')
    }
    else {
      Date<-format(dates[i],'%Y%m%d')
    }
    path_CMO<-paste(WorkingDirectory,"/CMORPH 3-hour/cmorph.3hr-025deg.",Date,".nc",sep="")
    data<-nc_open(path_CMO)
    if (j<5){
      k<-j+4
    }
    else{k<-j-4}
    P.CMO_3hour<- ncvar_get(data,"cmorph_precip",c(ny0_CMO,nx0_CMO,k),c(ny,nx,1)) #读取降雨数据  P1 <- ncvar_get(data,"precipitation",c(325,1129),c(24,40))
    #P.CMO_3hour0<- ncvar_get(data,"cmorph_precip",c(1,1,k),c(1440,480,1))
    write.csv(P.CMO_3hour,paste(Filepath,"/Output data/Processed CMORPH data/3-hour data/CMO",Date," ",Hour[j],"时.csv",sep=""))
    #write.csv(P.CMO_3hour0,paste(Filepath,"/Output data/Processed CMORPH data/3-hour data/CMO",Date," ",Hour[j],"时0.csv",sep=""))
    #P.CMO_3hour0[is.na(P.CMO_3hour0)]=-9999 #将NA替换为-9999，要不然asc数据读入到GIS会出错
    P.CMO_3hour<-ifelse(is.na(P.CMO_3hour),0,P.CMO_3hour)
    P<-P+P.CMO_3hour*3
    nc_close(data)
  }
  write.csv(P,paste(Filepath,"/Output data/Processed CMORPH data/Daily  data converting from 3-hour data/CMO",Date,".csv",sep=""))   #P为40行24列
  #P<-t(P)   #因为有下步，故无需转置
  dim(P) <- c(nrow(P)*ncol(P),1)      #将多行多列数据转换为一列,原P为40行（经度）24列（纬度），转换的时候按列读取，相当于完成了转置。写入ascii文件时和用P转置（24行40列）写入结果一样。
  P[is.na(P)]=-9999 #将NA替换为-9999，否则ascii数据读入到GIS会出错。
  P_CMO[,i-1]<-P[,1]
  fileCon.CMO <- paste(Filepath,"/Output data/Processed CMORPH data/Daily  data converting from 3-hour data/ascii/CMO",Date,".asc",sep="")
  writeLines('ncols\t40\nnrows\t23\nxllcorner\t102\nyllcorner\t21.25\ncellsize\t0.25\nNODATA_value\t-9999', fileCon.CMO)
  write.table(P_CMO[,i-1], fileCon.CMO,append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)
}
write.csv(P_CMO,paste(Filepath,"/Output data/Processed CMORPH data/P_CMO.csv",sep=""))


####################################
#           读取CMORPH日数据         #
####################################
###读取、存储日数据，最终得到的栅格在数据中和TRMM一一对应，坐标一致。
# P_CMO=matrix(0,nx*ny,ndays-1)   #存放最终的P数据，一列存放一个文件
# for (i in 2:ndays){
#   Date<-format(dates[i],'%Y%m%d')
#   path_CMO<-paste(WorkingDirectory,"/CMORPH 2017/CMORPH",Date,".txt",sep="")
#   P_CMO[,i-1]<-as.matrix(read.table(path_CMO))
#   }

####################################
#          读取PERSIANN数据        #
####################################
P_PER=matrix(0,nx*ny,ndays-1)   #存放最终的P数据，一列存放一个文件
for (i in 2:ndays){
  P<-matrix(0,nx*ny,1)
  for (j in 1:nHours){
    if (Hour[j]>11){
      Date<-format(dates[i-1],'%Y%m%d')
    }
    else {
      Date<-format(dates[i],'%Y%m%d')
    }
    path_PER<-paste(WorkingDirectory,"/PERSIANN/PERSIANN_3h",Date,Hour[j],".asc",sep="")
    data.grid<-read.asciigrid(path_PER)    #SpatialGridDataFrame格式  
    data.raster <- raster(data.grid,layer=1,values=T)    #raster格式
   # data.raster <- mask(data.raster,PERBase.raster)   #掩模提取，会报错。此处掩模为将,若掩模对象某处值为NA，则将被掩模对象该处的值也设为NA（和gis相同），但二者范围和分辨率需一致（gis则无此要求）
    data.raster<-crop(data.raster,SatBase.raster)   #crop裁剪
    P.PER_3hour<-as.data.frame(data.raster)
    colnames(P.PER_3hour)<-"p"
    P.PER_3hour$p<-ifelse(is.na(P.PER_3hour$p),0,P.PER_3hour$p)
    plot(data.raster)
    P<-P+P.PER_3hour*3
  }
  P_PER[,i-1]<-P[,1]
}

####################################
#            读取GSMaP数据         #   ######需输入相应范围及经纬度上栅格的个数，注意：循环里的28*28也需要输入
####################################
GridTopology10km<-GridTopology(cellcentre.offset = c(102.05, 21.25), cellsize = c(0.1,0.1),cells.dim = c(100,58))#cellcentre.offset：左下角栅格中心坐标（经度、纬度）；cellsize：栅格尺寸；dim：列数*行数；ascii为左下角栅格的左下角顶点坐标（最小经度最小纬度）
# P_GSM=matrix(0,nx*ny,ndays-1)   #存放最终的P数据，一列存放一个文件
P_GSM=matrix(0,58*100,ndays-1)   #存放最终的P数据，一列存放一个文件
HourG<-c("12","13","14","15","16","17","18","19","20","21","22","23","00","01","02","03","04","05","06","07","08","09","10","11")
nHoursG<-length(HourG)
for (i in 2:ndays){
  P<-matrix(0,100*58,1)   #####10km*10km的数据个数
  for (j in 1:nHoursG){
    if (HourG[j]>11){
      Date<-format(dates[i-1],'%Y%m%d')
    }
    else {
      Date<-format(dates[i],'%Y%m%d')
    }
    path_GSM<-paste(Filepath,"/Input data/GSMaP_data/GSMaP_NRT",Date,HourG[j],".txt",sep="")
    P.GSM_1hour<-read.table(path_GSM)
    P.GSM_1hour$V1<-ifelse(P.GSM_1hour$V1<0,0,P.GSM_1hour$V1)   ####不加v1不行(数据框文件用ifelse要加变量名称，矩阵则不需要)
    P.GSM_1hour$V1<-ifelse(is.na(P.GSM_1hour$V1),0,P.GSM_1hour$V1) 
    #P.GSM_1hour <- sapply(P.GSM_1hour,function(x){if(x < 0){0} else {x}})   #一循环就报错
    P<-P+P.GSM_1hour*1
  }
  #对日降雨数据进行重采样、裁剪
  P.GSM_10km <- SpatialGridDataFrame(grid = GridTopology10km, data = as.data.frame(P))#data需是数据框数据,填充进grid里
  P.GSM_10km.raster<-raster(P.GSM_10km,layer=1,values=T)
 # P.GSM_25km.raster<-resample(P.GSM_10km.raster,SatBase.raster,method="bilinear")
  data.raster<-crop(P.GSM_10km.raster,SatBase.raster)   #crop裁剪
  plot(data.raster)
  data<-as.data.frame(P.GSM_10km.raster)
  P_GSM[,i-1]<-as.matrix(data)[,1]
}
###############GSMaP的10km分辨率################
GSM.grid<-SpatialGridDataFrame(grid = GridTopology10km, data = as.data.frame(matrix(1,58*100)))
# write.asciigrid(GSM.grid,paste(Filepath,"/Input data/GSMBase.txt",sep = ""))
GSM.cor<-as.data.frame(GSM.grid)
GSM.cor<-GSM.cor[,-1]
colnames(GSM.cor)<-c("Xlon","Ylat")
GSM.raster<-raster(GSM.grid,layer=1,values=T)
a1.raster<-resample(SatBase.raster,GSM.raster,method="bilinear")
a1.grid<-SpatialGridDataFrame(grid = GridTopology10km, data =as.data.frame(a1.raster))  #raster转grid
write.asciigrid(a1.grid,paste(Filepath,"/Input data/a1.txt",sep = ""))
GSMBase.grid<-read.asciigrid(paste(Filepath,"/Input data/GSMBase.txt",sep = "")) 
GSMBase.raster<-raster(GSMBase.grid,layer=1,values=T)
plot(GSMBase.grid)
plot(SatBase.grid)
GSM.base<-as.data.frame(GSMBase.grid)
GSM.base<-GSM.base[,-1]
ngrids_GSM<-nrow(GSM.base)
colnames(GSM.base)<-c("Xlon","Ylat")
dem_10km.raster<-resample(dem_1km.raster,GSMBase.raster,method="bilinear")
dem_10km<-as.data.frame(rasterToPointsNA(dem_10km.raster) )
colnames(dem_10km)<-c("Xlon","Ylat","dem") 
slope_10km.raster<-resample(slope_1km.raster,GSMBase.raster,method="bilinear")
slope_10km<-as.data.frame(rasterToPointsNA(slope_10km.raster) )
colnames(slope_10km)<-c("Xlon","Ylat","slope") 
###################################
#        遥测精度评价及融合       #
###################################
raingauges<-read.delim(paste(Filepath,"/Input data/Raingauge data/气象站点坐标.txt",sep=""))      #读取雨量站（坐标）信息
#head(raingauges)
#combined <- sort(union(levels(pdata$STCD), levels(raingauges$STCD)))
#pdata1 <- left_join(mutate(pdata, STCD=factor(STCD, levels=combined)),
               #mutate(raingauges, STCD=factor(STCD, levels=combined)),by="STCD")
pdata<-left_join(pdata,raingauges,by="STCD")
Nstations<-length(raingauges$STCD)
 #pdata2<-inner_join(pdata,raingauges,by="STCD")
Stas<-matrix(NA,ndays-1,7)
rnames<-dates[2:ndays]
cnames<-c(Nstations)
lambda.TRR<-matrix(NA,ndays-1,ngrids)
lambda.CMO<-matrix(NA,ndays-1,ngrids)
lambda.PER<-matrix(NA,ndays-1,ngrids)
lambda.GSM<-matrix(NA,ndays-1,ngrids_GSM)
cn.TRR<-matrix(NA,ndays-1,ngrids)
cn.CMO<-matrix(NA,ndays-1,ngrids)
cn.PER<-matrix(NA,ndays-1,ngrids)
cn.GSM<-matrix(NA,ndays-1,ngrids_GSM)

#######################2014年7月-8月只有27个雨量站########
# lambda<-matrix(NA,ndays-1,27)
# cn<-matrix(NA,ndays-1,27)
##########################################################
lambda<-matrix(NA,ndays-1,Nstations)
cn<-matrix(NA,ndays-1,Nstations)
rsquare1<-vector(mode = "numeric",length = ndays-1)
rsquare2<-vector(mode = "numeric",length = ndays-1)
# Stas.TRM<-matrix(NA,ndays-1,8)   #8是统计指标的个数，此处用了RMSE,MAE,BIAS,CC,ME,POD,POFD和ETS共8个统计指标
Stas.TRR<-matrix(NA,ndays-1,8) 
Stas.CMO<-matrix(NA,ndays-1,8) 
Stas.PER<-matrix(NA,ndays-1,8) 
Stas.GSM<-matrix(NA,ndays-1,8) 
Stas.GWRR<-matrix(NA,ndays-1,8) 
Stas.GWRRK<-matrix(NA,ndays-1,8) 

RMSE<-matrix(NA,ndays-1,5) 
CVscore<-matrix(NA,ndays-1,5) 
rsquare<-matrix(NA,ndays-1,5) 
MAE<-matrix(NA,ndays-1,5) 

bw.SRmerge<-vector(mode = "numeric",length = ndays-1)
CV.score <- vector(mode = "numeric",length = ndays-1)
rmse<-vector(mode = "numeric",length = ndays-1)
rmspe<-vector(mode = "numeric",length = ndays-1)
#lambda<-vector(mode = "numeric",length = ndays-1)

############################################
#      数据点和预测点距离矩阵求取公式      #    
############################################
#A和B均为数据点坐标矩阵，行为点，列为经度和纬度（两列）-----可加高程，加高程则为3列
#A为data point，B为predict point
#输出矩阵的转置是gwrr的输入
distAB <- function(A,B){
  smatAA<-apply(A,1,crossprod)
  matAA<-matrix(smatAA,nrow=nrow(A),ncol=nrow(B))
  smatBB<-apply(B,1,crossprod)
  matBB<-t(matrix(smatBB,nrow=nrow(B),ncol=nrow(A)))
  matAB<-tcrossprod(A,B)
  distAB<-sqrt(matAA+matBB-2*matAB)
  distAB
}

#######################################
#              降尺度公式            #    
#######################################

####dat为数据框文件，第一列至第七列，变量名称分别为Xlon、Ylat、ndvi、dem、aspect、slope和p
downscale <- function(dat,Sat.cor,Sat.base,GridTopology25km,GridTopology1km,SatBase.raster,HydBase.raster){
  #Sat.base_1km为SatBase.raster重采样至HydBase.raster像元位置处（且相同栅格尺寸）转换成dataframe格式，25km卫星数据经转换后，有值部分其isin属性值为1，否则为n
  cat("Downscale: ",substitute(dat))
  time<-Sys.time()
  Sat.base_1km<-resample(SatBase.raster,HydBase.raster,method="ngb")
  # plot(Sat.base_1km)
  Sat.base_1km<-rasterToPointsNA(Sat.base_1km)
  Sat.base_1km<-as.data.frame(Sat.base_1km)
  colnames(Sat.base_1km)<-c("Xlon","Ylat","isin")
  ngrids<-NROW(Sat.base)
  SAT<-filter(dat,p>0)
  # SAT<-filter(TRR,p>0)
  if (nrow(SAT)>=15){
    locs<-cbind(SAT$Xlon,SAT$Ylat)    ####wet condition栅格坐标（25km）
    SAT.BC<-SAT
    SAT.BC$p<-((SAT.BC$p)^0.25-1)/0.25
    if (nrow(SAT)<=2000){
      SAT.est0 <- gwrr.est(p ~ ndvi + Xlon + Ylat + dem + slope,locs, SAT.BC, kernel = "bis", bw = T)
    }else{
      SAT.est0 <- gwrr.est(p ~ ndvi + Xlon + Ylat + dem + slope,locs, SAT.BC, kernel = "bis", bw = 0.7)
    }
    
    # yhat_25km<-cbind(SAT.BC,SAT.est0$yhat)[,-3:-5]
    y<-(0.25*SAT.est0$yhat+1)^(1/0.25)
    yhat_25km<-cbind(SAT,y)[,-3:-5]
    colnames(yhat_25km)<-c("Xlon","Ylat","p","yhat")
    CN<-cbind(locs,SAT.est0$cn)
    LAMBDA<-cbind(locs,SAT.est0$lambda)
    CN<-as.data.frame(CN)
    LAMBDA<-as.data.frame(LAMBDA)
    colnames(CN)<-c("Xlon","Ylat","cn")
    colnames(LAMBDA)<-c("Xlon","Ylat","lambda")
    # CN<-left_join(Sat.base,CN,by=c("Xlon","Ylat"))
    # LAMBDA<-left_join(Sat.base,LAMBDA,by=c("Xlon","Ylat"))
    CN<-merge(Sat.base,CN,by.x=c("Xlon","Ylat"),by.y=c("Xlon","Ylat"),all.x=T)   ####all.x=T是”=“而非”==“
    LAMBDA<-merge(Sat.base,LAMBDA,by.x=c("Xlon","Ylat"),by.y=c("Xlon","Ylat"),all.x=T)
    CN<-CN[order(-CN$Ylat,CN$Xlon),]   ####merge改变了原有的排序方式（left_join则不改变）,需转成与ascii一致的排序方式 
    LAMBDA<-LAMBDA[order(-CN$Ylat,LAMBDA$Xlon),]
    cn<-CN$cn
    lambda<-LAMBDA$lambda
    beta_SAT<-cbind(locs,SAT.est0$beta[1,],SAT.est0$beta[2,],SAT.est0$beta[3,],SAT.est0$beta[4,],SAT.est0$beta[5,],SAT.est0$beta[6,])
    beta_SAT<-as.data.frame(beta_SAT)
    colnames(beta_SAT)<-c("Xlon","Ylat","beta0","beta.ndvi","beta.Xlon","beta.Ylat","beta.dem","beta.slope")
    beta_25km<-left_join(Sat.cor,beta_SAT,by=c("Xlon","Ylat"))
    beta0_25km <- SpatialGridDataFrame(grid = GridTopology25km, data =as.data.frame(beta_25km$beta0))#data需是数据框数据,填充进grid里
    beta.ndvi_25km <- SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(beta_25km$beta.ndvi))
    beta.Xlon_25km <- SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(beta_25km$beta.Xlon))
    beta.Ylat_25km <- SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(beta_25km$beta.Ylat))
    beta.dem_25km <- SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(beta_25km$beta.dem))
    beta.slope_25km <- SpatialGridDataFrame(grid = GridTopology25km, data = as.data.frame(beta_25km$beta.slope))
    cat("Downscale: ",paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_0_",substitute(dat),".txt",sep=""))
    write.asciigrid(beta0_25km,paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_0_",substitute(dat),".txt",sep=""))
    write.asciigrid(beta.ndvi_25km,paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_NDVI_",substitute(dat),".txt",sep=""))
    write.asciigrid(beta.Xlon_25km,paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_Xlon_",substitute(dat),".txt",sep=""))
    write.asciigrid(beta.Ylat_25km,paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_Ylat_",substitute(dat),".txt",sep=""))
    write.asciigrid(beta.dem_25km,paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_dem_",substitute(dat),".txt",sep=""))
    write.asciigrid(beta.slope_25km,paste(Filepath,"/Output data/Merged precipitation/Beta/Downscaled/",substitute(dat),"/beta_slope_",substitute(dat),".txt",sep=""))
    beta0_25km.raster <- raster(beta0_25km,layer=1,values=T)   #SpatialGridDataframe转成raster
    beta.ndvi_25km.raster <- raster(beta.ndvi_25km,layer=1,values=T)   #SpatialGridDataframe转成raster
    beta.Xlon_25km.raster <- raster(beta.Xlon_25km,layer=1,values=T)   #SpatialGridDataframe转成raster
    beta.Ylat_25km.raster <- raster(beta.Ylat_25km,layer=1,values=T)   #SpatialGridDataframe转成raster
    beta.dem_25km.raster <- raster(beta.dem_25km,layer=1,values=T)   #SpatialGridDataframe转成raster
    beta.slope_25km.raster <- raster(beta.slope_25km,layer=1,values=T)   #SpatialGridDataframe转成raster
    beta0_1km <- resample(beta0_25km.raster,HydBase.raster,method="bilinear")  #栅格重采样 re_raster <- resample ("geotiff",“参照的raster”,"方法”（“bilinear"或者”ngb“）)
    beta1_1km <- resample(beta.ndvi_25km.raster,HydBase.raster,method="bilinear")
    beta2_1km <- resample(beta.Xlon_25km.raster,HydBase.raster,method="bilinear")
    beta3_1km <- resample(beta.Ylat_25km.raster,HydBase.raster,method="bilinear")
    beta4_1km <- resample(beta.dem_25km.raster,HydBase.raster,method="bilinear")
    beta5_1km <- resample(beta.slope_25km.raster,HydBase.raster,method="bilinear")
    P_Hyd.raster_GWRR<-beta0_1km+beta1_1km*NDVI_1km.raster+beta2_1km*Xlon_1km.raster+beta3_1km*Ylat_1km.raster+beta4_1km*dem_1km.raster+beta5_1km*slope_1km.raster
# P_Hyd.raster_GWRR<-(0.25*P_Hyd.raster_GWRR+1)^(1/0.25) #BOX-COX逆变换
# plot(P_Hyd.raster_GWRR)
    points<-rasterToPoints(P_Hyd.raster_GWRR)  #仅wet condition,GWR预测降雨
    points<-as.data.frame(points)
    colnames(points)<-c("Xlon","Ylat","p")
    # max(points$p)
    coordinates(points)<-~Xlon+Ylat
    P_Hyd.raster_25km<-resample(P_Hyd.raster_GWRR,SatBase.raster,method="bilinear")  #wet condition的GWRR预测（1km）
    yhat<-rasterToPoints(P_Hyd.raster_25km)  #wet condition的1km预测提取至25km
    yhat<-as.data.frame(yhat)
    colnames(yhat)<-c("Xlon","Ylat","yhat")
    Res<-left_join(SAT.BC,yhat,by=c("Xlon","Ylat"))     #wet condition栅格（25km）预测值
    res<-Res$p-Res$yhat      #wet condition的GWRR预测残差
    Res<-cbind(Res,res)
    ######GSM可能会有个别对不上
    Res<-filter(Res,is.na(Res$res)==F)
    
    plot_den<-ggplot(Res,aes(x=res,alpha = 1/10))+geom_density()
    coordinates(Res)<-~Xlon+Ylat
    v<-variogram(Res$res~1,Res)
    v.fit<-fit.variogram(v,vgm("Sph"))
    if (v.fit[2,]$psill<=0) {
      v.fit<-vgm(mean(order(v$gamma,decreasing = FALSE)[1:5]),"Sph",max(v$dist),mean(order(v$gamma,decreasing = TRUE)))
    }
    if (v.fit[2,]$range<0){
      v.fit[2,]$range<-1
    }
    jpeg(paste(substitute(dat),Date,".jpeg",sep=""))
    plot(v,model=v.fit)
    dev.off()
    plot_vgm<-plot(v,model=v.fit)
    res_ok<-krige(res~1,Res,points,model=v.fit)    #wet condition的GWRR残（25km）差插值至1km网格
    res_ok<-as.data.frame(res_ok)
    points<-as.data.frame(points)
    yhat.GWRRK<-left_join(points,res_ok,by=c("Xlon","Ylat")) #points含GWRR预测降雨，res_ok含kriging插值残差
    yhat<-yhat.GWRRK$p+yhat.GWRRK$var1.pred    #wet condition的1kmGWRRK结果
    yhat.GWRRK$p<-(0.25* yhat.GWRRK$p+1)^(1/0.25) #BOX-COX逆变换,GWRR预测
    yhat<-(0.25*yhat+1)^(1/0.25) #BOX-COX逆变换，GWRRK预测
    yhat<-ifelse(yhat<0.1,0,yhat)   #GWRRK预测小于0.1的设为0
    #####前0.2%的最大值用第0.2%最大替代
    m<-floor(NROW(yhat)*0.005)
    yhat_threshold<-yhat[order(yhat,decreasing=TRUE)[1:m]][m]
    yhat.GWRR_threshold<-yhat.GWRRK$p[order(yhat.GWRRK$p,decreasing=TRUE)[1:m]][m]
    yhat<-ifelse(yhat >yhat_threshold, yhat_threshold,yhat)
    yhat.GWRRK$p<-ifelse(yhat.GWRRK$p > yhat.GWRR_threshold, yhat.GWRR_threshold,yhat.GWRRK$p)
    yhat<-ifelse(yhat > 600, yhat_threshold,yhat)
    yhat.GWRRK$p<-ifelse(yhat.GWRRK$p > 600, yhat.GWRR_threshold,yhat.GWRRK$p)
    yhat.GWRRK<-cbind(yhat.GWRRK,yhat)
    yhat.GWRRK$p<-ifelse(yhat.GWRRK$p<0.1,0,yhat.GWRRK$p)   #GWR预测小于0.1的设为0
    yhat1<-left_join(Sat.base_1km,yhat.GWRRK,by=c("Xlon","Ylat"))  #yhat1为HydBase.raster同区域（矩形，为栅格输出）同分辨率（1km）---用Sat.base_1km为了引入卫星有效数据区（isin属性）
    yhat1$p<-ifelse(is.na(yhat1$p)&yhat1$isin==1,0,yhat1$p)  #卫星有效区域内dry condition栅格GWRr预测值均为0
    yhat1$yhat<-ifelse(is.na(yhat1$yhat)&yhat1$isin==1,0,yhat1$yhat)  #卫星有效区域内dry condition栅格GWRRK预测值均为0
    yhat0.grid<-SpatialGridDataFrame(GridTopology1km,data = as.data.frame(yhat1$p))
    yhat0.raster<-raster(yhat0.grid,layer=1,values=T)  #GWRR结果
    yhat.grid<-SpatialGridDataFrame(GridTopology1km,data = as.data.frame(yhat1$yhat))
    yhat.raster<-raster(yhat.grid,layer=1,values=T)  #GWRRK结果
  } else {
    yhat_25km<-NULL
    plot_den<-NULL
    plot_vgm<-NULL
    yhat<-left_join(Sat.cor,dat,by=c("Xlon","Ylat"))
    yhat.grid<-SpatialGridDataFrame(GridTopology25km,data=as.data.frame(yhat$p))
    yhat.raster<-raster(yhat.grid,layer=1,values=T)
    yhat.raster<-resample(yhat.raster,HydBase.raster,method="bilinear")
    yhat0.raster<- yhat.raster
    cn<-c(rep(NA,ngrids))
    lambda<-c(rep(NA,ngrids))
  }
#plot(yhat.raster)
  yhat0<-rasterToPointsNA(yhat0.raster)
  yhat0<-as.data.frame(yhat0)
  names(yhat0)<-c("Xlon","Ylat","p")
  yhat<-rasterToPointsNA(yhat.raster)
  yhat<-as.data.frame(yhat)
  names(yhat)<-c("Xlon","Ylat","p")
  # yhat.grid<-SpatialGridDataFrame(grid = GridTopology1km, data =as.data.frame(yhat$p))
  # path_yhat<-paste(Filepath,"/Output data/P_hyd",Date,".asc",sep="")
  # write.asciigrid(P_Hyd.grid,path_yhat)
  cat(ElapsedTime(time))
  params <- list(yhat_25km,yhat0,yhat,cn,lambda,yhat.grid,plot_den,plot_vgm)
  # yhat25km为GWRR预测结果（25km分辨率），yhat0为GWRR预测结果，yhat为GWRRK预测结果，yhat.grid为GWRRK预测的ASCII格式结果
  names(params) <- c("yhat_25km","yhat0","yhat", "cn","lambda","yhat.grid","plot_den","plot_vgm")
  params
}


#######################################
#            精度评价指标             #    
#######################################
#dat为包含y和yhat数据的矩阵或数据框文件，y为观测，yhat为模型预测
statistic <- function(dat, y, yhat){
  N <- length(y)
  rmse<-sqrt(sum((y - yhat)^2) / N)  #Root-mean-square error
  mae<- sum(abs(y - yhat)) / N  #Mean absolute error
  bias<-(sum(yhat)-sum(y))/sum(y)    #此为相对偏差，实际为bias<-(sum(yhat)-sum(y))/N，即me
  cc<-cor(y,yhat,method = "pearson")
  me<-mean(y-yhat)   #Mean error
  h<-nrow(subset(dat,y>0.5&yhat>0.5))
  f<-nrow(subset(dat,y<=0.5&yhat>0.5))
  m<-nrow(subset(dat,y>0.5&yhat<=0.5))
  c<-nrow(subset(dat,y<=0.5&yhat<=0.5))
  pod<-h/(h+m)
  pofd<-f/(f+c)
  r<-(h+m)*(h+f)/(h+f+m+c)
  ets<-(h-r)/(h+f+m-r)
  c(rmse,mae,bias,cc,me,pod,pofd,ets)
}

#######################################
#  1、寻找雨量站所在TRMM.base的网格点 #     ########1、每站-每日，主要目的导出遥感的对应站点数据##########
#######################################
# for (j in 1:Nstations){
#   print(raingauges$STCD[j])
#   cnames[j]<-raingauges$STCD[j]
#   pdata$STCD
#   b<-filter(pdata,STCD==raingauges$STCD[j])
#   Dis[,1]<-sqrt((b$Xlon[1]-TRMM.coord[,1])^2+(b$Ylat[1]-TRMM.coord[,2])^2)
#   Dis.min[j,1]<-which.min(Dis[,1])      #寻找最小值所在TRMM网格中的排序,该位置和TRMM降雨位置对应
#   for (i in 2:ndays){
#     P_GAUGES.stations[i-1,j]<-b$p[i-1]
#     P_TRMM.stations[i-1,j]<-P_TRMM[Dis.min[j,1],i-1]
#     P_CMO.stations[i-1,j]<-P_CMO[Dis.min[j,1],i-1]
#   }
# }
# colnames(P_GAUGES.stations)<-cnames
# colnames(P_TRMM.stations)<-cnames
# write.csv(P_GAUGES.stations,paste(Filepath,"/Output data/Merged precipitation/P.Gauges.csv",sep=""),row.names=FALSE,col.names=T)  #矩阵输出为表格，以csv格式保存
# write.csv(P_TRMM.stations,paste(Filepath,"/Output data/Merged precipitation/P.TRMM.csv",sep=""),row.names=FALSE,col.names=T)
# write.csv(P_CMO.stations,paste(Filepath,"/Output data/Merged precipitation/P.CMO.csv",sep=""),row.names=FALSE,col.names=T)
a<-filter(pdata,Date==dates[i])
for (j in 2:ndays){
  j<-2
  time2<-Sys.time()
  
  Date<-format(dates[j],'%Y%m%d')
  aa<-filter(pdata,Date==dates[j])
  if (j<3){
    aaa<-aa
  }else{
    aaa<-rbind(aaa,aa)
  }
}
 for (i in 2:ndays){
   i<-2
   time2<-Sys.time()
   
  Date<-format(dates[i],'%Y%m%d')
  Date_time<-as.Date(Date,format = "%Y%m%d")
  Date0<-as.Date(paste(substr(Date,1,4),"01","01",sep="/"),format = "%Y/%m/%d")  #一年第一天日期
  Year<-year(Date_time)
  Date1<-as.Date(paste(substr(Date,1,6),"01",sep=""),format = "%Y%m%d")  #月初日期
  # DayOfYear<-Date_time-Date0+1    #计算日序数
  DayOfYear<-Date1-Date0+1    #月初日序数  
  #########################读取土地利用类型######################
  if (Year>2012){
    path_LandUse<-paste(WorkingDirectory,"/LandUse/LandUse2010.tif",sep="")
  }else if (Year>2007){
    path_LandUse<-paste(WorkingDirectory,"/LandUse/LandUse2010.tif",sep="")
  }else if(Year>2002){
    path_LandUse<-paste(WorkingDirectory,"/LandUse/LandUse2005.tif",sep="") 
  }else{
    path_LandUse<-paste(WorkingDirectory,"/LandUse/LandUse2000.tif",sep="")  
  }
  LandUse.raster<-crop(raster(path_LandUse),HydBase.raster)
  plot(LandUse.raster)
  area(LandUse.raster)     #####可显示栅格尺寸大小
  #######################根据上式显示的resolution输入至下式,需手动#############################
  GridTopology_LandUse<-GridTopology(cellcentre.offset = c(LandUse.raster@extent@xmin+0.5*0.01005009, LandUse.raster@extent@ymin+0.5*0.01005009), cellsize = c(0.01005009,0.01005009),cells.dim = c(LandUse.raster@ncols,LandUse.raster@nrows))#cellcentre.offset：左下角栅格中心点坐标（经度、纬度）；cellsize：栅格尺寸；dim：列数*行数
  LuBase.grid <- SpatialGridDataFrame( GridTopology_LandUse, data = as.data.frame(vector(mode = "numeric",length = LandUse.raster@ncols*LandUse.raster@nrows)))#data需是数据框数据,填充进grid里
  Lu.base<-as.data.frame(LuBase.grid)
  
  LandUse<-rasterToPoints(LandUse.raster)
  colnames(LandUse)<-c("Xlon","Ylat","landuse")
  LandUse<-as.data.frame(LandUse)
  
  #########################读取NDVI数据#########################
  #####R的else if一定要这样写（中括号及位置）
  if (DayOfYear<10) {
    DayOfYear<-paste("00",DayOfYear,sep="")
    } else if (DayOfYear<100) DayOfYear<-paste("0",DayOfYear,sep = "")
  path_NDVI<-paste(WorkingDirectory,"/NDVI_tif/MOD13A3_",Year,DayOfYear,".1_km_monthly_NDVI.tif",sep = "")
  if (file.exists(path_NDVI)==TRUE) {
  NDVI<-raster(path_NDVI)   #利用raster包读取tif数据
  NDVI.raster<-0.0001*crop(NDVI,HydBase.raster)   #crop裁剪
  NDVI.raster<-resample(NDVI.raster,LandUse.raster)
  NDVI0.raster<- NDVI.raster
  # plot(NDVI0.raster)
  NDVI<-rasterToPoints(NDVI.raster)
  colnames(NDVI)<-c("Xlon","Ylat","ndvi")
  NDVI<-as.data.frame(NDVI)
  NDVI<-left_join(NDVI,LandUse,by=c("Xlon","Ylat"))
  NDVI$ndvi<-ifelse(NDVI$landuse==1|is.na(NDVI$landuse),NDVI$ndvi,NA)  #利用ifelse进行条件赋值
  NDVI$ndvi<-ifelse(NDVI$ndvi<0,NA,NDVI$ndvi)  #利用ifelse进行条件赋值
  NDVI<-NDVI[,-4]     #删除landuse列
  NDVI.grid<-SpatialGridDataFrame( GridTopology_LandUse, data = as.data.frame(NDVI$ndvi))
  NDVI.raster<-raster(NDVI.grid,layer=1,values=T)
  # plot(NDVI.raster)
  aaa<-as.data.frame(NDVI.raster)
  # gf<-focalWeight(NDVI.raster,2,"Gauss")
  # is.matrix(gf)
  # NDVI1.raster<-focal(NDVI.raster,gf)
  NDVI1.raster<-focal(NDVI.raster,matrix(1/25,nrow=5,ncol=5),na.rm=T,pad=FALSE,padValue=NA,NAonly=T)
  # AAA<-as.data.frame(NDVI1.raster)
  # plot(NDVI1.raster)
  # bb<-NDVI1.raster-NDVI0.raster
  bbb<-as.data.frame(NDVI1.raster)
  bbb$layer<-ifelse(bbb$layer<0.1,0.1,bbb$layer)
  NDVI1.grid<-SpatialGridDataFrame(GridTopology_LandUse,data=as.data.frame(bbb$layer))
  NDVI1.raster<-raster(NDVI1.grid,layer=1,values=T)
  NDVI_1km.raster<-resample(NDVI1.raster,HydBase.raster,method="bilinear")
#############################################平滑化  
  NDVI_2km.raster<-resample(NDVI_1km.raster,r,method="bilinear")
  NDVI_1km.raster<-resample(NDVI_2km.raster,NDVI_1km.raster,method="bilinear")
  # plot(NDVI_1km.raster)
  # plot(slope_1km.raster)
  # plot(dem_1km.raster)
  # plot(Ylat_1km.raster)
  # plot(Xlon_1km.raster)
  
  NDVI_1km<-as.data.frame(rasterToPoints(NDVI_1km.raster))
  names(NDVI_1km)<-c("Xlon","Ylat","ndvi")
  NDVI_25km.raster<-resample(NDVI_1km.raster,SatBase.raster,method="bilinear")
  NDVI_25km<-as.data.frame(rasterToPoints(NDVI_25km.raster))
  names(NDVI_25km)<-c("Xlon","Ylat","ndvi")
  plot(NDVI_1km.raster)
  NDVI_10km.raster<-resample(NDVI_1km.raster,GSMBase.raster,method="bilinear")
  NDVI_10km<-as.data.frame(rasterToPoints(NDVI_10km.raster))
  names(NDVI_10km)<-c("Xlon","Ylat","ndvi")
  }
  Sat.top<-inner_join(Sat.base,dem_25km,by=c("Xlon","Ylat"))   #只剩能匹配到的
  # Sat.base<-merge(Sat.base,dem_25km,by.x = c("Xlon","Ylat"),by.y = c("Xlon","Ylat")) #只剩能匹配到的,同上
  # Sat.top<-inner_join(Sat.top,aspect_25km,by=c("Xlon","Ylat"))
  Sat.top<-inner_join(Sat.top,slope_25km,by=c("Xlon","Ylat"))
  min(Sat.top$slope)
  max(Sat.top$slope)
  Sat.top<-right_join(NDVI_25km,Sat.top,by=c("Xlon","Ylat"))
  ########################################################
  ####################     易出问题      ###################
  ########################################################
  # GSM.top0<-left_join(GSM.base,dem_10km,by=c("Xlon","Ylat"))   #只剩能匹配到的
  # head(GSM.base)
  # head(dem_10km)
  # GSM.top<-inner_join(GSM.base,dem_10km,by=c("Xlon","Ylat"))   #只剩能匹配到的
  # GSM.top<-inner_join(GSM.top,slope_10km,by=c("Xlon","Ylat"))
  # GSM.top<-right_join(NDVI_10km,GSM.top,by=c("Xlon","Ylat"))
  GSM.top<-merge(GSM.base,dem_10km,by.x=c("Xlon","Ylat"),by.y=c("Xlon","Ylat"))  #raster中的merge
  GSM.top<-merge(GSM.top,slope_10km,by.x=c("Xlon","Ylat"),by.y=c("Xlon","Ylat"))
  GSM.top<-merge(NDVI_10km,GSM.top,by.x=c("Xlon","Ylat"),by.y=c("Xlon","Ylat"),all.y=T)
  GSM.top<-GSM.top[order(-GSM.top$Ylat,GSM.top$Xlon),]
  # Sat.top<-Sat.top[,-3]
  TRR<-cbind(Sat.cor,P_TRMM.RT[,i-1])
  CMO<-cbind(Sat.cor,P_CMO[,i-1])
  PER<-cbind(Sat.cor,P_PER[,i-1])
  GSM<-cbind(GSM.cor,P_GSM[,i-1])
  colnames(TRR)<-c("Xlon","Ylat","p")
  colnames(CMO)<-c("Xlon","Ylat","p")
  colnames(PER)<-c("Xlon","Ylat","p")
  colnames(GSM)<-c("Xlon","Ylat","p")
  TRR<-left_join(Sat.top,TRR,by=c("Xlon","Ylat"))    #Sat.top为550个25km分辨率栅格，具有Xlon、Ylat、dem、aspect和slope等属性数据
  CMO<-left_join(Sat.top,CMO,by=c("Xlon","Ylat"))
  PER<-left_join(Sat.top,PER,by=c("Xlon","Ylat"))

  # GSM<-left_join(GSM.top,GSM,by=c("Xlon","Ylat"))
  GSM<-merge(GSM.top,GSM,by.x=c("Xlon","Ylat"),by.y=c("Xlon","Ylat"),all.x=T)
  GSM<-GSM[order(-GSM$Ylat,GSM$Xlon),]
  TRR.DS<-downscale(TRR,Sat.cor,Sat.base, GridTopology25km,GridTopology1km,SatBase.raster,HydBase.raster)
  CMO.DS<-downscale(CMO,Sat.cor,Sat.base, GridTopology25km,GridTopology1km,SatBase.raster,HydBase.raster)
  PER.DS<-downscale(PER,Sat.cor,Sat.base, GridTopology25km,GridTopology1km,SatBase.raster,HydBase.raster)
  GSM.DS<-downscale(GSM,GSM.cor,GSM.base, GridTopology10km,GridTopology1km,GSMBase.raster,HydBase.raster)
  write.asciigrid(TRR.DS$yhat.grid,paste(Filepath,"/Output data/Merged precipitation/Rainfall (Downscaled)/TRR/DSRainfall_",Date,".txt",sep=""))
  write.asciigrid(CMO.DS$yhat.grid,paste(Filepath,"/Output data/Merged precipitation/Rainfall (Downscaled)/CMO/DSRainfall_",Date,".txt",sep=""))
  write.asciigrid(PER.DS$yhat.grid,paste(Filepath,"/Output data/Merged precipitation/Rainfall (Downscaled)/PER/DSRainfall_",Date,".txt",sep=""))
  write.asciigrid(GSM.DS$yhat.grid,paste(Filepath,"/Output data/Merged precipitation/Rainfall (Downscaled)/GSM/DSRainfall_",Date,".txt",sep=""))
  jpeg(paste("Den",Date,".jpeg",sep=""))
  multiplot(TRR.DS$plot_den,CMO.DS$plot_den,PER.DS$plot_den,GSM.DS$plot_den,cols=2)
  dev.off()
  jpeg(paste("TRR",Date,".jpeg",sep=""))
  par(frow=c(2,2))
  CMO.DS$plot_vgm
  PER.DS$plot_vgm
  GSM.DS$plot_vgm
  TRR.DS$plot_vgm
  dev.off()
  jpeg(paste("CMO",Date,".jpeg",sep=""))
  par(frow=c(2,2))
  TRR.DS$plot_vgm
  PER.DS$plot_vgm
  GSM.DS$plot_vgm
  CMO.DS$plot_vgm
  dev.off()
  jpeg(paste("PER",Date,".jpeg",sep=""))
  par(frow=c(2,2))
  TRR.DS$plot_vgm
  CMO.DS$plot_vgm
  GSM.DS$plot_vgm
  PER.DS$plot_vgm
  dev.off()
  jpeg(paste("GSM",Date,".jpeg",sep=""))
  par(frow=c(2,2))
  TRR.DS$plot_vgm
  CMO.DS$plot_vgm
  PER.DS$plot_vgm
  GSM.DS$plot_vgm
  dev.off()
  yhat_25km.TRR<-TRR.DS$yhat_25km
  yhat_25km.CMO<-CMO.DS$yhat_25km
  yhat_25km.PER<-PER.DS$yhat_25km
  yhat_10km.GSM<-GSM.DS$yhat_25km
  if (is.null(yhat_25km.TRR)==FALSE){
    yhat_25km.TRR<-cbind(Date=c(as.character(Date_time) ),yhat_25km.TRR)   #add a column of Date 
  }
  if (is.null(yhat_25km.CMO)==FALSE){
    yhat_25km.CMO<-cbind(Date=c(as.character(Date_time) ),yhat_25km.CMO)   #add a column of Date 
  }
  if (is.null(yhat_25km.PER)==FALSE){
    yhat_25km.PER<-cbind(Date=c(as.character(Date_time) ),yhat_25km.PER)   #add a column of Date 
  }
  if (is.null(yhat_10km.GSM)==FALSE){
    yhat_10km.GSM<-cbind(Date=c(as.character(Date_time) ),yhat_10km.GSM)   #add a column of Date
  }
 
  if (i==2) {
    Allyhat_25km.TRR<-yhat_25km.TRR
    Allyhat_25km.CMO<-yhat_25km.CMO
    Allyhat_25km.PER<-yhat_25km.PER
    Allyhat_10km.GSM<-yhat_10km.GSM
  }else {
    Allyhat_25km.TRR<-rbind(Allyhat_25km.TRR,yhat_25km.TRR)
    Allyhat_25km.CMO<-rbind(Allyhat_25km.CMO,yhat_25km.CMO)
    Allyhat_25km.PER<-rbind(Allyhat_25km.PER,yhat_25km.PER)
    Allyhat_10km.GSM<-rbind(Allyhat_10km.GSM,yhat_10km.GSM)
  }
  cn.TRR[i-1,]<-TRR.DS$cn
  cn.CMO[i-1,]<-CMO.DS$cn
  cn.PER[i-1,]<-PER.DS$cn
  cn.GSM[i-1,]<-GSM.DS$cn
  lambda.TRR[i-1,]<-TRR.DS$lambda
  lambda.CMO[i-1,]<-CMO.DS$lambda
  lambda.PER[i-1,]<-PER.DS$lambda
  lambda.GSM[i-1,]<-GSM.DS$lambda
  
  #######################################
  #  2、寻找雨量站所在TRMM.base的网格点 #     ########2、每日-每站，主要目的是为每天的数据融合提供数据源##########
  #######################################
  a<-filter(pdata,Date==dates[i])
  nstations<-length(a$STCD)
  dis<-matrix(0,nrow(TRR.DS$yhat),1)
  disGSM<-matrix(0,nrow(GSM.DS$yhat),1)
  dis.min<-matrix(NA,nstations,1)  #存放最小值所在的距离
  disGSM.min<-matrix(NA,nstations,1)  #存放最小值所在的距离
  P_SR<-matrix(NA,nstations,13)    #Satellite rainfall and raingauge data
  for (j in 1:nstations){
    dis[,1]<-sqrt((a$Xlon[j]-TRR.DS$yhat[,1])^2+(a$Ylat[j]-TRR.DS$yhat[,2])^2)
    
    ##############################################################################
    ###############                                               ################
    ###############                                               ################
    ##############################################################################
    #### disGSM[,1]<-sqrt((a$Xlon[j]-TRR.DS$yhat[,1])^2+(a$Ylat[j]-TRR.DS$yhat[,2])^2)
    
    disGSM[,1]<-sqrt((a$Xlon[j]-GSM.DS$yhat[,1])^2+(a$Ylat[j]-GSM.DS$yhat[,2])^2)
    dis.min[j,1]<-which.min(dis[,1])      #寻找最小值所在TRMM网格中的排序,该位置和TRMM降雨位置对应
    disGSM.min[j,1]<-which.min(disGSM[,1])      #寻找最小值所在TRMM网格中的排序,该位置和TRMM降雨位置对应
    P_SR[j,2]<-a$STCD[j]
    P_SR[j,3]<-a$Xlon[j]
    P_SR[j,4]<-a$Ylat[j]
    P_SR[j,5]<-a$p[j]   #P_error[j,1]本应写成P_error[j,i],但2014年7月1日-2014年8月不少站点数据缺失，各时间节点雨量站数目不同
    P_SR[j,6]<-TRR.DS$yhat0$p[dis.min[j,1]]   #GWRR结果
    P_SR[j,7]<-CMO.DS$yhat0$p[dis.min[j,1]]
    P_SR[j,8]<-PER.DS$yhat0$p[dis.min[j,1]]
    P_SR[j,9]<-GSM.DS$yhat0$p[disGSM.min[j,1]]
    P_SR[j,10]<-TRR.DS$yhat$p[dis.min[j,1]]   #GWRRK结果
    P_SR[j,11]<-CMO.DS$yhat$p[dis.min[j,1]]
    P_SR[j,12]<-PER.DS$yhat$p[dis.min[j,1]]
    P_SR[j,13]<-GSM.DS$yhat$p[disGSM.min[j,1]]
  }
  P_SR<-as.data.frame(P_SR)
  names(P_SR)[1:13]<-c("Date","STCD","Xlon","Ylat","P_raingauge","P_TRR.GWRR","P_CMO.GWRR","P_PER.GWRR","P_GSM.GWRR","P_TRMM.RT","P_CMO","P_PER","P_GSM")    #对data数据框的第一二三列变量进行重命名
  P_SR[,1]<-as.character(Date_time)   
  write.csv(P_SR,paste(Filepath,"/Output data/Merged precipitation/P_SR.csv",sep=""))
 
    ###################################
    #            遥感产品融合         #
    ###################################
  Sat.base_1km<-resample(SatBase.raster,HydBase.raster,method="ngb")
  plot(Sat.base_1km)
  plot(HydBase.raster)
  Sat.cor_1km<-rasterToPointsNA(Sat.base_1km)
  Sat.cor_1km<-as.data.frame(Sat.cor_1km)
  colnames(Sat.cor_1km)<-c("Xlon","Ylat","isin")
  Sat.base1_1km<-rasterToPoints(Sat.base_1km)
  Sat.base1_1km<-as.data.frame(Sat.base1_1km)
  colnames(Sat.base1_1km)<-c("Xlon","Ylat","isin")
  Sat.merge<-rasterToPoints(Sat.base_1km)
  Sat.merge<-as.data.frame(Sat.merge)
  colnames(Sat.merge)<-c("Xlon","Ylat","isin")
  
    locs<- cbind(P_SR$Xlon,P_SR$Ylat)
    P_SR.BC<-P_SR
    P_SR.BC$P_raingauge<-((P_SR.BC$P_raingauge)^0.25-1)/0.25
    P_SR.BC$P_TRMM.RT<-((P_SR.BC$P_TRMM.RT)^0.25-1)/0.25
    P_SR.BC$P_CMO<-((P_SR.BC$P_CMO)^0.25-1)/0.25
    P_SR.BC$P_PER<-((P_SR.BC$P_PER)^0.25-1)/0.25
    P_SR.BC$P_GSM<-((P_SR.BC$P_GSM)^0.25-1)/0.25
    A<-cbind(P_SR$Xlon,P_SR$Ylat)
    B<-cbind(Sat.merge$Xlon,Sat.merge$Ylat)
    dp.dist<-distAB(A,B)
    #dp.dist1<-distAB(A,A)
    dis1<-matrix(0,nrow(Sat.merge),1)
    dis1.min<-matrix(NA,nstations,1)  #存放最小值所在的距离
    S1<-matrix(NA,nstations,nstations)
    for (j in 1:nstations){
      dis1[,1]<-sqrt((a$Xlon[j]-Sat.merge$Xlon)^2+(a$Ylat[j]-Sat.merge$Ylat)^2)
      dis1.min[j,1]<-which.min(dis1[,1])   
      S1[,j]<-dp.dist[,dis1.min[j,1]] #S1为雨量站点与最近1km网格中心点的距离矩阵，同雨量站间距离矩阵有微小差别
    }
    P_SR.est <- gwrr.est(P_raingauge ~ P_TRMM.RT  + P_CMO + P_PER + P_GSM,locs, P_SR.BC, kernel = "bis", bw = TRUE)
    P_SR.est1 <- gwrr.est(P_raingauge ~ P_TRMM.RT,locs, P_SR.BC, kernel = "bis", bw = TRUE)
    P_SR.est2 <- gwrr.est(P_raingauge ~ P_CMO,locs, P_SR.BC, kernel = "bis", bw = TRUE)
    P_SR.est3 <- gwrr.est(P_raingauge ~ P_PER,locs, P_SR.BC, kernel = "bis", bw = TRUE)
    P_SR.est4 <- gwrr.est(P_raingauge ~ P_GSM,locs, P_SR.BC, kernel = "bis", bw = TRUE)
    bw.SRmerge[i-1]<-P_SR.est$phi
    lambda[i-1,]<- P_SR.est$lambda
    cn[i-1,]<- P_SR.est$cn
    RMSE[i-1,]<-c(P_SR.est$RMSE,P_SR.est1$RMSE,P_SR.est2$RMSE,P_SR.est3$RMSE,P_SR.est4$RMSE)
    CVscore[i-1,]<-c(P_SR.est$cvResults,P_SR.est1$cvResults,P_SR.est2$cvResults,P_SR.est3$cvResults,P_SR.est4$cvResults)
    rsquare[i-1,]<-c(P_SR.est$rsquare,P_SR.est1$rsquare,P_SR.est2$rsquare,P_SR.est3$rsquare,P_SR.est4$rsquare)
    MAE[i-1,]<-c(P_SR.est$MAE,P_SR.est1$MAE,P_SR.est2$MAE,P_SR.est3$MAE,P_SR.est4$MAE)
    P_SR.est <- gwrr.predict(P_raingauge ~ P_TRMM.RT  + P_CMO + P_PER + P_GSM, P_SR.BC, P_SR.BC, kernel = "bis", bw = P_SR.est$phi, S1)
    P_SR.est1 <- gwrr.predict(P_raingauge ~ P_TRMM.RT , P_SR.BC, P_SR.BC, kernel = "bis", bw = P_SR.est1$phi, S1)
    P_SR.est2 <- gwrr.predict(P_raingauge ~ P_CMO , P_SR.BC, P_SR.BC, kernel = "bis", bw = P_SR.est2$phi, S1)
    P_SR.est3 <- gwrr.predict(P_raingauge ~ P_PER, P_SR.BC, P_SR.BC, kernel = "bis", bw = P_SR.est3$phi, S1)
    P_SR.est4 <- gwrr.predict(P_raingauge ~ P_GSM, P_SR.BC, P_SR.BC, kernel = "bis", bw = P_SR.est4$phi, S1)
    P_GWRR<-(0.25*P_SR.est$yhat+1)^(1/0.25)
    P_GWRR1<-(0.25*P_SR.est1$yhat+1)^(1/0.25)
    P_GWRR2<-(0.25*P_SR.est2$yhat+1)^(1/0.25)
    P_GWRR3<-(0.25*P_SR.est3$yhat+1)^(1/0.25)
    P_GWRR4<-(0.25*P_SR.est4$yhat+1)^(1/0.25)
    P_SR<-cbind(P_SR,P_GWRR,P_GWRR1,P_GWRR2,P_GWRR3,P_GWRR4)
    P_SR$P_GWRR<-ifelse( P_SR$P_GWRR<0,0, P_SR$P_GWRR)
    P_SR$P_GWRR1<-ifelse( P_SR$P_GWRR1<0,0, P_SR$P_GWRR1)
    P_SR$P_GWRR2<-ifelse( P_SR$P_GWRR2<0,0, P_SR$P_GWRR2)
    P_SR$P_GWRR3<-ifelse( P_SR$P_GWRR3<0,0, P_SR$P_GWRR3)
    P_SR$P_GWRR4<-ifelse( P_SR$P_GWRR4<0,0, P_SR$P_GWRR4)
    # P_SR<-cbind(P_SR,P_GWRR)
    coordinates(P_SR)<-~Xlon+Ylat
    Sat.data<-cbind(TRR.DS$yhat$Xlon,TRR.DS$yhat$Ylat,TRR.DS$yhat$p,CMO.DS$yhat$p,PER.DS$yhat$p,GSM.DS$yhat$p)
    Sat.data<-as.data.frame(Sat.data)
    colnames(Sat.data)<-c("Xlon","Ylat","P_TRMM.RT","P_CMO","P_PER","P_GSM")
    Sat.merge<-left_join(Sat.merge,Sat.data,by=c("Xlon","Ylat"))
    Sat.merge$P_TRMM.RT<-((Sat.merge$P_TRMM.RT)^0.25-1)/0.25
    Sat.merge$P_CMO<-((Sat.merge$P_CMO)^0.25-1)/0.25
    Sat.merge$P_PER<-((Sat.merge$P_PER)^0.25-1)/0.25
    Sat.merge$P_GSM<-((Sat.merge$P_GSM)^0.25-1)/0.25
    
    # dis1<-matrix(0,nrow(Sat.merge),1)
    # dis1.min<-matrix(NA,nstations,1)  #存放最小值所在的距离
    # P_SR1<-matrix(NA,nstations,8)    #Satellite rainfall and raingauge data
    # for (j in 1:nstations){
    #   dis1[,1]<-sqrt((a$Xlon[j]-Sat.merge$Xlon)^2+(a$Ylat[j]-Sat.merge$Ylat)^2)
    #   dis1.min[j,1]<-which.min(dis1[,1])      #寻找最小值所在TRMM网格中的排序,该位置和TRMM降雨位置对应
    #   P_SR1[j,1]<-a$STCD[j]
    #   P_SR1[j,2]<-a$Xlon[j]
    #   P_SR1[j,3]<-a$Ylat[j]
    #   P_SR1[j,4]<-a$p[j]   #P_error[j,1]本应写成P_error[j,i],但2014年7月1日-2014年8月不少站点数据缺失，各时间节点雨量站数目不同
    #   P_SR1[j,5]<-Sat.merge$P_TRMM.RT[dis1.min[j,1]]   #GWRR结果
    #   P_SR1[j,6]<-Sat.merge$P_CMO[dis1.min[j,1]]
    #   P_SR1[j,7]<-Sat.merge$P_PER[dis1.min[j,1]]
    #   P_SR1[j,8]<-Sat.merge$P_GSM[dis1.min[j,1]]
    # }
    # coordinates(Sat.merge)<-~Xlon+Ylat
    # dp.dist <- gwrr.dist(dp.locat=coordinates(P_SR),rp.locat = coordinates(Sat.merge))    #运行速度过慢
    if (is(P_SR, "Spatial")) {
      P_SR <- as(P_SR, "data.frame")
    }
    if (is(Sat.merge, "Spatial")) {
      Sat.merge <- as(Sat.merge, "data.frame")
    }
    max(Sat.merge$P_TRMM.RT)
    max(Sat.merge$P_CMO)
    max(Sat.merge$P_PER)
    max(Sat.merge$P_GSM)
    P_SR.pre0 <- gwrr.predict(P_raingauge ~ P_TRMM.RT  + P_CMO + P_PER + P_GSM, P_SR.BC, Sat.merge, kernel = "bis", bw = bw.SRmerge[i-1], dp.dist)
    # yhat1.merge<-as.vector(P_SR.pre1$yhat)
    # yhat1.merge_VBC<-(0.25*yhat1.merge+1)^(1/0.25)
    # write.csv( yhat1.merge_VBC,paste(Filepath,"/Output data/Merged precipitation/Statistics/pre.csv",sep = ""))
    yhat.merge<-as.vector(P_SR.pre0$yhat)
    max(yhat.merge)
    cn.merge<-P_SR.pre0$cn
    lambda.merge<-P_SR.pre0$lambda
    beta.merge<-t(P_SR.pre0$beta)
    P_SR.pre<-as.data.frame(cbind(Sat.merge$Xlon,Sat.merge$Ylat,yhat.merge,beta.merge,cn.merge,lambda.merge))
    colnames(P_SR.pre)<-c("Xlon","Ylat","yhat","beta0","beta_TRR","beta_CMO","beta_PER","beta_GSM","cn.merge","lamda.merge")
    P_SR.merge<-left_join(Sat.cor_1km,P_SR.pre,by=c("Xlon","Ylat"))
    beta0.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$beta0))
    beta_TRR.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$beta_TRR))
    beta_CMO.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$beta_CMO))
    beta_PER.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$beta_PER))
    beta_GSM.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$beta_GSM))
    write.asciigrid(beta0.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/beta0",dates[i],".txt",sep = ""))
    write.asciigrid(beta_TRR.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/beta_TRR",dates[i],".txt",sep = ""))
    write.asciigrid(beta_CMO.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/beta_CMO",dates[i],".txt",sep = ""))
    write.asciigrid(beta_PER.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/beta_PER",dates[i],".txt",sep = ""))
    write.asciigrid(beta_GSM.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/beta_GSM",dates[i],".txt",sep = ""))
    cn_merge.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$cn.merge))
    lambda_merge.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$lamda.merge))
    write.asciigrid(cn_merge.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/cn_merge",dates[i],".txt",sep = ""))
    write.asciigrid(lambda_merge.grid,paste(Filepath,"/Output data/Merged precipitation/Statistics/日参数/lambda_merge",dates[i],".txt",sep = ""))
    P_SR.merge$yhat<-(0.25*P_SR.merge$yhat+1)^(1/0.25)
    P_merge.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_SR.merge$yhat))
    is.vector(P_merge.grid@data[[1]])
    P_merge.raster<- raster(P_merge.grid,layer=1,values=T)
    P_merge.raster1<-P_merge.raster
# plot(P_merge.raster1)
# a1<-rasterToPoints(P_merge.raster1)
    # P_GWRR<-vector(mode = "numeric",length = nstations)
    # for (j in 1:nstations){
    #   P_GWRR[j]<-P_SR.merge$yhat[dis.min[j,1]]
    # }
    # max(P_SR.merge$yhat)
    # P_SR<-cbind(P_SR,P_GWRR)
    
    # fileCon <- paste(Filepath,"/Output data/Merged precipitation/Rainfall (SR merged)/Rainfall_",Date,".txt",sep="")
    # writeLines('ncols\t40\nnrows\t23\nxllcorner\t102\nyllcorner\t21.25\ncellsize\t0.25\nNODATA_value\t-9999', fileCon)#ascii格式中坐标为左下角栅格单元的左下角坐标，即左下角栅格中心点坐标-0.5*栅格大小
    # write.table(P_SR.pre,fileCon,append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)#将数据框写成ASCII文件，是按行，读的第一个数据为左上角栅格（最小经度最大纬度

    #####################################################
    #########              残差处理          ############
    #####################################################
    res<-P_SR.BC$P_raingauge - P_SR.est$yhat
    den<-as.data.frame(res)
    jpeg(paste("DEN_M",Date,".jpeg",sep=""))
    plot_den<-ggplot(den,aes(x=res),alpha = 1/10)+geom_density()
    multiplot(plot_den)
    dev.off()
    
    Res<-cbind(locs,res)
    Res<-as.data.frame(Res)
    names(Res)<-c("Xlon","Ylat","res")
    coordinates(Res)<-~Xlon+Ylat
    v<-variogram(Res$res~1,Res)
    v.fit<-fit.variogram(v,vgm("Sph"))
    if (v.fit[2,]$psill<=0) {
      v.fit<-vgm(mean(order(v$gamma,decreasing = FALSE)[1:5]),"Sph",max(v$dist),mean(order(v$gamma,decreasing = TRUE)))
    }
    if (v.fit[2,]$range<0){
      v.fit[2,]$range<-1
    }
    #dev.new()
    jpeg(paste("VGM_M",Date,".jpeg",sep=""))
    plot_vgm<-plot(v,model=v.fit)
    plot_vgm
    dev.off()
   #  jpeg(paste("M",Date,".jpeg",sep=""))
   # # par(mfrow=c(1,2))
   #  split.screen(c(1,2)) 
   #  screen(1)
   #  plot(v,v.fit)
   #  title("变异函数")
   #  screen(2)
   #  title("概率密度分布图")
   #  dev.off()
   #  close.screen(all.screens = T)
   #  erase.screen()
   #  
   #  jpeg(paste("M",Date,".jpeg",sep=""))
   #  par(mfrow=c(2,2))
   #  plot.new()
   #  plot(v,v.fit,main="变异函数")
   #  plot.new()
   #  plot(v,v.fit,main="变异函数")
   #  plot.new()
   #  plot(v,v.fit,main="变异函数")
   #  plot.new()
   #  plot(v,v.fit,main="变异函数")
   #  dev.off()
   #  plot(v,add=T)
   #  plot(v.fit)
   #  
   #  jpeg(paste("M",Date,".jpeg",sep=""))
   #  par(mfrow=c(2,2))
   #  plot(v)
   #  plot(v)
   #  plot(v)
   #  plot(v)
   #  dev.off()
   #  plot.vgram(v.fit)
   #  plot(v,model=v.fit)
   #  
   # 
   #  jpeg(paste("M",Date,".jpeg",sep=""))
   #  par(mfrow=c(2,2))
   #  x<-seq(-pi,pi,by=0.01)
   #  plot(x,sin(x))
   #  plot(x,cos(x))
   #  plot(x,2*sin(x)*cos(x))
   #  plot(x,tan(x))
   #  dev.off()
   #  
   #  ggplot(v,aes(x=dist,y=gamma))+geom_abline(v.fit)
   #  plot(v.fit)
   #  
   #  jpeg(paste("M",Date,".jpeg",sep=""))
   #  par(mfrow=c(2,2))
   #  ggplot(den,aes(x=res),alpha = 1/10)+geom_density()
   #  ggplot(den,aes(x=res),alpha = 1/10)+geom_density()
   #  ggplot(den,aes(x=res),alpha = 1/10)+geom_density()
   #  # plot(v,v.fit,main="变异函数",xlab="滞后距",ylab="半方差")
   #  plot(v,v.fit,main="变异函数")
   #  #title("变异函数")
   #  ggplot(den,aes(x=res),alpha = 1/10)+geom_density()+geom_histogram(fill="navy") #alpha设置透明度
   #  # p<-ggplot(as.data.frame(res),aes(x=res,y=..density..))
   #  # p<-p+geom_density(colour="green")
   #  # p<-p+geom_histogram(fill="navy")
   # 
   #  
   #  ggplot(den,aes(x=res),alpha = 1)+geom_histogram()
   #  ggplot(count,aes(x=res),alpha = 1)
   #  ggplot(as.data.frame(res))+geom_density(aes(x=res))+geom_histogram(aes(x=res))
   #  ggplot(as.data.frame(res))+geom_histogram(aes(x=res))
   #  ggplot(as.data.frame(res))+geom_histogram(aes(x=res),stat="identity")
    coordinates(Sat.merge)<-~Xlon+Ylat
    res_ok<-krige(res~1,Res,Sat.merge,model=v.fit)
    # max(InvBoxCox(res_ok$var1.pred,lam))
    # yhat.mergek<-yhat.merge+InvBoxCox(res_ok$var1.pred,lam)
    yhat.mergek<-yhat.merge+res_ok$var1.pred
    P_Merge<-cbind(Sat.base1_1km,yhat.merge,res_ok$var1.pred,yhat.mergek)
    colnames(P_Merge)<-c("Xlon","Ylat","isin","yhat.GWRR","res","yhat.GWRRK")
    ######box-cox逆变换
    P_Merge$yhat.GWRR<-(0.25*P_Merge$yhat.GWRR+1)^(1/0.25)
    P_Merge$yhat.GWRRK<-(0.25*P_Merge$yhat.GWRRK+1)^(1/0.25)
    P_Merge$res<-P_Merge$yhat.GWRRK-P_Merge$yhat.GWRR
    P_Merge<-left_join(Sat.cor_1km,P_Merge,by=c("Xlon","Ylat"))
    P_Merge$yhat.GWRR<-ifelse(P_Merge$yhat.GWRR<0.01,0,P_Merge$yhat.GWRR)
    P_Merge$yhat.GWRRK<-ifelse(P_Merge$yhat.GWRRK<0.01,0,P_Merge$yhat.GWRRK)
    #####前2%的最大值用第0.2%最大替代(相当于500个0.25度栅格去除一个最大异常，即整个流域去除一个0.25度异常)
    n<-floor(NROW(P_Merge$yhat.GWRR[!is.na(P_Merge$yhat.GWRR)])*0.002)
    max(P_Merge$yhat.GWRR[!is.na(P_Merge$yhat.GWRR)])
    max(P_Merge$yhat.GWRRK[!is.na(P_Merge$yhat.GWRRK)])
    yhat.GWRR_threshold<-P_Merge$yhat.GWRR[!is.na(P_Merge$yhat.GWRR)][order(P_Merge$yhat.GWRR[!is.na(P_Merge$yhat.GWRR)],decreasing=TRUE)[1:n]][n]
    yhat.GWRRK_threshold<-P_Merge$yhat.GWRRK[!is.na(P_Merge$yhat.GWRR)][order(P_Merge$yhat.GWRRK[!is.na(P_Merge$yhat.GWRR)],decreasing=TRUE)[1:n]][n]
    P_Merge$yhat.GWRR<-ifelse(P_Merge$yhat.GWRR > yhat.GWRR_threshold, yhat.GWRR_threshold,P_Merge$yhat.GWRR)
    P_Merge$yhat.GWRRK<-ifelse(P_Merge$yhat.GWRRK > yhat.GWRRK_threshold, yhat.GWRRK_threshold,P_Merge$yhat.GWRRK)
    
    P_merge.grid1 <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_Merge$yhat.GWRR))
    P_merge.raster1<- raster(P_merge.grid1,layer=1,values=T)
    P_merge.grid <- SpatialGridDataFrame(grid = GridTopology1km, data = as.data.frame(P_Merge$yhat.GWRRK))
    P_merge.raster<- raster(P_merge.grid,layer=1,values=T)
    #######################################
    #                动态地图             #    
    #######################################
   # help(package="leaflet")
   # crs(P_merge.raster)<-CRS("+proj=longlat +datum=WGS84")
   # pal<-colorNumeric(c("transparent",topo.colors(100,alpha = NULL)),values(P_merge.raster),na.color = "transparent")
   # leaflet() %>% addTiles() %>%
   # addRasterImage(x=P_merge.raster, colors=pal, opacity = 0.8, project = T) %>%
   # addLegend(pal=pal,values=values(P_merge.raster),title = "merge")
   # 

   
   # jpeg(paste(Date,".jpeg",sep=""))
   # split.screen(c(1,2)) 
   #  screen(1)
   #  plot(raster(TRR.DS$yhat.grid,layer=1,values=T))
   #  title("TRR.DS")
   #  screen(2)
   #  plot(raster(CMO.DS$yhat.grid,layer=1,values=T))
   #  title("CMO.DS")
   #  # screen(3)
   #  # plot(PER.DS$yhat.grid)
   #  # title("PER.DS")
   #  # screen(4)
   #  # plot(TRR.DS$yhat.grid)
   #  # title("PER.DS")
   #  # screen(5)
   #  # plot(P_merge.raster1)
   #  # title("GWRR-based merged precipitation")
   #  # screen(6)
   #  # plot(P_merge.raster)
   #  # title("GWRRK-based merged precipitation")
   #  dev.off()
   #  close.screen(all.screens = T)
   #  erase.screen()
   #  
   #  jpeg(paste(Date,".jpeg",sep=""))
   #  par(mfrow=c(3,2))
   #  plot(TRR.DS$yhat.grid)
   #  title("TRR.DS")
   #  plot(CMO.DS$yhat.grid)
   #  title("CMO.DS")
   #  plot(PER.DS$yhat.grid)
   #  title("PER.DS")
   #  plot(TRR.DS$yhat.grid)
   #  title("PER.DS")
   #  plot(P_merge.raster1)
   #  title("GWRR-based merged precipitation")
   #  plot(P_merge.raster)
   #  title("GWRRK-based merged precipitation")
   #  dev.off()
   #  
   #######################################
   #              输出一张图             #    
   #######################################
    jpeg(paste(Date,".jpeg",sep=""))
    par(mfrow=c(3,2))
    plot(raster(TRR.DS$yhat.grid,layer=1,values=T))
    title("TRR.DS")
    plot(raster(CMO.DS$yhat.grid,layer=1,values=T))
    title("CMO.DS")
    plot(raster(PER.DS$yhat.grid,layer=1,values=T))
    title("PER.DS")
    plot(raster(GSM.DS$yhat.grid,layer=1,values=T))
    title("GSM.DS")
    plot(P_merge.raster1)
    title("GWRR-based merged precipitation")
    plot(P_merge.raster)
    title("GWRRK-based merged precipitation")
    dev.off()
    
    fileCon <- paste(Filepath,"/Output data/Merged precipitation/Rainfall (SR merged)/Rainfall_",Date,".txt",sep="")
    write.asciigrid(P_merge.grid,fileCon)
    write.asciigrid(P_merge.grid1,paste(Filepath,"/Output data/Merged precipitation/Rainfall (SR merged_GWR)/Rainfall_",Date,".txt",sep=""))
    g<-vector(mode="numeric",length=nstations)
    P_GWRRK<-vector(mode="numeric",length=nstations)
    g1<-vector(mode="numeric",length=nstations)
    dis1<-matrix(0,nrow(Sat.merge),1)
    dis1.min<-matrix(NA,nstations,1)  #存放最小值所在的距离
    P_SR1<-matrix(NA,nstations,8)    #Satellite rainfall and raingauge data
    for (j in 1:nstations){
      P_GWRRK[j]<-P_Merge$yhat.GWRRK[dis.min[j,1]]
      g1[j]<-P_Merge$yhat.GWRR[dis.min[j,1]]
      dis1[,1]<-sqrt((a$Xlon[j]-Sat.merge$Xlon)^2+(a$Ylat[j]-Sat.merge$Ylat)^2)
      dis1.min[j,1]<-which.min(dis1[,1])      #寻找最小值所在TRMM网格中的排序,该位置和TRMM降雨位置对应
      g[j]<-P_SR.pre0$yhat[dis1.min[j,1]]
    }
    g<-(0.25*g+1)^(1/0.25)
    P_SR<-cbind(P_SR,P_GWRRK,g1)
    P_SR$P_GWRR<-ifelse(P_SR$P_GWRR<0,0,P_SR$P_GWRR)     ###克里金插值之前不强制为0，插值之后强制为0，尽量保证残差正态
    ########################各卫星产品及融合产品精度评价结果###########################
    Stas.TRR[i-1,]<-statistic(P_SR,P_SR$P_raingauge,P_SR$P_TRMM.RT)
    Stas.CMO[i-1,]<-statistic(P_SR,P_SR$P_raingauge,P_SR$P_CMO)
    Stas.PER[i-1,]<-statistic(P_SR,P_SR$P_raingauge,P_SR$P_PER)
    Stas.GSM[i-1,]<-statistic(P_SR,P_SR$P_raingauge,P_SR$P_GSM)
    Stas.GWRR[i-1,]<-statistic(P_SR,P_SR$P_raingauge,P_SR$P_GWRR)
    Stas.GWRRK[i-1,]<-statistic(P_SR,P_SR$P_raingauge,P_SR$P_GWRRK)
    #######################卫星产品及融合产品站点雨量对比（每天及总表）#######################
    path_Station.P<-paste(Filepath,"/Output data/Merged precipitation/数据对比/Station_rainfall",Date,".xlsx",sep = "")
    write.xlsx(P_SR,path_Station.P)
    if (i==2) {
      AllP_SR<-P_SR
    }
    else {AllP_SR<-rbind(AllP_SR,P_SR)}
    # Sat.top<-Sat.top[,-3]   #
    # GSM.top<-GSM.top[,-3]
    cat("k =", i,"is done\tElapsedTime =", ElapsedTime(time2),"\n")
    bw.SRmerge[i-1]
 }

# stopCluster(cl)
path_Station.AllP<-paste(Filepath,"/Output data/Merged precipitation/Statistics/Station_rainfall总表.csv",sep = "")
write.csv(AllP_SR,path_Station.AllP)
write.csv(Allyhat_25km.TRR,paste(Filepath,"/Output data/Merged precipitation/Statistics/TRR降尺度总表.csv",sep = ""))
write.csv(Allyhat_25km.CMO,paste(Filepath,"/Output data/Merged precipitation/Statistics/CMO降尺度总表.csv",sep = ""))
write.csv(Allyhat_25km.PER,paste(Filepath,"/Output data/Merged precipitation/Statistics/PER降尺度总表.csv",sep = ""))
write.csv(Allyhat_10km.GSM,paste(Filepath,"/Output data/Merged precipitation/Statistics/GSM降尺度总表.csv",sep = ""))
ls()    #显示变量

###################################
#         精度评价结果汇总        #
###################################
write.csv(cn,paste(Filepath,"/Output data/Merged precipitation/Statistics/多卫星融合与单卫星融合/cn.csv",sep = ""))
write.csv(lambda,paste(Filepath,"/Output data/Merged precipitation/Statistics/多卫星融合与单卫星融合/lambda.csv",sep = ""))
write.csv(RMSE,paste(Filepath,"/Output data/Merged precipitation/Statistics/多卫星融合与单卫星融合/RMSE.csv",sep = ""))
write.csv(MAE,paste(Filepath,"/Output data/Merged precipitation/Statistics/多卫星融合与单卫星融合/MAE.csv",sep = ""))
write.csv(rsquare,paste(Filepath,"/Output data/Merged precipitation/Statistics/多卫星融合与单卫星融合/R2.csv",sep = ""))
write.csv(CVscore,paste(Filepath,"/Output data/Merged precipitation/Statistics/多卫星融合与单卫星融合/CVscore.csv",sep = ""))
#########################           统计总表---卫星总评（总时间序列）          ###########################
para<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
trr<-statistic(AllP_SR,AllP_SR$P_raingauge,AllP_SR$P_TRMM.RT)
cmo<-statistic(AllP_SR,AllP_SR$P_raingauge,AllP_SR$P_CMO)
per<-statistic(AllP_SR,AllP_SR$P_raingauge,AllP_SR$P_PER)
gsm<-statistic(AllP_SR,AllP_SR$P_raingauge,AllP_SR$P_GSM)
gwrr<-statistic(AllP_SR,AllP_SR$P_raingauge,AllP_SR$P_GWRR)
gwrrk<-statistic(AllP_SR,AllP_SR$P_raingauge,AllP_SR$P_GWRRK)
Tot.Sta<-cbind(para,trr,cmo,per,gsm,gwrr,gwrrk)
names(Tot.Sta)<-c("Para","TRR","CMO","PER","GSM","GWRR","GWRRK")
path_TOt.Sta<- paste(Filepath,"/Output data/Merged precipitation/Statistics/总表.csv",sep = "")
write.csv(Tot.Sta,path_TOt.Sta)

##########################      各参数分站点统计表--空间分析（总时间序列）     ###########################
for (k in 1:Nstations){
  c<-filter(AllP_SR,STCD==raingauges$STCD[k])
  trr.sta<-statistic(c,c$P_raingauge,c$P_TRMM.RT)
  cmo.sta<-statistic(c,c$P_raingauge,c$P_CMO)
  per.sta<-statistic(c,c$P_raingauge,c$P_PER)
  gsm.sta<-statistic(c,c$P_raingauge,c$P_GSM)
  gwrr.sta<-statistic(c,c$P_raingauge,c$P_GWRR)
  gwrrk.sta<-statistic(c,c$P_raingauge,c$P_GWRRK)
  rmse.sta<-c(c$STCD[1],trr.sta[1],cmo.sta[1],per.sta[1],gsm.sta[1],gwrr.sta[1],gwrrk.sta[1])
  mae.sta<-c(c$STCD[1],trr.sta[2],cmo.sta[2],per.sta[2],gsm.sta[2],gwrr.sta[2],gwrrk.sta[2])
  bias.sta<-c(c$STCD[1],trr.sta[3],cmo.sta[3],per.sta[3],gsm.sta[3],gwrr.sta[3],gwrrk.sta[3])
  cc.sta<-c(c$STCD[1],trr.sta[4],cmo.sta[4],per.sta[4],gsm.sta[4],gwrr.sta[4],gwrrk.sta[4])
  me.sta<-c(c$STCD[1],trr.sta[5],cmo.sta[5],per.sta[5],gsm.sta[5],gwrr.sta[5],gwrrk.sta[5])
  pod.sta<-c(c$STCD[1],trr.sta[6],cmo.sta[6],per.sta[6],gsm.sta[6],gwrr.sta[6],gwrrk.sta[6])
  pofd.sta<-c(c$STCD[1],trr.sta[7],cmo.sta[7],per.sta[7],gsm.sta[7],gwrr.sta[7],gwrrk.sta[7])
  ets.sta<-c(c$STCD[1],trr.sta[8],cmo.sta[8],per.sta[8],gsm.sta[8],gwrr.sta[8],gwrrk.sta[8])
  if (k==1) {
   rmse.stas<-rmse.sta
   mae.stas<-mae.sta
   bias.stas<-bias.sta
   cc.stas<-cc.sta
   me.stas<-me.sta
   pod.stas<-pod.sta
   pofd.stas<-pofd.sta
   ets.stas<-ets.sta
  }
  else {
    rmse.stas<-rbind(rmse.stas,rmse.sta)
    mae.stas<-rbind(mae.stas,mae.sta)
    bias.stas<-rbind(bias.stas,bias.sta)
    cc.stas<-rbind(cc.stas,cc.sta)
    me.stas<-rbind(me.stas,me.sta)
    pod.stas<-rbind(pod.stas,pod.sta)
    pofd.stas<-rbind(pofd.stas,pofd.sta)
    ets.stas<-rbind(ets.stas,ets.sta)
  }
}
row.names(rmse.stas)<-rmse.stas[,1]
rmse.stas<-rmse.stas[,-1]
row.names(mae.stas)<-mae.stas[,1]
mae.stas<-mae.stas[,-1]
row.names(bias.stas)<-bias.stas[,1]
bias.stas<-bias.stas[,-1]
row.names(cc.stas)<-cc.stas[,1]
cc.stas<-cc.stas[,-1]
row.names(me.stas)<-me.stas[,1]
me.stas<-me.stas[,-1]
row.names(pod.stas)<-pod.stas[,1]
pod.stas<-pod.stas[,-1]
row.names(pofd.stas)<-pofd.stas[,1]
pofd.stas<-pofd.stas[,-1]
row.names(ets.stas)<-ets.stas[,1]
ets.stas<-ets.stas[,-1]
rmse.stas<-as.data.frame(rmse.stas)
mae.stas<-as.data.frame(mae.stas)
bias.stas<-as.data.frame(bias.stas)
cc.stas<-as.data.frame(cc.stas)
me.stas<-as.data.frame(me.stas)
pod.stas<-as.data.frame(pod.stas)
pofd.stas<-as.data.frame(pofd.stas)
ets.stas<-as.data.frame(ets.stas)
names(rmse.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(mae.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(bias.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(cc.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(me.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(pod.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(pofd.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
names(ets.stas)<-c("TRR","CMO","PER","GSM","GWRR","GWRRK")
path_rmse.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/rmse_stas.csv",sep = "")
path_mae.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/mae_stas.csv",sep = "")
path_bias.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/bias_stas.csv",sep = "")
path_cc.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/cc_stas.csv",sep = "")
path_me.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/me_stas.csv",sep = "")
path_pod.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/pod_stas.csv",sep = "")
path_pofd.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/pofd_stas.csv",sep = "")
path_ets.stas<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分站点统计/ets_stas.csv",sep = "")
write.csv(rmse.stas,path_rmse.stas)
write.csv(mae.stas,path_mae.stas)
write.csv(bias.stas,path_bias.stas)
write.csv(cc.stas,path_cc.stas)
write.csv(me.stas,path_me.stas)
write.csv(pod.stas,path_pod.stas)
write.csv(pofd.stas,path_pofd.stas)
write.csv(ets.stas,path_ets.stas)

#########################         分卫星统计（每天数据）    ##############################
colnames(Stas.TRR)<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
row.names(Stas.TRR)<-c(as.character(dates[2:length(dates)]))
colnames(Stas.CMO)<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
row.names(Stas.CMO)<-c(as.character(dates[2:length(dates)]))
colnames(Stas.PER)<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
row.names(Stas.PER)<-c(as.character(dates[2:length(dates)]))
colnames(Stas.GSM)<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
row.names(Stas.GSM)<-c(as.character(dates[2:length(dates)]))
colnames(Stas.GWRR)<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
row.names(Stas.GWRR)<-c(as.character(dates[2:length(dates)]))
colnames(Stas.GWRRK)<-c("RMSE","MAE","BIAS","CC","ME","POD","POFD","ETS")
row.names(Stas.GWRRK)<-c(as.character(dates[2:length(dates)]))
write.csv(Stas.TRR,paste(Filepath,"/Output data/Merged precipitation/Statistics/各卫星参数统计/TRMM_RT Statistics.csv",sep=""),row.names=T,col.names=T)
write.csv(Stas.CMO,paste(Filepath,"/Output data/Merged precipitation/Statistics/各卫星参数统计/CMORPH Statistics.csv",sep=""),row.names=T,col.names=T)
write.csv(Stas.PER,paste(Filepath,"/Output data/Merged precipitation/Statistics/各卫星参数统计/PERSIANN Statistics.csv",sep=""),row.names=T,col.names=T)
write.csv(Stas.GSM,paste(Filepath,"/Output data/Merged precipitation/Statistics/各卫星参数统计/GSMaP Statistics.csv",sep=""),row.names=T,col.names=T)
write.csv(Stas.GWRR,paste(Filepath,"/Output data/Merged precipitation/Statistics/各卫星参数统计/GWRR Statistics.csv",sep=""),row.names=T,col.names=T)
write.csv(Stas.GWRRK,paste(Filepath,"/Output data/Merged precipitation/Statistics/各卫星参数统计/GWRRK Statistics.csv",sep=""),row.names=T,col.names=T)
###########################        分参数统计    #############################
RMSE.dat<-cbind(Stas.TRR[,1],Stas.CMO[,1],Stas.PER[,1],Stas.CMO[,1],Stas.GWRR[,1],Stas.GWRRK[,1])
MAE.dat<-cbind(Stas.TRR[,2],Stas.CMO[,2],Stas.PER[,2],Stas.CMO[,2],Stas.GWRR[,2],Stas.GWRRK[,2])
BIAS.dat<-cbind(Stas.TRR[,3],Stas.CMO[,3],Stas.PER[,3],Stas.CMO[,3],Stas.GWRR[,3],Stas.GWRRK[,3])
CC.dat<-cbind(Stas.TRR[,4],Stas.CMO[,4],Stas.PER[,4],Stas.CMO[,4],Stas.GWRR[,4],Stas.GWRRK[,4])
ME.dat<-cbind(Stas.TRR[,5],Stas.CMO[,5],Stas.PER[,5],Stas.CMO[,5],Stas.GWRR[,5],Stas.GWRRK[,5])
POD.dat<-cbind(Stas.TRR[,6],Stas.CMO[,6],Stas.PER[,6],Stas.CMO[,6],Stas.GWRR[,6],Stas.GWRRK[,6])
POFD.dat<-cbind(Stas.TRR[,7],Stas.CMO[,7],Stas.PER[,7],Stas.CMO[,7],Stas.GWRR[,7],Stas.GWRRK[,7])
ETS.dat<-cbind(Stas.TRR[,8],Stas.CMO[,8],Stas.PER[,8],Stas.CMO[,8],Stas.GWRR[,8],Stas.GWRRK[,8])
row.names(RMSE.dat)<-as.character(dates[2:length(dates)])
row.names(MAE.dat)<-as.character(dates[2:length(dates)])
row.names(RMSE.dat)<-as.character(dates[2:length(dates)])
row.names(CC.dat)<-as.character(dates[2:length(dates)])
row.names(ME.dat)<-as.character(dates[2:length(dates)])
row.names(POD.dat)<-as.character(dates[2:length(dates)])
row.names(POFD.dat)<-as.character(dates[2:length(dates)])
row.names(ETS.dat)<-as.character(dates[2:length(dates)])
colnames(RMSE.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(MAE.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(BIAS.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(CC.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(ME.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(POD.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(POFD.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
colnames(ETS.dat)<-c("TRMM.RT","CMO","PER","GSM","GWRR","GWRRK")
# path_RMSE<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/RMSE.xlsx",sep = "")
# path_MAE<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/MAE.xlsx",sep = "")
# path_BIAS<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/BIAS.xlsx",sep = "")
# path_CC<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/CC.xlsx",sep = "")
# path_ME<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/ME.xlsx",sep = "")
# path_POD<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/POD.xlsx",sep = "")
# path_POFD<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/POFD.xlsx",sep = "")
# path_ETS<-paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/ETS.xlsx",sep = "")
# write.xlsx(RMSE.dat,path_RMSE)
# write.xlsx(MAE.dat,path_MAE)
# write.xlsx(BIAS.dat,path_BIAS)
# write.xlsx(CC.dat,path_CC)
# write.xlsx(ME.dat,path_ME)
# write.xlsx(POD.dat,path_POD)
# write.xlsx(POFD.dat,path_POFD)
# write.xlsx(ETS.dat,path_ETS)
write.csv(RMSE.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/RMSE.csv",sep = ""))
write.csv(MAE.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/MAE.csv",sep = ""))
write.csv(BIAS.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/BIAS.csv",sep = ""))
write.csv(CC.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/CC.csv",sep = ""))
write.csv(ME.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/ME.csv",sep = ""))
write.csv(POD.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/POD.csv",sep = ""))
write.csv(POFD.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/POFD.csv",sep = ""))
write.csv(ETS.dat,paste(Filepath,"/Output data/Merged precipitation/Statistics/各参数分天统计/ETS.csv",sep = ""))
path_cn.TRR<-paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/TRR cn.xlsx",sep = "")
path_cn.CMO<-paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/CMO cn.xlsx",sep = "")
path_cn.PER<-paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/PER cn.xlsx",sep = "")
path_cn.GSM<-paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/GSM cn.xlsx",sep = "")
path_cn.est<-paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/EST cn.xlsx",sep = "")
path_lambda.TRR<-paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/TRR lambda.xlsx",sep = "")
path_lambda.CMO<-paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/CMO lambda.xlsx",sep = "")
path_lambda.PER<-paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/PER lambda.xlsx",sep = "")
path_lambda.GSM<-paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/GSM lambda.xlsx",sep = "")
path_lambda.est<-paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/EST lambda.xlsx",sep = "")
# write.xlsx(cn.TRR,path_cn.TRR)
# write.xlsx(cn.CMO,path_cn.CMO)
# write.xlsx(cn.PER,path_cn.PER)
# write.xlsx(cn.GSM,path_cn.GSM)
# write.xlsx(cn,path_cn.est)
write.csv(cn.TRR,paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/TRR cn.csv",sep = ""))
write.csv(cn.CMO,paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/CMO cn.csv",sep = ""))
write.csv(cn.PER,paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/PER cn.csv",sep = ""))
write.csv(cn.GSM,paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/GSM cn.csv",sep = ""))
write.csv(cn,paste(Filepath,"/Output data/Merged precipitation/Statistics/CN/EST cn.csv",sep = ""))

write.csv(lambda.TRR,paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/TRR lambda.csv",sep = ""))
write.csv(lambda.CMO,paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/CMO lambda.csv",sep = ""))
write.csv(lambda.PER,paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/PER lambda.csv",sep = ""))
write.csv(lambda.GSM,paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/GSM lambda.csv",sep = ""))
write.csv(lambda,paste(Filepath,"/Output data/Merged precipitation/Statistics/LAMBDA/EST lambda.csv",sep = ""))
# write.xlsx(lambda.TRR,path_lambda.TRR)
# write.xlsx(lambda.CMO,path_lambda.CMO)
# write.xlsx(lambda.PER,path_lambda.PER)
# write.xlsx(lambda.GSM,path_lambda.GSM)
# write.xlsx(lambda,path_lambda.est)
###################################
#               GWR               #
###################################
#GWR假设回归参数唯一连续表面（空间平滑），位置相邻的回归参数相似，以实测点i及其若干最近邻观测值为样本，建立局域回归模型
#GWR的综合考虑了因变量的局部自相关性（权重因子，通过空间权函数计算）和因变量与其他相关变量之间的互相关性（回归方程）
#GWR误差项正态独立假设
# coordinates(rainfall)=~Xcor+Ycor
# rain.dist<-gw.dist(dp.locat=coord_obj[,c(2,3)])    #生成距离矩阵
# rain.b<-bw.gwr(P ~ VX+VY , data=rain.m,approach="CV",adaptive=TRUE, kernel = "bisquare", dMat= rain.dist)  #返回带宽。adaptive=TRUE，则返回自适应带宽，否则返回固定带宽。确定最优带宽的方法有（广义）交叉验证法（CV)和AIC准则
# rain.gwr<-gwr.basic(P ~ VX+VY , data=rain.m, bw=rain.b,adaptive=TRUE, kernel = "bisquare", dMat=rain.dist)
# rain.pre<-gwr.predict(P ~ VX+VY , data=rain.m, predictdata=rain.P, bw=9, kernel="bisquare", adaptive=T,  dMat2=rain.dist)    #空间预测

ElapsedTime(time1)