setwd("D:\datasets\meteo0718\t")
datfiles=list.files(pattern = "TXT")
N=length(datfiles)
n=julian(as.Date("2018-12-31"),origin=as.Date("2007-01-01"))[[1]]+1  #返回时间长度
ST=read.table("D:\rice\meteosta\meteosta.txt", head=TRUE)
Num_Staion=length(ST$ID)

for(i in 1:Num_Staion)
{
name=as.character(ST$ID[i])
Data_Frame=data.frame(date=as.Date(1:n, origin="2006-12-31"),mean=array(dim=n,-999),max=array(dim=n,-999),min=array(dim=n,-999))
for(j in 1:N)
{
data <- array(scan(datfiles[j]), dim=c(13, length(scan(datfiles[j]))/13)) 
date=as.Date(paste(as.character(data[5,]),"-",as.character(data[6,]),"-",as.character(data[7,]),sep=""))
New_data=data.frame(ID=data[1,],Day=julian(date,origin=as.Date("2007-01-01"))[1:length(data[1,])]+1,mean=data[8,],max=data[9,],min=data[10,])
stations=unique(New_data$ID)

NN=which(New_data$ID==ST$ID[i])
if(length(NN)!=0)
{
  
  for (k in 1:length(NN)) {
    
mean_Sub=subset(New_data,New_data$ID==ST$ID[i])
max_Sub=subset(New_data,New_data$ID==ST$ID[i])
min_Sub=subset(New_data,New_data$ID==ST$ID[i])
Data_Frame$mean[mean_Sub$Day]=mean_Sub$mean
Data_Frame$max[max_Sub$Day]=max_Sub$max
Data_Frame$min[min_Sub$Day]=min_Sub$min
}
  print(j)
}
}
output=paste(name,".dat",sep="")
write.table(Data_Frame,output)

cat("The program is calculating station: ",ST$ID[i],"\n")
}
