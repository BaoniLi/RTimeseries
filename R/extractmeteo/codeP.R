setwd("E:/data/全国气象数据19652018/气象数据19652018/PRE")

datfiles=list.files(pattern = "TXT")
N=length(datfiles)
n=julian(as.Date("2017-12-31"),origin=as.Date("1965-01-01"))[[1]]+1  #返回时间长度
ST=read.table("E:/YANGTZE/YZstation/urbancountryside.txt", head=T)
Num_Staion=length(ST$ID)

for(i in 1:16)
{
name=as.character(ST$ID[i])
Data_Frame=data.frame(date=as.Date(1:n, origin="1964-12-31"),P=array(dim=n,-999))
for(j in 1:N)
{
data <- array(scan(datfiles[j]), dim=c(13, length(scan(datfiles[j]))/13)) 
date=as.Date(paste(as.character(data[5,]),"-",as.character(data[6,]),"-",as.character(data[7,]),sep=""))
New_data=data.frame(ID=data[1,],Day=julian(date,origin=as.Date("1965-01-01"))[1:length(data[1,])]+1,P=data[10,])
stations=unique(New_data$ID)

NN=which(New_data$ID==ST$ID[i])
if(length(NN)!=0)
{
   
for (k in 1:length(NN)) {
  

  if(New_data$P[NN[k]]>10000)
  {
    
    if(New_data$P[NN[k]]==32700)
    {
      New_data$P[NN[k]]=0
    }
    
    else if(trunc(New_data$P[NN[k]]/1000)==32)
    {New_data$P[NN[k]]=New_data$P[NN[k]]%%1000}
    else if(trunc(New_data$P[NN[k]]/1000)==31)
    {New_data$P[NN[k]]=New_data$P[NN[k]]%%1000}
    else if(trunc(New_data$P[NN[k]]/1000)==30)
    {New_data$P[NN[k]]=New_data$P[NN[k]]%%1000}
  }
P_Sub=subset(New_data,New_data$ID==ST$ID[i])
Data_Frame$P[P_Sub$Day]=P_Sub$P
}
  print(j)
}
}
output=paste(name,".dat",sep="")
write.table(Data_Frame,output)

cat("The program is calculating station: ",ST$ID[i],"\n")
}