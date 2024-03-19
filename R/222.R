setwd("E:/PCA/PCA_input")

data_file=list.files(pattern = "csv")
N=length(data_file)

SBC=array(dim=c(5,N))
p_KS=array(dim=c(5,N))
v1=array(dim=c(1,N))
v2=array(dim=c(1,N))
v3=array(dim=c(1,N))

for(i in 1:1)
{
  data<-read.csv(data_file[i],head=T)
  Y=data$Y
  TEM=data$T
  IM=data$I
  LAI=data$L
  n=length(Y)
  t=c(1:n)
  
  FIRST <- lm(TEM~IM)
  SECOND <- lm(TEM~LAI)
  THIRD <- lm(TEM~IM+LAI)
  
}