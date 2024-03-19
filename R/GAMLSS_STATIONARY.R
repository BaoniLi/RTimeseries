library(maxLik)
library(gamlss)
library(Kendall)
library(moments)
library(PearsonDS)

#change the work direction##########################################
setwd("E:/dissertation/3/GAMLSS")

data_file=list.files(pattern = "csv")
N=length(data_file)
SBC=array(dim=c(5,N))
p_KS=array(dim=c(5,N))
v1=array(dim=c(1,N))
v2=array(dim=c(1,N))
v3=array(dim=c(1,N))


for(i in 1:N)
{
  data<-read.csv(data_file[i],head=T)
  Y=data$Y
  RH=data$RH
  Q=data$Q
  EA=data$EA
  VPD=data$VPD
  Te=data$T
  P=data$P
  PET=data$PET
  n=length(RH)
  t=c(1:n)
  
  #change the variable########################################
  test_data=PET
  

  #stationary=====================================
  
  GU0<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=GU(mu.link = "identity",sigma.link = "identity"))
  WEI0<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=WEI(mu.link = "identity",sigma.link = "log"))
  GA0<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=GA(mu.link = "identity",sigma.link = "identity"))
  N0<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
  LN0<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
  
  
  SBC[1,i]=GU0$sbc
  SBC[2,i]=WEI0$sbc
  SBC[3,i]=GA0$sbc
  SBC[4,i]=N0$sbc
  SBC[5,i]=LN0$sbc
  
  p_KS[1,i]=ks.test(pnorm(GU0$residuals),"punif",0,1)$p.value
  p_KS[2,i]=ks.test(pnorm(WEI0$residuals),"punif",0,1)$p.value
  p_KS[3,i]=ks.test(pnorm(GA0$residuals),"punif",0,1)$p.value
  p_KS[4,i]=ks.test(pnorm(N0$residuals),"punif",0,1)$p.value
  p_KS[5,i]=ks.test(pnorm(LN0$residuals),"punif",0,1)$p.value
  
  mn=which(SBC[,i]==min(SBC[,i]))
  
  if(mn==1)
  {
    print("GU")
    print(GU0)
  }
  if(mn==2)
  {
    print("WEI")
    print(WEI0)
  }
  if(mn==3)
  {
    print("GA")
    print(GA0)
  }
  if(mn==4)
  {
    print("NO")
    print(N0)
  }
  if(mn==5)
  {
    print("LNO")
    print(LN0)
  }
  
  print(p_KS[mn,i])
  
  print("*************Following is the next station***********************")
  
  v1[1,i]=SBC[mn,i]
  v2[1,i]=p_KS[mn,i]
  v3[1,i]=mn
}
result=rbind(v1,v2,v3)
write.csv(result,"E:/dissertation/3/111/PETSTA.csv")