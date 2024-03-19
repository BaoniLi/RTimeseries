library(maxLik)
library(gamlss)
library(Kendall)
library(moments)
library(PearsonDS)

#change the work direction##########################################
setwd("E:/dissertation/3/GAMLSS")
data_file=list.files(pattern = "csv")
N=length(data_file)
v1=array(dim=c(1,N))#########################################################
v2=array(dim=c(1,N))
v3=array(dim=c(1,N))
v4=array(dim=c(1,N))
musigma=c(1:53)

for(i in 1:N)
{
  data<-read.csv(data_file[i],head=T)
  
  Y=data$Y
  RH=data$RH
  Q=data$Q
  EA=data$EA
  VPD=data$VPD
  Te=data$T
  n=length(RH)
  t=c(1:n)
  
  ##change the variable#########################################
  test_data=Q
  #=========================================================================================================================
  type=1
  cp=0
  NO1<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
  SBC=NO1$sbc
  
  for(j in 10:(n-10))
  {
    tt=array(dim=n,1)
    tt[1:j]=0
    NO2<-gamlss(test_data~I(tt^1),sigma.formula=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
    if(SBC>NO2$sbc)
    {
      SBC=NO2$sbc
      type=2
      cp=j
    }
  }
  
  for(j in 10:(n-10))
  {
    tt=array(dim=n,1)
    tt[1:j]=0
    NO3<-gamlss(test_data~I(t^0),sigma.formula=~I(tt^1),family=NO(mu.link = "identity",sigma.link = "identity"))
    if(SBC>NO3$sbc)
    {
      SBC=NO3$sbc
      type=3
      cp=j
    }
  }
  
  
  for(j in 10:(n-10))
  {
    tt=array(dim=n,1)
    tt[1:j]=0
    NO4<-gamlss(test_data~I(tt^1),sigma.formula=~I(tt^1),family=NO(mu.link = "identity",sigma.link = "identity"))
    if(SBC>NO4$sbc)
    {
      SBC=NO4$sbc
      type=4
      cp=j
    }
  }
  
  
  #=========================================================================================================================
  
  LNO1<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
  if(SBC>LNO1$sbc)
  {
    SBC=LNO1$sbc
    type=5
  }
  
  for(j in 10:(n-10))
  {
    tt=array(dim=n,1)
    tt[1:j]=0
    LNO2<-gamlss(test_data~I(tt^1),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
    if(SBC>LNO2$sbc)
    {
      SBC=LNO2$sbc
      type=6
      cp=j
    }
  }
  
  for(j in 10:(n-10))
  {
    tt=array(dim=n,1)
    tt[1:j]=0
    LNO3<-gamlss(test_data~I(t^0),sigma.formula=~I(tt^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
    if(SBC>LNO3$sbc)
    {
      SBC=LNO3$sbc
      type=7
      cp=j
    }
  }
  
  
  for(j in 10:(n-10))
  {
    tt=array(dim=n,1)
    tt[1:j]=0
    LNO4<-gamlss(test_data~I(tt^1),sigma.formula=~I(tt^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
    if(SBC>LNO4$sbc)
    {
      SBC=LNO4$sbc
      type=8
      cp=j
    }
  }
  
  
  
  qQ=array(dim=c(n,6))
  #=====================================
  if(type==1)
  {
    print(NO1)
    p_KS=ks.test(pnorm(NO1$residuals),"punif",0,1)$p.value
    print(p_KS)
    mu=array(dim=n,NO1$mu.coefficients[1])
    sigma=array(dim=n,NO1$sigma.coefficients[1])
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  if(type==2)
  {
    tt=array(dim=n,1)
    tt[1:cp]=0
    NO2<-gamlss(test_data~I(tt^1),sigma.formula=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
    print(cp)
    print(NO2)
    p_KS=ks.test(pnorm(NO2$residuals),"punif",0,1)$p.value
    print(p_KS)
    
    
    mu=array(dim=n,NO2$mu.coefficients[1]+NO2$mu.coefficients[2]*tt)
    sigma=array(dim=n,NO2$sigma.coefficients[1])
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
    
  }
  
  if(type==3)
  {
    tt=array(dim=n,1)
    tt[1:cp]=0
    NO3<-gamlss(test_data~I(t^0),sigma.formula=~I(tt^1),family=NO(mu.link = "identity",sigma.link = "identity"))
    print(cp)
    print(NO3)
    p_KS=ks.test(pnorm(NO3$residuals),"punif",0,1)$p.value
    print(p_KS)
    
    mu=array(dim=n,NO3$mu.coefficients[1])
    sigma=array(dim=n,NO3$sigma.coefficients[1]+NO3$sigma.coefficients[2]*tt)
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
    
  }
  
  
  if(type==4)
  {
    tt=array(dim=n,1)
    tt[1:cp]=0
    NO4<-gamlss(test_data~I(tt^1),sigma.formula=~I(tt^1),family=NO(mu.link = "identity",sigma.link = "identity"))
    print(cp)
    print(NO4)
    p_KS=ks.test(pnorm(NO4$residuals),"punif",0,1)$p.value
    print(p_KS)
    
    mu=array(dim=n,NO4$mu.coefficients[1]+NO4$mu.coefficients[2]*tt)
    sigma=array(dim=n,NO4$sigma.coefficients[1]+NO4$sigma.coefficients[2]*tt)
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  
  if(type==5)
  {
    print(LNO1)
    p_KS=ks.test(pnorm(LNO1$residuals),"punif",0,1)$p.value
    print(p_KS)
    mu=array(dim=n,LNO1$mu.coefficients[1])
    sigma=array(dim=n,LNO1$sigma.coefficients[1])
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  if(type==6)
  {
    tt=array(dim=n,1)
    tt[1:cp]=0
    LNO2<-gamlss(test_data~I(tt^1),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
    print(cp)
    print(LNO2)
    p_KS=ks.test(pnorm(LNO2$residuals),"punif",0,1)$p.value
    print(p_KS)
    
    
    mu=array(dim=n,LNO2$mu.coefficients[1]+LNO2$mu.coefficients[2]*tt)
    sigma=array(dim=n,LNO2$sigma.coefficients[1])
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
    
  }
  
  if(type==7)
  {
    tt=array(dim=n,1)
    tt[1:cp]=0
    LNO3<-gamlss(test_data~I(t^0),sigma.formula=~I(tt^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
    print(cp)
    print(LNO3)
    p_KS=ks.test(pnorm(LNO3$residuals),"punif",0,1)$p.value
    print(p_KS)
    
    mu=array(dim=n,LNO3$mu.coefficients[1])
    sigma=array(dim=n,LNO3$sigma.coefficients[1]+LNO3$sigma.coefficients[2]*tt)
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
    
  }
  
  
  if(type==8)
  {
    tt=array(dim=n,1)
    tt[1:cp]=0
    LNO4<-gamlss(test_data~I(tt^1),sigma.formula=~I(tt^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
    print(cp)
    print(LNO4)
    p_KS=ks.test(pnorm(LNO4$residuals),"punif",0,1)$p.value
    print(p_KS)
    
    mu=array(dim=n,LNO4$mu.coefficients[1]+LNO4$mu.coefficients[2]*tt)
    sigma=array(dim=n,LNO4$sigma.coefficients[1]+LNO4$sigma.coefficients[2]*tt)
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  
  #change the output filename###############################
  output=paste("E:/dissertation/3/QUANTILES1/","Q_",data_file[i],sep="")
  write.csv(qQ,output)
  musigma=cbind(musigma,mu,sigma)
  
  
  v1[1,i]=SBC
  v2[1,i]=p_KS
  v3[1,i]=cp
  v4[1,i]=type
  
  
  print("=================Following is the next station============================================")
}
output=paste("E:/dissertation/3/Q_musigma.csv")
write.csv(musigma,output)
result=rbind(v1,v2,v3,v4)
write.csv(result,"E:/dissertation/3/Q_CP.csv")