library(maxLik)
library(gamlss)
library(Kendall)
library(moments)
library(PearsonDS)

#change the work direction##########################################
setwd("E:/dissertation/3/GAMLSS")

data_file=list.files(pattern = "csv")
N=length(data_file)
SBC=array(dim=c(21,N))#5 distributions multiplied by 7 scenarios
p_KS=array(dim=c(21,N))
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
  
  #change the variable######################################
  test_data=PET
  
  
  GA1<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=GA(mu.link = "identity",sigma.link = "identity"))
  GA2<-gamlss(test_data~I(t^1),sigma.formula=~I(t^0),family=GA(mu.link = "identity",sigma.link = "identity"))
  GA3<-gamlss(test_data~I(t^0),sigma.formula=~I(t^1),family=GA(mu.link = "identity",sigma.link = "identity"))
  GA4<-gamlss(test_data~I(t^1),sigma.formula=~I(t^1),family=GA(mu.link = "identity",sigma.link = "identity"))
  GA5<-gamlss(test_data~I(t^1),sigma.formula=~I(t^0),family=GA(mu.link = "log",sigma.link = "log"))
  GA6<-gamlss(test_data~I(t^0),sigma.formula=~I(t^1),family=GA(mu.link = "log",sigma.link = "log"))
  GA7<-gamlss(test_data~I(t^1),sigma.formula=~I(t^1),family=GA(mu.link = "log",sigma.link = "log"))
  
  
  #===========================================
  NO1<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
  NO2<-gamlss(test_data~I(t^1),sigma.formula=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
  NO3<-gamlss(test_data~I(t^0),sigma.formula=~I(t^1),family=NO(mu.link = "identity",sigma.link = "identity"))
  NO4<-gamlss(test_data~I(t^1),sigma.formula=~I(t^1),family=NO(mu.link = "identity",sigma.link = "identity"))
  NO5<-gamlss(test_data~I(t^1),sigma.formula=~I(t^0),family=NO(mu.link = "log",sigma.link = "log"))
  NO6<-gamlss(test_data~I(t^0),sigma.formula=~I(t^1),family=NO(mu.link = "log",sigma.link = "log"))
  NO7<-gamlss(test_data~I(t^1),sigma.formula=~I(t^1),family=NO(mu.link = "log",sigma.link = "log"))
  
  #===========================================
  LNO1<-gamlss(test_data~I(t^0),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
  LNO2<-gamlss(test_data~I(t^1),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
  LNO3<-gamlss(test_data~I(t^0),sigma.formula=~I(t^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
  LNO4<-gamlss(test_data~I(t^1),sigma.formula=~I(t^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
  LNO5<-gamlss(test_data~I(t^1),sigma.formula=~I(t^0),family=LNO(mu.link = "log",sigma.link = "log"))
  LNO6<-gamlss(test_data~I(t^0),sigma.formula=~I(t^1),family=LNO(mu.link = "log",sigma.link = "log"))
  LNO7<-gamlss(test_data~I(t^1),sigma.formula=~I(t^1),family=LNO(mu.link = "log",sigma.link = "log"))
  
  
  
  
  SBC[1,i]=GA1$sbc
  SBC[2,i]=GA2$sbc
  SBC[3,i]=GA3$sbc
  SBC[4,i]=GA4$sbc
  SBC[5,i]=GA5$sbc
  SBC[6,i]=GA6$sbc
  SBC[7,i]=GA7$sbc
  
  
  
  SBC[8,i]=NO1$sbc
  SBC[9,i]=NO2$sbc
  SBC[10,i]=NO3$sbc
  SBC[11,i]=NO4$sbc
  SBC[12,i]=NO5$sbc
  SBC[13,i]=NO6$sbc
  SBC[14,i]=NO7$sbc
  
  SBC[15,i]=LNO1$sbc
  SBC[16,i]=LNO2$sbc
  SBC[17,i]=LNO3$sbc
  SBC[18,i]=LNO4$sbc
  SBC[19,i]=LNO5$sbc
  SBC[20,i]=LNO6$sbc
  SBC[21,i]=LNO7$sbc
  
  
 
  
  p_KS[1,i]=ks.test(pnorm(GA1$residuals),"punif",0,1)$p.value
  p_KS[2,i]=ks.test(pnorm(GA2$residuals),"punif",0,1)$p.value
  p_KS[3,i]=ks.test(pnorm(GA3$residuals),"punif",0,1)$p.value
  p_KS[4,i]=ks.test(pnorm(GA4$residuals),"punif",0,1)$p.value
  p_KS[5,i]=ks.test(pnorm(GA5$residuals),"punif",0,1)$p.value
  p_KS[6,i]=ks.test(pnorm(GA6$residuals),"punif",0,1)$p.value
  p_KS[7,i]=ks.test(pnorm(GA7$residuals),"punif",0,1)$p.value
  
 
  p_KS[8,i]=ks.test(pnorm(NO1$residuals),"punif",0,1)$p.value
  p_KS[9,i]=ks.test(pnorm(NO2$residuals),"punif",0,1)$p.value
  p_KS[10,i]=ks.test(pnorm(NO3$residuals),"punif",0,1)$p.value
  p_KS[11,i]=ks.test(pnorm(NO4$residuals),"punif",0,1)$p.value
  p_KS[12,i]=ks.test(pnorm(NO5$residuals),"punif",0,1)$p.value
  p_KS[13,i]=ks.test(pnorm(NO6$residuals),"punif",0,1)$p.value
  p_KS[14,i]=ks.test(pnorm(NO7$residuals),"punif",0,1)$p.value
  
  p_KS[15,i]=ks.test(pnorm(LNO1$residuals),"punif",0,1)$p.value
  p_KS[16,i]=ks.test(pnorm(LNO2$residuals),"punif",0,1)$p.value
  p_KS[17,i]=ks.test(pnorm(LNO3$residuals),"punif",0,1)$p.value
  p_KS[18,i]=ks.test(pnorm(LNO4$residuals),"punif",0,1)$p.value
  p_KS[19,i]=ks.test(pnorm(LNO5$residuals),"punif",0,1)$p.value
  p_KS[20,i]=ks.test(pnorm(LNO6$residuals),"punif",0,1)$p.value
  p_KS[21,i]=ks.test(pnorm(LNO7$residuals),"punif",0,1)$p.value
  
  
  mn=which(SBC[,i]==min(SBC[,i]))
  
  # if(mn<=7)
  # {
  #   print("GU")
  # }
  # if(mn>7&mn<=14)
  # {
  #   print("WEI")
  # }
  # if(mn>14&mn<=21)
  # {
  #   print("GA")
  # }
  # if(mn>21&mn<=28)
  # {
  #   print("NO")
  # }
  # if(mn>28)
  # {
  #   print("LNO")
  # }
  #=====================================
  
  qQ=array(dim=c(n,6))
  
  

  
  if(mn==1)
  {
    print(GA1)
    mu=array(dim=n,GA1$mu.coefficients[1])
    sigma=array(dim=n,GA1$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==2)
  {
    print(GA2)
    mu=array(dim=n,GA2$mu.coefficients[1])
    sigma=array(dim=n,GA2$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==3)
  {
    print(GA3)
    mu=array(dim=n,GA3$mu.coefficients[1])
    sigma=array(dim=n,GA3$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==4)
  {
    print(GA4)
    mu=array(dim=n,GA4$mu.coefficients[1])
    sigma=array(dim=n,GA4$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==5)
  {
    print(GA5)
    mu=array(dim=n,GA5$mu.coefficients[1])
    sigma=array(dim=n,GA5$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==6)
  {
    print(GA6)
    mu=array(dim=n,GA6$mu.coefficients[1])
    sigma=array(dim=n,GA6$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==7)
  {
    print(GA7)
    mu=array(dim=n,GA7$mu.coefficients[1])
    sigma=array(dim=n,GA7$sigma.coefficients[1])
    qQ[,1]=qGA(0.05,mu,sigma)
    qQ[,2]=qGA(0.25,mu,sigma)
    qQ[,3]=qGA(0.5,mu,sigma)
    qQ[,4]=qGA(0.75,mu,sigma)
    qQ[,5]=qGA(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  
  if(mn==8)
  {
    print(NO1)
    mu=array(dim=n,NO1$mu.coefficients[1])
    sigma=array(dim=n,NO1$sigma.coefficients[1])
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==9)
  {
    print(NO2)
    mu=array(dim=n,NO2$mu.coefficients[1]+NO2$mu.coefficients[2]*t)
    sigma=array(dim=n,NO2$sigma.coefficients[1])
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==10)
  {
    print(NO3)
    mu=array(dim=n,NO3$mu.coefficients[1])
    sigma=array(dim=n,NO3$sigma.coefficients[1]+NO3$sigma.coefficients[2]*t)
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==11)
  {
    print(NO4)
    mu=array(dim=n,NO4$mu.coefficients[1]+NO4$mu.coefficients[2]*t)
    sigma=array(dim=n,NO4$sigma.coefficients[1]+NO4$sigma.coefficients[2]*t)
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==12)
  {
    print(NO5)
    mu=array(dim=n,exp(NO5$mu.coefficients[1]+NO5$mu.coefficients[2]*t))
    sigma=array(dim=n,exp(NO5$sigma.coefficients[1]))
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==13)
  {
    print(NO6)
    mu=array(dim=n,exp(NO6$mu.coefficients[1]))
    sigma=array(dim=n,exp(NO6$sigma.coefficients[1]+NO6$sigma.coefficients[2]*t))
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==14)
  {
    print(NO7)
    mu=array(dim=n,exp(NO7$mu.coefficients[1]+NO7$mu.coefficients[2]*t))
    sigma=array(dim=n,exp(NO7$sigma.coefficients[1]+NO7$sigma.coefficients[2]*t))
    qQ[,1]=qNO(0.05,mu,sigma)
    qQ[,2]=qNO(0.25,mu,sigma)
    qQ[,3]=qNO(0.5,mu,sigma)
    qQ[,4]=qNO(0.75,mu,sigma)
    qQ[,5]=qNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  if(mn==15)
  {
    print(LNO1)
    mu=array(dim=n,LNO1$mu.coefficients[1])
    sigma=array(dim=n,LNO1$sigma.coefficients[1])
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==16)
  {
    print(LNO2)
    mu=array(dim=n,LNO2$mu.coefficients[1]+LNO2$mu.coefficients[2]*t)
    sigma=array(dim=n,LNO2$sigma.coefficients[1])
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==17)
  {
    print(LNO3)
    mu=array(dim=n,LNO3$mu.coefficients[1])
    sigma=array(dim=n,LNO3$sigma.coefficients[1]+LNO3$sigma.coefficients[2]*t)
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==18)
  {
    print(LNO4)
    mu=array(dim=n,LNO4$mu.coefficients[1]+LNO4$mu.coefficients[2]*t)
    sigma=array(dim=n,LNO4$sigma.coefficients[1]+LNO4$sigma.coefficients[2]*t)
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==19)
  {
    print(LNO5)
    mu=array(dim=n,exp(LNO5$mu.coefficients[1]+LNO5$mu.coefficients[2]*t))
    sigma=array(dim=n,exp(LNO5$sigma.coefficients[1]))
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==20)
  {
    print(LNO6)
    mu=array(dim=n,exp(LNO6$mu.coefficients[1]))
    sigma=array(dim=n,exp(LNO6$sigma.coefficients[1]+LNO6$sigma.coefficients[2]*t))
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  if(mn==21)
  {
    print(LNO7)
    mu=array(dim=n,exp(LNO7$mu.coefficients[1]+LNO7$mu.coefficients[2]*t))
    sigma=array(dim=n,exp(LNO7$sigma.coefficients[1]+LNO7$sigma.coefficients[2]*t))
    qQ[,1]=qLNO(0.05,mu,sigma)
    qQ[,2]=qLNO(0.25,mu,sigma)
    qQ[,3]=qLNO(0.5,mu,sigma)
    qQ[,4]=qLNO(0.75,mu,sigma)
    qQ[,5]=qLNO(0.95,mu,sigma)
    qQ[,6]=test_data
  }
  
  print(p_KS[mn,i])
  
  #change the output filename###################################### 
  output=paste("E:/dissertation/3/QUANTILES_TR/","PET_",data_file[i],sep="")
  write.csv(qQ,output)
  v1[1,i]=SBC[mn,i]
  v2[1,i]=p_KS[mn,i]
  v3[1,i]=mn
  print("====================Following is the next station=============================")
}
result=rbind(v1,v2,v3)
write.csv(result,"E:/dissertation/3/PET_TR.csv")