library(maxLik)
library(gamlss)
library(copula)
library(moments)
library(PearsonDS)
library(CDVine)
library(trend)
library(VineCopula)
library(maxLik)
library(graphics)
library(nlme)


#读取数据
setwd("D:/长江流域/MKTEST")
data<-read.table(file="E.txt",head=T)
E=data$E

n=length(E)
t=c(1:n)

mk.test(E)
pettitt.test(E)
#基于时变矩模型的序列趋势性诊断
#选择分布类型
#1
BIC_EG=c(-999,-999,-999,-999,-999,-999)
EG<-gamlss(E~I(t^0),sigma.fo=~I(t^0),family=WEI(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
BIC_EG[1]=EG$sbc
#2
EG<-gamlss(E~I(t^0),sigma.fo=~I(t^0),family=GU(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
BIC_EG[2]=EG$sbc
#3
EG<-gamlss(E~I(t^0),sigma.fo=~I(t^0),family=GA(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
BIC_EG[3]=EG$sbc
#4
EG<-gamlss(E~I(t^0),sigma.fo=~I(t^0),family=LO(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
BIC_EG[4]=EG$sbc
#5
EG<-gamlss(E~I(t^0),sigma.fo=~I(t^0),family=NO(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
BIC_EG[5]=EG$sbc
#6
EG<-gamlss(E~I(t^0),sigma.fo=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
BIC_EG[6]=EG$sbc

min(BIC_EG)
uE=pLNO(E,muE,sigmaE)
ks.test(uE,"punif",0,1)

#趋势性分析
BIC_EG1=c(-999,-999)
EG<-gamlss(E~I(t^1),sigma.formula=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
BIC_EG1[1]=EG$sbc
muE=array(dim=n,EG$mu.coefficients[2]*t+EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])

EG<-gamlss(E~I(t^0),sigma.formula=~I(t^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
BIC_EG1[2]=EG$sbc
muE=array(dim=n,EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1]+EG$sigma.coefficients[2]*t)

uE=pLNO(E,muE,sigmaE)
ks.test(uE,"punif",0,1)

min(BIC_PG1)

#基于时变矩模型的变点诊断
BIC_CE=array(dim=n,-999)
BIC0=669.3723
  
for (i in 10:40) {
  
tte1=array(dim=n,1)
tte1[1:i]=0
EG<-gamlss(E~I(tte1^1),sigma.fo=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
BIC_CE[i]=EG$sbc
if(BIC0>BIC_CE[i])
{
  BIC0=BIC_CE[i]
  print(i)
}
}

for (i in 10:40) {
  
  tte2=array(dim=n,1)
  tte2[1:i]=0
  EG<-gamlss(E~I(t^0),sigma.fo=~I(tte2^1),family=LNO(mu.link = "identity",sigma.link = "identity"))
  BIC_CE[i]=EG$sbc
  if(BIC0>BIC_CE[i])
  {
    BIC0=BIC_CE[i]
    print(i)
  }
}
BIC_CE[26]
tte1=array(dim=n,1)
tte1[1:26]=0
EG<-gamlss(E~I(tte1^1),sigma.fo=~I(t^0),family=LNO(mu.link = "identity",sigma.link = "identity"))
muE=array(dim=n,EG$mu.coefficients[2]*tte1+EG$mu.coefficients[1])
sigmaE=array(dim=n,EG$sigma.coefficients[1])
uE=pLNO(E,muE,sigmaE)
BIC_CE1=EG$sbc
ks.test(uE,"punif",0,1)
E0.5=qLO(0.5,muE,sigmaE)
write.csv(E0.5,"E0.5.csv")

ttq2=array(dim=n,1)
ttq2[1:32]=0
QG<-gamlss(Q~I(t^0),sigma.fo=~I(ttq2^1),family=LO(mu.link = "identity",sigma.link = "identity"))
muQ=array(dim=n,QG$mu.coefficients[2]*ttq+QG$mu.coefficients[1])
sigmaQ=array(dim=n,QG$sigma.coefficients[1])
uQ=pLO(Q,muQ,sigmaQ)
ks.test(uQ,"punif",0,1)