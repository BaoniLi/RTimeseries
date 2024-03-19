library(utils)


stationlevel <- read.csv("E:/dissertation/4/IM.csv", head=T)
co = array(0,8)
pvalue = array(0,8)
pv = array(0,8)
for (i in seq(1,8,1))
{
  tem <- stationlevel[,i]
  im <- stationlevel[,i+1]
  
  aa <- cor.test(tem,im)
  co[i] <- aa$estimate
  pvalue[i] <- aa$p.value
  
  if (pvalue[i]<0.05){
    pv[i] = 0
  } else if (pvalue[i]<0.1){
    pv[i] = 1
  } else {
    pv[i] = 2
  }
}
write.csv(co,"E:/dissertation/CO.csv")
write.csv(pvalue,"E:/dissertation/PV.csv")