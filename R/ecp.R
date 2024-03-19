library("ecp")

setwd("E:/YANGTZE")

stationlevel <- read.csv("CP.csv")
cp = array(0,16)
cp1 = array(0,16)
cp2 = array(0,16)
cp3 = array(0,16)

for (i in seq(1,ncol(stationlevel),1))
{
       rh <- stationlevel[,i]
       # tem <- stationlevel[,i+1]
       # q <- stationlevel[,i+2]
       # ea <- stationlevel[,i+3]
       # vpd <- stationlevel[,i+4]
       
       RH <- matrix(c(rh), ncol = 1)
       # TEM <- matrix(c(tem), ncol = 1)
       # Q <- matrix(c(q), ncol = 1)
       # EA <- matrix(c(ea), ncol = 1)
       # VPD <- matrix(c(vpd), ncol = 1)
      
       output1 <- e.divisive(RH, sig.lvl = 0.05, R = 499, k = NULL, min.size = 2, alpha = 1)
       # output2 <- e.divisive(TEM, sig.lvl = 0.05, R = 499, k = NULL, min.size = 2, alpha = 1)
       # output3 <- e.divisive(Q, sig.lvl = 0.05, R = 499, k = NULL, min.size = 2, alpha = 1)
       # output4 <- e.divisive(EA, sig.lvl = 0.05, R = 499, k = NULL, min.size = 2, alpha = 1)
       # output5 <- e.divisive(VPD, sig.lvl = 0.05, R = 499, k = NULL, min.size = 2, alpha = 1)
       
       cp[i] = output1$estimates[2] + 1964
       cp1[i] = output1$estimates[3] + 1964
       cp2[i] = output1$estimates[4] + 1964
       cp3[i] = output1$estimates[5] + 1964
       # cp[i+3] = output4$estimates[2] + 1964
       # cp[i+4] = output5$estimates[2] + 1964
       # ts.plot(RH, ylab = "Annual mean temperature", main = "Change in Annual mean temperature")
       # abline(v = output1$estimates, col = "red", lty = 2)
}