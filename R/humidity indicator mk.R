library(trend)
library(utils)
library(ecp)

setwd("D:/202111/NH/STATION-LEVEL")

stationlevel <- read.table("te.txt", header = F)
cp = array(0,16)
pv = array(0,16)
slope = array(0,16)
x = c(1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017)
cp0 = array(0,16)
cp1 = array(0,16)
cp2 = array(0,16)
cp3 = array(0,16)

for (i in seq(1,nrow(stationlevel),1))
{
  rh <- stationlevel[i,]
  rh2 <- as.vector(as.matrix(rh))
  
  mk <- mk.test(rh2)
  cp[i] <- mk$p.value
  if (cp[i]<0.1){
    pv[i] = 1
  } else {
    pv[i] = 0
  }
  
  fit <- lm(rh2~x)
  slope[i] <- fit$coefficients[2]
  
  RH <- matrix(c(rh), ncol = 1)
  output1 <- e.divisive(RH, sig.lvl = 0.05, R = 499, k = NULL, min.size = 2, alpha = 1)
  cp0[i] = output1$estimates[2] + 1964
  cp1[i] = output1$estimates[3] + 1964
  cp2[i] = output1$estimates[4] + 1964
  cp3[i] = output1$estimates[5] + 1964
}
# write.csv(slope,"111.csv")