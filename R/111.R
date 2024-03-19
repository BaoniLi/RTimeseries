library(utils)

setwd("E:/YANGTZE/scatter")

stationlevel <- read.csv("YRDRH.csv")
vrh <- c(1:17)
vet <- c(1:17)

for (i in seq(1,ncol(stationlevel),2))
{
  rh <- stationlevel[,i]
  et <- stationlevel[,i+1]
  vrh <- c(vrh,rh)
  vet <- c(vet,et)
}
write.csv(vrh,"E:/YANGTZE/scatter/vrh.csv")