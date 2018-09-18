#############################################
##                                         ##
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 July 2017                            ##
##                                         ##
## Brook trout eDNA and habitat covariates ##
##                                         ##
#############################################

data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")  

par(mfrow=c(1,3))
par(mar=c(4,4,0.5,0), las=1,cex=1.2, lwd=2, lend=2)
plot(BRK2_Quant ~ MS_Hist, data, xlab="Summer discharge (CFS)", ylab="mtDNA copies/PCR",
     pch=19, col=rgb(0,0,0, alpha=0.2))

par(mar=c(4,2,0.5,2))
plot(BRK2_Quant ~ S1_93_11, data, xlab=expression(paste("August stream temp (C", degree,")")), ylab="", axes=F,
     pch=19, col=rgb(0,0,0, alpha=0.2))
axis(1, at=c(6,8,10,12))
box()

par(mar=c(4,0,0.5,4))
plot(BRK2_Quant ~ SLOPE, data, xlab="Stream gradient", ylab="", axes=F,
     pch=19, col=rgb(0,0,0, alpha=0.2))
axis(1)
box()
