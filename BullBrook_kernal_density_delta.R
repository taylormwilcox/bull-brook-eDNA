#############################################
##                                         ##
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 July 2017                            ##
##                                         ##
## Kernal density plot                     ##
##                                         ##
#############################################

### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### some required packages
library("MASS")
library("RColorBrewer")

### ranges for calculations/plotting
xrange <- range(data$S1_93_11)
yrange <- range(data$MS_Hist)
kd_xrange <- xrange + c(-0.5,0.5)
kd_yrange <- yrange + c(-0.5,0.5)

### catageories of detections
BUTnoBRK <- subset(data, data$BUT1_P == 1 & data$BRK2_P == 0)
BUTandBRK <- subset(data, data$BUT1_P == 1 & data$BRK2_P == 1)
BRKnoBUT <- subset(data, data$BUT1_P == 0 & data$BRK2_P == 1)

### densities
z1 <- kde2d(x=data$S1_93_11, y=data$MS_Hist, n=100,lims=c(kd_xrange, kd_yrange))
z2 <- kde2d(x=BUTnoBRK$S1_93_11, y=BUTnoBRK$MS_Hist, n=100,lims=c(kd_xrange, kd_yrange))
z3 <- kde2d(x=BUTandBRK$S1_93_11, BUTandBRK$MS_Hist, n=100,lims=c(kd_xrange, kd_yrange))
z4 <- kde2d(x=BRKnoBUT$S1_93_11, y=BRKnoBUT$MS_Hist, n=100,lims=c(kd_xrange, kd_yrange))

### standardize kernal densities
z1$z <- z1$z/sum(z1$z)
z2$z <- z2$z/sum(z2$z)
z3$z <- z3$z/sum(z3$z)
z4$z <- z4$z/sum(z4$z)

### make delta z (standardized)
z2delta <- (z2$z - z1$z)
z3delta <- (z3$z - z1$z)
z4delta <- (z4$z - z1$z)

### colors for delta plots
rf <- colorRampPalette(c(rev(brewer.pal(9,'Blues')),rep("white",50),brewer.pal(9,'Reds')))
r <- rf(100)

### bull trout only
par(mfrow = c(1,3))
par(mar=c(4,4,0.5,0), las=1,cex=1.2, lwd=2, lend=2)
plot(NULL, NULL, xlim=xrange, ylim=yrange, axes=F, xlab="", ylab="Discharge (CFS)")
image(list(x=z1$x, y=z1$y, z=z2delta), col=r, add=T)
points(MS_Hist ~ S1_93_11, BUTnoBRK, pch=19, cex = 0.2)
abline(h=10, lty=2, lwd=1)
axis(1, at=c(6,8,10,12), lwd=2)
axis(2, at=c(10,50,100), lwd=2)
box()
text(5.7,87,"Bull trout only", adj=c(0,0), cex=1)

### bull and brook trout
par(mar=c(4,2,0.5,2))
plot(NULL, NULL, xlim=xrange, ylim=yrange, xlab=expression(paste("Temperature (",degree,"C)")), ylab="", axes=F)
image(list(x=z1$x, y=z1$y, z=z3delta), col=r, add=T)
points(MS_Hist ~ S1_93_11, BUTandBRK, pch=19, cex = 0.2)
axis(1, at=seq(from=6,to=12,by=2), lwd=2)
abline(h=10, lty=2, lwd=1)
box()
text(5.7,87,"Bull & brook trout", adj=c(0,0))

### brook trout only
par(mar=c(4,0,0.5,4))
plot(NULL, NULL, xlim=xrange, ylim=yrange, xlab="",axes=F)
image(list(x=z1$x, y=z1$y, z=z4delta), col=r, add=T)
points(MS_Hist ~ S1_93_11, BRKnoBUT, pch=19, cex = 0.2)
axis(1, at=seq(from=6,to=12,by=2), lwd=2)
abline(h=10, lty=2, lwd=1)
box()
text(5.7,87,"Brook trout only", adj=c(0,0))


