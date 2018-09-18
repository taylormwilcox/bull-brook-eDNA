#############################################
##                                         ##
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 July 2017                            ##
##                                         ##
## Model response curves                   ##
##                                         ##
#############################################

### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### some required packages
require(lme4)
require(grDevices)

### bring in fit model and patch influence object
load("Q_model_fits.rda")     # previously-fit models
model_Q <- Q_model_fits[[3]] # top model
load("Q_influence.rda")      # inf object

### user-defined function to create response curve given a matrix "x" of made-up covariate values
pred_BULL <- function(x){
  pred <- c()
  for(i in 1:nrow(x)){
    p <-model_Q@beta[1] +
        model_Q@beta[2]*log(x$flow[i]) + 
        model_Q@beta[3]*x$temp[i] +
        model_Q@beta[4]*log(x$brk[i]+1) +
        model_Q@beta[5]*log(x$brk[i]+1)^2 +
        model_Q@beta[6]*log(x$brk[i]+1)*log(x$flow[i])
    pred[i] <- exp(p)/(1+exp(p))
  }
  print(pred)
}

### Pr(bull trout) ~ FLOW
### user-defined function to create response curves with FLOW as the x-axis; different curve for each patch left out in model fit
pred_BULL_inf_FLOW <- function(x, c){
  
  for(j in 1:nrow(Q_influence$alt.fixed)){
    betas <- as.vector(Q_influence$alt.fixed[j,])
    pred <- c()
    for(i in 1:nrow(x)){
      p <- betas[1] +
        betas[2]*log(x$flow[i]) + 
        betas[3]*x$temp[i] +
        betas[4]*log(x$brk[i]+1) +
        betas[5]*log(x$brk[i]+1)^2 +
        betas[6]*log(x$brk[i]+1)*log(x$flow[i])
      pred[i] <- exp(p)/(1+exp(p))
    }
    lines(pred ~ x$flow, col=c)
  }
}

### plotting area   
par(mfrow=c(1,2),mar=c(4,4,0.5,0.5),las=1, cex=1.2, lend=2)
plot(NULL, NULL, xlim=range(data$MS_Hist), ylim=c(0,1),
     xlab=paste("Mean summer discharge (CFS)"), ylab="Pr(Bull trout)", axes=F)
axis(1, lwd=2, cex=1.5)
axis(2, at=seq(from=0, to=1, by=0.2), lwd=2)
box(,lwd=2)

### label for subplot
text(5,.95,"a.", cex=1.2)

### dataframe for median temp and no brook v flow
flow <- seq(from=min(data$MS_Hist), to=max(data$MS_Hist), by = 0.1)
temp <- rep(median(data$S1_93_11), length(flow))
brk  <- rep(0, length(flow))
medTemp_noBRK <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_FLOW(x=medTemp_noBRK, c=adjustcolor("black", alpha=0.1)) # dropping patches
medTemp_noBRK$bull <- pred_BULL(medTemp_noBRK) # full data set
lines(bull ~ flow, medTemp_noBRK, lwd=2)

### add a few brook trout
data_10 <- subset(data, data$BRK2_Quant > 0 & data$BRK2_Quant <= 100) # limit lines to observed range
flow <- seq(from=min(data_10$MS_Hist), to=max(data_10$MS_Hist), by = 0.1)
temp <- rep(median(data$S1_93_11), length(flow))
brk  <- rep(10, length(flow))
medTemp_10BRK <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_FLOW(x=medTemp_10BRK, c=adjustcolor("darkorange", alpha=0.1)) 
medTemp_10BRK$bull <- pred_BULL(medTemp_10BRK)
lines(bull ~ flow, medTemp_10BRK, lwd=2, col="darkorange")

### add a bunch of brook trout
data_100 <- subset(data, data$BRK2_Quant > 100 & data$BRK2_Quant <= 200) 
flow <- seq(from=min(data_100$MS_Hist), to=max(data_100$MS_Hist), by = 0.1)
temp <- rep(median(data$S1_93_11), length(flow))
brk  <- rep(100, length(flow))
medTemp_100BRK <- as.data.frame(cbind(flow,temp,brk))
#pred_BULL_inf_FLOW(x=medTemp_100BRK, c=adjustcolor("darkblue", alpha=0.2)) 
medTemp_100BRK$bull <- pred_BULL(medTemp_100BRK)
lines(bull ~ flow, medTemp_100BRK, lwd=2, col="darkblue")

###EXTRAPLOATED###
### add a bunch of brook trout
flow <- seq(from=min(data$MS_Hist), to=max(data$MS_Hist), by = 0.1)
brk  <- rep(100, length(flow))
medTemp_100BRK <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_FLOW(x=medTemp_100BRK, c=adjustcolor("darkblue", alpha=0.1)) 
medTemp_100BRK$bull <- pred_BULL(medTemp_100BRK)
lines(bull ~ flow, medTemp_100BRK, lwd=2, col="darkblue",lty=2)

### add a terrible number of brook trout
data_200 <- subset(data, data$BRK2_Quant > 200 & data$BRK2_Quant <= 500)
flow <- seq(from=min(data_200$MS_Hist), to=max(data_200$MS_Hist), by = 0.1)
temp <- rep(median(data$S1_93_11), length(flow))
brk  <- rep(200, length(flow))
medTemp_200BRK <- as.data.frame(cbind(flow,temp,brk))
#pred_BULL_inf_FLOW(x=medTemp_200BRK, c=adjustcolor("darkred", alpha=0.2)) 
medTemp_200BRK$bull <- pred_BULL(medTemp_200BRK)
lines(bull ~ flow, medTemp_200BRK, lwd=2, col="darkred")

###EXTRAPLOATED###
### add a bunch of brook trout
flow <- seq(from=min(data$MS_Hist), to=max(data$MS_Hist), by = 0.1)
brk  <- rep(200, length(flow))
temp <- rep(median(data$S1_93_11), length(flow))
medTemp_200BRK <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_FLOW(x=medTemp_200BRK, c=adjustcolor("darkred", alpha=0.1)) 
medTemp_200BRK$bull <- pred_BULL(medTemp_200BRK)
lines(bull ~ flow, medTemp_200BRK, lwd=2, col="darkred",lty=2)

### add a horrendous number of brook trout
data_500 <- subset(data, data$BRK2_Quant > 500)
flow <- seq(from=min(data_500$MS_Hist), to=max(data_500$MS_Hist), by = 0.1)
temp <- rep(median(data$S1_93_11), length(flow))
brk  <- rep(500, length(flow))
medTemp_500BRK <- as.data.frame(cbind(flow,temp,brk))
#pred_BULL_inf_FLOW(x=medTemp_500BRK, c=adjustcolor("purple", alpha=0.2)) 
medTemp_500BRK$bull <- pred_BULL(medTemp_500BRK)
lines(bull ~ flow, medTemp_500BRK, lwd=2, col="purple")

###EXTRAPLOATED###
### add a bunch of brook trout
flow <- seq(from=min(data$MS_Hist), to=max(data$MS_Hist), by = 0.1)
brk  <- rep(500, length(flow))
temp <- rep(median(data$S1_93_11), length(flow))
medTemp_500BRK <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_FLOW(x=medTemp_500BRK, c=adjustcolor("purple", alpha=0.1)) 
medTemp_500BRK$bull <- pred_BULL(medTemp_500BRK)
lines(bull ~ flow, medTemp_500BRK, lwd=2, col="purple",lty=2)

### legend
#text(80,0.58,"Brook tr DNA")
#text(79,.5,"copies/PCR", cex=1)
legend(x=70, y=.4, pch=19, col=c("black","darkorange","darkblue","darkred","purple"),
legend=c("0  ","10 ","100","200","500"),
cex=0.5,
title="copies/PCR")

### Pr(bull tr) ~ TEMP
### user-defined function for influence object and temp response
pred_BULL_inf_TEMP <- function(x, c){
  
  for(j in 1:nrow(Q_influence$alt.fixed)){
    betas <- as.vector(Q_influence$alt.fixed[j,])
    pred <- c()
    for(i in 1:nrow(x)){
      p <- betas[1] +
        betas[2]*log(x$flow[i]) + 
        betas[3]*x$temp[i] +
        betas[4]*log(x$brk[i]+1) +
        betas[5]*log(x$brk[i]+1)^2 +
        betas[6]*log(x$brk[i]+1)*log(x$flow[i])
      pred[i] <- exp(p)/(1+exp(p))
    }
    lines(pred ~ temp, col=c)
  }
}

### plot area
par(mar=c(4,2,0.5,2.5))
plot(NULL, NULL, xlim=range(data$S1_93_11), ylim=c(0,1), # limit to 10 CFS
     xlab=expression(paste("Mean Aug temperature (",degree,"C)")), ylab="", axes=F)
axis(1, lwd=2, at=seq(from=6,to=12,by=2))
#axis(2, at=seq(from=0, to=1, by=0.2), lwd=2)
box(lwd=2)
text(5.9,.95,"b.", cex=1.2)

### no brook and 10 CFS
data_0 <- subset(data, data$BRK2_Quant == 0) # range with no brook
temp <- seq(from=min(data_0$S1_93_11), to=max(data_0$S1_93_11), by = 0.1)
flow <- rep(10, length(temp))
brk  <- rep(0, length(flow))
FLOW10_BRK0 <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_TEMP(FLOW10_BRK0, c=adjustcolor("black", alpha=0.2))
FLOW10_BRK0$bull <- pred_BULL(FLOW10_BRK0)
lines(bull ~ temp, FLOW10_BRK0, lwd=2, col="black")

### 10 copies BRK and 10 CFS
temp <- seq(from=min(data_10$S1_93_11), to=max(data_10$S1_93_11), by = 0.1)
flow <- rep(10, length(temp))
brk  <- rep(10, length(temp))
FLOW10_BRK10 <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_TEMP(FLOW10_BRK10, c=adjustcolor("darkorange", alpha=0.2))
FLOW10_BRK10$bull <- pred_BULL(FLOW10_BRK10)
lines(bull ~ temp, FLOW10_BRK10, lwd=2, col="darkorange")

### 100 copies BRK and 10 CFS
temp <- seq(from=min(data_100$S1_93_11), to=max(data_100$S1_93_11), by = 0.1)
range(temp)
flow <- rep(10, length(temp))
brk  <- rep(100, length(temp))
FLOW10_BRK100 <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_TEMP(FLOW10_BRK100, c=adjustcolor("darkblue", alpha=0.2))
FLOW10_BRK100$bull <- pred_BULL(FLOW100_BRK10)
lines(bull ~ temp, FLOW10_BRK100, lwd=2, col="darkblue")

### 200 copies BRK and 10 CFS
temp <- seq(from=min(data_200$S1_93_11), to=max(data_200$S1_93_11), by = 0.1)
flow <- rep(10, length(temp))
brk  <- rep(200, length(temp))
FLOW10_BRK200 <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_TEMP(FLOW10_BRK200, c=adjustcolor("darkred", alpha=0.2))
FLOW10_BRK200$bull <- pred_BULL(FLOW10_BRK200)
lines(bull ~ temp, FLOW10_BRK200, lwd=2, col="darkred")

### 500 copies BRK and 10 CFS
temp <- seq(from=min(data_500$S1_93_11), to=max(data_500$S1_93_11), by = 0.1)
flow <- rep(10, length(temp))
brk  <- rep(500, length(temp))
FLOW10_BRK500 <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_TEMP(FLOW10_BRK500, c=adjustcolor("purple", alpha=0.2))
FLOW10_BRK500$bull <- pred_BULL(FLOW10_BRK500)
lines(bull ~ temp, FLOW10_BRK500, lwd=2, col="purple")

