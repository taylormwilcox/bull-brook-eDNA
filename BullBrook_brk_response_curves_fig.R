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
load("Q_model_fits.rda") # previously-fit models
model_Q <- Q_model_fits[[3]] # top model
load("Q_influence.rda") # inf object

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

### Pr(bull tr) ~ brook trout eDNA
### user-defined function for influence object curves
pred_BULL_inf_BRK <- function(x, c){
  
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
    lines(pred ~ x$brk, col=c)
  }
}

### plot space
plot(NULL, NULL, xlim=c(0,max(data$BRK2_Quant)), ylim=c(0,1),
     xlab="Brook trout eDNA copies/PCR", ylab="Pr(Bull trout)", axes=F)
axis(1, at=seq(from=0, to=800, by=200), lwd=2)
axis(2, at=seq(from=0, to=1, by=0.2), lwd=2)
box(lwd=2)

### low flow median temp
data_1FLOW <- subset(data, data$MS_Hist <= 5)
brk  <- seq(from=min(data_1FLOW$BRK2_Quant), to=max(data_1FLOW$BRK2_Quant), by = 1)
flow <- rep(1, length(brk))
temp <- rep(median(data$S1_93_11), length(brk))
medTemp_1FLOW <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_BRK(medTemp_1FLOW, c=adjustcolor("purple", alpha=0.2))
medTemp_1FLOW$bull <- pred_BULL(medTemp_1FLOW)
lines(bull ~ brk, medTemp_1FLOW, lwd=2, col="purple")

### medflow at 5 CFS
data_5FLOW <- subset(data, data$MS_Hist >5 & data$MS_Hist <= 10) 
brk  <- seq(from=min(data_10FLOW$BRK2_Quant), to=max(data_10FLOW$BRK2_Quant), by = 1)
flow <- rep(5, length(brk))
temp <- rep(median(data$S1_93_11), length(brk))
medTemp_5FLOW <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_BRK(medTemp_5FLOW, c=adjustcolor("darkred", alpha=0.2))
medTemp_5FLOW$bull <- pred_BULL(medTemp_5FLOW)
lines(bull ~ brk, medTemp_5FLOW, lwd=2, col="darkred")

### highish flow at 10 CFS
data_10FLOW <- subset(data, data$MS_Hist >10 & data$MS_Hist <= 25) 
brk  <- seq(from=min(data_10FLOW$BRK2_Quant), to=max(data_10FLOW$BRK2_Quant), by = 1)
flow <- rep(10, length(brk))
temp <- rep(median(data$S1_93_11), length(brk))
medTemp_10FLOW <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_BRK(medTemp_10FLOW, c=adjustcolor("darkblue", alpha=0.2))
medTemp_10FLOW$bull <- pred_BULL(medTemp_10FLOW)
lines(bull ~ brk, medTemp_10FLOW, lwd=2, col="darkblue")

### higher flow 25 CFS
data_25FLOW <- subset(data, data$MS_Hist >25 & data$MS_Hist <= 50) 
brk  <- seq(from=min(data_25FLOW$BRK2_Quant), to=max(data_25FLOW$BRK2_Quant), by = 1)
flow <- rep(25, length(brk))
temp <- rep(median(data$S1_93_11), length(brk))
medTemp_25FLOW <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_BRK(medTemp_25FLOW, c=adjustcolor("darkgreen", alpha=0.2))
medTemp_25FLOW$bull <- pred_BULL(medTemp_25FLOW)
lines(bull ~ brk, medTemp_25FLOW, lwd=2, col="darkgreen")

### highest flow
data_50FLOW <- subset(data, data$MS_Hist >50) 
brk  <- seq(from=min(data_50FLOW$BRK2_Quant), to=max(data_50FLOW$BRK2_Quant), by = 1)
flow <- rep(50, length(brk))
temp <- rep(median(data$S1_93_11), length(brk))
medTemp_50FLOW <- as.data.frame(cbind(flow,temp,brk))
pred_BULL_inf_BRK(medTemp_50FLOW, c=adjustcolor("black", alpha=0.2))
medTemp_50FLOW$bull <- pred_BULL(medTemp_50FLOW)
lines(bull ~ brk, medTemp_50FLOW, lwd=2, col="black")

### legend
legend("topright", bty="n", pch=19, col=c("black","darkgreen","darkblue","darkred","purple"),
       legend=c("50 CFS","25 CFS","10 CFS","5 CFS","1 CFS"))

