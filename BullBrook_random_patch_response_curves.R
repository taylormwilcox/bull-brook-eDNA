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

### Some figures showing the random intercepts and slopes
patch_intercepts <- ranef(model_Q)[[1]][,1] # patch intercepts
patch_flow <- ranef(model_Q)[[1]][,2] # patch slopes for flow
patch_brk <- ranef(model_Q)[[1]][,3] # for brook trout

### brook trout quant on x axis
par(mfrow = c(1,3))
par(mar=c(4,4,0.5,0), las=1,cex=1.2, lwd=2, lend=2)
plot(NULL, NULL, xlim=range(data$BRK2_Quant), ylim=c(0,1),
     xlab="BRK mtDNA copies/PCR", ylab="Pr(Bull trout)", axes=F)
axis(1, at=seq(from=0, to=800, by=200))
axis(2, at=seq(from=0, to=1, by=0.2))
box()

t <- median(data$S1_93_11) # median temperature (9.1 C)
f <- median(data$MS_Hist)  # median flow (3.8 CFS)
for(j in 1:length(patch_intercepts)){
  
  patch_data <- subset(data, data$Patch == levels(data$Patch)[j]) # patches must be in same order
  brk <- seq(from=min(patch_data$BRK2_Quant), to=max(patch_data$BRK2_Quant), by=1)
  
  pred <- c()
  for(i in 1:length(brk)){
    p <-model_Q@beta[1] + patch_intercepts[j] +
        model_Q@beta[2]*log(f) + patch_flow[j]*log(f) +
        model_Q@beta[3]*t +
        model_Q@beta[4]*log(brk[i]+1) + patch_brk[j]*log(brk[i]+1) +
        model_Q@beta[5]*log(brk[i]+1)^2 +
        model_Q@beta[6]*log(brk[i]+1)*log(f)
    pred[i] <- exp(p)/(1+exp(p))
  }
  
  lines(pred ~ brk, col=adjustcolor("black", alpha=0.2))
}

### line for the mean
brk <- seq(from=0, to=max(data$BRK2_Quant), by=1)
pred <- c()
for(i in 1:length(brk)){
  p <-model_Q@beta[1] +
    model_Q@beta[2]*log(f) + 
    model_Q@beta[3]*t +
    model_Q@beta[4]*log(brk[i]+1) + 
    model_Q@beta[5]*log(brk[i]+1)^2 +
    model_Q@beta[6]*log(brk[i]+1)*log(f)
  pred[i] <- exp(p)/(1+exp(p))
}
lines(pred ~ brk, lwd=2, col="black")

### versus flow
par(mar=c(4,2,0.5,2))
plot(NULL, NULL, xlim=range(data$MS_Hist), ylim=c(0,1),
     xlab="Summer discharge (CFS)", ylab="Pr(Bull trout)", axes=F)
axis(1)
box()

t <- median(data$S1_93_11) # median temperature (9.1 C)
brk <- 0

for(j in 1:length(patch_intercepts)){
  
  patch_data <- subset(data, data$Patch == levels(data$Patch)[j]) # patches must be in same order
  f <- seq(from=min(patch_data$MS_Hist), to=max(patch_data$MS_Hist), by=0.1)
  
  pred <- c()
  for(i in 1:length(f)){
    p <-model_Q@beta[1] + patch_intercepts[j] +
      model_Q@beta[2]*log(f[i]) + patch_flow[j]*log(f[i]) +
      model_Q@beta[3]*t +
      model_Q@beta[4]*log(brk+1) + patch_brk[j]*log(brk+1) +
      model_Q@beta[5]*log(brk+1)^2 +
      model_Q@beta[6]*log(brk+1)*log(f[i])
    pred[i] <- exp(p)/(1+exp(p))
  }
  
  lines(pred ~ f, col=adjustcolor("black", alpha=0.2))
}

### line for the mean
f <- seq(from=min(data$MS_Hist), to=max(data$MS_Hist), by=0.1)
pred <- c()
for(i in 1:length(f)){
  p <-model_Q@beta[1] +
      model_Q@beta[2]*log(f[i]) + 
      model_Q@beta[3]*t +
      model_Q@beta[4]*log(brk+1) + 
      model_Q@beta[5]*log(brk+1)^2 +
      model_Q@beta[6]*log(brk+1)*log(f[i])
  pred[i] <- exp(p)/(1+exp(p))
}
lines(pred ~ f, lwd=2, col="black")

### versus temp
par(mar=c(4,0,0.5,4))
plot(NULL, NULL, xlim=range(data$S1_93_11), ylim=c(0,1),
     xlab="August stream temp (C)", ylab="Pr(Bull trout)", axes=F)
axis(1, at=c(6,8,10,12))
box()

f <- median(data$MS_Hist) # median flow
brk <- 0

for(j in 1:length(patch_intercepts)){
  
  patch_data <- subset(data, data$Patch == levels(data$Patch)[j]) # patches must be in same order
  t <- seq(from=min(patch_data$S1_93_11), to=max(patch_data$S1_93_11), by=0.1)
  
  pred <- c()
  for(i in 1:length(t)){
    p <-model_Q@beta[1] + patch_intercepts[j] +
        model_Q@beta[2]*log(f) + patch_flow[j]*log(f) +
        model_Q@beta[3]*t[i] +
        model_Q@beta[4]*log(brk+1) + patch_brk[j]*log(brk+1) +
        model_Q@beta[5]*log(brk+1)^2 +
        model_Q@beta[6]*log(brk+1)*log(f)
    pred[i] <- exp(p)/(1+exp(p))
  }
  
  lines(pred ~ t, col=adjustcolor("black", alpha=0.2))
}

### line for the mean
t <- seq(from=min(data$S1_93_11), to=max(data$S1_93_11), by=0.1)
pred <- c()
for(i in 1:length(t)){
  p <-model_Q@beta[1] +
      model_Q@beta[2]*log(f) + 
      model_Q@beta[3]*t[i] +
      model_Q@beta[4]*log(brk+1) + 
      model_Q@beta[5]*log(brk+1)^2 +
      model_Q@beta[6]*log(brk+1)*log(f)
  pred[i] <- exp(p)/(1+exp(p))
}
lines(pred ~ t, lwd=2, col="black")
