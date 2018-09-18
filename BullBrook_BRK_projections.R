#############################################
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 Mar 18                               ##
##                                         ##
## Brook trout quant model and climate proj##       
##                                         ##
#############################################

### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### some required packages and functions
require("lme4") # for model fitting
require("merTools") # model assessment
require("MuMIn") # model assessment
require("RColorBrewer") # colors
require("pROC") # for model assessment
antilogit <- function(x) { exp(x) / (1 + exp(x) ) } # back transformation

### Exploratory analysis
plot(log(BRK2_Quant + 1) ~ SLOPE, data)
lines(lowess(log(data$BRK2_Quant + 1) ~ data$SLOPE), col='red')

plot(log(BRK2_Quant + 1) ~ S1_93_11, data)
lines(lowess(log(data$BRK2_Quant + 1) ~ data$S1_93_11), col='red') 

plot(log(BRK2_Quant + 1) ~ log(MS_Hist), data)
lines(lowess(log(data$BRK2_Quant + 1) ~ log(data$MS_Hist)), col='red')

### FIT AND PROJECT BROOK TROUT ################################################
### Code for sites within patches with brook trout
patch_occ <- aggregate(data$BRK2_P, by = list(data$Patch), FUN = sum)           
patch_occ <- subset(patch_occ, patch_occ$x > 0)$Group.1                         # patch has at least one positive
data$patch_code <- ifelse(data$Patch %in% patch_occ, 1, 0)
sum(data$patch_code)                                                            # 476 sites within these patches
data_brk <- subset(data, data$patch_code == 1)                                  # subsetted data
data_brk <- droplevels(data_brk)
data_no_brk <- subset(data, data$patch_code == 0)
nrow(data_brk)

### Fit a brook trout model; only wth sites within patches with brook trout
brk1 <- lmer(log(BRK2_Quant + 1) ~ S1_93_11 + MS_Hist + SLOPE + (1 | Patch), 
             data = data_brk)
RMSE.merMod(brk1) # RMSE = 1.464
r.squaredGLMM(brk1) # conditional Rsq = 0.521
plot(log(data_brk$BRK2_Quant + 1) ~ fitted(brk1), # figure of predicted versus observed
     xlab="predicted ln(mtDNA copies)", ylab="observed ln(mtDNA copies)")
lines(lowess(log(data_brk$BRK2_Quant + 1) ~ fitted(brk1)), col='red', lwd=3)
abline(0,1)

### Assess brook trout model
brk_model_ci <- confint.merMod(brk1, 
                             level=0.95,
                             method="boot",
                             nsim=200)
summary(brk1)
brk_model_ci

### Use model to predict scenarios: Brook trout stay within patches (historical, 2040, 2080)
### Use same model to predict future scenarios

### Function to predict brook trout
brk_pred <- function(temp, flow, slope){                                        
  brk1@beta[1] +                                                                
    brk1@beta[2] * temp[i] +
    brk1@beta[3] * flow[i] +
    brk1@beta[4] * slope[i]
}

brk_stay_hist <- c()                                                            # to fill
brk_stay_2040 <- c()
brk_stay_2080 <- c()

for(i in 1:nrow(data)){
  if(data$patch_code[i] == 1){
    brk_stay_hist[i] <- brk_pred(temp = data$S1_93_11, flow = data$MS_Hist, slope = data$SLOPE)
    brk_stay_2040[i] <- brk_pred(temp = data$S30_2040D.N.19.11, flow = data$MS_2040.N.16.6, slope = data$SLOPE)
    brk_stay_2080[i] <- brk_pred(temp = data$S32_2080D.N.19.11, flow = data$MS_2080.N.16.6, slope = data$SLOPE)
  } else{
    brk_stay_hist[i] <- 0
    brk_stay_2040[i] <- 0
    brk_stay_2080[i] <- 0
  }
}

data$brk_stay_hist <- ifelse(exp(brk_stay_hist) - 1 >= 0, brk_stay_hist, 0)     # convert values less than 1 to zero
data$brk_stay_2040 <- ifelse(exp(brk_stay_2040) - 1 >= 0, brk_stay_2040, 0)
data$brk_stay_2080 <- ifelse(exp(brk_stay_2080) - 1 >= 0, brk_stay_2080, 0)

### SIMULATIONS OF BROOK TROUT #################################################
### Prediction intervals to estimate standard deviation
brk1_interval <- predictInterval(brk1, newdata = data_brk, which = "fixed", level = 0.68)
brk1_sd <- brk1_interval[,2] - brk1_interval[,1] # aprox. standard deviation
brk1_sd_frame <- data.frame(ID_Tag = data_brk$ID_Tag, BRK_SD = brk1_sd)         # attach to dataframe
data <- merge(data,brk1_sd_frame, by = "ID_Tag", all = T)

### Simulate brook trout 
sim_brk_stay_hist <- matrix(, ncol = 500, nrow = nrow(data))
sim_brk_stay_2040 <- matrix(, ncol = 500, nrow = nrow(data))
sim_brk_stay_2080 <- matrix(, ncol = 500, nrow = nrow(data))

for(i in 1:nrow(data)){
  if(data$patch_code[i] == 1){
    sim_brk_stay_hist[i,] <- rnorm(n = 500, 
                                       mean = data$brk_stay_hist[i], 
                                       sd = data$BRK_SD[i])
    sim_brk_stay_2040[i,] <- rnorm(n = 500, 
                                   mean = data$brk_stay_2040[i], 
                                   sd = data$BRK_SD[i])
    sim_brk_stay_2080[i,] <- rnorm(n = 500, 
                                   mean = data$brk_stay_2080[i], 
                                   sd = data$BRK_SD[i])
  } else{
    sim_brk_stay_hist[i,] <- rep(0, 500)
    sim_brk_stay_2040[i,] <- rep(0, 500)
    sim_brk_stay_2080[i,] <- rep(0, 500)
  }
}

sim_brk_stay_hist <- ifelse(exp(sim_brk_stay_hist) - 1 < 0, 0, sim_brk_stay_hist) # remove less than zero values
sim_brk_stay_2040 <- ifelse(exp(sim_brk_stay_2040) - 1 < 0, 0, sim_brk_stay_2040)
sim_brk_stay_hist <- ifelse(exp(sim_brk_stay_2040) - 1 < 0, 0, sim_brk_stay_2040) 

### Now use these to project bull trout
load("Q_model_fits.rda")                                                        # top bull trout model
Q_model <- Q_model_fits[[3]] 

pred_fun <- function(flow, temp, brk){                                          # prediction function
  print(Q_model@beta[1] + 
          Q_model@beta[2] * log(flow[i]) +
          Q_model@beta[3] * temp[i] +
          Q_model@beta[4] * brk[i] + 
          Q_model@beta[5] * brk[i]^2 + 
          Q_model@beta[6] * log(flow[i]) * brk[i])
}

sim_but_hist <- matrix(, ncol = 500, nrow = nrow(data))
sim_but_2040 <- matrix(, ncol = 500, nrow = nrow(data))
sim_but_2080 <- matrix(, ncol = 500, nrow = nrow(data))
for(j in 1:500){
  for(i in 1:nrow(data)){
    sim_but_hist[i,j] <- antilogit(pred_fun(flow = data$MS_Hist, temp = data$S1_93_11, brk = sim_brk_stay_hist[,j]))
    sim_but_2040[i,j] <- antilogit(pred_fun(flow = data$MS_2040.N.16.6, temp = data$S30_2040D.N.19.11, brk = sim_brk_stay_2040[,j]))
    sim_but_2080[i,j] <- antilogit(pred_fun(flow = data$MS_2080.N.16.6, temp = data$S32_2080D.N.19.11, brk = sim_brk_stay_2080[,j]))
  }
}

### Zero brook trout scenario
sim_brk_stay_2040_nobrk <- apply(sim_brk_stay_2040, 1, mean)
sim_brk_stay_2080_nobrk <- apply(sim_brk_stay_2080, 1, mean)

nobrk_2040 <- ifelse(sim_brk_stay_2040_nobrk > 5, 5, sim_brk_stay_2040_nobrk) # If brook trout, drop down to 5
nobrk_2080 <- ifelse(sim_brk_stay_2080_nobrk > 5, 5, sim_brk_stay_2080_nobrk)

but_nobrk_2040 <- c()
but_nobrk_2080 <- c()

for(i in 1:length(nobrk_2040)){
  but_nobrk_2040[i] <- pred_fun(flow = data$MS_2040.N.16.6, temp = data$S30_2040D.N.19.11, brk = nobrk_2040)
  but_nobrk_2080[i] <- pred_fun(flow = data$MS_2080.N.16.6, temp = data$S32_2080D.N.19.11, brk = nobrk_2080)
}

### Save and load
save(sim_but_hist, file = "sim_but_hist.Rda")
save(sim_but_2040, file = "sim_but_2040.Rda")
save(sim_but_2080, file = "sim_but_2080.Rda")

load("sim_but_hist.RDA")
load("sim_but_2040.RDA")
load("sim_but_2080.RDA")

### Summarize these
sim_but_hist_mean <- c()
for(i in 1:nrow(data)){
  sim_but_hist_mean[i] <- antilogit(pred_fun(flow = data$MS_Hist, temp = data$S1_93_11, brk = brk_stay_hist))
}


data$sim_but_hist_mean <- sim_but_hist_mean
data$sim_but_2040_mean <- apply(sim_but_2040, 1, mean)
data$sim_but_2040_min <-  apply(sim_but_2040, 1, function(x){quantile(x, 0.05)})
data$sim_but_2040_max <-  apply(sim_but_2040, 1, function(x){quantile(x, 0.95)})

data$sim_but_2080_mean <- apply(sim_but_2080, 1, mean)
data$sim_but_2080_min <-  apply(sim_but_2080, 1, function(x){quantile(x, 0.05)})
data$sim_but_2080_max <-  apply(sim_but_2080, 1, function(x){quantile(x, 0.95)})

data$but_nobrk_2040 <- antilogit(but_nobrk_2040)
data$but_nobrk_2080 <- antilogit(but_nobrk_2080)

### > 12C to zero
data$sim_but_hist_mean <- ifelse(data$S1_93_11 > 12, 0, data$sim_but_hist_mean)
data$sim_but_hist_min  <- ifelse(data$S1_93_11 > 12, 0, data$sim_but_hist_min)
data$sim_but_hist_max  <- ifelse(data$S1_93_11 > 12, 0, data$sim_but_hist_max)

data$sim_but_2040_mean <- ifelse(data$S30_2040D.N.19.11 > 12, 0, data$sim_but_2040_mean)
data$sim_but_2040_min  <- ifelse(data$S30_2040D.N.19.11 > 12, 0, data$sim_but_2040_min)
data$sim_but_2040_max  <- ifelse(data$S30_2040D.N.19.11 > 12, 0, data$sim_but_2040_max)

data$sim_but_2080_mean <- ifelse(data$S32_2080D.N.19.11 > 12, 0, data$sim_but_2080_mean)
data$sim_but_2080_min  <- ifelse(data$S32_2080D.N.19.11 > 12, 0, data$sim_but_2080_min)
data$sim_but_2080_max  <- ifelse(data$S32_2080D.N.19.11 > 12, 0, data$sim_but_2080_max)

data$but_nobrk_2040 <- ifelse(data$S30_2040D.N.19.11 > 12, NA, data$but_nobrk_2040)
data$but_nobrk_2080 <- ifelse(data$S32_2080D.N.19.11 > 12, NA, data$but_nobrk_2080)

### Plotting
data$seg_col <- brewer.pal(10, "RdYlBu")[cut(data$sim_but_hist_mean, 
                                            breaks = seq(from = 0, to = 1, by = 0.1))]
data$pch_2040 <- ifelse(data$sim_but_2040_mean == 0, NA, 21) # allow jittering of zeros
data$pch_2040b <- ifelse(data$sim_but_2040_mean == 0, 21, NA)
data$pch_2080 <- ifelse(data$sim_but_2080_mean == 0, NA, 21)
data$pch_2080b <- ifelse(data$sim_but_2080_mean == 0, 21, NA)

sort_d <- data[order(data$sim_but_hist_mean),] 
par(mfrow=c(1,2), mar = c(4,4,1,1), las = 1, lend = 2)
plot(sort_d$sim_but_hist_mean, pch = ".",
     xlab = "Site", ylab = "Pr(bull trout)")
segments(1:630,
         sort_d$sim_but_2040_min,
         1:630,
         sort_d$sim_but_2040_max,
         col = "darkgrey")
segments(1:630,
         sort_d$sim_but_2040_mean,
         1:630,
         sort_d$but_nobrk_2040, col = "red")
points(sort_d$sim_but_2040_mean, pch = sort_d$pch_2040, 
       col = "black", bg = sort_d$seg_col)
points(jitter(sort_d$sim_but_2040_mean, amount = 0.01), pch = sort_d$pch_2040b, 
       col = "black", bg = sort_d$seg_col)
lines(sort_d$sim_but_hist_mean, lwd = 2)
box(lwd = 2)
text(300, 0.9, "2040", cex = 1.5)

legend("topleft", bty = "n", lwd = 4, col = rev(brewer.pal(9, "RdYlBu")),
       legend = c(rev(seq(from = 0.1, to = 0.9, by = 0.1))), adj = 0, inset = 0.01)

par(mar = c(4,1,1,4))
plot(sort_d$sim_but_hist_mean, pch = ".",
     xlab = "Site", ylab = "",
     axes = F)
axis(1)
segments(1:630,
         sort_d$sim_but_2080_min,
         1:630,
         sort_d$sim_but_2080_max,
         col = "darkgrey")
segments(1:630,
         sort_d$sim_but_2080_mean,
         1:630,
         sort_d$but_nobrk_2080, col = "red")
points(sort_d$sim_but_2080_mean, pch = sort_d$pch_2080,
       col = "black", bg = sort_d$seg_col)
points(jitter(sort_d$sim_but_2080_mean, amount = 0.01), pch = sort_d$pch_2080b, 
       col = "black", bg = sort_d$seg_col)
lines(sort_d$sim_but_hist_mean, lwd = 2)
box(lwd = 2)
text(300, 0.9, "2080", cex = 1.5)

### Scatter plot
par(mar = c(4,4,1,1))
plot(log(MS_Hist) ~  S1_93_11, sort_d,                                          # plotting
     col = sort_d$seg_col, pch=19,
     xlim = c(5.06, 13.22),                                                     # axis extents consistent with kernal figures
     ylim = c(-1.46, log(99.20)),
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "Discharge (CFS)",
     cex = 1, axes = F)
axis(1, at = c(6,8,10,12))
axis(2, at = log(c(1,10,50,100)), labels = c(1,10,50,100))

segments(x0 = sort_d$S1_93_11, 
         y0 = log(sort_d$MS_Hist),
         x1 = sort_d$S30_2040D.N.19.11, 
         y1 = log(sort_d$MS_2040.N.16.6),
         col = sort_d$seg_col, lwd=1)
abline(v = 12, lty =2)
text(6, 4.2, "2040", cex = 1.5)
box(lwd = 2)

par(mar = c(4,1,1,4))
plot(log(MS_Hist) ~  S1_93_11, sort_d,                                          # plotting
     col = sort_d$seg_col, pch=19,
     xlim = c(5.06, 13.22),                                                     # axis extents consistent with kernal figures
     ylim = c(-1.46, log(99.20)),
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "Discharge (CFS)",
     cex = 1, axes = F)
axis(1, at = c(6,8,10,12))

segments(x0 = sort_d$S1_93_11, 
         y0 = log(sort_d$MS_Hist),
         x1 = sort_d$S32_2080D.N.19.11, 
         y1 = log(sort_d$MS_2080.N.16.6),
         col = sort_d$seg_col, lwd=1)
abline(v = 12, lty =2)
text(6, 4.2, "2080", cex = 1.5)
box(lwd = 2)

### SOME MODEL ASSESSMENT #######################################################
### modified confusion matrix and anitlogit and AUC functions
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }

confmat <- function(fits, dataset, threshold){ 
  fitted <- fits
  pred <- ifelse(fitted < threshold, 0, 1)
  obs_pred <- cbind(dataset$BUT1_P, pred)
  TP <- nrow(subset(obs_pred, obs_pred[,1]==1 & obs_pred[,2]==1))
  TN <- nrow(subset(obs_pred, obs_pred[,1]==0 & obs_pred[,2]==0))
  FP <- nrow(subset(obs_pred, obs_pred[,1]==0 & obs_pred[,2]==1))
  FN <- nrow(subset(obs_pred, obs_pred[,1]==1 & obs_pred[,2]==0))
  accuracy <- (TP+TN)/nrow(obs_pred)
  print(accuracy)
}

auc_ <- function(fits, dataset){
  fitted <- fits                                                      
  AUC <- roc(dataset$BUT1_P, fitted, auc=T)                                     
  AUC_CI <- ci(AUC, of="auc")                                                  
  print(AUC$auc)
  #print(AUC_CI)
}

confmat(data$sim_but_hist_mean, data, 0.5) # accuracy = 0.759
auc_(fits = data$sim_but_hist_mean, dataset = data) # AUC = 0.801

### Some summary statistics
summary(data$S30_2040D.N.19.11 - data$S1_93_11)
summary(data$S32_2080D.N.19.11 - data$S1_93_11)

summary(data$MS_2040.N.16.6 - data$MS_Hist)

summary(exp(data$brk_stay_2040) - exp(data$brk_stay_hist))
summary(exp(data$brk_stay_2080) - exp(data$brk_stay_hist))

summary(data$sim_but_2040_mean - data$sim_but_hist_mean)
summary(data$sim_but_2080_mean - data$sim_but_hist_mean)

sum(data$sim_but_2040_mean - data$sim_but_hist_mean <= -0.25)/630
sum(data$sim_but_2080_mean - data$sim_but_hist_mean <= -0.25)/630

sum(data$sim_but_2040_mean - data$sim_but_hist_mean <= -0.25 &
      data$sim_but_2040_mean == 0)/sum(data$sim_but_2040_mean - data$sim_but_hist_mean <= -0.25)
sum(data$sim_but_2080_mean - data$sim_but_hist_mean <= -0.25 &
      data$sim_but_2080_mean == 0)/sum(data$sim_but_2080_mean - data$sim_but_hist_mean <= -0.25)

sum(data$sim_but_2040_mean > data$sim_but_hist_mean)/630
sum(data$sim_but_2080_mean > data$sim_but_hist_mean)/630

summary(data$sim_but_2040_max - data$sim_but_2040_min)
summary(data$sim_but_2080_max - data$sim_but_2080_min)

sum(data$sim_but_2040_max < data$sim_but_hist_mean &
      data$S30_2040D.N.19.11 <=12)/630
sum(data$sim_but_2080_max < data$sim_but_hist_mean &
      data$S32_2080D.N.19.11 <=12)/630

sum(data$sim_but_2040_max < data$sim_but_hist_mean)/630
sum(data$sim_but_2080_max < data$sim_but_hist_mean)/630

sum(data$sim_but_2040_min > data$sim_but_hist_mean)/630
sum(data$sim_but_2080_min > data$sim_but_hist_mean)/630

subset(data, data$sim_but_2040_min > data$sim_but_hist_mean)$BRK2_Quant
subset(data, data$sim_but_2040_min > data$sim_but_hist_mean)$brk_stay_2040

sum(data$S30_2040D.N.19.11 > 12 & data$MS_Hist > 10)/sum(data$MS_Hist > 10)
sum(data$S32_2080D.N.19.11 > 12 & data$MS_Hist > 10)/sum(data$MS_Hist > 10)

sum(data$sim_but_2040_max > data$sim_but_hist_mean)
sum(data$sim_but_2080_max > data$sim_but_hist_mean)/630

sim_brk_stay_hist > sim_brk_stay_2040

### Comparing no brook versus increased brook trout
hist(sort_d$but_nobrk_2040 - sort_d$sim_but_2040_mean)
hist(sort_d$but_nobrk_2080 - sort_d$sim_but_2080_mean)

mean(sort_d$but_nobrk_2040 - sort_d$sim_but_2040_mean, na.rm = T) #4.08%
mean(sort_d$but_nobrk_2080 - sort_d$sim_but_2080_mean, na.rm = T) #3.99%

range(sort_d$but_nobrk_2040 - sort_d$sim_but_2040_mean, na.rm = T)
range(sort_d$but_nobrk_2080 - sort_d$sim_but_2080_mean, na.rm = T)

high_hist_occ <- subset(sort_d, sort_d$sim_but_hist_mean >= 0.75) # High historic probability of bull trout
mean(high_hist_occ$but_nobrk_2040 - high_hist_occ $sim_but_2040_mean, na.rm = T) #8.86%
mean(high_hist_occ$but_nobrk_2080 - high_hist_occ $sim_but_2080_mean, na.rm = T) #9.27%

plot(c(data$but_nobrk_2040 - data$sim_but_2040_mean) ~ 
       data$MS_Hist)
plot(c(data$but_nobrk_2040 - data$sim_but_2040_mean) ~ 
       data$S1_93_11)



-
mean((data$MS_Hist - data$MS_2040.N.16.6)/data$MS_Hist)
mean((data$MS_Hist - data$MS_2080.N.16.6)/data$MS_Hist)

mean((data$BRK2_Quant - data$MS_2080.N.16.6)/data$BRK2_Quant)