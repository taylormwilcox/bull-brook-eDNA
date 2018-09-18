#############################################
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 21 June 2018                            ##
##                                         ##
## Concordance with electrofishing         ##       
##                                         ##
#############################################

### bring in data 
data <- read.csv(file.choose("UpperClarkFork_eDNA-electrofishing-concordance.csv"))

par(mfrow = c(1,1))
par(mar = c(5, 5, 0.5, 0.5), las = 1, cex = 1.2)
plot(log(BRK2_Quant + 1) ~ log(BRK + 1), data,
     xlab = "log(Brook trout captured + 1)", ylab = "log(Brook trout mtDNA copies/rxn + 1)",
     pch = 19, cex = 1, axes = F, cex.lab = 1.2)

axis(1, at = log(c(1, 6, 11, 51, 101)), labels = c(0, 5, 10, 50, 100))
axis(2, at = log(c(1, 11, 101, 201, 501)), labels = c(0, 10, 100, 200, 500))
box(lwd = 2)

cor.test(data$BRK2_Quant, data$BRK)
cor.test(log(data$BRK2_Quant + 1), log(data$BRK + 1)) # r = 0.759, P < 0.001
summary(lm(log(BRK2_Quant + 1) ~ log(BRK + 1), data)) # Rsq = 0.576, P < 0.001
abline(0.4649, 0.8818, lwd = 2)

### Bull trout eDNA v electrofishing
eDNA_BULL <- ifelse(data$BUT1_Wells > 0, 1, 0)
elec_BULL <- ifelse(data$BULL > 0 | data$BULL_BRK > 0, 1, 0)
sum(eDNA_BULL == elec_BULL)/length(eDNA_BULL) # 91.7% concordance
sum(eDNA_BULL == 1 & elec_BULL == 1) #TP
sum(eDNA_BULL == 0 & elec_BULL == 0) #TN
sum(eDNA_BULL == 1 & elec_BULL == 0) #FP - 4 eDNA detections where no fish with electrofishing
sum(eDNA_BULL == 0 & elec_BULL == 1) #FN


data[c(19, 34, 42, 43), ] # exploring mismatches
# Little Blackfoot = repeated temp/spatial detections
# MF Warm Springs = bull trout were detected in Warm Springs and WF Warm Springs < 3 km away (nearest sampled sites)
# Twin Lakes = bull trout detected downstream and there is a known population described in the report upstream

### Brook trout eDNA v electrofishing
eDNA_BRK <- ifelse(data$BRK2_Quant > 0, 1, 0)
elec_BRK <- ifelse(data$BRK > 0, 1, 0)
sum(eDNA_BRK == elec_BRK)/length(eDNA_BRK) # 93.8% concordance
sum(eDNA_BRK == 1 & elec_BRK == 1) #TP
sum(eDNA_BRK == 0 & elec_BRK == 0) #TN
sum(eDNA_BRK == 1 & elec_BRK == 0) #FP - 3 eDNA without electrofishing
sum(eDNA_BRK == 0 & elec_BRK == 1) #FN 

data[c(5, 11, 12), ] # exploring mismatches
# Elk Creek = brook trout electrofished < 2 km downstream (nearest site)
# Boulder = brook trout electrofishd both upstream and downstream < 4 km (nearest sites)

