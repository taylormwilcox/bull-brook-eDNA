#############################################
##                                         ##
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 July 2017                            ##
##                                         ##
## Data summaries                          ##
##                                         ##
#############################################

### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")                         # Raw data
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### covariate summaries table
data_cov <- data.frame(data$S1_93_11, data$MS_Hist, data$SLOPE, data$BRK2_Quant)
data_cov_summaries <- matrix(, nrow=ncol(data_cov), ncol=4)
for(i in 1:ncol(data_cov)){
  data_cov_summaries[i,] <-
    c(min(data_cov[,i]),max(data_cov[,i]),mean(data_cov[,i]),median(data_cov[,i]))
}
rownames(data_cov_summaries) <- c("Temp","Summer flow","Slope","BRK eDNA copies/rxn")
colnames(data_cov_summaries) <- c("Min","Max","Mean","Median")
data_cov_summaries <- round(data_cov_summaries,3)
data_cov_summaries
write.csv(data_cov_summaries, "tables/cov_summaries.csv")

### covariate correlation matrix table
data_cov2 <- data.frame(data$S1_93_11, data$MS_Hist, data$SLOPE, data$BRK2_Quant, data$BRK2_P, data$BUT1_P)
cor_matrix <- round(cor(data_cov2),3)
rownames(cor_matrix) <- c("Temp","Summer flow","Slope","BRK copies", "BRK eDNA", "BUT eDNA")
colnames(cor_matrix) <- c("Temp","Summer flow","Slope","BRK copies", "BRK eDNA", "BUT eDNA")
write.csv(cor_matrix, "tables/cor_matrix.csv")

### detection types
BULL <- subset(data, data$BUT1_P == 1)                                          # bull trout detected
BULLnoBRK <- subset(data, data$BUT1_P == 1 & data$BRK2_P == 0)                  # bull no brook trout
BULLandBRK <- subset(data, data$BUT1_P == 1 & data$BRK2_P == 1)                 # both species
BRK <- subset(data, data$BRK2_P == 1)                                           # brook trout detected
BRKnoBULL <- subset(data, data$BUT1_P == 0 & data$BRK2_P == 1)                  # brook and no bull
NO <- subset(data, data$BUT1_P == 0 & data$BRK2_P == 0)                         # neither species
BULLorBRK <- subset(data, data$BUT1_P == 1 | data$BRK2_P == 1)                  # either species
HI_BRK <- subset(data, data$BRK2_Quant >= 100)                                  # high brook trout eDNA concentration

### how many detections of each kind
nrow(BULL) # 262
nrow(BRK) # 293
nrow(BULLandBRK) # 121
nrow(BRK) - nrow(HI_BRK) # 222 of brook trout detections < 100 copies/rxn

### more summaries
median(BULL$MS_Hist) # 9.32 
nrow(subset(BULL, BULL$MS_Hist < 5))/nrow(BULL) # 27.9
nrow(subset(data, data$MS_Hist < 5))/nrow(data) # 57.0
median(BULLandBRK$MS_Hist) # 15.977
nrow(subset(BRK, BRK$MS_Hist < 5))/nrow(BRK) # 49.8
median(BRKnoBULL$MS_Hist) # 2.59
median(HI_BRK$MS_Hist) # 4.51
median(subset(BRK, BRK$BRK2_Quant < 100)$MS_Hist) #5.47
min(BULLorBRK$S1_93_11) # min temp Salvelinus detection 6.69 C
summary(subset(data, data$BRK2_Quant > 100)$MS_Hist)
summary(subset(data, data$BRK2_Quant > 500)$MS_Hist)
