#############################################
##                                         ##
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 July 2017                            ##
##                                         ##
## Hierarchical model fitting for bull tr  ##
##                                         ##
#############################################

### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### some required packages
require("lme4") # for model fitting
require("pROC") # AUC estimation

### some user-defined functions to use below ###
### AUC estimation
auc_lme4 <- function(model, dataset){
  fitted <- fitted(model)                                                       # fitted are back-transformed already!
  AUC <- roc(dataset$BUT1_P, fitted, auc=T)                                     # specific to these models/dataset
  AUC_CI <- ci(AUC, of="auc")                                                   # 95% CI
  print(AUC$auc)
  #print(AUC_CI)
}

### confusion matrix stats
confmat_lme4 <- function(model, dataset, threshold){ 
  fitted <- fitted(model)
  pred <- ifelse(fitted < threshold, 0, 1)
  obs_pred <- cbind(dataset$BUT1_P, pred)
  TP <- nrow(subset(obs_pred, obs_pred[,1]==1 & obs_pred[,2]==1))
  TN <- nrow(subset(obs_pred, obs_pred[,1]==0 & obs_pred[,2]==0))
  FP <- nrow(subset(obs_pred, obs_pred[,1]==0 & obs_pred[,2]==1))
  FN <- nrow(subset(obs_pred, obs_pred[,1]==1 & obs_pred[,2]==0))
  accuracy <- (TP+TN)/nrow(obs_pred)
  print(accuracy)
}

### brook trout P/A as candidate covariate models first - NOT included in final manuscript due to weaker model performance
### model formulas for different random effects with same "full" fixed effects structure
sat_PA_models <- list(
  M1 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 | Patch)",
  M2 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + S1_93_11 | Patch)",
  M3 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + log(MS_Hist) | Patch)",
  M4 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + BRK2_P | Patch)",
  M5 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + S1_93_11 + log(MS_Hist) | Patch)",
  M6 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + S1_93_11 + BRK2_P | Patch)",
  M7 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + log(MS_Hist) + BRK2_P | Patch)",
  M8 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + S1_93_11 + log(MS_Hist) + BRK2_P | Patch)")

### loop through list of model formulas and fit each
sat_PA_model_summaries <- matrix(, nrow=length(sat_PA_models), ncol=3)          # AIC, AUC, and accuracy summaries
sat_PA_model_fits <- list()                                                     # list of fitted models

for(i in 1:length(sat_PA_models)){
  the_model <- glmer(paste(sat_PA_models[i]),
                      data = data, family = "binomial")

  the_AIC <- AIC(the_model)
  the_AUC <- auc_lme4(the_model, data)[1:1]
  the_acc <- confmat_lme4(the_model, data, threshold=0.5)

  sat_PA_model_summaries[i,] <- c(the_AIC, the_AUC, the_acc)
  sat_PA_model_fits[i] <- the_model
}

rownames(sat_PA_model_summaries) <- sat_PA_models
colnames(sat_PA_model_summaries) <- c("AIC","AUC","Accuracy")
save(sat_PA_model_fits, file="sat_PA_model_fits.rda")
write.csv(sat_PA_model_summaries, "tables/sat_PA_model_summaries.csv")

### now selection on fixed effects (PA models)
PA_models <- list(
  M2 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + BRK2_P + BRK2_P*log(MS_Hist) + BRK2_P*S1_93_11 + (1 + log(MS_Hist) + BRK2_P | Patch)",
  M3 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*log(MS_Hist) + (1 + log(MS_Hist) + BRK2_P | Patch)",
  M4 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + BRK2_P + BRK2_P*log(MS_Hist) + (1 + log(MS_Hist) + BRK2_P | Patch)",
  M5 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + BRK2_P*S1_93_11 + (1 + log(MS_Hist) + BRK2_P | Patch)",
  M6 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + BRK2_P + BRK2_P*S1_93_11 + (1 + log(MS_Hist) + BRK2_P | Patch)",
  M7 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + BRK2_P + (1 + BRK2_P + log(MS_Hist) | Patch)",
  M8 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + BRK2_P + (1 + log(MS_Hist) + BRK2_P | Patch)")

### loop through list of model formulas and fit each
PA_model_summaries <- matrix(, nrow=length(PA_models), ncol=3)                  # AIC, AUC, and accuracy summaries
PA_model_fits <- list()                                                         # list of fitted models

for(i in 1:length(PA_models)){
  the_model <- glmer(paste(PA_models[i]),
                     data = data, family = "binomial")

  the_AIC <- AIC(the_model)
  the_AUC <- auc_lme4(the_model, data)[1:1]
  the_acc <- confmat_lme4(the_model, data, threshold=0.5)
  
  PA_model_summaries[i,] <- c(the_AIC, the_AUC, the_acc)
  PA_model_fits[i] <- the_model
}

rownames(PA_model_summaries) <- PA_models
colnames(PA_model_summaries) <- c("AIC","AUC","Accuracy")
save(PA_model_fits, file="PA_model_fits.rda")
write.csv(PA_model_summaries, "tables/PA_model_summaries.csv")

### same, but for quant models - reported in the main manuscript
### random effect selection
sat_Q_models <- c(
  M1 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 | Patch)",
  M2 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + S1_93_11 | Patch)",
  M3 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + log(MS_Hist) | Patch)",
  M4 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + log(BRK2_Quant+1) | Patch)",
  M5 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + S1_93_11 + log(MS_Hist) | Patch)",
  M6 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + S1_93_11 + log(BRK2_Quant+1) | Patch)",
  M7 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M8 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + log(BRK2_Quant+1)*S1_93_11 + I(log(BRK2_Quant+1)^2) + (1 + S1_93_11 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)")

### loop through list of model formulas and fit each
sat_Q_model_summaries <- matrix(, nrow=length(sat_Q_models), ncol=3)            # AIC, AUC, and accuracy summaries
sat_Q_model_fits <- list()                                                      # list of fitted models

for(i in 1:length(sat_Q_models)){
  the_model <- glmer(paste(sat_Q_models[i]),
                     data = data, family = "binomial")
  
  the_AIC <- AIC(the_model)
  the_AUC <- auc_lme4(the_model, data)[1:1]
  the_acc <- confmat_lme4(the_model, data, threshold=0.5)
  
  sat_Q_model_summaries[i,] <- c(the_AIC, the_AUC, the_acc)
  sat_Q_model_fits[i] <- the_model
}
rownames(sat_Q_model_summaries) <- sat_Q_models
colnames(sat_Q_model_summaries) <- c("AIC","AUC","Accuracy")
save(sat_Q_model_fits, file="sat_Q_model_fits.rda")
write.csv(sat_Q_model_summaries, "tables/sat_Q_model_summaries.csv")

### now select on fixed effects (quant models)
Q_models <- c(
  M1 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + log(BRK2_Quant+1)*S1_93_11 + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M2 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + log(BRK2_Quant+1)*S1_93_11 + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M3 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M4 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",

  M5 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + log(BRK2_Quant+1)*S1_93_11 + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M6 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + log(BRK2_Quant+1)*S1_93_11 + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M7 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M8 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + log(BRK2_Quant+1) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  
  M9 <- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + log(BRK2_Quant+1)*S1_93_11 + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M10<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + log(BRK2_Quant+1)*S1_93_11 + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M11<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M12<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + I(log(BRK2_Quant+1)^2) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  
  M13<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*S1_93_11 + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M14<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*S1_93_11 + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M15<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + log(BRK2_Quant+1)*log(MS_Hist) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)",
  M16<- "BUT1_P ~ log(MS_Hist) + S1_93_11 + SLOPE + log(BRK2_Quant+1) + (1 + log(MS_Hist) + log(BRK2_Quant+1) | Patch)")

### loop through list of model formulas and fit each
Q_model_summaries <- matrix(, nrow=length(Q_models), ncol=3)                    # AIC, AUC, and accuracy summaries
Q_model_fits <- list()                                                          # list of fitted models

for(i in 1:length(Q_models)){
  the_model <- glmer(paste(Q_models[i]),
                     data = data, family = "binomial")
  
  the_AIC <- AIC(the_model)
  the_AUC <- auc_lme4(the_model, data)[1:1]
  the_acc <- confmat_lme4(the_model, data, threshold=0.5)
  
  Q_model_summaries[i,] <- c(the_AIC, the_AUC, the_acc)
  Q_model_fits[i] <- the_model
}
rownames(Q_model_summaries) <- Q_models
colnames(Q_model_summaries) <- c("AIC","AUC","Accuracy")
save(Q_model_fits, file="Q_model_fits.rda")
write.csv(Q_model_summaries, "tables/Q_model_summaries.csv")

