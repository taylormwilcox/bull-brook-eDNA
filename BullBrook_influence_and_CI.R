#############################################
##                                         ##
## Bull trout/brook trout interactions     ##
## Taylor Wilcox                           ##
## 20 July 2017                            ##
##                                         ##
## Patch influence and fixed effect CIs    ##
##                                         ##
#############################################

### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### some packages
require('lme4') 
require("pROC")
require("influence.ME") # to estimate Patch influence

### model fits from before
load("Q_model_fits.rda")

### assign models and create influence objects
Q_model <- Q_model_fits[[3]]                                                    # top model
Q_influence <- influence(Q_model, group="Patch")                                # hold each patch out and re-fit the model with same form (COMPUTATIONALLY INTENSIVE)
save(Q_influence, file="Q_influence.rda")                                       # save object

### cooks distance
load("Q_influence.rda") # load in the influence object that was fit
cooks <- cooks.distance.estex(Q_influence)
patch_counts <- aggregate(data$Patch, by=list(data$Patch), FUN=length)

### cooks distance versus sites/patch (figure for supplement)
par(mfrow=c(1,1))
plot(cooks ~ patch_counts$x, pch=".",
     xlab="sites/patch", ylab="Cooks distance",
     axes=F)
text(cooks ~ patch_counts$x, labels=levels(data$Patch))
axis(1)
axis(2)
box()

### bootstrapping confidence intervals for betas
Q_model_ci <- confint.merMod(Q_model, 
                     level=0.95,
                     method="boot",
                     nsim=200)
save(Q_model_ci, file="Q_model_ci.rda") # save
Q_model_ci <- load("Q_model_ci.rda")
Q_model_ci
