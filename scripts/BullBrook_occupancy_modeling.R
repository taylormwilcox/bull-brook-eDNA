### bring in data and setup some fields
data <- read.csv("FlintRock_UClarkFork_Bull_Brook.csv")                         # Raw data
data$BUT1_P <- ifelse(data$BUT1_Wells > 0, 1, 0)                                # P/A bull trout eDNA
data$BRK2_P <- ifelse(data$BRK2_Wells > 0, 1, 0)                                # P/A brook trout eDNA
data$Patch <- as.factor(data$Patch)                                             # make sure "Patch" is treated as a factor

### packages
require("unmarked")
require("tidyverse")

### for back transformations
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }

### munge data (bull)
bull <- data$BUT1_Wells %>% 
  sapply(., function(x) c(rep(1,x), rep(0,3-x))) %>%
  t(.) %>%
  unmarkedFrameOccu(., siteCovs=data[,]) 

### fit a simple occupancy model (bull)
occ.1 <- occu(~ 1 ~ 1, data=bull)
summary(occ.1) # simple model
antilogit(2.52) # 0.926 per PCR replicate
1 - (1-antilogit(2.52))^3 # per triplicate PCR

#### munge data (brook)
brook <- data$BRK2_Wells %>% 
  sapply(., function(x) c(rep(1,x), rep(0,3-x))) %>%
  t(.) %>%
  unmarkedFrameOccu(., siteCovs=data[,]) 

### fit a simple occupancy model (brook)
occ.2 <- occu(~ 1 ~ 1, data=brook)
summary(occ.2) # simple model
antilogit(2.8) # 0.943 per PCR replicate
1 - (1-antilogit(2.8))^3 # per triplicate PCR



