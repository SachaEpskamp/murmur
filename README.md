# MUltivaRiate MUlti-level Regression (murmur)

Example:

```r
# Install package:
library("devtools")
install_github("sachaepskamp/murmur")

# Load package:
library("murmur")

### Setup:
nVar <- 3
nPerson <- 25
nTime <- 25

### Simulate model:
Mod <- mlVARsim(nPerson, nVar, nTime)
Data <- Mod$Data

# Setting no X runs VAR:
Res1 <- murmur(Y = c("V1","V2","V3"),ID = "ID", data = Data, 
               n.chain = 1, n.iter = 500, n.burnin = 100)

# Compare fixed effects:
cor(c(t(Res1$Beta_fixed)),Mod$parameters$fixed)

# Specifying L1_.. in X uses lagged indicators (here lag1 and 2 for V1 and V2, V3 is a non-lagged co-variate):
Res2 <- murmur(Y = c("V1","V2"),X = c("L1_V1","L1_V2","L2_V1","L2_V3","V3"),
               ID = "ID", data = Data, 
               n.chain = 1, n.iter = 500, n.burnin = 100)

# Fixed effects:
Res2$Beta_fixed
```