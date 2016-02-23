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
Res1 <- murmur(Y = c("V1","V2","V3"),ID = "ID", data = Data)

# Compare fixed effects:
cor(c(t(Res1$Beta_fixed)),Mod$parameters$fixed)

# Specifying L1_.. in X uses lagged indicators:
Res2 <- murmur(Y = c("V1","V2"),X = c("L1_V1","L1_V2","L1_V3","L2_V1","L2_V3"),ID = "ID", data = Data)

# Fixed effects:
Res2$Beta_fixed


```