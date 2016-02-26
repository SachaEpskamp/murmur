### JAGS MODEL ###
JAGSModel_Bayes_seq <- function(){
  # Y
  # X
  # nX = number of X variables
  # nY = number of Y variables
  # nT = number of time points
  # id = id vector
  # nP = number of persons
  
  ### Identity matrices:
  ### Priors: ###
  
  ## Fixed effects:
  for (i in 1:nY){
    beta0_fixed[i] ~ dnorm(0,0.001)
    for (j in 1:nX){
      Beta_fixed[i,j] ~ dnorm(0,0.001)
    }
  }
  
  ## Random effects:
  # Vector of random effects
  
  
  # Standard deviations of error:
  for (i in 1:nY){
    Tau_error[i,i] ~ dgamma(0.001,0.001)
    Sigma_error[i,i] <- 1/Tau_error[i,i]
  }
  
  for (i in 1:nY){
    Tau_randomEffects[1:(nX+1),1:(nX+1),i] ~ dwish(R_Tau_randomEffects[((i-1)*(nX+1))+(1:(nX+1)),((i-1)*(nX+1))+(1:(nX+1))], (nX+1)) 
    Sigma_randomEffects[1:(nX+1),1:(nX+1),i] <- inverse(Tau_randomEffects[,,i])
    
    for (p in 1:nP){
      Randoms[1:(nX+1),i,p] ~ dmnorm(randMean[1:(nX+1)], Tau_randomEffects[,,i])
      beta0_randomEffects[i,p] <- Randoms[1,i,p]
      Beta_randomEffects[i,1:nX,p]  <- Randoms[1+1:nX,i,p]
    }
  }
  
  # likelihood:
  for (t in 1:nT){
    for (y in 1:nY){
      Y[t,y] ~ dnorm(beta0_fixed[y] + beta0_randomEffects[y,id[t]] + 
                       (Beta_fixed[y,1:nX] +  Beta_randomEffects[y,1:nX,id[t]]) %*% X[t,1:nX], Tau_error[y,y])      
    }
  }
}

### JAGS MODEL ###
JAGSModel_Bayes_orth <- function(){
  # Y
  # X
  # nX = number of X variables
  # nY = number of Y variables
  # nT = number of time points
  # id = id vector
  # nP = number of persons
  
  ### Identity matrices:
  ### Priors: ###
  
  ## Fixed effects:
  for (i in 1:nY){
    beta0_fixed[i] ~ dnorm(0,0.001)
    for (j in 1:nX){
      Beta_fixed[i,j] ~ dnorm(0,0.001)
    }
  }
  
  ## Random effects:
  # Vector of random effects
  
  
  # Standard deviations of error:
  for (i in 1:nY){
    Tau_error[i,i] ~ dgamma(0.001,0.001)
    Sigma_error[i,i] <- 1/Tau_error[i,i]
  }
  
  for (i in 1:nY){
    for (j in 1:(nX+1)){
      Tau_randomEffects[i,j] ~ dgamma(1/2, R_Tau_randomEffects[((i-1)*(nX+1))+j,((i-1)*(nX+1))+j]) 
      # Tau_randomEffects[j,j,i] ~ dgamma(0.00001,0.00001) 
      Sigma_randomEffects[i,j] <- 1/Tau_randomEffects[i,j]
    }
    
    for (p in 1:nP){
      for (j in 1:(nX+1)){
        Randoms[i,j,p] ~ dnorm(0, Tau_randomEffects[i,j]) 
      }
      beta0_randomEffects[i,p] <- Randoms[i,1,p]
      Beta_randomEffects[i,1:nX,p]  <- Randoms[i,1+(1:nX),p]
    }
  }
  
#   for (i in 1:nY){
#     for (j in 1:(nX+1)){
#       Tau_randomEffects[j,j,i] ~ dgamma(0.001,0.001)
#       Sigma_randomEffects[j,j,i] <- 1/Tau_randomEffects[j,j,i]
#     }
#     
#     
#     for (p in 1:nP){
#       for (j in 1:(nX+1)){
#         Randoms[j,i,p] ~ dnorm(0, Tau_randomEffects[j,j,i])
#       }
#       beta0_randomEffects[i,p] <- Randoms[1,i,p]
#       Beta_randomEffects[i,1:nX,p]  <- Randoms[1+1:nX,i,p]
#     }
#   }
  
  # likelihood:
  for (t in 1:nT){
    for (y in 1:nY){
      Y[t,y] ~ dnorm(beta0_fixed[y] + beta0_randomEffects[y,id[t]] + 
                       (Beta_fixed[y,1:nX] +  Beta_randomEffects[y,1:nX,id[t]]) %*% X[t,1:nX], Tau_error[y,y])      
    }
  }
}

### JAGS MODEL ###
JAGSModel_Bayes_full <- function(){
  # Y
  # X
  # nX = number of X variables
  # nY = number of Y variables
  # nT = number of time points
  # id = id vector
  # nP = number of persons
  
  ### Identity matrices:
  ### Priors: ###
  
  ## Fixed effects:
  for (i in 1:nY){
    beta0_fixed[i] ~ dnorm(0,0.001)
    for (j in 1:nX){
      Beta_fixed[i,j] ~ dnorm(0,0.001)
    }
  }
  
  ## Random effects:
  # Vector of random effects
  Tau_randomEffects[1:((nX+1)*nY),1:((nX+1)*nY)] ~ dwish(R_Tau_randomEffects, ((nX+1)*nY)) 
  Sigma_randomEffects[1:((nX+1)*nY),1:((nX+1)*nY)] <- inverse(Tau_randomEffects)  
  
  # General error structure:
  Tau_error[1:nY,1:nY] ~ dwish(R_Tau_error[1:nY,1:nY], nY)
  Sigma_error[1:nY,1:nY] <- inverse(Tau_error[1:nY,1:nY])
  
  for (p in 1:nP){
    Randoms[1:((nX+1)*nY),p] ~ dmnorm(randMean, Tau_randomEffects)
    
    for (i in 1:nY){
      beta0_randomEffects[i,p] <- Randoms[(i-1)*(nX+1)+1,p]
      for (j in 1:nX){
        Beta_randomEffects[i,j,p]  <- Randoms[(i-1)*(nX+1)+1+j,p]
      }
    }
  }
  
  # likelihood:
  for (t in 1:nT){
    Y[t,1:nY] ~ dmnorm(beta0_fixed[1:nY] + beta0_randomEffects[1:nY,id[t]] + 
                         (Beta_fixed[1:nY,1:nX] +  Beta_randomEffects[1:nY,1:nX,id[t]]) %*% X[t,1:nX], Tau_error[1:nY,1:nY])
  }
}
