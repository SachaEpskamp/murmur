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
    beta0_fixed[i] ~ dnorm(0,1)
    for (j in 1:nX){
      Beta_fixed[i,j] ~ dnorm(0,1)
    }
  }
  
  ## Random effects:
  # Vector of random effects
  Tau_randomEffects[1:nRandom,1:nRandom] ~ dwish(R_Tau_randomEffects, nRandom) 
  Sigma_randomEffects[1:nRandom,1:nRandom] <- inverse(Tau_randomEffects)  
  
  # General error structure:
  Tau_error[1:nY,1:nY] ~ dwish(R_Tau_error[1:nY,1:nY], nY)
  Sigma_error[1:nY,1:nY] <- inverse(Tau_error[1:nY,1:nY])
  
  for (p in 1:nP){
    Randoms[1:nRandom,p] ~ dmnorm(randMean, Tau_randomEffects)

    
    for (i in 1:nY){
      beta0_randomEffects[i,p] <- Randoms[i,p]
      for (j in 1:nX){
        Beta_randomEffects[i,j,p]  <- Randoms[(nY + (i-1)*nX + j ),p]
      }
    }
  }
  
  # likelihood:
  for (t in 1:nT){
    Y[t,1:nY] ~ dmnorm(beta0_fixed[1:nY] + beta0_randomEffects[1:nY,id[t]] + 
                         (Beta_fixed[1:nY,1:nX] +  Beta_randomEffects[1:nY,1:nX,id[t]]) %*% X[t,1:nX], Tau_error[1:nY,1:nY])
  }
}

Bayes_full <- function(Y,X,ID,data,...){
  data <- data %>% select_(.dots = c(X,Y,ID)) %>% na.omit
  
  Y <- as.matrix(data[,Y])
  X <- as.matrix(data[,X])
  ID <- data[[ID]]

  # obtain data:
  jagsData <- list(
    Y=Y,
    X=X,
    nX = ncol(X), 
    nY = ncol(Y), 
    nT = nrow(Y), 
    id = ID, 
    nP = max(ID),
    R_Tau_error = diag(ncol(Y)),
    R_Tau_randomEffects = diag(ncol(Y) + ncol(Y) * ncol(X)),
    nRandom = ncol(Y) + ncol(Y) * ncol(X),
    randMean = rep(0, ncol(Y) + ncol(Y) * ncol(X))
  )
  
  ### RUN JAGS ###
  samples <- jags(jagsData, parameters.to.save = 
                    c("beta0_fixed","Beta_fixed",
                  "beta0_randomEffects",
                  "Beta_randomEffects","Sigma_error",
                  "Sigma_randomEffects"),
                  model.file ="mlRegression.txt", ...)
  
  Results <- list(output = samples)
  
  ### Fixed effects:
  Results$Beta_fixed <-  samples$BUGSoutput$sims.list$Beta_fixed %>% 
    apply(2:3, mean)
  Results$Beta_random <-  lapply(1:jagsData$nP,
                                 function(i){
                                   samples$BUGSoutput$sims.list$Beta_randomEffects[,,,i] %>% 
                                     apply(2:3, mean)
                                 })
  Results$beta0_fixed <-  samples$BUGSoutput$sims.list$beta0_fixed %>% 
    colMeans
  Results$beta0_random <-  lapply(1:jagsData$nP,
                                  function(i){
                                    samples$BUGSoutput$sims.list$beta0_randomEffects[,,i] %>% 
                                      colMeans
                                  })
  
  Results$Sigma_error <-  lapply(1:jagsData$nP,
                                 function(i){
                                   samples$BUGSoutput$sims.list$Sigma_error[,,,i] %>% 
                                     apply(2:3, mean)
                                 })
  Results$Sigma_randomEffects <-  samples$BUGSoutput$sims.list$Sigma_randomEffects %>% 
    apply(2:3, mean)
  
  
  return(Results)
  
}