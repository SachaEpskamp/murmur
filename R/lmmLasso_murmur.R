lmmlasso_murmur <- function(Y,X,ID,data,orthogonal = FALSE, verbose=TRUE,...){
  data <- data %>% select_(.dots = c(X,Y,ID)) %>% na.omit
  
  nY <- length(Y)
  nX <- length(X)
  nP <- length(unique(data[[ID]]))
  Y <- as.matrix(data[,Y])
  X <- as.matrix(data[,X])
  ID <- data[[ID]]
  
  # Add intercept:
  Xint <- cbind(1,X)
  
  # Output list:
  lmmmlassoResults <- list()
 
  if (verbose){
    pb <- txtProgressBar(min = 0, max = length(Y), style = 3)
  }
  # start nodewise loop:
  for (i in 1:ncol(Y)){
    
    # Rune lmer:
    lmmmlassoResults[[i]] <- BIClmmlasso(Xint, Y[,i], Xint, ID, pdMat = ifelse(orthogonal,"pdDiag","pdSym"))
    
    if (verbose){
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose){
    close(pb)
  }

  
  # Results:
  Results <- list(output = lmmmlassoResults)
  
  # Setup matrices to fill:
  Results$Beta_fixed <-  matrix(NA, nY, nX)
  Results$Beta_random <-  lapply(1:nP,
                                 function(i){
                                   matrix(NA, nY, nX)
                                 })
  Results$beta0_fixed <-  rep(NA, nY)
  Results$beta0_random <- lapply(1:nP,function(i){
                                 rep(NA, nY)
                                 })
  
  Results$Sigma_error <-  matrix(NA, nY, nY)

  Results$Sigma_randomEffects <-  matrix(NA, nY * (nX+1), nY * (nX + 1))
  
  # Populate matrices:
  for (i in 1:nY){
    
    Results$Beta_fixed[i,] <- coef(lmmmlassoResults[[i]])[-1]
    Results$beta0_fixed[i] <- coef(lmmmlassoResults[[i]])[1]
    
    randMat <- matrix(lmmmlassoResults[[i]]$random,nP,nX+1,byrow=TRUE)
    
    for (p in 1:nP){
      Results$Beta_random[[p]][i,] <- randMat[p,-1]
      Results$beta0_random[[p]][i] <- randMat[p,1]
    }
    
    Results$Sigma_error[i,i] <- lmmmlassoResults[[i]]$sigma^2

    Results$Sigma_randomEffects[ (i-1) * (nX+1) + (1:(nX+1)) , (i-1) * (nX+1) + (1:(nX+1)) ] <- lmmmlassoResults[[i]]$Psi
  }
  
  
  return(Results)
  
}