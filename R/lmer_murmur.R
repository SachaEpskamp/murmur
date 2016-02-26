lmer_murmur <- function(Y,X,ID,data,orthogonal = FALSE, verbose=TRUE,...){
  data <- data %>% select_(.dots = c(X,Y,ID)) %>% na.omit
  
  nY <- length(Y)
  nX <- length(X)
  nP <- length(unique(data[[ID]]))
  
  # Output list:
  lmerResults <- list()
 
  if (verbose){
    pb <- txtProgressBar(min = 0, max = length(Y), style = 3)
  }
  # start nodewise loop:
  for (i in seq_along(Y)){
    
    # Setup model:
    mod <- paste0(
      Y[i],
      " ~ ",
      paste(X,collapse="+"),
      " + (",
      paste(X,collapse="+"),
      ifelse(orthogonal,"||","|"),
      ID,
      ")"
    )
    
    # Formula:
    formula <- as.formula(mod)
    
    # Rune lmer:
    lmerResults[[i]] <- suppressWarnings(lmer(formula, data = data, ...))

    if (verbose){
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose){
    close(pb)
  }
  

  
  # Results:
  Results <- list(output = lmerResults)
  
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
  
  Results$beta0_SD <- rep(NA,nY)
  Results$Beta_SD <- matrix(NA, nY, nX)
  
  # Populate matrices:
  for (i in 1:nY){
    Results$Beta_fixed[i,] <- fixef(lmerResults[[i]])[-1]
    Results$beta0_fixed[i] <- fixef(lmerResults[[i]])[1]
    
    for (p in 1:nP){

      Results$Beta_random[[p]][i,] <- unlist(ranef(lmerResults[[i]])[[ID]][p,-1])
      Results$beta0_random[[p]][i] <- unlist(ranef(lmerResults[[i]])[[ID]][p,1]) 
    }
 
    Results$Sigma_error[i,i] <- sigma(lmerResults[[i]])^2
    
    if (orthogonal){
      vars <- VarCorr(lmerResults[[i]])
      mat <- diag(unlist(vars[grepl(ID,names(vars))]))
    } else {
      mat <- VarCorr(lmerResults[[i]])[[ID]][,]      
    }

    # Vector of relevant random effects:
    Results$Sigma_randomEffects[ (i-1) * (nX+1) + (1:(nX+1)) , (i-1) * (nX+1) + (1:(nX+1)) ] <- mat
    
    Results$beta0_SD[i] <- sqrt(diag(mat)[1])
    Results$Beta_SD[i,] <- sqrt(diag(mat)[-1])
  }
  
  
  return(Results)
  
}