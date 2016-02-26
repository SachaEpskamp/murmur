
Bayes_murmur <- function(Y,X,ID,data,verbose=TRUE,prior = c("lmer_orth","identity","lmer_seq"),
                       type = c("full","sequential","orthogonal"),...){
  prior <- match.arg(prior)

  # Random effects prior:
  if (prior == "identity"){

    REprior <- diag(length(Y) + length(Y) * length(X))
  } else 
  {
    if (verbose)
    {
      message("computing prior")
    }
    if (prior == "lmer_seq"){
      
      lmerRes <- lmer_murmur(Y=Y,X=X,ID=ID,data=data,orthogonal = FALSE,verbose=verbose)
      REprior <- diag(diag(lmerRes$Sigma_randomEffects))
    } else if (prior == "lmer_orth"){
      lmerRes <- lmer_murmur(Y=Y,X=X,ID=ID,data=data,orthogonal = FALSE,verbose=verbose)
      REprior <- diag(diag(lmerRes$Sigma_randomEffects))
    } else stop("Prior method not implemented yet.")
    diag(REprior) <- pmax(diag(REprior),0.0001)
    
    # REprior <- solve(REprior)
  }
  
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
      R_Tau_randomEffects = REprior,
      # nRandom = ncol(Y) + ncol(Y) * ncol(X),
      randMean = rep(0, ncol(Y) + ncol(Y) * ncol(X))
    )
  
    ### RUN JAGS ###
    if (verbose){
      message("Estimating model")
    }

    samples <- jags(jagsData, parameters.to.save = 
                      c("beta0_fixed","Beta_fixed",
                        "beta0_randomEffects",
                        "Beta_randomEffects","Sigma_error",
                        "Sigma_randomEffects"),
                    model.file = switch(type,
                                full = JAGSModel_Bayes_full,
                                sequential = JAGSModel_Bayes_seq,
                                orthogonal = JAGSModel_Bayes_orth,
                                ),...)
    
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
    
    Results$Sigma_error <-  samples$BUGSoutput$sims.list$Sigma_error %>% 
      apply(2:3, mean)

    if (type == "full"){
      
      Results$Sigma_randomEffects <-  samples$BUGSoutput$sims.list$Sigma_randomEffects %>% 
        apply(2:3, mean)
      
    } else if (type == "sequential"){
      
      Results$Sigma_randomEffects <-  matrix(NA, ncol(Y) * (ncol(X)+1), ncol(Y) * (ncol(X) + 1))
    
      for (i in seq_len(ncol(Y))){
            Results$Sigma_randomEffects[ (i-1) * (ncol(X)+1) + (1:(ncol(X)+1)) , (i-1) * (ncol(X)+1) + (1:(ncol(X)+1)) ] <- 
            samples$BUGSoutput$sims.list$Sigma_randomEffects[,,,i] %>% 
              apply(2:3, mean)
      }  
      } else {
        
        Results$Sigma_randomEffects <-  matrix(NA, ncol(Y) * (ncol(X)+1), ncol(Y) * (ncol(X) + 1))
        
        for (i in seq_len(ncol(Y))){

          Results$Sigma_randomEffects[ (i-1) * (ncol(X)+1) + (1:(ncol(X)+1)) , (i-1) * (ncol(X)+1) + (1:(ncol(X)+1)) ] <- 
            samples$BUGSoutput$sims.list$Sigma_randomEffects[,i,] %>% 
            colSums %>% diag
        }  
        
      }
      
    Results$beta0_SD <- sqrt(diag(Results$Sigma_randomEffects)[((1:ncol(Y))-1)*(ncol(X)+1)+1])
    Results$Beta_SD <- matrix(sqrt(diag(Results$Sigma_randomEffects)[-(((1:ncol(Y))-1)*(ncol(X)+1)+1)]),ncol(Y),ncol(X),byrow=TRUE)
    
    return(Results)
    
  }