BIClmmlasso <- function(x,y,z,grp,..., startLambda = 0.1, lambdaIncrease = 5){
  library("lmmlasso")
  
  Results <- list()
  lambdaCur <- startLambda
  BICs <- numeric(0)
  lambdas <- numeric(0)
  nVar <- numeric(0)
  c <- 1
  repeat{
    lambdas <- c(lambdas,lambdaCur)
    sink(tempfile())
    Results[[c]] <- suppressMessages(suppressWarnings(lmmlasso::lmmlasso(x,y,z,grp,..., lambda = lambdaCur, nonpen = 1)))
    sink()
    BICs <- c(BICs,Results[[c]]$bic)
    nVar <- c(nVar,sum(coef(Results[[c]])[-1]!=0))
  
    if (nVar[c]==0){
      break
    } else {
      c <- c + 1
      lambdaCur <- lambdaCur * lambdaIncrease
    }
  }
  
  return(Results[[which.min(BICs)]])
}