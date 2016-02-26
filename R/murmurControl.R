# Setup consistency:
# FUNCTIONNAME_ARGUMENT = ....
# some defaults set...

murmurControl <- function(
# jags:
  jags_n.iter = 2000,
  jags_n.burnin = 500,
  jags_n.chain = 3,
  
# lmer:
  lmer_control = list(optimizer = "bobyqa"),
  ...
){

  # Contorol list:
  Controls <- list(
    jags = list(
      n.iter = jags_n.iter,
      n.burnin = jags_n.burnin,
      n.chain = jags_n.chain
    ),
    lmer = list(
      control = do.call('lmerControl',lmer_control)
    )
  )
  class(Controls) <- "murmurControl"
  
  # add dots:
  dots <- list(...)
  if (length(dots)>1){
    names <- stringr::str_split_fixed(names(dots), pattern = '_', n = 2)
    for (i in 1:nrow(names)){
      Controls[[names[i,1]]][[names[i,2]]] <- dots[[i]]
    }
    
  }
  

  return(Controls)
}