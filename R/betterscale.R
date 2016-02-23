betterscale <- function(x,center=TRUE,scale=TRUE){
  
  if (sd(x,na.rm=TRUE)==0){
    return(0) 
  } else {
    if (center){
      x <- x - mean(x,na.rm=TRUE)
    }
    if (scale){
      x <- x / sd(x, na.rm=TRUE)
    }
    return(x)
  }
}
