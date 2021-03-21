returnP <- function(wOptim, data, j, ....){
  rp <- sum(wOptim*data[j+1,])
  return(rp)
}

