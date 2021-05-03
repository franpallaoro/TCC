
#função -------------
metrics <- function(x){
  meanR <- mean(x, na.rm = T)*sqrt(12)
  sdR <- sd(x, na.rm = T)*sqrt(12) 
  SR <- meanR/sdR
  return(list(meanR = meanR, sdR = sdR, SR = SR))
}

turnover <- function(w, l, N){
  difW <- matrix(NA, ncol = N, nrow = l-1)
  for (i in 1:(l - 1)) {
    difW[i,] <- abs(w[i+1,] - w[i,])
  }
  turnover <- sum(difW)/(l-1)
  return(turnover)
}


