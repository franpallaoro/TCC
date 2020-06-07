optimW <- function(N, mean, cov, v, ...){
  library(nloptr)
  w0 <- rep(0, N)
  wOptim <- auglag(x0 = w0, 
                   fn = CRRA, heq = constraint, 
                   localsolver = "LBFGS")$par
  return(wOptim)
} 

