optimWCRRA <- function(N, mean, cov, v, ...){
  library(CVXR)
  
  w <- Variable(N)
  r <- quad_form(w, cov)
  objective <- Maximize( ( 1/(1-v) ) + 
                           (t(w)%*%mean) - (v * ( (r) + (t(w)%*%mean)^2 )/2) )
  
  constraints <- list(w >= 0, sum(w) == 1)
  problem <- Problem(objective, constraints)
  result <- solve(problem)
  
  wOptim <- result$getValue(w)
  
  return(wOptim)
} 


optimWMeanVar <- function(N, mean, cov, v, ...){
  library(CVXR)
  
  w <- Variable(N)
  r <- quad_form(w, cov)
  objective <- Minimize( r - ((1/v) * t(mean)%*%w))
  
  constraints <- list(w >= 0, sum(w) == 1)
  problem <- Problem(objective, constraints)
  result <- solve(problem)
  
  wOptim <- result$getValue(w)
  
  return(wOptim)
} 

optimWMinVar <- function(N, mean, cov, v, ...){
  library(CVXR)
  
  w <- Variable(N)
  r <- quad_form(w, cov)
  objective <- Minimize(r)
  
  constraints <- list(w >= 0, sum(w) == 1)
  problem <- Problem(objective, constraints)
  result <- solve(problem)
  
  wOptim <- result$getValue(w)
  
  return(wOptim)
} 