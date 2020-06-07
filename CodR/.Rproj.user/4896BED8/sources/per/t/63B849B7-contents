CRRA <-function(w, ....){ 
  ut <- - (( 1/(1-v) ) + 
             (t(w)%*%mean) - (v * ( (t(w)%*%cov%*%w) + 
                                      (t(w)%*%mean)^2 )/2)) 
  return(ut)
}

constraint <- function(w,...){
  sum(w) - 1
}
