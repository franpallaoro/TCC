estGAS <- function(specf, data, i, j, ...){
  fit <- MultiGASFit(GASSpec = specf, data = data[i:j,])
  n <- nrow(fit@Estimates$Moments$mean)
  meanGAS <- fit@Estimates$Moments$mean[n,]
  matrixCovGAS <- fit@Estimates$Moments$cov[,,n]  
  return(list(mean = meanGAS, matrixCov = matrixCovGAS))
}
