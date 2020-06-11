estGAS <- function(specf, data, i, j, nWin, ...){
  fit <- MultiGASFit(GASSpec = specf, data = data[i:j,])
  meanGAS <- fit@Estimates$Moments$mean[nWin,]
  matrixCovGAS <- fit@Estimates$Moments$cov[,,nWin]  
  return(list(mean = meanGAS, matrixCov = matrixCovGAS))
}
