estGAS <- function(specf, data, i, j, nWin, ...){
  fit <- MultiGASFit(GASSpec = specf, data = data[i:j,])
  meanGAS <- fit@Estimates$Moments$mean[nWin+1,]
  matrixCovGAS <- fit@Estimates$Moments$cov[,,nWin+1]  
  return(list(mean = meanGAS, matrixCov = matrixCovGAS))
}
