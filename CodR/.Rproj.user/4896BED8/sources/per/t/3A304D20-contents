estSample <- function(data,i,j,...){
  meanSample <- apply(data[i:j,], 2, mean)
  matrixCovSample <- cov(data[i:j,])
  return(list(mean = meanSample, matrixCov = matrixCovSample))
}
