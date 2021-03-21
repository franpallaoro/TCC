
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


#portfolio 6------------------

dataPort6 <- data.frame(returnPortEW6, returnPortGASnormal6CRRA, 
                        returnPortGASnormal6MeanVar, 
                        returnPortGASnormal6MinVar,
                        #returnPortGASt6CRRA, returnPortGASt6MeanVar,
                        #returnPortGASt6MinVar,
                        returnPortSample6CRRA,
                        returnPortSample6MeanVar, returnPortSample6MinVar)

write.csv(dataPort6, file = "dataPort6")

dataWPort6 <- list(wOptimGASnormal6CRRA, wOptimGASnormal6MeanVar, 
                   wOptimGASnormal6MinVar, 
                   #wOptimGASt6CRRA, wOptimGASt6MeanVar, wOptimGASt6MinVar,
                   wOptimSample6CRRA,
                   wOptimSample6MeanVar, wOptimSample6MinVar)

write.csv(dataWPort6, file = "dataWPort6")

metricsPort6 <- as.data.frame(unlist(apply(dataPort6, 2, metrics)))
write.csv(metricsPort6, file = "metricsPort6")

turnoverPort6 <- vector()
for (i in 1:6) {
  turnoverPort6[i] <- turnover(dataWPort6[[i]], l, N) 
}
