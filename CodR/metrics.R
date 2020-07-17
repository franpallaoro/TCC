
#função -------------
metrics <- function(x){
  meanR <- mean(x, na.rm = T)*sqrt(12)
  sdR <- sd(x, na.rm = T)*sqrt(12) 
  SR <- meanR/sdR
  return(list(meanR = meanR, sdR = sdR, SR = SR))
}




#portfolio 6------------------

dataPort6 <- data.frame(returnPortEW6, returnPortGASnormal6CRRA, 
                        returnPortGASnormal6MeanVar, 
                        returnPortGASnormal6MinVar,
                        returnPortGASt6CRRA, returnPortGASt6MeanVar,
                        returnPortGASt6MinVar, returnPortSample6CRRA,
                        returnPortSample6MeanVar, returnPortSample6MinVar)

write.csv(dataPort6, file = "dataPort6")

metricsPort6 <- as.data.frame(unlist(apply(dataPort6, 2, metrics)))
write.csv(metricsPort6, file = "metricsPort6")
