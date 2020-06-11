
data <- read.csv('6_Portfolios_2x3.CSV', header = T)
data <- data[,-1]
window <- 120
N <- ncol(data)
returnPortEW <- vector()


for (i in 1: window) {
  i = i
  j = i + window
  returnPortEW <- sum(data[j+1,])/N
}