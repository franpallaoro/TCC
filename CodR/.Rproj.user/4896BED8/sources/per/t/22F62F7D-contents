library(GAS)

data <- read.csv('6_Portfolios_2x3.CSV', header = T)
data <- data[,-1]
window <- 120
N <- ncol(data)
v <- 5
returnPortGASnormal <- vector()
specf <- MultiGASSpec() 

for (i in 1: window) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = 1, 
                j = 120, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  returnPortGASnormal <- c(returnPortGASnormal, 
                          returnP(wOptim = optimW(N = N, 
                                                  mean = mean, 
                                                  cov = cov, v = v),
                                  data = data, j = j))
}

system("xdg-open 'https://www.youtube.com/watch?v=TOamHghGCfg'")