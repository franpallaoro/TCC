library(GAS)

data <- read.csv('6_Portfolios_2x3.CSV', header = T)
data <- data[,-1]
window <- 120
N <- ncol(data)
v <- 5
returnPortGASt <- vector()
specf <- MultiGASSpec(Dist = "mvt") 

for (i in 1: window) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  returnPortGASt <- c(returnPortGASt, 
                      returnP(wOptim = optimW(N = N, 
                                              mean = mean, 
                                              cov = cov, v = v),
                                   data = data, j = j))
}

system("xdg-open 'https://www.youtube.com/watch?v=TOamHghGCfg'")