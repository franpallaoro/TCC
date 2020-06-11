data <- read.csv('6_Portfolios_2x3.CSV', header = T)
data <- data[,-1]
window <- 120
N <- ncol(data)
v <- 5
returnPortSample <- vector()

for (i in 1: window) {
  i = i
  j = i + window
  est <- estSample(data = data, i = i, j = j)
  mean <- est[[1]]
  cov <- est[[2]]
  returnPortSample <- c(returnPortSample, 
                        returnP(wOptim = optimW(N = N, 
                                                mean = mean, 
                                                cov = cov, v = v),
                                data = data, j = j))
}

system("xdg-open 'https://www.youtube.com/watch?v=TOamHghGCfg'")