library(GAS)

#fixed objects ------------------

window <- 200
v <- 5


# 6 portfolio ----------------

data <- read.csv('data6.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "1998-01-01" &  data$date <="2019-12-01", -1]
N <- ncol(data)
returnPortSample6CRRA <- vector()
returnPortSample6MeanVar <- vector()
returnPortSample6MinVar <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estSample(data = data, i = i, j = j)
  mean <- est[[1]]
  cov <- est[[2]]
  
  returnPortSample6CRRA <- c(returnPortSample6CRRA, 
                            returnP(wOptim = optimWCRRA(N = N, 
                                                        mean = mean, 
                                                        cov = cov, v = v),
                                    data = data, j = j))
  
  returnPortSample6MeanVar <- c(returnPortSample6MeanVar, 
                                returnP(wOptim = optimWMeanVar(N = N, 
                                                              mean = mean, 
                                                              cov = cov, 
                                                              v = v),
                                        data = data, j = j))
  
  returnPortSample6MinVar <- c(returnPortSample6MinVar, 
                              returnP(wOptim = optimWMinVar(N = N, 
                                                            mean = mean, 
                                                            cov = cov, 
                                                            v = v),
                                      data = data, j = j))
}


# 25 portfolio -----------------------------

data <- read.csv('data25.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "2000-01-01" &  data$date <="2019-12-01", -1]
returnPortSample25CRRA <- vector()
returnPortSample25MeanVar <- vector()
returnPortSample25MinVar <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estSample(data = data, i = i, j = j)
  mean <- est[[1]]
  cov <- est[[2]]
  
  returnPortSample25CRRA <- c(returnPortSample25CRRA, 
                              returnP(wOptim = optimWCRRA(N = N, 
                                                          mean = mean, 
                                                          cov = cov, v = v),
                                      data = data, j = j))
  
  returnPortSample25MeanVar <- c(returnPortSample25MeanVar, 
                                 returnP(wOptim = optimWMeanVar(N = N, 
                                                                mean = mean, 
                                                                cov = cov, 
                                                                v = v),
                                        data = data, j = j))
  
  returnPortSample25MinVar <- c(returnPortSample25MinVar, 
                                returnP(wOptim = optimWMinVar(N = N, 
                                                              mean = mean, 
                                                              cov = cov, 
                                                              v = v),
                                        data = data, j = j))
}



#100 portfolio -------------------------
data <- read.csv('data100.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "2000-01-01" &  data$date <="2019-12-01", -1]
returnPortSample100CRRA <- vector()
returnPortSample100MeanVar <- vector()
returnPortSample100MinVar <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estSample(data = data, i = i, j = j)
  mean <- est[[1]]
  cov <- est[[2]]
  returnPortSample100CRRA <- c(returnPortSample100CRRA, 
                               returnP(wOptim = optimWCRRA(N = N, 
                                                           mean = mean, 
                                                           cov = cov, v = v),
                                      data = data, j = j))
  
  returnPortSample100MeanVar <- c(returnPortSample100MeanVar, 
                                  returnP(wOptim = optimWMeanVar(N = N, 
                                                                mean = mean, 
                                                                cov = cov, 
                                                                v = v),
                                          data = data, j = j))
  
  returnPortSample100MinVar <- c(returnPortSample100MinVar, 
                                returnP(wOptim = optimWMinVar(N = N, 
                                                              mean = mean, 
                                                              cov = cov, 
                                                              v = v),
                                        data = data, j = j))
}

system("xdg-open 'https://www.youtube.com/watch?v=TOamHghGCfg'")