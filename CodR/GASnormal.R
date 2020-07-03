library(GAS)

#fixed objects ------------------

window <- 200
v <- 5
specf <- MultiGASSpec() 

# 6 portfolio ----------------

data <- read.csv('data6.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "1998-01-01" &  data$date <="2019-12-01", -1]
N <- ncol(data)
returnPortGASnormal6CRRA <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  
  returnPortGASnormal6CRRA <- c(returnPortGASnormal6CRRA, 
                                returnP(wOptim = optimWCRRA(N = N, 
                                                          mean = mean, 
                                                          cov = cov, v = v),
                                        data = data, j = j))
  
  returnPortGASnormal6MeanVar <- c(returnPortGASnormal6MeanVar, 
                                returnP(wOptim = optimWMeanVar(N = N, 
                                                              mean = mean, 
                                                              cov = cov, 
                                                              v = v),
                                        data = data, j = j))
  
  returnPortGASnormal6MinVar <- c(returnPortGASnormal6MinVar, 
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
returnPortGASnormal25CRRA <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  
  returnPortGASnormal25CRRA <- c(returnPortGASnormal25CRRA, 
                                 returnP(wOptim = optimWCRRA(N = N, 
                                                             mean = mean, 
                                                             cov = cov, v = v),
                                         data = data, j = j))
  
  returnPortGASnormal25MeanVar <- c(returnPortGASnormal25MeanVar, 
                                    returnP(wOptim = optimWMeanVar(N = N, 
                                                                   mean = mean, 
                                                                   cov = cov, 
                                                                   v = v),
                                            data = data, j = j))
  
  returnPortGASnormal25MinVar <- c(returnPortGASnormal25MinVar, 
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
returnPortGASnormal100CRRA <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  returnPortGASnormal100CRRA <- c(returnPortGASnormal100CRRA, 
                                  returnP(wOptim = optimWCRRA(N = N, 
                                                              mean = mean, 
                                                              cov = cov, v = v),
                                          data = data, j = j))
  
  returnPortGASnormal100MeanVar <- c(returnPortGASnormal100MeanVar, 
                                     returnP(wOptim = optimWMeanVar(N = N, 
                                                                    mean = mean, 
                                                                    cov = cov, 
                                                                    v = v),
                                              data = data, j = j))
  
  returnPortGASnormal100MinVar <- c(returnPortGASnormal100MinVar, 
                                    returnP(wOptim = optimWMinVar(N = N, 
                                                                  mean = mean, 
                                                                  cov = cov, 
                                                                  v = v),
                                            data = data, j = j))
}



system("xdg-open 'https://www.youtube.com/watch?v=TOamHghGCfg'")