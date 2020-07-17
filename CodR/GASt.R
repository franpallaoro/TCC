library(GAS)

#fixed objects ------------------

window <- 200
v <- 5
specf <- MultiGASSpec(Dist = "mvt") 

# 6 portfolio ----------------

data <- read.csv('data6.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "1998-01-01" &  data$date <="2019-12-01", -1]
N <- ncol(data)
returnPortGASt6CRRA <- vector()
returnPortGASt6MeanVar <- vector()
returnPortGASt6MinVar <- vector()

for (i in 38:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  
  returnPortGASt6CRRA <- c(returnPortGASt6CRRA, 
                                returnP(wOptim = optimWCRRA(N = N, 
                                                            mean = mean, 
                                                            cov = cov, v = v),
                                        data = data, j = j))
  
  returnPortGASt6MeanVar <- c(returnPortGASt6MeanVar, 
                                   returnP(wOptim = optimWMeanVar(N = N, 
                                                                  mean = mean, 
                                                                  cov = cov, 
                                                                  v = v),
                                           data = data, j = j))
  
  returnPortGASt6MinVar <- c(returnPortGASt6MinVar, 
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
returnPortGASt25CRRA <- vector()
returnPortGASt25MeanVar <- vector()
returnPortGASt25MinVar <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  
  returnPortGASt25CRRA <- c(returnPortGASt25CRRA, 
                                 returnP(wOptim = optimWCRRA(N = N, 
                                                             mean = mean, 
                                                             cov = cov, v = v),
                                         data = data, j = j))
  
  returnPortGASt25MeanVar <- c(returnPortGASt25MeanVar, 
                                    returnP(wOptim = optimWMeanVar(N = N, 
                                                                   mean = mean, 
                                                                   cov = cov, 
                                                                   v = v),
                                            data = data, j = j))
  
  returnPortGASt25MinVar <- c(returnPortGASt25MinVar, 
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
returnPortGASt100CRRA <- vector()
returnPortGASt100MeanVar <- vector()
returnPortGASt100MinVar <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  est <- estGAS(specf = specf, data = data, i = i, 
                j = j, nWin = window)
  mean <- est[[1]]
  cov <- est[[2]]
  returnPortGASt100CRRA <- c(returnPortGASt100CRRA, 
                                  returnP(wOptim = optimWCRRA(N = N, 
                                                              mean = mean, 
                                                              cov = cov, v = v),
                                          data = data, j = j))
  
  returnPortGASt100MeanVar <- c(returnPortGASt100MeanVar, 
                                     returnP(wOptim = optimWMeanVar(N = N, 
                                                                    mean = mean, 
                                                                    cov = cov, 
                                                                    v = v),
                                             data = data, j = j))
  
  returnPortGASt100MinVar <- c(returnPortGASt100MinVar, 
                                    returnP(wOptim = optimWMinVar(N = N, 
                                                                  mean = mean, 
                                                                  cov = cov, 
                                                                  v = v),
                                            data = data, j = j))
}



system("xdg-open 'https://www.youtube.com/watch?v=TOamHghGCfg'")