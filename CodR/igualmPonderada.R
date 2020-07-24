

window <- 200

#6portfolio --------------------------------

data <- read.csv('data6.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "1998-01-01" &  data$date <="2019-12-01", -1]
N <- ncol(data)
returnPortEW6 <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window-1
  returnPortEW6 <- c(returnPortEW6, sum(data[j+1,])/N)
}

#25 portfolio ---------------------------------

data <- read.csv('data25.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "1998-01-01" &  data$date <="2019-12-01", -1]
N <- ncol(data)
returnPortEW25 <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  returnPortEW25 <- c(returnPortEW25, sum(data[j+1,])/N)
}

#100 portfolio ------------------------------

data <- read.csv('data100.csv', header = T, na.strings = "NA")
data$date <- as.Date(paste0(as.character(data$date), '01'), format='%Y%m%d')
data <- data[data$date >= "1998-01-01" &  data$date <="2019-12-01", -1]
N <- ncol(data)
returnPortEW100 <- vector()

for (i in 1:(nrow(data) - window)) {
  i = i
  j = i + window
  returnPortEW100 <- c(returnPortEW100, sum(data[j+1,])/N)
}