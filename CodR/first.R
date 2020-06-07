library(GAS)

dados <- read.csv('6_Portfolios_2x3.CSV', header = T)

inSample <- dados[1:190,]
outSample <- dados[190:nrow(dados),]

#Media e matriz de cov amostral--------
estSample <- function(data,i,j,...){
  meanSample <- apply(data[i:j,], 2, mean)
  matrixCovSample <- cov(data[i:j,])
  return(list(mean = meanSample, matrixCov = matrixCovSample))
}


#GAS ---------------------------
specf <- MultiGASSpec()
estGAS <- function(specf, data, i, j, nW, ...){
  fit <- MultiGASFit(GASSpec = specf, data = data[i:j,])
  meanGAS <- as.matrix(fit@Estimates$Moments$mean[nW,])
  matrixCovGAS <- fit@Estimates$Moments$cov[,,nW]  
  return(list(mean = meanGAS, matrixCov = matrixCov))
}




#otim--------

w <- as.matrix(rep(0, ncol(inSample)-1))
muP <- t(w)%*%M1
sig2P <- t(w)%*%M2%*%w
m1P <- muP
m2P <- (t(w)%*%M2%*%w) + (t(w)%*%M1)^2

otim <- function(w, M1, M2, v,...){
 #ut <- 
  w <- Variable(1,6)
  constraint <- list(sum(w) == 1)
  objective <- Minimize()
  problem <- Problem(objective)
  result <- solve(problem)
  return(result)
}

w <- rep(0,6)

otim(w = w, M1 = meanGAS, M2 = covGAS, v = 5)

wOtim <- optim(par = w, 
               fn = CRRA, 
               M1 = meanGAS, 
               M2 = covGAS,
               v = 5)$par



library(quadprog)

uns <-as.matrix(rep(1, 6)) 
otim <- solve.QP(covGAS, meanGAS, uns, 1 )$solution
sum(otim)

sum(sol*outSample[1, 2:7])
sum(outSample[1, 2:7])/6


# funcionou para otimização ---------------------

library(nloptr)

ut <-function(w){ - (( 1/(1-5) ) + 
           (t(w)%*%meanGAS) - 
           (5 * ( (t(w)%*%covGAS%*%w) + (t(w)%*%meanGAS)^2 )/2)) }

w0 <- rep(0,6)

constraint <- function(w){
  sum(w) - 1
}

sol <- auglag(x0 = w0, fn = ut, heq = constraint, localsolver = "LBFGS")$par
sum(sol)
