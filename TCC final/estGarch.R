estGarch = function(data,i,j,k,...){
  if(!require(fGarch)){install.packages("fGarch")}
  library(fGarch)
  
  fit = garchFit(~arma(2,0)+garch(1,1), data = data[i:j,k],
                 cond.dist = "sstd", include.mean = F, trace = F)
  
  
  var = as.numeric((predict(fit, 1)[3])^2)
  
  u_i = as.vector(psstd(residuals(fit, standard= T),
              nu = coef(fit)[7], 
              xi = coef(fit)[6]))

  return(list(u_i, var))
}



