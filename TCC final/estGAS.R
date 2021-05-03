estGAS <- function(data, i, j, k, Tns, ...){
  #i-j janela
  #k é o ativo 1...N
  if(!require(GAS)){install.packages("GAS")}
  library(GAS)
  GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                       GASPar = list(location = FALSE, scale = TRUE,
                                     shape = FALSE))
  fit = UniGASFit(GASSpec, data[i:j,k])
  sim = UniGASSim(fit = fit, T.sim = Tns)
  var_sim = var(sim@Data$vY)
  u_i = as.vector(fit@Estimates$vU)
  return(list(u_i, var_sim))
}





