estGAS <- function(data, i, j, k, ...){
  #i-j janela
  #k é o ativo 1...N
  GASSpec = UniGASSpec(Dist = "std", ScalingType = "Identity",
                       GASPar = list(location = FALSE, scale = TRUE,
                                     shape = FALSE))
  fit = UniGASFit(GASSpec, data[i:j,k])
  sim = UniGASSim(fit = fit, T.sim = 1000)
  var_sim = var(sim@Data$vY)
  u_i = fit@Estimates$vU
  return(list(u_i, var_sim))
}





