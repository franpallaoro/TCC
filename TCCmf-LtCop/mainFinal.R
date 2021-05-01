source("Compute_BLOCK_correlation_matrix.R")
source("Compute_R_t_inv_and_det_R_t_L_matrix.R")
source("Compute_S_l_matrix.R")
source("LogLik_Copula_LT_factor_given_omega.R")
source("Moment_based_omega_LT_FC.R")
source("Sim_u_mat.R")
source("Sym2Vech.R")
source("LogLik_Copula_1f_eq.R")
source("Compute_R_t.R")
source("Compute_R_t_inv_and_det_R_t.R")
source("Sim_x_mat_1f_eq.R")
source("estGas.R")
source("estCop.R")

window <- 200
l <- nrow(data) - window
T = nrow(data)
N = ncol(data)
u_mat = matrix(0, T, N)
var_vec = rep(0, N)
asset_group_vec #vetor dos grupos
Tns = 1000 #simulações

for (i in 1:l) {
  i = i
  j = i + window-1
  
  #estimação dos modelos marginais
  for (k in 1:N) {
    aux = estGAS(data, i, j, k)
    u_mat[,k]  = aux[[1]]
    var_vec[k] = aux[[2]]
  }
  
  #estimação da cópula
  cov_mat = estcop(u_mat, asset_group_vec, Tns)
  
  diag(cov_mat) = var_vec
  
  wOptimGASt6MinVar[i,] <- optimWMinVar(N = N, mean = mean, cov = cov, v = v)
  returnPortGASt6MinVar <- c(returnPortGASt6MinVar, 
                             returnP(wOptim = wOptimGASt6MinVar[i,],
                                     data = data, j = j))
}