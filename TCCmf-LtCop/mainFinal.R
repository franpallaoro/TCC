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

u_mat = matrix(0, T, N)#matriz PITs
var_vec_gas = rep(0, N) #variancia modelos marginais
weightsCop = matrix(0, l, N) #matriz guarda os pesos
asset_group_vec #vetor dos grupos
Tns = 1000 #simulações
u_mat_list = list()
cov_mat_cop = list()
params_cop  = list()

for (i in 1:l) {
  i = i
  j = i + window-1
  
  #estimação dos modelos marginais
  for (k in 1:N) {
    aux            = estGAS(data, i, j, k, Tns)
    u_mat[,k]      = aux[[1]]
    var_vec_gas[k] = aux[[2]]
  }
  
  #estimação da cópula
  u_mat_list[[i]]        = u_mat
  auxCop                 = estcop(u_mat, asset_group_vec, Tns, 1)
  cov_mat_cop[[i]]       = auxCop[[1]]
  diag(cov_mat_cop[[i]]) = var_vec_gas
  params_cop[[i]]        = auxCop[[2]]
  
  #amostral
  cov_mat_sample[[i]] = cov(data[i:j,])
  
  #minima variancia
  weightsCop[i,]    = minVar(cov_mat_cop[[i]])
  weightsSample[i,] = minVar(cov_mat_sample[[i]])
  
  #retornos 
  returnCop[i]    = returnP(wOptim = weightsCop[i,], data = data, j = j)
  returnSample[i] = returnP(wOptim = weightsSample[i,], data = data, j = j)
}