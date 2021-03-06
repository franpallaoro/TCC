estCop = function(u_mat, asset_group_vec, Tns, ind_t_dist){
  
  T          = dim(u_mat)[1]
  N          = dim(u_mat)[2]
  G = max(asset_group_vec)
  n_vec = unique(asset_group_vec)*matrix(1,G,1)
  
  if(!require(pracma)){install.packages("pracma")}
  if(!require(NlcOptim)){install.packages("NlcOptim")}
  library(pracma)
  library(NlcOptim)
  
  nr_l  = G*(G+1)/2
  ind_pos =  c(1, cumsum(seq(G, 2, by = -1)) + 1)
  
  Lb_Mf_LT_1 = -3*matrix(1,nr_l,1)
  Lb_Mf_LT_1[ind_pos] = .Machine$double.eps
  Rb_Mf_LT_1 =  3*matrix(1,nr_l,1)
  para_start_Mf_LT_1 = 0.2 * matrix(1,nr_l,1)
  
  Lb_Mf_LT_2 = c(.Machine$double.eps*matrix(1,1,2), 2.1)
  Rb_Mf_LT_2 = c(5*matrix(1,1,1), 0.99999, 5000)
  para_start_Mf_LT_2 = c(0.01, 0.96, 20)
  
  x_mat    = qnorm(u_mat)#norminv(u_mat);
  Rt_Block = Compute_BLOCK_correlation_matrix(x_mat,n_vec)[[1]]
  Rt_sample_vech = Sym2Vech(G,t(Rt_Block))
  

  step1 = fmincon(x0 = para_start_Mf_LT_1, fn = Moment_based_omega_LT_FC, 
                  lb = Lb_Mf_LT_1, ub = Rb_Mf_LT_1,
                  G = G, rho_vec_sample_bl_vech = Rt_sample_vech)
  f_bar = step1$par

  step2 = fmincon(x0 = para_start_Mf_LT_2, fn = LogLik_Copula_LT_factor_given_omega,
                  lb = Lb_Mf_LT_2, ub = Rb_Mf_LT_2, N = N, T = T, 
                  f_hat_vec = f_bar, u_mat = u_mat, 
                  n_vec = n_vec, ind_t_dist = ind_t_dist, ind_Rt = 0, ind_sim = 0)
  
  theta_opt = step2$par
  
  u_mat_sim = Sim_u_mat(N, G, Tns, theta_opt, f_bar, n_vec, 1, 0)
  
  return(list(cov(u_mat_sim), theta_opt))
}

