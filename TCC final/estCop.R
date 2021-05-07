estCop = function(u_mat, asset_group_vec, n_vec, Tns, ind_t_dist){
  
  T          = dim(u_mat)[1]
  N          = dim(u_mat)[2]
  G = max(asset_group_vec)
  n_vec = n_vec*matrix(1,G,1)
  
  if(!require(nloptr)){install.packages("nloptr")}
  library(nloptr)
  
  nr_l  = G*(G+1)/2
  ind_pos =  c(1, cumsum(seq(G, 2, by = -1)) + 1)
  
  Lb_Mf_LT_1 = -3*matrix(1,nr_l,1)
  Lb_Mf_LT_1[ind_pos] = .Machine$double.eps
  Rb_Mf_LT_1 =  3*matrix(1,nr_l,1)
  para_start_Mf_LT_1 = 0.2 * matrix(1,nr_l,1)
  
  Lb_Mf_LT_2 = c(.Machine$double.eps*matrix(1,1,2), 2.1)
  Rb_Mf_LT_2 = c(5*matrix(1,1,1), 0.99999, 5000)
  para_start_Mf_LT_2 = c(0.01, 0.96, 20)
  
  x_mat          = qnorm(u_mat)#norminv(u_mat);
  Rt_Block       = Compute_BLOCK_correlation_matrix(x_mat,n_vec)[[1]]
  Rt_sample_vech = Sym2Vech(G,t(Rt_Block))
  

  step1 = slsqp(x0 = para_start_Mf_LT_1, fn = Moment_based_omega_LT_FC, 
                lower = Lb_Mf_LT_1, upper = Rb_Mf_LT_1,
                G = G, rho_vec_sample_bl_vech = Rt_sample_vech)
  f_bar = step1$par

  step2 = slsqp(x0 = para_start_Mf_LT_2, fn = LogLik_Copula_LT_factor_given_omega,
                lower = Lb_Mf_LT_2, upper = Rb_Mf_LT_2, N = N, T = T, 
                f_hat_vec = f_bar, u_mat = u_mat, 
                n_vec = n_vec, ind_t_dist = 1, 
                ind_Rt = 0, asset_group_vec = asset_group_vec)
  
  theta_opt = step2$par
  
  sim       = Sim_u_mat(N, G, Tns, theta_opt, f_bar, n_vec, 1, 0)
  u_mat_sim = sim[[1]]
  x_mat_sim = sim[[2]]
  
  return(list(cov(u_mat_sim), cor(u_mat_sim), cov(x_mat_sim), cor(x_mat_sim), theta_opt))
  
}

