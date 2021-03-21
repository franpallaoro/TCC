Moment_based_omega_LT_FC = function(G, params, rho_vec_sample_bl_vech){ # function [dist_sq,L_tilde_prime_L_tilde] = Moment_based_omega_LT_FC(G, params, rho_vec_sample_bl_vech);
  # sa?das  dist_sq, L_tilde_prime_L_tilde
  # This function contains the first step of the two-step
  # approach when estimating the MF-LT model with a scaled loading
  # matrix defined in eq(6). This is a moment based estimator of f_bar.
  
  # The function computes L_m of eq (12), using eq(11).
  
  # Input
  # rho_vec_sample_bl_vech:  a (vech) of the G times G matrix R_M defined below eq(11),
  #                          containing the average within/between industry correlation.
  # G: number of groups (industries)
  # params: parameters
  
  # Output
  # L_m: the squared distance  dist_sq
  # L_tilde_prime_L_tilde: matrix of off-diagonal correlations,
  # see eq(12)################################
  ##########################################
  
  g_vec_cum            = cumsum(seq(G, 1, by = -1)) #cumsum(G:-1:1);
  omega_vec            = params #qual o tamanho desse vetor?
  
  f_t_vec       = t(omega_vec) #omega_vec'; 
  f_t_vec_2     = f_t_vec^2
  
  # build the model implied correlation matrix Rm of eq (12)
  teller_i = 1
  denom_f_vec_g = rep(1, G) #ones(G,1)
  f_prime_mat_t = matrix(0, ncol = G, nrow = G) #zeros(G,G)
  for(i in 1:G){
    g_i_c              = g_vec_cum[i]
    denom_f_vec_g      = denom_f_vec_g + c(rep(0, i-1), t(f_t_vec_2[teller_i:g_i_c]))#denom_f_vec_g + [zeros(i-1,1); f_t_vec_2(teller_i:g_i_c)'];
    f_prime_mat_t[,i]  = c(rep(0, i-1), t(f_t_vec[teller_i:g_i_c]))#f_prime_mat_t(:,i) = #[zeros(i-1,1); f_t_vec(teller_i:g_i_c)'];
    teller_i = g_i_c + 1
  }
  
  
  denom_f_mat_g            = denom_f_vec_g%*%t(rep(1, G))#denom_f_vec_g*ones(1,G);
  L_tilde_prime_mat_t      = f_prime_mat_t/sqrt(denom_f_mat_g)#f_prime_mat_t./sqrt(denom_f_mat_g);
  L_tilde_prime_L_tilde    = L_tilde_prime_mat_t %*% t(L_tilde_prime_mat_t)#L_tilde_prime_mat_t * L_tilde_prime_mat_t';
  vech_FC                  = Sym2Vech(G,L_tilde_prime_L_tilde)#Admin.Sym2Vech(G,L_tilde_prime_L_tilde); #fun?ao do admin
  
  # compute the Frobenius norm of the difference between the sample
  # moments and the model-implied moments
  dist_vec = vech_FC - rho_vec_sample_bl_vech #vech_FC - rho_vec_sample_bl_vech;
  dist_sq  = t(dist_vec)%*%dist_vec #dist_vec'*dist_vec;
  
  return(list("dist_sq" = dist_sq, "L_tilde_prime_L_tilde" = L_tilde_prime_L_tilde))         
}
