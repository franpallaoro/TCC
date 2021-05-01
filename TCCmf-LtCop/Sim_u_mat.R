Sim_u_mat = function(N, G, T, params, mu_moment, n_vec, ind_t_dist, ind_Rt){
#entrar com A B omega e nu
#entrar com x inicial sefor 0 não precisa
#entrar com indicadora se é t ou gaussian
#entrar com n_vec que é o número de ativos por grupo
#entrar com T, G e N
#sair com PITs


g_vec_cum            = cumsum(seq(G, 1, by = -1))
k                    = G
p                    = G*(G+1)/2

A_vec     = params[1]
B_vec     = params[2]
if(ind_t_dist == 1){
  end       = length(params)
  nu        = params[end]
}
#omega_vec = (1-B_vec)%*%t(f_hat_vec)
omega_vec  =  t(mu_moment)*(1-B)


f_mat               = matrix(0, nrow = T, ncol = p)
s_mat               = matrix(0, nrow = T, ncol = p)
x_mat               = matrix(0, nrow = T, ncol = N)
u_mat               = matrix(0, nrow = T, ncol = N)

if(ind_Rt==1){
  R_mat   = array(0, dim = c(G, G, T)) 
}else{
  R_mat   = vector() 
}



f_mat[1,] = mu_moment


teller_0 = 1
asset_group_vec_new = matrix(0, nrow = N, ncol = 1)
for(m in 1:G){
  asset_group_vec_new[teller_0:(teller_0+n_vec[m]-1)] = m
  teller_0                                          = teller_0 + n_vec[m]
}


S_l_mat    = Compute_S_l_matrix(asset_group_vec_new)[[1]] 

Ind_lambda_to_group_vec = vector()
for(i in 1:G){
  Ind_lambda_to_group_vec = c(Ind_lambda_to_group_vec, t(i:G))
}
Ind_lambda_to_group_mat      = cbind(1:p, Ind_lambda_to_group_vec)

# we are now ready to build "dau_L_tilde_dau_lambda_tilde
dau_L_tilde_dau_lambda_tilde = matrix(0, nrow = N*k, ncol = p)

ind_match_lambda_group_cell = matrix(list(), G, 1) 


for(i in 1:G){
  # step 2
  if(i == 1){
    aux = seq(from = 1, to = sum(n_vec[1:i])*G, by = G)
    dau_L_tilde_dau_lambda_tilde[aux,i] = 1 
  }else{
    aux = seq(from = sum(n_vec[1:(i-1)])*G+1, to = sum(n_vec[1:i])*G, by = G)
    dau_L_tilde_dau_lambda_tilde[aux,i] = 1 
  }
  
  
  # step 3
  ind_lambda_g_group_i           = Ind_lambda_to_group_mat[Ind_lambda_to_group_mat[ ,2] ==i, 1]
  ind_match_lambda_group_cell[[i]] = ind_lambda_g_group_i
  
  # step 4
  for(j in 1:(i-1)){
    lambda_j                   = ind_lambda_g_group_i[j+1]
    aux1 = seq(sum(n_vec[1:i-1])*G+1+j, sum(n_vec[1:i])*G+j, by = G)
    dau_L_tilde_dau_lambda_tilde[aux1,lambda_j] =1
  }
}

for(j in 1:T){
  # Step 0: compute the time varying coefficient and the
  #         data vector at this time
  
  f_t_vec   = f_mat[j,]
  f_t_vec_2 = f_t_vec^2
  x_j       = t(x_mat[j,])
  
  # Build the matrix L_til_prime_mat_t. See (F.19) of the web
  # appendix.                 
  # We also build a second matrix L_tilde_prime_mat_tg (F.19)
  # with N_g = 1 for all values of g)                 
  # This G x G matrix is based on only unique elements of f_j,t.
  # We (could) use this later on to build the G x G block correlation
  # matrix with within and between group correlations.
  
  teller_i        = 1
  denom_f_un      = matrix(1, nrow = G, ncol = 1)
  denom_f_vec     = matrix(1, nrow = N, ncol = 1)
  f_prime_mat_t   = matrix(0, nrow = N, ncol = G)
  f_prime_mat_tg  = matrix(0, nrow = G, ncol = G)
  
  for(i in 1:k){ 
    if(i != G){
      g_i_c           = g_vec_cum[i]
      denom_f_vec     = denom_f_vec + S_l_mat[,i:G]%*%f_t_vec_2[teller_i:g_i_c]  # Based on (F.19) and used repeatedly. 
      denom_f_un      = denom_f_un + c(matrix(0, nrow = i-1, ncol = 1), f_t_vec_2[teller_i:g_i_c])  # This is (F.21) and used repeatedly.
      f_prime_mat_t[,i]  = S_l_mat[,i:G]%*%f_t_vec[teller_i:g_i_c]
      f_prime_mat_tg[,i] = c(matrix(0, nrow = i-1, ncol = 1), f_t_vec[teller_i:g_i_c])
      teller_i            = g_i_c+1
    }else{
      g_i_c           = g_vec_cum[i]
      denom_f_vec     = denom_f_vec + S_l_mat[,i:G]*f_t_vec_2[teller_i:g_i_c]  # Based on (F.19) and used repeatedly. 
      denom_f_un      = denom_f_un + c(matrix(0, nrow = i-1, ncol = 1), f_t_vec_2[teller_i:g_i_c])  # This is (F.21) and used repeatedly.
      f_prime_mat_t[,i]  = S_l_mat[,i:G]*f_t_vec[teller_i:g_i_c]
      f_prime_mat_tg[,i] = c(matrix(0, nrow = i-1, ncol = 1), f_t_vec[teller_i:g_i_c])
      teller_i            = g_i_c+1
    }
    
  }
  
  # Step 1b: compute the correlation matrix (if ind_Rt=1) and its
  #          inverse and determinant
  denom_f_mat             = denom_f_vec%*%matrix(1, nrow = 1,ncol = k)
  lambda_til_prime_mat_t  = f_prime_mat_t/sqrt(denom_f_mat)
  sigma_2_vec_t           = (1/denom_f_vec)
  aux_Compute_R_t_inv_and_det_R_t_L_matrix = Compute_R_t_inv_and_det_R_t_L_matrix(sigma_2_vec_t,lambda_til_prime_mat_t)
  R_inv_t                 = as.matrix(aux_Compute_R_t_inv_and_det_R_t_L_matrix[[1]])
  det_R_t                 = aux_Compute_R_t_inv_and_det_R_t_L_matrix[[2]]
  Flag                    = aux_Compute_R_t_inv_and_det_R_t_L_matrix[[3]]
  if(ind_Rt==1){
    denom_f_mat_tg          = denom_f_un%*%matrix(1,1,G)
    L_tilde_prime_mat_tg    = f_prime_mat_tg/sqrt(denom_f_mat_tg)
    Rt_Block                = L_tilde_prime_mat_tg%*%t(L_tilde_prime_mat_tg)
    R_mat[,,j] = Rt_Block
  }
  
  if(Flag==1 || is.nan(det_R_t)){
    loglike_vec = 1e14
    break 
  }else{
    
    #Step 2: Compute the likelihood for time t=j and store it,
    #         compute the score function and update the time-varying
    #         parameters
    
    # Compute the score is factor-model-specific.
    
    # Step 2a: compute R^{-1} x(t) and x(t)' R^{-1} x(t)                    
    coef_R_1    = -0.5
    R_inv_x     = R_inv_t %*% t(x_j)
    x_R_inv_x   = x_j%*%R_inv_x
    # Step 2b: compute the derivative of the copula density with respect to
    #          R(t) as in (A8)-(A9)                    
    if(ind_t_dist==1){
      coef_R_2 = 0.5*(nu + N)/(nu-2+ x_R_inv_x ) 
    }else{
      coef_R_2 = 0.5 
    }
    R_2 = R_inv_x%*% t(R_inv_x)  # outer product: R^{-1} x(t) x(t)' R^{-1}
    
    # dau_log_c_dau_vec_R
    A_mat = as.numeric(coef_R_2) * R_2 + coef_R_1 * R_inv_t 
    
    # Step 2c: compute the score, using  (F.2),(F.3) (F.20)-(F.25)
    # of the online appendix.
    
    # auxilary term used repeatedly 
    denom_f_un_3_2           = denom_f_un^(3/2)
    
    # (F.22) for j = 1,2,..., G(G+1/2)
    dau_lambda_tilde_dau_f_vec  = 1/sqrt(denom_f_un[Ind_lambda_to_group_mat[,2]]) - t(f_t_vec_2)/denom_f_un_3_2[Ind_lambda_to_group_mat[,2]]
    dau_lambda_tilde_dau_f_mat  = diag(length(dau_lambda_tilde_dau_f_vec))
    diag(dau_lambda_tilde_dau_f_mat) = dau_lambda_tilde_dau_f_vec             
    
    dau_sigma2_dau_f_mat         = matrix(0, G, p)
    
    # fill dau_lambda_tilde_dau_f and dau_sigma2_dau_f
    
    # The loop below computes (F.23) for all values of j, looping over all G columns. 
    # That is: check which values of f_(t,j) are in column g. Then compute the cross deratives.                     
    
    # The loop also computes (F.25) for each g (g = 1,...,G) 
    # That is: check which values of f_(t,j) are in column g. 
    # Then compute (F.25) w.r.t. to all these f_(t,j) values.                      
    for(g in 1:G){
      # find out which f_(t,j) are in column g
      ind_cross_lambda_g =  ind_match_lambda_group_cell[[g]]
      
      if(g>1){
        for(m in 1:g){
          ind_2   = ind_cross_lambda_g
          ind_m1  = ind_cross_lambda_g[m]
          ind_2 = ind_2[-m]
          
          #(F.23)
          dau_lambda_tilde_i_dau_m_vec             = -(f_t_vec[ind_m1]*f_t_vec[ind_2])/denom_f_un_3_2[g]
          dau_lambda_tilde_dau_f_mat[ind_m1,ind_2] = dau_lambda_tilde_i_dau_m_vec
        }
      }
      
      #(F.25)
      dau_sigma2_dau_f_mat[g,ind_cross_lambda_g] = -2*f_t_vec[ind_cross_lambda_g]/(denom_f_un[g]^2)
    }
    
    
    #(F.20)
    dau_L_tilde_dau_f_mat =  dau_L_tilde_dau_lambda_tilde%*%dau_lambda_tilde_dau_f_mat
    
    
    # (F.24)                                         
    dau_D_dau_f_mat = S_l_mat%*%dau_sigma2_dau_f_mat
    
    # combining above with (F.2) and (F.3)
    L_A         = t(lambda_til_prime_mat_t)%*%A_mat
    s_vec_a     = 2* t(c(L_A))%*%dau_L_tilde_dau_f_mat
    s_vec_b     = t(diag(A_mat))%*%dau_D_dau_f_mat
    s_mat[j,]  = s_vec_a + s_vec_b
    
    # Step 2d: update the time-varying parameters using the score for the next
    #          iteration
    if(j<T){
      f_mat[j+1,] = omega_vec +  A_vec*s_mat[j,] + B_vec*f_mat[j,] #código da verossimilhança
      }
    if(!require(invgamma)){install.packages("invgamma")}
    if(!require(mvtnorm)){install.packages("mvtnorm")}
    library(invgamma)
    library(mvtnorm)
    if(ind_t_dist == 1){
      ginv = rinvgamma(1, 1/2*nu, 1/2*nu)#comum pra todos os i 
      zt = rnorm(1)#comum pra todos os i
      for (i in 1:N) {
        lambda_tilde_x_i = as.matrix(lambda_til_prime_vec_t[i])
        sigma_x_i = sqrt(sigma_2_vec_t[i])
        eps_i = rnorm(1)#específico para x_i
        x_mat[j,i] = sqrt(ginv)*(lambda_tilde_x_i*zt + sigma_x_i*eps_i)
        x_j_scaled = x_mat[j,i] * 1/ sqrt((nu-2)/nu);
        u_mat[j,i] = pt(x_j_scaled,nu)
      }
    }else{
      zt = rnorm(1)#comum pra todos os i
      for (i in 1:N) {
        lambda_tilde_x_i = as.matrix(lambda_til_prime_vec_t[i])
        sigma_x_i = sqrt(sigma_2_vec_t[i])
        eps_i = rnorm(1)#específico para x_i
        x_mat[j,i] = (lambda_tilde_x_i*zt + sigma_x_i*eps_i)
        u_mat[j,i] = pnorm(x_mat[j,i])
        }
      }   
    }
  }


return(u_mat)

}
  

