Sim_x_mat_1f_eq = function(N, T, params, ind_t_dist){
  # N = N
  # T = T
  # params = LogLik_Copula_1f_eq_optim$par
  # ind_t_dist = 1
  omega_l = params[1]
  A    = params[2]
  B    = params[3]
  
  if(ind_t_dist == 1){
    end       = length(params)
    nu        = params[end]
  }
  
  #Initialize the storage objects for the output
  loglike_vec = matrix(0, T, 1)
  f_vec       = matrix(0, T, 1)
  s_vec       = matrix(0, T, 1)
  R_vec       = matrix(0, T, 1)
  x_mat       = matrix(0, T, N)
  u_mat       = matrix(0, T, N)
  #initialize fixed matrices required during computations
  i_vec = matrix(1, N, 1)
  
  # Compute the starting value of f(t)
  #For the one-factor-equi model, we target the average
  # sample correlation and then use the parameterization
  # of the loading as rho = f / sqrt(1+f*f)..
  #RC_eq       = 1/(N*(N-1))*(t(i_vec)%*%cor(u_mat)%*%i_vec-N)
  #f_vec[1]    = sqrt(RC_eq/(1-RC_eq))
  f_vec[1]    = (1-B)^(-1)*omega_l
  # Loop over all observations to compute the time varying
  # parameter f(t) and implied time-varying correlation
  # parameter.
  for(j in 1:T){
    # Step 0: compute the time varying coefficient and the
    #         data vector at this time
    f_t     = f_vec[j]
    f_2_t   = f_t^2
    x_j     = x_mat[j,]
    
    # Step 1a: compute the factor loading, the idiosyncratic
    #          variance, and some auxiliary quantities
    denom_t                 = 1 + f_2_t	#denominator value, used repeatedly
    lambda_til_prime_vec_t  = f_t/sqrt(denom_t)*i_vec
    sigma_2_vec_t           = i_vec%*%(1/denom_t)
    
    
    # Step 1b: compute the correlation matrix and its
    #          inverse and determinant
    R_inv_t = Compute_R_t_inv_and_det_R_t(sigma_2_vec_t,lambda_til_prime_vec_t)[[1]]
    det_R_t = Compute_R_t_inv_and_det_R_t(sigma_2_vec_t,lambda_til_prime_vec_t)[[2]]
    R_t               = Compute_R_t(sigma_2_vec_t,lambda_til_prime_vec_t)
    # Step 1c: store the time-varying correlation for this model
    R_vec[j,1]        = R_t[1,2]
    
    # CHECKUP: if something wrong with the correlation matrix (not PDS),
    # then return likelihood failure by huge -loglik value.
    if (det_R_t <=0 | is.nan(det_R_t)){
      loglike_vec = 1e14
      break
    }else{
      # Step 2: Compute the likelihood for time t=j and store it,
      #         compute the score function and update the time-varying
      #         parameter
      
      # Step 2a: compute R^{-1} x(t) and x(t)' R^{-1} x(t)
      coef_R_1    = -0.5
      R_inv_x     = R_inv_t%*%x_j
      x_R_inv_x   = x_j%*%t(R_inv_x)#t(x_j)%*%R_inv_x
      
      # Step 2b: compute the derivative of the copula density with respect to
      #          R(t) as in (A8)-(A9)
      if(ind_t_dist==1){
        coef_R_2 = 0.5*(nu + N)/(nu-2+ x_R_inv_x )
      }else{
        coef_R_2 = 0.5
      }
      R_2 = R_inv_x%*%t(R_inv_x)    # outer product: R^{-1} x(t) x(t)' R^{-1}
      # dau log c/ dau vec R_t
      A_mat    = coef_R_1 * R_inv_t + coef_R_2 * R_2
      
      # Step 2c: compute the score, using vec(A_mat)' * vec(B_mat) =
      #          trace(A_mat' * Bmat) = trace(A_mat * Bmat) given
      #the symmetry of A_mat. See (F.3) in the online appendix.
      g_f = 2 * f_t / (denom_t ^ 2)
      s_vec[j,1] = g_f * (sum(A_mat) - sum(diag(A_mat)))
      
      # Step 2d: update the time-varying parameter using the score for the next
      #          iteration
      if(j<T){
        f_vec[j+1] = omega_l +  A * s_vec[j,1] + B * f_vec[j,1]
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
