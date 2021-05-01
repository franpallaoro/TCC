LogLik_Copula_1f_eq = function(params, N, T, u_mat, ind_t_dist, ind_sim){

# 1-factor Copula: 1 equi loading with Normal or t_distribution,
# see line 12-28 of this function for the general input and output.

# Get out the parameters for the model from the parameter vector params
# Order: omega_l, A, B [, nu]  (nu only there if ind_t_dist == 1)
  omega_l = params[1]
  A    = params[2]
  B    = params[3]

# Compute the inv cdf (Gaussian or Student's t) of the PITs
            # as input for the copula density later on
  if(ind_sim == 1){
    x_mat = u_mat
    if(ind_t_dist == 1){
      end = length(params)
      nu      = params[end]
    }
  }else{
    if (ind_t_dist == 1){
      end = length(params)
      nu      = params[end]
      x_mat_0 = qt(u_mat, nu) 
      x_mat   = sqrt((nu-2)/nu) * x_mat_0
    }else{
      x_mat   = qnorm(u_mat)
    }
  }
    

                       
  #Initialize the storage objects for the output
  loglike_vec = matrix(0, T, 1)
  f_vec       = matrix(0, T, 1)
  s_vec       = matrix(0, T, 1)
  R_vec       = matrix(0, T, 1)
                       
  #initialize fixed matrices required during computations
  i_vec = matrix(1, N, 1)
                       
  # Compute the starting value of f(t)
  #For the one-factor-equi model, we target the average
  # sample correlation and then use the parameterization
   # of the loading as rho = f / sqrt(1+f*f)..
  RC_eq       = 1/(N*(N-1))*(t(i_vec)%*%cor(x_mat)%*%i_vec-N)
  f_vec[1]    = sqrt(RC_eq/(1-RC_eq))
                        
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
      
      # Step 2e: compute the core of the log likelihood for this observation,
      #          with additional constants added outside the loop, and the
      #          "marginal" part of the copula density also added outside the
      #          loop.
      if(ind_t_dist == 1){
        # Student's t copula core
        loglike_vec[j,1] = -0.5 *(nu+N)*log(1 + x_R_inv_x/(nu-2)) - 0.5*log(det_R_t)
      }else{
        # Gaussian copula core
        loglike_vec[j,1] = -0.5 * x_R_inv_x - 0.5*log(det_R_t)
      }
    }
  }
                
  # CHECKUP: if the likelihood computations were successful, proceed
  # and return the likelihood value; otherwise, return a failure by
  # returning a high -likelihood value.
  if(!is.complex(loglike_vec) && sum(loglike_vec) != 1e14){
    # complete the likelihood
    if(ind_t_dist==1){
      # Student's t case: add the marginal Student t density parts
        # and the integrating constants
      
      loglike_vec = loglike_vec + log(gamma(0.5*(nu+N))) + (N-1)*log(gamma(nu/2)) - N * log(gamma(0.5*(nu+1))) + 0.5*(nu+1)*rowSums(log(matrix(1,T,N)+(1/(nu-2))*x_mat^2))
    }else{
      # Gaussian case: add the marginal Gaussian density parts
      loglike_vec = loglike_vec + 0.5*rowSums(x_mat^2)
      }
    
    # store the -log-likelihood value
    LLF = -sum(loglike_vec)                   
                       
  }else{
    LLF = 1e14
  }
  return(LLF)
  
}                                                  
                                                  
                                                  
                                                  
                                                  
                                                  
                                                    
                                                 
                                                  
                                                  
                                                  
                                                  
                                                 
                       
                       
                       
                       
            
            