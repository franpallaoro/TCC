Compute_R_t_inv_and_det_R_t_L_matrix <- function(sigma_2_vec_t,L_til_prime_t){
# Compute inverse of R_t and det R_t using eq (5) of Opschoor et al (2020)

  Flag    = 0
  if(is.vector(L_til_prime_t)){
    m = length(L_til_prime_t)
  }else{
    m = ncol(L_til_prime_t)
  }#size(L_til_prime_t,2)

  sigma_2_inv           = 1/sigma_2_vec_t
  L_til_t_times_D_inv_t = L_til_prime_t*(sigma_2_inv%*%matrix(1,1,m))
  temp_0                = diag(m) +  t(L_til_t_times_D_inv_t)%*%L_til_prime_t
            
  ind_pd = try(chol(temp_0), TRUE)
  if(!is.matrix(ind_pd)){
    Flag = 1
    R_inv_t = NaN
    det_R_t = NaN
  }else{
    #D_inv_t = diag(1/sigma_2_vec_t)
    D_inv_t = diag(length(sigma_2_vec_t))
    diag(D_inv_t) = 1/sigma_2_vec_t
    R_inv_t = D_inv_t - (L_til_t_times_D_inv_t%*%solve(temp_0))%*%t(L_til_t_times_D_inv_t)
    det_D_t = prod(sigma_2_vec_t)
    det_R_t = det(temp_0)*det_D_t
  }
  
  return(list("R_inv_t" = R_inv_t, "det_R_t" = det_R_t, "Flag" = Flag))
}             
