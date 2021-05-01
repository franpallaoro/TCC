Compute_R_t_inv_and_det_R_t = function(sigma_2_vec_t,L_til_prime_t){
  # Compute inverse of R_t and det R_t using eq (5) of Opschoor et al (2020)

  sigma_2_inv           = 1/sigma_2_vec_t
  L_til_t_times_D_inv_t = L_til_prime_t*sigma_2_inv
  temp_0                = 1 + t(L_til_t_times_D_inv_t)%*%L_til_prime_t
  D_inv_t               = diag(sigma_2_inv)
            
  R_inv_t = D_inv_t - (1/temp_0)*(t(L_til_t_times_D_inv_t)%*%L_til_t_times_D_inv_t)
  det_D_t =  prod(sigma_2_vec_t)
  det_R_t = temp_0*det_D_t
  
  return(list(R_inv_t,det_R_t,temp_0,det_D_t))
}
