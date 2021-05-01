Compute_R_t = function(sigma_2_vec_t,L_til_prime_t){
  #Compute R_t following eq (4) of Opschoor et al. 2020

  R_t = L_til_prime_t%*%t(L_til_prime_t) + diag(sigma_2_vec_t)
  return(R_t)        
}