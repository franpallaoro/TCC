LogLik_Copula_LT_factor_given_omega = function(N,T,params,f_hat_vec,u_mat,asset_group_vec,n_vec,ind_t_dist,ind_Rt,ind_sim,...){
  G = max(asset_group_vec)
  g_vec_cum            = cumsum(seq(G, 1, by = -1))
  k                    = G
  p                    = G*(G+1)/2

  A_vec  = params[1]
  B_vec  = params[2]
  omega_vec = (1-B_vec)%*%t(f_hat_vec)

  if(ind_sim == 1){
    x_mat = u_mat
    if(ind_t_dist == 1){
      end = length(params)
      nu      = params[end]
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

  

  loglike_vec         = matrix(0, nrow = T, ncol = 1)            
  f_mat               = matrix(0, nrow = T, ncol = p)
  s_mat               = matrix(0, nrow = T,ncol = p)

  if(ind_Rt==1){
    R_mat   = array(0, dim = c(G, G, T)) 
  }else{
    R_mat   = vector() #nao sei 
  }


  f_mat[1,] = f_hat_vec


  teller_0 = 1
  asset_group_vec_new = matrix(0, nrow = N, ncol = 1)
  for(m in 1:G){
    asset_group_vec_new[teller_0:(teller_0+n_vec[m]-1)] = m
    teller_0                                          = teller_0 + n_vec[m]
  }


  S_l_mat    = Compute_S_l_matrix(asset_group_vec_new)[[1]] #função que retorna duas matrizes mas aqui ta pedindo só uma??

  Ind_lambda_to_group_vec = vector()
  for(i in 1:G){
    Ind_lambda_to_group_vec = c(Ind_lambda_to_group_vec, t(i:G))
  }
  Ind_lambda_to_group_mat      = cbind(1:p, Ind_lambda_to_group_vec)

  # we are now ready to build "dau_L_tilde_dau_lambda_tilde"
  dau_L_tilde_dau_lambda_tilde = matrix(0, nrow = N*k, ncol = p)

  ind_match_lambda_group_cell = matrix(list(), G, 1) #oq substituir


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

    for(i in 1:k){ # erro: i = 10, último dimensao da sl_ma
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
      dau_L_tilde_dau_f_mat =  dau_L_tilde_dau_lambda_tilde%*%dau_lambda_tilde_dau_f_mat;


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
        f_mat[j+1,] = omega_vec +  A_vec*s_mat[j,] + B_vec*f_mat[j,] 
      }

      # Step 2e: compute the core of the log likelihood for this observation,
      #          with additional constants added outside the loop, and the
      #          "marginal" part of the copula density also added outside the
      #          loop.                    
      if(ind_t_dist ==1){
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
    # complete likelihood
    if(ind_t_dist==1){

      loglike_vec = loglike_vec + log(gamma(0.5*(nu+N))) + (N-1)*log(gamma(nu/2)) - N * log(gamma(0.5*(nu+1))) + 0.5*(nu+1)*rowSums(log(matrix(1,T,N)+(1/(nu-2))*x_mat^2))

    }else{
      loglike_vec = loglike_vec + 0.5*rowSums(x_mat^2)
    }

    LLF = -sum(loglike_vec)
  }else{
    LLF = 1e14
  }

  return(list("LLF" = LLF, "loglike_vec" = loglike_vec, 
              "R_mat" = R_mat, "s_mat" = s_mat, "f_mat" = f_mat))
}








# function [LLF,loglike_vec,R_mat,s_mat,f_mat] = LogLik_Copula_LT_factor_given_omega(N,T,params,f_hat_vec,u_mat,asset_group_vec,n_vec,ind_t_dist,ind_Rt);
# 
# # MF_LT copula, (lower triangular matrix) g factors, see eq (6)
# # see line 12-28 of this function for the general input and output.
# 
# # Extra input:
#   # ind_Rt: 1 if you would like to compute and store the block correlation matrix, 0 otherwise
# # f_hat_vec: estimated f_bar using the first step moment
# #            estimator using eq(12).
# # n_vec    : G by 1 vector where entry g (g = 1,..,G) denotes
# # the number of assets in group g.
# 
# #################################################################################
# #################################################################################
# #################################################################################
# 
# # IMPORTANT NOTE:
#   # asset_group_vec and u_mat needs to be SORTED before running this function (!!!)
#     
#     #################################################################################
# #################################################################################
# #################################################################################
# 
# 
# # specify number of factors (k), number of groups (G) and the number
# # of unique factor loadings (p)            
# G = max(asset_group_vec);
# g_vec_cum            = cumsum(G:-1:1);
# k                    = G;
# p                    = G*(G+1)/2;
# 
# # Get out the parameters for the model from the parameter vector params
# # Order: A_vec, B_vec [, nu]  (nu only there if ind_t_dist == 1)           
# A_vec  = params(1);
# B_vec  = params(2);
# omega_vec = (1-B_vec)*f_hat_vec';
#             
#             # Compute the inv cdf (Gaussian or Student's t) of the PITs
# # as input for the copula density later on
# if ind_t_dist ==1
# nu      = params(end);
# x_mat_0 = tinv(u_mat,nu);
# x_mat   = sqrt((nu-2)/nu) * x_mat_0;
# else
#   x_mat   = norminv(u_mat);
# end
# 
# # Initialize the storage objects for the output
# loglike_vec         = zeros(T,1);            
# f_mat               = zeros(T,p);
# s_mat               = zeros(T,p);
# if ind_Rt==1
# R_mat   = zeros(G,G,T);
# else
#   R_mat   = [];
# end
# 
# # Compute the starting value of f(t)
# # For the MF_LT model, we set f(1,:)
# # f_bar, estimated by moment-based estimator, see eq(11) and eq(12)            
# f_mat(1,:) = f_hat_vec;
# 
# # Maka a new asset_group_vec that is ascending from 1 to G,
# # where 1 is your first sorted industry, 2 the second etc. 
# 
# teller_0 = 1;
# asset_group_vec_new = zeros(N,1);
# for m = 1:G
# asset_group_vec_new(teller_0:teller_0+n_vec(m)-1) = m;
# teller_0                                          = teller_0 + n_vec(m);
# end
# 
# # initialize fixed matrices required during computations
# S_l_mat    = Admin.Compute_S_l_matrix(asset_group_vec_new);
# 
# # Part of the score dau_c/dau_vec_f_t is build outside the time-loop
# # this is the selection matrix "dau_L_tilde_dau_lambda_tilde".
# 
# # Consider eq(6); we have G(G+1)/2 different f_values.
# # Since lambda_tilde_1,1,t and lambda_tilde_2,2,t are in the
# # same row, they will both appear in the non-linear scaling
# # formulas of (3).
# 
# # Here we assign which lambda_tilde values belong to the same
# # row of L_tilde_prime (or COLUMN of L_tilde)
# Ind_lambda_to_group_vec = [];
# for i = 1:G
# Ind_lambda_to_group_vec = [Ind_lambda_to_group_vec; [i:G]'];
# end                                   
# Ind_lambda_to_group_mat      = [[1:p]' Ind_lambda_to_group_vec];
# 
# # we are now ready to build "dau_L_tilde_dau_lambda_tilde"
# dau_L_tilde_dau_lambda_tilde = zeros(N*k,p);
# 
# 
# ind_match_lambda_group_cell = cell(G,1);
# 
# #     Building dau_L_tilde_dau_lambda_tilde consists of the following steps:            
#   #     1) Loop over each column i of L_tilde (i = 1:G)
# #     (note that lambda_tilde_i is the first element in each colomn)
# #     2) compute dau_L_t_lambda_tilde_i;
# #     3) check which lambda_tilde_j (j = 1,...,p) are in the same column i of L_tilde
# #        store these into a vector ind_lambda_group_i (and in a
#                                                        #        cell ind_match_lambda_group_cell)
# #     4) if i > 1 we have more than on lambda_tilde_j in the
# #        i-th column of L_tilde. More speficically, we have exactly (i-1) other
# #        lambda_tilde_j parameters.
# #        4a) loop over these (i-1) other lambda_tilde_j values (using the vector ind_lambda_g_group_i) 
# #            and compute dau_L_t_lambda_tilde_j
# #     5) return to step 1).
# 
# for i =  1:G
# # step 2
# dau_L_tilde_dau_lambda_tilde(sum(n_vec(1:i-1))*G+1:G:sum(n_vec(1:i))*G,i) =1;
# 
# # step 3
# ind_lambda_g_group_i           = Ind_lambda_to_group_mat(Ind_lambda_to_group_mat(:,2)==i,1);
# ind_match_lambda_group_cell{i} = ind_lambda_g_group_i;
# 
# # step 4
# for j = 1:i-1
# lambda_j                   = ind_lambda_g_group_i(j+1);
# dau_L_tilde_dau_lambda_tilde(sum(n_vec(1:i-1))*G+1+j:G:sum(n_vec(1:i))*G+j,lambda_j) =1;
# end
# end
# 
# 
# # Loop over all observations to compute the time varying
# # parameters f(t) and implied time-varying correlation
# # parameters.
# 
# 
# for j = 1:T
# 
# # Step 0: compute the time varying coefficient and the
# #         data vector at this time
# 
# f_t_vec   = f_mat(j,:);
# f_t_vec_2 = f_t_vec.^2;
# x_j       = x_mat(j,:)';
#                                 
#                 # Build the matrix L_til_prime_mat_t. See (F.19) of the web
#                 # appendix.                 
#                 # We also build a second matrix L_tilde_prime_mat_tg (F.19)
#                 # with N_g = 1 for all values of g)                 
#                 # This G x G matrix is based on only unique elements of f_j,t.
#                 # We (could) use this later on to build the G x G block correlation
#                 # matrix with within and between group correlations.
#                 
#                 teller_i        = 1;
#                 denom_f_un      = ones(G,1);
#                 denom_f_vec     = ones(N,1);
#                 f_prime_mat_t   = zeros(N,G);
#                 f_prime_mat_tg  = zeros(G,G);
#                 for i = 1:k
#                     g_i_c           = g_vec_cum(i);
#                     denom_f_vec     = denom_f_vec + S_l_mat(:,i:G)*f_t_vec_2(teller_i:g_i_c)';  # Based on (F.19) and used repeatedly. 
# denom_f_un      = denom_f_un + [zeros(i-1,1); f_t_vec_2(teller_i:g_i_c)'];  # This is (F.21) and used repeatedly.
#                     f_prime_mat_t(:,i)  = S_l_mat(:,i:G)*f_t_vec(teller_i:g_i_c)';
#                                 f_prime_mat_tg(:,i) = [zeros(i-1,1); f_t_vec(teller_i:g_i_c)'];
#                     teller_i            = g_i_c+1;
#                 end
#                 
#                 # Step 1b: compute the correlation matrix (if ind_Rt=1) and its
#                 #          inverse and determinant
#                 denom_f_mat             = denom_f_vec*ones(1,k);
#                 lambda_til_prime_mat_t  = f_prime_mat_t./ sqrt(denom_f_mat);
#                 sigma_2_vec_t           = (1./denom_f_vec);
#                 [R_inv_t,det_R_t,Flag] = Admin.Compute_R_t_inv_and_det_R_t_L_matrix(sigma_2_vec_t,lambda_til_prime_mat_t);
#                 if ind_Rt==1
#                     denom_f_mat_tg          = denom_f_un*ones(1,G);
#                     L_tilde_prime_mat_tg    = f_prime_mat_tg./sqrt(denom_f_mat_tg);
#                     Rt_Block                = L_tilde_prime_mat_tg* L_tilde_prime_mat_tg';
#                                                        R_mat(:,:,j) = Rt_Block;
#                                                        end
#                                                        
#                                                        # CHECKUP: if something wrong with the correlation matrix (not PDS),
#                                                        # then return likelihood failure by huge -loglik value.                
#                                                        
#                                                        if Flag==1 || isnan(det_R_t)
#                                                        loglike_vec = 1e14;
#                                                        break
#                                                        else
#                                                          #Step 2: Compute the likelihood for time t=j and store it,
#                                                        #         compute the score function and update the time-varying
#                                                        #         parameters
#                                                        
#                                                        # Compute the score is factor-model-specific.
#                                                        
#                                                        # Step 2a: compute R^{-1} x(t) and x(t)' R^{-1} x(t)                    
#                     coef_R_1    = -0.5;
#                     R_inv_x     = (R_inv_t * x_j);
#                     x_R_inv_x   =  x_j' * R_inv_x;
#                                                        
#                                                        # Step 2b: compute the derivative of the copula density with respect to
#                                                        #          R(t) as in (A8)-(A9)                    
#                                                        if ind_t_dist==1
#                                                        coef_R_2 = 0.5*(nu + N)/(nu-2+ x_R_inv_x );
#                                                        else
#                                                          coef_R_2 = 0.5;
#                                                        end
#                                                        
#                                                        R_2 = R_inv_x * R_inv_x';  # outer product: R^{-1} x(t) x(t)' R^{-1}
#                                                        
#                                                        # dau_log_c_dau_vec_R
#                                                        A_mat = coef_R_1 * R_inv_t + coef_R_2 *R_2;
#                                                        
#                                                        # Step 2c: compute the score, using  (F.2),(F.3) (F.20)-(F.25)
#                                                        # of the online appendix.
#                                                        
#                                                        # auxilary term used repeatedly 
#                                                        denom_f_un_3_2           = denom_f_un.^(3/2);
#                                                        
#                                                        # (F.22) for j = 1,2,..., G(G+1/2)
#                                                        dau_lambda_tilde_dau_f_vec  = 1./sqrt(denom_f_un(Ind_lambda_to_group_mat(:,2))) - f_t_vec_2'./denom_f_un_3_2(Ind_lambda_to_group_mat(:,2));
#                     dau_lambda_tilde_dau_f_mat  = diag(dau_lambda_tilde_dau_f_vec);
#                     
#                     dau_sigma2_dau_f_mat         = zeros(G,p);
#                     
#                     # fill dau_lambda_tilde_dau_f and dau_sigma2_dau_f
#                                         
#                     # The loop below computes (F.23) for all values of j, looping over all G columns. 
#                     # That is: check which values of f_(t,j) are in column g. Then compute the cross deratives.                     
#                      
#                     # The loop also computes (F.25) for each g (g = 1,...,G) 
#                     # That is: check which values of f_(t,j) are in column g. 
#                     # Then compute (F.25) w.r.t. to all these f_(t,j) values.                      
#                     for g = 1:G
#                         # find out which f_(t,j) are in column g
#                         ind_cross_lambda_g =  ind_match_lambda_group_cell{g};
#                         if g>1
#                             for m = 1:g
#                                 ind_2   = ind_cross_lambda_g;
#                                 ind_m1  = ind_cross_lambda_g(m);
#                                 ind_2(m)= [];
#                                 
#                                 #(F.23)
#                                 dau_lambda_tilde_i_dau_m_vec             = -(f_t_vec(ind_m1).*f_t_vec(ind_2))./denom_f_un_3_2(g);
#                                 dau_lambda_tilde_dau_f_mat(ind_m1,ind_2) = dau_lambda_tilde_i_dau_m_vec;
#                             end
#                         end
#                         #(F.25)
#                         dau_sigma2_dau_f_mat(g,ind_cross_lambda_g) = -2*f_t_vec(ind_cross_lambda_g)./ (denom_f_un(g)^2);
#                     end
#                     
#                     #(F.20)
#                     dau_L_tilde_dau_f_mat =  dau_L_tilde_dau_lambda_tilde * dau_lambda_tilde_dau_f_mat;
#                     
#                     
#                     # (F.24)                                         
#                     dau_D_dau_f_mat = S_l_mat *  dau_sigma2_dau_f_mat;
#                     
#                     # combining above with (F.2) and (F.3)
#                     L_A         = lambda_til_prime_mat_t' * A_mat;
#                                                        s_vec_a     = 2* (L_A(:))'* dau_L_tilde_dau_f_mat;
#                     s_vec_b     = diag(A_mat)' * dau_D_dau_f_mat;
#                                                        s_mat(j,:)  = s_vec_a + s_vec_b;
#                                                        
#                                                        # Step 2d: update the time-varying parameters using the score for the next
#                                                        #          iteration
#                                                        if j<T
#                                                        f_mat(j+1,:) = omega_vec +  A_vec.*s_mat(j,:) + B_vec.*f_mat(j,:);
#                                                        end
#                                                        
#                                                        # Step 2e: compute the core of the log likelihood for this observation,
#                                                        #          with additional constants added outside the loop, and the
#                                                        #          "marginal" part of the copula density also added outside the
#                                                        #          loop.                    
#                                                        if ind_t_dist ==1
#                                                        # Student's t copula core
#                         loglike_vec(j,1) = -0.5 *(nu+N)*log(1 + x_R_inv_x/(nu-2)) - 0.5*log(det_R_t);
#                     else
#                         # Gaussian copula core
#                         loglike_vec(j,1) = -0.5 * x_R_inv_x - 0.5*log(det_R_t);
#                     end
#                 end
#                 
#                 
#             end
#             # CHECKUP: if the likelihood computations were successful, proceed
#             # and return the likelihood value; otherwise, return a failure by
#             # returning a high -likelihood value.
#             
#             
#             if isreal(loglike_vec) && sum(loglike_vec) ~=1e14
#                 # complete likelihood
#                 if ind_t_dist==1
#                     # Student's t case: add the marginal Student t density parts
#                                                        # and the integrating constants                    
#                                                        loglike_vec = loglike_vec + gammaln(0.5*(nu+N)) + (N-1)*gammaln(nu/2) - N * gammaln(0.5*(nu+1)) ...
#                                                        + 0.5*(nu+1)*sum(log(ones(T,N)+(1/(nu-2))*x_mat.^2),2);
#                                                        else
#                                                          # Gaussian case: add the marginal Gaussian density parts
#                                                        loglike_vec = loglike_vec + 0.5*sum(x_mat.^2,2);
#                                                        end
#                                                        # store the -log-likelihood value
#                                                        LLF = -sum(loglike_vec);
#                                                        else
#                                                          LLF = 1e14;
#                                                        end
#                                                        
#                                                        
#                                                        end        
#                                                        