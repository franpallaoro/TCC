Compute_BLOCK_correlation_matrix = function(x_mat, n_vec){ #function [Rho_block_mat,Rt_Block] = Compute_BLOCK_correlation_matrix(x_mat, n_vec)
# build a G times G block correlatoin matrix where
# entry (i,j) (i not equal to j) means the average beween
# correlation between industry i and j
# and entry (i,i) indicates the average WITHIN correlation of
# industry i.

T = dim(x_mat)[1] 
N = dim(x_mat)[1] #[T,N] = size(x_mat);
G = length(n_vec) #length(n_vec);
Rt_mat = cor(x_mat) #corr(x_mat);

teller_i = 1 #;
Rho_block_mat = matrix(0, ncol = G, nrow = G) #zeros(G,G);
Rt_0 = matrix(0, ncol = N, nrow = N) #zeros(N,N);
Rt_1 = matrix(0, ncol = N, nrow = N) #zeros(N,N);
# loop over all G industries and compute the average WITHIN and BETWEEN
# correlatinos.
for (i in 1:G) {
  n_i = n_vec[i]
  Rt_mat_ii = Rt_mat[teller_i:(teller_i+n_i-1),teller_i:(teller_i+n_i-1)]
  if (n_i>1){
    rho_ii = as.numeric((t(matrix(1, n_i, 1))%*%Rt_mat_ii%*%matrix(1, n_i, 1) - n_i)/(n_i*(n_i-1)))
    
  }else{
    rho_ii = 0
  }
  Rho_block_mat[i,i] = rho_ii
  
  Rt_0[teller_i:(teller_i+n_i-1), teller_i:(teller_i+n_i-1)] = (1- rho_ii)*diag(n_i)
  
  Rt_1[teller_i:(teller_i+n_i-1), teller_i:(teller_i+n_i-1)] = rho_ii
  
  teller_j = sum(n_vec[1:i])+1
  if(i != G){
    for (b in (i+1):G) {
      Rt_mat_ij =  Rt_mat[teller_i:(teller_i+n_i-1), teller_j:(teller_j+n_vec[b]-1)]
      rho_ij = mean(Rt_mat_ij)
      Rho_block_mat[i, b] = rho_ij
      Rho_block_mat[b,i] = rho_ij
      Rt_1[teller_i:(teller_i+n_i-1), teller_j:(teller_j+n_vec[b]-1)] = rho_ij
      Rt_1[teller_j:(teller_j+n_vec[b]-1), teller_i:(teller_i+n_i-1)] = rho_ij
      
      teller_j = teller_j+n_vec[b]
    }
  }

  teller_i = teller_i+n_vec[i]
}

Rt_Block  = Rt_0 + Rt_1

return(list("Rho_block_mat" = Rho_block_mat, "Rt_Block" = Rt_Block))
}
# for i = 1:G
# n_i = n_vec(i);
# Rt_mat_ii = Rt_mat(teller_i:teller_i+n_i-1,teller_i:teller_i+n_i-1);
# if n_i>1
# rho_ii = (ones(n_i,1)'*Rt_mat_ii*ones(n_i,1) - n_i)/(n_i*(n_i-1));
#                 else
#                     rho_ii = 0;
#                 end
#                 Rho_block_mat(i,i) = rho_ii;
#                 
#                 Rt_0(teller_i:teller_i+n_i-1,teller_i:teller_i+n_i-1) = (1- rho_ii)*eye(n_i);
#                 Rt_1(teller_i:teller_i+n_i-1,teller_i:teller_i+n_i-1) = rho_ii;
#                 
#                 teller_j = sum(n_vec(1:i))+1;
#                 for b = i+1:G
#                     Rt_mat_ij = Rt_mat(teller_i:teller_i+n_i-1,teller_j:teller_j+n_vec(b)-1);
#                     rho_ij = mean(Rt_mat_ij(:));
#                     Rho_block_mat(i,b) = rho_ij;
#                     Rho_block_mat(b,i) = rho_ij;
#                     
#                     Rt_1(teller_i:teller_i+n_i-1,teller_j:teller_j+n_vec(b)-1) = rho_ij;
#                     Rt_1(teller_j:teller_j+n_vec(b)-1,teller_i:teller_i+n_i-1) = rho_ij;
#                     
#                     teller_j = teller_j+n_vec(b);
#                 end
#                 
#                 teller_i = teller_i+n_vec(i);
#             end
#             
#             Rt_Block  = Rt_0 + Rt_1;
#         end