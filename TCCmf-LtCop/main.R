#import data and gather T and N
u_mat      = read.csv('Testdata.csv', header = F)
u_mat[,1]  = NULL
u_mat      = as.matrix(u_mat)
T          = dim(u_mat)[1]
N          = dim(u_mat)[2]

# put to if you would like to compute standard errors based on the sandwich
# covariance matrix
ind_se     = 1

# information about group size and group allocation
n_vec           = 10*matrix(1,10,1)            # number of assets per group
asset_group_vec = kronecker(t(1:10),matrix(1,10,1))
G               = max(asset_group_vec)     #number of groups


##########################################################################
##########################################################################
############### Matlab routine to set up constrained optimization ########

#install.packages("pracma")
#install.packages("NlcOptim")
library(pracma)
library(NlcOptim)


#options = optimset('fmincon');
#options = optimset(options, 'Display', 'off');
#options = optimset(options, 'Diagnostics', 'off');
#options = optimset(options, 'Algorithm', 'sqp');
MaxIter     = 30000
MaxFunEvals = 30000

# Create starting values and lower/upper bounds for the contrained optimization

#para_start_1f_eq  = [0.001 0.01 0.92 30]';
#para_start_1f_gr  = [0.001*ones(1,G) 0.02 0.9 30]';
#para_start_2f = [0.001*ones(1,G) 0.02 0.01*ones(1,2) 0.9 30]';
para_start_Mf = c(0.05, 0.02*matrix(1,1,G), 0.01*matrix(1,1,2), 0.9, 30)

# lower and upper bounds
#Lb_1f_eq = [-1 eps*ones(1,2) 2.1]';
#Rb_1f_eq = [ 0.1 0.10 0.99999 5000]';
#Lb_1f_gr = [-1*ones(1,G) eps*ones(1,2) 2.1]';
#Rb_1f_gr = [0.1*ones(1,G) 0.10 0.99999 5000]';
#Lb_2f = [-1*ones(1,G+1) eps*ones(1,3) 2.1]';
#Rb_2f = [0.1*ones(1,G+1) 0.1 0.1 0.99999 5000]';
Lb_Mf = c(.Machine$double.eps*matrix(1, 1, 4+10), 2.1) #[ eps*ones(1,4+G) 2.1]';
Rb_Mf = c(0.5*matrix(1,1,G+1), 0.2, 0.2, 0.99999, 5000)

nr_l  = G*(G+1)/2
ind_pos =  c(1, cumsum(seq(G, 2, by = -1)) + 1)#[1; cumsum(G:-1:2)' + 1];

# Lower and upper bounds for the first step moment estimator of MF_LT
Lb_Mf_LT_1 = -3*matrix(1,nr_l,1)
Lb_Mf_LT_1[ind_pos] = .Machine$double.eps
Rb_Mf_LT_1 =  3*matrix(1,nr_l,1)
para_start_Mf_LT_1 = 0.2 * matrix(1,nr_l,1)

# Lower and upper bounds for the second step ML estimator of MF_LT
Lb_Mf_LT_2 = c(.Machine$double.eps*matrix(1,1,2), 2.1)
Rb_Mf_LT_2 = c(5*matrix(1,1,1), 0.99999, 5000)
para_start_Mf_LT_2 = c(0.01, 0.96, 20)
            
#create a empty cell array for parameters, standard errors (optional) and strings (will be used later on to print the results)
# para_cell = cell(10,1);
# if ind_se==1
# se_para_cell = cell(10,1);
# end
# para_cell_str = cell(10,1);
# #para_1f_eq_str = {'omega', 'A','B', 'nu Cop'};
# #para_1f_gr_str = {'omega_1','omega_2','omega_3','omega_4','omega_5','omega_6','omega_7','omega_8','omega_9','omega_10', 'A','B', 'nu Cop'};
# #para_2f_str = {'omega_eq','omega_1','omega_2','omega_3','omega_4','omega_5','omega_6','omega_7','omega_8','omega_9','omega_10', 'A_eq','A_gr','B', 'nu Cop'};
# para_Mf_str = para_2f_str;
# para_Mf_LT_str = {'A','B', 'nu Cop'};
#             
# para_cell_str{1} = para_1f_eq_str(1:end-1);
# para_cell_str{2} = para_1f_eq_str;
# para_cell_str{3} = para_1f_gr_str(1:end-1);
# para_cell_str{4} = para_1f_gr_str;
# para_cell_str{5} = para_2f_str(1:end-1);
# para_cell_str{6} = para_2f_str;
# para_cell_str{7} = para_Mf_str(1:end-1);
# para_cell_str{8} = para_Mf_str;
# para_cell_str{9} = para_Mf_LT_str(1:end-1);
# para_cell_str{10} = para_Mf_LT_str;
#             

# step 1: moment estimator of f_bar (see eq 11/12 of Opschoor et al. (2020))
# NOTE: DATA IS ALREADY BEEN SORTED ACROSS INDUSTRY (from 1 to ..... 10)
# (Note also that the moment estimator is the same as in the Gaussian case)
# Compute matrix R_m of eq(12)

x_mat    = qnorm(u_mat)#norminv(u_mat);
Rt_Block = Compute_BLOCK_correlation_matrix(x_mat,n_vec)[[1]]
Rt_sample_vech = Sym2Vech(G,t(Rt_Block))
            

Moment_based_omega_LT_FC_Optim <- function(G, params, rho_vec_sample_bl_vech) {
  Moment_based_omega_LT_FC(G = G, params = params, rho_vec_sample_bl_vech = Rt_sample_vech)[[1]]
}
step1 = fmincon(x0 = para_start_Mf_LT_1, fn = Moment_based_omega_LT_FC_Optim, 
                lb = Lb_Mf_LT_1, ub = Rb_Mf_LT_1, maxfeval = MaxFunEvals, maxiter = MaxIter,
                G = G, rho_vec_sample_bl_vech = Rt_sample_vech)
f_bar = step1$par

# step 2: ML (t copula) given f_bar
LogLik_Copula_LT_factor_given_omega_Optim <- function(N,T,params,f_hat_vec,u_mat,asset_group_vec,n_vec,ind_t_dist,ind_Rt){
  LogLik_Copula_LT_factor_given_omega(N,T,params,f_bar,u_mat,asset_group_vec,n_vec,1,0)[[1]]
}
step2 = fmincon(x0 = para_start_Mf_LT_2, fn = LogLik_Copula_LT_factor_given_omega_Optim, 
                lb = Lb_Mf_LT_2, ub = Rb_Mf_LT_2, maxfeval = MaxFunEvals, maxiter = MaxIter,
                N = N, T = T, f_hat_vec = f_bar, u_mat = u_mat, asset_group_vec = asset_group_vec, 
                n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0)



# used to print the results
theta_opt = step2[[1]]
max_LogLik = -step2[[2]]
            