rm(list = ls())

a <- Sys.time()

source("Compute_BLOCK_correlation_matrix.R")
source("Compute_R_t_inv_and_det_R_t_L_matrix.R")
source("Compute_S_l_matrix.R")
source("LogLik_Copula_LT_factor_given_omega.R")
source("Moment_based_omega_LT_FC.R")
source("Sim_u_mat.R")
source("Sym2Vech.R")
source("LogLik_Copula_1f_eq.R")
source("Compute_R_t.R")
source("Compute_R_t_inv_and_det_R_t.R")
source("Sim_x_mat_1f_eq.R")

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




#install.packages("pracma")
#install.packages("NlcOptim")
library(pracma)
library(NlcOptim)
library(nloptr)


#options = optimset('fmincon');
#options = optimset(options, 'Display', 'off');
#options = optimset(options, 'Diagnostics', 'off');
#options = optimset(options, 'Algorithm', 'sqp');
MaxIter     = 30000
MaxFunEvals = 30000

# Create starting values and lower/upper bounds for the contrained optimization

para_start_Mf = c(0.05, 0.02*matrix(1,1,G), 0.01*matrix(1,1,2), 0.9, 30)

# lower and upper bounds

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
#Lb_Mf_LT_2 = c(.Machine$double.eps*matrix(1,1,2), 2.1)#certo
Rb_Mf_LT_2 = c(5*matrix(1,1,1), 0.99999, 5000)
#Rb_Mf_LT_2 = c(5*matrix(1,1,1), 0.99999, 5000)#certo
#para_start_Mf_LT_2 = c(0, 0.9, 30)
para_start_Mf_LT_2 = c(0.01, 0.96, 20)#certo





# step 1: moment estimator of f_bar (see eq 11/12 of Opschoor et al. (2020))
# NOTE: DATA IS ALREADY BEEN SORTED ACROSS INDUSTRY (from 1 to ..... 10)
# (Note also that the moment estimator is the same as in the Gaussian case)
# Compute matrix R_m of eq(12)

x_mat    = qnorm(u_mat)#norminv(u_mat);
Rt_Block = Compute_BLOCK_correlation_matrix(x_mat,n_vec)[[1]]
Rt_sample_vech = Sym2Vech(G,t(Rt_Block))


# Moment_based_omega_LT_FC_Optim <- function(G, params, rho_vec_sample_bl_vech) {
#   Moment_based_omega_LT_FC(G = G, params = params, rho_vec_sample_bl_vech = Rt_sample_vech)[[1]]
# }
step1 = slsqp(x0 = para_start_Mf_LT_1, fn = Moment_based_omega_LT_FC, 
                lower = Lb_Mf_LT_1, upper = Rb_Mf_LT_1,
                G = G, rho_vec_sample_bl_vech = Rt_sample_vech)
f_bar = step1$par

# step 2: ML (t copula) given f_bar
# LogLik_Copula_LT_factor_given_omega_Optim <- function(N,T,params,f_hat_vec,u_mat,asset_group_vec,n_vec,ind_t_dist,ind_Rt,ind_sim){
#   LogLik_Copula_LT_factor_given_omega(N,T,params,f_bar,u_mat,asset_group_vec,n_vec,1,0,0)[[1]]
# }
step2 = slsqp(x0 = para_start_Mf_LT_2, fn = LogLik_Copula_LT_factor_given_omega, 
                lower = Lb_Mf_LT_2, upper = Rb_Mf_LT_2,
                N = N, T = T, f_hat_vec = f_bar, u_mat = u_mat, 
                n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0, ind_sim = 0)
step2 = fmincon(x0 = para_start_Mf_LT_2,fn = LogLik_Copula_LT_factor_given_omega,
        lb = Lb_Mf_LT_2, ub = Rb_Mf_LT_2, N = N, T = T, f_hat_vec = f_bar, u_mat = u_mat, 
        n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0, ind_sim = 0)

theta_opt = step2$par

# used to print the results

#max_LogLik = -step2[[2]]

u_mat_sim = Sim_u_mat(N, G, T, theta_opt, f_bar, n_vec, 1, 0)


b <- Sys.time()

print(b - a)

#simulação 
N = 100
G = 10
T = 1000
A = 0.015
B = 0.97
nu = 35
omega = seq(-0.10,0.9, length.out = 55)
#omegas = (1-B_vec)%*%t(f_hat_vec)
theta = c(A, B, nu)
n_vec           = matrix(10,10,1)  
x_mat_sim = Sim_x_mat(N, G, T, theta, n_vec, 1, 0, omega)
       
x_mat_sim = as.matrix(x_mat_sim)

T          = dim(x_mat_sim)[1]
N          = dim(x_mat_sim)[2]

n_vec           = matrix(10,10,1)            
asset_group_vec = kronecker(t(1:10),matrix(1,10,1))
G               = max(asset_group_vec)    


Rt_Block_sim = Compute_BLOCK_correlation_matrix(x_mat_sim,n_vec)[[1]]
Rt_sample_vech_sim = Sym2Vech(G,t(Rt_Block_sim))



step1_sim = slsqp(x0 = para_start_Mf_LT_1, fn = Moment_based_omega_LT_FC, 
                lower = Lb_Mf_LT_1, upper = Rb_Mf_LT_1,
                G = G, rho_vec_sample_bl_vech = Rt_sample_vech_sim)
f_bar_sim = step1_sim$par



step2_sim = slsqp(x0 = para_start_Mf_LT_2, fn = LogLik_Copula_LT_factor_given_omega, 
                lower = Lb_Mf_LT_2, upper = Rb_Mf_LT_2,
                N = N, T = T, f_hat_vec = f_bar_sim, u_mat = x_mat_sim, asset_group_vec = asset_group_vec, 
                n_vec = n_vec, ind_t_dist = 1, ind_Rt = 0, ind_sim = 1)

step2_sim$par

# copula 1 fator 
para_start_1f_eq  = c(0.001, 0.01, 0.92, 30)
Lb_1f_eq = c(-1, rep(.Machine$double.eps, 2), 2.1)
Rb_1f_eq = c(0.1, 0.10, 0.99999, 5000)



LogLik_Copula_1f_eq_optim = slsqp(x0 = para_start_1f_eq, fn = LogLik_Copula_1f_eq, 
                                  lower = Lb_1f_eq, upper = Rb_1f_eq, N = N, 
                                  T = T, u_mat = u_mat, ind_t_dist = 1, ind_sim = 0)
T = 1000
N = 100
ns = 500
par = matrix(0, ns, 4)
for (i in 1:ns) {
  x_mat_sim = Sim_x_mat_1f_eq(N = N, T = T, c(0.07, 0.0085, 0.92, 35), ind_t_dist = 1)
  aux = nloptr(x0 = para_start_1f_eq, eval_f  = LogLik_Copula_1f_eq, 
               lb = Lb_1f_eq, ub = Rb_1f_eq, 
               N = N, T = T, u_mat = x_mat_sim, ind_t_dist = 1, ind_sim = 0,
               opts = list("algorithm"="NLOPT_LN_COBYLA",
                           "xtol_rel"=1.0e-8))
  par[i,] = aux$solution
  print(i)
}
para_start_1f_eq  = c(0.001, 0.01, 0.92, 30)
Lb_1f_eq = c(-1, rep(.Machine$double.eps, 2), 2.1)
Rb_1f_eq = c(0.1, 0.10, 0.99999, 5000)

Lb_1f_eq = c(0,0,0,20)
Rb_1f_eq = c(0.10000,   0.10000,    0.99999, 100)
para_start_1f_eq = c(0.05, 0.01, 0.92, 34)
x_mat_sim = Sim_x_mat_1f_eq(N = N, T = T, c(0.07, 0.0085, 0.92, 35), ind_t_dist = 1)
aux = nloptr(x0 = para_start_1f_eq, eval_f  = LogLik_Copula_1f_eq, 
             lb = Lb_1f_eq, ub = Rb_1f_eq, 
             N = N, T = T, u_mat = x_mat_sim, ind_t_dist = 1, ind_sim = 0,
             opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND",
                         "xtol_rel"=1.0e-8))
print(aux)

LogLik_Copula_1f_eq_optim_sim = slsqp(x0 = par1, fn = LogLik_Copula_1f_eq, 
                                      lower = Lb_1f_eq, upper = Rb_1f_eq, N = N, 
                                      T = T, u_mat = x_mat_sim, ind_t_dist = 1, ind_sim = 0,
                                      control = list(maxeval = 3000, xtol_rel= 1e-08))

par1 = c(LogLik_Copula_1f_eq_optim_sim$par[1:3]+runif(3,-0.01,0.01), LogLik_Copula_1f_eq_optim_sim$par[4]+runif(1,1,2))

aux = optim(par = para_start_1f_eq, fn = LogLik_Copula_1f_eq, 
            lower = Lb_1f_eq, upper = Rb_1f_eq, method = "L-BFGS-B",
            N = N, T = T, u_mat = x_mat_sim, ind_t_dist = 1, ind_sim = 0)
x = aux$par + runif(4)

aux = optim(par = x, fn = LogLik_Copula_1f_eq, 
            lower = Lb_1f_eq, upper = Rb_1f_eq, method = "L-BFGS-B",
            N = N, T = T, u_mat = x_mat_sim, ind_t_dist = 1, ind_sim = 0)

