source("Compute_BLOCK_correlation_matrix.R")
source("Compute_R_t_inv_and_det_R_t_L_matrix.R")
source("Compute_S_l_matrix.R")
source("LogLik_Copula_LT_factor_given_omega.R")
source("Moment_based_omega_LT_FC.R")
source("Sim_u_mat.R")
source("Sym2Vech.R")
source("estGAS.R")
source("estCop.R")
source("estSample.R")
source("estGarch.R")
source("returnP.R")
source("minVar.R")
source("getData.R")

# 1. (3)Comunicações: OIBR4 TIMS3 VIVT3
# 2. (15)Consumo cíclico: GRND3 LAME4 BTOW3 COGN3 CVCB3 CYRE3 EZTC3 JHSF3 LCAM3 LREN3 MGLU3 MRVE3 RENT3 VVAR3 YDUQ3
# 3. (5)Consumo não cíclico: ABEV3 BRFS3 BEEF3 JBSS3 MRFG3 
# 4. (6)Bens industriais: CCRO3 EMBR3 GOLL4 POMO4 WEGE3 ECOR3
# 5. (14)Financeiro: BBAS3 BBDC3 BBDC4 ITSA4 ITUB4 PSSA3 B3SA3 BBSE3 BRML3 CIEL3 MULT3 SANB11 SULA11
# 6. (11)Materiais básicos: BRAP4 BRKM5 CSNA3 DTEX3 GGBR3 GGBR4 GOAU4 KLBN4 UNIP6 USIM5 VALE3 
# 7. (16)Utilidade pública: CESP6 CMIG3 CMIG4 CPLE6 EGIE3 ELET3 ELET6 LIGT3 SBSP3 TRPL4 CPFE3 ENBR3 ENEV3 ENGI11 EQTL3 TAEE11
# 8. (5)Petróleo, gás e biocombustíveis: PETR3 PETR4 UGPA3 CSAN3 PRIO3
# 9. (4)Saúde: FLRY3 HYPE3 QUAL3 RADL3


data = getData('2014-01-01', '2019-12-31')
data$HGTX3.SA <- NULL
data$PCAR3.SA <- NULL
data$SUZB3.SA <- NULL
window = 1000
asset_group_vec = c(rep(1,3), rep(2,15), rep(3,5),
                    rep(4,6), rep(5,14), rep(6,11),
                    rep(7,16), rep(8,5), rep(9,4)) #vetor dos grupos
n_vec = c(3, 15, 5, 6, 14, 11, 16, 5, 4)

l = nrow(data) - window
T = nrow(data)
N = ncol(data)
Tns = 20000 #simulações

weightsCop_cov_u = matrix(0, l, N)#matriz guarda os pesos
weightsCop_cov_x = matrix(0, l, N)
weightsCop_cor_u = matrix(0, l, N)
weightsCop_cor_x = matrix(0, l, N)
weightsSample    = matrix(0, l, N)
returnCop_cov_u  = rep(0, l)
returnCop_cov_x  = rep(0, l)
returnCop_cor_u  = rep(0, l)
returnCop_cor_x  = rep(0, l)
returnSample     = rep(0, l)

u_mat_list     = list()
cov_mat_cop_u  = list()
cor_mat_cop_u  = list()
cov_mat_cop_x  = list()
cor_mat_cop_x  = list()
params_cop     = list()
cov_mat_sample = list()


count = 1
##############
for (i in 1:200) {
  #fran de 1:200
  #juju de 201:300
  #hudson 301:l
  a = Sys.time()
  i = i
  j = i + window-1
  
  u_mat = matrix(0, window, N)#matriz PITs
  var_vec_garch = rep(0, N) #variancia modelos marginais
  #estimação dos modelos marginais --------------------
  for (k in 1:N) {
    aux            = estGarch(data, i, j, k)
    u_mat[,k]      = aux[[1]]
    var_vec_garch[k] = aux[[2]]
  }
  
  #estimação da cópula --------------------
  u_mat_list[[i]]        = u_mat
  auxCop                 = estCop(u_mat, asset_group_vec, n_vec, Tns, 1)
  params_cop[[i]]        = auxCop[[5]]
  
   ##u_mat
  cov_mat_cop_u[[i]]       = auxCop[[1]]
  diag(cov_mat_cop_u[[i]]) = var_vec_garch
  cor_mat_cop_u[[i]]       = auxCop[[2]]
  diag(cor_mat_cop_u[[i]]) = var_vec_garch
  
  
   ##x_mat
  cov_mat_cop_x[[i]]       = auxCop[[3]]
  diag(cov_mat_cop_x[[i]]) = var_vec_garch
  cor_mat_cop_x[[i]]       = auxCop[[4]]
  diag(cor_mat_cop_x[[i]]) = var_vec_garch
  
   ##amostral
  cov_mat_sample[[i]] = cov(data[i:j,])
  
  #minima variancia
  weightsSample[i,] = minVar(cov_mat_sample[[i]])
  
  weightsCop_cov_u[i,]    = minVar(cov_mat_cop_u[[i]])
  weightsCop_cov_x[i,]    = minVar(cov_mat_cop_x[[i]])
  
  weightsCop_cor_u[i,]    = minVar(cor_mat_cop_u[[i]])
  weightsCop_cor_x[i,]    = minVar(cor_mat_cop_x[[i]])
  
  #retornos 
  returnSample[i] = returnP(wOptim = weightsSample[i,], data = data, j = j)
  
  returnCop_cov_u[i]    = returnP(wOptim = weightsCop_cov_u[i,], data = data, j = j)
  returnCop_cov_x[i]    = returnP(wOptim = weightsCop_cov_x[i,], data = data, j = j)
  
  returnCop_cor_u[i]    = returnP(wOptim = weightsCop_cor_u[i,], data = data, j = j)
  returnCop_cor_x[i]    = returnP(wOptim = weightsCop_cor_x[i,], data = data, j = j)
  
  
  b = Sys.time()
  
  print(count)
  count = count + 1
  print(b-a)
}
 