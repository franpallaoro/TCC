minVar = function(cov){
  if(!require(slam)){install.packages("slam")}
  cov_ret_1 = solve(cov)
  weigths = slam::row_sums(cov_ret_1) / sum(cov_ret_1)
  return(weigths)
}