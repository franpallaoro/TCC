Compute_S_l_matrix <- function (asset_group_vec){

  nr_groups = max(asset_group_vec)
  n   = length(asset_group_vec)

  S_l = matrix(0,n,nr_groups)

  for(i in 1:n){
    S_l[i,asset_group_vec[i]] = 1
  }

  S_l_c = matrix(1,n,nr_groups) - S_l

  return(list("S_l" = S_l, "S_l_c" = S_l_c))
}


# compute selection matrix N x g
# put a 1 when asset i (i = 1,..N) belongs to group j (j =
                                                         # 1,...g)