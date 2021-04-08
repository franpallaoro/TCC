Sym2Vech = function(k, A){
  vechA = matrix(0, nrow = k*(k+1)/2, ncol = 1)
  l = 1
  for (j in 1:k) {
   for (i in j:k) {
     vechA[l] = A[i,j]
     l = l+1
   } 
  }
  return(vechA)
}



# function [ vechA ] = Sym2Vech(k, A )
# # this function puts A into vech of A
# 
# vechA= zeros(k*(k+1)/2, 1);
# l=1;
# for j=1:k
# for i=j:k
# vechA(l) = A(i,j);
# l=l+1;
# end
# end
# 
# end