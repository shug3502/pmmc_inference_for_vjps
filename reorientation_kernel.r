reorientation_kernel <- function(theta2){
#assume a circular uniform distribution for reorientations
#takes theta2 (previous reorientation angle as argument, but no dependence on this in this case)
#################### 
  theta1 = pi*(2*runif(1)-1)
  return(theta1)
}
