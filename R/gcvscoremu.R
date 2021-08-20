# Compute GCV scores for mean function estimates with various roughness parameters
gcvscoremu <- function(Xtilde, B, B_der, penpar){

  sum_m <- nrow(B)
  hat_m <- B %*% solve(t(B) %*% B + sum_m * penpar * B_der) %*% t(B)
  return(sum((Xtilde - hat_m %*% Xtilde)^2)/(1 - mean(diag(hat_m)))^2)

}
