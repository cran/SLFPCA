# Some functions for binary case
linkfun <- function(x){
  1 - 1/(exp(x) + 1)
}

linkinv <- function(p){
  log(p/(1-p))
}

linkder <- function(x){
  t <- exp(x)
  t/(t + 1)^2
}
