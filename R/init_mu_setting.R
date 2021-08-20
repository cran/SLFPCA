# Getting initial value for coefficients of the mean functions
init_mu_setting <- function(Ly, Lt, bw, xout, basis){
  Lq <- lapply(Ly, function(x){2 * x - 1})
  mufun_est <- fdapace::Lwls1D(bw = bw, kernel_type = "epan", xin = sort(unlist(Lt)),
                               yin = unlist(Lq)[order(unlist(Lt))], xout = xout)
  mufd_est <- fda::smooth.basis(xout, mufun_est, basis)$fd
  return(mufd_est$coefs)
}
