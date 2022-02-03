# Getting initial value for coefficients of the eigenfunctions and scores
init_eig_setting <- function(Ly, Lt, bwmu, bwcov, xout, npc, basis){

  n <- length(Lt)
  Lq <- lapply(Ly, function(x){2 * x - 1})

  gridobs <- sort(unique(unlist(Lt)))
  mufine_est <- fdapace::Lwls1D(bw = bwmu, kernel_type = "epan", xin = sort(unlist(Lt)),
                                yin = unlist(Lq)[order(unlist(Lt))], xout = gridobs)

  Lq_cen <- list()
  for(i in 1:n){
    Lq_cen[[i]] <- Lq[[i]] - mufine_est[match(Lt[[i]], gridobs)]
  }
  xin2D <- NULL
  yin2D <- NULL
  for(i in 1:n){
    xin2D <- rbind(xin2D, t(utils::combn(Lt[[i]], 2)))
    yin2D <- rbind(yin2D, t(utils::combn(Lq_cen[[i]], 2)))
  }
  xin_pair <- rbind(xin2D, cbind(xin2D[,2], xin2D[,1]))
  yin_pair <- yin2D[,1] * yin2D[, 2]
  yin_pair <- c(yin_pair, yin_pair)
  covfun_est <- fdapace::Lwls2D(bw = bwcov, kern = "epan", xin = xin_pair, yin = yin_pair,
                                xout1 = xout, xout2 = xout)

  eigfun_est <- eigen(covfun_est)$vectors[,1:npc]
  eigfd_est <- fda::smooth.basis(xout, eigfun_est, basis)$fd
  Theta_est <- t(eigfd_est$coefs)

  eigfd_adj <- list()
  eigval <- eigen(covfun_est)$values[1:npc]
  for(i in 1:npc){
    multi_eig <- as.numeric(fda::inprod(eigfd_est[i], eigfd_est[i]))
    eigfd_adj[[i]] <- eigfd_est[i] * multi_eig^(-1/2)
    Theta_est[i,] <- Theta_est[i,] * multi_eig^(-1/2)
    eigval[i] <- eigval[i] * multi_eig
  }

  score_est <- matrix(0, nrow = n, ncol = npc)

  init <- list()
  init$Theta_est <- Theta_est
  init$score_est <- score_est

  return(init)
}
