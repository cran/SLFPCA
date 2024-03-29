SLFPCA_sub <- function(Ly, Lt, interval, npc, nknots, norder, kappa_theta, sparse_pen,
                       nRegGrid = 51, bwmu_init = 0.5, bwcov_init = 1, kappa_mu,
                       itermax = 100, tol = 10){

  n <- length(Ly)
  start_time <- interval[1]
  end_time <- interval[2]
  gridequal <- seq(start_time, end_time, length.out = nRegGrid)
  Lq <- lapply(Ly, function(x){2 * x - 1})

  #B-spline
  basis_mod <- fda::create.bspline.basis(interval, nbasis = nknots + norder, norder = norder,
                                         breaks = seq(start_time, end_time, length.out = nknots + 2))
  B_der <- fda::bsplinepen(basis_mod)
  B <- splines::bs(x = unlist(Lt), degree = norder - 1, knots = seq(start_time, end_time, length.out = nknots + 2)[-c(1, nknots + 2)],
                   intercept = T)
  L <- ncol(B)
  BBinv <- solve(t(B) %*% B)

  Blist <- list()
  index <- 1
  for(i in 1:n){
    Blist[[i]] <- B[index:(index + length(Lt[[i]]) - 1),]
    index <- index + length(Lt[[i]])
  }

  ####### initial value settings####################
  # mean function
  mu_ini <- init_mu_setting(Ly, Lt, bw = bwmu_init, xout = gridequal, basis = basis_mod)

  # eigenfunctions and scores
  init_eig <- init_eig_setting(Ly, Lt, bwmu = bwmu_init, bwcov = bwcov_init, xout = gridequal,
                               npc = npc, basis = basis_mod)
  Theta_ini <- init_eig$Theta_est
  score_ini <- init_eig$score_est


  #################Estimation###################
  BICscore <- array(0, c(length(kappa_mu), length(sparse_pen), length(kappa_theta)))
  mu_res <- array(0, c(L, length(kappa_mu), length(sparse_pen), length(kappa_theta)))
  score_res <- array(0, c(n, npc, length(kappa_mu), length(sparse_pen), length(kappa_theta)))
  Theta_res <- array(0, c(npc, L, length(kappa_mu), length(sparse_pen), length(kappa_theta)))

  for(mm in 1:length(kappa_mu)){
    for(ii in 1:length(sparse_pen)){
      for(jj in 1:length(kappa_theta)){

        print(paste("Tuning parameters: kappa_theta = ", kappa_theta[jj],
                    ",  sparse_pen = ", sparse_pen[ii], ", kappa_mu = ",
                    kappa_mu[mm], sep = ""))

        mu_est <- mu_ini
        Theta_est <- Theta_ini
        score_est <- score_ini

        loglike_loop <- NULL
        loop <- 1
        while(loop <= itermax){ # outer loop

          ################estimation of mean function#######################
          # Compute Xtilde
          Btheta <- NULL
          Bmu <- NULL
          for(i in 1:n){
            Bmu <- c(Bmu, Blist[[i]] %*% mu_est)
            Btheta <- c(Btheta, Blist[[i]] %*% t(Theta_est) %*% score_est[i,])
          }
          Delta <- Bmu + Btheta

          loglike <- -sum(log(linkfun(unlist(Lq) * Delta)))

          XX <- Delta + 4 * unlist(Lq) * (1 - linkfun(unlist(Lq) * Delta))

          Xtilde <- XX - Btheta

          # estimate
          sum_m <- nrow(B)
          mu_est <- solve(t(B) %*% B + sum_m * kappa_mu[mm] * B_der) %*% t(B) %*% Xtilde

          #################estimation of scores and eigenfunctions###################
          df <- NULL

          for(k in 1:npc){

            Bthetasub <- NULL
            for(i in 1:n){
              Bthetasub <- c(Bthetasub, Blist[[i]] %*% t(Theta_est)[,-k] %*% score_est[i, -k])
            }
            Xk <- XX - B %*% mu_est - Bthetasub

            m <- 1
            while(m <= itermax){

              ## score estimation
              phi_est <- B %*% Theta_est[k,]

              index <- 1
              for(i in 1:n){
                range <- index:(index + length(Lt[[i]]) - 1)
                score_est[i,k] <- sum(phi_est[range] * Xk[range])/sum(phi_est[range]^2)
                index <- index + length(Lt[[i]])
              }
              nanid <- which(is.nan(score_est[,k]))
              score_est[nanid, k] <- mean(score_est[,k], na.rm = T)

              ## eigenfunction estimation
              U <- matrix(0, nrow = 1, ncol = L)
              Ulist <- list()
              for(i in 1:n){
                Ulist[[i]] <- score_est[i, k] * Blist[[i]]
                U <- rbind(U, Ulist[[i]])
              }
              U <- U[-1,]
              Theta_est[k,] <- tryCatch(slos_temp(Xk, U, Maxiter = 100, lambda = sparse_pen[ii],
                                                  gamma = kappa_theta[jj], beta.basis = basis_mod,
                                                  absTol = 0.0004, Cutoff = 0),
                                        error = function(err){return(rep(Inf, nknots + norder))})
              if(sum(Theta_est[k,]) == Inf){
                Theta_est[k,] <- rep(0.5, nknots + norder)
                break
              }

              ## negative loglikelihood
              loglike_old <- loglike

              Btheta <- NULL
              Bmu <- NULL
              for(i in 1:n){
                Bmu <- c(Bmu, Blist[[i]] %*% mu_est)
                Btheta <- c(Btheta, Blist[[i]] %*% t(Theta_est) %*% score_est[i,])
              }
              Delta <- Bmu + Btheta
              loglike <- -sum(log(linkfun(unlist(Lq) * Delta)))

              # cat(ii, " ", jj, "  ", k, "  ", m, "  ", loglike, "\n")

              if(abs(loglike - loglike_old) < tol){
                break
              }

              m <- m + 1

            }

            ####df
            if(sum(abs(diff(c(Theta_est[k,])))) == 0){
              df[k] <- Inf
              break
            }else{
              sparse.idx = which(Theta_est[k,] == 0)
              if(length(sparse.idx) == 0)
              {
                ula = U
                vla = B_der
              }
              else{
                ula  = U[, -sparse.idx]
                vla  = B_der[-sparse.idx, -sparse.idx]
              }
              hat2 = ula%*%solve(t(ula)%*%ula + sum_m * kappa_theta[jj]*vla)%*%t(ula)
              df[k]  = psych::tr(hat2)
            }

          }

          if(sum(df) == Inf){
            break
          }

          loglike_loop[loop] <- loglike
          if(loop != 1){

            if((abs(loglike_loop[loop] - loglike_loop[loop - 1]) < tol)|(loglike_loop[loop] > loglike_loop[loop - 1])){
              break
            }

          }

          print(loop)
          print(loglike_loop[loop])
          loop <- loop + 1

        }

        mu_res[, mm, ii, jj] <- mu_est
        score_res[, , mm, ii, jj] <- score_est
        Theta_res[, , mm, ii, jj] <- Theta_est

        #BIC
        BICscore[mm, ii, jj] <- 2 * loglike + log(length(unlist(Lt))) * sum(df) +
          0.5 * sum(df) * log(npc * (nknots + norder) + n * npc)


      }

    }

  }

  ######choose tuning parameters that achieve lowest BIC score#############
  mm <- which(BICscore == min(BICscore), arr.ind = T)[1]
  ii <- which(BICscore == min(BICscore), arr.ind = T)[2]
  jj <- which(BICscore == min(BICscore), arr.ind = T)[3]
  mu_est <- mu_res[, mm, ii, jj]
  Theta_est <- as.matrix(Theta_res[, , mm, ii, jj])
  score_est <- as.matrix(score_res[, , mm, ii, jj])

  B_reg <- splines::bs(x = gridequal, degree = norder - 1,
                       knots = seq(start_time, end_time, length.out = nknots + 2)[-c(1, nknots + 2)],
                       intercept = T)
  mufd_est <- fda::smooth.basis(gridequal, B_reg %*% mu_est, basis_mod)$fd
  if(npc == 1){
    eigfd_est <- fda::smooth.basis(gridequal, B_reg %*% Theta_est, basis_mod)$fd
  }else{
    eigfd_est <- fda::smooth.basis(gridequal, B_reg %*% t(Theta_est), basis_mod)$fd
  }
  multi <- NULL
  eigfd_est_st <- list()
  for(i in 1:npc){
    multi[i] <- sqrt(as.numeric(fda::inprod(eigfd_est[i], eigfd_est[i])))
    eigfd_est_st[[i]] <- eigfd_est[i] * (1/multi[i])
    score_est[,i] <- score_est[,i] * multi[i]
  }

  ret <- list()

  ret$mufd <- mufd_est
  ret$eigfd_list <- eigfd_est_st
  ret$score <- score_est
  ret$kappa_mu <- kappa_mu[mm]
  ret$kappa_theta <- kappa_theta[jj]
  ret$sparse_pen <- sparse_pen[ii]
  ret$EBICscore <- BICscore[mm, ii, jj]

  return(ret)

}
