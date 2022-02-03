#' @title Sparse logistic functional principal component analysis
#'
#' @description Sparse logistic functional principal component analysis (SLFPCA) for binary data. The estimated eigenfunctions from SLFPCA can be strictly zero on some sub-intervals, which is helpful for interpretation.
#'
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the binary measurements of each subject at the observation time correspond to \code{Lt}.
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param interval A \code{vector} of length two denoting the supporting interval.
#' @param npc An integer denoting the number of FPCs.
#' @param L_list A \code{vector} denoting the candidates for the number of B-spline basis functions.
#' @param norder An integer denoting the order of the using B-spline basis, which is one higher than their degree.
#' @param kappa_theta A \code{vector} denoting the smoothing parameters for eigenfunctions, the optimal tuning parameter is chosen from them.
#' @param sparse_pen A \code{vector} denoting the sparseness parameters for eigenfunctions, the optimal tuning parameter is chosen from them.
#' @param nRegGrid An integer denoting the number of equally spaced time points in the supporting interval. The eigenfunctions and mean function are estimated at these equally spaced time points first, before transforming into functional data object. (default: 51)
#' @param bwmu_init A scalar denoting the bandwidth for mean function estimation in the setting of initial values. (default: 0.5)
#' @param bwcov_init A scalar denoting the bandwidth for covariance function estimation in the setting of initial values. (default: 1)
#' @param stepmu A scalar denoting the length between each considered smoothing parameter for mean function. For selection of smoothing parameter for mean function, we start from zero and increase the value until GCV score begins increasing.
#' @param mucand_num An integer denoting the maximum number of the considered smoothing parameter for mean function. (default: 100)
#' @param itermax An integer denoting the maximum number of iterations. (default: 100)
#' @param tol A scalar. When difference of the loglikelihood functions between two consecutive iteration is less than \code{tol}, the convergence is supposed to be reached. (default: 0.5)
#'
#' @return A \code{list} containing the following components:
#' \item{mufd}{A functional data object for the mean function estimate.}
#' \item{eigfd_list}{A \code{list} containing \code{npc} functional data objects, which are the eigenfunction estimates.}
#' \item{score}{A \emph{n} by \code{npc} \code{matrix} containing the estimates of the FPC scores, where \emph{n} is the sample size.}
#' \item{kappa_mu}{A scalar denoting the selected smoothing parameter for mean function.}
#' \item{kappa_theta}{A scalar denoting the selected smoothing parameter for eigenfunctions.}
#' \item{sparse_pen}{A scalar denoting the selected sparseness parameter for eigenfunctions.}
#' \item{L_select}{A scalar denoting the selected number of B-spline basis functions.}
#' \item{EBICscore}{A \code{vector} denoting the selected EBIC scores corresponding to different numbers of B-spline basis functions in \code{L_list}.}
#' @export
#'
#' @references
#' \cite{Rou Zhong, Shishi Liu, Haocheng Li, Jingxiao Zhang (2021). "Sparse logistic functional principal component analysis for binary data." <arXiv: https://arxiv.org/abs/2109.08009>.}
#'
#' @examples
#' #Generate data
#' n <- 100
#' npc <- 1
#' interval <- c(0, 10)
#' gridequal <- seq(0, 10, length.out = 51)
#' basis <- fda::create.bspline.basis(c(0, 10), nbasis = 13, norder = 4,
#'          breaks = seq(0, 10, length.out = 11))
#' meanfun <- function(t){2 * sin(pi * t/5)/sqrt(5)}
#' lambda_1 <- 3^2 #the first eigenvalue
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)))
#' eigfun <- list()
#' eigfun[[1]] <- function(t){cos(pi * t/5)/sqrt(5)}
#' eigfd <- list()
#' for(i in 1:npc){
#'   eigfd[[i]] <- fda::smooth.basis(gridequal, eigfun[[i]](gridequal), basis)$fd
#' }
#' DataNew <- GenBinaryFD(n, interval, sparse = 8:12, regular = FALSE,
#'            meanfun = meanfun, score, eigfd)
#' SLFPCA_list <- SLFPCA(DataNew$Ly, DataNew$Lt, interval, npc, L_list = 13,
#'                norder = 4, kappa_theta = 0.2, sparse_pen = 0,
#'                nRegGrid = 51, stepmu = 0.005)
#' plot(SLFPCA_list$eigfd_list[[1]])
#'
SLFPCA <- function(Ly, Lt, interval, npc, L_list, norder, kappa_theta, sparse_pen,
                   nRegGrid = 51, bwmu_init = 0.5, bwcov_init = 1, stepmu, mucand_num = 100,
                   itermax = 100, tol = 0.5){

  EBICscore <- NULL
  SLFPCA_l <- list()
  for(l_id in 1:length(L_list)){

    nknots <- L_list[l_id] - norder
    SLFPCA_l[[l_id]] <- SLFPCA_sub(Ly, Lt, interval, npc, nknots = nknots, norder,
                                   kappa_theta, sparse_pen, nRegGrid, bwmu_init, bwcov_init,
                                   stepmu, mucand_num, itermax, tol)
    EBICscore[l_id] <- SLFPCA_l[[l_id]]$EBICscore

    print(paste("L =", L_list[l_id]))

  }

  l_id_select <- which.min(EBICscore)
  L_select <- L_list[l_id_select]
  SLFPCA_ret <- SLFPCA_l[[l_id_select]]

  ret <- list()

  ret$mufd <- SLFPCA_ret$mufd
  ret$eigfd_list <- SLFPCA_ret$eigfd_list
  ret$score <- SLFPCA_ret$score
  ret$kappa_mu <- SLFPCA_ret$kappa_mu
  ret$kappa_theta <- SLFPCA_ret$kappa_theta
  ret$sparse_pen <- SLFPCA_ret$sparse_pen
  ret$L_select <- L_select
  ret$EBICscore <- EBICscore

  return(ret)

}
