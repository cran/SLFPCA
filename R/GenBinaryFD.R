#' @title Generate binary functional data
#'
#' @description Generate binary functional data through latent process.
#'
#' @param n An integer denoting the number of sample size.
#' @param interval A \code{vector} of length two denoting the supporting interval.
#' @param sparse A \code{vector} denoting the possible numbers of observation size. The elements are chosen with equal chance. The length of \code{sparse} must be one if \code{regular = TRUE}.
#' @param regular Logical; If \code{TRUE}, the observation grids are equally-spaced.
#' @param meanfun A function for the mean.
#' @param score A \emph{n} by \code{npc} \code{matrix} containing the FPC scores, where \code{npc} is the number of FPCs.
#' @param eigfd A \code{list} containing functional objects for the eigenfunctions.
#'
#' @return A \code{list} containing the following components:
#' \item{Lt}{A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.}
#' \item{Lx}{A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains vales of the latent process of each subject at the observation time correspond to \code{Lt}.}
#' \item{Ly}{A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the binary measurements of each subject at the observation time correspond to \code{Lt}.}
#' @export
#'
#' @examples
#' n <- 100
#' npc <- 2
#' interval <- c(0, 10)
#' gridequal <- seq(0, 10, length.out = 51)
#' basis <- fda::create.bspline.basis(c(0, 10), nbasis = 13, norder = 4,
#'          breaks = seq(0, 10, length.out = 11))
#' meanfun <- function(t){2 * sin(pi * t/5)/sqrt(5)}
#' lambda_1 <- 3^2 #the first eigenvalue
#' lambda_2 <- 2^2 #the second eigenvalue
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)), rnorm(n, 0, sqrt(lambda_2)))
#' eigfun <- list()
#' eigfun[[1]] <- function(t){cos(pi * t/5)/sqrt(5)}
#' eigfun[[2]] <- function(t){sin(pi * t/5)/sqrt(5)}
#' eigfd <- list()
#' for(i in 1:npc){
#'   eigfd[[i]] <- fda::smooth.basis(gridequal, eigfun[[i]](gridequal), basis)$fd
#' }
#' DataNew <- GenBinaryFD(n, interval, sparse = 8:12, regular = FALSE,
#'            meanfun = meanfun, score, eigfd)
#'
GenBinaryFD <- function(n, interval, sparse, regular, meanfun, score, eigfd){

  npc <- ncol(score)
  gridequal <- seq(interval[1], interval[2], length.out = sparse[1])

  Lt <- list()
  Lx <- list()
  Ly <- list()
  Lx_reg <- list()
  Lq <- list()
  for(i in 1:n){
    if(regular){
      Lt[[i]] <- gridequal
    }else{
      num <- sample(sparse, 1) #observation size for the i-th subject
      Lt[[i]] <- sort(unique(stats::runif(num, min = interval[1], max = interval[2])))
    }

    scorei <- score[i,]
    x <- meanfun(Lt[[i]])
    for(j in 1:npc){
      x <- x + scorei[j] * fda::eval.fd(Lt[[i]], eigfd[[j]])
    }
    Lx[[i]] <- x
    Ly[[i]] <- sapply(linkfun(Lx[[i]]), function(x){stats::rbinom(1, 1, x)})
  }

  ret <- list()

  ret$Lt <- Lt
  ret$Lx <- Lx
  ret$Ly <- Ly

  return(ret)

}
