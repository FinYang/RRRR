#' Reduced-Rank Regression using Gaussian MLE
#'
#' Gaussian Maximum Likelihood Estimation method for Reduced-Rank Regression.
#' This method is not robust in the sense that it assumes a Gaussian distribution
#' for the innovations
#' which does not take into account the heavy-tailedness of the true distribution and
#' outliers.
#'
#' The formulation of the reduced-rank regression is as follow:
#' \deqn{y = \mu +AB'  x + D z+innov,}
#' where for each realization \eqn{y} is a vector of dimension \eqn{P} for the \eqn{P} response variables,
#' \eqn{x} is a vector of dimension \eqn{Q} for the \eqn{Q} explanatory variables that will be projected to
#' reduce the rank,
#' \eqn{z} is a vector of dimension \eqn{R} for the \eqn{R} explanatory variables
#' that will not be projected,
#' \eqn{\mu} is the constant vector of dimension \eqn{P},
#' \eqn{innov} is the innovation vector of dimension \eqn{P},
#' \eqn{D} is a coefficient matrix for \eqn{z} with dimension \eqn{P*R},
#' \eqn{A} is the so called exposure matrix with dimension \eqn{P*r}, and
#' \eqn{B} is the so called factor matrix with dimension \eqn{Q*r}.
#' The matrix resulted from \eqn{AB'} will be a reduced rank coefficient matrix with rank of \eqn{r}.
#' The function estimates parameters \eqn{\mu}, \eqn{A}, \eqn{B}, \eqn{D}, and \eqn{Sigma}, the covariance matrix of
#' the innovation's distribution, assuming the innovation has a Gaussian distribution.
#'
#' @param y Matrix of dimension N*P. The matrix for the response variables. See \code{Detail}.
#' @param x Matrix of dimension N*Q. The matrix for the explanatory variables to be projected. See \code{Detail}.
#' @param z Matrix of dimension N*R. The matrix for the explanatory variables not to be projected. See \code{Detail}.
#' @param mu Logical. Indicating if a constant term is included.
#' @param r Integer. The rank for the reduced-rank matrix \eqn{AB'}. See \code{Detail}.
#'
#' @return A list of the estimated parameters of class \code{RRR}.
#' \describe{
#' \item{spec}{The input specifications. \eqn{N} is the sample size.}
#' \item{mu}{The estimated constant vector. Can be \code{NULL}.}
#' \item{A}{The estimated exposure matrix.}
#' \item{B}{The estimated factor matrix.}
#' \item{D}{The estimated coefficient matrix of \code{z}. Can be \code{NULL}.}
#' \item{Sigma}{The estimated covariance matrix of the innovation distribution.}
#' }
#'
#' @examples
#' set.seed(2222)
#' data <- RRR_sim()
#' res <- RRR(y=data$y, x=data$x, z = data$z)
#' res
#' @seealso For robust reduced-rank regression estimation see function \code{\link{RRRR}}.
#' @author Yangzhuoran Yang
#' @references S. Johansen, "Estimation and Hypothesis Testing of Cointegration Vectors in Gaussian Vector Autoregressive Models,"Econometrica, vol. 59,p. 1551, Nov. 1991.
#' @export
RRR <- function(y, x, z = NULL, mu = TRUE, r=1){

  N <- nrow(y)
  if(!is.null(z)){
    z <- as.matrix(z)
    R <- ncol(z)
    z <- t(z)
    if(mu)  z <- rbind(z, 1)
    znull <- FALSE
  } else {
    R <- 0
    znull <- TRUE
    if (mu){
      z <- matrix(rep(1, N), nrow = 1)
    }
  }

  if(nrow(y) != nrow(x)){
    if(!is.null(z) && nrow(y) != ncol(z))
      stop("The numbers of realizations are not consistant in the inputs.")
  }


  P <- ncol(y)
  Q <- ncol(x)
  y <- t(y)
  x <- t(x)
  if(mu || !znull){
    M <- diag(1, N) - t(z) %*% solve(z %*% t(z)) %*% z
    R_0 <- y %*% M
    R_1 <- x %*% M
  } else {
    R_0 <- y
    R_1 <- x
  }
  S_01 <- 1/(N) * R_0 %*% t(R_1)
  S_10 <- 1/(N) * R_1 %*% t(R_0)
  S_00 <- 1/(N) * R_0 %*% t(R_0)
  S_11 <- 1/(N) * R_1 %*% t(R_1)

  SSSSS <- solve(expm::sqrtm(S_11)) %*% S_10 %*% solve(S_00) %*% S_01 %*% solve(expm::sqrtm(S_11))
  V <- eigen(SSSSS)$vectors[, seq_len(r)]

  B_hat <- solve(expm::sqrtm(S_11)) %*% V
  # beta_hat <- beta_hat * (1/beta_hat[[1]])

  # alpha_hat <- S_01 %*% (beta_hat) %*% solve(t(beta_hat) %*% S_11 %*% (beta_hat))
  A_hat <- S_01 %*% B_hat
  Pi_hat <- A_hat %*% t(B_hat)

  if(mu || !znull){
    Gamma_hat <- (y - A_hat %*% t(B_hat) %*% x) %*% t(z) %*% solve(z %*% t(z))
    Sigma_hat <- 1/(N) * (y - A_hat %*% t(B_hat) %*% x - Gamma_hat %*% z) %*% t(y - A_hat %*% t(B_hat) %*% x - Gamma_hat %*% z)
  } else {
    # Gamma_hat <- (y - A_hat %*% t(B_hat) %*% x)
    Gamma_hat <- NULL
    Sigma_hat <- 1/(N) * (y - A_hat %*% t(B_hat) %*% x) %*% t(y - A_hat %*% t(B_hat) %*% x )
  } #else {
  #   Gamma_hat <- (y - A_hat %*% t(B_hat) %*% x)
  #   Sigma_hat <- 1/(N) * (y -Gamma_hat - A_hat %*% t(B_hat) %*% x) %*% t(y - Gamma_hat- A_hat %*% t(B_hat) %*% x )
  # }

  # out$spec$Gamma
  D <- NULL
  if(!znull){
    if(mu){
      D <- Gamma_hat[,-ncol(Gamma_hat), drop = FALSE]
    }  else {
      D <- Gamma_hat
    }
  }
  if(mu)
    mu <- Gamma_hat[,ncol(Gamma_hat), drop = FALSE]
  else
    mu <- NULL
  output <- list(spec = list(N = N, P = P, Q = Q,R=R, r = r),
                 mu = mu , A = A_hat, B = B_hat,
                 D = D, Sigma = Sigma_hat)
  return(new_RRR(output))

}

