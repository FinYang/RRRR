
#' Simulating data for Reduced-Rank Regression
#'
#' Simulate data for Reduced-rank regression. See \code{Detail} for the formulation
#' of the simulated data.
#'
#' The data simulated can be used for the standard reduced-rank regression testing
#' with the following formulation
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
#' The function simulates \eqn{x}, \eqn{z} from multivariate normal distribution and \eqn{y} by specifying
#' parameters \eqn{\mu}, \eqn{A}, \eqn{B}, \eqn{D}, and \eqn{varcov}, the covariance matrix of
#' the innovation's distribution.  The constant \eqn{\mu} and the term \eqn{Dz} can be
#' dropped by setting \code{NULL} for arguments \code{mu} and \code{D}. The \code{innov} in the argument is
#' the collection of innovations of all the realizations.
#'
#' @param N Integer. The total number of simulated realizations.
#' @param P Integer. The dimension of the response variable matrix. See \code{Detail}.
#' @param Q Integer. The dimension of the explanatory variable matrix to be projected. See \code{Detail}.
#' @param R Integer. The dimension of the explanatory variable matrix not to be projected. See \code{Detail}.
#' @param r Integer. The rank of the reduced rank coefficient matrix. See \code{Detail}.
#' @param mu Vector with length P. The constants. Can be \code{NULL} to drop the term. See \code{Detail}.
#' @param A Matrix with dimension P*r. The exposure matrix. See \code{Detail}.
#' @param B Matrix with dimension Q*r. The factor matrix. See \code{Detail}.
#' @param D Matrix with dimension P*R. The coefficient matrix for \code{z}. Can be \code{NULL} to drop the term. See \code{Detail}.
#' @param varcov Matrix with dimension P*P. The covariance matrix of the innovation. See \code{Detail}.
#' @param innov Matrix with dimension N*P. The innovations. Default to be simulated from a Student t distribution, See \code{Detail}.
#' @param mean_x Integer. The mean of the normal distribution \eqn{x} is simulated from.
#' @param mean_z Integer. The mean of the normal distribution \eqn{z} is simulated from.
#' @param x Matrix with dimension N*Q. Can be used to specify \eqn{x} instead of simulating form a normal distribution.
#' @param z Matrix with dimension N*R. Can be used to specify \eqn{z} instead of simulating form a normal distribution.
#'
#' @return A list of the input specifications and the data \eqn{y}, \eqn{x}, and \eqn{z}, of class \code{RRR_data}.
#' \describe{
#' \item{y}{Matrix of dimension N*P}
#' \item{x}{Matrix of dimension N*Q}
#' \item{z}{Matrix of dimension N*R}
#' }
#'
#' @examples
#' set.seed(2222)
#' data <- RRR_sim()
#'
#' @author Yangzhuoran Yang
#' @importFrom stats rnorm
#' @export
RRR_sim <- function(N = 1000, P = 3, Q = 3, R = 1, r = 1,
                    mu = rep(0.1, P),
                    A = matrix(rnorm(P*r), ncol = r),
                    B = matrix(rnorm(Q*r), ncol = r),
                    D = matrix(rnorm(P*R), ncol = R),
                    varcov = diag(P),
                    innov = mvtnorm::rmvt(N, sigma = varcov, df = 3),
                    mean_x = 0, mean_z = 0,
                    x = NULL,
                    z = NULL
) {
  if(P != NROW(A)) stop("P is not the same with the row number of matrix A.")
  if(Q != NROW(B)) stop("Q is not the same with the row number of matrix B.")
  if(R==0 || is.null(D)){
    R <- 0
  } else {
    if(R != NCOL(D)) stop("R is not the same with the column number of matrix D.")
  }

  if(is.null(x))  x <- mvtnorm::rmvnorm(N, mean = rep(mean_x, Q))

  if((!is.null(D)) && is.null(z)){
    if (R > 1) {
      z <- mvtnorm::rmvnorm(N, mean = rep(mean_z, R))
    } else {
      z <- stats::rnorm(N, mean = mean_z)
    }
  }
  C <- A %*% t(B)
  # r <- matrixcalc::matrix.rank(Pi)
  if(is.null(mu)){
    constant <- 0
  } else {
    constant <- do.call(rbind, rep(list(mu), N))
  }
  if(is.null(D)){
    y <- constant + x %*% t(C) + innov
    # output <- list(spec = list(
    #   P = P, N = N, Q = Q, R = R, A = A, B = B, mu = mu, r = r,
    #   Sigma = varcov
    # ), y = y, x = x)
  } else {
    y <- constant + x %*% t(C) + z %*% t(D) + innov
  }
  output <- list(spec = list(
    P = P, N = N, Q = Q, R = R, A = A, B = B, mu = mu, D = D,  r = r,
    Sigma = varcov, innov = innov
  ), y = y, x = x, z = z)
  class(output) <- c("RRR_data", "list")
  return(output)
}
