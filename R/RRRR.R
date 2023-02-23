#' Robust Reduced-Rank Regression using Majorisation-Minimisation
#'
#' Majorisation-Minimisation based Estimation for Reduced-Rank Regression with a Cauchy Distribution Assumption.
#' This method is robust in the sense that it assumes a heavy-tailed Cauchy distribution
#' for the innovations. This method is an iterative optimization algorithm. See \code{References} for a similar setting.
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
#' the innovation's distribution, assuming the innovation has a Cauchy distribution.
#'
#' @param y Matrix of dimension N*P. The matrix for the response variables. See \code{Detail}.
#' @param x Matrix of dimension N*Q. The matrix for the explanatory variables to be projected. See \code{Detail}.
#' @param z Matrix of dimension N*R. The matrix for the explanatory variables not to be projected. See \code{Detail}.
#' @param mu Logical. Indicating if a constant term is included.
#' @param r Integer. The rank for the reduced-rank matrix \eqn{AB'}. See \code{Detail}.
#' @param itr Integer. The maximum number of iteration.
#' @param earlystop Scalar. The criteria to stop the algorithm early. The algorithm will stop if the improvement
#'     on objective function is small than \eqn{earlystop * objective_from_last_iteration}.
#' @param initial_mu Vector of length P. The initial value for constant \eqn{\mu} See \code{Detail}.
#' @param initial_A Matrix of dimension P*r. The initial value for matrix \eqn{A}. See \code{Detail}.
#' @param initial_B Matrix of dimension Q*r. The initial value for matrix \eqn{B}. See \code{Detail}.
#' @param initial_D Matrix of dimension P*R. The initial value for matrix \eqn{D}. See \code{Detail}.
#' @param initial_mu Matrix of dimension P*1. The initial value for the constant \eqn{mu}. See \code{Detail}.
#' @param initial_Sigma Matrix of dimension P*P. The initial value for matrix Sigma. See \code{Detail}.
#' @param return_data Logical. Indicating if the data used is return in the output.
#' If set to \code{TRUE}, \code{update.RRRR} can update the model by simply provide new data.
#' Set to \code{FALSE} to save output size.
#'
#' @return A list of the estimated parameters of class \code{RRRR}.
#' \describe{
#' \item{spec}{The input specifications. \eqn{N} is the sample size.}
#' \item{history}{The path of all the parameters during optimization and the path of the objective value.}
#' \item{mu}{The estimated constant vector. Can be \code{NULL}.}
#' \item{A}{The estimated exposure matrix.}
#' \item{B}{The estimated factor matrix.}
#' \item{D}{The estimated coefficient matrix of \code{z}.}
#' \item{Sigma}{The estimated covariance matrix of the innovation distribution.}
#' \item{obj}{The final objective value.}
#' \item{data}{The data used in estimation if \code{return_data} is set to \code{TRUE}. \code{NULL} otherwise.}
#' }
#'
#' @examples
#' set.seed(2222)
#' data <- RRR_sim()
#' res <- RRRR(y=data$y, x=data$x, z = data$z)
#' res
#'
#' @author Yangzhuoran Yang
#' @references Z. Zhao and D. P. Palomar, "Robust maximum likelihood estimation of sparse vector error correction model," in2017 IEEE Global Conference on Signal and Information Processing (GlobalSIP),  pp. 913--917,IEEE, 2017.
#' @importFrom magrittr %>%
#' @export
RRRR <- function(y, x, z = NULL, mu = TRUE, r=1,
                 itr = 100, earlystop = 1e-4,
                 initial_A = matrix(rnorm(P*r), ncol =  r),
                 initial_B = matrix(rnorm(Q*r), ncol =  r),
                 initial_D = matrix(rnorm(P*R), ncol =  R),
                 initial_mu = matrix(rnorm(P)),
                 initial_Sigma = diag(P),
                 return_data = TRUE){


  # warning(str(initial_D))
  if(return_data){
    returned_data <- list(y=y, x=x, z=z)
  } else {
    returned_data <- NULL
  }
  N <- nrow(y)
  P <- ncol(y)
  Q <- ncol(x)
  if(NROW(initial_A) != P)
    stop("Mismatched dimension. The row number of initial_A is not the same as P.")
  if(NROW(initial_mu) != P)
    stop("Mismatched dimension. The row number (length) of initial_mu is not the same as P.")
  if(NROW(initial_B) != Q)
    stop("Mismatched dimension. The row number of initial_B is not the same as Q.")
  if(NCOL(initial_A) != r)
    stop("Mismatched dimension. The column number of initial_A is not the same as r.")
  if(NCOL(initial_B) != r)
    stop("Mismatched dimension. The column number of initial_B is not the same as r.")


  # check if z is NULL
  # add column of ones for mu
  if(!is.null(z)){
    z <- as.matrix(z)
    R <- ncol(z)
    if(NCOL(initial_D) != R)
      stop("Mismatched dimension. The column number of initial_D is not the same as the column number of variable z.")
    if(NROW(initial_D) != P)
      stop("Mismatched dimension. The row number of initial_D is not the same as P.")

    z <- t(z)
    if(mu){
      z <- rbind(z, 1)
      initial_D <- cbind(initial_D, initial_mu)
    }

    znull <- FALSE
  } else {
    R <- 0
    znull <- TRUE
    if (mu){
      z <- matrix(rep(1, N), nrow = 1)
      initial_D <- as.matrix(initial_mu)
    }
  }
  muorz <- mu || !znull

  if(nrow(y) != nrow(x)){
    if(!is.null(z) && nrow(y) != ncol(z))
      stop("The numbers of realizations are not consistant in the inputs.")
  }




  ##
  A <- vector("list", itr + 1)
  B <- vector("list", itr + 1)
  Pi <- vector("list", itr + 1)
  D <- vector("list", itr + 1)
  Sigma <- vector("list", itr + 1)
  A[[1]] <- initial_A
  B[[1]] <- initial_B
  Pi[[1]] <- A[[1]] %*% t(B[[1]])
  D[[1]] <- initial_D
  Sigma[[1]] <- initial_Sigma

  obj <- numeric(itr+1)

  ybar <- vector("list", itr)
  xbar <- vector("list", itr)
  zbar <- vector("list", itr)
  Mbar <- vector("list", itr)

  xk <- vector("list", itr)
  wk <- vector("list", itr)
  # pb <- progress_estimated(itr)
  runtime <- vector("list", itr+1)
  runtime[[1]] <- Sys.time()
  y <- t(y)
  x <- t(x)
  # z <- t(z)
  # z <- rbind(z, 1)
  # warning(str(D))
  # warning(str(z))
  for (k in seq_len(itr)) {
    if(muorz){
      temp <- t(y - Pi[[k]] %*% x - D[[k]] %*% z) %>% split(seq_len(nrow(.)))
    } else {
      temp <- t(y - Pi[[k]] %*% x) %>% split(seq_len(nrow(.)))
    }
    xk[[k]] <- sapply(temp, function(tem) t(tem) %*% solve(Sigma[[k]]) %*% tem)
    wk[[k]] <- 1/(1 + xk[[k]])

    dwk <- diag(sqrt(wk[[k]]))
    ybar[[k]] <- y %*% dwk
    xbar[[k]] <- x %*% dwk
    if(muorz){
      zbar[[k]] <- z %*% dwk

      Mbar[[k]] <- diag(1, N) - t(zbar[[k]]) %*% solve(zbar[[k]] %*% t(zbar[[k]])) %*% zbar[[k]]
      ##
      R_0 <- ybar[[k]] %*% Mbar[[k]]
      R_1 <- xbar[[k]] %*% Mbar[[k]]
    } else {
      R_0 <- ybar[[k]]
      R_1 <- xbar[[k]]
    }

    S_01 <- 1/N * R_0 %*% t(R_1)
    S_10 <- 1/N * R_1 %*% t(R_0)
    S_00 <- 1/N * R_0 %*% t(R_0)
    S_11 <- 1/N * R_1 %*% t(R_1)

    SSSSS <- solve(expm::sqrtm(S_11)) %*% S_10 %*% solve(S_00) %*% S_01 %*% solve(expm::sqrtm(S_11))
    V <- eigen(SSSSS)$vectors[, seq_len(r)]

    B[[k + 1]] <- solve(expm::sqrtm(S_11)) %*% V
    # beta_hat <- beta_hat * (1/beta_hat[[1]])

    # alpha_hat <- S_01 %*% (beta_hat) %*% solve(t(beta_hat) %*% S_11 %*% (beta_hat))
    A[[k + 1]] <- S_01 %*% B[[k + 1]]
    Pi[[k + 1]] <- A[[k + 1]] %*% t(B[[k + 1]])

    if(muorz){
      D[[k + 1]] <- (ybar[[k]] - A[[k + 1]] %*% t(B[[k + 1]]) %*% xbar[[k]]) %*% t(zbar[[k]]) %*% solve(zbar[[k]] %*% t(zbar[[k]]))
      # out$spec$Gamma
      Sigma[[k + 1]] <- (P + 1)/(N - 2) *
        (ybar[[k]] - A[[k + 1]] %*% t(B[[k + 1]]) %*%
           xbar[[k]] - D[[k + 1]] %*% zbar[[k]]) %*%
        t(ybar[[k]] - A[[k + 1]] %*% t(B[[k + 1]]) %*%
            xbar[[k]] - D[[k + 1]] %*% zbar[[k]])
    } else {
      Sigma[[k + 1]] <- (P + 1)/(N - 2) *
        (ybar[[k]] - A[[k + 1]] %*% t(B[[k + 1]]) %*%
           xbar[[k]]) %*%
        t(ybar[[k]] - A[[k + 1]] %*% t(B[[k + 1]]) %*%
            xbar[[k]])

    }
    obj[[k]] <- 1/2 * log(det(Sigma[[k]])) +(1+P)/(2*(N)) * sum(log(1+xk[[k]]))
    # obj[[k]] <- 1/2 * log(det(Sigma[[k+1]])) +(1+P)/(2*(N)) * sum(log(1+xk[[k]]))
    last <- k+1
    # pb$tick()$print()
    runtime[[k+1]] <- Sys.time()
    if(k>1){
      if(abs((obj[[k]]-obj[[k-1]])/obj[[k-1]]) <= 1e-5){
        break()
      }
    }
  }


  if(muorz){
    temp <- t(y - Pi[[k+1]] %*% x - D[[k+1]] %*% z) %>% split(seq_len(nrow(.)))
  } else {
    temp <- t(y - Pi[[k+1]] %*% x) %>% split(seq_len(nrow(.)))
  }
  xkk <- sapply(temp, function(tem) t(tem) %*% solve(Sigma[[k+1]]) %*% tem)
  obj[[k+1]] <- 1/2 * log(det(Sigma[[k+1]])) +(1+P)/(2*(N)) * sum(log(1+xkk))
  # warning(str(D))
  if(mu){
    mu <- lapply(D[sapply(D, function(x) !is.null(x))], function(x) x[,ncol(x)])
    if(!znull)
      D <- lapply(D[sapply(D, function(x) !is.null(x))], function(x) x[,seq_len(ncol(x)-1)])
  } else {
    mu <- NULL
  }
  if(znull){
    D <- NULL
  }
  history <- list(mu = mu, A = A, B = B, D = D, Sigma = Sigma, obj = obj, runtime = c(0,diff(do.call(base::c,runtime)))) %>%
    # lapply(function(x) x[sapply(x, function(z) !is.null(z))])
    lapply(function(x) x[seq_len(last)])
  output <- list(spec = list(N = N, P = P,Q=Q, R=R,  r = r),
                 history = history,
                 mu = mu[[last]],
                 A = A[[last]],
                 B = B[[last]],
                 D = D[[last]],
                 Sigma = Sigma[[last]],
                 obj = obj[[last]],
                 data = returned_data)
  return(new_RRRR(output))
}


