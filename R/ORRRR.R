#' Online Robust Reduced-Rank Regression
#'
#'
#'
#' Majorization-Minimization based Estimation for Reduced-Rank Regression with a Cauchy Distribution Assumption
#' This method is robust in the sense that it assumes a heavy-tailed Cauchy distribution
#' for the innovations. This method is an iterative optimization algorithm. See \code{source} for a similar setting.
#'
#' The fomulation of the reduced-rank regression is as follow:
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
#' @param itr Interger. The maximum number of iteration.
#' @param earlystop Scalar. The criteria to stop the algorithm early. The algorithm will stop if the improvement
#'     on objective function is small than \eqn{earlystop * objective_from_last_iteration}.
#' @param initial_mu Vector of length P. The initial value for constant \eqn{\mu} See \code{Detail}.
#' @param initial_A Matrix of dimension P*r. The initial value for matrix \eqn{A}. See \code{Detail}.
#' @param initial_B Matrix of dimension Q*r. The initial value for matrix \eqn{B}. See \code{Detail}.
#' @param initial_D Matrix of dimension P*R. The initial value for matrix \eqn{D}. See \code{Detail}.
#' @param initial_Sigma Matrix of dimension P*P. The initial value for matrix Sigma. See \code{Detail}.
#'
#' @return A list of the estimated parameters of class \code{RRRR}.
#' \describe{
#' \item{spec}{The input specifications. \eqn{N} is the sample size.}
#' \item{history}{The path of all the parameters during optimization and the path of the objective value.}
#' \item{mu}{The estimated constant vector. Can be \code{NULL}.}
#' \item{A}{The estimated exposure matrix.}
#' \item{B}{The estimated factor matrix.}
#' \item{D}{The estimated coefficient matrix of \code{z}.}
#' \item{Sigma}{The estimated covariance matrix of the innovarion distribution.}
#' \item{obj}{The final objective value.}
#' }
#'
#' @examples
#' set.seed(2222)
#' data <- RRR_sim()
#' res <- RRRR(y=data$y, x=data$x, z = data$z)
#' res
#'
#' @author Yangzhuoran Fin Yang
#' @source Z. Zhao and D. P. Palomar, "Robust maximum likelihood estimation of sparse vector error correction model," in2017 IEEE Global Conferenceon Signal and Information Processing (GlobalSIP),  pp. 913--917,IEEE, 2017.
#' @importFrom magrittr %>%
#' @export
ORRRR <- function(y, x, z = NULL, mu = TRUE, r = 1,
                  initial_size = 100, addon = 10,
                  method = c("SMM", "SAA"),
                  submethod = c("optim", "MM", "GMLE"),
                  ...,
                  initial_A = matrix(rnorm(P*r), ncol =  r),
                  initial_B = matrix(rnorm(Q*r), ncol =  r),
                  initial_D = matrix(rnorm(P*R), ncol =  R),
                  initial_Sigma = diag(P),
                  ProgressBar = requireNamespace("dplyr")){
  if (ProgressBar && !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package \"dplyr\" needed for progress bar to work. Please install it.",
         call. = FALSE)
  }
  method <- method[[1]]
  N <- nrow(y)

  P <- ncol(y)
  Q <- ncol(x)
  # check if z is NULL
  # add column of ones for mu
  if(!is.null(z)){
    z <- as.matrix(z)
    R <- ncol(z)
    if(mu){
      z <- cbind(z, 1)
      initial_D <- cbind(initial_D, rnorm(P))
    }
    znull <- FALSE
  } else {
    R <- 0
    znull <- TRUE
    if (mu){
      z <- matrix(rep(1, N))
      initial_D <- matrix(rnorm(P))
    }
  }
  muorz <- mu || !znull

  if(nrow(y) != nrow(x)){
    if(!is.null(z) && nrow(y) != ncol(z))
      stop("The numbers of realizations are not consistant in the inputs.")
  }



  yy <- y
  xx <- x
  zz <- z

  A <- list()
  B <- list()
  Pi <- list()
  D <- list()
  Sigma <- list()
  A[[1]] <- initial_A
  B[[1]] <- initial_B
  Pi[[1]] <- A[[1]] %*% t(B[[1]])
  if(muorz)
    D[[1]] <- initial_D
  Sigma[[1]] <- initial_Sigma

  itr <- (N-initial_size)/addon +1
  if(ProgressBar)
    pb <- dplyr::progress_estimated(itr)
  obj <- numeric(itr+1)

  ybar <- list()
  xbar <- list()
  zbar <- list()
  Mbar <- list()

  xk <- list()
  wk <- list()

  runtime <- vector("list",itr+1)
  runtime[[1]] <- Sys.time()
  for (k in seq_len(itr)) {
    if(k==1){
      y <- t(yy[seq(1, k*initial_size), ])
      x <- t(xx[seq(1, k*initial_size), ])
      if(muorz)
        z <- t(zz[seq(1, k*initial_size), ])
    } else {
      if(floor(itr) < itr && k==floor(itr)){
        y <- t(yy)
        x <- t(xx)
        if(muorz)
          z <- t(zz)
      } else {
        y <- t(yy[seq(1, initial_size+(k-1)*addon),])
        x <- t(xx[seq(1, initial_size+(k-1)*addon),])
        if(muorz)
          z <- t(zz[seq(1, initial_size+(k-1)*addon),])
      }
    }
    N <- ncol(y)
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
    # obj
    obj[[k]] <- 1/2 * log(det(Sigma[[k]])) +(1+P)/(2*(N)) * sum(log(1+xk[[k]]))
    if(ProgressBar)
      pb$tick()$print()
    runtime[[k+1]] <- Sys.time()
  }

  if(muorz){
    temp <- t(y - Pi[[k+1]] %*% x - D[[k+1]] %*% z) %>% split(seq_len(nrow(.)))
  } else {
    temp <- t(y - Pi[[k+1]] %*% x) %>% split(seq_len(nrow(.)))
  }
  xkk <- sapply(temp, function(tem) t(tem) %*% solve(Sigma[[k+1]]) %*% tem)
  obj[[k+1]] <- 1/2 * log(det(Sigma[[k+1]])) +(1+P)/(2*(N)) * sum(log(1+xkk))

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

  history <- list(mu = mu, A = A, B = B, D = D, Sigma = Sigma, obj = obj, runtime = c(0,diff(do.call(base::c,runtime))))
  output <- list(method = method,
                 spec = list(N = N, P = P, R = R, r = r),
                 history = history,
                 mu = mu[[length(mu)]],
                 A = A[[length(A)]],
                 B = B[[length(B)]],
                 D = D[[length(D)]],
                 Sigma = Sigma[[length(Sigma)]],
                 obj = obj[[length(obj)]])

  return(new_ORRRR(output))
}
