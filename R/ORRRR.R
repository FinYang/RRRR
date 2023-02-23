#' Online Robust Reduced-Rank Regression
#'
#'
#' Online robust reduced-rank regression with two major estimation methods:
#' \describe{
#' \item{SMM}{Stochastic Majorisation-Minimisation}
#' \item{SAA}{Sample Average Approximation}
#' }
#'
#'
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
#' the innovation's distribution.
#'
#' The algorithm is online in the sense that the data is continuously incorporated
#' and the algorithm can update the parameters accordingly. See \code{?update.RRRR} for more details.
#'
#' At each iteration of SAA, a new realisation of the parameters is achieved by
#' solving the minimisation problem of the sample average of
#' the desired objective function using the data currently incorporated.
#' This can be computationally expensive when the objective function is highly nonconvex.
#' The SMM method overcomes this difficulty by replacing the objective function
#' by a well-chosen majorising surrogate function which can be much easier to optimise.
#'
#' SMM method is robust in the sense that it assumes a heavy-tailed Cauchy distribution
#' for the innovations.
#'
#'
#' @param y Matrix of dimension N*P. The matrix for the response variables. See \code{Detail}.
#' @param x Matrix of dimension N*Q. The matrix for the explanatory variables to be projected. See \code{Detail}.
#' @param z Matrix of dimension N*R. The matrix for the explanatory variables not to be projected. See \code{Detail}.
#' @param mu Logical. Indicating if a constant term is included.
#' @param r Integer. The rank for the reduced-rank matrix \eqn{AB'}. See \code{Detail}.
#' @param initial_size Integer. The number of data points to be used in the first iteration.
#' @param addon Integer. The number of data points to be added in the algorithm in each iteration after the first.
#' @param method Character. The estimation method. Either "SMM" or "SAA". See \code{Description} and \code{Detail}.
#' @param SAAmethod Character. The sub solver used in each iteration when the \code{method} is chosen to be "SAA". See \code{Detail}.
#' @param ... Additional arguments to function
#' \describe{
#' \item{\code{optim}}{when the \code{method} is "SAA" and the \code{SAAmethod} is "optim"}
#' \item{\code{RRRR}}{when the \code{method} is "SAA" and the \code{SAAmethod} is "MM"}
#' }
#' @param initial_mu Vector of length P. The initial value for constant \eqn{\mu} See \code{Detail}.
#' @param initial_A Matrix of dimension P*r. The initial value for matrix \eqn{A}. See \code{Detail}.
#' @param initial_B Matrix of dimension Q*r. The initial value for matrix \eqn{B}. See \code{Detail}.
#' @param initial_D Matrix of dimension P*R. The initial value for matrix \eqn{D}. See \code{Detail}.
#' @param initial_mu Matrix of dimension P*1. The initial value for the constant \eqn{mu}. See \code{Detail}.
#' @param initial_Sigma Matrix of dimension P*P. The initial value for matrix Sigma. See \code{Detail}.
#' @param ProgressBar Logical. Indicating if a progress bar is shown during the estimation process.
#' The progress bar requires package \code{lazybar} to work.
#' @param return_data Logical. Indicating if the data used is return in the output.
#' If set to \code{TRUE}, \code{update.RRRR} can update the model by simply provide new data.
#' Set to \code{FALSE} to save output size.
#'
#' @return A list of the estimated parameters of class \code{ORRRR}.
#' \describe{
#' \item{method}{The estimation method being used}
#' \item{SAAmethod}{If SAA is the major estimation method, what is the sub solver in each iteration.}
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
#' @seealso \code{update.RRRR}, \code{RRRR}, \code{RRR}
#' @examples
#' \donttest{
#' set.seed(2222)
#' data <- RRR_sim()
#' res <- ORRRR(y=data$y, x=data$x, z = data$z)
#' res
#' }
#' @author Yangzhuoran Yang
#' @importFrom magrittr %>%
#' @export
ORRRR <- function(y, x, z = NULL, mu = TRUE, r = 1,
                  initial_size = 100, addon = 10,
                  method = c("SMM", "SAA"),
                  SAAmethod = c("optim", "MM"),
                  ...,
                  initial_A = matrix(rnorm(P*r), ncol =  r),
                  initial_B = matrix(rnorm(Q*r), ncol =  r),
                  initial_D = matrix(rnorm(P*R), ncol =  R),
                  initial_mu = matrix(rnorm(P)),
                  initial_Sigma = diag(P),
                  ProgressBar = requireNamespace("lazybar"),
                  return_data = TRUE){

  if (ProgressBar && !requireNamespace("lazybar", quietly = TRUE)) {
    stop("Package \"lazybar\" needed for progress bar to work. Please install it.",
         call. = FALSE)
  }
  method <- method[[1]]
  if(!method %in% c("SMM", "SAA")) stop("Unrecognised method")
  if(method == "SAA")  SAAmethod <- SAAmethod[[1]] else SAAmethod <- "NULL"
  if(method == "SAA" && !SAAmethod %in% c("optim", "MM")) stop("Unrecognised SAAmethod")

  if(SAAmethod == "MM"){
    RRRR_argument <- list(...)
    if(is.null(RRRR_argument$itr)) RRRR_argument$itr <- 10
    if(is.null(RRRR_argument$earlystop)) RRRR_argument$earlystop <- 1e-4

  }
  if(return_data){
    returned_data <- list(y=y, x=x, z=z)
  } else {
    returned_data <- NULL
  }

  N <- nrow(y)

  P <- ncol(y)
  Q <- ncol(x)
  if(NCOL(initial_A) != r)
    stop("Mismatched dimension. The column number of initial_A is not the same as r.")
  if(NCOL(initial_B) != r)
    stop("Mismatched dimension. The column number of initial_B is not the same as r.")
  if(NROW(initial_A) != P)
    stop("Mismatched dimension. The row number of initial_A is not the same as P.")
  if(NROW(initial_mu) != P)
    stop("Mismatched dimension. The row number (length) of initial_mu is not the same as P.")
  if(NROW(initial_B) != Q)
    stop("Mismatched dimension. The row number of initial_B is not the same as Q.")
  # check if z is NULL
  # add column of ones for mu
  # if SAA-MM, pass z as it is
  if(!is.null(z)){
    z <- as.matrix(z)
    R <- ncol(z)
    if(NCOL(initial_D) != R)
      stop("Mismatched dimension. The column number of initial_D is not the same as the column number of variable z.")
    if(NROW(initial_D) != P)
      stop("Mismatched dimension. The row number of initial_D is not the same as P.")

    if(mu){
      if(SAAmethod != "MM"){
        z <- cbind(z, 1)
        initial_D <- cbind(initial_D, initial_mu)
      }
    }
    znull <- FALSE
  } else {
    R <- 0
    znull <- TRUE
    if (mu){
      if(SAAmethod != "MM")
        z <- matrix(rep(1, N))
      initial_D <- initial_mu
    }
  }
  muorz <- mu || !znull

  if(nrow(y) != nrow(x)){
    if(!is.null(z) && nrow(y) != ncol(z))
      stop("The numbers of realizations are not consistant in the inputs.")
  }



  yy <- y
  xx <- x
  if(SAAmethod != "MM"|| (SAAmethod=="MM" && !znull))
    zz <- z

  # initialise loop
  # SMM and MM are the same
  if(method=="SMM" || SAAmethod=="MM"){
    A <- list()
    B <- list()
    Pi <- list()
    D <- list()
    MM_mu <- list()
    Sigma <- list()
    A[[1]] <- initial_A
    B[[1]] <- initial_B
    Pi[[1]] <- A[[1]] %*% t(B[[1]])
    if(muorz)
      D[[1]] <- initial_D
    if(SAAmethod == "MM")
      MM_mu[[1]] <- initial_mu
    Sigma[[1]] <- initial_Sigma


    ybar <- list()
    xbar <- list()
    zbar <- list()
    Mbar <- list()
  } else if(method=="SAA"){
    if(SAAmethod == "optim"){
      make_symm <- function(m) {
        m[upper.tri(m)] <- t(m)[upper.tri(m)]
        return(m)
      }

      para <- list()
      if(muorz){
        para[[1]] <- c(initial_A, initial_B, initial_D,
                       initial_Sigma[lower.tri(initial_Sigma, diag = TRUE)])
      } else {
        para[[1]] <- c(initial_A, initial_B,
                       initial_Sigma[lower.tri(initial_Sigma, diag = TRUE)])
      }

    }
  }
  xk <- list()
  wk <- list()

  # progress bar and run time
  itr <- (N-initial_size)/addon +1
  obj <- numeric(itr+1)

  runtime <- vector("list",itr+1)
  if(ProgressBar)
    pb <- lazybar::lazyProgressBar(itr, method = "drift")
  runtime[[1]] <- Sys.time()
  # loop
  for (k in seq_len(itr)) {
    # assign different data each itreation
    if(k==1){
      y <- t(yy[seq(1, k*initial_size), ])
      x <- t(xx[seq(1, k*initial_size), ])
      if(muorz && (SAAmethod != "MM" || (SAAmethod=="MM" && !znull)))
        z <- t(zz[seq(1, k*initial_size), ])
    } else {
      if(floor(itr) < itr && k==floor(itr)){
        y <- t(yy)
        x <- t(xx)
        if(muorz && (SAAmethod != "MM" || (SAAmethod=="MM" && !znull)))
          z <- t(zz)
      } else {
        y <- t(yy[seq(1, initial_size+(k-1)*addon),])
        x <- t(xx[seq(1, initial_size+(k-1)*addon),])
        if(muorz && (SAAmethod != "MM" || (SAAmethod=="MM" && !znull)))
          z <- t(zz[seq(1, initial_size+(k-1)*addon),])
      }
    }
    N <- ncol(y)

    # SMM core
    if(method == "SMM"){
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

    } else if(method=="SAA"){
      # SAA-optim core
      if(SAAmethod == "optim"){
        if(muorz){
          ne_log_likihood_loss <- function(para){
            A <- matrix(para[1:(P*r)], nrow = P)
            B <- matrix(para[(P*r+1):(2*P*r)], nrow = P)
            D <- matrix(para[(2*P*r+1):(P*r*2+length(initial_D))], nrow = P)
            Sigma <- matrix(nrow = P, ncol = P)
            Sigma[lower.tri(Sigma,diag=TRUE)] <- para[(P*r*2+length(initial_D)+1):(length(para))]
            Sigma <- make_symm(Sigma)
            if(!matrixcalc::is.positive.definite(Sigma))
              return(Inf)
            Pi <- A %*% t(B)
            temp <- t(y - Pi %*% x - D %*% z) %>% split(seq_len(nrow(.)))
            xk <- sapply(temp, function(tem) t(tem) %*% solve(Sigma) %*% tem)

            return(1/2 * log(det(Sigma)) +(1+P)/(2*(N-2)) * sum(log(1+xk)))

          }
        } else {
          ne_log_likihood_loss <- function(para){
            A <- matrix(para[1:(P*r)], nrow = P)
            B <- matrix(para[(P*r+1):(2*P*r)], nrow = P)
            Sigma <- matrix(nrow = P, ncol = P)
            Sigma[lower.tri(Sigma,diag=TRUE)] <- para[(2*P*r+1):(length(para))]
            Sigma <- make_symm(Sigma)
            if(!matrixcalc::is.positive.definite(Sigma))
              return(Inf)
            Pi <- A %*% t(B)
            temp <- t(y - Pi %*% x ) %>% split(seq_len(nrow(.)))
            xk <- sapply(temp, function(tem) t(tem) %*% solve(Sigma) %*% tem)

            return(1/2 * log(det(Sigma)) +(1+P)/(2*(N)) * sum(log(1+xk)))

          }
        }


        sub_res <- stats::optim(para[[k]], ne_log_likihood_loss, ...)
        para[[k+1]] <- sub_res$par


        obj[[k+1]] <- sub_res$value
      } else if(SAAmethod == "MM"){
        # SAA-MM core
        if(!znull){
          # warning(str(D[[k]]))
          # warning(str(z))
          # warning(str(x))
          sub_res <- RRRR(y=t(y), x=t(x), z = t(z), mu = mu, r=r,
                          initial_A = A[[k]],
                          initial_B = B[[k]],
                          initial_D = D[[k]],
                          initial_mu = MM_mu[[k]],
                          initial_Sigma = Sigma[[k]],
                          itr = RRRR_argument$itr,
                          earlystop = RRRR_argument$earlystop)
          D[[k+1]] <- sub_res$D

        } else {
          sub_res <- RRRR(y=t(y), x=t(x), mu = mu, r=r,
                          initial_A = A[[k]],
                          initial_B = B[[k]],
                          initial_D = NULL,
                          initial_mu = MM_mu[[k]],
                          initial_Sigma = Sigma[[k]],
                          itr = RRRR_argument$itr,
                          earlystop = RRRR_argument$earlystop)
        }
        MM_mu[[k+1]] <- sub_res$mu
        A[[k+1]] <- sub_res$A
        B[[k+1]] <- sub_res$B
        # D[[k+1]] <- sub_res$D
        Sigma[[k+1]] <- sub_res$Sigma

        obj[[k+1]] <- sub_res$obj
      }

      # for SAA
      # assign the first obj value
      if(k==1){
        if(SAAmethod != "MM"){
          initial_Pi <- initial_A %*% t(initial_B)
          if(muorz){
            temp <- t(y - initial_Pi %*% x - initial_D %*% z) %>% split(seq_len(nrow(.)))
          } else {
            temp <- t(y - initial_Pi %*% x ) %>% split(seq_len(nrow(.)))
          }
          xk <- sapply(temp, function(tem) t(tem) %*% solve(initial_Sigma) %*% tem)

          obj[[1]] <- 1/2 * log(det(initial_Sigma)) +(1+P)/(2*(N)) * sum(log(1+xk))
        } else {
          initial_Pi <- initial_A %*% t(initial_B)
          if(mu && znull){
            # stop(paste("y",dim(y),"PI", dim(initial_Pi),"x", dim(x), "mu", dim(initial_mu)))
            temp <- t(y - initial_Pi %*% x - initial_mu %*% matrix(rep(1, ncol(x)), nrow = 1)) %>% split(seq_len(nrow(.)))
          } else if(!mu && !znull){
            temp <- t(y - initial_Pi %*% x - initial_D %*% z) %>% split(seq_len(nrow(.)))

          } else {
            temp <- t(y - initial_Pi %*% x ) %>% split(seq_len(nrow(.)))
          }
          xk <- sapply(temp, function(tem) t(tem) %*% solve(initial_Sigma) %*% tem)

          obj[[1]] <- 1/2 * log(det(initial_Sigma)) +(1+P)/(2*(N)) * sum(log(1+xk))
        }
      }
    }
    if(ProgressBar)
      pb$tick()$print()
    runtime[[k+1]] <- Sys.time()
  }

  # calculate obj value for SMM
  if(method == "SMM"){
    if(muorz){
      temp <- t(y - Pi[[k+1]] %*% x - D[[k+1]] %*% z) %>% split(seq_len(nrow(.)))
    } else {
      temp <- t(y - Pi[[k+1]] %*% x) %>% split(seq_len(nrow(.)))
    }
    xkk <- sapply(temp, function(tem) t(tem) %*% solve(Sigma[[k+1]]) %*% tem)
    obj[[k+1]] <- 1/2 * log(det(Sigma[[k+1]])) +(1+P)/(2*(N)) * sum(log(1+xkk))

  } else if(method == "SAA"){
    # recover parameters from SAA-optim
    if(SAAmethod == "optim"){
      A <- lapply(para, function(para) matrix(para[1:(P*r)], nrow = P))
      B <- lapply(para, function(para) matrix(para[(P*r+1):(2*P*r)], nrow = P))
      if(muorz){

        D <- lapply(para, function(para) matrix(para[(2*P*r+1):(P*r*2+length(initial_D))], nrow = P))
        Sigma <- lapply(para,
                        function(para){
                          Sigma <- matrix(nrow = P, ncol = P)
                          Sigma[lower.tri(Sigma,diag=TRUE)] <- para[(P*r*2+length(initial_D)+1):(length(para))]
                          Sigma <- make_symm(Sigma)
                          return(Sigma)
                        })
      } else {
        Sigma <- lapply(para,
                        function(para){
                          Sigma <- matrix(nrow = P, ncol = P)
                          Sigma[lower.tri(Sigma,diag=TRUE)] <- para[(2*P*r+1):(length(para))]
                          Sigma <- make_symm(Sigma)
                          return(Sigma)
                        })
      }
    }
  }
  # seperate mu and D
  # when it's not SAA-MM
  if(SAAmethod != "MM"){
    if(mu){
      mu <- lapply(D[sapply(D, function(x) !is.null(x))], function(x) x[,ncol(x)])
      if(!znull)
        D <- lapply(D[sapply(D, function(x) !is.null(x))], function(x) x[,seq_len(ncol(x)-1)])
    } else {
      mu <- NULL
    }
  } else {
    # recover mu for SAA-MM
    if(mu){
      mu <- MM_mu
    } else {
      mu <- NULL
    }
  }
  if(znull){
    D <- NULL
  }

  # arrange output
  history <- list(mu = mu, A = A, B = B, D = D, Sigma = Sigma, obj = obj, runtime = c(0,diff(do.call(base::c,runtime))))
  output <- list(method = method,
                 SAAmethod = SAAmethod,
                 spec = list(N = N, P = P, R = R, r = r, initial_size = initial_size, addon = addon),
                 history = history,
                 mu = mu[[length(mu)]],
                 A = A[[length(A)]],
                 B = B[[length(B)]],
                 D = D[[length(D)]],
                 Sigma = Sigma[[length(Sigma)]],
                 obj = obj[[length(obj)]],
                 data = returned_data)

  return(new_ORRRR(output))
}
