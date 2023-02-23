#' Update the RRRR/ORRRR type model with addition data
#'
#' \code{update.RRRR} will update online robust reduced-rank regression model with class \code{RRRR}(\code{ORRRR}) using newly added data
#' to achieve online estimation.
#' Estimation methods:
#' \describe{
#' \item{SMM}{Stochastic Majorisation-Minimisation}
#' \item{SAA}{Sample Average Approximation}
#' }
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
#' See \code{?ORRRR} for details about the estimation methods.
#'
#'
#' @param object A model with class \code{RRRR}(\code{ORRRR})
#' @param newy Matrix of dimension N*P, the new data y. The matrix for the response variables. See \code{Detail}.
#' @param newx Matrix of dimension N*Q, the new data x. The matrix for the explanatory variables to be projected. See \code{Detail}.
#' @param newz Matrix of dimension N*R, the new data z. The matrix for the explanatory variables not to be projected. See \code{Detail}.
#' @param addon Integer. The number of data points to be added in the algorithm in each iteration after the first.
#' @param method Character. The estimation method. Either "SMM" or "SAA". See \code{Description}.
#' @param SAAmethod Character. The sub solver used in each iteration when the \code{methid} is chosen to be "SAA". See \code{Detail}.
#' @param ... Additional arguments to function
#' \describe{
#' \item{\code{optim}}{when the \code{method} is "SAA" and the \code{SAAmethod} is "optim"}
#' \item{\code{RRRR}}{when the \code{method} is "SAA" and the \code{SAAmethod} is "MM"}
#' }
#' @param ProgressBar Logical. Indicating if a progress bar is shown during the estimation process.
#' The progress bar requires package \code{lazybar} to work.
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
#' \item{data}{The data used in estimation.}
#' }
#'
#'
#' @seealso \code{ORRRR}, \code{RRRR}, \code{RRR}
#' @examples
#' \donttest{
#' set.seed(2222)
#' data <- RRR_sim()
#' newdata <- RRR_sim(A = data$spec$A,
#'                    B = data$spec$B,
#'                    D = data$spec$D)
#' res <- ORRRR(y=data$y, x=data$x, z = data$z)
#' res <- update(res, newy=newdata$y, newx=newdata$x, newz=newdata$z)
#' res
#' }
#' @author Yangzhuoran Yang
#' @importFrom magrittr %>%
#' @importFrom stats update
#' @export
update.RRRR <- function(object, newy, newx, newz=NULL,
                        addon = object$spec$addon,
                        method = object$method,
                        SAAmethod = object$SAAmethod,
                        ...,
                        ProgressBar = requireNamespace("lazybar")){
  if(!identical(ncol(object$data$y), ncol(newy)))
    stop("The dimension of the new data y is not the same as the data stored in model")
  if(!identical(ncol(object$data$x), ncol(newx)))
    stop("The dimension of the new data x is not the same as the data stored in model")
  if(!is.null(object$data$z) && is.null(newz))
    stop("New data z not supplied.")
  if(is.null(object$data$z) && !is.null(newz))
    stop("New data z ignored.")
  if(!is.null(object$data$z) && !identical(ncol(object$data$z), ncol(newz)))
    stop("The dimension of the new data z is not the same as the data stored in model")
  newy <- rbind(object$data$y, newy)
  newx <- rbind(object$data$x, newx)
  if(!is.null(object$data$z)){
    newz <- rbind(as.matrix(object$data$z), as.matrix(newz))
  } else {
    newz <- NULL
  }

  if(!"ORRRR" %in% class(object)){
    if(is.null(addon)) addon <- 10
    if(is.null(method)) method <- "SMM"
  }
  res <- ORRRR(y=newy, x=newx, z=newz,
               mu = "mu" %in% colnames(coef(object)),
               r = object$spec$r,
               initial_size = object$spec$N+addon,
               addon = addon,
               method = method,
               SAAmethod = SAAmethod,
               ...,
               initial_A = object$A,
               initial_B = object$B,
               initial_D = object$D,
               initial_mu = object$mu,
               initial_Sigma = object$Sigma,
               ProgressBar = ProgressBar,
               return_data = TRUE)
  # res$data$y <- rbind(object$data$y, y)
  # res$data$x <- rbind(object$data$x, x)
  # if(!is.null(object$data$z)){
  #   res$data$z <- rbind(object$data$z, z)
  # } else {
  #   res$data$z <- NULL
  # }

  res$history$mu <- c(object$history$mu, res$history$mu)
  res$history$A <- c(object$history$A, res$history$A)
  res$history$B <- c(object$history$B, res$history$B)
  res$history$D <- c(object$history$D, res$history$D)
  res$history$Sigma <- c(object$history$Sigma, res$history$Sigma)
  res$history$obj <- c(object$history$obj, res$history$obj)
  res$history$runtime <- c(object$history$runtime, res$history$runtime+object$history$runtime[[length(object$history$runtime)]])


  return(new_ORRRR(res))
}
