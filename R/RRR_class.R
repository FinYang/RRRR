new_RRR <- function(x = list()){
  stopifnot(is.list(x))
  stopifnot(all(c("spec", "mu", "A", "B", "D", "Sigma") %in% names(x)))
  structure(x, class = "RRR")
}

#' @importFrom stats coef
#' @export
coef.RRR <- function(object, ...){
  coefficient <- with(object, cbind(mu, A, B, D, Sigma))
  r <- object$spec$r
  P <- object$spec$P
  R <- object$spec$R
  name <- with(object, c(rep("mu", !is.null(mu)),
                         paste0("A", ifelse(r>1, seq(r), "")),
                         paste0("B", ifelse(r>1, seq(r), "")),
                         rep(paste0("D", ifelse(R>1, seq(R), "")), !is.null(D)),
                         paste0("Sigma", seq(P))
  )
  )
  colnames(coefficient) <- name
  return(coefficient)
}

#' @importFrom stats coef
#' @export
print.RRR <- function(x,  digits = max(3L, getOption("digits") - 2L), ...){
  cat("Reduced-Rank Regression\n")
  cat("------------\n")
  cat("Specifications:\n")
  print(do.call(base::c, x$spec))
  cat("\nCoefficients:\n")
  print.default(coef(x), digits = digits)
}

# summary.RRR <- function(object, ...){
#   out <- list(spec = do.call(base::c, object$spec), coef = coef(object))
#   class(out) <- "summary.RRR"
#   return(out)
# }
#
# print.summary.RRR <- function(x, digits = max(3L, getOption("digits") - 2L), ...){
#   cat("Reduced-Rank Regression\n")
#   cat("------------\n")
#   cat("Specifications:\n")
#   print(x$spec)
#   cat("\nCoefficients:\n")
#   print.default(x$coef, digits = digits)
# }

#' @export
print.RRR_data <- function(x, digits = max(3L, getOption("digits") - 2L), ...){
  coefficient <- with(x$spec, cbind(mu, A, B, D, Sigma))
  r <- x$spec$r
  P <- x$spec$P
  R <- x$spec$R
  name <- with(x$spec, c(rep("mu", !is.null(mu)),
                         paste0("A", ifelse(r>1, seq(r), "")),
                         paste0("B", ifelse(r>1, seq(r), "")),
                         rep(paste0("D", ifelse(R>1, seq(R), "")), !is.null(D)),
                         paste0("Sigma", seq(P))
  )
  )
  spec <- with(x$spec, c(N = N, P =P, Q = Q,R = R,r= r))
  colnames(coefficient) <- name
  cat("Simulated Data for Reduced-Rank Regression\n")
  cat("------------\n")
  cat("Specifications:\n")
  print(spec)
  print.default(coefficient, digits = digits)
}