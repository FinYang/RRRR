new_RRRR <- function(x = list()){
  stopifnot(is.list(x))
  stopifnot(all(c("spec","history", "mu", "A", "B", "D", "Sigma", "obj") %in% names(x)))
  structure(x, class = c("RRRR","RRR"))
}

#' @importFrom stats coef
#' @export
print.RRRR <- function(x,  digits = max(3L, getOption("digits") - 2L), ...){
  cat("Robust Reduced-Rank Regression\n")
  cat("Majorization-Minimization\n")
  cat("------------\n")
  cat("Specifications:\n")
  print(do.call(base::c, x$spec))
  cat("\nCoefficients:\n")
  print.default(coef(x), digits = digits)
}
