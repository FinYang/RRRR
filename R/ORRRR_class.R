new_ORRRR <- function(x = list()){
  stopifnot(is.list(x))
  stopifnot(all(c("method","spec","history", "mu", "A", "B", "D", "Sigma", "obj") %in% names(x)))
  structure(x, class = c("ORRRR", "RRRR","RRR"))
}

#' @importFrom stats coef
#' @export
print.ORRRR <- function(x,  digits = max(3L, getOption("digits") - 2L), ...){
  cat("Online Robust Reduced-Rank Regression\n")
  switch(x$method,
         "SMM" = cat("Stochastic Majorization Minimization"),
         "SAA" = cat("Sample Average Approximation"))
  cat("\n------------\n")
  cat("Specifications:\n")
  print(do.call(base::c, x$spec))
  cat("\nCoefficients:\n")
  print.default(coef(x), digits = digits)
}

#' @importFrom stats update
#' @export
update.ORRRR <- function(x){

}
