new_ORRRR <- function(x = list()){
  stopifnot(is.list(x))
  stopifnot(all(c("method","SAAmethod","spec","history", "mu", "A", "B", "D", "Sigma", "obj") %in% names(x)))
  structure(x, class = c("ORRRR", "RRRR","RRR"))
}

#' @importFrom stats coef
#' @export
print.ORRRR <- function(x,  digits = max(3L, getOption("digits") - 2L), ...){
  cat("Online Robust Reduced-Rank Regression")
  cat("\n------\n")
  switch(x$method,
         "SMM" = cat("Stochastic Majorisation-Minimisation"),
         "SAA" = cat("Sample Average Approximation"))
  if(x$method == "SAA"){
    cat("\nSub solver: ")
    switch(x$SAAmethod,
           "optim" = cat("stats::optim"),
           "MM" = cat("Majorisation Minimisation"))
  }
  cat("\n------------\n")
  cat("Specifications:\n")
  print(do.call(base::c, x$spec))
  cat("\nCoefficients:\n")
  print(coef(x), digits = digits)
}


