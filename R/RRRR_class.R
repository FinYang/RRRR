new_RRRR <- function(x = list()){
  stopifnot(is.list(x))
  stopifnot(all(c("spec","history", "mu", "A", "B", "D", "Sigma", "obj") %in% names(x)))
  structure(x, class = c("RRRR","RRR"))
}
