new_RRRR <- function(x = list()){
  stopifnot(is.list(x))
  stopifnot(all(c("spec","history", "mu", "A", "B", "D", "Sigma", "obj") %in% names(x)))
  structure(x, class = c("RRRR","RRR"))
}

#' @importFrom stats coef
#' @export
print.RRRR <- function(x,  digits = max(3L, getOption("digits") - 2L), ...){
  cat("Robust Reduced-Rank Regression\n")
  cat("------\n")
  cat("Majorisation-Minimisation\n")
  cat("------------\n")
  cat("Specifications:\n")
  print(do.call(base::c, x$spec))
  cat("\nCoefficients:\n")
  print(coef(x), digits = digits)
}

#' Plot Objective value of a Robust Reduced-Rank Regression
#'
#' @param x An RRRR object.
#' @param aes_x Either "iteration" or "runtime". The x axis in the plot.
#' @param xlog10 Logical, indicates whether the scale of x axis is log 10 transformed.
#' @param ... Additional argument to \code{ggplot2}.
#' @return An ggplot2 object
#' @author Yangzhuoran Fin Yang
#' @examples
#' set.seed(2222)
#' data <- RRR_sim()
#' res <- RRRR(y=data$y, x=data$x, z = data$z)
#' plot(res)
#' @export
plot.RRRR <- function(x,
                      aes_x = c("iteration", "runtime"),
                      xlog10 = TRUE, ...){

  plot_data <- data.frame(runtime = cumsum(x$history$runtime),
                          obj = x$history$obj,
                          iteration = seq_along(x$history$obj))
  if(xlog10) plot_data$runtime[[1]] <- 0.001
  output <- ggplot2::ggplot(plot_data, ...) +
    ggplot2::geom_line(ggplot2::aes_string(x=aes_x[[1]], y="obj")) +
    ggplot2::ylab("Objective value")
  if(xlog10) output <- output + ggplot2::coord_trans(x="log10")
  return(output)
}
