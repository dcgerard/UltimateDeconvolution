## Some plotting functions

#' Basic line plot for \code{\link{em_cpp}}.
#'
#' @param llike_vec The y-values
#' @param itermax A positive integer. The x-values are \code{1:itermax}.
#'
#' @author David Gerard
#'
plot_llike <- function(llike_vec, itermax) {
  y <- rep(NA, itermax)
  y[1:length(llike_vec)] <- llike_vec
  graphics::plot(1:itermax, y, xlab = "Iteration", ylab = "Log-likelihood",
                 type = "l", main = "Log-likelihood vs EM Iteration")
  if (length(llike_vec) / itermax >= 0.9) {
    graphics::mtext("Almost There!", at = length(llike_vec))
  }
}
