## ultimate_deconvolution

#' Multivariate normal density at mean 0 and a
#' covariance that is the sum of a rank-1 matrix and a diagonal matrix.
#'
#' @param x The data.
#' @param v The square-root of the rank-1 matrix.
#' @param s_diag The diagonal elements of the diagonal component of the
#'     covariance matrix.
#' @param mu The mean. Defaults to 0.
#' @param log A logical. Should we return the log-density (\code{TRUE}) or the
#'     density (\code{FALSE})?
#' @author David Gerard
#'
dr1_norm <- function(x, v, s_diag, mu = rep(0, length(x)), log = FALSE) {

  ## Check input
  R <- length(x)
  assertthat::are_equal(length(R), 1)
  assertthat::assert_that(!is.null(R))
  assertthat::assert_that(is.numeric(x))
  assertthat::assert_that(is.numeric(v))
  assertthat::assert_that(is.numeric(s_diag))
  assertthat::assert_that(is.numeric(mu))
  assertthat::assert_that(is.logical(log))
  assertthat::are_equal(R, length(v))
  assertthat::are_equal(R, length(s_diag))
  assertthat::are_equal(R, length(mu))


  vsv <- sum(v ^ 2 / s_diag)
  xmusv <- sum((x - mu) * v / s_diag)
  xmusxmu <- sum((x - mu) ^ 2 / s_diag)

  llike <- (-R * log(2 * pi) - log(1 + vsv) - sum(log(s_diag)) -
              xmusxmu + (xmusv ^ 2) / (1 + vsv)) / 2

  if (log) {
    return(llike)
  } else {
    return(exp(llike))
  }
}
