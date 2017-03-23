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
#' @export
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


#' Calculates the gaussian mixture density.
#'
#' This function assumes that the mixing means are all zeros and
#' the mixing covariances are rank-1 matrices. Each observation
#' has its own independent noise (variances collected in
#' \code{s_mat}).
#'
#' @param x_mat A matrix of data. The rows index the observations and
#'     the columns index the variables.
#' @param pi_vec A vector of mixture proportions.
#' @param v_mat A matrix of square-roots of the rank-1 mixing
#'     covariance matrices. The columns index the mixing
#'     components, the rows index the variables.
#' @param s_mat A matrix of variances (NOT standard deviations).
#'     The rows index the observations and the columns index the
#'     variables.
#' @param log A logical. Should we return the log-density
#'     (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @export
#'
dmixlike <- function(x_mat, s_mat, v_mat, pi_vec, log = FALSE) {

  ## Test input -------------------------------------------
  assertthat::assert_that(is.matrix(x_mat))
  assertthat::assert_that(is.matrix(v_mat))
  assertthat::assert_that(is.matrix(s_mat))
  assertthat::assert_that(is.logical(log))

  R <- ncol(x_mat)
  N <- nrow(x_mat)
  K <- ncol(v_mat)

  assertthat::are_equal(R, nrow(v_mat))
  assertthat::are_equal(N, nrow(s_mat))
  assertthat::are_equal(R, ncol(s_mat))
  assertthat::are_equal(K, length(pi_vec))
  assertthat::assert_that(abs(sum(pi_vec) - 1) < 10 ^ -12)

  ## Get matrix of log-likelihood values ------------------
  llike_mat <- get_llike_mat(x_mat = x_mat,
                             s_mat = s_mat,
                             v_mat = v_mat,
                             pi_vec = pi_vec)

  ## log-sum-exponential trick for each row, then sum.
  rowmax <- apply(llike_mat, 1, max)
  llike_rowtot <- log(rowSums(exp(llike_mat - rowmax))) + rowmax
  llike_tot <- sum(llike_rowtot)

  if (log) {
    return(llike_tot)
  } else {
    return(exp(llike_tot))
  }
}


#' Get's a matrix of likeklihood values.
#'
#' Element (i, j) of the returned matrix is the log of
#' pi_j N(x_i | 0, v_j v_j^T + S_i)
#'
#' @inheritParams dmixlike
#'
#' @author David Gerard
#'
get_llike_mat <- function(x_mat, s_mat, v_mat, pi_vec) {

  ## Test input -------------------------------------------
  assertthat::assert_that(is.matrix(x_mat))
  assertthat::assert_that(is.matrix(v_mat))
  assertthat::assert_that(is.matrix(s_mat))

  R <- ncol(x_mat)
  N <- nrow(x_mat)
  K <- ncol(v_mat)

  assertthat::are_equal(R, nrow(v_mat))
  assertthat::are_equal(N, nrow(s_mat))
  assertthat::are_equal(R, ncol(s_mat))
  assertthat::are_equal(K, length(pi_vec))
  assertthat::assert_that(abs(sum(pi_vec) - 1) < 10 ^ -12)

  ## Get matrix of log-likelihood values ------------------
  llike_mat <- matrix(NA, nrow = N, ncol = K)

  for (obs_index in 1:N) {
    x_current <- x_mat[obs_index, ]
    s_current <- s_mat[obs_index, ]
    for (mix_index in 1:K) {
      v_current <- v_mat[, mix_index]
      pi_current <- pi_vec[mix_index]
      llike <- dr1_norm(x = x_current, v = v_current,
                        s_diag = s_current, log = TRUE)
      llike_mat[obs_index, mix_index] <- llike + log(pi_current)
    }
  }

  return(llike_mat)
}



#' Random draw from mixture of multivariate normals.
#'
#' It is assumed that the mixing covariances are rank-1 (with their square roots in
#' \code{v_mat}) and the observation covariances are diagonal (with their variances in s_mat).
#'
#' The number of observations drawn is equal to \code{nrow(s_mat)}.
#'
#' @inheritParams dmixlike
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso dmixlike
rmixmultnorm <- function(s_mat, v_mat, pi_vec) {
  ## Check input -------------------------------------------------------------
  assertthat::assert_that(is.matrix(v_mat))
  assertthat::assert_that(is.matrix(s_mat))

  R <- nrow(v_mat)
  K <- ncol(v_mat)
  N <- nrow(s_mat)

  assertthat::are_equal(length(pi_vec), K)
  assertthat::are_equal(ncol(s_mat), R)

  x_mat <- matrix(NA, nrow = N, ncol = R)

  for (index in 1:N) {
    aj <- stats::rnorm(1)
    zj <- sample(x = 1:K, prob = pi_vec, size = 1)
    x_mat[index, ] <- stats::rnorm(n = R, mean = aj * v_mat[, zj], sd = sqrt(s_mat[index, ]))
  }

  return(x_mat)
}
