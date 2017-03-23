## EM algorithm here

#' EM algorithm to fit multivariate Guassian convolution problem.
#'
#' @param x_mat A matrix of data. The rows index the observations
#'     and the columns index the variables.
#' @param s_mat A matrix of variances (NOT standard deviations).
#'     The rows index the observations and the columns index the
#'     variables.
#' @param v_mat A matrix of initial values for the low-rank mixture
#'     covariances.
#' @param pi_vec A vector of initial values of the mixing proportions.
#' @param itermax The maximum number of fixed-point iterations for
#'     the EM to run.
#' @param tol The tolerance for the stopping criterion. The current
#'     stopping criterion is the ratio of successive
#'     likelihoods minus 1.
#' @param print_iter A logical. Should we print the likelihood at
#'     each step (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @export
#'
ultimate_deconvolution <- function(x_mat, s_mat, v_mat, pi_vec,
                                   itermax = 500, tol = 10 ^ -5,
                                   print_iter = FALSE) {

  ## Test input -------------------------------------------
  assertthat::assert_that(is.matrix(x_mat))
  assertthat::assert_that(is.matrix(v_mat))
  assertthat::assert_that(is.matrix(s_mat))
  assertthat::assert_that(is.logical(print_iter))
  assertthat::assert_that(tol >= 0)
  assertthat::assert_that(itermax >= 1)

  R <- ncol(x_mat)
  N <- nrow(x_mat)
  K <- ncol(v_mat)

  assertthat::are_equal(R, nrow(v_mat))
  assertthat::are_equal(N, nrow(s_mat))
  assertthat::are_equal(R, ncol(s_mat))
  assertthat::are_equal(K, length(pi_vec))
  assertthat::assert_that(abs(sum(pi_vec) - 1) < 10 ^ -12)

  llike_current <- dmixlike_cpp(x_mat = x_mat,
                                s_mat = s_mat,
                                v_mat = v_mat,
                                pi_vec = pi_vec,
                                return_log = TRUE)

  iterindex <- 1
  err <- tol + 1
  llike_vec <- llike_current

  while(iterindex <= itermax & err > tol) {

    llike_old <- llike_current

    fout <- em_fix(x_mat = x_mat,
                       s_mat = s_mat,
                       v_mat = v_mat,
                       pi_vec = pi_vec)

    v_mat <- fout$v_mat
    pi_vec <- fout$pi_vec

    llike_current <- dmixlike_cpp(x_mat = x_mat,
                                  s_mat = s_mat,
                                  v_mat = v_mat,
                                  pi_vec = pi_vec,
                                  return_log = TRUE)

    llike_vec <- c(llike_vec, llike_current)
    ## Make sure likelihood increases (within some tolerance)
    assertthat::assert_that(llike_current - llike_old > -10 ^ -12)

    err <- abs(exp(llike_current - llike_old) - 1)

    iterindex <- iterindex + 1
  }

  return(list(pi_vec = pi_vec, v_mat = v_mat, llike_vec = llike_vec))
}

#' Fixed iteration from the EM algorithm.
#'
#' @inheritParams dmixlike
#'
#' @author David Gerard
#'
em_fix <- function(x_mat, s_mat, v_mat, pi_vec) {

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

  ## Get the W matrix (the element-specific mixing proportions)
  llike_mat <- get_llike_mat(x_mat = x_mat,
                             s_mat = s_mat,
                             v_mat = v_mat,
                             pi_vec = pi_vec)

  rowmax <- apply(llike_mat, 1, max)
  ldenom <- log(rowSums(exp(llike_mat - rowmax))) + rowmax

  wmat <- exp(llike_mat - ldenom)
  ## assertthat::assert_that(all(abs(rowSums(wmat) - 1) < 10 ^ -12))

  ## ------------------------------------------------------------
  ## get sigma_kj (the element-specific mixing variances) and
  ## mu_kj (the element-specific mixing means).
  ## Then combine them with the w_kj's to finish up the E-step
  ## ------------------------------------------------------------

  theta_mat <- matrix(NA, nrow = N, ncol = K)
  eta_mat <- matrix(NA, nrow = N, ncol = K)
  for (obs_index in 1:N) {
    x_current <- x_mat[obs_index, ]
    s_current <- s_mat[obs_index, ]
    for (mix_index in 1:K) {
      v_current <- v_mat[, mix_index]
      w_current <- wmat[obs_index, mix_index]
      vsv <- sum(v_current ^ 2 / s_current)
      xsv <- sum(v_current * x_current / s_current)

      sigma2_kj <- 1 / (vsv + 1)
      mu_kj <- sigma2_kj * xsv

      theta_mat[obs_index, mix_index] <- w_current * mu_kj
      eta_mat[obs_index, mix_index] <-
        w_current * (mu_kj ^ 2 + sigma2_kj)
    }
  }

  ## Update the mixing proportions ------------------------------
  pi_new <- colSums(wmat)
  pi_new <- pi_new / sum(pi_new)

  ## Update the rank-1 matrices ---------------------------------
  lincom_s <- 1 / crossprod(eta_mat, 1 / s_mat)
  lincom_xs <- crossprod(theta_mat, x_mat / s_mat)
  v_new <- t(lincom_s * lincom_xs)

  return(list(pi_vec = pi_new, v_mat = v_new))
}
