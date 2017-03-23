context("Densities")

test_that("Normal Density is Accurate", {

  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    set.seed(17)
    R <- 11
    x <- stats::rnorm(R)
    s_diag <- stats::rchisq(R, df = 5) / 5
    v <- stats::rnorm(R)
    mu <- rep(0, length(x))
    log <- TRUE

  medense <- dr1_norm(x = x, v = v, s_diag = s_diag,
                      mu = mu, log = log)

  true_cov <- tcrossprod(v) + diag(s_diag)
  themdense <- mvtnorm::dmvnorm(x = x, mean = mu,
                                sigma = true_cov, log = log)

  testthat::expect_equal(medense, themdense)

  } else {
    skip("mvtnorm not installed")
  }
}
)

test_that("dmixlike runs OK", {

  if (requireNamespace("mixAK", quietly = TRUE)) {
    set.seed(45)
    K <- 3
    R <- 11
    N <- 7

    x_mat <- matrix(stats::rnorm(N * R), nrow = N)
    s_mat <- matrix(stats::rchisq(N * R, df = 5), nrow = N) / 5
    pi_vec <- stats::runif(K)
    pi_vec <- pi_vec / sum(pi_vec)
    v_mat <- matrix(stats::rnorm(R * K), nrow = R)

    me_like <- dmixlike(x_mat = x_mat, s_mat = s_mat, pi_vec = pi_vec,
                        v_mat = v_mat, log = TRUE)

    ## Set input to be same as in mixAK package
    llike_vec <- rep(NA, length = N)
    mean_mat <- matrix(0, nrow = K, ncol = R)
    for (obs_index in 1:N) {
      Sigma_list <- list()
      for (mix_index in 1:K) {
        Sigma_list[[mix_index]] <- diag(s_mat[obs_index, ]) +
          tcrossprod(v_mat[, mix_index])
      }
      llike <- mixAK::dMVNmixture(x = x_mat[obs_index, ],
                                  weight = pi_vec,
                                  mean = mean_mat,
                                  Sigma = Sigma_list,
                                  log = TRUE)
      llike_vec[obs_index] <- llike
    }
    them_like <- sum(llike_vec)
    expect_equal(me_like, them_like)
  } else {
    skip("mixAK not installed")
  }
}
)

test_that("rmixmultnorm works", {
  set.seed(3456)
  K <- 1
  R <- 11
  N <- 10000
  s_mat <- matrix(0.1 , nrow = N, ncol = R)
  pi_vec <- stats::runif(K)
  pi_vec <- pi_vec / sum(pi_vec)
  v_mat <- matrix(stats::rnorm(R * K), nrow = R)

  xout <- rmixmultnorm(s_mat = s_mat, v_mat = v_mat, pi_vec = pi_vec)

  sample_cov <- crossprod(xout) / N

  ## This should be small
  expect_true(max(sample_cov - (tcrossprod(v_mat) + diag(0.1, nrow = R))) < 0.02)
}
)


test_that("dr1_norm (R version) and dnorm_rank1 (C++ version) return same llike", {
  set.seed(17)
  R <- 11
  x <- stats::rnorm(R)
  s_diag <- stats::rchisq(R, df = 5) / 5
  v <- stats::rnorm(R)
  mu <- rep(0, length(x))
  log <- TRUE

  rllike <- dr1_norm(x = x, v = v, s_diag = s_diag,
                     mu = mu, log = log)
  cppllike <- dnorm_rank1(x = x, v = v, s_diag = s_diag, mu = mu, return_log = log)

  expect_equal(rllike, cppllike)
}
)
