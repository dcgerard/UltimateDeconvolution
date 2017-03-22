context("Densities")

test_that("Normal Density is Accurate", {

  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    set.seed(17)
    R <- 11
    x <- stats::rnorm(R)
    s_diag <- stats::rchisq(R, df = 5)
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
