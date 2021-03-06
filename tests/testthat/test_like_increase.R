context("Likelihood Increases")

test_that("Likelihood Increases", {

  ## Generate fake data --------------------------------
  set.seed(823)
  K <- 3
  R <- 11
  N <- 7

  x_mat <- matrix(stats::rnorm(N * R), nrow = N)
  s_mat <- matrix(stats::rchisq(N * R, df = 5), nrow = N)
  pi_vec <- stats::runif(K)
  pi_vec <- pi_vec / sum(pi_vec)
  v_mat <- matrix(stats::rnorm(R * K), nrow = R)

  llike_current <- dmixlike(x_mat = x_mat,
                            s_mat = s_mat,
                            v_mat = v_mat,
                            pi_vec = pi_vec,
                            log = TRUE)

  llike_vec <- llike_current
  itermax <- 10
  for (index in 1:itermax) {
    fout <- em_fix(x_mat = x_mat,
                   s_mat = s_mat,
                   v_mat = v_mat,
                   pi_vec = pi_vec)

    v_mat <- fout$v_mat
    pi_vec <- fout$pi_vec

    llike_current <- dmixlike(x_mat = x_mat,
                              s_mat = s_mat,
                              v_mat = v_mat,
                              pi_vec = pi_vec,
                              log = TRUE)

    llike_vec <- c(llike_vec, llike_current)
  }

  expect_true(all(llike_vec[1:(itermax - 1)] <= llike_vec[2:itermax]))

}
)

test_that("ultimate_deconvolution will run", {
  set.seed(83)
  K <- 3
  R <- 11
  N <- 7

  x_mat <- matrix(stats::rnorm(N * R), nrow = N)
  s_mat <- matrix(stats::rchisq(N * R, df = 5), nrow = N)
  pi_vec <- stats::runif(K)
  pi_vec <- pi_vec / sum(pi_vec)
  v_mat <- matrix(stats::rnorm(R * K), nrow = R)

  rout <- em_r(x_mat = x_mat, s_mat = s_mat, v_mat = v_mat, pi_vec = pi_vec, itermax = 10)

  cppout <- em_cpp(x_mat = x_mat, s_mat = s_mat, v_mat = v_mat, pi_vec = pi_vec,
                   itermax = 10)

  expect_equal(rout[[1]], cppout[[1]])
  expect_equal(rout[[2]], cppout[[2]])
  expect_equal(rout[[3]], cppout[[3]])
}
)

# test_that("em_fix_cpp and em_fix return same values", {
#   ## Generate fake data --------------------------------
#   set.seed(823)
#   K <- 3
#   R <- 11
#   N <- 7
#
#   x_mat <- matrix(stats::rnorm(N * R), nrow = N)
#   s_mat <- matrix(stats::rchisq(N * R, df = 5), nrow = N)
#   pi_vec <- stats::runif(K)
#   pi_vec <- pi_vec / sum(pi_vec)
#   v_mat <- matrix(stats::rnorm(R * K), nrow = R)
#
#   for (index in 1:5) {
#     r_fout <- em_fix(x_mat = x_mat, s_mat = s_mat, v_mat = v_mat, pi_vec = pi_vec)
#     v_r <- r_fout$v_mat
#     pi_r <- r_fout$pi_vec
#     cpp_fout <- em_fix(x_mat = x_mat, s_mat = s_mat, v_mat = v_mat, pi_vec = pi_vec)
#     v_r
#     v_mat
#
#     expect_equal(r_fout, cpp_fout)
#   }
# }
# )

