context("Plotting")

test_that("plot_llike works", {
  set.seed(634)
  llike_vec <- sort(stats::rnorm(100))
  itermax <- 500
  plot_llike(llike_vec = llike_vec, itermax = itermax)
}
)
