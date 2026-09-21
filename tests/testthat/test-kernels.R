test_that("the Gram matrix equals the squared L2 norm of the estimator", {
  set.seed(1)
  x <- rnorm(12)
  h <- 0.4
  p <- runif(12); p <- p / sum(p)

  quadratic <- as.numeric(t(p) %*% tiltdens:::self_conv_gram(x, h) %*% p)

  grid <- seq(-12, 12, length.out = 40001)
  y <- as.vector(tiltdens:::kernel_matrix(grid, x, h) %*% p)
  numeric_val <- sum(diff(grid) * (y[-1]^2 + y[-length(y)]^2) / 2)

  expect_equal(quadratic, numeric_val, tolerance = 1e-6)
})

test_that("the Gram matrix is symmetric and numerically positive semidefinite", {
  # A Gaussian kernel matrix is positive definite in exact arithmetic, but it is
  # severely ill-conditioned: for n = 20 the smallest eigenvalue is already at
  # the level of rounding error. That is why solve_simplex_qp() ridges the
  # quadratic form before handing it to quadprog, and why it keeps a
  # projected-gradient fallback.
  set.seed(2)
  A <- tiltdens:::self_conv_gram(rnorm(20), 0.5)
  expect_equal(A, t(A))

  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), -1e-12 * max(ev))
})
