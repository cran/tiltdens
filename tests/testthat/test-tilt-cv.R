test_that("the CV criterion matches its definition computed directly", {
  set.seed(8)
  x <- rnorm(15)
  h <- 0.4
  p <- runif(15); p <- p / sum(p)
  n <- length(x)

  grid <- seq(-15, 15, length.out = 100001)
  y <- as.vector(tiltdens:::kernel_matrix(grid, x, h) %*% p)
  first <- sum(diff(grid) * (y[-1]^2 + y[-length(y)]^2) / 2)

  second <- sum(vapply(seq_len(n), function(i) {
    keep <- setdiff(seq_len(n), i)
    sum(p[keep] * dnorm((x[i] - x[keep]) / h) / h)
  }, numeric(1)))

  expect_equal(as.numeric(tilt_cv(x, h, p)), first - 2 * second / n,
               tolerance = 1e-6)
})

test_that("CV selection returns a valid fit that improves on uniform weights", {
  set.seed(9)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  fit <- tilt_density_cv(x, max_groups = 3)

  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  expect_true(all(fit$weights >= -1e-9))
  expect_gt(fit$bw, 0)
  expect_equal(fit$cv, as.numeric(tilt_cv(x, fit$bw, fit$weights)),
               tolerance = 1e-8)
  expect_lte(fit$cv, fit$trace$cv[1] + 1e-10)
})

test_that("the trace records one row per step", {
  set.seed(10)
  x <- rnorm(40)
  fit <- tilt_density_cv(x, max_groups = 3)
  expect_equal(nrow(fit$trace), 3L)
  expect_equal(fit$trace$n_groups, 1:3)
})
