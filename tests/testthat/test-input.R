test_that("bad input is rejected with a clear message", {
  expect_error(tilt_density("a"), "numeric")
  expect_error(tilt_density(c(1, 2)), "at least")
  expect_error(tilt_density(c(1, 2, Inf, 4)), "finite")
  expect_error(tilt_density(rnorm(20), comparator = "gaussian"), "sinc")
  expect_error(tilt_density(rnorm(20), bw = -1), "positive")
  expect_error(tilt_cv(rnorm(20), 0.5, weights = rep(1, 5)), "same length")
})

test_that("missing values are dropped with a warning", {
  set.seed(16)
  x <- c(rnorm(30), NA)
  expect_warning(fit <- tilt_density(x, m = 3), "Missing")
  expect_equal(fit$n, 30L)
})
