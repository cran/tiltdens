test_that("tilt weights are a proper probability distribution", {
  set.seed(5)
  x <- c(rnorm(25, -1.5), rnorm(25, 1.5))

  for (m in list(3, Inf)) {
    fit <- tilt_density(x, m = m)
    expect_true(all(fit$weights >= -1e-9))
    expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  }
})

test_that("m = 3 really uses only three distinct weights", {
  set.seed(5)
  x <- c(rnorm(25, -1.5), rnorm(25, 1.5))
  fit <- tilt_density(x, m = 3)
  expect_lte(length(unique(round(fit$weights, 8))), 3L)
})

test_that("tilting never does worse than uniform weights on its own criterion", {
  set.seed(6)
  x  <- c(rnorm(20, -1.5), rnorm(20, 1.5))
  hc <- as.numeric(bw_comparator_cv(x))
  fit <- tilt_density(x, m = Inf, comparator_bw = hc)

  n <- length(x)
  A <- tiltdens:::self_conv_gram(x, fit$bw)
  b <- tiltdens:::comparator_cross(x, x, fit$bw, hc, "sinc")
  const <- tiltdens:::comparator_normsq(x, hc, "sinc")$value

  uniform <- rep(1 / n, n)
  loss_uniform <- as.numeric(t(uniform) %*% A %*% uniform) -
    2 * sum(b * uniform) / n + const

  expect_lte(fit$distance2, loss_uniform + 1e-10)
  expect_gte(fit$distance2, -1e-8)
})

test_that("the tilted estimate is a proper density", {
  set.seed(7)
  x <- rnorm(30)
  fit <- tilt_density(x, m = 3, from = -15, to = 15, n = 40001)
  expect_true(all(fit$y >= 0))
  mass <- sum(diff(fit$x) * (fit$y[-1] + fit$y[-length(fit$y)]) / 2)
  expect_equal(mass, 1, tolerance = 1e-5)
})

test_that("both comparators are accepted and give sane fits", {
  set.seed(8)
  x <- c(rnorm(25, -1.5), rnorm(25, 1.5))
  for (comparator in c("sinc", "trapezoid")) {
    fit <- tilt_density(x, m = 3, comparator = comparator)
    expect_equal(fit$comparator, comparator)
    expect_true(all(fit$y >= 0))
  }
})

test_that("a bandwidth grid picks the best of the grid", {
  set.seed(9)
  x <- c(rnorm(25, -1.5), rnorm(25, 1.5))
  grid <- c(0.3, 0.5, 0.8)
  fit <- tilt_density(x, m = 3, bw_grid = grid)
  expect_true(fit$bw %in% grid)
  expect_equal(fit$distance2, min(fit$distance_grid))
})
