## Regression tests for the issues raised in the JSS editorial review.

test_that("multivariate data frames with large-scale columns work (faithful)", {
  # bw_flattop() used to size its grid by 50/scale; for waiting (scale ~ 14)
  # the grid was shorter than the run it searched for and seq_len() failed.
  fit <- tilt_density(datasets::faithful)
  expect_s3_class(fit, "tiltdens_md")
  expect_equal(fit$d, 2L)
  expect_true(all(fit$weights >= -1e-9))
  expect_true(is.finite(bw_flattop(datasets::faithful$waiting)))
  expect_true(is.finite(bw_flattop(datasets::faithful$waiting * 1000)))
  expect_s3_class(tilt_density_cv(datasets::faithful), "tiltdens_md")
})

test_that("sharpen_density refuses multivariate input with a clear message", {
  expect_error(sharpen_density(datasets::faithful), "univariate")
})

test_that("the multistart optimizer matches the genetic algorithm", {
  set.seed(1)
  x <- c(rnorm(40, -1.5), rnorm(40, 1.5))
  ms <- sharpen_density(x, m = 3, optimizer = "multistart")
  ga <- sharpen_density(x, m = 3, optimizer = "ga",
                        control = list(max_iterations = 60, population = 30))
  expect_equal(ms$distance2, ga$distance2, tolerance = 1e-3)
  expect_equal(ms$n_groups, 3L)
})

test_that("a user-supplied optimizer is honoured and bounds are respected", {
  set.seed(2)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  calls <- 0L
  probe <- function(fn, lower, upper) {
    calls <<- calls + 1L
    expect_true(all(lower < 0) && all(upper > 0) && length(lower) == 3L)
    list(par = rep(0, length(lower)), value = fn(rep(0, length(lower))))
  }
  fit <- sharpen_density(x, m = 3, optimizer = probe, polish = FALSE)
  expect_equal(calls, 1L)
  expect_equal(fit$shifts, rep(0, length(x)))
})

test_that("polish uses a bounded method and cannot leave the box", {
  set.seed(3)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  fit <- sharpen_density(x, m = 3, max_shift = 0.05)
  span <- 0.05 * diff(range(x))
  expect_true(all(abs(fit$shifts) <= span + 1e-8))
})

test_that("make_kernel reproduces known constants", {
  logistic <- make_kernel(function(u) exp(-u) / (1 + exp(-u))^2, support = Inf)
  expect_equal(logistic$moment2, pi^2 / 3, tolerance = 1e-4)
  expect_equal(logistic$roughness, 1 / 6, tolerance = 1e-4)

  g <- make_kernel(dnorm, support = Inf)
  expect_equal(g$canonical, bw_canonical_factor("gaussian"), tolerance = 1e-4)
  expect_equal(g$ft(c(0.5, 2)), exp(-c(0.5, 2)^2 / 2), tolerance = 1e-5)
  expect_equal(g$self_conv(0.7), dnorm(0.7, sd = sqrt(2)), tolerance = 1e-5)

  epa <- make_kernel(function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0), support = 1)
  expect_equal(epa$canonical, bw_canonical_factor("epanechnikov"), tolerance = 1e-4)
})

test_that("check_kernel passes good kernels and flags bad ones", {
  for (name in c("gaussian", "epanechnikov", "laplace2")) {
    r <- suppressMessages(capture.output(out <- check_kernel(name)))
    expect_true(out$integrates_to_one && out$non_negative && out$symmetric,
                label = name)
  }
  bad <- make_kernel(function(u) ifelse(abs(u) <= 1, 1 - u^2, 0), support = 1)
  invisible(capture.output(out <- check_kernel(bad)))
  expect_false(out$integrates_to_one)
  asym <- make_kernel(function(u) ifelse(u >= 0 & u <= 1, 2 * (1 - u), 0), support = 1)
  invisible(capture.output(out <- check_kernel(asym)))
  expect_false(out$symmetric)
})

test_that("a custom kernel from make_kernel works end to end", {
  set.seed(4)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  logistic <- make_kernel(function(u) exp(-u) / (1 + exp(-u))^2, support = Inf)
  fit <- tilt_density(x, m = 3, kernel = logistic)
  expect_true(all(fit$y >= 0))
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
})

test_that("conventional_density is a proper density and matches uniform tilting", {
  set.seed(5)
  x <- rnorm(50)
  cd <- conventional_density(x, bw = 0.4, from = -6, to = 6, n = 2001)
  mass <- sum(diff(cd$x) * (cd$y[-1] + cd$y[-length(cd$y)]) / 2)
  expect_equal(mass, 1, tolerance = 1e-5)
  # identical to the tilted estimator with uniform weights
  expect_equal(cd$y, tiltdens:::new_tiltdens(x, 0.4, rep(1/50, 50), 2001, -6, 6,
                                             "x", NULL, "t")$y, tolerance = 1e-12)
})

test_that("the tilt_kernel dispatcher returns kernel objects unchanged", {
  k <- tilt_kernel("biweight")
  expect_identical(tilt_kernel(k), k)
})
