test_that("sharpening never does worse than leaving the data alone", {
  set.seed(10)
  x  <- c(rnorm(15, -1.5), rnorm(15, 1.5))
  hc <- as.numeric(bw_comparator_cv(x))

  fit <- sharpen_density(x, m = 3, comparator_bw = hc,
                         control = list(max_iterations = 15, population = 15))

  const <- tiltdens:::comparator_normsq(x, hc, "sinc")$value
  no_shift <- tiltdens:::sharpen_objective(rep(0, length(x)), x, fit$bw, hc,
                                           "sinc", const)

  expect_lte(fit$distance2, no_shift + 1e-10)
  expect_equal(fit$sharpened, x + fit$shifts, tolerance = 1e-12)
})

test_that("the sharpened estimate stays non-negative", {
  set.seed(11)
  x <- c(rnorm(15, -1.5), rnorm(15, 1.5))
  fit <- sharpen_density(x, m = 3,
                         control = list(max_iterations = 10, population = 10))
  expect_true(all(fit$y >= 0))
})

test_that("the genetic algorithm minimises a simple function", {
  set.seed(12)
  fit <- real_ga(function(v) sum((v - c(1, -2))^2), c(-5, -5), c(5, 5),
                 max_iterations = 60L, population = 30L)
  expect_equal(fit$par, c(1, -2), tolerance = 0.15)
  expect_lt(fit$value, 0.05)
})
