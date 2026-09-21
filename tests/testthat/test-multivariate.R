make_bivariate <- function(n = 60) {
  cbind(c(rnorm(n, -1.5), rnorm(n, 1.5)), c(rnorm(n, -1), rnorm(n, 1)))
}

test_that("multivariate quantities factorise across dimensions", {
  set.seed(50)
  X <- make_bivariate(15)
  h <- 0.5

  A <- tiltdens:::self_conv_gram_md(X, h)
  expect_equal(A, tiltdens:::self_conv_gram(X[, 1], h) *
                  tiltdens:::self_conv_gram(X[, 2], h))
  expect_equal(A, t(A))

  b <- tiltdens:::comparator_cross_md(X, X, h, 0.6, "sinc")
  expect_equal(b, rowSums(
    tiltdens:::comparator_cross_matrix(X[, 1], X[, 1], h, 0.6, "sinc") *
    tiltdens:::comparator_cross_matrix(X[, 2], X[, 2], h, 0.6, "sinc")))
})

test_that("a one-column matrix reproduces the univariate quantities", {
  set.seed(51)
  x <- rnorm(12)
  h <- 0.4
  expect_equal(tiltdens:::self_conv_gram_md(matrix(x), h),
               tiltdens:::self_conv_gram(x, h))
  expect_equal(tiltdens:::comparator_cross_md(matrix(x), matrix(x), h, 0.5, "sinc"),
               tiltdens:::comparator_cross(x, x, h, 0.5, "sinc"))
})

test_that("the bivariate estimate is a proper density", {
  set.seed(52)
  X <- make_bivariate(40)
  fit <- tilt_density(X)

  expect_s3_class(fit, "tiltdens_md")
  expect_equal(fit$d, 2L)
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  expect_true(all(fit$weights >= -1e-9))

  grid <- as.matrix(expand.grid(seq(-5, 5, length.out = 30),
                                seq(-5, 5, length.out = 30)))
  expect_true(all(predict(fit, grid) >= 0))
})

test_that("the bivariate comparator goes negative where the tilted fit cannot", {
  set.seed(53)
  X <- make_bivariate(40)
  fit <- tilt_density(X)
  grid <- as.matrix(expand.grid(seq(-5, 5, length.out = 30),
                                seq(-5, 5, length.out = 30)))

  comparator <- tiltdens:::comparator_values_md(X, fit$comparator_bw, grid, "sinc")
  expect_lt(min(comparator), 0)
  expect_gte(min(predict(fit, grid)), 0)
})

test_that("the multivariate CV criterion matches a direct computation", {
  set.seed(54)
  X <- make_bivariate(10)
  h <- 0.5
  n <- nrow(X)
  p <- rep(1 / n, n)

  M <- tiltdens:::kernel_matrix_md(X, X, h)
  diag(M) <- 0
  A <- tiltdens:::self_conv_gram_md(X, h)
  direct <- as.numeric(t(p) %*% A %*% p) - 2 * sum(colSums(M) * p) / n

  expect_equal(as.numeric(tiltdens:::tilt_cv_md(X, h, p)), direct, tolerance = 1e-12)
})

test_that("three dimensions and the CV front end both work", {
  set.seed(55)
  X <- cbind(make_bivariate(25), rnorm(50))
  expect_equal(tilt_density(X)$d, 3L)
  expect_s3_class(tilt_density_cv(make_bivariate(25)), "tiltdens_md")
})

test_that("unsupported multivariate options are refused clearly", {
  set.seed(56)
  X <- make_bivariate(20)
  expect_error(tilt_density(X, constraint = "unimodal"), "univariate")
  expect_warning(tilt_density(X, m = 3), "ignored")
  expect_error(plot(tilt_density(cbind(X, rnorm(40)))), "two-dimensional")
})
