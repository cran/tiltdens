trapz <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
## Wide enough for the heaviest-tailed kernel here: laplace4 has variance 8, so
## a range of a few units truncates its Fourier integral noticeably.
fine  <- seq(-60, 60, length.out = 400001)

test_that("every kernel is a density with the right Fourier transform at zero", {
  for (name in tilt_kernels()) {
    K <- tilt_kernel(name)
    expect_equal(trapz(fine, K$dens(fine)), 1, tolerance = 1e-5, label = name)
    expect_true(all(K$dens(fine) >= 0), label = name)
    expect_equal(K$ft(0), 1, tolerance = 1e-6, label = name)
  }
})

test_that("each Fourier transform matches numerical integration of its kernel", {
  tt <- c(0.5, 1.5, 3)
  for (name in tilt_kernels()) {
    K <- tilt_kernel(name)
    numeric_ft <- vapply(tt, function(t) trapz(fine, cos(t * fine) * K$dens(fine)),
                         numeric(1))
    expect_equal(K$ft(tt), numeric_ft, tolerance = 1e-5, label = name)
  }
})

test_that("each self-convolution matches a direct convolution", {
  for (name in tilt_kernels()) {
    K <- tilt_kernel(name)
    for (v in c(-1.3, -0.3, 0, 0.55, 1.7)) {
      direct <- trapz(fine, K$dens(fine) * K$dens(v - fine))
      expect_equal(K$self_conv(v), direct, tolerance = 1e-6,
                   label = paste(name, "at", v))
    }
  }
})

test_that("the Laplace family matches the examples in the 2016 paper", {
  # Section 3.1 gives K(x) = exp(-|x|)/2 for k = 1 and (1+|x|)exp(-|x|)/4 for k = 2.
  u <- c(-2, -0.4, 0, 0.7, 3)
  expect_equal(tilt_kernel("laplace")$dens(u), exp(-abs(u)) / 2)
  expect_equal(tilt_kernel("laplace2")$dens(u), (1 + abs(u)) * exp(-abs(u)) / 4)
  # The self-convolution of the k = 1 kernel is the k = 2 kernel.
  expect_equal(tilt_kernel("laplace")$self_conv(u),
               tilt_kernel("laplace2")$dens(u), tolerance = 1e-8)
})

test_that("the Gram matrix is positive semidefinite for every kernel", {
  # A has Fourier transform phi^2 >= 0, so this must hold whatever K is.
  set.seed(70)
  x <- rnorm(25)
  for (name in tilt_kernels()) {
    A <- tiltdens:::self_conv_gram(x, 0.6, tilt_kernel(name))
    ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(ev), -1e-8 * max(ev), label = name)
  }
})

test_that("the Gram matrix equals the squared L2 norm for a non-Gaussian kernel", {
  set.seed(71)
  x <- rnorm(10)
  h <- 0.5
  p <- runif(10); p <- p / sum(p)
  K <- tilt_kernel("epanechnikov")

  quadratic <- as.numeric(t(p) %*% tiltdens:::self_conv_gram(x, h, K) %*% p)
  grid <- seq(-8, 8, length.out = 200001)
  y <- as.vector(tiltdens:::kernel_matrix(grid, x, h, K) %*% p)
  expect_equal(quadratic, trapz(grid, y^2), tolerance = 1e-5)
})

test_that("fits with any kernel are proper densities", {
  set.seed(72)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  for (name in c("gaussian", "laplace2", "epanechnikov", "triangular")) {
    fit <- tilt_density(x, m = 3, kernel = name)
    expect_equal(sum(fit$weights), 1, tolerance = 1e-8, label = name)
    expect_true(all(fit$weights >= -1e-9), label = name)
    expect_true(all(fit$y >= 0), label = name)
    expect_equal(fit$kernel$name, name)
  }
})

test_that("kernels reach the CV, sharpening and multivariate front ends", {
  set.seed(73)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  expect_true(all(tilt_density_cv(x, kernel = "laplace2")$y >= 0))
  expect_true(all(sharpen_density(x, m = 3, kernel = "epanechnikov",
                                  control = list(max_iterations = 5,
                                                 population = 10))$y >= 0))
  X <- cbind(x, c(rnorm(30, -1), rnorm(30, 1)))
  fit <- tilt_density(X, kernel = "laplace2")
  expect_equal(fit$kernel$name, "laplace2")
  expect_true(all(predict(fit, X) >= 0))
})

test_that("a custom kernel can be supplied and bad ones are refused", {
  logistic <- list(
    name = "logistic",
    dens = function(u) exp(-u) / (1 + exp(-u))^2,
    ft   = function(t) ifelse(abs(t) < 1e-8, 1, pi * t / sinh(pi * t)),
    self_conv = NULL,
    support = Inf
  )
  expect_error(tilt_kernel(logistic), "self_conv")

  set.seed(74)
  bump <- list(dens = function(u) dnorm(u), ft = function(t) exp(-t^2 / 2),
               self_conv = function(v) dnorm(v, sd = sqrt(2)), support = Inf)
  fit <- tilt_density(rnorm(30), m = 3, kernel = bump)
  expect_equal(fit$kernel$name, "custom")

  expect_error(tilt_kernel(list(dens = dnorm)), "ft")
  expect_error(tilt_kernel(list(dens = dnorm, ft = function(t) rep(2, length(t)))),
               "ft\\(0\\)")
  expect_error(tilt_kernel("cauchy"), "Unknown kernel")
})

test_that("canonical factors match the published Marron-Nolan values", {
  published <- c(gaussian = 0.7764, epanechnikov = 1.7188,
                 biweight = 2.0362, triweight = 2.3122, triangular = 1.8882)
  for (name in names(published)) {
    expect_equal(bw_canonical_factor(name), published[[name]],
                 tolerance = 1e-4, label = name)
  }
})

test_that("roughness and second moment agree with numerical integration", {
  for (name in tilt_kernels()) {
    K <- tilt_kernel(name)
    expect_equal(K$roughness, trapz(fine, K$dens(fine)^2),
                 tolerance = 1e-5, label = name)
    expect_equal(K$moment2, trapz(fine, fine^2 * K$dens(fine)),
                 tolerance = 1e-4, label = name)
    # R(K) is the self-convolution at zero, for a symmetric kernel.
    expect_equal(K$roughness, K$self_conv(0), tolerance = 1e-12, label = name)
  }
})

test_that("bandwidth conversion is consistent and reversible", {
  expect_equal(bw_convert(0.5, "gaussian", "gaussian"), 0.5)
  expect_equal(bw_convert(bw_convert(0.5, "gaussian", "biweight"),
                          "biweight", "gaussian"), 0.5, tolerance = 1e-12)
  expect_gt(bw_convert(0.5, "gaussian", "epanechnikov"), 0.5)  # wider support
  expect_error(bw_convert(-1, "gaussian", "biweight"), "positive")
})

test_that("the rule of thumb is rescaled for the kernel", {
  set.seed(80)
  x <- rnorm(100)
  ratio <- bw_nrd_robust(x, "epanechnikov") / bw_nrd_robust(x, "gaussian")
  expect_equal(ratio, bw_canonical_factor("epanechnikov") /
                      bw_canonical_factor("gaussian"), tolerance = 1e-10)
})

test_that("canonical rescaling makes kernels smooth by comparable amounts", {
  # The point of the rescaling: with it, the choice of kernel changes the shape
  # of the fit but not how much it is smoothed, so accuracy varies far less.
  set.seed(2016)
  x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
  truth <- function(t) 0.5 * dnorm(t, -1.5) + 0.5 * dnorm(t, 1.5)
  kernels <- c("gaussian", "epanechnikov", "biweight", "triangular")

  rescaled <- vapply(kernels, function(k) {
    ise(tilt_density(x, m = 3, kernel = k), truth)
  }, numeric(1))

  h_gauss <- tilt_density(x, m = 3, kernel = "gaussian")$bw
  fixed <- vapply(kernels, function(k) {
    ise(tilt_density(x, m = 3, kernel = k, bw = h_gauss), truth)
  }, numeric(1))

  expect_lt(max(rescaled) / min(rescaled), max(fixed) / min(fixed))
  expect_lt(max(rescaled) / min(rescaled), 1.5)
})

test_that("default bandwidths follow the canonical ratio", {
  set.seed(81)
  x <- c(rnorm(40, -1.5), rnorm(40, 1.5))
  g <- tilt_density(x, m = 3, kernel = "gaussian")
  e <- tilt_density(x, m = 3, kernel = "epanechnikov")
  expect_equal(e$bw / g$bw, bw_canonical_factor("epanechnikov") /
                            bw_canonical_factor("gaussian"), tolerance = 1e-10)
})
