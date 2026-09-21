## Count modes by prominence: how far a local maximum rises above the lower of
## the two troughs flanking it, relative to the global peak. Counting sign
## changes of the derivative instead would report every ripple, including ones
## of relative size 1e-5 left over between constraint grid points, which are an
## artefact of the discretisation rather than a feature of the estimate.
n_modes <- function(fit, x, min_prominence = 0.01) {
  g <- seq(min(x), max(x), length.out = 800)
  y <- predict(fit, g)
  dy <- diff(y)
  peaks <- which(dy[-length(dy)] > 0 & dy[-1] <= 0) + 1L
  if (!length(peaks)) return(0L)

  ## Topographic prominence: walk out from the peak in each direction until the
  ## curve rises above it, and take the higher of the two intervening minima.
  ## Using the global minima instead would credit a ripple sitting on the
  ## shoulder of a real mode with the whole height of that shoulder.
  prominence <- vapply(peaks, function(k) {
    left <- k
    while (left > 1L && y[left - 1L] <= y[k]) left <- left - 1L
    right <- k
    while (right < length(y) && y[right + 1L] <= y[k]) right <- right + 1L
    (y[k] - max(min(y[left:k]), min(y[k:right]))) / max(y)
  }, numeric(1))

  sum(prominence > min_prominence)
}

test_that("the unimodal constraint holds exactly on its grid", {
  set.seed(7)
  x <- c(rnorm(40, 0, 1), rnorm(6, 2.6, 0.12))   # a spurious secondary bump
  fit <- tilt_density(x, m = Inf, constraint = "unimodal")

  grid <- seq(min(x), max(x), length.out = 100)
  y <- predict(fit, grid)
  dy <- diff(y)
  k <- which.min(abs(grid - fit$mode_location))

  expect_gte(min(dy[seq_len(k - 1L)]), -1e-9)     # increasing before the mode
  expect_lte(max(dy[k:length(dy)]), 1e-9)         # decreasing after it
  expect_equal(n_modes(fit, x), 1L)
})

test_that("the constraint removes a mode the unconstrained fit keeps", {
  set.seed(7)
  x <- c(rnorm(40, 0, 1), rnorm(6, 2.6, 0.12))
  expect_gt(n_modes(tilt_density(x, m = Inf), x), 1L)
  expect_equal(n_modes(tilt_density(x, m = Inf, constraint = "unimodal"), x), 1L)
})

test_that("constraining cannot improve the unconstrained criterion", {
  set.seed(60)
  x <- c(rnorm(40, 0, 1), rnorm(6, 2.6, 0.12))
  free <- tilt_density(x, m = Inf)
  cons <- tilt_density(x, m = Inf, constraint = "unimodal")
  expect_gte(cons$distance2, free$distance2 - 1e-10)
})

test_that("monotone constraints are enforced", {
  set.seed(61)
  x <- rexp(100)

  dec <- tilt_density(x, m = Inf, constraint = "decreasing")
  grid <- seq(min(x), max(x), length.out = 100)
  expect_lte(max(diff(predict(dec, grid))), 1e-9)

  inc <- tilt_density(-x, m = Inf, constraint = "increasing")
  grid <- seq(min(-x), max(-x), length.out = 100)
  expect_gte(min(diff(predict(inc, grid))), -1e-9)
})

test_that("constrained fits are still proper densities", {
  set.seed(62)
  x <- c(rnorm(40, 0, 1), rnorm(6, 2.6, 0.12))
  fit <- tilt_density(x, m = Inf, constraint = "unimodal")
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  expect_true(all(fit$weights >= -1e-9))
  expect_true(all(fit$y >= 0))
})

test_that("an infeasible constraint gives an explanatory error", {
  set.seed(63)
  x <- c(rnorm(40, 0, 1), rnorm(6, 2.6, 0.12))
  expect_error(tilt_density_cv(x, max_groups = 2, constraint = "unimodal"),
               "distinct weight")
  expect_error(tilt_density(x, constraint = "banana"), "constraint")
})

test_that("the constrained solver returns a feasible point", {
  set.seed(64)
  x <- rnorm(30)
  h <- 0.5
  n <- length(x)
  A <- tiltdens:::self_conv_gram(x, h)
  b <- tiltdens:::comparator_cross(x, x, h, h, "sinc")
  grid <- seq(min(x), max(x), length.out = 40)
  C <- tiltdens:::constraint_matrix(x, h, "unimodal", grid, mode_index = 20)

  fit <- tiltdens:::solve_constrained_simplex_qp(A, b / n, rep(1, n), C)
  expect_true(fit$feasible)
  expect_gte(min(C %*% fit$q), -1e-9)
  expect_equal(sum(fit$q), 1, tolerance = 1e-9)
})
