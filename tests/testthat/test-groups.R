test_that("blocks partition the sample and are contiguous in rank order", {
  set.seed(4)
  x <- rnorm(30)
  blocks <- tiltdens:::blocks_from_breaks(x, c(10, 20))

  expect_equal(sum(blocks$G), 30)
  expect_equal(as.numeric(blocks$group_size), c(10, 10, 10))
  expect_equal(length(unique(blocks$group_of)), 3L)

  labels <- blocks$group_of[order(x)]
  expect_equal(sum(diff(labels) != 0), 2L)
  expect_true(all(diff(labels) >= 0))
})

test_that("equal breaks give near-equal blocks", {
  expect_equal(tiltdens:::equal_breaks(30, 3), c(10L, 20L))
  expect_equal(tiltdens:::equal_breaks(30, 2), 15L)
  expect_equal(tiltdens:::equal_breaks(100, 1), integer(0))
})

test_that("the minimum block size is respected", {
  br <- tiltdens:::enforce_min_block(c(1, 2, 3), n = 100, min_block = 10, m = 4)
  edges <- c(0, br, 100)
  expect_true(all(diff(edges) >= 10))
})

test_that("breakpoint search beats equal spacing on its own criterion", {
  set.seed(30)
  x <- c(rnorm(30, -1.5), rnorm(30, 1.5))
  h <- 0.5
  A <- tiltdens:::self_conv_gram(x, h)
  b <- tiltdens:::comparator_cross(x, x, h, h, "sinc")
  min_block <- tiltdens:::default_min_block(length(x), 3)

  optimal <- tiltdens:::fit_blocked_weights(A, b, x, 3, "optimal", min_block)
  equal   <- tiltdens:::fit_blocked_weights(A, b, x, 3, "equal", min_block)

  expect_lte(optimal$value, equal$value + 1e-12)
  expect_length(optimal$breaks, 2L)
})

test_that("modal breaks land near the trough of a bimodal density", {
  set.seed(31)
  x <- c(rnorm(60, -3), rnorm(60, 3))
  br <- tiltdens:::modal_breaks(x, 2, min_block = 5)
  expect_length(br, 1L)
  expect_lt(abs(sort(x)[br]), 1.5)     # the trough sits near zero
})

test_that("m = Inf gives one weight per observation", {
  set.seed(4)
  x <- rnorm(30)
  fit <- tilt_density(x, m = Inf)
  expect_equal(fit$n_groups, 30L)
  expect_length(fit$breaks, 0L)
})

test_that("user-supplied breakpoints are honoured", {
  set.seed(32)
  x <- rnorm(40)
  fit <- tilt_density(x, m = 3, breaks = c(12, 25))
  expect_equal(fit$breaks, c(12L, 25L))
  expect_equal(fit$breaks_method, "user")
  expect_lte(length(unique(round(fit$weights, 8))), 3L)
})
