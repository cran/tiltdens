test_that("comparator estimators integrate to one", {
  set.seed(1)
  x <- rnorm(8)
  grid <- seq(-400, 400, length.out = 200001)

  for (comparator in c("sinc", "trapezoid")) {
    y <- tiltdens:::comparator_values(x, 0.5, grid, comparator)
    mass <- sum(diff(grid) * (y[-1] + y[-length(y)]) / 2)
    expect_equal(mass, 1, tolerance = 2e-3, label = comparator)
  }
})

test_that("the sinc estimator takes negative values, as it should", {
  set.seed(3)
  x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
  expect_lt(min(sinc_density(x)$y), 0)
})

test_that("the comparator norm matches numerical integration", {
  set.seed(3)
  x  <- rnorm(6)
  hc <- 0.5
  grid <- seq(-200, 200, length.out = 200001)

  for (comparator in c("sinc", "trapezoid")) {
    analytic <- tiltdens:::comparator_normsq(x, hc, comparator)$value
    y <- tiltdens:::comparator_values(x, hc, grid, comparator)
    numeric_val <- sum(diff(grid) * (y[-1]^2 + y[-length(y)]^2) / 2)
    expect_equal(analytic, numeric_val, tolerance = 1e-3, label = comparator)
  }
})

test_that("the cross vector matches numerical integration", {
  set.seed(2)
  x  <- rnorm(6)
  h  <- 0.35
  hc <- 0.5
  n  <- length(x)
  grid <- seq(-60, 60, length.out = 200001)

  for (comparator in c("sinc", "trapezoid")) {
    b <- tiltdens:::comparator_cross(x, x, h, hc, comparator)
    fcheck <- tiltdens:::comparator_values(x, hc, grid, comparator)

    numeric_val <- vapply(seq_len(n), function(i) {
      k <- dnorm(grid, mean = x[i], sd = h)
      prod_vals <- k * fcheck
      n * sum(diff(grid) * (prod_vals[-1] + prod_vals[-length(prod_vals)]) / 2)
    }, numeric(1))

    expect_equal(b, numeric_val, tolerance = 1e-4, label = comparator)
  }
})

test_that("the closed-form trapezoid self-convolution matches quadrature", {
  set.seed(20)
  x  <- rnorm(7)
  hc <- 0.45
  t  <- outer(x, x, "-") / hc

  closed <- tiltdens:::trapezoid_self_conv(t) / hc
  quad   <- tiltdens:::comparator_normsq(x, hc, "trapezoid")$self_conv

  expect_equal(closed, quad, tolerance = 1e-8)
})
