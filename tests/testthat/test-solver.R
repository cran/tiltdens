test_that("the simplex QP solves a problem with a known answer", {
  # Minimising ||q - target||^2 over the simplex returns the projection.
  target <- c(0.7, 0.5, -0.2)
  fit <- tiltdens:::solve_simplex_qp(diag(3), target)
  expect_equal(fit$q, c(0.6, 0.4, 0), tolerance = 1e-6)
  expect_equal(sum(fit$q), 1, tolerance = 1e-10)
})

test_that("the exact small solver agrees with the iterative one", {
  set.seed(40)
  for (m in 2:5) {
    M <- matrix(rnorm(m * m), m)
    H <- crossprod(M)
    f <- rnorm(m) / m
    s <- sample(1:5, m, TRUE)

    exact <- tiltdens:::solve_simplex_qp_small(H, f, s)
    iter  <- tiltdens:::fista_simplex(H / outer(s, s), f / s, s / sum(s),
                                      50000L, 1e-14) / s
    iter_value <- as.numeric(t(iter) %*% H %*% iter - 2 * sum(f * iter))

    expect_equal(exact$value, iter_value, tolerance = 1e-8)
    expect_equal(sum(s * exact$q), 1, tolerance = 1e-10)
    expect_true(all(exact$q >= -1e-12))
  }
})

test_that("the fallback solver agrees with quadprog", {
  set.seed(5)
  A <- tiltdens:::self_conv_gram(rnorm(15), 0.5)
  f <- runif(15) / 15

  via_qp    <- tiltdens:::solve_simplex_qp(A, f)
  u0        <- rep(1 / 15, 15)
  via_fista <- tiltdens:::fista_simplex(A, f, u0, 20000L, 1e-14)

  expect_equal(sum(via_fista), 1, tolerance = 1e-8)
  expect_equal(as.numeric(t(via_fista) %*% A %*% via_fista - 2 * sum(f * via_fista)),
               via_qp$value, tolerance = 1e-6)
})

test_that("simplex projection returns a probability vector", {
  set.seed(6)
  for (i in 1:20) {
    p <- tiltdens:::project_simplex(rnorm(10))
    expect_equal(sum(p), 1, tolerance = 1e-12)
    expect_true(all(p >= 0))
  }
})
