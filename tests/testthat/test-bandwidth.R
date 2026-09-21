test_that("bandwidth selectors return single positive numbers", {
  set.seed(17)
  x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
  for (f in list(bw_nrd_robust, bw_flattop, bw_comparator_cv)) {
    bw <- as.numeric(f(x))
    expect_length(bw, 1L)
    expect_gt(bw, 0)
    expect_true(is.finite(bw))
  }
})

test_that("the cross-validation criterion is attached and lies on a 1/h grid", {
  set.seed(18)
  x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
  bw <- bw_comparator_cv(x, q_step = 0.05)
  cv <- attr(bw, "cv")

  expect_s3_class(cv, "data.frame")
  expect_equal(as.numeric(bw), cv$h[which.min(cv$cv)])
  # `%%` binds tighter than `/`, hence the parentheses.
  expect_equal((1 / as.numeric(bw)) %% 0.05, 0, tolerance = 1e-8)
})

test_that("the criterion matches a direct evaluation of its definition", {
  set.seed(19)
  x <- rnorm(20)
  n <- length(x)

  for (comparator in c("sinc", "trapezoid")) {
    bw <- bw_comparator_cv(x, comparator = comparator,
                           q_grid = c(1.5, 2.5))
    cv <- attr(bw, "cv")

    for (k in 1:2) {
      h <- cv$h[k]
      d <- outer(x, x, "-")
      diag(d) <- 1
      M <- if (comparator == "sinc") {
        v <- sin(d / h) / (pi * d); diag(v) <- 0; v
      } else {
        v <- h * (cos(d / h) - cos(2 * d / h)) / (pi * d^2); diag(v) <- 0; v
      }
      first <- tiltdens:::comparator_normsq(x, h, comparator)$value
      expect_equal(cv$cv[k], first - 2 * sum(M) / (n * (n - 1)),
                   tolerance = 1e-8, label = paste(comparator, k))
    }
  }
})

test_that("tied observations do not break the cross-validation bandwidth", {
  # Real data are often recorded to limited precision. faithful$eruptions has
  # 146 ties among 272 values, which used to give 0/0 in the criterion.
  x <- datasets::faithful$eruptions
  bw <- bw_comparator_cv(x)
  expect_true(is.finite(bw) && bw > 0)
  expect_true(all(is.finite(attr(bw, "cv")$cv)))

  # Adding an exact duplicate must give a finite, nearby answer, not NaN.
  set.seed(90)
  y <- rnorm(40)
  b1 <- as.numeric(bw_comparator_cv(y))
  b2 <- as.numeric(bw_comparator_cv(c(y, y[1])))
  expect_true(is.finite(b2))
  expect_lt(abs(log(b2 / b1)), 0.5)

  for (comparator in c("sinc", "trapezoid")) {
    expect_true(is.finite(bw_comparator_cv(c(1, 1, 2, 2, 3, 3, 4, 5, 6),
                                           comparator = comparator)))
  }
})
