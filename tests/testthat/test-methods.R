test_that("fits inherit from density so base plotting works", {
  set.seed(13)
  x <- rnorm(40)
  fit <- tilt_density(x, m = 3)

  expect_s3_class(fit, "tiltdens")
  expect_s3_class(fit, "density")
  expect_true(all(c("x", "y", "bw", "n") %in% names(fit)))
  expect_silent(invisible(capture.output(print(fit))))
})

test_that("predict matches the stored grid values", {
  set.seed(14)
  x <- rnorm(40)
  fit <- tilt_density(x, m = 3)
  expect_equal(predict(fit), fit$y)
  expect_equal(predict(fit, fit$x[1:5]), fit$y[1:5], tolerance = 1e-12)
})

test_that("ise is zero when the estimate is compared with itself", {
  set.seed(15)
  fit <- tilt_density(rnorm(40), m = 3)
  f <- stats::approxfun(fit$x, fit$y, rule = 2)
  expect_equal(ise(fit, f), 0, tolerance = 1e-12)
})
