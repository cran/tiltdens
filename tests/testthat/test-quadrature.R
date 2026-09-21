test_that("Gauss-Legendre integrates polynomials exactly", {
  rule <- tiltdens:::gauss_legendre(0, 3, 2L, 8L)
  # An 8-point rule per panel is exact up to degree 15.
  expect_equal(sum(rule$weights * rule$nodes^7), 3^8 / 8, tolerance = 1e-8)
  expect_equal(sum(rule$weights), 3, tolerance = 1e-12)
})

test_that("the cosine integral matches its closed form at high frequency", {
  # int_0^1 cos(t z) dz = sin(t) / t.  Checked well past the point where a
  # fixed coarse rule loses all accuracy.
  t <- c(0.5, 5, 25, 60)
  got <- tiltdens:::cos_integral(t, function(z) rep(1, length(z)), 0, 1)
  expect_equal(got, sin(t) / t, tolerance = 1e-9)
})

test_that("the cosine integral preserves matrix shape", {
  t <- matrix(seq(0.1, 1.2, length.out = 6), nrow = 2)
  got <- tiltdens:::cos_integral(t, function(z) rep(1, length(z)), 0, 1)
  expect_true(is.matrix(got))
  expect_equal(dim(got), dim(t))
})
