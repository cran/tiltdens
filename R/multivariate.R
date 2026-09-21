#' Multivariate building blocks
#'
#' Section 3.3 of Doosti and Hall (2016) extends the method to l-variate data by
#' replacing the univariate kernel with a product kernel
#' \eqn{K_{\mathrm{pr}}(x) = \prod_j K(x_j)} and taking a scalar bandwidth,
#' \eqn{H = \mathrm{diag}(h, \dots, h)}.
#'
#' That choice makes the extension almost free. Every quantity the method needs
#' is an integral of a product of two product kernels, so it factorises across
#' dimensions:
#'
#' \deqn{\int K_{\mathrm{pr},h}(x - u) L_{\mathrm{pr},g}(x - v)\,dx
#'       = \prod_{j=1}^{l} \int K_h(x_j - u_j) L_g(x_j - v_j)\,dx_j.}
#'
#' So the l-variate matrix is the elementwise product of the l univariate
#' matrices, each computed by the code that already exists. The weights are
#' still an n-vector on the probability simplex, so the optimisation problem is
#' unchanged: same convex quadratic program, same solver.
#'
#' @param y,x Numeric matrices with one row per observation.
#' @param h,hc Positive scalar bandwidths.
#' @param comparator `"sinc"` or `"trapezoid"`.
#'
#' @name multivariate_internals
#' @keywords internal
#' @noRd
NULL

#' @keywords internal
#' @noRd
self_conv_gram_md <- function(y, h, kernel = NULL) {
  y <- as.matrix(y)
  out <- self_conv_gram(y[, 1L], h, kernel)
  for (j in seq_len(ncol(y))[-1L]) out <- out * self_conv_gram(y[, j], h, kernel)
  (out + t(out)) / 2
}

#' @keywords internal
#' @noRd
comparator_cross_md <- function(y, x, h, hc, comparator, kernel = NULL) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  out <- comparator_cross_matrix(y[, 1L], x[, 1L], h, hc, comparator, kernel)
  for (j in seq_len(ncol(y))[-1L]) {
    out <- out * comparator_cross_matrix(y[, j], x[, j], h, hc, comparator, kernel)
  }
  rowSums(out)
}

#' @keywords internal
#' @noRd
comparator_normsq_md <- function(x, hc, comparator) {
  x <- as.matrix(x)
  out <- comparator_normsq(x[, 1L], hc, comparator)$self_conv
  for (j in seq_len(ncol(x))[-1L]) {
    out <- out * comparator_normsq(x[, j], hc, comparator)$self_conv
  }
  list(value = sum(out) / nrow(x)^2, self_conv = out)
}

#' @keywords internal
#' @noRd
kernel_matrix_md <- function(u, v, h, kernel = NULL) {
  u <- as.matrix(u)
  v <- as.matrix(v)
  out <- kernel_matrix(u[, 1L], v[, 1L], h, kernel)
  for (j in seq_len(ncol(u))[-1L]) out <- out * kernel_matrix(u[, j], v[, j], h, kernel)
  out
}

#' Product-kernel comparator estimator on a set of points
#' @keywords internal
#' @noRd
comparator_values_md <- function(x, h, points, comparator) {
  x <- as.matrix(x)
  points <- as.matrix(points)
  out <- comparator_values_matrix(points[, 1L], x[, 1L], h, comparator)
  for (j in seq_len(ncol(x))[-1L]) {
    out <- out * comparator_values_matrix(points[, j], x[, j], h, comparator)
  }
  rowMeans(out)
}

#' @keywords internal
#' @noRd
comparator_values_matrix <- function(grid, x, h, comparator) {
  d <- outer(as.vector(grid), as.vector(x), "-")
  zero <- d == 0
  d[zero] <- 1

  switch(
    match_comparator(comparator),
    sinc = { v <- sin(d / h) / (pi * d); v[zero] <- 1 / (pi * h); v },
    trapezoid = {
      v <- h * (cos(d / h) - cos(2 * d / h)) / (pi * d^2)
      v[zero] <- 3 / (2 * pi * h)
      v
    }
  )
}

#' Bandwidth for multivariate data
#'
#' Applies the univariate cross-validation rule to each coordinate and takes the
#' geometric mean, then rescales by \eqn{n^{1/5 - 1/(l+4)}} so that the exponent
#' matches the l-variate optimal rate.
#'
#' @keywords internal
#' @noRd
bw_comparator_cv_md <- function(x, comparator) {
  x <- as.matrix(x)
  n <- nrow(x)
  l <- ncol(x)
  per_dim <- vapply(seq_len(l), function(j) {
    as.numeric(bw_comparator_cv(x[, j], comparator = comparator))
  }, numeric(1))
  exp(mean(log(per_dim))) * n^(1 / 5 - 1 / (l + 4))
}

#' Is this input multivariate?
#' @keywords internal
#' @noRd
is_multivariate <- function(x) {
  (is.matrix(x) || is.data.frame(x)) && ncol(as.matrix(x)) > 1L
}

#' Validate a multivariate sample
#' @keywords internal
#' @noRd
check_sample_md <- function(x, min_n = 3L) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  keep <- stats::complete.cases(x)
  if (!all(keep)) {
    x <- x[keep, , drop = FALSE]
    warning("Rows with missing values were removed.", call. = FALSE)
  }
  if (!all(is.finite(x))) stop("`x` must contain only finite values.", call. = FALSE)
  if (nrow(x) < min_n) {
    stop("`x` must have at least ", min_n, " rows.", call. = FALSE)
  }
  if (is.null(colnames(x))) colnames(x) <- paste0("V", seq_len(ncol(x)))
  x
}
