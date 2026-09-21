#' Scaled Gaussian kernel matrix
#'
#' @param u,v Numeric vectors.
#' @param h Positive bandwidth.
#'
#' @return A `length(u)` by `length(v)` matrix with entries
#'   \eqn{K_h(u_i - v_j)}, where
#'   \eqn{K_h(z) = \exp(-z^2/(2h^2))/(h\sqrt{2\pi})}.
#'
#' @keywords internal
#' @noRd
kernel_matrix <- function(u, v, h, kernel = NULL) {
  d <- outer(as.vector(u), as.vector(v), "-")
  if (is.null(kernel)) return(exp(-0.5 * (d / h)^2) / (h * sqrt(2 * pi)))
  matrix(kernel$dens(as.vector(d / h)) / h, nrow = nrow(d))
}

#' Gram matrix of the squared L2 norm of a Gaussian kernel estimator
#'
#' Returns the matrix `A` with
#' \eqn{A_{ij} = \int K_h(x - y_i) K_h(x - y_j)\,dx = K_{\sqrt{2}h}(y_i - y_j)},
#' using the fact that the convolution of two Gaussian kernels of bandwidth `h`
#' is a Gaussian kernel of bandwidth `sqrt(2) * h`. Hence for weights `p`,
#' \eqn{\int (\sum_i p_i K_h(x - y_i))^2 dx = p^\top A p}.
#'
#' `A` is symmetric positive definite, which is what makes weight selection a
#' convex quadratic program rather than a general nonlinear search.
#'
#' @param y Numeric vector of (possibly shifted) observations.
#' @param h Positive bandwidth.
#'
#' @return A symmetric `length(y)` by `length(y)` matrix.
#'
#' @keywords internal
#' @noRd
self_conv_gram <- function(y, h, kernel = NULL) {
  d <- outer(as.vector(y), as.vector(y), "-")
  a <- if (is.null(kernel)) {
    exp(-0.25 * (d / h)^2) / (2 * h * sqrt(pi))
  } else {
    matrix(kernel$self_conv(as.vector(d / h)) / h, nrow = nrow(d))
  }
  (a + t(a)) / 2
}
