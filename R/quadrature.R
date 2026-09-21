#' Composite Gauss-Legendre nodes and weights
#'
#' Splits `[lo, hi]` into `n_panels` equal panels and places an `order`-point
#' Gauss-Legendre rule on each. The reference nodes come from the Golub-Welsch
#' eigenvalue problem, so no external quadrature routine is needed.
#'
#' @param lo,hi Numeric endpoints of the interval.
#' @param n_panels Number of equal panels.
#' @param order Points per panel.
#'
#' @return A list with components `nodes` and `weights`, both numeric vectors of
#'   length `n_panels * order`, such that `sum(weights * g(nodes))` approximates
#'   the integral of `g` over `[lo, hi]`.
#'
#' @keywords internal
#' @noRd
gauss_legendre <- function(lo, hi, n_panels, order = 8L) {
  order    <- max(2L, as.integer(order))
  n_panels <- max(1L, as.integer(n_panels))

  k    <- seq_len(order - 1L)
  beta <- k / sqrt(4 * k^2 - 1)
  jacobi <- matrix(0, order, order)
  jacobi[cbind(k, k + 1L)] <- beta
  jacobi[cbind(k + 1L, k)] <- beta

  eig    <- eigen(jacobi, symmetric = TRUE)
  ord    <- order(eig$values)
  x_ref  <- eig$values[ord]
  w_ref  <- 2 * eig$vectors[1L, ord]^2

  edges    <- seq(lo, hi, length.out = n_panels + 1L)
  half_len <- (edges[2L] - edges[1L]) / 2
  centres  <- (edges[-length(edges)] + edges[-1L]) / 2

  list(
    nodes   = as.vector(outer(half_len * x_ref, centres, "+")),
    weights = rep(half_len * w_ref, times = n_panels)
  )
}

#' Oscillatory cosine integral
#'
#' Evaluates `integral from lo to hi of weight_fun(z) * cos(t * z) dz` for every
#' element of `t`.
#'
#' Every inner product between two kernel estimators in this package reduces to
#' an integral of this form by Parseval's identity: if `K` and `L` have Fourier
#' transforms `phiK` and `phiL`, then
#'
#' \deqn{\int K_h(x-u) L_g(x-v)\,dx =
#'       \frac{1}{\pi g}\int_0^\infty \phi_K(wz)\phi_L(z)\cos(z(u-v)/g)\,dz}
#'
#' with \eqn{w = h/g}. Choosing `weight_fun` and the limits selects the kernel
#' pair.
#'
#' The number of quadrature panels is scaled with `max(abs(t))` so that the
#' integrand stays resolved when the cosine oscillates rapidly, which happens
#' whenever the data range is large relative to the bandwidth. A fixed-node rule
#' silently loses all accuracy in that regime.
#'
#' @param t Numeric vector or matrix of frequencies.
#' @param weight_fun Function of a numeric vector, returning a vector of the
#'   same length.
#' @param lo,hi Integration limits.
#' @param panels_per_cycle Quadrature panels allocated per oscillation.
#'
#' @return A numeric object with the same shape as `t`.
#'
#' @keywords internal
#' @noRd
cos_integral <- function(t, weight_fun, lo, hi, panels_per_cycle = 4) {
  max_t    <- max(abs(t))
  n_cycles <- max_t * (hi - lo) / (2 * pi)
  n_panels <- min(max(8, ceiling(panels_per_cycle * n_cycles)), 20000)

  rule <- gauss_legendre(lo, hi, n_panels, 8L)
  g    <- weight_fun(rule$nodes) * rule$weights

  t_vec  <- as.vector(t)
  result <- numeric(length(t_vec))

  ## Chunked so the (length(t) x length(nodes)) cosine matrix stays bounded.
  chunk <- max(1L, floor(4e6 / length(rule$nodes)))
  starts <- seq(1L, length(t_vec), by = chunk)
  for (i0 in starts) {
    idx <- i0:min(i0 + chunk - 1L, length(t_vec))
    result[idx] <- as.vector(cos(outer(t_vec[idx], rule$nodes)) %*% g)
  }

  if (is.matrix(t)) matrix(result, nrow = nrow(t)) else result
}
