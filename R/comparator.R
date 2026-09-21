#' Infinite-order kernel density estimators
#'
#' Evaluate the two comparator estimators used by the tilting and data
#' sharpening methods. Both are built from kernels of unlimited order: they
#' converge very fast when the underlying density is smooth, but they can take
#' negative values and they oscillate in the tails. Removing those defects
#' without giving up the convergence rate is the point of [tilt_density()].
#'
#' `sinc_density()` uses \eqn{L(u) = \sin(u)/(\pi u)}, whose Fourier transform
#' is flat on \eqn{[-1, 1]} and zero outside.
#'
#' `trapezoid_density()` uses \eqn{L(u) = (\cos u - \cos 2u)/(\pi u^2)}, whose
#' Fourier transform is 1 on \eqn{[-1, 1]} and tapers linearly to 0 at
#' \eqn{|z| = 2} (Politis 2003; Politis and Romano 1999).
#'
#' @param x Numeric vector of observations.
#' @param bw Positive bandwidth. Defaults to [bw_comparator_cv()] for the
#'   matching kernel.
#' @param n Number of grid points at which to evaluate.
#' @param from,to Range of the grid. Defaults to the data range extended by
#'   three bandwidths.
#'
#' @return An object of class `"density"`, so that `plot()`, `lines()` and
#'   `print()` work as they do for [stats::density()].
#'
#' @references
#' Politis, D. N. (2003). Adaptive bandwidth choice.
#' *Journal of Nonparametric Statistics* **15**, 517-533.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
#' plot(sinc_density(x))
#' abline(h = 0, col = "grey")   # the estimate dips below zero
#'
#' @name comparator_density
NULL

#' @rdname comparator_density
#' @export
sinc_density <- function(x, bw = NULL, n = 512, from = NULL, to = NULL) {
  comparator_density_impl(x, bw, n, from, to, "sinc", deparse(substitute(x)))
}

#' @rdname comparator_density
#' @export
trapezoid_density <- function(x, bw = NULL, n = 512, from = NULL, to = NULL) {
  comparator_density_impl(x, bw, n, from, to, "trapezoid", deparse(substitute(x)))
}

comparator_density_impl <- function(x, bw, n, from, to, comparator, data_name) {
  x <- check_sample(x)
  if (is.null(bw)) bw <- bw_comparator_cv(x, comparator = comparator)
  check_bw(bw)

  if (is.null(from)) from <- min(x) - 3 * bw
  if (is.null(to))   to   <- max(x) + 3 * bw
  grid <- seq(from, to, length.out = n)

  structure(
    list(
      x = grid,
      y = comparator_values(x, bw, grid, comparator),
      bw = bw,
      n = length(x),
      call = match.call(),
      data.name = data_name,
      has.na = FALSE,
      comparator = comparator
    ),
    class = "density"
  )
}

#' Values of a comparator estimator on a grid
#'
#' @keywords internal
#' @noRd
comparator_values <- function(x, h, grid, comparator) {
  d <- outer(as.vector(grid), as.vector(x), "-")
  zero <- d == 0
  d[zero] <- 1  # placeholder; overwritten below

  vals <- switch(
    match_comparator(comparator),
    sinc = {
      v <- sin(d / h) / (pi * d)
      v[zero] <- 1 / (pi * h)
      v
    },
    trapezoid = {
      v <- h * (cos(d / h) - cos(2 * d / h)) / (pi * d^2)
      v[zero] <- 3 / (2 * pi * h)
      v
    }
  )
  as.vector(rowMeans(vals))
}

#' Cross term between a Gaussian estimator and a comparator
#'
#' Returns the vector `b` with
#' \eqn{b_i = \sum_j \int K_h(x - y_i) L_{h_c}(x - x_j)\,dx}.
#'
#' With `fcheck` the comparator estimator and
#' \eqn{\hat f(x) = \sum_i p_i K_h(x - y_i)}, the squared distance minimised by
#' the 2016 method is
#' \eqn{\|\hat f - \check f\|^2 = p^\top A p - (2/n) p^\top b + \mathrm{const}},
#' where the constant involves neither `p` nor `h`.
#'
#' @keywords internal
#' @noRd
comparator_cross <- function(y, x, h, hc, comparator, kernel = NULL) {
  rowSums(comparator_cross_matrix(y, x, h, hc, comparator, kernel))
}

#' Matrix of pairwise cross terms
#'
#' The entry \eqn{(i,j)} is \eqn{\int K_h(x - y_i) L_{h_c}(x - x_j)\,dx}. The
#' multivariate case needs the individual entries rather than their row sums,
#' because a product kernel makes the integral factorise across dimensions: the
#' l-variate cross term is the elementwise product of the l univariate ones.
#'
#' @keywords internal
#' @noRd
comparator_cross_matrix <- function(y, x, h, hc, comparator, kernel = NULL) {
  w <- h / hc
  t <- outer(as.vector(y), as.vector(x), "-") / hc

  ## phi_K enters through the transfer function of the estimator's kernel.
  phi <- if (is.null(kernel)) function(z) exp(-0.5 * z^2) else kernel$ft
  base <- function(z) phi(w * z)
  g <- cos_integral(t, base, 0, 1)

  if (match_comparator(comparator) == "trapezoid") {
    g <- g + cos_integral(t, function(z) phi(w * z) * (2 - z), 1, 2)
  }

  g / (pi * hc)
}

#' Squared L2 norm of a comparator estimator
#'
#' @return A list with the scalar `value` and the matrix `self_conv` whose
#'   entries are \eqn{\int L_{h_c}(x - x_i) L_{h_c}(x - x_j)\,dx}.
#'
#' @keywords internal
#' @noRd
comparator_normsq <- function(x, hc, comparator) {
  n <- length(x)
  t <- outer(as.vector(x), as.vector(x), "-") / hc

  self_conv <- switch(
    match_comparator(comparator),
    ## The sinc kernel is idempotent under convolution, since phi^2 = phi.
    sinc = cos_integral(t, function(z) rep(1, length(z)), 0, 1) / (pi * hc),
    trapezoid = (cos_integral(t, function(z) rep(1, length(z)), 0, 1) +
                 cos_integral(t, function(z) (2 - z)^2, 1, 2)) / (pi * hc)
  )

  list(value = sum(self_conv) / n^2, self_conv = self_conv)
}

#' Closed-form self-convolution of the trapezoidal kernel
#'
#' Returns \eqn{(1/\pi)\int_0^2 \phi(z)^2 \cos(tz)\,dz} for the trapezoidal
#' transfer function. The caller divides by the bandwidth. Having this in closed
#' form keeps quadrature out of the bandwidth search loop.
#'
#' @keywords internal
#' @noRd
trapezoid_self_conv <- function(t) {
  small <- abs(t) < 1e-6
  tt <- t
  tt[small] <- 1e-6

  s <- sin(tt)
  cc <- cos(tt)

  ## int_0^1 u^2 cos(t u) du and int_0^1 u^2 sin(t u) du
  a <-  s / tt + 2 * cc / tt^2 - 2 * s / tt^3
  b <- -cc / tt + 2 * s / tt^2 + 2 * cc / tt^3 - 2 / tt^3

  value <- s / tt + cos(2 * tt) * a + sin(2 * tt) * b
  value[small] <- 4 / 3            # limit as t -> 0: 1 + 1/3
  value / pi
}
