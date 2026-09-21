#' Bandwidth selectors
#'
#' Three rules, for three different jobs.
#'
#' `bw_comparator_cv()` chooses the bandwidth of an infinite-order comparator
#' estimator by minimising least-squares cross-validation,
#' \deqn{CV(h) = \int \check f(x|h)^2 dx - \frac{2}{n}\sum_i \check f_{-i}(x_i|h),}
#' equation (2.7) of Doosti and Hall (2016). It is the default throughout the
#' package, and it is the rule the bandwidths used for the published simulations
#' were computed with.
#'
#' `bw_flattop()` implements the empirical rule of Politis (2003): take
#' \eqn{\hat q} to be the smallest frequency beyond which the modulus of the
#' empirical characteristic function stays below \eqn{c\sqrt{\log_{10}(n)/n}}
#' over a stretch of length `kn`, and set \eqn{h = 1/\hat q}. It is also used
#' internally to cap the cross-validation search.
#'
#' `bw_nrd_robust()` is a rule of thumb for the *conventional* estimator:
#' \eqn{(4/3)^{1/5} n^{-1/5}} times a robust scale, rescaled for the kernel in
#' use (see below).
#'
#' `bw_canonical_factor()` and `bw_convert()` handle the fact that bandwidths
#' mean different things for different kernels.
#'
#' @section Why the cross-validation search is bounded:
#' Least-squares cross-validation for an infinite-order kernel is multimodal in
#' `h`, and its deepest minimum is often a spurious one at a very small
#' bandwidth. Past the frequency at which the empirical characteristic function
#' decays into sampling noise, the criterion is chasing noise, so its global
#' minimum badly undersmooths. The search is therefore capped at `q_max_factor`
#' times the frequency identified by `bw_flattop()`.
#'
#' @section Canonical bandwidths:
#' A bandwidth is only meaningful relative to the kernel it scales: `bw = 0.5`
#' smooths far less with the Epanechnikov kernel, supported on \eqn{[-1,1]},
#' than with the Gaussian, which has unit variance. Marron and Nolan (1988)
#' resolve this with the canonical factor
#' \deqn{\delta_K = \left(R(K)/\mu_2(K)^2\right)^{1/5},}
#' where \eqn{R(K) = \int K^2} and \eqn{\mu_2(K) = \int u^2K(u)\,du}.
#' Bandwidths in units of \eqn{\delta_K} are comparable across kernels: the
#' same number produces the same amount of smoothing, and the asymptotically
#' optimal bandwidths coincide.
#'
#' `bw_canonical_factor()` returns \eqn{\delta_K}. `bw_convert()` maps a
#' bandwidth from one kernel to the equivalent for another, multiplying by the
#' ratio of their factors. The estimators use this automatically: when you do
#' not supply `bw`, the bandwidth chosen for the comparator is rescaled to the
#' kernel in use, so switching kernels changes the shape of the fit without
#' changing how much it is smoothed.
#'
#' @param x Numeric vector of observations.
#' @param kernel Kernel the bandwidth is for; see [tilt_kernels()].
#' @param from,to Kernels to convert between, for `bw_convert()`.
#' @param bw Bandwidth to convert.
#' @param comparator `"sinc"` (default) or `"trapezoid"`.
#' @param q_step Spacing of the frequency grid \eqn{q = 1/h} (default 0.05).
#' @param q_max Largest frequency searched. Defaults to `q_max_factor` divided
#'   by `bw_flattop(x)`.
#' @param q_max_factor Multiplier for the automatic cutoff (default 2).
#' @param q_grid Supply the frequency grid outright, overriding the above.
#' @param c_thresh Threshold constant in the flat-top rule (default 2).
#' @param kn Length of the stretch that must stay below threshold. Defaults to
#'   `max(5, sqrt(log10(n)))`.
#' @param t_max Largest frequency searched by the flat-top rule.
#' @param n_grid Number of frequencies searched by the flat-top rule.
#'
#' @return A single positive number. `bw_comparator_cv()` attaches an attribute
#'   `"cv"`, a data frame of the frequency grid and criterion, which is worth
#'   plotting: the criterion is multimodal, and seeing that is more informative
#'   than the single number.
#'
#' @references
#' Doosti, H. and Hall, P. (2016). Making a non-parametric density estimator
#' more attractive, and more accurate, by data perturbation.
#' *Journal of the Royal Statistical Society B* **78**, 445-462.
#'
#' Politis, D. N. (2003). Adaptive bandwidth choice.
#' *Journal of Nonparametric Statistics* **15**, 517-533.
#'
#' Marron, J. S. and Nolan, D. (1988). Canonical kernels for density
#' estimation. *Statistics and Probability Letters* **7**, 195-199.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
#'
#' bw_comparator_cv(x)
#' bw_flattop(x)
#' bw_nrd_robust(x)
#'
#' ## The criterion is multimodal; the search is capped to avoid the
#' ## spurious minima at small bandwidths.
#' cv <- attr(bw_comparator_cv(x), "cv")
#' plot(cv$q, cv$cv, type = "l", xlab = "frequency 1/h", ylab = "CV")
#'
#' ## Bandwidths are not comparable across kernels until rescaled.
#' bw_canonical_factor("gaussian")
#' bw_canonical_factor("epanechnikov")
#' bw_convert(0.5, from = "gaussian", to = "epanechnikov")
#'
#' @name bandwidth
NULL

#' @rdname bandwidth
#' @export
bw_comparator_cv <- function(x, comparator = c("sinc", "trapezoid"),
                             q_step = 0.05, q_max = NULL, q_max_factor = 2,
                             q_grid = NULL) {
  x <- check_sample(x)
  comparator <- match_comparator(match.arg(comparator))
  n <- length(x)

  if (is.null(q_grid)) {
    if (is.null(q_max)) q_max <- q_max_factor / bw_flattop(x)
    q_grid <- seq(q_step, max(q_max, 10 * q_step), by = q_step)
  }
  h_grid <- 1 / q_grid

  ## Both terms depend on the data only through sums over pairwise differences.
  ## Tied observations give a zero difference, at which both kernels have a
  ## finite limit rather than the 0/0 the formulas below would produce; those
  ## pairs are counted separately at the limiting value.
  d     <- outer(x, x, "-")
  d_off <- d[lower.tri(d) | upper.tri(d)]
  tied  <- d_off == 0
  n_tied <- sum(tied)
  d_off <- d_off[!tied]

  cv    <- numeric(length(h_grid))
  chunk <- max(1L, floor(4e6 / max(1L, length(d_off))))

  for (i0 in seq(1L, length(h_grid), by = chunk)) {
    idx <- i0:min(i0 + chunk - 1L, length(h_grid))
    hh  <- h_grid[idx]
    tt  <- outer(d_off, hh, "/")

    if (comparator == "sinc") {
      kern <- sin(tt) / (pi * d_off)
      sum_kern <- colSums(kern) + n_tied / (pi * hh)
      ## The sinc kernel is its own self-convolution.
      sum_self <- sum_kern + n / (pi * hh)
    } else {
      kern <- sweep(cos(tt) - cos(2 * tt), 2L, hh, "*") / (pi * d_off^2)
      sum_kern <- colSums(kern) + n_tied * 3 / (2 * pi * hh)
      sum_self <- colSums(sweep(trapezoid_self_conv(tt), 2L, hh, "/")) +
        (n + n_tied) * (4 / 3) / (pi * hh)
    }

    cv[idx] <- sum_self / n^2 - 2 * sum_kern / (n * (n - 1))
  }

  best <- which.min(cv)
  structure(
    h_grid[best],
    cv = data.frame(q = q_grid, h = h_grid, cv = cv),
    comparator = comparator
  )
}

#' @rdname bandwidth
#' @export
bw_flattop <- function(x, c_thresh = 2, kn = NULL, t_max = NULL, n_grid = 4000) {
  x <- check_sample(x)
  n <- length(x)
  scale <- robust_scale(x)

  if (is.null(kn))    kn    <- max(5, sqrt(log10(n)))
  if (is.null(t_max)) t_max <- 50 / scale

  t <- seq(0, t_max, length.out = n_grid)
  phi_abs <- abs(as.vector(exp(1i * outer(t, x)) %*% rep(1 / n, n)))

  below <- phi_abs < c_thresh * sqrt(log10(n) / n)
  dt      <- t[2L] - t[1L]
  run_len <- max(1L, ceiling(kn / dt))

  q_hat <- NA_real_
  for (k in seq_len(length(t) - run_len)) {
    if (all(below[k:(k + run_len)])) {
      q_hat <- t[k]
      break
    }
  }

  if (is.na(q_hat) || q_hat <= 0) bw_nrd_robust(x) else 1 / q_hat
}

#' @rdname bandwidth
#' @export
bw_nrd_robust <- function(x, kernel = "gaussian") {
  x <- check_sample(x)
  gaussian_rule <- (4 / 3)^(1 / 5) * length(x)^(-1 / 5) * robust_scale(x)
  gaussian_rule * canonical_ratio(kernel, "gaussian")
}

#' @rdname bandwidth
#' @export
bw_canonical_factor <- function(kernel = "gaussian") {
  tilt_kernel(kernel)$canonical
}

#' @rdname bandwidth
#' @export
bw_convert <- function(bw, from, to) {
  if (!is.numeric(bw) || any(!is.finite(bw)) || any(bw <= 0)) {
    stop("`bw` must be positive and finite.", call. = FALSE)
  }
  bw * canonical_ratio(to, from)
}

#' Ratio of canonical factors, for rescaling a bandwidth between kernels
#' @keywords internal
#' @noRd
canonical_ratio <- function(to, from) {
  tilt_kernel(to)$canonical / tilt_kernel(from)$canonical
}
