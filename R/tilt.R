#' Tilted density estimation by distance to a comparator
#'
#' Chooses non-negative weights `p` summing to one for the estimator
#' \deqn{\hat f(x \mid h, p) = \sum_i p_i K_h(x - x_i)}
#' by minimising its L2 distance to an infinite-order comparator estimator.
#' This is the method of Doosti and Hall (2016): it inherits the fast
#' convergence rate of the comparator while staying a proper, non-negative
#' density with no oscillatory tails.
#'
#' The problem is a convex quadratic program over the probability simplex, so it
#' has a unique solution.
#'
#' @section Choosing `m`:
#' `m = Inf` lets every observation take its own weight, which is the most
#' flexible option and the one the paper calls \eqn{T_n}. `m = 3` restricts the
#' weights to three distinct values (\eqn{T_3}); it is faster, and in the
#' paper's simulations it was often more accurate, because fewer free parameters
#' means less overfitting.
#'
#' When `m` is finite the block boundaries matter as much as the number of
#' blocks, and `breaks` controls them. The default `"optimal"` treats the
#' boundaries as part of the optimisation, which is what Section 4.1 of the 2016
#' paper specifies. On bimodal data it typically places them at the trough
#' between the modes without being told to, which is the same answer `"modal"`
#' arrives at by construction from a pilot estimate.
#'
#' @param x Numeric vector of observations.
#' @param m Number of distinct weight values. `Inf` (default) or 3 are the cases
#'   studied in the paper.
#' @param comparator `"sinc"` (default) or `"trapezoid"`.
#' @param kernel Kernel of the perturbed estimator: `"gaussian"` (default),
#'   `"laplace"`, `"laplace2"`, `"epanechnikov"`, `"biweight"` and others; see
#'   [tilt_kernels()] for the full list, or supply a custom kernel as a list.
#'   Bandwidths are on each kernel's own scale and are not comparable across
#'   kernels.
#' @param bw Bandwidth of the tilted estimator. Defaults to the comparator
#'   bandwidth, as in the papers.
#' @param comparator_bw Bandwidth of the comparator. Defaults to
#'   [bw_comparator_cv()].
#' @param bw_grid Optional grid of candidate bandwidths. When supplied, the
#'   weights are re-optimised at each one and the pair `(h, p)` with the
#'   smallest distance is returned, which is the joint minimisation over
#'   \eqn{\theta = (h, p)} that the paper defines.
#' @param breaks How the block boundaries are chosen when `m` is finite. One of:
#'   \describe{
#'     \item{`"optimal"`}{(default) The boundaries are chosen jointly with the
#'       weights, minimising the same criterion. This is what Section 4.1 of the
#'       2016 paper specifies: the breakpoints \eqn{r_1, r_2} are part of the
#'       optimisation. Exhaustive for `m <= 3`, coordinate descent above.}
#'     \item{`"equal"`}{Blocks of near-equal size, Algorithm A of the 2018
#'       paper. Much cheaper, and often nearly as good.}
#'     \item{`"modal"`}{Boundaries placed at the troughs of a pilot density
#'       estimate, following the practical suggestion in Section 2 of the 2018
#'       paper.}
#'   }
#'   A numeric vector of `m - 1` breakpoints, given as ranks in the sorted
#'   sample, may be supplied instead.
#' @param min_block Smallest number of observations allowed in a block, used by
#'   `"optimal"` and `"modal"`. Defaults to `max(2, floor(0.02 * n))`.
#' @param constraint Shape constraint imposed on the fitted density: `"none"`
#'   (default), `"unimodal"`, `"increasing"` or `"decreasing"`. The estimator is
#'   linear in the weights, so a shape requirement evaluated on a grid becomes a
#'   set of linear inequalities and the problem stays a convex quadratic
#'   program. Univariate only. See Details.
#' @param constraint_range Interval on which the shape constraint is enforced.
#'   Defaults to `range(x)`. It cannot usefully extend past the observations: a
#'   Gaussian kernel mixture must rise from zero on the left and fall to zero on
#'   the right, so a monotonicity requirement outside the data range has no
#'   feasible solution.
#' @param constraint_grid_n Number of points at which the shape constraint is
#'   enforced (default 100). The constraint holds exactly on this grid; between
#'   grid points a negligible violation is possible, so raise this if you need a
#'   tighter guarantee.
#' @param n,from,to Grid on which to evaluate the fitted density.
#'
#' @return An object of class `"tiltdens"`, which inherits from `"density"`, so
#'   `plot()`, `lines()` and `points()` work directly. Components beyond those
#'   of [stats::density()] include `weights`, `distance2` (the attained squared
#'   L2 distance), `comparator`, `comparator_bw`, `n_groups` and `group_of`.
#'
#' @references
#' Doosti, H. and Hall, P. (2016). Making a non-parametric density estimator
#' more attractive, and more accurate, by data perturbation.
#' *Journal of the Royal Statistical Society B* **78**, 445-462.
#'
#' @seealso [tilt_density_cv()] for the faster cross-validation variant,
#'   [sharpen_density()] for perturbing the observations instead of their
#'   weights.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(25, -1.5), rnorm(25, 1.5))
#'
#' fit <- tilt_density(x, m = 3)
#' fit
#' plot(fit)
#' lines(sinc_density(x), col = "red", lty = 2)
#' abline(h = 0, col = "grey")
#'
#' ## The comparator dips below zero; the tilted estimate cannot.
#' min(sinc_density(x)$y)
#' min(fit$y)
#'
#' \donttest{
#' ## Choosing the block boundaries, rather than fixing them, lowers the
#' ## criterion the method is minimising. Searching them is the slow option:
#' ## it examines every admissible pair, so its cost grows with n^2.
#' c(equal   = tilt_density(x, m = 3, breaks = "equal")$distance2,
#'   modal   = tilt_density(x, m = 3, breaks = "modal")$distance2,
#'   optimal = tilt_density(x, m = 3, breaks = "optimal")$distance2)
#' }
#'
#' @export
tilt_density <- function(x, m = Inf, comparator = c("sinc", "trapezoid"),
                         bw = NULL, comparator_bw = NULL, bw_grid = NULL,
                         kernel = "gaussian",
                         breaks = c("optimal", "equal", "modal"),
                         min_block = NULL, constraint = "none",
                         constraint_range = NULL, constraint_grid_n = 100L,
                         n = 512, from = NULL, to = NULL) {
  data_name  <- deparse(substitute(x))
  comparator <- match_comparator(match.arg(comparator))

  if (is_multivariate(x)) {
    if (!identical(constraint, "none")) {
      stop("Shape constraints are available for univariate data only.",
           call. = FALSE)
    }
    if (is.finite(m)) {
      warning("Multivariate fits use one weight per observation; `m` ignored.",
              call. = FALSE)
    }
    return(tilt_density_md(x, comparator, bw, comparator_bw, bw_grid,
                           data_name, match.call(), n,
                           criterion = "comparator",
                           kernel = tilt_kernel(kernel)))
  }

  x          <- check_sample(x)
  constraint <- match_constraint(constraint)
  kernel     <- tilt_kernel(kernel)
  if (!is.numeric(breaks)) breaks <- match.arg(breaks)
  n_obs      <- length(x)
  if (is.null(min_block)) min_block <- default_min_block(n_obs, m)

  if (is.null(comparator_bw)) {
    comparator_bw <- as.numeric(bw_comparator_cv(x, comparator = comparator))
  }
  check_bw(comparator_bw)

  h_candidates <- if (!is.null(bw_grid)) {
    as.vector(bw_grid)
  } else if (!is.null(bw)) {
    bw
  } else {
    ## The comparator bandwidth is on the comparator's scale, which is the
    ## Gaussian's for the papers' setting. Rescale it canonically so that
    ## changing the kernel changes the shape of the fit, not how much it is
    ## smoothed. For the Gaussian the factor is one, so the papers' behaviour
    ## is unchanged.
    comparator_bw * canonical_ratio(kernel, "gaussian")
  }
  for (h in h_candidates) check_bw(h)

  const_term <- comparator_normsq(x, comparator_bw, comparator)$value

  distances <- numeric(length(h_candidates))
  weights   <- matrix(0, n_obs, length(h_candidates))
  fits      <- vector("list", length(h_candidates))

  for (k in seq_along(h_candidates)) {
    h <- h_candidates[k]
    A <- self_conv_gram(x, h, kernel)
    b <- comparator_cross(x, x, h, comparator_bw, comparator, kernel)

    fit <- fit_blocked_weights(A, b, x, m, breaks, min_block, constraint, h,
                                 constraint_range, constraint_grid_n, kernel)

    fits[[k]]    <- fit
    weights[, k] <- fit$weights
    distances[k] <- fit$value + const_term
  }

  best <- which.min(distances)
  chosen <- fits[[best]]
  new_tiltdens(
    x = x, bw = h_candidates[best], weights = weights[, best],
    n = n, from = from, to = to, data_name = data_name,
    call = match.call(), method = "tilt_density", kernel = kernel,
    extra = list(
      distance2     = distances[best],
      comparator    = comparator,
      comparator_bw = comparator_bw,
      n_groups      = chosen$n_groups,
      group_of      = chosen$group_of,
      breaks        = chosen$breaks,
      break_values  = if (length(chosen$breaks)) sort(x)[chosen$breaks] else numeric(0),
      breaks_method = if (is.numeric(breaks)) "user" else breaks,
      constraint    = constraint,
      kernel        = kernel,
      mode_location = chosen$mode,
      bw_grid       = if (length(h_candidates) > 1L) h_candidates else NULL,
      distance_grid = if (length(h_candidates) > 1L) distances else NULL
    )
  )
}
