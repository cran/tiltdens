#' Objective for data sharpening
#'
#' Returns \eqn{\|\hat f(\cdot \mid h, q) - \check f\|^2}, where
#' \eqn{\hat f(x \mid h,q) = (nh)^{-1}\sum_i K((x - x_i - q_i)/h)} is the
#' data-sharpened estimator, equation (2.2) of Doosti and Hall (2016), and
#' \eqn{\check f} is the comparator.
#'
#' Unlike the tilting problem this objective is not convex in `q`, which is why
#' [sharpen_density()] uses a stochastic search.
#'
#' @param shifts Numeric vector of shifts, one per observation.
#' @param x Numeric vector of observations.
#' @param bw,comparator_bw,comparator As in [sharpen_density()].
#' @param const_term Optional constant added to the result; pass the squared
#'   norm of the comparator to obtain the full squared distance.
#'
#' @return A single number.
#'
#' @keywords internal
#' @noRd
sharpen_objective <- function(shifts, x, bw, comparator_bw, comparator,
                              const_term = 0, kernel = NULL) {
  n <- length(x)
  y <- x + shifts
  A <- self_conv_gram(y, bw, kernel)
  b <- comparator_cross(y, x, bw, comparator_bw, comparator, kernel)
  sum(A) / n^2 - 2 * sum(b) / n^2 + const_term
}

#' Data-sharpened density estimation
#'
#' Finds shifts `q` for the estimator
#' \deqn{\hat f(x \mid h, q) = \frac{1}{nh}\sum_i K\!\left(\frac{x - x_i - q_i}{h}\right)}
#' by minimising its L2 distance to an infinite-order comparator. This is the
#' second form of data perturbation in Doosti and Hall (2016): the observations
#' are moved rather than re-weighted.
#'
#' The objective is not convex, so the result depends on the random seed. Set
#' one before calling if you need reproducibility. With `m = Inf` the search is
#' over `length(x)` variables and is slow; the paper reports run times of a few
#' minutes per sample of size 100. In the published simulations tilting was
#' usually at least as accurate and far cheaper, so [tilt_density()] and
#' [tilt_density_cv()] are the better default.
#'
#' @param x Numeric vector of observations.
#' @param m Number of distinct shift values. `Inf` (default) or 3.
#' @param comparator `"sinc"` (default) or `"trapezoid"`.
#' @param kernel Kernel of the estimator; see [tilt_kernels()].
#' @param bw,comparator_bw,breaks,min_block As in [tilt_density()]. Breakpoint
#'   search is not available here, so `breaks` defaults to `"equal"`; supply
#'   `"modal"` or an explicit numeric vector if you want different boundaries.
#' @param max_shift Half-width of the box searched, as a multiple of the data
#'   range (default 1).
#' @param polish Refine the search result with [stats::optim()] (default `TRUE`).
#' @param control A list of options passed to the genetic algorithm; see
#'   [real_ga()].
#' @param n,from,to Grid on which to evaluate the fitted density.
#'
#' @return An object of class `"tiltdens"`; see [tilt_density()]. The `shifts`
#'   and `sharpened` components hold `q` and `x + q`.
#'
#' @references
#' Doosti, H. and Hall, P. (2016). Making a non-parametric density estimator
#' more attractive, and more accurate, by data perturbation.
#' *Journal of the Royal Statistical Society B* **78**, 445-462.
#'
#' @seealso [tilt_density()], [tilt_density_cv()]
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(15, -1.5), rnorm(15, 1.5))
#'
#' ## A deliberately tiny search budget, so that this runs in about a second.
#' ## It is far too small for real use; see below for a realistic call.
#' fit <- sharpen_density(x, m = 3, polish = FALSE,
#'                        control = list(max_iterations = 10, population = 10))
#' plot(fit)
#' round(unique(fit$shifts), 3)
#'
#' \donttest{
#' ## Realistic settings. Expect this to take a few seconds, and much longer
#' ## with m = Inf, where the search ranges over one shift per observation.
#' x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
#' fit <- sharpen_density(x, m = 3)
#' fit
#' }
#'
#' @export
sharpen_density <- function(x, m = Inf, comparator = c("sinc", "trapezoid"),
                            bw = NULL, comparator_bw = NULL, kernel = "gaussian",
                            breaks = c("equal", "modal"), min_block = NULL,
                            max_shift = 1, polish = TRUE, control = list(),
                            n = 512, from = NULL, to = NULL) {
  data_name  <- deparse(substitute(x))
  x          <- check_sample(x)
  comparator <- match_comparator(match.arg(comparator))
  if (!is.numeric(breaks)) breaks <- match.arg(breaks)
  kernel <- tilt_kernel(kernel)
  if (is.null(min_block)) min_block <- default_min_block(length(x), m)

  if (is.null(comparator_bw)) {
    comparator_bw <- as.numeric(bw_comparator_cv(x, comparator = comparator))
  }
  check_bw(comparator_bw)
  if (is.null(bw)) bw <- comparator_bw * canonical_ratio(kernel, "gaussian")
  check_bw(bw)

  groups <- sharpen_blocks(x, m, breaks, min_block)
  n_var  <- ncol(groups$G)

  const_term <- comparator_normsq(x, comparator_bw, comparator)$value
  objective <- function(v) {
    sharpen_objective(as.vector(groups$G %*% v), x, bw, comparator_bw,
                      comparator, const_term, kernel)
  }

  span  <- max_shift * diff(range(x))
  lower <- rep(-span, n_var)
  upper <- rep(span, n_var)

  ga <- do.call(real_ga, c(list(objective, lower, upper), control))
  best  <- ga$par
  value <- ga$value

  if (polish) {
    fit <- try(stats::optim(best, objective, method = "Nelder-Mead",
                            control = list(maxit = 200 * n_var)), silent = TRUE)
    if (!inherits(fit, "try-error") && fit$value < value) {
      best  <- fit$par
      value <- fit$value
    }
  }

  shifts <- as.vector(groups$G %*% best)

  new_tiltdens(
    x = x + shifts, bw = bw, weights = rep(1 / length(x), length(x)),
    n = n, from = from, to = to, data_name = data_name,
    call = match.call(), method = "sharpen_density", kernel = kernel,
    extra = list(
      shifts        = shifts,
      sharpened     = x + shifts,
      original      = x,
      distance2     = value,
      comparator    = comparator,
      comparator_bw = comparator_bw,
      n_groups      = n_var,
      group_of      = groups$group_of,
      breaks        = groups$breaks,
      breaks_method = if (is.numeric(breaks)) "user" else breaks,
      history       = ga$history
    )
  )
}

#' Real-coded genetic algorithm for small box-constrained problems
#'
#' Minimises `fn` over the box `[lower, upper]`. This is a tidied,
#' dependency-free version of the algorithm used for the data sharpening results
#' of Doosti and Hall (2016).
#'
#' @param fn Function of a numeric vector, returning a scalar.
#' @param lower,upper Numeric vectors of bounds; their length sets the dimension.
#' @param max_iterations Generations (default 130).
#' @param population Population size (default 50).
#' @param crossover_rate Fraction of the population replaced by offspring (0.8).
#' @param mutation_rate Initial fraction of mutants (0.4, decayed after a third
#'   of the run).
#' @param mutation_step Mutation scale as a fraction of the box width (0.02).
#' @param blend_gamma Blend-crossover expansion factor (0.3).
#' @param selection_pressure Exponent of the fitness-proportional selection (8).
#' @param trace Print the best cost each generation.
#'
#' @return A list with `par`, `value` and `history`.
#'
#' @export
real_ga <- function(fn, lower, upper, max_iterations = 130L, population = 50L,
                    crossover_rate = 0.8, mutation_rate = 0.4,
                    mutation_step = 0.02, blend_gamma = 0.3,
                    selection_pressure = 8, trace = FALSE) {
  lower <- as.vector(lower)
  upper <- as.vector(upper)
  n_var <- length(lower)

  n_pop <- as.integer(population)
  n_off <- 2L * round(crossover_rate * n_pop / 2)

  clamp <- function(v) pmin(pmax(v, lower), upper)

  pos  <- matrix(stats::runif(n_pop * n_var), n_pop, n_var)
  pos  <- sweep(sweep(pos, 2L, upper - lower, "*"), 2L, lower, "+")
  cost <- apply(pos, 1L, fn)
  ord  <- order(cost)
  pos  <- pos[ord, , drop = FALSE]
  cost <- cost[ord]

  worst   <- cost[n_pop]
  history <- numeric(max_iterations)
  mut     <- mutation_rate

  for (it in seq_len(max_iterations)) {
    denom <- if (worst == 0 || !is.finite(worst)) 1 else abs(worst)
    prob  <- exp(-selection_pressure * cost / denom)
    if (!all(is.finite(prob)) || sum(prob) <= 0) prob <- rep(1, n_pop)
    prob <- prob / sum(prob)

    if (it > max_iterations / 3) mut <- max(0.05, mut - 0.005)
    n_mut <- max(1L, round(mut * n_pop))

    ## Crossover
    off <- matrix(0, n_off, n_var)
    for (k in seq(1L, n_off, by = 2L)) {
      if (it < max_iterations / 3) {
        i1 <- sample.int(n_pop, 1L)
        i2 <- sample.int(n_pop, 1L)
      } else {
        i1 <- sample.int(n_pop, 1L, prob = prob)
        i2 <- sample.int(n_pop, 1L, prob = prob)
      }
      alpha <- stats::runif(n_var, -blend_gamma, 1 + blend_gamma)
      off[k, ]      <- clamp(alpha * pos[i1, ] + (1 - alpha) * pos[i2, ])
      off[k + 1L, ] <- clamp(alpha * pos[i2, ] + (1 - alpha) * pos[i1, ])
    }
    off_cost <- apply(off, 1L, fn)

    ## Mutation
    mutants <- matrix(0, n_mut, n_var)
    for (k in seq_len(n_mut)) {
      parent <- pos[sample.int(n_pop, 1L), ]
      j <- sample.int(n_var, max(1L, ceiling(mutation_step * n_var)))
      child <- parent
      child[j] <- parent[j] + 0.1 * (upper[j] - lower[j]) * stats::rnorm(length(j))
      mutants[k, ] <- clamp(child)
    }
    mut_cost <- apply(mutants, 1L, fn)

    all_pos  <- rbind(pos, off, mutants)
    all_cost <- c(cost, off_cost, mut_cost)
    ord <- order(all_cost)
    all_pos  <- all_pos[ord, , drop = FALSE]
    all_cost <- all_cost[ord]

    worst <- max(worst, all_cost[length(all_cost)])
    pos   <- all_pos[seq_len(n_pop), , drop = FALSE]
    cost  <- all_cost[seq_len(n_pop)]

    history[it] <- cost[1L]
    if (trace) {
      message(sprintf("generation %4d   best cost = %.6g", it, cost[1L]))
    }
  }

  list(par = pos[1L, ], value = cost[1L], history = history)
}
