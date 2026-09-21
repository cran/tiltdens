#' Cross-validation criterion for a tilted estimator
#'
#' Evaluates
#' \deqn{CV(h, p) = \int \hat f(x \mid h,p)^2 dx -
#'                  \frac{2}{n}\sum_i \hat f_{-i}(x_i \mid h,p),}
#' equation (2.1) of Doosti, Hall and Mateu (2018), where \eqn{\hat f_{-i}} is
#' the estimator with the term in \eqn{x_i} deleted. As in the paper, `p` is
#' *not* re-standardised after deletion; that step turns out to be unnecessary
#' and skipping it saves computation.
#'
#' The first term is \eqn{p^\top A p} and the second is linear in `p`, so this
#' criterion has the same quadratic form as the 2016 distance criterion. The
#' only difference is the linear term: the 2016 method compares against a
#' comparator estimator, the 2018 method against the leave-one-out fit. That is
#' why one solver serves both.
#'
#' @param x Numeric vector of observations.
#' @param bw Positive bandwidth.
#' @param kernel Kernel of the estimator; see [tilt_kernels()].
#' @param weights Tilt weights. Defaults to the uniform weights `1/n`, in which
#'   case this reduces to ordinary least-squares cross-validation for the
#'   bandwidth.
#'
#' @return A single number, with attributes `"A"` and `"C"` holding the
#'   quadratic and linear parts, which [tilt_density_cv()] reuses.
#'
#' @references
#' Doosti, H., Hall, P. and Mateu, J. (2018). Nonparametric tilted density
#' function estimation: a cross-validation criterion.
#' *Journal of Statistical Planning and Inference* **197**, 51-68.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(60)
#' tilt_cv(x, bw = 0.4)
#'
#' @export
tilt_cv <- function(x, bw, weights = NULL, kernel = "gaussian") {
  x <- check_sample(x)
  check_bw(bw)
  kernel <- tilt_kernel(kernel)
  n <- length(x)
  if (is.null(weights)) weights <- rep(1 / n, n)
  weights <- as.vector(weights)
  if (length(weights) != n) {
    stop("`weights` must have the same length as `x`.", call. = FALSE)
  }

  A <- self_conv_gram(x, bw, kernel)

  M <- kernel_matrix(x, x, bw, kernel)
  diag(M) <- 0                      # omit the j = i terms
  cvec <- colSums(M)                # c_j = sum over i != j of K_h(x_i - x_j)

  value <- as.numeric(t(weights) %*% A %*% weights) - 2 * sum(cvec * weights) / n
  structure(value, A = A, C = cvec)
}

#' Tilted density estimation by cross-validation
#'
#' Chooses the bandwidth and the tilt weights together by minimising
#' [tilt_cv()], rather than by comparing against an infinite-order estimator.
#' This is method (III) of Doosti, Hall and Mateu (2018). It needs no comparator
#' and is much cheaper than [tilt_density()], which is the point: the 2016
#' construction is accurate but slow, and this makes it practical.
#'
#' The sequential scheme of the paper is used. Start with uniform weights and
#' choose the bandwidth by ordinary least-squares cross-validation. At step `r`,
#' minimise `CV(h, p)` jointly over `h` and over weights taking `r + 1` distinct
#' values, subject to `h >= rho * h_previous`; the constraint reflects the
#' optimal bandwidth growing as the bias shrinks. Stop at `max_groups` blocks.
#'
#' @param x Numeric vector of observations.
#' @param max_groups Largest number of distinct weight values (default 3).
#' @param bw_grid Candidate bandwidths. Defaults to 40 points spaced
#'   logarithmically over \eqn{[n^{-1/5-c_1}, n^{-1/5+c_1}]} times a robust
#'   scale, matching the class `H` of equation (3.1) with \eqn{c_1 = 1/30}.
#' @param rho Constant \eqn{\rho \in (0,1)} of equation (2.6) (default 0.8).
#' @param kernel Kernel of the estimator; see [tilt_kernels()].
#' @param breaks How the block boundaries are chosen. `"equal"` (the default)
#'   gives blocks of near-equal size, which is Algorithm A of the paper.
#'   `"modal"` places them at the troughs of a pilot density estimate, the
#'   refinement the paper suggests for multimodal densities. `"optimal"` chooses
#'   them jointly with the weights, which is more expensive but usually gives a
#'   smaller criterion. A numeric vector of breakpoints may be supplied instead.
#'   See [tilt_density()] for the full description.
#' @param min_block Smallest number of observations allowed in a block.
#' @param constraint_range,constraint_grid_n As in [tilt_density()].
#' @param constraint Shape constraint, as in [tilt_density()]. Applying a shape
#'   constraint to the cross-validation criterion follows the extension
#'   suggested in the discussion of the 2018 paper.
#' @param n,from,to Grid on which to evaluate the fitted density.
#'
#' @return An object of class `"tiltdens"`; see [tilt_density()]. The `trace`
#'   component records the number of blocks, bandwidth and criterion value at
#'   each step, which is useful for seeing whether extra blocks bought anything.
#'
#' @references
#' Doosti, H., Hall, P. and Mateu, J. (2018). Nonparametric tilted density
#' function estimation: a cross-validation criterion.
#' *Journal of Statistical Planning and Inference* **197**, 51-68.
#'
#' @seealso [tilt_density()], [tilt_cv()]
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
#'
#' fit <- tilt_density_cv(x)
#' fit
#' fit$trace
#' plot(fit)
#'
#' @export
tilt_density_cv <- function(x, max_groups = 3, bw_grid = NULL, rho = 0.8,
                            kernel = "gaussian",
                            breaks = c("equal", "modal", "optimal"),
                            min_block = NULL, constraint = "none",
                            constraint_range = NULL, constraint_grid_n = 100L,
                            n = 512, from = NULL, to = NULL) {
  data_name <- deparse(substitute(x))

  if (is_multivariate(x)) {
    if (!identical(constraint, "none")) {
      stop("Shape constraints are available for univariate data only.",
           call. = FALSE)
    }
    return(tilt_density_md(x, "sinc", NULL, NULL, bw_grid, data_name,
                           match.call(), n, criterion = "cv", rho = rho,
                           kernel = tilt_kernel(kernel)))
  }

  x         <- check_sample(x)
  constraint <- match_constraint(constraint)
  kernel    <- tilt_kernel(kernel)
  if (!is.numeric(breaks)) breaks <- match.arg(breaks)
  n_obs     <- length(x)
  if (is.null(min_block)) min_block <- default_min_block(n_obs, max_groups)

  if (is.null(bw_grid)) {
    c1 <- 1 / 30
    bw_grid <- robust_scale(x) * canonical_ratio(kernel, "gaussian") *
      exp(seq(log(n_obs^(-1 / 5 - c1)), log(n_obs^(-1 / 5 + c1)), length.out = 40))
  }
  bw_grid <- as.vector(bw_grid)

  ## Step 1: uniform weights, ordinary least-squares cross-validation.
  uniform  <- rep(1 / n_obs, n_obs)
  cv_grid  <- vapply(bw_grid, function(h) as.numeric(tilt_cv(x, h, uniform)), numeric(1))
  best_ix  <- which.min(cv_grid)

  cv_best  <- cv_grid[best_ix]
  bw_best  <- bw_grid[best_ix]
  p_best   <- uniform
  m_best   <- 1L
  group_of <- rep(1L, n_obs)
  breaks_best <- integer(0)

  trace <- data.frame(n_groups = 1L, bw = bw_best, cv = cv_best)

  ## Steps 2, 3, ...: more distinct weights, bandwidth restricted from below.
  for (m in seq_len(max_groups)[-1L]) {
    allowed <- bw_grid[bw_grid >= rho * bw_best]
    if (!length(allowed)) allowed <- bw_best

    step_cv <- Inf
    step_bw <- bw_best
    step_p  <- p_best
    step_fit <- NULL

    for (h in allowed) {
      parts <- tilt_cv(x, h, uniform, kernel)
      fit <- fit_blocked_weights(attr(parts, "A"), attr(parts, "C"),
                                 x, m, breaks, min_block, constraint, h,
                                 constraint_range, constraint_grid_n, kernel)
      if (fit$value < step_cv) {
        step_cv  <- fit$value
        step_bw  <- h
        step_p   <- fit$weights
        step_fit <- fit
      }
    }

    trace <- rbind(trace, data.frame(n_groups = m, bw = step_bw, cv = step_cv))

    if (step_cv < cv_best) {
      cv_best     <- step_cv
      bw_best     <- step_bw
      p_best      <- step_p
      m_best      <- m
      group_of    <- step_fit$group_of
      breaks_best <- step_fit$breaks
    }
  }

  new_tiltdens(
    x = x, bw = bw_best, weights = p_best,
    n = n, from = from, to = to, data_name = data_name,
    call = match.call(), method = "tilt_density_cv", kernel = kernel,
    extra = list(
      cv       = cv_best,
      n_groups      = m_best,
      group_of      = group_of,
      breaks        = breaks_best,
      break_values  = if (length(breaks_best)) sort(x)[breaks_best] else numeric(0),
      breaks_method = if (is.numeric(breaks)) "user" else breaks,
      constraint    = constraint,
      kernel        = kernel,
      trace         = trace
    )
  )
}
