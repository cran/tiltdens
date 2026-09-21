#' Shape constraints for tilted density estimation
#'
#' The tilted estimator \eqn{\hat f(x) = \sum_i p_i K_h(x - x_i)} is *linear* in
#' the weights. So a shape requirement evaluated on a grid becomes a set of
#' linear inequalities in `p`, and imposing it leaves the optimisation a convex
#' quadratic program. Nothing about the method changes; the feasible set simply
#' gains some rows.
#'
#' Writing \eqn{\Phi_{kj} = K_h(t_k - x_j)} for the grid \eqn{t_1 < \dots < t_M},
#' the fitted curve is \eqn{\Phi p}, so
#'
#' * increasing on the grid is \eqn{(\Phi_{k+1,\cdot} - \Phi_{k,\cdot})\,p \ge 0};
#' * decreasing is the same with the sign reversed;
#' * unimodal with the mode at grid point \eqn{k_0} is increasing up to
#'   \eqn{k_0} and decreasing after it.
#'
#' The mode location is not known in advance, so for `"unimodal"` the problem is
#' solved once for each candidate mode on a coarse grid and the best solution
#' kept. That is the same device Hall and Huang (2002) use for unimodal density
#' estimation by tilting, which is the published antecedent of this feature.
#'
#' @section What the constraint guarantees:
#' The shape requirement is imposed at the `constraint_grid_n` grid points, and
#' it holds there to solver tolerance. Between grid points a violation is
#' possible in principle. In practice it is negligible: on a fit where the
#' unconstrained estimate has a clear spurious mode, the constrained fit is flat
#' across that region and ripples by a relative amount of order 1e-5. Raise
#' `constraint_grid_n` if you need a tighter guarantee; the cost is linear in
#' the number of grid points.
#'
#' @section Provenance:
#' Shape-constrained tilting is due to Hall and Huang (2002), which both of the
#' papers implemented here cite. Applying the same device to the
#' cross-validation criterion of Doosti, Hall and Mateu (2018) follows the
#' extension suggested in the discussion of that paper.
#'
#' @references
#' Hall, P. and Huang, L. S. (2002). Unimodal density estimation using kernel
#' methods. *Statistica Sinica* **12**, 965-990.
#'
#' @name constraints
#' @keywords internal
#' @noRd
NULL

#' Build the constraint matrix for a shape requirement
#'
#' @param x Numeric sample.
#' @param h Bandwidth of the tilted estimator.
#' @param constraint `"increasing"`, `"decreasing"`, or `"unimodal"` with the
#'   mode fixed at index `mode_index` of the grid.
#' @param grid Grid on which the shape is enforced.
#' @param mode_index Index of the mode within `grid`, for `"unimodal"`.
#'
#' @return A matrix `C` such that the constraint is `C %*% p >= 0`.
#'
#' @keywords internal
#' @noRd
constraint_matrix <- function(x, h, constraint, grid, mode_index = NULL,
                              kernel = NULL) {
  Phi <- kernel_matrix(grid, x, h, kernel)
  D <- Phi[-1L, , drop = FALSE] - Phi[-nrow(Phi), , drop = FALSE]

  switch(
    constraint,
    increasing = D,
    decreasing = -D,
    unimodal = {
      if (is.null(mode_index)) stop("`mode_index` is required.", call. = FALSE)
      up   <- seq_len(nrow(D)) < mode_index
      rbind(D[up, , drop = FALSE], -D[!up, , drop = FALSE])
    },
    stop("Unknown constraint \"", constraint, "\".", call. = FALSE)
  )
}

#' Minimise a quadratic form over the simplex subject to linear inequalities
#'
#' Solves \eqn{\min_q\ q^\top H q - 2 f^\top q} subject to \eqn{s^\top q = 1},
#' \eqn{q \ge 0} and \eqn{C q \ge 0}. The extra rows keep the problem convex, so
#' this is still a quadratic program with a unique solution.
#'
#' @return A list with `q`, `value` and `feasible`. When the constraint set is
#'   empty, `feasible` is `FALSE` and `q` is `NULL`.
#'
#' @keywords internal
#' @noRd
solve_constrained_simplex_qp <- function(H, f, s, C) {
  m <- length(f)
  H <- (H + t(H)) / 2
  Hs <- H / outer(s, s)
  Hs <- (Hs + t(Hs)) / 2
  fs <- f / s

  ## The constraint acts on q = u / s, so rescale its columns to match.
  Cs <- sweep(C, 2L, s, "/")

  Amat <- cbind(rep(1, m), diag(m), t(Cs))
  bvec <- c(1, rep(0, m), rep(0, nrow(Cs)))

  ridge <- 1e-10 * max(mean(diag(Hs)), .Machine$double.eps)
  for (attempt in seq_len(12L)) {
    fit <- try(
      quadprog::solve.QP(Dmat = 2 * (Hs + diag(ridge, m)), dvec = 2 * fs,
                         Amat = Amat, bvec = bvec, meq = 1L),
      silent = TRUE
    )
    if (!inherits(fit, "try-error") && all(is.finite(fit$solution))) {
      q <- pmax(fit$solution, 0) / s
      q <- q / sum(s * q)
      return(list(q = q, feasible = TRUE,
                  value = as.numeric(t(q) %*% H %*% q - 2 * sum(f * q))))
    }
    ridge <- ridge * 100
  }

  list(q = NULL, value = Inf, feasible = FALSE)
}

#' Fit weights under a shape constraint
#'
#' @keywords internal
#' @noRd
fit_constrained_weights <- function(A, linear, x, h, constraint, G,
                                    group_size, constraint_range = NULL,
                                    constraint_grid_n = 100L,
                                    mode_candidates = 25L, kernel = NULL) {
  n <- length(x)
  ## The constraint is enforced on the data range, not beyond it. A Gaussian
  ## kernel mixture necessarily rises from zero as the data are approached from
  ## the left and falls to zero on the right, so asking for a monotone fit on a
  ## region extending past the observations is infeasible by construction.
  if (is.null(constraint_range)) constraint_range <- range(x)
  grid <- seq(constraint_range[1L], constraint_range[2L],
              length.out = constraint_grid_n)

  H <- t(G) %*% A %*% G
  f <- as.vector(t(G) %*% linear) / n

  if (constraint %in% c("increasing", "decreasing")) {
    C <- constraint_matrix(x, h, constraint, grid, kernel = kernel) %*% G
    fit <- solve_constrained_simplex_qp(H, f, group_size, C)
    if (!fit$feasible) stop(infeasible_message(constraint, ncol(G)), call. = FALSE)
    return(list(q = fit$q, value = fit$value, mode = NA_real_))
  }

  ## Unimodal: the mode location is unknown, so try candidates and keep the best.
  candidates <- unique(round(seq(3, constraint_grid_n - 2L,
                                 length.out = min(mode_candidates,
                                                  constraint_grid_n - 4L))))
  best <- list(q = NULL, value = Inf, mode = NA_real_)
  for (k in candidates) {
    C <- constraint_matrix(x, h, "unimodal", grid, mode_index = k,
                           kernel = kernel) %*% G
    fit <- solve_constrained_simplex_qp(H, f, group_size, C)
    if (fit$feasible && fit$value < best$value) {
      best <- list(q = fit$q, value = fit$value, mode = grid[k])
    }
  }
  if (is.null(best$q)) {
    stop(infeasible_message("unimodal", ncol(G)), call. = FALSE)
  }
  best
}

#' Normalise a constraint argument
#' @keywords internal
#' @noRd
match_constraint <- function(constraint) {
  constraint <- tolower(as.character(constraint)[1L])
  valid <- c("none", "unimodal", "increasing", "decreasing")
  if (!constraint %in% valid) {
    stop("`constraint` must be one of ", paste0("\"", valid, "\"", collapse = ", "),
         ".", call. = FALSE)
  }
  constraint
}

#' Explain an infeasible constraint
#'
#' A shape constraint has to be satisfiable within the weights available. With
#' `m` distinct weights there are only `m` degrees of freedom, and a shape
#' requirement evaluated on a grid of sixty points may simply have no solution.
#' Saying so is more useful than reporting that the solver failed.
#'
#' @keywords internal
#' @noRd
infeasible_message <- function(constraint, n_groups) {
  paste0(
    "No weights satisfy the \"", constraint, "\" constraint with ", n_groups,
    if (n_groups == 1L) " distinct weight." else " distinct weights.",
    if (n_groups < 10L)
      " Shape constraints need enough freedom to work with; try m = Inf."
    else
      " Try widening `constraint_range` or loosening the constraint."
  )
}
