#' Minimise a quadratic form over a weighted probability simplex
#'
#' Solves \eqn{\min_q\; q^\top H q - 2 f^\top q} subject to
#' \eqn{s^\top q = 1} and \eqn{q \ge 0}, for symmetric positive semidefinite
#' `H`. With `s` the vector of block sizes this is exactly the weight-selection
#' problem of both papers.
#'
#' The problem is convex, so it has a unique solution and there is no need for a
#' derivative-free search. [quadprog::solve.QP()] is tried first; if the
#' quadratic form is too ill-conditioned for it even after ridging, an
#' accelerated projected-gradient method takes over, which is slower but always
#' applicable.
#'
#' @param H Symmetric positive semidefinite matrix.
#' @param f Numeric vector, the linear coefficients.
#' @param s Positive vector of constraint weights. Defaults to all ones.
#' @param max_iter,tol Controls for the fallback solver.
#'
#' @return A list with the minimiser `q` and the attained `value`.
#'
#' @keywords internal
#' @noRd
solve_simplex_qp <- function(H, f, s = NULL, max_iter = 20000L, tol = 1e-12) {
  m <- length(f)
  f <- as.vector(f)
  if (is.null(s)) s <- rep(1, m)
  s <- as.vector(s)

  H <- (H + t(H)) / 2

  ## For a handful of variables an exact enumeration of the supports is both
  ## faster and more accurate than a general QP routine.
  if (m <= 6L) return(solve_simplex_qp_small(H, f, s))

  ## Rescale by the block sizes so the constraint becomes the standard simplex.
  Hs <- H / outer(s, s)
  fs <- f / s
  Hs <- (Hs + t(Hs)) / 2

  u <- s / sum(s)          # uniform weights p_i = 1/n
  solved <- FALSE

  if (requireNamespace("quadprog", quietly = TRUE)) {
    Amat <- cbind(rep(1, m), diag(m))
    bvec <- c(1, rep(0, m))
    ridge <- 1e-10 * max(mean(diag(Hs)), .Machine$double.eps)
    for (attempt in seq_len(12L)) {
      fit <- try(
        quadprog::solve.QP(
          Dmat = 2 * (Hs + diag(ridge, m)),
          dvec = 2 * fs, Amat = Amat, bvec = bvec, meq = 1L
        ),
        silent = TRUE
      )
      if (!inherits(fit, "try-error") && all(is.finite(fit$solution))) {
        u <- fit$solution
        solved <- TRUE
        break
      }
      ridge <- ridge * 100
    }
  }

  if (!solved) u <- fista_simplex(Hs, fs, u, max_iter, tol)

  u <- project_simplex(u)
  q <- u / s
  list(q = q, value = as.numeric(t(q) %*% H %*% q - 2 * sum(f * q)))
}

#' Accelerated projected gradient on the standard simplex
#' @keywords internal
#' @noRd
fista_simplex <- function(H, f, u0, max_iter, tol) {
  lipschitz <- 2 * norm(H, type = "2")
  if (!is.finite(lipschitz) || lipschitz <= 0) lipschitz <- 1
  step <- 1 / lipschitz

  u <- u0
  z <- u
  t_prev <- 1
  for (it in seq_len(max_iter)) {
    grad  <- 2 * as.vector(H %*% z) - 2 * f
    u_new <- project_simplex(z - step * grad)
    t_new <- (1 + sqrt(1 + 4 * t_prev^2)) / 2
    z <- u_new + ((t_prev - 1) / t_new) * (u_new - u)
    if (sum(abs(u_new - u)) < tol) return(u_new)
    u <- u_new
    t_prev <- t_new
  }
  u
}

#' Euclidean projection onto the standard simplex
#' @keywords internal
#' @noRd
project_simplex <- function(v) {
  v <- as.vector(v)
  m <- length(v)
  u <- sort(v, decreasing = TRUE)
  css <- cumsum(u)
  idx <- seq_len(m)
  rho <- which(u > (css - 1) / idx)
  if (!length(rho)) return(rep(1 / m, m))
  rho <- max(rho)
  pmax(v - (css[rho] - 1) / rho, 0)
}

#' Exact solver for a simplex QP in a handful of variables
#'
#' Enumerates the possible supports. For a convex quadratic program over
#' \eqn{\{q \ge 0,\ s^\top q = 1\}} the optimum has some set of strictly
#' positive coordinates; solving the equality-constrained problem on every
#' candidate support, discarding infeasible solutions and keeping the best gives
#' the global optimum exactly.
#'
#' With `m` at most six this costs at most 63 tiny linear solves, which is far
#' quicker than calling a general quadratic programming routine. That matters
#' because [search_breaks()] solves one of these for every candidate set of
#' breakpoints, of which there are order \eqn{n^2} when `m = 3`.
#'
#' @keywords internal
#' @noRd
solve_simplex_qp_small <- function(H, f, s) {
  m <- length(f)
  H <- (H + t(H)) / 2

  best_q <- NULL
  best_value <- Inf

  for (code in seq_len(2^m - 1L)) {
    support <- which(bitwAnd(code, 2^(seq_len(m) - 1L)) > 0L)
    k <- length(support)

    Hs <- H[support, support, drop = FALSE]
    fs <- f[support]
    ss <- s[support]

    kkt <- rbind(cbind(2 * Hs, ss), c(ss, 0))
    rhs <- c(2 * fs, 1)

    sol <- try(solve(kkt, rhs), silent = TRUE)
    if (inherits(sol, "try-error")) next

    q_sup <- sol[seq_len(k)]
    if (any(!is.finite(q_sup)) || any(q_sup < -1e-12)) next

    q <- numeric(m)
    q[support] <- pmax(q_sup, 0)
    value <- as.numeric(t(q) %*% H %*% q - 2 * sum(f * q))

    if (value < best_value) {
      best_value <- value
      best_q <- q
    }
  }

  if (is.null(best_q)) {
    best_q <- s / sum(s^2)
    best_value <- as.numeric(t(best_q) %*% H %*% best_q - 2 * sum(f * best_q))
  }

  list(q = best_q, value = best_value)
}
