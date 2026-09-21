#' Block structure from a set of breakpoints
#'
#' Given breakpoints \eqn{1 \le i_1 < \dots < i_{m-1} < n} in rank order, returns
#' the block index of each observation. Blocks are
#' \eqn{[1, i_1], [i_1+1, i_2], \dots, [i_{m-1}+1, n]} in the ranked sample,
#' which is the structure both papers use.
#'
#' @param x Numeric vector of observations.
#' @param breaks Integer vector of `m - 1` breakpoints, in rank order.
#'
#' @return A list with `G`, `group_size`, `group_of` and `breaks`, as
#'   as the other block helpers return.
#'
#' @keywords internal
#' @noRd
blocks_from_breaks <- function(x, breaks) {
  n <- length(x)
  breaks <- sort(unique(as.integer(breaks)))
  breaks <- breaks[breaks >= 1L & breaks < n]

  edges <- c(0L, breaks, n)
  label <- rep(seq_along(diff(edges)), times = diff(edges))

  group_of <- integer(n)
  group_of[order(x)] <- label

  m <- length(edges) - 1L
  G <- matrix(0, n, m)
  G[cbind(seq_len(n), group_of)] <- 1

  list(G = G, group_size = colSums(G), group_of = group_of, breaks = breaks)
}

#' Equally spaced breakpoints
#' @keywords internal
#' @noRd
equal_breaks <- function(n, m) {
  if (m <= 1L) return(integer(0))
  as.integer(round(seq(0, n, length.out = m + 1L))[-c(1L, m + 1L)])
}

#' Breakpoints placed near the modes of a pilot estimate
#'
#' Implements the practical suggestion in Section 2 of Doosti, Hall and Mateu
#' (2018): build a pilot kernel estimate, find its modes, and put breakpoints at
#' the ranks of the observations closest to them. Remaining breakpoints, if `m`
#' exceeds the number of modes plus one, are filled in at equal spacing.
#'
#' The rationale is that the weights need to change where the density changes
#' character. A block boundary in the middle of a mode wastes a degree of
#' freedom; one at the foot of a mode does not.
#'
#' @keywords internal
#' @noRd
modal_breaks <- function(x, m, min_block) {
  n <- length(x)
  if (m <= 1L) return(integer(0))

  pilot <- stats::density(x, bw = bw_nrd_robust(x), n = 512)
  y <- pilot$y
  is_mode <- rep(FALSE, length(y))
  if (length(y) > 2L) {
    inner <- 2:(length(y) - 1L)
    is_mode[inner] <- y[inner] > y[inner - 1L] & y[inner] >= y[inner + 1L]
  }
  mode_x <- pilot$x[is_mode]

  ## Put boundaries at the troughs between consecutive modes: those are the
  ## points where the shape of the density changes most.
  candidates <- numeric(0)
  if (length(mode_x) >= 2L) {
    for (k in seq_len(length(mode_x) - 1L)) {
      window <- pilot$x >= mode_x[k] & pilot$x <= mode_x[k + 1L]
      candidates <- c(candidates, pilot$x[window][which.min(y[window])])
    }
  }

  xs <- sort(x)
  ranks <- vapply(candidates, function(cx) which.min(abs(xs - cx)), integer(1))
  ranks <- sort(unique(ranks))

  ## Keep the modal boundaries, then fill any remaining slots by splitting the
  ## widest blocks. Filling with a separate equally spaced grid would let the
  ## fillers displace the modal boundaries, which defeats the purpose.
  if (length(ranks) > m - 1L) {
    ranks <- ranks[seq_len(m - 1L)]
  }
  while (length(ranks) < m - 1L) {
    edges <- c(0L, ranks, n)
    widest <- which.max(diff(edges))
    new_break <- as.integer(round((edges[widest] + edges[widest + 1L]) / 2))
    if (new_break %in% ranks || new_break <= 0L || new_break >= n) break
    ranks <- sort(c(ranks, new_break))
  }

  enforce_min_block(ranks, n, min_block, m)
}

#' Push breakpoints apart so every block has at least `min_block` observations
#' @keywords internal
#' @noRd
enforce_min_block <- function(breaks, n, min_block, m) {
  if (!length(breaks)) return(integer(0))
  breaks <- sort(unique(pmin(pmax(as.integer(breaks), min_block), n - min_block)))
  for (k in seq_along(breaks)) {
    lower <- if (k == 1L) min_block else breaks[k - 1L] + min_block
    breaks[k] <- min(max(breaks[k], lower), n - (length(breaks) - k + 1L) * min_block)
  }
  breaks <- sort(unique(breaks))
  breaks[breaks >= min_block & breaks <= n - min_block]
}

#' Search over breakpoints for the best block structure
#'
#' Chooses the breakpoints and the weights together, minimising the same
#' criterion. This is what Section 4.1 of Doosti and Hall (2016) describes for
#' the m = 3 estimator: the breakpoints \eqn{r_1, r_2} are part of the
#' optimisation, not fixed in advance.
#'
#' For `m <= 3` every admissible set of breakpoints is examined. That is `n - 1`
#' candidates for `m = 2` and about `n^2/2` for `m = 3`, which is affordable
#' because the block sums of `A` and of the linear term are read off
#' two-dimensional cumulative sums in constant time, and the resulting quadratic
#' program has only `m` variables.
#'
#' For `m > 3` exhaustive search is not affordable, so breakpoints are refined
#' by coordinate descent from an equally spaced start: each breakpoint in turn is
#' moved to its best position with the others held fixed, sweeping until nothing
#' improves.
#'
#' @param A Quadratic form in the original observation order.
#' @param linear Linear coefficient vector in the original observation order.
#' @param x The sample, used only for its ranks.
#' @param m Number of blocks.
#' @param min_block Smallest number of observations allowed in a block.
#' @param max_sweeps Coordinate-descent sweeps when `m > 3`.
#'
#' @return A list with `breaks`, `q` (the distinct weight values) and `value`.
#'
#' @keywords internal
#' @noRd
search_breaks <- function(A, linear, x, m, min_block, max_sweeps = 10L) {
  n <- length(x)
  ord <- order(x)

  ## Rank-ordered cumulative sums make any block sum an O(1) lookup.
  As <- A[ord, ord, drop = FALSE]
  cum_A <- rbind(0, cbind(0, apply(apply(As, 1L, cumsum), 1L, cumsum)))
  cum_b <- c(0, cumsum(linear[ord]))

  block_sum_A <- function(r1, r2, c1, c2) {
    cum_A[r2 + 1L, c2 + 1L] - cum_A[r1, c2 + 1L] -
      cum_A[r2 + 1L, c1] + cum_A[r1, c1]
  }

  evaluate <- function(breaks) {
    edges <- c(0L, breaks, n)
    k <- length(edges) - 1L
    lo <- edges[-(k + 1L)] + 1L
    hi <- edges[-1L]
    if (any(hi - lo + 1L < 1L)) return(list(value = Inf))

    Sa <- matrix(0, k, k)
    for (i in seq_len(k)) {
      for (j in i:k) {
        Sa[i, j] <- block_sum_A(lo[i], hi[i], lo[j], hi[j])
        Sa[j, i] <- Sa[i, j]
      }
    }
    sb <- cum_b[hi + 1L] - cum_b[lo]
    sizes <- hi - lo + 1L

    fit <- solve_simplex_qp(Sa, sb / n, sizes)
    list(value = fit$value, q = fit$q, sizes = sizes)
  }

  if (m <= 1L) {
    res <- evaluate(integer(0))
    return(list(breaks = integer(0), q = res$q, value = res$value))
  }

  best <- list(value = Inf)

  if (m == 2L) {
    for (r1 in seq(min_block, n - min_block)) {
      res <- evaluate(r1)
      if (res$value < best$value) best <- c(res, list(breaks = r1))
    }
  } else if (m == 3L) {
    for (r1 in seq(min_block, n - 2L * min_block)) {
      for (r2 in seq(r1 + min_block, n - min_block)) {
        res <- evaluate(c(r1, r2))
        if (res$value < best$value) best <- c(res, list(breaks = c(r1, r2)))
      }
    }
  } else {
    breaks <- enforce_min_block(equal_breaks(n, m), n, min_block, m)
    best <- c(evaluate(breaks), list(breaks = breaks))

    for (sweep in seq_len(max_sweeps)) {
      improved <- FALSE
      for (k in seq_along(breaks)) {
        lower <- if (k == 1L) min_block else breaks[k - 1L] + min_block
        upper <- if (k == length(breaks)) n - min_block else breaks[k + 1L] - min_block
        if (upper < lower) next
        for (candidate in seq(lower, upper)) {
          trial <- breaks
          trial[k] <- candidate
          res <- evaluate(trial)
          if (res$value < best$value - 1e-14) {
            best <- c(res, list(breaks = trial))
            breaks <- trial
            improved <- TRUE
          }
        }
      }
      if (!improved) break
    }
  }

  list(breaks = best$breaks, q = best$q, value = best$value)
}

#' Resolve a breakpoint specification into fitted weights
#'
#' Shared back end for the `breaks` argument of [tilt_density()] and
#' [tilt_density_cv()]. Handles the three named strategies and a user-supplied
#' numeric vector.
#'
#' @keywords internal
#' @noRd
fit_blocked_weights <- function(A, linear, x, m, breaks, min_block,
                                constraint = "none", h = NULL,
                                constraint_range = NULL,
                                constraint_grid_n = 100L, kernel = NULL) {
  n <- length(x)
  constrained <- !identical(constraint, "none")

  if (is.null(m) || !is.finite(m) || m >= n) {
    if (constrained) {
      fit <- fit_constrained_weights(A, linear, x, h, constraint,
                                     diag(n), rep(1, n), constraint_range,
                                     constraint_grid_n, kernel = kernel)
      return(list(weights = fit$q, value = fit$value, breaks = integer(0),
                  group_of = seq_len(n), n_groups = n, mode = fit$mode))
    }
    fit <- solve_simplex_qp(A, linear / n, rep(1, n))
    return(list(weights = fit$q, value = fit$value,
                breaks = integer(0), group_of = seq_len(n), n_groups = n))
  }
  m <- max(1L, as.integer(round(m)))

  ## Searching breakpoints inside a constrained fit would mean re-solving a
  ## constrained program for every candidate, which is not worth the cost; fall
  ## back to equal blocks and say so.
  if (constrained && identical(breaks, "optimal")) breaks <- "equal"

  if (is.numeric(breaks)) {
    blocks <- blocks_from_breaks(x, breaks)
  } else {
    breaks <- match.arg(breaks, c("optimal", "equal", "modal"))
    if (breaks == "optimal") {
      found <- search_breaks(A, linear, x, m, min_block)
      blocks <- blocks_from_breaks(x, found$breaks)
    } else {
      chosen <- if (breaks == "equal") {
        enforce_min_block(equal_breaks(n, m), n, min_block, m)
      } else {
        modal_breaks(x, m, min_block)
      }
      blocks <- blocks_from_breaks(x, chosen)
    }
  }

  G <- blocks$G
  if (constrained) {
    fit <- fit_constrained_weights(A, linear, x, h, constraint, G,
                                   blocks$group_size, constraint_range,
                                   constraint_grid_n, kernel = kernel)
    return(list(weights = as.vector(G %*% fit$q), value = fit$value,
                breaks = blocks$breaks, group_of = blocks$group_of,
                n_groups = ncol(G), mode = fit$mode))
  }

  fit <- solve_simplex_qp(
    H = t(G) %*% A %*% G,
    f = as.vector(t(G) %*% linear) / n,
    s = blocks$group_size
  )

  list(weights = as.vector(G %*% fit$q), value = fit$value,
       breaks = blocks$breaks, group_of = blocks$group_of,
       n_groups = ncol(G))
}

#' Default minimum block size
#' @keywords internal
#' @noRd
default_min_block <- function(n, m) {
  if (is.null(m) || !is.finite(m)) return(1L)
  max(2L, as.integer(floor(0.02 * n)))
}

#' Block structure for data sharpening
#'
#' Data sharpening has no closed-form criterion to search breakpoints against,
#' so only the fixed strategies are offered here.
#'
#' @keywords internal
#' @noRd
sharpen_blocks <- function(x, m, breaks, min_block) {
  n <- length(x)
  if (is.null(m) || !is.finite(m) || m >= n) {
    return(list(G = diag(n), group_size = rep(1, n),
                group_of = seq_len(n), breaks = integer(0)))
  }
  m <- max(1L, as.integer(round(m)))

  chosen <- if (is.numeric(breaks)) {
    breaks
  } else if (identical(breaks, "modal")) {
    modal_breaks(x, m, min_block)
  } else {
    enforce_min_block(equal_breaks(n, m), n, min_block, m)
  }
  blocks_from_breaks(x, chosen)
}
