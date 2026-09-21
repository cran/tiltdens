#' Multivariate tilted density estimation
#'
#' Back end for [tilt_density()] and [tilt_density_cv()] when the sample is a
#' matrix. Only `m = Inf` (one weight per observation) is supported: the block
#' structure of the univariate case rests on the order statistics, which have no
#' natural counterpart in more than one dimension.
#'
#' Data sharpening is deliberately not offered here. The search would range over
#' `n * l` variables and the run times reported in the paper for the univariate
#' case, minutes per sample of size 100, make that impractical rather than
#' merely slow.
#'
#' @keywords internal
#' @noRd
tilt_density_md <- function(x, comparator, bw, comparator_bw, bw_grid,
                            data_name, call, grid_n, criterion = "comparator",
                            max_groups = 1L, rho = 0.8, kernel = NULL) {
  x <- check_sample_md(x)
  n <- nrow(x)
  l <- ncol(x)

  if (criterion == "comparator") {
    if (is.null(comparator_bw)) comparator_bw <- bw_comparator_cv_md(x, comparator)
    check_bw(comparator_bw)
    h_candidates <- if (!is.null(bw_grid)) as.vector(bw_grid) else
      if (!is.null(bw)) bw else comparator_bw
    for (h in h_candidates) check_bw(h)

    const_term <- comparator_normsq_md(x, comparator_bw, comparator)$value

    values  <- numeric(length(h_candidates))
    weights <- matrix(0, n, length(h_candidates))
    for (k in seq_along(h_candidates)) {
      h <- h_candidates[k]
      A <- self_conv_gram_md(x, h, kernel)
      b <- comparator_cross_md(x, x, h, comparator_bw, comparator, kernel)
      fit <- solve_simplex_qp(A, b / n, rep(1, n))
      weights[, k] <- fit$q
      values[k]    <- fit$value + const_term
    }
    best <- which.min(values)
    extra <- list(distance2 = values[best], comparator = comparator,
                  comparator_bw = comparator_bw)
    chosen_bw <- h_candidates[best]
    chosen_p  <- weights[, best]

  } else {
    ## Cross-validation criterion, the multivariate form of tilt_cv().
    if (is.null(bw_grid)) {
      scale <- exp(mean(log(apply(x, 2L, robust_scale))))
      c1 <- 1 / 30
      bw_grid <- scale *
        exp(seq(log(n^(-1 / (l + 4) - c1)), log(n^(-1 / (l + 4) + c1)),
                length.out = 25))
    }
    bw_grid <- as.vector(bw_grid)

    uniform <- rep(1 / n, n)
    cv_grid <- vapply(bw_grid, function(h) {
      as.numeric(tilt_cv_md(x, h, uniform, kernel))
    }, numeric(1))
    ix <- which.min(cv_grid)

    bw_uniform <- bw_grid[ix]
    cv_uniform <- cv_grid[ix]

    ## With every weight free there is nothing to gain from the sequential
    ## block scheme, so a single pass over the bandwidth grid is the whole job.
    best_cv <- Inf
    chosen_bw <- bw_uniform
    chosen_p  <- uniform
    for (h in bw_grid[bw_grid >= rho * bw_uniform]) {
      parts <- tilt_cv_md(x, h, uniform, kernel)
      fit <- solve_simplex_qp(attr(parts, "A"), attr(parts, "C") / n, rep(1, n))
      if (fit$value < best_cv) {
        best_cv   <- fit$value
        chosen_bw <- h
        chosen_p  <- fit$q
      }
    }
    extra <- list(cv = best_cv,
                  trace = data.frame(n_groups = c(1L, n),
                                     bw = c(bw_uniform, chosen_bw),
                                     cv = c(cv_uniform, best_cv)))
  }

  structure(
    c(list(
      data = x, weights = chosen_p, bw = chosen_bw, n = n, d = l,
      call = call, data.name = data_name,
      method = if (criterion == "comparator") "tilt_density" else "tilt_density_cv",
      n_groups = n, group_of = seq_len(n), breaks = integer(0),
      breaks_method = "none", kernel = kernel
    ), extra),
    class = c("tiltdens_md", "tiltdens")
  )
}

#' Multivariate cross-validation criterion
#' @keywords internal
#' @noRd
tilt_cv_md <- function(x, bw, weights, kernel = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)

  A <- self_conv_gram_md(x, bw, kernel)
  M <- kernel_matrix_md(x, x, bw, kernel)
  diag(M) <- 0
  cvec <- colSums(M)

  value <- as.numeric(t(weights) %*% A %*% weights) - 2 * sum(cvec * weights) / n
  structure(value, A = A, C = cvec)
}

#' @rdname predict.tiltdens
#' @export
predict.tiltdens_md <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(as.vector(object$fitted))
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != object$d) {
    stop("`newdata` must have ", object$d, " columns.", call. = FALSE)
  }
  as.vector(kernel_matrix_md(newdata, object$data, object$bw,
                             object$kernel) %*% object$weights)
}

#' Print a multivariate fit
#'
#' @param x A `"tiltdens_md"` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.tiltdens_md <- function(x, digits = getOption("digits") - 2L, ...) {
  label <- if (identical(x$method, "tilt_density")) {
    "Tilted density estimate (Doosti & Hall 2016, Section 3.3)"
  } else {
    "Tilted density estimate, cross-validated (Doosti, Hall & Mateu 2018)"
  }
  cat("\n", label, "\n\n", sep = "")
  cat("Call:      "); print(x$call)
  cat("Data:      ", x$data.name, " (", x$n, " obs., ", x$d, " dimensions)\n", sep = "")
  cat("Bandwidth: ", format(x$bw, digits = digits), " (product kernel)\n", sep = "")
  if (!is.null(x$comparator)) {
    cat("Comparator:", x$comparator, "with bandwidth",
        format(x$comparator_bw, digits = digits), "\n")
  }
  if (!is.null(x$distance2)) {
    cat("Distance:  ", format(x$distance2, digits = digits),
        " (squared L2 to the comparator)\n", sep = "")
  }
  if (!is.null(x$cv)) cat("CV:        ", format(x$cv, digits = digits), "\n", sep = "")
  cat("\nUse predict() to evaluate, or plot() for a two-dimensional contour.\n\n")
  invisible(x)
}

#' Contour plot of a bivariate fit
#'
#' @param x A `"tiltdens_md"` object with two dimensions.
#' @param n Grid points per axis.
#' @param points Draw the observations over the contours.
#' @param ... Passed to [graphics::contour()].
#'
#' @return `x`, invisibly.
#'
#' @export
plot.tiltdens_md <- function(x, n = 60, points = TRUE, ...) {
  if (x$d != 2L) {
    stop("Only two-dimensional fits can be plotted; use predict() instead.",
         call. = FALSE)
  }
  rng <- apply(x$data, 2L, range)
  pad <- 3 * x$bw
  gx <- seq(rng[1L, 1L] - pad, rng[2L, 1L] + pad, length.out = n)
  gy <- seq(rng[1L, 2L] - pad, rng[2L, 2L] + pad, length.out = n)

  grid <- as.matrix(expand.grid(gx, gy))
  z <- matrix(predict(x, grid), n, n)

  graphics::contour(gx, gy, z,
                    xlab = colnames(x$data)[1L], ylab = colnames(x$data)[2L], ...)
  if (points) graphics::points(x$data, pch = 16, cex = 0.5, col = "grey40")
  invisible(x)
}
