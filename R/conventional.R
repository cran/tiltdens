#' Conventional kernel density estimator
#'
#' The ordinary estimator with uniform weights, evaluated with the same kernel
#' code as the tilted estimators. It is provided so that comparisons within the
#' package, and the replication materials of the accompanying article, do not
#' depend on [stats::density()], whose output changed slightly between
#' \R versions 4.3 and 4.4.
#'
#' @param x Numeric vector of observations.
#' @param bw Positive bandwidth. Defaults to [bw_nrd_robust()] for the kernel.
#' @param kernel Kernel name or object; see [tilt_kernels()].
#' @param n,from,to Grid on which to evaluate.
#'
#' @return An object of class `"density"`.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' plot(conventional_density(x))
#' lines(conventional_density(x, kernel = "epanechnikov"), col = "red")
#'
#' @export
conventional_density <- function(x, bw = NULL, kernel = "gaussian",
                                 n = 512, from = NULL, to = NULL) {
  data_name <- deparse(substitute(x))
  x <- check_sample(x)
  k <- tilt_kernel(kernel)
  if (is.null(bw)) bw <- bw_nrd_robust(x, k)
  check_bw(bw)
  if (is.null(from)) from <- min(x) - 3 * bw
  if (is.null(to))   to   <- max(x) + 3 * bw
  grid <- seq(from, to, length.out = n)

  structure(
    list(x = grid,
         y = as.vector(kernel_matrix(grid, x, bw, k) %*% rep(1 / length(x), length(x))),
         bw = bw, n = length(x), call = match.call(), data.name = data_name,
         has.na = FALSE, kernel = k$name),
    class = "density"
  )
}
