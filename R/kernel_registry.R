#' Kernels for the perturbed estimator
#'
#' The estimator \eqn{\hat f(x) = \sum_i p_i K_h(x - x_i)} needs a kernel `K`.
#' The papers use the standard normal throughout, noting in Section 4.1 of
#' Doosti and Hall (2016) that it "simplified numerical work"; the theory in
#' Section 3.1 is in fact stated for kernels satisfying condition (3.1), which
#' holds when `K` is a k-fold convolution of a Laplace density. Both families
#' are available here, along with the usual compactly supported kernels.
#'
#' @section Available kernels:
#' \describe{
#'   \item{`"gaussian"`}{The standard normal density. The default, and the one
#'     used for every numerical result in both papers.}
#'   \item{`"laplace"`}{\eqn{\tfrac12 e^{-|u|}}. The k = 1 case of condition
#'     (3.1) of the 2016 paper.}
#'   \item{`"laplace2"`}{\eqn{\tfrac14 (1+|u|)e^{-|u|}}, the twofold convolution
#'     of a Laplace density and the k = 2 case of condition (3.1). Given as an
#'     example in the paper.}
#'   \item{`"laplace3"`, `"laplace4"`}{Higher-order members of the same family,
#'     progressively smoother.}
#'   \item{`"epanechnikov"`}{\eqn{\tfrac34 (1-u^2)} on \eqn{[-1,1]}. The
#'     minimum-variance choice under the usual asymptotics.}
#'   \item{`"biweight"`}{\eqn{\tfrac{15}{16}(1-u^2)^2} on \eqn{[-1,1]}.}
#'   \item{`"triweight"`}{\eqn{\tfrac{35}{32}(1-u^2)^3} on \eqn{[-1,1]}.}
#'   \item{`"triangular"`}{\eqn{1-|u|} on \eqn{[-1,1]}.}
#' }
#'
#' @section Why any kernel keeps the problem convex:
#' The quadratic form is \eqn{A_{ij} = (K\ast K)(x_i - x_j)}, whose Fourier
#' transform is \eqn{\phi_K^2 \ge 0}. A function with a non-negative Fourier
#' transform has a positive semidefinite Gram matrix, so `A` is positive
#' semidefinite whatever `K` is, and weight selection remains a convex quadratic
#' program. Equivalently, \eqn{p^\top A p = \int \hat f^2 \ge 0} by
#' construction.
#'
#' @section A caution about bandwidths:
#' Bandwidths are on each kernel's own scale, so they are not comparable across
#' kernels: `bw = 0.5` means something different for the Gaussian, which has
#' unit variance, than for the Epanechnikov, which is supported on
#' \eqn{[-1,1]}. If you switch kernels, let the bandwidth be chosen afresh
#' rather than carrying a number over.
#'
#' @param name Kernel name, or a list describing a custom kernel (see Details).
#'
#' @return `tilt_kernel()` returns a list with components `name`, `dens` (the
#'   kernel), `ft` (its Fourier transform), `self_conv` (the self-convolution at
#'   bandwidth 1) and `support` (`Inf` for the unbounded kernels).
#'   `tilt_kernels()` returns the available names.
#'
#' @details
#' A custom kernel may be supplied as a list with elements `dens`, `ft` and
#' `support`. `dens` and `ft` must be vectorised, `ft` must satisfy
#' `ft(0) == 1`, and `support` is the half-width of the support or `Inf`. The
#' self-convolution is then computed numerically.
#'
#' @examples
#' tilt_kernels()
#'
#' set.seed(1)
#' x <- c(rnorm(50, -1.5), rnorm(50, 1.5))
#'
#' ## The Laplace-convolution kernel of the paper's condition (3.1).
#' fit <- tilt_density(x, m = 3, kernel = "laplace2")
#' fit
#'
#' @name kernels
#' @export
tilt_kernels <- function() {
  c("gaussian", "laplace", "laplace2", "laplace3", "laplace4",
    "epanechnikov", "biweight", "triweight", "triangular")
}

#' @rdname kernels
#' @export
tilt_kernel <- function(name = "gaussian") {
  if (inherits(name, "tilt_kernel")) return(name)
  if (is.list(name)) return(custom_kernel(name))

  name <- tolower(as.character(name)[1L])
  if (name == "normal") name <- "gaussian"

  kern <- switch(
    name,
    gaussian = list(
      dens = function(u) exp(-u^2 / 2) / sqrt(2 * pi),
      ft   = function(t) exp(-t^2 / 2),
      ## Gaussian convolved with itself is Gaussian with variance 2.
      self_conv = function(v) exp(-v^2 / 4) / (2 * sqrt(pi)),
      moment2 = 1,
      support = Inf
    ),
    laplace  = laplace_kernel(1L),
    laplace2 = laplace_kernel(2L),
    laplace3 = laplace_kernel(3L),
    laplace4 = laplace_kernel(4L),
    epanechnikov = compact_kernel(
      function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0),
      function(t) {
        out <- 3 * (sin(t) - t * cos(t)) / t^3
        out[abs(t) < 1e-4] <- 1
        out
      }, moment2 = 1 / 5),
    biweight = compact_kernel(
      function(u) ifelse(abs(u) <= 1, (15 / 16) * (1 - u^2)^2, 0),
      function(t) {
        out <- 15 * (3 * sin(t) - 3 * t * cos(t) - t^2 * sin(t)) / t^5
        out[abs(t) < 1e-2] <- 1 - t[abs(t) < 1e-2]^2 / 14
        out
      }, moment2 = 1 / 7),
    triweight = compact_kernel(
      function(u) ifelse(abs(u) <= 1, (35 / 32) * (1 - u^2)^3, 0),
      NULL, moment2 = 1 / 9),
    triangular = compact_kernel(
      function(u) ifelse(abs(u) <= 1, 1 - abs(u), 0),
      function(t) {
        out <- 2 * (1 - cos(t)) / t^2
        out[abs(t) < 1e-4] <- 1
        out
      }, moment2 = 1 / 6),
    stop("Unknown kernel \"", name, "\". Available: ",
         paste(tilt_kernels(), collapse = ", "), ".", call. = FALSE)
  )

  kern$name <- name
  kern <- finalise_kernel(kern)
  kern$validated <- TRUE
  kern
}

#' k-fold convolution of a Laplace density
#'
#' The density with Fourier transform \eqn{(1+t^2)^{-k}} is
#' \deqn{f_k(u) = \frac{e^{-|u|}}{2^k (k-1)!}
#'       \sum_{j=0}^{k-1} \frac{(k-1+j)!}{j!\,(k-1-j)!}\frac{|u|^{k-1-j}}{2^j},}
#' which for k = 1 and k = 2 gives the two examples in Section 3.1 of Doosti and
#' Hall (2016). Its self-convolution is \eqn{f_{2k}}, since squaring the Fourier
#' transform doubles the order, so no numerical convolution is needed.
#'
#' @keywords internal
#' @noRd
laplace_kernel <- function(k) {
  dens_k <- function(order) {
    function(u) {
      a <- abs(u)
      j <- 0:(order - 1L)
      coef <- factorial(order - 1L + j) /
        (factorial(j) * factorial(order - 1L - j) * 2^j)
      powers <- outer(a, order - 1L - j, "^")
      exp(-a) * as.vector(powers %*% coef) / (2^order * factorial(order - 1L))
    }
  }
  list(
    dens = dens_k(k),
    ft = function(t) (1 + t^2)^(-k),
    self_conv = dens_k(2L * k),
    ## A Laplace density has variance 2, so the k-fold convolution has 2k.
    moment2 = 2 * k,
    support = Inf
  )
}

#' Compactly supported kernel with a numerically computed self-convolution
#' @keywords internal
#' @noRd
compact_kernel <- function(dens, ft, moment2 = NULL) {
  if (is.null(ft)) ft <- numeric_ft(dens, 1)
  if (is.null(moment2)) moment2 <- numeric_moment2(dens, 1)
  list(dens = dens, ft = ft, self_conv = numeric_self_conv(dens, 1),
       moment2 = moment2, support = 1)
}

#' Self-convolution computed once on a grid and interpolated
#'
#' \eqn{(K\ast K)(v) = \int K(w) K(v-w)\,dw}. For a compactly supported kernel
#' this is a finite integral, so it can be evaluated accurately on a grid and
#' splined. Doing it once at construction keeps it out of the inner loops.
#'
#' @keywords internal
#' @noRd
numeric_self_conv <- function(dens, support) {
  ## The self-convolution of a compactly supported kernel is a piecewise
  ## polynomial whose pieces meet at multiples of the support half-width, and a
  ## single spline laid across those joints smooths the kinks away. Fitting one
  ## spline per piece keeps the interpolation exact to rounding error.
  knots <- support * c(-2, -1, 0, 1, 2)

  conv_at <- function(v) {
    lo <- max(-support, v - support)
    hi <- min(support, v + support)
    if (hi <= lo) return(0)

    ## The integrand K(w) K(v - w) can be non-smooth where either factor is:
    ## at w = 0 and w = v for a kernel with a corner at the origin, such as the
    ## triangular. Gauss-Legendre across a corner loses most of its accuracy, so
    ## the range is split there before integrating.
    cuts <- sort(unique(c(lo, hi, c(0, v)[c(0, v) > lo & c(0, v) < hi])))
    total <- 0
    for (k in seq_len(length(cuts) - 1L)) {
      rule <- gauss_legendre(cuts[k], cuts[k + 1L], 20L, 8L)
      total <- total + sum(rule$weights * dens(rule$nodes) * dens(v - rule$nodes))
    }
    total
  }

  pieces <- lapply(seq_len(4L), function(k) {
    grid <- seq(knots[k], knots[k + 1L], length.out = 401L)
    stats::splinefun(grid, vapply(grid, conv_at, numeric(1)), method = "fmm")
  })

  function(v) {
    out <- numeric(length(v))
    for (k in seq_len(4L)) {
      inside <- v >= knots[k] & v <= knots[k + 1L]
      if (any(inside)) out[inside] <- pieces[[k]](v[inside])
    }
    out[abs(v) >= 2 * support] <- 0
    pmax(out, 0)
  }
}

#' Fourier transform computed by quadrature
#' @keywords internal
#' @noRd
numeric_ft <- function(dens, support) {
  function(t) {
    rule <- gauss_legendre(-support, support, 200L, 8L)
    w <- rule$weights * dens(rule$nodes)
    as.vector(cos(outer(as.vector(t), rule$nodes)) %*% w)
  }
}

#' Validate a user-supplied kernel
#' @keywords internal
#' @noRd
custom_kernel <- function(spec) {
  required <- c("dens", "ft")
  missing <- setdiff(required, names(spec))
  if (length(missing)) {
    stop("A custom kernel needs components: ", paste(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  if (abs(spec$ft(0) - 1) > 1e-6) {
    stop("`ft(0)` must equal 1 for a kernel that integrates to one.",
         call. = FALSE)
  }
  support <- if (is.null(spec$support)) Inf else spec$support
  if (is.null(spec$self_conv)) {
    if (!is.finite(support)) {
      stop("A custom kernel with unbounded support must supply `self_conv`.",
           call. = FALSE)
    }
    spec$self_conv <- numeric_self_conv(spec$dens, support)
  }
  spec$support <- support
  if (is.null(spec$moment2)) {
    spec$moment2 <- numeric_moment2(spec$dens, if (is.finite(support)) support else 50)
  }
  spec$name <- if (is.null(spec$name)) "custom" else spec$name
  finalise_kernel(spec)
}

#' Attach the canonical bandwidth factor
#'
#' The canonical factor of Marron and Nolan (1988) is
#' \deqn{\delta_K = \left(R(K)/\mu_2(K)^2\right)^{1/5},
#'       \quad R(K)=\int K^2, \quad \mu_2(K)=\int u^2 K(u)\,du.}
#' Rescaling a kernel to \eqn{K(u/\delta_K)/\delta_K} equalises the
#' asymptotically optimal bandwidth across kernels, so that the same numerical
#' bandwidth produces the same amount of smoothing whichever kernel is used.
#'
#' \eqn{R(K)} is read off the self-convolution: for a symmetric kernel
#' \eqn{(K\ast K)(0) = \int K(u)K(-u)\,du = \int K^2}.
#'
#' @keywords internal
#' @noRd
finalise_kernel <- function(kern) {
  kern$roughness <- kern$self_conv(0)
  kern$canonical <- (kern$roughness / kern$moment2^2)^(1 / 5)
  structure(kern, class = "tilt_kernel")
}

#' Second moment by quadrature
#' @keywords internal
#' @noRd
numeric_moment2 <- function(dens, support) {
  rule <- gauss_legendre(-support, support, 400L, 8L)
  sum(rule$weights * rule$nodes^2 * dens(rule$nodes))
}

#' @export
print.tilt_kernel <- function(x, ...) {
  cat("<tilt_kernel> ", x$name,
      if (is.finite(x$support)) paste0(", support [-", x$support, ", ",
                                       x$support, "]") else ", unbounded support",
      "\n", sep = "")
  cat("  roughness R(K) = ", format(x$roughness, digits = 5),
      ",  second moment = ", format(x$moment2, digits = 5),
      ",  canonical factor = ", format(x$canonical, digits = 5), "\n", sep = "")
  invisible(x)
}

#' Construct and assess a custom kernel
#'
#' `make_kernel()` builds a kernel object from a density function alone,
#' computing everything else the package needs -- the Fourier transform, the
#' self-convolution, the second moment and the canonical bandwidth factor -- by
#' quadrature. By default it then checks that the result is suitable and stops
#' with an explanation if it is not. `check_kernel()` runs the same checks on
#' any kernel, built-in or custom, and reports them; `summary()` on a kernel
#' object returns the same information as an object.
#'
#' A kernel is suitable for tilting if it is a symmetric probability density.
#' Non-negativity is what guarantees that the tilted estimate is a density;
#' symmetry is assumed by the Fourier representation of the cross terms; and
#' unit mass is what keeps the weights on the simplex meaningful.
#'
#' The estimators [tilt_density()], [tilt_density_cv()] and [sharpen_density()]
#' refuse a kernel that fails these checks, so an unsuitable kernel built with
#' `check = FALSE` can be inspected but not used.
#'
#' @param dens A vectorised function giving the kernel density.
#' @param support Half-width of the support, or `Inf`. For an unbounded kernel
#'   the numerical integrals are truncated at `support_eff` (default 40), which
#'   is ample for any kernel with tails at least as light as the Laplace.
#' @param name A name for the kernel.
#' @param ft,self_conv Optional closed forms for the Fourier transform and the
#'   self-convolution; computed numerically when omitted.
#' @param support_eff Truncation point used for unbounded kernels.
#' @param check If `TRUE` (the default), stop with an informative error when
#'   the kernel is not a symmetric probability density. Set to `FALSE` to
#'   construct an unsuitable kernel for inspection.
#' @param kernel A kernel name or object, for `check_kernel()`.
#' @param tolerance Tolerance for the mass, symmetry and positivity checks.
#'
#' @return `make_kernel()` returns a `"tilt_kernel"` object usable wherever a
#'   kernel name is accepted. `check_kernel()` prints a report and returns,
#'   invisibly, the `"summary.tilt_kernel"` object described in
#'   [summary.tilt_kernel()].
#'
#' @examples
#' ## The logistic kernel, from its density alone.
#' logistic <- make_kernel(function(u) exp(-u) / (1 + exp(-u))^2,
#'                         support = Inf, name = "logistic")
#' logistic
#' summary(logistic)
#'
#' set.seed(1)
#' x <- c(rnorm(40, -1.5), rnorm(40, 1.5))
#' tilt_density(x, m = 3, breaks = "equal", kernel = logistic)$distance2
#'
#' ## A function that is not a density is refused at construction ...
#' try(make_kernel(function(u) ifelse(abs(u) <= 1, 1 - u^2, 0), support = 1))
#'
#' ## ... and if built anyway for inspection, cannot be used for estimation.
#' bad <- make_kernel(function(u) ifelse(abs(u) <= 1, 1 - u^2, 0),
#'                    support = 1, check = FALSE)
#' check_kernel(bad)
#' try(tilt_density(x, kernel = bad))
#'
#' @name make_kernel
#' @export
make_kernel <- function(dens, support = Inf, name = "custom", ft = NULL,
                        self_conv = NULL, support_eff = 40, check = TRUE) {
  if (!is.function(dens)) stop("`dens` must be a function.", call. = FALSE)
  lim <- if (is.finite(support)) support else support_eff

  if (is.null(ft)) ft <- numeric_ft(dens, lim)
  if (is.null(self_conv)) {
    self_conv <- if (is.finite(support)) numeric_self_conv(dens, support)
                 else numeric_self_conv_unbounded(dens, lim)
  }
  kern <- finalise_kernel(list(name = name, dens = dens, ft = ft,
                               self_conv = self_conv,
                               moment2 = numeric_moment2(dens, lim),
                               support = support))
  if (check) assert_kernel(kern, context = "make_kernel()")
  kern
}

#' @rdname make_kernel
#' @export
check_kernel <- function(kernel, tolerance = 1e-3) {
  out <- summary(tilt_kernel(kernel), tolerance = tolerance)
  print(out)
  invisible(out)
}

#' Summarise a kernel
#'
#' Reports the properties the estimators rely on -- unit mass, non-negativity,
#' symmetry, and \eqn{\phi_K(0) = 1} -- together with the second moment, the
#' roughness \eqn{R(K) = \int K^2} and the canonical bandwidth factor of Marron
#' and Nolan (1988).
#'
#' @param object A `"tilt_kernel"` object, from [tilt_kernel()] or
#'   [make_kernel()].
#' @param tolerance Tolerance for the mass, symmetry and positivity checks.
#' @param ... Ignored.
#'
#' @return An object of class `"summary.tilt_kernel"`: a list with the kernel
#'   `name` and `support`; logical elements `integrates_to_one`,
#'   `non_negative`, `symmetric`, `ft_at_zero_is_one` and their conjunction
#'   `suitable`; and numeric elements `mass`, `moment2`, `roughness` and
#'   `canonical`. It has a `print` method.
#'
#' @examples
#' summary(tilt_kernel("epanechnikov"))
#' summary(tilt_kernel("laplace2"))$canonical
#'
#' @export
summary.tilt_kernel <- function(object, tolerance = 1e-3, ...) {
  k   <- object
  lim <- if (is.finite(k$support)) k$support else 40
  rule <- gauss_legendre(-lim, lim, 400L, 8L)
  y    <- k$dens(rule$nodes)
  mass <- sum(rule$weights * y)
  asym <- max(abs(y - k$dens(-rule$nodes)))

  out <- list(
    name              = k$name,
    support           = k$support,
    integrates_to_one = abs(mass - 1) < tolerance,
    non_negative      = all(y >= -tolerance),
    symmetric         = asym < tolerance,
    ft_at_zero_is_one = abs(k$ft(0) - 1) < tolerance,
    mass              = mass,
    moment2           = k$moment2,
    roughness         = k$roughness,
    canonical         = k$canonical
  )
  out$suitable <- out$integrates_to_one && out$non_negative &&
                  out$symmetric && out$ft_at_zero_is_one
  structure(out, class = "summary.tilt_kernel")
}

#' @rdname summary.tilt_kernel
#' @param x A `"summary.tilt_kernel"` object.
#' @param digits Number of significant digits.
#' @export
print.summary.tilt_kernel <- function(x, digits = 5L, ...) {
  yn <- function(v) if (v) "yes" else "NO "
  cat("Kernel:", x$name,
      if (is.finite(x$support)) paste0("(support [-", x$support, ", ", x$support, "])")
      else "(unbounded support)", "\n")
  cat(sprintf("  integrates to one  %s  (mass = %.6f)\n", yn(x$integrates_to_one), x$mass))
  cat(sprintf("  non-negative       %s\n", yn(x$non_negative)))
  cat(sprintf("  symmetric          %s\n", yn(x$symmetric)))
  cat(sprintf("  Fourier transform  %s  (phi(0) = 1)\n", yn(x$ft_at_zero_is_one)))
  cat("  second moment     ", format(x$moment2, digits = digits), "\n")
  cat("  roughness R(K)    ", format(x$roughness, digits = digits), "\n")
  cat("  canonical factor  ", format(x$canonical, digits = digits), "\n")
  cat(if (x$suitable) "Suitable for tilting.\n" else
      "NOT suitable: the tilted estimate would not be a proper density.\n")
  invisible(x)
}

#' Plot a kernel
#'
#' @param x A `"tilt_kernel"` object.
#' @param from,to Range of the plot. Defaults to the support, or to four
#'   standard deviations for an unbounded kernel.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' plot(tilt_kernel("biweight"))
#'
#' @export
plot.tilt_kernel <- function(x, from = NULL, to = NULL, ...) {
  half <- if (is.finite(x$support)) 1.1 * x$support else 4 * sqrt(x$moment2)
  if (is.null(from)) from <- -half
  if (is.null(to))   to   <- half
  u <- seq(from, to, length.out = 401)
  graphics::plot(u, x$dens(u), type = "l", xlab = "u", ylab = "K(u)",
                 main = paste("Kernel:", x$name), ...)
  graphics::abline(h = 0, col = "grey70")
  invisible(x)
}

#' Stop if a kernel is not a symmetric probability density
#' @keywords internal
#' @noRd
assert_kernel <- function(kern, context) {
  if (isTRUE(kern$validated)) return(invisible(kern))
  s <- summary(kern)
  if (s$suitable) return(invisible(kern))
  failed <- c(if (!s$integrates_to_one) sprintf("integrates to %.4g, not 1", s$mass),
              if (!s$non_negative) "takes negative values",
              if (!s$symmetric) "is not symmetric",
              if (!s$ft_at_zero_is_one && s$integrates_to_one) "has phi(0) != 1")
  stop(context, ": kernel \"", kern$name, "\" is not a symmetric probability ",
       "density -- it ", paste(failed, collapse = ", "), ". ",
       "See check_kernel().", call. = FALSE)
}

#' Self-convolution of an unbounded kernel by quadrature and interpolation
#' @keywords internal
#' @noRd
numeric_self_conv_unbounded <- function(dens, lim) {
  grid <- seq(-2 * lim, 2 * lim, length.out = 2001L)
  rule <- gauss_legendre(-lim, lim, 200L, 8L)
  vals <- vapply(grid, function(v) {
    sum(rule$weights * dens(rule$nodes) * dens(v - rule$nodes))
  }, numeric(1))
  interp <- stats::splinefun(grid, vals, method = "natural")
  function(v) pmax(interp(v), 0)
}
