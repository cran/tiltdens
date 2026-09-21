#' Validate a sample vector
#' @keywords internal
#' @noRd
check_sample <- function(x, min_n = 3L) {
  if (!is.numeric(x)) stop("`x` must be a numeric vector.", call. = FALSE)
  x <- as.vector(x)
  if (anyNA(x)) {
    x <- x[!is.na(x)]
    warning("Missing values in `x` were removed.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("`x` must contain only finite values.", call. = FALSE)
  }
  if (length(x) < min_n) {
    stop("`x` must have at least ", min_n, " observations.", call. = FALSE)
  }
  x
}

#' Validate a bandwidth
#' @keywords internal
#' @noRd
check_bw <- function(bw) {
  if (!is.numeric(bw) || length(bw) != 1L || !is.finite(bw) || bw <= 0) {
    stop("Bandwidth must be a single positive finite number.", call. = FALSE)
  }
  invisible(bw)
}

#' Normalise a comparator name
#' @keywords internal
#' @noRd
match_comparator <- function(comparator) {
  comparator <- tolower(as.character(comparator)[1L])
  if (comparator %in% c("sinc")) return("sinc")
  if (comparator %in% c("trapezoid", "trapezoidal", "flattop", "flat-top")) {
    return("trapezoid")
  }
  stop("`comparator` must be \"sinc\" or \"trapezoid\", not \"", comparator, "\".",
       call. = FALSE)
}

#' Robust scale estimate used by several bandwidth rules
#' @keywords internal
#' @noRd
robust_scale <- function(x) {
  scale <- min(mad(x, constant = 1) / 0.6745, sd(x))
  if (!is.finite(scale) || scale <= 0) scale <- sd(x)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  scale
}
