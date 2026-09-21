## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.2)

## ----setup--------------------------------------------------------------------
library(tiltdens)
set.seed(2016)

## ----problem------------------------------------------------------------------
x <- c(rnorm(50, -1.5), rnorm(50, 1.5))

conventional <- density(x)
sinc <- sinc_density(x)

plot(sinc, main = "Sinc kernel estimator", ylim = range(0, sinc$y))
lines(conventional, col = "grey40", lty = 2)
abline(h = 0, col = "red")
legend("topright", c("sinc", "conventional"), col = c("black", "grey40"),
       lty = c(1, 2), bty = "n")

min(sinc$y)

## ----tilting------------------------------------------------------------------
fit <- tilt_density(x, m = 3)
fit

## ----tilting-plot-------------------------------------------------------------
plot(fit, main = "Tilted vs sinc")
lines(sinc, col = "red", lty = 2)
abline(h = 0, col = "grey")
legend("topright", c("tilted", "sinc"), col = c("black", "red"),
       lty = c(1, 2), bty = "n")

## ----m-comparison-------------------------------------------------------------
fit_n <- tilt_density(x, m = Inf)
fit_3 <- tilt_density(x, m = 3)

c(m_n = fit_n$distance2, m_3 = fit_3$distance2)
length(unique(round(fit_3$weights, 8)))

## ----breaks-------------------------------------------------------------------
c(equal   = tilt_density(x, m = 3, breaks = "equal")$distance2,
  modal   = tilt_density(x, m = 3, breaks = "modal")$distance2,
  optimal = tilt_density(x, m = 3, breaks = "optimal")$distance2)

## ----breaks-where-------------------------------------------------------------
fit_opt <- tilt_density(x, m = 3, breaks = "optimal")
fit_opt$breaks          # ranks in the sorted sample
fit_opt$break_values    # the corresponding data values

## ----breaks-manual------------------------------------------------------------
tilt_density(x, m = 3, breaks = c(30, 70))$distance2

## ----weights------------------------------------------------------------------
plot(x, fit_n$weights, xlab = "observation", ylab = "tilt weight",
     main = "Fitted weights, m = Inf")
abline(h = 1 / length(x), lty = 2, col = "grey40")

## ----cv-----------------------------------------------------------------------
fit_cv <- tilt_density_cv(x)
fit_cv
fit_cv$trace

## ----sharpen------------------------------------------------------------------
fit_s <- sharpen_density(x, m = 3,
                         control = list(max_iterations = 30, population = 25))
round(unique(fit_s$shifts), 3)

## ----bandwidths---------------------------------------------------------------
c(conventional = bw_nrd_robust(x),
  flat_top     = bw_flattop(x),
  comparator   = as.numeric(bw_comparator_cv(x)))

## ----bw-curve-----------------------------------------------------------------
cv <- attr(bw_comparator_cv(x), "cv")
plot(cv$q, cv$cv, type = "l", xlab = "frequency 1/h", ylab = "CV criterion")
abline(v = 1 / as.numeric(bw_comparator_cv(x)), col = "red", lty = 2)

## ----methods------------------------------------------------------------------
predict(fit, newdata = c(-2, 0, 2))

## Integrated squared error against a known truth
truth <- function(t) 0.5 * dnorm(t, -1.5) + 0.5 * dnorm(t, 1.5)
c(conventional = sum(diff(conventional$x) *
                     ((conventional$y - truth(conventional$x))^2)[-1]),
  tilted_3 = ise(fit_3, truth),
  tilted_cv = ise(fit_cv, truth))

## ----kernels------------------------------------------------------------------
tilt_kernels()

## ----kernel-compare-----------------------------------------------------------
sapply(c("gaussian", "laplace2", "epanechnikov", "biweight"), function(k) {
  tilt_density(x, m = 3, kernel = k)$distance2
})

## ----canonical----------------------------------------------------------------
sapply(c("gaussian", "epanechnikov", "biweight"), bw_canonical_factor)
bw_convert(0.5, from = "gaussian", to = "epanechnikov")

## ----canonical-effect---------------------------------------------------------
truth <- function(t) 0.5 * dnorm(t, -1.5) + 0.5 * dnorm(t, 1.5)
kernels <- c("gaussian", "epanechnikov", "biweight", "triangular")

rescaled <- sapply(kernels, function(k) ise(tilt_density(x, m = 3, kernel = k), truth))
h_gauss  <- tilt_density(x, m = 3, kernel = "gaussian")$bw
fixed    <- sapply(kernels, function(k)
                   ise(tilt_density(x, m = 3, kernel = k, bw = h_gauss), truth))

c(rescaled = max(rescaled) / min(rescaled),
  fixed    = max(fixed) / min(fixed))

