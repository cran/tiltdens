# tiltdens 0.1.1

First release.

* `tilt_density()` implements the tilted estimator of Doosti and Hall (2016),
  with weights chosen by minimising the L2 distance to a sinc or trapezoidal
  comparator estimator.
* When the number of distinct weights `m` is finite, the block boundaries can be
  chosen in three ways through the `breaks` argument: `"optimal"` selects them
  jointly with the weights, as Section 4.1 of the 2016 paper specifies;
  `"equal"` uses near-equal blocks, Algorithm A of the 2018 paper; `"modal"`
  places them at the troughs of a pilot estimate, the refinement suggested in
  Section 2 of the 2018 paper. Explicit boundaries may also be supplied.
* `tilt_density_cv()` implements method (III) of Doosti, Hall and Mateu (2018),
  choosing the bandwidth and the weights together by cross-validation.
* `sharpen_density()` implements the data-sharpened estimator of Doosti and
  Hall (2016).
* The kernel of the perturbed estimator is selectable through the `kernel`
  argument: `"gaussian"` (the default, and the one used for every numerical
  result in both papers), the Laplace-convolution family `"laplace"` through
  `"laplace4"` that Section 3.1 of the 2016 paper states its condition (3.1)
  for, and `"epanechnikov"`, `"biweight"`, `"triweight"` and `"triangular"`.
  Custom kernels may be supplied as a list. `tilt_kernels()` lists them.
* Bandwidths are rescaled across kernels using the canonical factor of Marron
  and Nolan (1988), so that the same numerical bandwidth smooths by the same
  amount whichever kernel is used. `bw_canonical_factor()` returns the factor
  and `bw_convert()` maps a bandwidth between kernels; `bw_nrd_robust()` gained
  a `kernel` argument. The estimators apply the rescaling to their default
  bandwidth, so changing the kernel changes the shape of the fit rather than
  how much it is smoothed. For the Gaussian the factor is one, so nothing about
  the papers' setting changes.
  Weight selection stays a convex quadratic program whatever the kernel,
  because the Gram matrix has Fourier transform `phi^2 >= 0` and so is positive
  semidefinite by construction.
* Multivariate data are supported by `tilt_density()` and `tilt_density_cv()`:
  pass a matrix. This is Section 3.3 of the 2016 paper. The product kernel makes
  every quantity factorise across dimensions, so the l-variate matrices are
  elementwise products of the univariate ones and the optimisation problem is
  unchanged. `plot()` draws a contour for bivariate fits. Data sharpening
  remains univariate.
* Shape constraints are available through the `constraint` argument of
  `tilt_density()` and `tilt_density_cv()`: `"unimodal"`, `"increasing"` or
  `"decreasing"`. The estimator is linear in the weights, so a shape
  requirement on a grid becomes linear inequalities and the problem stays a
  convex quadratic program. Shape-constrained tilting is due to Hall and Huang
  (2002); applying it to the cross-validation criterion follows the extension
  suggested in the discussion of the 2018 paper.
* `sinc_density()` and `trapezoid_density()` provide the infinite-order
  comparator estimators on their own.
* `bw_comparator_cv()`, `bw_flattop()` and `bw_nrd_robust()` select bandwidths.
* Fitted objects inherit from `"density"`, so base plotting works unchanged,
  and have `print()`, `summary()` and `predict()` methods.
* `bw_comparator_cv()` handles tied observations, which real data recorded to
  limited precision routinely contain.
