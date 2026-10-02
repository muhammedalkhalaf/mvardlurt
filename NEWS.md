# mvardlurt 1.1.0

This release corrects the bootstrap and the interpretation of the test so
that they follow Sam, McNown, Goh and Goh (2025), Studies in Economics and
Econometrics 49(1), 17-33, doi:10.1080/03796205.2024.2439101. The version on
CRAN is 1.0.2.

## Breaking changes

* `level` is now the significance level of the two tests (default 0.05, the
  level used in the paper's simulations); in 1.0.2 it was a confidence level
  (default 0.95) that was stored but never used, the decision always being
  made at 10 percent. Values of 0.5 or more are rejected with a message
  pointing to `1 - level`.
* The lag convention follows equation (1) of the paper: an ARDL(p, q) model
  has p - 1 lagged differences of y and q - 1 lagged differences of x
  (p, q at least 1). `fixlag = c(p, q)` in 1.0.2 corresponds to
  `fixlag = c(p + 1, q + 1)` in 1.1.0 (plus the contemporaneous difference
  of x, see below). The paper itself is inconsistent: Section 4.1 and
  Table 4 count lagged differences directly (the 1.0.2 convention), so the
  paper's ARDL(0, 2) is `fixlag = c(1, 3)` here. `maxlag` must be between
  1 and 12.
* `reps` was renamed to `nboot` (default 999, minimum 99); `reps` still works
  with a deprecation warning. `pi_coef`, `pi_se`, `delta_coef` and
  `delta_se` were renamed to `b1`, `b1_se`, `beta2` and `beta2_se`;
  `fstat_p` was removed (see below).

## Corrections

* Bootstrap (Section 4.2 of the paper). Version 1.0.2 simulated y* and x* as
  independent Gaussian random walks without any estimated coefficient and
  without resampling residuals, so the critical values did not depend on the
  data (Monte Carlo size of the t test with drift was 0.010 at the nominal
  5 percent level). The t and F statistics are now bootstrapped separately
  with the respective null imposed: the restricted regression (without
  y[t-1] for the t test, without x[t-1] for the F test) is estimated, its
  recentred residuals are resampled with replacement, y* is generated
  recursively from the restricted estimated equation with the observed x
  held fixed, and the unrestricted regression is re-estimated on (y*, x).
  The contemporaneous difference of x, which equation (24) of the paper
  omits although equation (1) has it, is included in the restricted
  regressions and in the bootstrap, following equation (1). Bootstrap
  critical values at 10, 5, 2.5 and 1 percent and bootstrap p-values
  (`t_pval`, `f_pval`) are reported; the bootstrap distributions are
  returned in `t_boot` and `f_boot`.
* Four cases (Section 3.2 of the paper). The labels were wrong in 1.0.2. They
  are now: Case I, neither test rejects, y is I(1) with no cointegration;
  Case II, only the t test rejects, y is I(0); Case III, only the F test
  rejects, y is I(2); Case IV, both reject, y is I(1) through cointegration
  with x. The order of integration is stated in the output
  (`decision$integration`). The wording "degenerate case", "spurious" and
  "H0: delta = 0, no cointegration" was removed; the F test is for
  H0: beta2 = 0 (equation (12)).
* Test regression (equation (1) of the paper). The contemporaneous difference
  of x is now included in the regression, in both restricted regressions and
  in the bootstrap (equation (1); see Breaking changes for the lag
  convention).
* `x` may be a matrix with several covariates; the F statistic is then the
  joint Wald F test on all lagged levels. `beta2`, `beta2_se` and `lr_mult`
  have one element per covariate, and `coef()` returns `b1`,
  `beta2.<name>` and `lr_mult.<name>`.
* The decision uses `level` (see Breaking changes). The significance codes
  for 10, 5, 2.5 and 1 percent are still reported.
* `seed` now defaults to `NULL`. When a seed is supplied the caller's RNG
  state is saved and restored on exit, so a call no longer resets the global
  random number stream (1.0.2 called `set.seed(12345)` unconditionally).
* `fstat_p`, a p-value from the standard F distribution that is invalid here
  (Section 3.3 of the paper), was removed in favour of the bootstrap
  p-value.
* Lag selection compares all candidate models on a common sample
  (t = maxlag + 1, ..., n); 1.0.2 compared AIC or BIC values computed on
  different samples. The IC table is indexed by p = 1, ..., maxlag and
  q = 1, ..., maxlag.
* The citation of the paper was corrected to volume 49(1), pages 17-33,
  2025 (it read 2024, 1-17).
* The data frame of the test regression is returned as `data`; the `model`
  component is for inspection, and can be refitted with
  `lm(formula(object$model), data = object$data)`.

# mvardlurt 1.0.0

* Initial CRAN release.
* Implements the multivariate ARDL unit root test (Sam, McNown, Goh and Goh, 2025) [citation corrected in 1.1.0].
* Features:
  - Automatic lag selection via AIC/BIC
  - Bootstrap critical values for correct size
  - Three deterministic cases (none, intercept, intercept + trend)
  - Four-case decision framework for inference
  - Comprehensive print and summary methods
  - Diagnostic plots
