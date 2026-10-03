#' Multivariate ARDL Unit Root Test
#'
#' @description
#' Implements the multivariate autoregressive distributed lag (ARDL) unit root
#' test of Sam, McNown, Goh and Goh (2025). The test augments the ADF
#' regression with the lagged level, the current difference and lagged
#' differences of one or more covariates, so that a possible cointegrating
#' relationship between the series under test and the covariates is taken into
#' account. Critical values and p-values are obtained by the residual
#' bootstrap of Section 4.2 of the paper.
#'
#' @details
#' The test regression is equation (1) of the paper,
#'
#' \deqn{\Delta y_t = c_1 + c_2 t + b_1 y_{t-1} + \beta_2' x_{t-1}
#'   + \sum_{i=1}^{p-1} \phi_i \Delta y_{t-i}
#'   + \sum_{j=1}^{q-1} \Phi_j' \Delta x_{t-j} + \omega' \Delta x_t + u_t,}
#'
#' estimated by OLS on t = max(p, q) + 1, ..., n. The orders follow equation
#' (1) of the paper: an ARDL(p, q) model contains p - 1 lagged differences of
#' y and q - 1 lagged differences of x in addition to the contemporaneous
#' difference of x, so p and q are at least 1 and ARDL(1, 1) contains no
#' lagged differences at all. The paper is not consistent on this point:
#' equations (1) and (24) sum to p - 1 and q - 1, whereas Section 4.1
#' describes the ARDL(1, 1) data generating process as having one lagged
#' difference of each variable and Table 4 reports an ARDL(0, 2) model with
#' no lagged difference of y. This package follows equation (1); the paper's
#' ARDL(0, 2) corresponds to \code{fixlag = c(1, 3)} here, and
#' \code{fixlag = c(p, q)} in versions before 1.1.0 corresponds to
#' \code{fixlag = c(p + 1, q + 1)} in 1.1.0 (apart from the contemporaneous
#' difference of x, which earlier versions omitted).
#'
#' Two statistics are computed (Section 3.3, equations (11) and (12)):
#' \itemize{
#'   \item the t statistic on \eqn{y_{t-1}} for \eqn{H_0: b_1 = 0} against
#'     \eqn{H_1: b_1 < 0} (left tailed);
#'   \item the F (Wald) statistic for the joint hypothesis
#'     \eqn{H_0: \beta_2 = 0} on the lagged levels of all covariates
#'     (right tailed). With one covariate this equals the squared t statistic
#'     on \eqn{x_{t-1}}.
#' }
#'
#' Their null distributions depend on nuisance parameters, so both are
#' bootstrapped separately with the respective null imposed: the restricted
#' regression (without \eqn{y_{t-1}} for the t test, without \eqn{x_{t-1}} for
#' the F test) is estimated, its recentred residuals are resampled with
#' replacement, a bootstrap series \eqn{y^*} is generated recursively from the
#' restricted estimated equation with the observed x held fixed, and the
#' unrestricted regression is re-estimated on \eqn{(y^*, x)}. The
#' contemporaneous difference of x, which equation (24) of the paper omits
#' although equation (1) has it, is included in the restricted regressions
#' and in the bootstrap, following equation (1). The initial values
#' \eqn{y^*_1, \ldots, y^*_{\max(p, q)}} are set to the observed y, so the
#' recursion starts at t = max(p, q) + 1 with the observed history. Critical
#' values
#' are the percentiles of the \code{nboot} bootstrap statistics, and the
#' bootstrap p-values are the proportions of bootstrap statistics at least as
#' extreme as the observed ones.
#'
#' The lag orders p and q of the bootstrap are those of the test regression,
#' either fixed through \code{fixlag} or selected by AIC or BIC on the observed
#' data. They are held fixed in every bootstrap replication: the restricted
#' regression is estimated with these orders and the unrestricted regression
#' is re-estimated on \eqn{(y^*, x)} with the same orders, and the lag
#' selection is not repeated on the bootstrap samples. This follows Section
#' 4.2 of the paper, where the restricted regression (24) is estimated with
#' the lag length of the data generating process and equation (1) is
#' re-estimated on \eqn{y^*_t} and \eqn{x_t}. When the orders are selected
#' from the data, the bootstrap distribution therefore does not reflect the
#' selection step. In a Monte Carlo experiment by the package author (n = 100,
#' y and x independent random walks, 200 replications, 199 bootstrap
#' replications, 5 percent level) the rejection frequencies of the bootstrap
#' tests were 0.055 to 0.070 with fixed lag orders and 0.080 to 0.100 with
#' AIC selection (Monte Carlo standard error about 0.02). Fixing the orders
#' through \code{fixlag} avoids this source of size distortion.
#'
#' The decision at significance level \code{level} classifies the series
#' according to the four cases of Section 3.2 of the paper:
#' \itemize{
#'   \item \strong{Case I} (t not rejected, F not rejected): nonstationary
#'     process, no cointegration; y is I(1).
#'   \item \strong{Case II} (t rejected, F not rejected): stationary process;
#'     y is I(0).
#'   \item \strong{Case III} (t not rejected, F rejected): degenerate lagged
#'     dependent variable; y is I(2).
#'   \item \strong{Case IV} (both rejected): nonstationary process with
#'     cointegration between y and x; y is I(1).
#' }
#' The Monte Carlo experiments of the paper use a 5 percent level for both
#' tests, which is the default here. In the empirical application of the paper
#' the F statistic is reported as significant at the 10 percent level and
#' cointegration is concluded (Section 6); set \code{level = 0.10} to
#' reproduce that convention. Before version 1.1.0 \code{level} was a
#' confidence level (default 0.95); values of 0.5 or more are now rejected.
#'
#' Lag orders are either fixed through \code{fixlag} or selected by AIC or BIC
#' over p = 1, ..., \code{maxlag} and q = 1, ..., \code{maxlag}. All candidate
#' models are estimated on the same observations (t = \code{maxlag} + 1, ...,
#' n) so that the criteria are comparable; the selected model is then
#' re-estimated on its own full sample.
#'
#' @param y A numeric vector or time series. The series under test.
#' @param x A numeric vector, time series or numeric matrix with one column
#'   per covariate (the I(1) forcing variables of the paper).
#' @param case Integer. Deterministic specification:
#'   \itemize{
#'     \item \code{1}: No deterministic terms
#'     \item \code{3}: Intercept only (default)
#'     \item \code{5}: Intercept and linear trend
#'   }
#' @param maxlag Integer. Maximum ARDL order for AIC/BIC selection. Default
#'   is 10. Must be between 1 and 12.
#' @param ic Character. Information criterion for lag selection: \code{"aic"}
#'   (default) or \code{"bic"}.
#' @param fixlag Optional numeric vector of length 2, the fixed orders
#'   \code{c(p, q)} of the ARDL(p, q) model (see Details). Both must be at
#'   least 1. If provided, overrides automatic lag selection.
#' @param nboot Integer. Number of bootstrap replications. Default is 999.
#'   Minimum is 99.
#' @param level Numeric. Significance level of the two tests used for the
#'   case classification (0 to 1). Default is 0.05.
#' @param seed Optional integer. If supplied, the random number generator is
#'   seeded with it for the bootstrap and the caller's RNG state is restored
#'   on exit. Default \code{NULL} leaves the RNG untouched.
#' @param boot Logical. Whether to compute bootstrap critical values and
#'   p-values. Default is \code{TRUE}.
#' @param reps Deprecated name for \code{nboot}, kept for backward
#'   compatibility.
#'
#' @return An object of class \code{"mvardlurt"} containing:
#'   \item{tstat}{t statistic for \eqn{H_0: b_1 = 0}}
#'   \item{fstat}{F statistic for \eqn{H_0: \beta_2 = 0}}
#'   \item{t_pval, f_pval}{Bootstrap p-values (NA when \code{boot = FALSE})}
#'   \item{t_cv, f_cv}{Bootstrap critical values at the 10, 5, 2.5 and 1
#'     percent levels}
#'   \item{t_boot, f_boot}{The bootstrap distributions}
#'   \item{b1, b1_se}{Estimate and standard error of \eqn{b_1}}
#'   \item{beta2, beta2_se}{Estimates and standard errors of \eqn{\beta_2}
#'     (one per covariate)}
#'   \item{lr_mult}{Long-run multipliers \eqn{-\beta_2 / b_1}}
#'   \item{opt_p, opt_q}{Orders of the ARDL(p, q) model used}
#'   \item{case, casename}{Deterministic case and its description}
#'   \item{nboot}{Number of bootstrap replications}
#'   \item{nobs}{Number of observations in the regression}
#'   \item{aic, bic, r_squared, adj_r_squared}{Fit statistics of the model}
#'   \item{ic_table}{Matrix of information criteria for all (p, q) candidates
#'     (NULL with \code{fixlag})}
#'   \item{decision}{List with the decision at \code{level}: \code{reject_t},
#'     \code{reject_f}, \code{case_num}, \code{case_result},
#'     \code{integration} (the implied order of integration of y) and the
#'     significance codes}
#'   \item{model}{The fitted \code{lm} object, for inspection (its call
#'     refers to an internal data frame; to refit use
#'     \code{lm(formula(object$model), data = object$data)})}
#'   \item{data}{The data frame of the test regression}
#'   \item{y, x}{The data used (x as a matrix)}
#'   \item{residuals}{Residuals from the fitted model}
#'
#' @references
#' Sam, C. Y., McNown, R., Goh, S. K. and Goh, K. L. (2025). A multivariate
#' autoregressive distributed lag unit root test. \emph{Studies in Economics
#' and Econometrics}, 49(1), 17-33.
#' \doi{10.1080/03796205.2024.2439101}
#'
#' @examples
#' # Cointegrated data: y is I(1) through x (Case IV)
#' set.seed(123)
#' n <- 200
#' x <- cumsum(rnorm(n))
#' y <- 0.5 * x + rnorm(n, sd = 0.5)
#'
#' result <- mvardlurt(y, x, case = 3, fixlag = c(2, 2), nboot = 199)
#' print(result)
#' summary(result)
#'
#' @export
mvardlurt <- function(y, x, case = 3L, maxlag = 10L, ic = "aic",
                      fixlag = NULL, nboot = 999L, level = 0.05,
                      seed = NULL, boot = TRUE, reps = NULL) {

  if (!is.null(reps)) {
    warning("'reps' is deprecated; use 'nboot' instead.")
    nboot <- reps
  }

  # Input validation
  if (!is.numeric(y) || !is.numeric(x)) {
    stop("'y' and 'x' must be numeric.")
  }
  y <- as.numeric(y)
  if (is.matrix(x)) {
    X <- x
  } else {
    X <- matrix(as.numeric(x), ncol = 1L)
  }
  storage.mode(X) <- "double"
  if (nrow(X) != length(y)) {
    stop("'y' and 'x' must have the same length.")
  }
  k <- ncol(X)
  xn <- colnames(X)
  if (is.null(xn) || any(xn == "") || anyDuplicated(xn) || any(xn == "y")) {
    xn <- if (k == 1L) "x" else paste0("x", seq_len(k))
  }
  colnames(X) <- make.names(xn)

  # Remove NAs (complete cases)
  complete <- complete.cases(y, X)
  y <- y[complete]
  X <- X[complete, , drop = FALSE]

  n <- length(y)
  if (n < 30) {
    stop(sprintf("Too few observations (%d). Need at least 30.", n))
  }

  case <- as.integer(case)
  if (!case %in% c(1L, 3L, 5L)) {
    stop("'case' must be 1 (none), 3 (intercept), or 5 (intercept + trend).")
  }

  maxlag <- as.integer(maxlag)
  if (is.na(maxlag) || maxlag < 1 || maxlag > 12) {
    stop("'maxlag' must be between 1 and 12.")
  }

  ic <- tolower(ic)
  if (!ic %in% c("aic", "bic")) {
    stop("'ic' must be 'aic' or 'bic'.")
  }

  nboot <- as.integer(nboot)
  if (is.na(nboot) || nboot < 99) {
    stop("'nboot' must be at least 99.")
  }

  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("'level' must be between 0 and 1 (exclusive).")
  }
  if (level >= 0.5) {
    stop("'level' is the significance level since mvardlurt 1.1.0 ",
         "(default 0.05); did you mean 1 - level?")
  }

  if (!is.null(fixlag)) {
    if (length(fixlag) != 2 || !is.numeric(fixlag)) {
      stop("'fixlag' must be a numeric vector of length 2: c(p, q).")
    }
    fixlag <- as.integer(fixlag)
    if (any(fixlag < 1)) {
      stop("'fixlag' values must be at least 1 (ARDL(p, q) with p, q >= 1).")
    }
  }

  # Seed handling: never leave the global RNG state modified
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  casename <- switch(as.character(case),
                     "1" = "No Deterministic Terms",
                     "3" = "Intercept Only",
                     "5" = "Intercept and Trend")

  # Lag selection or fixed lags
  if (!is.null(fixlag)) {
    opt_p <- fixlag[1]
    opt_q <- fixlag[2]
    ic_table <- NULL
    manual_lag <- TRUE
  } else {
    lag_result <- .select_lags(y, X, case, maxlag, ic)
    opt_p <- lag_result$opt_p
    opt_q <- lag_result$opt_q
    ic_table <- lag_result$ic_table
    manual_lag <- FALSE
  }

  # Estimate the test regression
  design <- .ardl_design(y, X, opt_p, opt_q, case)
  model <- lm(.ardl_formula(design), data = design$data)
  nobs <- nrow(design$data)
  msum <- summary(model)
  coefs <- coef(model)
  se <- msum$coefficients[, "Std. Error"]

  b1 <- as.numeric(coefs[design$ty])
  b1_se <- as.numeric(se[design$ty])
  tstat <- b1 / b1_se

  beta2 <- coefs[design$tx]
  beta2_se <- se[design$tx]
  names(beta2) <- names(beta2_se) <- colnames(X)

  # Joint Wald F statistic on the lagged levels of x
  V <- vcov(model)[design$tx, design$tx, drop = FALSE]
  fstat <- as.numeric(crossprod(beta2, solve(V, beta2))) / k

  lr_mult <- if (b1 != 0) -beta2 / b1 else rep(NA_real_, k)

  # Bootstrap
  if (boot) {
    br <- .bootstrap_cv(y, X, opt_p, opt_q, case, nboot)
    t_cv <- br$t_cv
    f_cv <- br$f_cv
    t_boot <- br$t_boot
    f_boot <- br$f_boot
  } else {
    t_cv <- f_cv <- c(cv10 = NA_real_, cv05 = NA_real_, cv025 = NA_real_,
                      cv01 = NA_real_)
    t_boot <- f_boot <- NULL
  }

  decision <- .make_decision(tstat, fstat, t_boot, f_boot, t_cv, f_cv, level,
                             boot)

  result <- list(
    tstat = tstat,
    fstat = fstat,
    t_pval = decision$t_pval,
    f_pval = decision$f_pval,
    b1 = b1,
    b1_se = b1_se,
    beta2 = beta2,
    beta2_se = beta2_se,
    lr_mult = lr_mult,
    opt_p = opt_p,
    opt_q = opt_q,
    case = case,
    casename = casename,
    nboot = nboot,
    nobs = nobs,
    aic = AIC(model),
    bic = BIC(model),
    r_squared = msum$r.squared,
    adj_r_squared = msum$adj.r.squared,
    t_cv = t_cv,
    f_cv = f_cv,
    t_boot = t_boot,
    f_boot = f_boot,
    ic_table = ic_table,
    ic = ic,
    decision = decision,
    model = model,
    data = design$data,
    y = y,
    x = X,
    residuals = residuals(model),
    manual_lag = manual_lag,
    boot = boot,
    level = level,
    seed = seed,
    call = match.call()
  )

  class(result) <- "mvardlurt"
  result
}
