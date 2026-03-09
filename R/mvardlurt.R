#' Multivariate ARDL Unit Root Test
#'
#' @description
#' Implements the multivariate autoregressive distributed lag (ARDL) unit root
#' test proposed by Sam, McNown, Goh, and Goh (2024). The test augments the
#' standard ADF regression with lagged levels of a covariate (independent
#' variable) to improve power, especially when cointegration exists. Bootstrap
#' critical values ensure correct size regardless of nuisance parameters.
#'
#' @details
#' The test estimates the following ARDL regression:
#'
#' \deqn{\Delta y_t = \pi y_{t-1} + \delta x_{t-1} + \sum_{j=1}^{p} \gamma_j
#' \Delta y_{t-j} + \sum_{j=1}^{q} \theta_j \Delta x_{t-j} + deterministics +
#' \varepsilon_t}
#'
#' The test produces two statistics:
#' \itemize{
#'   \item **t-statistic**: Tests \eqn{H_0: \pi = 0} (unit root in y)
#'   \item **F-statistic**: Tests \eqn{H_0: \delta = 0} (no cointegration)
#' }
#'
#' Based on the four-case framework (Sam et al., 2024):
#' \itemize{
#'   \item **Case I**: Reject both \eqn{\rightarrow} Cointegration
#'   \item **Case II**: Reject t, Accept F \eqn{\rightarrow} Degenerate case 1
#'     (y may be I(0))
#'   \item **Case III**: Accept t, Reject F \eqn{\rightarrow} Degenerate case 2
#'     (spurious)
#'   \item **Case IV**: Accept both \eqn{\rightarrow} No cointegration
#' }
#'
#' @param y A numeric vector or time series. The dependent variable.
#' @param x A numeric vector or time series. The independent variable
#'   (covariate).
#' @param case Integer. Deterministic specification:
#'   \itemize{
#'     \item \code{1}: No deterministic terms
#'     \item \code{3}: Intercept only (default)
#'     \item \code{5}: Intercept and linear trend
#'   }
#' @param maxlag Integer. Maximum lag order for AIC/BIC selection. Default is
#'   10. Must be between 0 and 10.
#' @param ic Character. Information criterion for lag selection: \code{"aic"}
#'   (default) or \code{"bic"}.
#' @param fixlag Optional numeric vector of length 2, specifying fixed lag
#'   orders \code{c(p, q)} for \eqn{\Delta y} and \eqn{\Delta x} respectively.
#'   If provided, overrides automatic lag selection.
#' @param reps Integer. Number of bootstrap replications. Default is 1000.
#'   Minimum is 100.
#' @param level Numeric. Confidence level for inference (0 to 1). Default is
#'   0.95.
#' @param seed Integer. Random seed for reproducibility. Default is 12345.
#' @param boot Logical. Whether to compute bootstrap critical values. Default
#'   is \code{TRUE}.
#'
#' @return An object of class \code{"mvardlurt"} containing:
#'   \item{tstat}{t-statistic for the unit root test (on \eqn{\pi})}
#'   \item{fstat}{F-statistic for the cointegration test (on \eqn{\delta})}
#'   \item{fstat_p}{Asymptotic p-value for the F-statistic}
#'   \item{pi_coef}{Coefficient estimate of \eqn{\pi} (lagged y)}
#'   \item{pi_se}{Standard error of \eqn{\pi}}
#'   \item{delta_coef}{Coefficient estimate of \eqn{\delta} (lagged x)}
#'   \item{delta_se}{Standard error of \eqn{\delta}}
#'   \item{lr_mult}{Long-run multiplier \eqn{-\delta/\pi} (if \eqn{\pi \neq 0})}
#'   \item{opt_p}{Selected lag order for \eqn{\Delta y}}
#'   \item{opt_q}{Selected lag order for \eqn{\Delta x}}
#'   \item{case}{Deterministic case used}
#'   \item{casename}{Description of the deterministic case}
#'   \item{reps}{Number of bootstrap replications}
#'   \item{nobs}{Number of observations used}
#'   \item{aic}{AIC value of the selected model}
#'   \item{bic}{BIC value of the selected model}
#'   \item{r_squared}{R-squared of the regression}
#'   \item{t_cv}{Bootstrap critical values for t-statistic (10\%, 5\%, 2.5\%,
#'     1\%)}
#'   \item{f_cv}{Bootstrap critical values for F-statistic (10\%, 5\%, 2.5\%,
#'     1\%)}
#'   \item{ic_table}{Matrix of IC values for all (p, q) combinations}
#'   \item{decision}{List containing test decisions and significance levels}
#'   \item{model}{The fitted \code{lm} object}
#'   \item{y}{Original dependent variable}
#'   \item{x}{Original independent variable}
#'   \item{residuals}{Residuals from the fitted model}
#'
#' @references
#' Sam, C. Y., McNown, R., Goh, S. K., & Goh, K. L. (2024). A multivariate
#' autoregressive distributed lag unit root test. \emph{Studies in Economics
#' and Econometrics}, 1-17.
#' \doi{10.1080/03796205.2024.2439101}
#'
#' @examples
#' # Generate example data with cointegration
#' set.seed(123)
#' n <- 200
#' x <- cumsum(rnorm(n))  # I(1) process
#' y <- 0.5 * x + rnorm(n, sd = 0.5)  # Cointegrated with x
#'
#' # Run the test
#' result <- mvardlurt(y, x, case = 3, reps = 500)
#' print(result)
#' summary(result)
#'
#' # With fixed lags
#' result2 <- mvardlurt(y, x, fixlag = c(2, 2))
#'
#' @export
mvardlurt <- function(y, x, case = 3L, maxlag = 10L, ic = "aic",
                      fixlag = NULL, reps = 1000L, level = 0.95,
                      seed = 12345L, boot = TRUE) {


  # Input validation
  if (!is.numeric(y) || !is.numeric(x)) {
    stop("'y' and 'x' must be numeric vectors.")
  }

  if (length(y) != length(x)) {
    stop("'y' and 'x' must have the same length.")
  }

  y <- as.numeric(y)
  x <- as.numeric(x)

  # Remove NAs (complete cases)
  complete <- complete.cases(y, x)
  y <- y[complete]
  x <- x[complete]

  n <- length(y)
  if (n < 30) {
    stop(sprintf("Too few observations (%d). Need at least 30.", n))
  }

  # Validate case
  case <- as.integer(case)
  if (!case %in% c(1L, 3L, 5L)) {
    stop("'case' must be 1 (none), 3 (intercept), or 5 (intercept + trend).")
  }

  # Validate maxlag
  maxlag <- as.integer(maxlag)
  if (maxlag < 0 || maxlag > 10) {
    stop("'maxlag' must be between 0 and 10.")
  }

  # Validate IC
  ic <- tolower(ic)
  if (!ic %in% c("aic", "bic")) {
    stop("'ic' must be 'aic' or 'bic'.")
  }

  # Validate reps
  reps <- as.integer(reps)
  if (reps < 100) {
    stop("'reps' must be at least 100.")
  }

  # Validate level
  if (level <= 0 || level >= 1) {
    stop("'level' must be between 0 and 1 (exclusive).")
  }

  # Validate fixlag
  if (!is.null(fixlag)) {
    if (length(fixlag) != 2 || !is.numeric(fixlag)) {
      stop("'fixlag' must be a numeric vector of length 2: c(p, q).")
    }
    fixlag <- as.integer(fixlag)
    if (any(fixlag < 0)) {
      stop("'fixlag' values must be non-negative.")
    }
  }

  # Set seed
  set.seed(seed)

  # Case labels
  casename <- switch(
    as.character(case),
    "1" = "No Deterministic Terms",
    "3" = "Intercept Only",
    "5" = "Intercept and Trend"
  )

  # Compute differences
  dy <- diff(y)
  dx <- diff(x)

  # Lag selection or use fixed lags
  if (!is.null(fixlag)) {
    opt_p <- fixlag[1]
    opt_q <- fixlag[2]
    ic_table <- NULL
    manual_lag <- TRUE
  } else {
    lag_result <- .select_lags(y, x, dy, dx, case, maxlag, ic)
    opt_p <- lag_result$opt_p
    opt_q <- lag_result$opt_q
    ic_table <- lag_result$ic_table
    manual_lag <- FALSE
  }

  # Estimate the optimal model
  model_result <- .estimate_ardl(y, x, dy, dx, opt_p, opt_q, case)

  # Extract results
  model <- model_result$model
  nobs <- model_result$nobs
  r2 <- summary(model)$r.squared
  adj_r2 <- summary(model)$adj.r.squared
  ll <- as.numeric(logLik(model))
  nparams <- length(coef(model))
  aic_val <- AIC(model)
  bic_val <- BIC(model)

  # Extract pi (coefficient on L.y) and delta (coefficient on L.x)
  coefs <- coef(model)
  se <- summary(model)$coefficients[, "Std. Error"]

  pi_coef <- coefs["L.y"]
  pi_se <- se["L.y"]
  tstat <- pi_coef / pi_se

  delta_coef <- coefs["L.x"]
  delta_se <- se["L.x"]

  # F-test on delta
  fstat_result <- .wald_test(model, "L.x")
  fstat <- fstat_result$fstat
  fstat_p <- fstat_result$pval

  # Long-run multiplier
  if (pi_coef != 0) {
    lr_mult <- -delta_coef / pi_coef
  } else {
    lr_mult <- NA_real_
  }

  # Bootstrap critical values
  if (boot) {
    boot_result <- .bootstrap_cv(y, x, opt_p, opt_q, case, reps)
    t_cv <- boot_result$t_cv
    f_cv <- boot_result$f_cv
  } else {
    t_cv <- c(cv10 = NA_real_, cv05 = NA_real_, cv025 = NA_real_,
              cv01 = NA_real_)
    f_cv <- c(cv10 = NA_real_, cv05 = NA_real_, cv025 = NA_real_,
              cv01 = NA_real_)
  }

  # Determine significance and decisions
  decision <- .make_decision(tstat, fstat, t_cv, f_cv, boot)

  # Build result object
  result <- list(
    tstat = as.numeric(tstat),
    fstat = fstat,
    fstat_p = fstat_p,
    pi_coef = as.numeric(pi_coef),
    pi_se = as.numeric(pi_se),
    delta_coef = as.numeric(delta_coef),
    delta_se = as.numeric(delta_se),
    lr_mult = lr_mult,
    opt_p = opt_p,
    opt_q = opt_q,
    case = case,
    casename = casename,
    reps = reps,
    nobs = nobs,
    aic = aic_val,
    bic = bic_val,
    r_squared = r2,
    adj_r_squared = adj_r2,
    t_cv = t_cv,
    f_cv = f_cv,
    ic_table = ic_table,
    ic = ic,
    decision = decision,
    model = model,
    y = y,
    x = x,
    residuals = residuals(model),
    manual_lag = manual_lag,
    boot = boot,
    level = level,
    seed = seed,
    call = match.call()
  )

  class(result) <- "mvardlurt"
  return(result)
}
