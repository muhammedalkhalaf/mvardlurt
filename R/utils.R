#' Internal utility functions for mvardlurt
#'
#' @name utils
#' @keywords internal
NULL

#' Build the ARDL test regression (equation (1) of Sam, McNown, Goh and Goh)
#'
#' Constructs the dependent variable and the design matrix of
#' \deqn{\Delta y_t = c_1 + c_2 t + b_1 y_{t-1} + \beta_2' x_{t-1}
#'   + \sum_{i=1}^{p-1} \phi_i \Delta y_{t-i}
#'   + \sum_{j=1}^{q-1} \Phi_j' \Delta x_{t-j} + \omega' \Delta x_t + u_t}
#' for t = start, ..., n, where start = max(p, q) + 1 unless a larger value is
#' supplied (used to put all candidate models on a common sample).
#'
#' @param y Numeric vector, the series under test.
#' @param X Numeric matrix (n x k) of covariates.
#' @param p Integer, ARDL order of y (p - 1 lagged differences of y).
#' @param q Integer, ARDL order of x (q - 1 lagged differences of x, plus the
#'   contemporaneous difference).
#' @param case Deterministic case (1, 3, or 5).
#' @param start Optional first time index of the regression sample.
#'
#' @return A list with \code{Y} (response), \code{Z} (design matrix with column
#'   names), \code{idx} (time indices used), \code{ty} (name of the y_{t-1}
#'   column), \code{tx} (names of the x_{t-1} columns) and \code{data} (the
#'   same information as a data frame, for \code{lm}).
#' @keywords internal
.ardl_design <- function(y, X, p, q, case, start = NULL) {
  n <- length(y)
  k <- ncol(X)
  xn <- colnames(X)
  minstart <- max(p, q) + 1L
  if (is.null(start)) start <- minstart
  if (start < minstart) stop("'start' is smaller than the lag structure allows.")
  if (start > n - 5L) {
    stop("Not enough observations for the specified lag structure.")
  }
  idx <- start:n

  dy <- c(NA_real_, diff(y))
  dX <- rbind(NA_real_, diff(X))
  dX <- matrix(dX, nrow = n, ncol = k)

  cols <- list()
  if (case %in% c(3L, 5L)) cols[["(Intercept)"]] <- rep(1, length(idx))
  if (case == 5L) cols[["trend"]] <- as.numeric(idx)
  cols[["L.y"]] <- y[idx - 1L]
  for (j in seq_len(k)) cols[[paste0("L.", xn[j])]] <- X[idx - 1L, j]
  if (p > 1L) {
    for (i in seq_len(p - 1L)) cols[[paste0("L", i, ".dy")]] <- dy[idx - i]
  }
  for (j in seq_len(k)) cols[[paste0("d.", xn[j])]] <- dX[idx, j]
  if (q > 1L) {
    for (i in seq_len(q - 1L)) {
      for (j in seq_len(k)) {
        cols[[paste0("L", i, ".d", xn[j])]] <- dX[idx - i, j]
      }
    }
  }
  Z <- do.call(cbind, cols)
  colnames(Z) <- names(cols)
  Y <- dy[idx]

  data <- as.data.frame(Z[, setdiff(colnames(Z), "(Intercept)"), drop = FALSE])
  data <- cbind(dy = Y, data)

  list(Y = Y, Z = Z, idx = idx, ty = "L.y", tx = paste0("L.", xn),
       data = data, has_intercept = case %in% c(3L, 5L))
}


#' Formula for the ARDL test regression
#'
#' @param design Output of \code{.ardl_design}.
#' @param drop Character vector of regressors to exclude (for restricted
#'   models).
#' @return A formula.
#' @keywords internal
.ardl_formula <- function(design, drop = character(0)) {
  regs <- setdiff(colnames(design$Z), c("(Intercept)", drop))
  rhs <- paste(paste0("`", regs, "`"), collapse = " + ")
  if (!design$has_intercept) rhs <- paste(rhs, "- 1")
  as.formula(paste("dy ~", rhs))
}


#' OLS t-statistic on y_{t-1} and joint F-statistic on x_{t-1}
#'
#' Fast QR-based computation used inside the bootstrap. Gives the same values
#' as \code{lm} and \code{anova}.
#'
#' @param Y Response vector.
#' @param Z Design matrix with named columns.
#' @param ty Name of the y_{t-1} column.
#' @param tx Names of the x_{t-1} columns.
#' @return Numeric vector \code{c(t = , F = )}.
#' @keywords internal
.ols_stats <- function(Y, Z, ty, tx) {
  fit <- lm.fit(Z, Y)
  b <- fit$coefficients
  df <- length(Y) - fit$rank
  s2 <- sum(fit$residuals^2) / df
  R <- chol2inv(fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE])
  pv <- fit$qr$pivot[seq_len(fit$rank)]
  V <- matrix(NA_real_, ncol(Z), ncol(Z), dimnames = list(colnames(Z), colnames(Z)))
  V[pv, pv] <- s2 * R
  tstat <- b[ty] / sqrt(V[ty, ty])
  bx <- b[tx]
  Vx <- V[tx, tx, drop = FALSE]
  fstat <- as.numeric(crossprod(bx, solve(Vx, bx))) / length(tx)
  c(t = as.numeric(tstat), F = fstat)
}


#' Select the ARDL orders (p, q) by AIC or BIC on a common sample
#'
#' All candidate models with p = 1, ..., maxlag and q = 1, ..., maxlag are
#' estimated on the same observations (t = maxlag + 1, ..., n), so that the
#' information criteria are comparable.
#'
#' @param y Series under test.
#' @param X Covariate matrix.
#' @param case Deterministic case.
#' @param maxlag Maximum order.
#' @param ic "aic" or "bic".
#' @return List with opt_p, opt_q and ic_table.
#' @keywords internal
.select_lags <- function(y, X, case, maxlag, ic) {
  ic_table <- matrix(NA_real_, nrow = maxlag, ncol = maxlag,
                     dimnames = list(paste0("p=", 1:maxlag),
                                     paste0("q=", 1:maxlag)))
  start <- maxlag + 1L
  best <- Inf
  opt_p <- opt_q <- 1L
  for (p in 1:maxlag) {
    for (q in 1:maxlag) {
      d <- .ardl_design(y, X, p, q, case, start = start)
      fit <- lm(.ardl_formula(d), data = d$data)
      val <- if (ic == "aic") AIC(fit) else BIC(fit)
      ic_table[p, q] <- val
      if (is.finite(val) && val < best) {
        best <- val
        opt_p <- p
        opt_q <- q
      }
    }
  }
  list(opt_p = as.integer(opt_p), opt_q = as.integer(opt_q),
       ic_table = ic_table)
}


#' Residual bootstrap of the t and F statistics (Section 4.2 of the paper)
#'
#' For each statistic the regression is estimated with its null imposed
#' (b1 = 0 for the t test, beta2 = 0 for the F test), the recentred restricted
#' residuals are resampled with replacement, the bootstrap series y* is built
#' recursively from the restricted estimated equation with the observed x held
#' fixed, and the unrestricted regression is re-estimated on (y*, x).
#'
#' @param y Series under test.
#' @param X Covariate matrix.
#' @param p,q ARDL orders.
#' @param case Deterministic case.
#' @param nboot Number of bootstrap replications.
#' @return List with \code{t_boot}, \code{f_boot} (bootstrap distributions),
#'   \code{t_cv}, \code{f_cv} (10, 5, 2.5 and 1 percent critical values).
#' @keywords internal
.bootstrap_cv <- function(y, X, p, q, case, nboot) {
  d <- .ardl_design(y, X, p, q, case)
  t_boot <- .boot_one(y, d, p, drop = d$ty, stat = "t", nboot)
  f_boot <- .boot_one(y, d, p, drop = d$tx, stat = "F", nboot)

  t_cv <- quantile(t_boot, probs = c(0.10, 0.05, 0.025, 0.01), names = FALSE)
  f_cv <- quantile(f_boot, probs = c(0.90, 0.95, 0.975, 0.99), names = FALSE)
  names(t_cv) <- names(f_cv) <- c("cv10", "cv05", "cv025", "cv01")
  list(t_boot = t_boot, f_boot = f_boot, t_cv = t_cv, f_cv = f_cv)
}


#' One bootstrap distribution (t or F) with a given null imposed
#'
#' @param y Series under test.
#' @param d Output of \code{.ardl_design} for the unrestricted model.
#' @param p ARDL order of y.
#' @param drop Columns removed to impose the null.
#' @param stat "t" or "F".
#' @param nboot Number of replications.
#' @return Numeric vector of bootstrap statistics.
#' @keywords internal
.boot_one <- function(y, d, p, drop, stat, nboot) {
  Z <- d$Z
  keep <- setdiff(colnames(Z), drop)
  Zr <- Z[, keep, drop = FALSE]
  rfit <- lm.fit(Zr, d$Y)
  b <- rfit$coefficients
  b[is.na(b)] <- 0
  e <- rfit$residuals
  e <- e - mean(e)

  # Split the restricted equation into the part driven by y* (its own lagged
  # level and lagged differences) and the part fixed in the bootstrap
  # (deterministics, x_{t-1}, current and lagged differences of x).
  ycols <- c("L.y", if (p > 1L) paste0("L", seq_len(p - 1L), ".dy"))
  ycols <- intersect(ycols, keep)
  xcols <- setdiff(keep, ycols)
  exog <- if (length(xcols)) as.numeric(Zr[, xcols, drop = FALSE] %*% b[xcols]) else
    numeric(length(d$Y))

  b1 <- if ("L.y" %in% ycols) b["L.y"] else 0
  phi <- if (p > 1L) b[paste0("L", seq_len(p - 1L), ".dy")] else numeric(0)
  # Levels recursion: y_t = a_1 y_{t-1} + ... + a_p y_{t-p} + w_t, where
  # a_1 = 1 + b1 + phi_1, a_i = phi_i - phi_{i-1}, a_p = -phi_{p-1}.
  if (p == 1L) {
    a <- 1 + b1
  } else {
    a <- c(1 + b1 + phi[1], if (p > 2L) diff(phi), -phi[p - 1L])
  }
  a <- as.numeric(a)
  start <- d$idx[1]
  init <- y[(start - 1L):(start - p)]
  m <- length(d$Y)

  out <- numeric(nboot)
  for (bb in seq_len(nboot)) {
    w <- exog + sample(e, m, replace = TRUE)
    ystar <- y
    ystar[d$idx] <- as.numeric(stats::filter(w, a, method = "recursive",
                                             init = init))
    # Re-estimate the unrestricted regression on (y*, x): only the columns
    # built from y change.
    Zb <- Z
    Zb[, "L.y"] <- ystar[d$idx - 1L]
    if (p > 1L) {
      dys <- c(NA_real_, diff(ystar))
      for (i in seq_len(p - 1L)) Zb[, paste0("L", i, ".dy")] <- dys[d$idx - i]
    }
    Yb <- ystar[d$idx] - ystar[d$idx - 1L]
    s <- .ols_stats(Yb, Zb, d$ty, d$tx)
    out[bb] <- s[stat]
  }
  out
}


#' Classify the outcome according to Section 3.2 of the paper
#'
#' @param tstat,fstat Observed statistics.
#' @param t_boot,f_boot Bootstrap distributions (may be NULL).
#' @param t_cv,f_cv Critical values at 10, 5, 2.5 and 1 percent.
#' @param level Significance level used for the decision.
#' @param boot Whether a bootstrap was run.
#' @return List with decision information.
#' @keywords internal
.make_decision <- function(tstat, fstat, t_boot, f_boot, t_cv, f_cv, level,
                           boot) {
  if (!boot || is.null(t_boot) || is.null(f_boot)) {
    return(list(level = level, t_pval = NA_real_, f_pval = NA_real_,
                t_cv_level = NA_real_, f_cv_level = NA_real_,
                t_sig = "", t_level = "N/A", f_sig = "", f_level = "N/A",
                reject_t = NA, reject_f = NA,
                case_num = NA_integer_, case_result = "N/A",
                integration = "N/A"))
  }

  t_pval <- mean(t_boot <= tstat)
  f_pval <- mean(f_boot >= fstat)
  t_cv_level <- quantile(t_boot, probs = level, names = FALSE)
  f_cv_level <- quantile(f_boot, probs = 1 - level, names = FALSE)
  reject_t <- as.logical(tstat < t_cv_level)
  reject_f <- as.logical(fstat > f_cv_level)

  sig <- function(reject_at) {
    if (reject_at[4]) return(c("***", "1%"))
    if (reject_at[3]) return(c("**", "2.5%"))
    if (reject_at[2]) return(c("*", "5%"))
    if (reject_at[1]) return(c("+", "10%"))
    c("", "n.s.")
  }
  ts <- sig(tstat < t_cv)
  fs <- sig(fstat > f_cv)

  if (reject_t && reject_f) {
    case_num <- 4L
    case_result <- "Nonstationary process, cointegration"
    integration <- "I(1)"
  } else if (reject_t && !reject_f) {
    case_num <- 2L
    case_result <- "Stationary process, degenerate lagged independent variable"
    integration <- "I(0)"
  } else if (!reject_t && reject_f) {
    case_num <- 3L
    case_result <- "Second order integration, degenerate lagged dependent variable"
    integration <- "I(2)"
  } else {
    case_num <- 1L
    case_result <- "Nonstationary process, no cointegration"
    integration <- "I(1)"
  }

  list(level = level, t_pval = t_pval, f_pval = f_pval,
       t_cv_level = t_cv_level, f_cv_level = f_cv_level,
       t_sig = ts[1], t_level = ts[2], f_sig = fs[1], f_level = fs[2],
       reject_t = reject_t, reject_f = reject_f,
       case_num = case_num, case_result = case_result,
       integration = integration)
}
