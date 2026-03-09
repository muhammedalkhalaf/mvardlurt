#' Internal utility functions for mvardlurt
#'
#' @name utils
#' @keywords internal
NULL

#' Select optimal lag orders using AIC or BIC
#'
#' @param y Original dependent variable
#' @param x Original independent variable
#' @param dy First difference of y
#' @param dx First difference of x
#' @param case Deterministic case (1, 3, or 5)
#' @param maxlag Maximum lag order
#' @param ic Information criterion ("aic" or "bic")
#'
#' @return List with opt_p, opt_q, and ic_table
#' @keywords internal
.select_lags <- function(y, x, dy, dx, case, maxlag, ic) {

  n_orig <- length(y)
  dim <- maxlag + 1
  ic_table <- matrix(NA_real_, nrow = dim, ncol = dim,
                     dimnames = list(paste0("p=", 0:maxlag),
                                     paste0("q=", 0:maxlag)))

  best_ic <- Inf
  opt_p <- 0
  opt_q <- 0

  for (p in 0:maxlag) {
    for (q in 0:maxlag) {

      # Build data frame for regression
      result <- tryCatch({
        model_data <- .build_model_data(y, x, dy, dx, p, q, case)
        if (nrow(model_data) < 15) {
          return(NULL)
        }

        # Fit model
        if (case == 1L) {
          formula <- .build_formula(p, q, case, intercept = FALSE)
          fit <- lm(formula, data = model_data)
        } else {
          formula <- .build_formula(p, q, case, intercept = TRUE)
          fit <- lm(formula, data = model_data)
        }

        # Compute IC
        if (ic == "aic") {
          ic_val <- AIC(fit)
        } else {
          ic_val <- BIC(fit)
        }

        list(fit = fit, ic_val = ic_val)
      }, error = function(e) NULL)

      if (is.null(result)) next

      ic_table[p + 1, q + 1] <- result$ic_val

      if (result$ic_val < best_ic) {
        best_ic <- result$ic_val
        opt_p <- p
        opt_q <- q
      }
    }
  }

  return(list(opt_p = opt_p, opt_q = opt_q, ic_table = ic_table))
}


#' Build model data frame for ARDL regression
#'
#' @param y Original dependent variable
#' @param x Original independent variable
#' @param dy First difference of y
#' @param dx First difference of x
#' @param p Lag order for dy
#' @param q Lag order for dx
#' @param case Deterministic case
#'
#' @return Data frame with all variables
#' @keywords internal
.build_model_data <- function(y, x, dy, dx, p, q, case) {

  n <- length(dy)

  # Start index to ensure all lags are available
  start <- max(p, q) + 1
  end <- n

  if (start > end) {
    stop("Not enough observations for the specified lag structure.")
  }

  # Dependent variable: dy[t]
  # We need dy from position 'start' to 'end'
  dy_dep <- dy[start:end]

  # Lagged levels: y[t-1] and x[t-1]
  # Since dy[t] = y[t] - y[t-1], and we're using indices in dy,

  # dy[start] corresponds to original y[start+1] - y[start]
  # So L.y for dy[start] is y[start], and L.x is x[start]
  L_y <- y[start:(end)]
  L_x <- x[start:(end)]

  df <- data.frame(
    dy = dy_dep,
    L.y = L_y,
    L.x = L_x
  )

  # Lagged differences of dy
  if (p > 0) {
    for (j in 1:p) {
      col_name <- paste0("L", j, ".dy")
      df[[col_name]] <- dy[(start - j):(end - j)]
    }
  }

  # Lagged differences of dx
  if (q > 0) {
    for (j in 1:q) {
      col_name <- paste0("L", j, ".dx")
      df[[col_name]] <- dx[(start - j):(end - j)]
    }
  }

  # Trend variable for case 5
  if (case == 5L) {
    df$trend <- seq_len(nrow(df))
  }

  return(df)
}


#' Build formula for ARDL regression
#'
#' @param p Lag order for dy
#' @param q Lag order for dx
#' @param case Deterministic case
#' @param intercept Whether to include intercept
#'
#' @return Formula object
#' @keywords internal
.build_formula <- function(p, q, case, intercept = TRUE) {

  # Base regressors
  regressors <- c("L.y", "L.x")

  # Lagged differences
  if (p > 0) {
    regressors <- c(regressors, paste0("L", 1:p, ".dy"))
  }
  if (q > 0) {
    regressors <- c(regressors, paste0("L", 1:q, ".dx"))
  }

  # Trend
  if (case == 5L) {
    regressors <- c(regressors, "trend")
  }

  # Build formula
  rhs <- paste(regressors, collapse = " + ")
  if (!intercept) {
    rhs <- paste(rhs, "- 1")
  }

  formula <- as.formula(paste("dy ~", rhs))
  return(formula)
}


#' Estimate the ARDL model with optimal lags
#'
#' @param y Original dependent variable
#' @param x Original independent variable
#' @param dy First difference of y
#' @param dx First difference of x
#' @param p Lag order for dy
#' @param q Lag order for dx
#' @param case Deterministic case
#'
#' @return List with model and nobs
#' @keywords internal
.estimate_ardl <- function(y, x, dy, dx, p, q, case) {

  model_data <- .build_model_data(y, x, dy, dx, p, q, case)

  if (case == 1L) {
    formula <- .build_formula(p, q, case, intercept = FALSE)
  } else {
    formula <- .build_formula(p, q, case, intercept = TRUE)
  }

  model <- lm(formula, data = model_data)

  return(list(model = model, nobs = nrow(model_data)))
}


#' Wald test for a single coefficient
#'
#' @param model lm object
#' @param varname Name of the variable to test
#'
#' @return List with fstat and pval
#' @keywords internal
.wald_test <- function(model, varname) {

  coefs <- coef(model)
  vcov_mat <- vcov(model)

  if (!varname %in% names(coefs)) {
    stop(sprintf("Variable '%s' not found in model.", varname))
  }

  idx <- which(names(coefs) == varname)
  beta <- coefs[idx]
  var_beta <- vcov_mat[idx, idx]

  # Wald F-statistic (chi-sq / 1 = F with df1=1)
  fstat <- (beta^2) / var_beta
  df2 <- model$df.residual
  pval <- pf(fstat, df1 = 1, df2 = df2, lower.tail = FALSE)

  return(list(fstat = as.numeric(fstat), pval = as.numeric(pval)))
}


#' Bootstrap critical values under the null hypothesis
#'
#' @param y Original dependent variable
#' @param x Original independent variable
#' @param p Lag order for dy
#' @param q Lag order for dx
#' @param case Deterministic case
#' @param reps Number of bootstrap replications
#'
#' @return List with t_cv and f_cv vectors
#' @keywords internal
.bootstrap_cv <- function(y, x, p, q, case, reps) {

  n <- length(y)

  # Compute differences
  dy <- diff(y)
  dx <- diff(x)

  # Estimate null model (imposing pi = 0, i.e., no error correction)
  # Under null: dy[t] = sum of lagged dy + sum of lagged dx + det + error
  # We need to generate data under the null of a unit root in y

  # Get residual variance from original model for scaling
  model_data <- .build_model_data(y, x, dy, dx, p, q, case)
  if (case == 1L) {
    formula <- .build_formula(p, q, case, intercept = FALSE)
  } else {
    formula <- .build_formula(p, q, case, intercept = TRUE)
  }
  model <- lm(formula, data = model_data)
  sigma_hat <- summary(model)$sigma

  # Storage for bootstrap statistics
  t_boot <- numeric(reps)
  f_boot <- numeric(reps)

  # Bootstrap loop
  for (b in 1:reps) {
    # Generate y under null (random walk)
    # y_boot[t] = y_boot[t-1] + epsilon
    eps_y <- rnorm(n, sd = sigma_hat)
    y_boot <- cumsum(eps_y)

    # x remains as observed (fixed in repeated samples) or generate as I(1)
    eps_x <- rnorm(n, sd = sd(dx))
    x_boot <- cumsum(eps_x)

    # Compute differences
    dy_boot <- diff(y_boot)
    dx_boot <- diff(x_boot)

    # Estimate model
    boot_result <- tryCatch({
      model_data_boot <- .build_model_data(y_boot, x_boot, dy_boot, dx_boot,
                                           p, q, case)

      if (case == 1L) {
        fit_boot <- lm(formula, data = model_data_boot)
      } else {
        fit_boot <- lm(formula, data = model_data_boot)
      }

      # Extract t-statistic on L.y
      coefs_boot <- coef(fit_boot)
      se_boot <- summary(fit_boot)$coefficients[, "Std. Error"]

      if (!"L.y" %in% names(coefs_boot)) {
        return(list(t = NA_real_, f = NA_real_))
      }

      t_stat <- coefs_boot["L.y"] / se_boot["L.y"]

      # F-statistic on L.x
      f_result <- .wald_test(fit_boot, "L.x")

      list(t = as.numeric(t_stat), f = f_result$fstat)
    }, error = function(e) {
      list(t = NA_real_, f = NA_real_)
    })

    t_boot[b] <- boot_result$t
    f_boot[b] <- boot_result$f
  }

  # Remove NAs
  t_boot <- t_boot[!is.na(t_boot)]
  f_boot <- f_boot[!is.na(f_boot)]

  # Compute critical values
  # For t-statistic (left-tailed test for unit root)
  t_cv <- quantile(t_boot, probs = c(0.10, 0.05, 0.025, 0.01), na.rm = TRUE)
  names(t_cv) <- c("cv10", "cv05", "cv025", "cv01")

  # For F-statistic (right-tailed test)
  f_cv <- quantile(f_boot, probs = c(0.90, 0.95, 0.975, 0.99), na.rm = TRUE)
  names(f_cv) <- c("cv10", "cv05", "cv025", "cv01")

  return(list(t_cv = t_cv, f_cv = f_cv))
}


#' Make decision based on test statistics and critical values
#'
#' @param tstat t-statistic
#' @param fstat F-statistic
#' @param t_cv Bootstrap critical values for t
#' @param f_cv Bootstrap critical values for F
#' @param boot Whether bootstrap was performed
#'
#' @return List with decision information
#' @keywords internal
.make_decision <- function(tstat, fstat, t_cv, f_cv, boot) {

  if (!boot || any(is.na(t_cv)) || any(is.na(f_cv))) {
    return(list(
      t_sig = "",
      t_level = "N/A",
      f_sig = "",
      f_level = "N/A",
      case_result = "N/A",
      case_num = NA_integer_,
      reject_t = NA,
      reject_f = NA
    ))
  }

  # t-test significance (left-tailed: reject if tstat < critical value)
  t_sig <- ""
  t_level <- "n.s."
  reject_t <- FALSE

  if (tstat < t_cv["cv01"]) {
    t_sig <- "***"
    t_level <- "1%"
    reject_t <- TRUE
  } else if (tstat < t_cv["cv025"]) {
    t_sig <- "**"
    t_level <- "2.5%"
    reject_t <- TRUE
  } else if (tstat < t_cv["cv05"]) {
    t_sig <- "*"
    t_level <- "5%"
    reject_t <- TRUE
  } else if (tstat < t_cv["cv10"]) {
    t_sig <- "+"
    t_level <- "10%"
    reject_t <- TRUE
  }

  # F-test significance (right-tailed: reject if fstat > critical value)
  f_sig <- ""
  f_level <- "n.s."
  reject_f <- FALSE

  if (fstat > f_cv["cv01"]) {
    f_sig <- "***"
    f_level <- "1%"
    reject_f <- TRUE
  } else if (fstat > f_cv["cv025"]) {
    f_sig <- "**"
    f_level <- "2.5%"
    reject_f <- TRUE
  } else if (fstat > f_cv["cv05"]) {
    f_sig <- "*"
    f_level <- "5%"
    reject_f <- TRUE
  } else if (fstat > f_cv["cv10"]) {
    f_sig <- "+"
    f_level <- "10%"
    reject_f <- TRUE
  }

  # Four-case framework
  if (reject_t && reject_f) {
    case_result <- "Cointegration"
    case_num <- 1L
  } else if (reject_t && !reject_f) {
    case_result <- "Degenerate case 1 (y may be I(0))"
    case_num <- 2L
  } else if (!reject_t && reject_f) {
    case_result <- "Degenerate case 2 (spurious)"
    case_num <- 3L
  } else {
    case_result <- "No cointegration"
    case_num <- 4L
  }

  return(list(
    t_sig = t_sig,
    t_level = t_level,
    f_sig = f_sig,
    f_level = f_level,
    case_result = case_result,
    case_num = case_num,
    reject_t = reject_t,
    reject_f = reject_f
  ))
}
