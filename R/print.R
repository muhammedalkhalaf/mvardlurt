#' Print method for mvardlurt objects
#'
#' @param x An object of class \code{"mvardlurt"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
#' @method print mvardlurt
print.mvardlurt <- function(x, ...) {

  cat("\n")
  cat("  Multivariate ARDL Unit Root Test -- Sam, McNown, Goh and Goh (2025)\n")
  cat("\n")

  cat(sprintf("  Model: ARDL(%d, %d)    Deterministic case %d (%s)\n",
              x$opt_p, x$opt_q, x$case, x$casename))
  cat(sprintf("  Covariates: %s\n", paste(colnames(x$x), collapse = ", ")))
  cat(sprintf("  Observations: %d    R-squared: %.6f\n",
              x$nobs, x$r_squared))
  cat("\n")

  line <- paste(rep("-", 66), collapse = "")
  cat(line, "\n")
  cat("  Test Statistics:\n")
  cat(line, "\n")
  cat(sprintf("    t-statistic: %10.4f  (H0: b1 = 0, coefficient on y[t-1])\n",
              x$tstat))
  cat(sprintf("    F-statistic: %10.4f  (H0: beta2 = 0, coefficients on x[t-1])\n",
              x$fstat))
  cat(line, "\n")

  if (x$boot && !any(is.na(x$t_cv))) {
    cat("\n")
    cat(sprintf("  Bootstrap Critical Values (%d replications):\n", x$nboot))
    cat(line, "\n")
    cat(sprintf("    Sig. Level   %10s %10s %10s %10s %10s\n",
                "10%", "5%", "2.5%", "1%", "p-value"))
    cat(line, "\n")
    cat(sprintf("    t-critical   %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                x$t_cv["cv10"], x$t_cv["cv05"], x$t_cv["cv025"],
                x$t_cv["cv01"], x$t_pval))
    cat(sprintf("    F-critical   %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                x$f_cv["cv10"], x$f_cv["cv05"], x$f_cv["cv025"],
                x$f_cv["cv01"], x$f_pval))
    cat(line, "\n")
  }

  if (x$boot && !is.na(x$decision$case_num)) {
    cat("\n")
    cat(sprintf("  Decision at the %s level:\n", .fmt_level(x$level)))
    cat(line, "\n")

    t_decision <- if (x$decision$reject_t) {
      paste0("Reject", x$decision$t_sig, " (", x$decision$t_level, ")")
    } else {
      "Fail to reject"
    }
    f_decision <- if (x$decision$reject_f) {
      paste0("Reject", x$decision$f_sig, " (", x$decision$f_level, ")")
    } else {
      "Fail to reject"
    }

    cat(sprintf("    t-test (H0: b1 = 0):     %s\n", t_decision))
    cat(sprintf("    F-test (H0: beta2 = 0):  %s\n", f_decision))
    cat("\n")
    cat(sprintf("  => CASE %s: %s\n",
                as.roman(x$decision$case_num), x$decision$case_result))
    cat(sprintf("     y is %s\n", x$decision$integration))
    cat(line, "\n")
  }

  cat("\n")
  cat("  Significance: *** 1%  ** 2.5%  * 5%  + 10%  n.s. not significant\n")
  cat("\n")

  invisible(x)
}


#' Format a significance level as a percentage
#' @param level Numeric level.
#' @return Character string.
#' @keywords internal
.fmt_level <- function(level) {
  paste0(format(100 * level, trim = TRUE), "%")
}


#' Summary method for mvardlurt objects
#'
#' @param object An object of class \code{"mvardlurt"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
#' @method summary mvardlurt
summary.mvardlurt <- function(object, ...) {

  x <- object
  line <- paste(rep("=", 78), collapse = "")
  line2 <- paste(rep("-", 78), collapse = "")

  cat("\n")
  cat(line, "\n")
  cat("  Table 1: Multivariate ARDL Unit Root Test\n")
  cat("  Sam, McNown, Goh and Goh (2025)\n")
  cat(line, "\n")

  cat(sprintf("  Model: ARDL(%d, %d)                Case: %d (%s)\n",
              x$opt_p, x$opt_q, x$case, x$casename))
  cat(sprintf("  Covariates: %s\n", paste(colnames(x$x), collapse = ", ")))
  cat(sprintf("  Observations: %d                  R-squared: %.6f\n",
              x$nobs, x$r_squared))
  cat(sprintf("  AIC: %.4f                        BIC: %.4f\n",
              x$aic, x$bic))
  cat(line, "\n")
  cat(sprintf("  t-statistic: %12.6f        (H0: b1 = 0)\n", x$tstat))
  cat(sprintf("  F-statistic: %12.6f        (H0: beta2 = 0)\n", x$fstat))
  cat(line, "\n")
  cat("\n")

  if (x$boot && !any(is.na(x$t_cv))) {
    cat(line, "\n")
    cat(sprintf("  Table 2: Bootstrap Critical Values and p-values (%d replications)\n",
                x$nboot))
    cat(line, "\n")
    cat(sprintf("  Sig. Level     %12s %12s %12s %12s %10s\n",
                "10%", "5%", "2.5%", "1%", "p-value"))
    cat(line2, "\n")
    cat(sprintf("  t-critical     %12.4f %12.4f %12.4f %12.4f %10.4f\n",
                x$t_cv["cv10"], x$t_cv["cv05"], x$t_cv["cv025"],
                x$t_cv["cv01"], x$t_pval))
    cat(sprintf("  F-critical     %12.4f %12.4f %12.4f %12.4f %10.4f\n",
                x$f_cv["cv10"], x$f_cv["cv05"], x$f_cv["cv025"],
                x$f_cv["cv01"], x$f_pval))
    cat(line, "\n")
    cat("\n")
  }

  cat(line, "\n")
  cat("  Table 3: Long-run Coefficients\n")
  cat(line, "\n")
  cat(sprintf("  %-22s %12s %12s %12s\n",
              "Parameter", "Coefficient", "Std. Error", "t-stat"))
  cat(line2, "\n")
  cat(sprintf("  %-22s %12.6f %12.6f %12.4f\n", "b1 (y[t-1])",
              x$b1, x$b1_se, x$tstat))
  for (j in seq_along(x$beta2)) {
    cat(sprintf("  %-22s %12.6f %12.6f %12.4f\n",
                paste0("beta2 (", names(x$beta2)[j], "[t-1])"),
                x$beta2[j], x$beta2_se[j], x$beta2[j] / x$beta2_se[j]))
  }
  if (all(!is.na(x$lr_mult))) {
    cat(line2, "\n")
    for (j in seq_along(x$lr_mult)) {
      cat(sprintf("  Long-run multiplier (-beta2/b1) for %s: %12.6f\n",
                  names(x$beta2)[j], x$lr_mult[j]))
    }
  }
  cat(line, "\n")
  cat("\n")

  if (x$boot && !is.na(x$decision$case_num)) {
    cat(line, "\n")
    cat(sprintf("  Table 4: Decision and Inference at the %s level\n",
                .fmt_level(x$level)))
    cat(line, "\n")
    cat("\n")

    cat("  A. Hypothesis Tests (bootstrap)\n")
    cat(line2, "\n")
    cat(sprintf("  %-10s %-16s %12s %10s %15s %8s\n",
                "Test", "Null Hypothesis", "Statistic", "p-value",
                "Decision", "Sig."))
    cat(line2, "\n")

    t_decision <- if (x$decision$reject_t) {
      paste0("Reject", x$decision$t_sig)
    } else {
      "Fail to reject"
    }
    f_decision <- if (x$decision$reject_f) {
      paste0("Reject", x$decision$f_sig)
    } else {
      "Fail to reject"
    }

    cat(sprintf("  %-10s %-16s %12.4f %10.4f %15s %8s\n",
                "t-test", "H0: b1 = 0", x$tstat, x$t_pval, t_decision,
                x$decision$t_level))
    cat(sprintf("  %-10s %-16s %12.4f %10.4f %15s %8s\n",
                "F-test", "H0: beta2 = 0", x$fstat, x$f_pval, f_decision,
                x$decision$f_level))
    cat(line2, "\n")
    cat("\n")

    cat("  B. Four Cases (Section 3.2 of Sam, McNown, Goh and Goh, 2025)\n")
    cat(line2, "\n")
    cat(sprintf("  %2s %-5s %-9s %-9s %-6s %-40s\n",
                "", "Case", "t-test", "F-test", "y is", "Interpretation"))
    cat(line2, "\n")

    cases <- list(
      list(num = "I", t = "Accept", f = "Accept", ord = "I(1)",
           interp = "Nonstationary, no cointegration"),
      list(num = "II", t = "Reject", f = "Accept", ord = "I(0)",
           interp = "Stationary process"),
      list(num = "III", t = "Accept", f = "Reject", ord = "I(2)",
           interp = "Degenerate lagged dependent variable"),
      list(num = "IV", t = "Reject", f = "Reject", ord = "I(1)",
           interp = "Nonstationary, cointegrated with x")
    )

    for (i in seq_along(cases)) {
      mark <- if (i == x$decision$case_num) "=>" else "  "
      cat(sprintf("  %2s %-5s %-9s %-9s %-6s %-40s\n",
                  mark, cases[[i]]$num, cases[[i]]$t, cases[[i]]$f,
                  cases[[i]]$ord, cases[[i]]$interp))
    }
    cat(line2, "\n")
    cat("\n")

    cat("  C. Conclusion\n")
    cat(line2, "\n")

    if (x$decision$case_num == 1L) {
      cat("  CASE I: Nonstationary process, no cointegration\n")
      cat("    b1 = 0 and beta2 = 0: the coefficient on y[t-1] is unity in levels\n")
      cat("    => y is I(1); no long-run relationship with x is detected\n")
    } else if (x$decision$case_num == 2L) {
      cat("  CASE II: Stationary process\n")
      cat("    b1 < 0 and beta2 = 0: y reverts to its mean without a long-run\n")
      cat("    contribution from the levels of x\n")
      cat("    => y is I(0)\n")
    } else if (x$decision$case_num == 3L) {
      cat("  CASE III: Degenerate lagged dependent variable\n")
      cat("    b1 = 0 and beta2 != 0: y inherits the unit root of x on top of\n")
      cat("    its own unit root\n")
      cat("    => y is I(2)\n")
    } else {
      cat("  CASE IV: Nonstationary process, cointegration\n")
      cat("    b1 < 0 and beta2 != 0: y is cointegrated with the I(1) covariates\n")
      cat("    and inherits their nonstationarity\n")
      cat("    => y is I(1)\n")
    }
    cat(line2, "\n")
    cat("\n")
  }

  cat("  Significance: *** 1%  ** 2.5%  * 5%  + 10%  n.s. not significant\n")
  cat(line, "\n")
  cat("\n")

  invisible(x)
}


#' Coef method for mvardlurt objects
#'
#' @param object An object of class \code{"mvardlurt"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector with \code{b1}, the \code{beta2} coefficients
#'   (one per covariate) and the long-run multipliers.
#'
#' @export
#' @method coef mvardlurt
coef.mvardlurt <- function(object, ...) {
  result <- c(object$b1, object$beta2, object$lr_mult)
  names(result) <- c("b1", paste0("beta2.", names(object$beta2)),
                     paste0("lr_mult.", names(object$beta2)))
  result
}


#' Extract residuals from mvardlurt objects
#'
#' @param object An object of class \code{"mvardlurt"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Numeric vector of residuals.
#'
#' @export
#' @method residuals mvardlurt
residuals.mvardlurt <- function(object, ...) {
  object$residuals
}


#' Extract fitted values from mvardlurt objects
#'
#' @param object An object of class \code{"mvardlurt"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Numeric vector of fitted values.
#'
#' @export
#' @method fitted mvardlurt
fitted.mvardlurt <- function(object, ...) {
  fitted(object$model)
}
