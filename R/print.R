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
  cat("  Multivariate ARDL Unit Root Test -- Sam, McNown, Goh and Goh (2024)\n")
  cat("\n")

  # Basic info
  cat(sprintf("  Optimal Model: ARDL(%d, %d)    Case %d (%s)\n",
              x$opt_p, x$opt_q, x$case, x$casename))
  cat(sprintf("  Observations: %d    R-squared: %.6f\n",
              x$nobs, x$r_squared))
  cat("\n")

  # Test statistics
  line <- paste(rep("-", 60), collapse = "")
  cat(line, "\n")
  cat("  Test Statistics:\n")
  cat(line, "\n")
  cat(sprintf("    t-statistic: %10.4f  (H0: pi = 0, unit root)\n", x$tstat))
  cat(sprintf("    F-statistic: %10.4f  (H0: delta = 0, no cointegration)\n",
              x$fstat))
  cat(line, "\n")

  # Bootstrap critical values
  if (x$boot && !any(is.na(x$t_cv))) {
    cat("\n")
    cat(sprintf("  Bootstrap Critical Values (%d replications):\n", x$reps))
    cat(line, "\n")
    cat(sprintf("    Sig. Level   %10s %10s %10s %10s\n",
                "10%", "5%", "2.5%", "1%"))
    cat(line, "\n")
    cat(sprintf("    t-critical   %10.4f %10.4f %10.4f %10.4f\n",
                x$t_cv["cv10"], x$t_cv["cv05"], x$t_cv["cv025"], x$t_cv["cv01"]))
    cat(sprintf("    F-critical   %10.4f %10.4f %10.4f %10.4f\n",
                x$f_cv["cv10"], x$f_cv["cv05"], x$f_cv["cv025"], x$f_cv["cv01"]))
    cat(line, "\n")
  }

  # Decision
  if (x$boot && !is.na(x$decision$case_num)) {
    cat("\n")
    cat("  Decision:\n")
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

    cat(sprintf("    t-test (H0: pi = 0):    %s\n", t_decision))
    cat(sprintf("    F-test (H0: delta = 0): %s\n", f_decision))
    cat("\n")
    cat(sprintf("  => CASE %s: %s\n",
                as.roman(x$decision$case_num), x$decision$case_result))
    cat(line, "\n")
  }

  cat("\n")
  cat("  Significance: *** 1%  ** 2.5%  * 5%  + 10%  n.s. not significant\n")
  cat("\n")

  invisible(x)
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
  cat("  Table 1: ARDL Unit Root Test\n")
  cat("  Sam, McNown, Goh and Goh (2024)\n")
  cat(line, "\n")

  cat(sprintf("  Optimal Model: ARDL(%d, %d)        Case: %d (%s)\n",
              x$opt_p, x$opt_q, x$case, x$casename))
  cat(sprintf("  Observations: %d                  R-squared: %.6f\n",
              x$nobs, x$r_squared))
  cat(sprintf("  AIC: %.4f                        BIC: %.4f\n",
              x$aic, x$bic))
  cat(line, "\n")
  cat(sprintf("  t-statistic: %12.6f        (H0: pi = 0, unit root)\n",
              x$tstat))
  cat(sprintf("  F-statistic: %12.6f        (H0: delta = 0, no cointegration)\n",
              x$fstat))
  cat(line, "\n")
  cat("\n")

  # Bootstrap critical values
  if (x$boot && !any(is.na(x$t_cv))) {
    cat(line, "\n")
    cat(sprintf("  Table 2: Bootstrap Critical Values (%d replications)\n",
                x$reps))
    cat(line, "\n")
    cat(sprintf("  Sig. Level     %12s %12s %12s %12s\n",
                "10%", "5%", "2.5%", "1%"))
    cat(line2, "\n")
    cat(sprintf("  t-critical     %12.4f %12.4f %12.4f %12.4f\n",
                x$t_cv["cv10"], x$t_cv["cv05"], x$t_cv["cv025"],
                x$t_cv["cv01"]))
    cat(sprintf("  F-critical     %12.4f %12.4f %12.4f %12.4f\n",
                x$f_cv["cv10"], x$f_cv["cv05"], x$f_cv["cv025"],
                x$f_cv["cv01"]))
    cat(line, "\n")
    cat("\n")
  }

  # Coefficient summary
  cat(line, "\n")
  cat("  Table 3: ARDL Coefficient Summary\n")
  cat(line, "\n")
  cat(sprintf("  %-20s %12s %12s %12s\n",
              "Parameter", "Coefficient", "Std. Error", "t-stat"))
  cat(line2, "\n")
  cat(sprintf("  pi (unit root)     %12.6f %12.6f %12.4f\n",
              x$pi_coef, x$pi_se, x$tstat))
  cat(sprintf("  delta (cointegr.)  %12.6f %12.6f\n",
              x$delta_coef, x$delta_se))
  if (!is.na(x$lr_mult) && x$pi_coef != 0) {
    cat(line2, "\n")
    cat(sprintf("  Long-run multiplier (-delta/pi): %12.6f\n", x$lr_mult))
  }
  cat(line, "\n")
  cat("\n")

  # Decision and inference
  if (x$boot && !is.na(x$decision$case_num)) {
    cat(line, "\n")
    cat("  Table 4: Decision and Inference\n")
    cat(line, "\n")
    cat("\n")

    # Part A: Hypothesis Tests
    cat("  A. Hypothesis Tests\n")
    cat(line2, "\n")
    cat(sprintf("  %-10s %-20s %12s %15s %8s\n",
                "Test", "Null Hypothesis", "Statistic", "Decision", "Sig."))
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

    cat(sprintf("  %-10s %-20s %12.4f %15s %8s\n",
                "t-test", "H0: pi = 0", x$tstat, t_decision, x$decision$t_level))
    cat(sprintf("  %-10s %-20s %12.4f %15s %8s\n",
                "F-test", "H0: delta = 0", x$fstat, f_decision,
                x$decision$f_level))
    cat(line2, "\n")
    cat("\n")

    # Part B: Four Cases Framework
    cat("  B. Four-Case Framework (Sam, McNown, Goh and Goh, 2024)\n")
    cat(line2, "\n")
    cat(sprintf("  %2s %-6s %-10s %-10s %-35s\n",
                "", "Case", "t-test", "F-test", "Interpretation"))
    cat(line2, "\n")

    cases <- list(
      list(num = "I", t = "Reject", f = "Reject", interp = "Cointegration"),
      list(num = "II", t = "Reject", f = "Accept",
           interp = "Degenerate case 1 (y may be I(0))"),
      list(num = "III", t = "Accept", f = "Reject",
           interp = "Degenerate case 2 (spurious)"),
      list(num = "IV", t = "Accept", f = "Accept",
           interp = "No cointegration")
    )

    for (i in seq_along(cases)) {
      mark <- if (i == x$decision$case_num) "=>" else "  "
      cat(sprintf("  %2s %-6s %-10s %-10s %-35s\n",
                  mark, cases[[i]]$num, cases[[i]]$t, cases[[i]]$f,
                  cases[[i]]$interp))
    }
    cat(line2, "\n")
    cat("\n")

    # Part C: Conclusion
    cat("  C. Conclusion\n")
    cat(line2, "\n")

    if (x$decision$case_num == 1L) {
      cat("  CASE I: Cointegration\n")
      cat("    pi != 0 : Error correction exists; y adjusts to equilibrium\n")
      cat("    delta != 0 : x enters the long-run equation\n")
      cat("    => y and x are cointegrated\n")
      cat("    => The ECM is valid; long-run equilibrium relationship exists\n")
    } else if (x$decision$case_num == 2L) {
      cat("  CASE II: Degenerate case 1\n")
      cat("    pi != 0 : y is stationary or error-correcting\n")
      cat("    delta = 0 : x has no long-run effect on y\n")
      cat("    => y may be I(0); x is not a cointegrator\n")
      cat("    => No long-run relationship via x\n")
    } else if (x$decision$case_num == 3L) {
      cat("  CASE III: Degenerate case 2\n")
      cat("    pi = 0 : y has a unit root (I(1))\n")
      cat("    delta != 0 : x coefficient is significant\n")
      cat("    => Spurious result; unit root invalidates the relationship\n")
      cat("    => The long-run relationship is not reliable\n")
    } else {
      cat("  CASE IV: No cointegration\n")
      cat("    pi = 0 : y has a unit root (I(1))\n")
      cat("    delta = 0 : No cointegrating relationship\n")
      cat("    => No evidence of long-run equilibrium\n")
      cat("    => y and x are not cointegrated\n")
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
#' @return Named numeric vector of key coefficients.
#'
#' @export
#' @method coef mvardlurt
coef.mvardlurt <- function(object, ...) {
  result <- c(object$pi_coef, object$delta_coef, object$lr_mult)
  names(result) <- c("pi", "delta", "lr_mult")
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
