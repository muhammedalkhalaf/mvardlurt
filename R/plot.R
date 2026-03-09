#' Plot method for mvardlurt objects
#'
#' Creates diagnostic plots for the multivariate ARDL unit root test.
#'
#' @param x An object of class \code{"mvardlurt"}.
#' @param which Integer vector indicating which plots to produce:
#'   \itemize{
#'     \item \code{1}: Residuals vs Fitted
#'     \item \code{2}: Q-Q plot of residuals
#'     \item \code{3}: Residuals over time
#'     \item \code{4}: ACF of residuals
#'     \item \code{5}: Time series of y and x
#'     \item \code{6}: Information criterion surface
#'   }
#'   Default is \code{c(1, 2, 3, 4)}.
#' @param ask Logical. If \code{TRUE}, prompt before each plot. Default is
#'   \code{TRUE} when multiple plots are requested in an interactive session.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
#' @method plot mvardlurt
plot.mvardlurt <- function(x, which = c(1, 2, 3, 4),
                           ask = (length(which) > 1 && dev.interactive()),
                           ...) {

  # Store original par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  if (ask) {
    par(ask = TRUE)
  }

  resid <- residuals(x)
  fitted_vals <- fitted(x)

  for (w in which) {
    if (w == 1) {
      # Residuals vs Fitted
      plot(fitted_vals, resid,
           xlab = "Fitted values",
           ylab = "Residuals",
           main = "Residuals vs Fitted",
           pch = 20, col = "steelblue", ...)
      abline(h = 0, lty = 2, col = "red")
      lines(lowess(fitted_vals, resid), col = "red", lwd = 2)

    } else if (w == 2) {
      # Q-Q plot
      qqnorm(resid, main = "Normal Q-Q Plot of Residuals",
             pch = 20, col = "steelblue", ...)
      qqline(resid, col = "red", lwd = 2)

    } else if (w == 3) {
      # Residuals over time
      plot(seq_along(resid), resid,
           type = "l",
           xlab = "Observation",
           ylab = "Residuals",
           main = "Residuals over Time",
           col = "steelblue", lwd = 1.5, ...)
      abline(h = 0, lty = 2, col = "red")

    } else if (w == 4) {
      # ACF of residuals
      acf(resid, main = "ACF of Residuals", ...)

    } else if (w == 5) {
      # Time series plot
      par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

      plot(x$y, type = "l",
           xlab = "Observation", ylab = "y",
           main = "Dependent Variable (y)",
           col = "steelblue", lwd = 1.5, ...)

      plot(x$x, type = "l",
           xlab = "Observation", ylab = "x",
           main = "Independent Variable (x)",
           col = "darkgreen", lwd = 1.5, ...)

      par(mfrow = c(1, 1))

    } else if (w == 6) {
      # IC surface
      if (!is.null(x$ic_table)) {
        ic_mat <- x$ic_table
        ic_name <- toupper(x$ic)

        # Create heatmap
        image(x = 0:(nrow(ic_mat) - 1),
              y = 0:(ncol(ic_mat) - 1),
              z = ic_mat,
              xlab = "p (lags of dy)",
              ylab = "q (lags of dx)",
              main = paste(ic_name, "Values for ARDL(p, q)"),
              col = heat.colors(12, rev = TRUE), ...)

        # Add grid lines
        abline(h = 0:(ncol(ic_mat) - 1) - 0.5, col = "gray80", lty = 3)
        abline(v = 0:(nrow(ic_mat) - 1) - 0.5, col = "gray80", lty = 3)

        # Mark optimal
        points(x$opt_p, x$opt_q, pch = 4, cex = 2, lwd = 3, col = "blue")
        legend("topright",
               legend = sprintf("Optimal: (%d, %d)", x$opt_p, x$opt_q),
               pch = 4, col = "blue", bg = "white")
      } else {
        plot.new()
        text(0.5, 0.5, "IC table not available\n(manual lag selection used)",
             cex = 1.2)
      }

    } else {
      warning(sprintf("Plot type %d not recognized.", w))
    }
  }

  invisible(x)
}


#' Create a combined diagnostic plot
#'
#' @param x An object of class \code{"mvardlurt"}.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
autoplot.mvardlurt <- function(x, ...) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  resid <- residuals(x)
  fitted_vals <- fitted(x)

  # 1. Residuals vs Fitted
  plot(fitted_vals, resid,
       xlab = "Fitted", ylab = "Residuals",
       main = "Residuals vs Fitted",
       pch = 20, col = "steelblue")
  abline(h = 0, lty = 2, col = "red")

  # 2. Q-Q plot
  qqnorm(resid, main = "Normal Q-Q",
         pch = 20, col = "steelblue")
  qqline(resid, col = "red", lwd = 2)

  # 3. Residuals over time
  plot(seq_along(resid), resid,
       type = "l",
       xlab = "Observation", ylab = "Residuals",
       main = "Residuals over Time",
       col = "steelblue", lwd = 1.5)
  abline(h = 0, lty = 2, col = "red")

  # 4. ACF
  acf(resid, main = "ACF of Residuals")

  par(mfrow = c(1, 1))

  invisible(x)
}
