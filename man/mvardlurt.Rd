\name{mvardlurt}
\alias{mvardlurt}
\title{Multivariate ARDL Unit Root Test}
\description{
Implements the multivariate autoregressive distributed lag (ARDL) unit root
test proposed by Sam, McNown, Goh, and Goh (2024). The test augments the
standard ADF regression with lagged levels of a covariate (independent
variable) to improve power, especially when cointegration exists. Bootstrap
critical values ensure correct size regardless of nuisance parameters.
}
\usage{
mvardlurt(y, x, case = 3L, maxlag = 10L, ic = "aic",
          fixlag = NULL, reps = 1000L, level = 0.95,
          seed = 12345L, boot = TRUE)
}
\arguments{
\item{y}{A numeric vector or time series. The dependent variable.}

\item{x}{A numeric vector or time series. The independent variable (covariate).}

\item{case}{Integer. Deterministic specification:
  \itemize{
    \item \code{1}: No deterministic terms
    \item \code{3}: Intercept only (default)
    \item \code{5}: Intercept and linear trend
  }}

\item{maxlag}{Integer. Maximum lag order for AIC/BIC selection. Default is
  10. Must be between 0 and 10.}

\item{ic}{Character. Information criterion for lag selection: \code{"aic"}
  (default) or \code{"bic"}.}

\item{fixlag}{Optional numeric vector of length 2, specifying fixed lag
  orders \code{c(p, q)} for \eqn{\Delta y} and \eqn{\Delta x} respectively.
  If provided, overrides automatic lag selection.}

\item{reps}{Integer. Number of bootstrap replications. Default is 1000.
  Minimum is 100.}

\item{level}{Numeric. Confidence level for inference (0 to 1). Default is
  0.95.}

\item{seed}{Integer. Random seed for reproducibility. Default is 12345.}

\item{boot}{Logical. Whether to compute bootstrap critical values. Default
  is \code{TRUE}.}
}
\details{
The test estimates the following ARDL regression:

\deqn{\Delta y_t = \pi y_{t-1} + \delta x_{t-1} + \sum_{j=1}^{p} \gamma_j
\Delta y_{t-j} + \sum_{j=1}^{q} \theta_j \Delta x_{t-j} + deterministics +
\varepsilon_t}

The test produces two statistics:
\itemize{
  \item \strong{t-statistic}: Tests \eqn{H_0: \pi = 0} (unit root in y)
  \item \strong{F-statistic}: Tests \eqn{H_0: \delta = 0} (no cointegration)
}

Based on the four-case framework (Sam et al., 2024):
\itemize{
  \item \strong{Case I}: Reject both \eqn{\rightarrow} Cointegration
  \item \strong{Case II}: Reject t, Accept F \eqn{\rightarrow} Degenerate case 1
    (y may be I(0))
  \item \strong{Case III}: Accept t, Reject F \eqn{\rightarrow} Degenerate case 2
    (spurious)
  \item \strong{Case IV}: Accept both \eqn{\rightarrow} No cointegration
}
}
\value{
An object of class \code{"mvardlurt"} containing:
\item{tstat}{t-statistic for the unit root test (on \eqn{\pi})}
\item{fstat}{F-statistic for the cointegration test (on \eqn{\delta})}
\item{fstat_p}{Asymptotic p-value for the F-statistic}
\item{pi_coef}{Coefficient estimate of \eqn{\pi} (lagged y)}
\item{pi_se}{Standard error of \eqn{\pi}}
\item{delta_coef}{Coefficient estimate of \eqn{\delta} (lagged x)}
\item{delta_se}{Standard error of \eqn{\delta}}
\item{lr_mult}{Long-run multiplier \eqn{-\delta/\pi} (if \eqn{\pi \neq 0})}
\item{opt_p}{Selected lag order for \eqn{\Delta y}}
\item{opt_q}{Selected lag order for \eqn{\Delta x}}
\item{case}{Deterministic case used}
\item{casename}{Description of the deterministic case}
\item{reps}{Number of bootstrap replications}
\item{nobs}{Number of observations used}
\item{aic}{AIC value of the selected model}
\item{bic}{BIC value of the selected model}
\item{r_squared}{R-squared of the regression}
\item{t_cv}{Bootstrap critical values for t-statistic (10\%, 5\%, 2.5\%, 1\%)}
\item{f_cv}{Bootstrap critical values for F-statistic (10\%, 5\%, 2.5\%, 1\%)}
\item{ic_table}{Matrix of IC values for all (p, q) combinations}
\item{decision}{List containing test decisions and significance levels}
\item{model}{The fitted \code{lm} object}
\item{y}{Original dependent variable}
\item{x}{Original independent variable}
\item{residuals}{Residuals from the fitted model}
}
\references{
Sam, C. Y., McNown, R., Goh, S. K., & Goh, K. L. (2024). A multivariate
autoregressive distributed lag unit root test. \emph{Studies in Economics
and Econometrics}, 1-17.
\doi{10.1080/03796205.2024.2439101}
}
\examples{
# Generate example data with cointegration
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.5)

# Run the test (reduced reps for speed)
result <- mvardlurt(y, x, case = 3, reps = 200)
print(result)

# With fixed lags
result2 <- mvardlurt(y, x, fixlag = c(2, 2), reps = 200)
}
\author{
Merwan Roudane
}
