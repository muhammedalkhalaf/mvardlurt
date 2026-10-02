\name{mvardlurt-package}
\alias{mvardlurt-package}
\docType{package}
\title{
Multivariate ARDL Unit Root Test
}
\description{
Implements the multivariate autoregressive distributed lag (ARDL) unit root
test of Sam, McNown, Goh and Goh (2025). The test augments the ADF
regression with the lagged level, the current difference and lagged
differences of one or more covariates, so that a possible cointegrating
relationship between the series under test and the covariates is taken into
account. Critical values and p-values come from the residual bootstrap
described in the paper.
}
\details{
The main function is \code{\link{mvardlurt}}, which performs the multivariate
ARDL unit root test. The package provides:

\itemize{
  \item Automatic lag selection via AIC or BIC on a common sample
  \item Residual bootstrap critical values and p-values with the null
    imposed (Section 4.2 of the paper)
  \item Three deterministic specifications (none, intercept, intercept + trend)
  \item One or several covariates (joint F test on their lagged levels)
  \item Four-case classification of the order of integration of y
  \item Comprehensive print and summary output
  \item Diagnostic plots
}

The test produces two statistics:
\itemize{
  \item \strong{t-statistic}: Tests \eqn{H_0: b_1 = 0} on the lagged
    level of y
  \item \strong{F-statistic}: Tests \eqn{H_0: \beta_2 = 0} jointly on the
    lagged levels of the covariates
}

The four cases of Section 3.2 of the paper are:
\itemize{
  \item Case I: neither rejected; y is I(1), no cointegration
  \item Case II: t rejected, F not rejected; y is I(0)
  \item Case III: t not rejected, F rejected; y is I(2)
  \item Case IV: both rejected; y is I(1) and cointegrated with x
}
}
\author{Muhammad Alkhalaf}
\references{
Sam, C. Y., McNown, R., Goh, S. K. and Goh, K. L. (2025). A multivariate
autoregressive distributed lag unit root test. \emph{Studies in Economics
and Econometrics}, 49(1), 17-33.
\doi{10.1080/03796205.2024.2439101}
}
\keyword{package}
\seealso{
\code{\link{mvardlurt}}
}
\examples{
# Generate cointegrated data
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.5)

# Run the test
result <- mvardlurt(y, x, case = 3, fixlag = c(2, 2), nboot = 199)
print(result)
}
