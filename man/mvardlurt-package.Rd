\name{mvardlurt-package}
\alias{mvardlurt-package}
\docType{package}
\title{
Multivariate ARDL Unit Root Test
}
\description{
Implements the multivariate autoregressive distributed lag (ARDL) unit root
test proposed by Sam, McNown, Goh, and Goh (2024). The test augments the
standard ADF regression with lagged levels of a covariate to improve power
when cointegration exists. Bootstrap critical values ensure correct size
regardless of nuisance parameters.
}
\details{
The main function is \code{\link{mvardlurt}}, which performs the multivariate
ARDL unit root test. The package provides:

\itemize{
  \item Automatic lag selection via AIC or BIC
  \item Bootstrap critical values for correct test size
  \item Three deterministic specifications (none, intercept, intercept + trend)
  \item Four-case decision framework for inference
  \item Comprehensive print and summary output
  \item Diagnostic plots
}

The test produces two statistics:
\itemize{
  \item \strong{t-statistic}: Tests \eqn{H_0: \pi = 0} (unit root)
  \item \strong{F-statistic}: Tests \eqn{H_0: \delta = 0} (no cointegration)
}

The four-case framework interprets results as:
\itemize{
  \item Case I: Both reject \eqn{\rightarrow} Cointegration
  \item Case II: Reject t, Accept F \eqn{\rightarrow} y may be I(0)
  \item Case III: Accept t, Reject F \eqn{\rightarrow} Spurious
  \item Case IV: Both accept \eqn{\rightarrow} No cointegration
}
}
\author{Muhammad Alkhalaf}
\references{
Sam, C. Y., McNown, R., Goh, S. K., & Goh, K. L. (2024). A multivariate
autoregressive distributed lag unit root test. \emph{Studies in Economics
and Econometrics}, 1-17.
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
result <- mvardlurt(y, x, case = 3, reps = 200)
print(result)
}
