\name{plot.mvardlurt}
\alias{plot.mvardlurt}
\title{Plot Method for mvardlurt Objects}
\description{
Creates diagnostic plots for the multivariate ARDL unit root test.
}
\usage{
\method{plot}{mvardlurt}(x, which = c(1, 2, 3, 4),
     ask = (length(which) > 1 && dev.interactive()), \dots)
}
\arguments{
\item{x}{An object of class \code{"mvardlurt"}.}

\item{which}{Integer vector indicating which plots to produce:
  \itemize{
    \item \code{1}: Residuals vs Fitted
    \item \code{2}: Q-Q plot of residuals
    \item \code{3}: Residuals over time
    \item \code{4}: ACF of residuals
    \item \code{5}: Time series of y and x
    \item \code{6}: Information criterion surface
  }
  Default is \code{c(1, 2, 3, 4)}.}

\item{ask}{Logical. If \code{TRUE}, prompt before each plot. Default is
  \code{TRUE} when multiple plots are requested in an interactive session.}

\item{\dots}{Additional arguments passed to plotting functions.}
}
\value{
Invisibly returns \code{x}.
}
\examples{
set.seed(123)
n <- 100
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.5)

result <- mvardlurt(y, x, reps = 100)

# Default diagnostic plots
plot(result, ask = FALSE)

# IC surface plot
plot(result, which = 6)
}
\seealso{
\code{\link{mvardlurt}}, \code{\link{autoplot.mvardlurt}}
}
\author{
Merwan Roudane
}
