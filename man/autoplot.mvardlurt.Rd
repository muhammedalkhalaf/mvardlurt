\name{autoplot.mvardlurt}
\alias{autoplot.mvardlurt}
\title{Create a Combined Diagnostic Plot for mvardlurt Objects}
\description{
Creates a 2x2 panel of diagnostic plots for the multivariate ARDL unit root
test including residuals vs fitted, Q-Q plot, residuals over time, and ACF.
}
\usage{
\method{autoplot}{mvardlurt}(x, \dots)
}
\arguments{
\item{x}{An object of class \code{"mvardlurt"}.}

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
autoplot.mvardlurt(result)
}
\seealso{
\code{\link{mvardlurt}}, \code{\link{plot.mvardlurt}}
}
\author{
Merwan Roudane
}
