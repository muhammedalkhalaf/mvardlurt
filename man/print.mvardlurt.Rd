\name{print.mvardlurt}
\alias{print.mvardlurt}
\alias{summary.mvardlurt}
\alias{coef.mvardlurt}
\alias{residuals.mvardlurt}
\alias{fitted.mvardlurt}
\title{Methods for mvardlurt Objects}
\description{
Print, summary, and accessor methods for objects of class \code{"mvardlurt"}.
}
\usage{
\method{print}{mvardlurt}(x, \dots)

\method{summary}{mvardlurt}(object, \dots)

\method{coef}{mvardlurt}(object, \dots)

\method{residuals}{mvardlurt}(object, \dots)

\method{fitted}{mvardlurt}(object, \dots)
}
\arguments{
\item{x}{An object of class \code{"mvardlurt"}.}

\item{object}{An object of class \code{"mvardlurt"}.}

\item{\dots}{Additional arguments (ignored).}
}
\value{
\code{print} and \code{summary} invisibly return \code{x}/\code{object}.

\code{coef} returns a named numeric vector with \code{pi}, \code{delta},
and \code{lr_mult}.

\code{residuals} returns the numeric vector of residuals from the fitted model.

\code{fitted} returns the numeric vector of fitted values.
}
\examples{
set.seed(123)
n <- 100
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.5)

result <- mvardlurt(y, x, reps = 100)

# Print method
print(result)

# Summary method
summary(result)

# Extract coefficients
coef(result)

# Extract residuals
head(residuals(result))

# Extract fitted values
head(fitted(result))
}
\seealso{
\code{\link{mvardlurt}}
}
\author{Muhammad Alkhalaf}
