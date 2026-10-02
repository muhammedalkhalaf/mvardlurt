\name{mvardlurt}
\alias{mvardlurt}
\title{Multivariate ARDL Unit Root Test}
\description{
Implements the multivariate autoregressive distributed lag (ARDL) unit root
test of Sam, McNown, Goh and Goh (2025). The test augments the ADF
regression with the lagged level, the current difference and lagged
differences of one or more covariates, so that a possible cointegrating
relationship between the series under test and the covariates is taken into
account. Critical values and p-values are obtained by the residual bootstrap
of Section 4.2 of the paper.
}
\usage{
mvardlurt(y, x, case = 3L, maxlag = 10L, ic = "aic",
          fixlag = NULL, nboot = 999L, level = 0.05,
          seed = NULL, boot = TRUE, reps = NULL)
}
\arguments{
\item{y}{A numeric vector or time series. The series under test.}

\item{x}{A numeric vector, time series or numeric matrix with one column per
  covariate (the I(1) forcing variables of the paper).}

\item{case}{Integer. Deterministic specification:
  \itemize{
    \item \code{1}: No deterministic terms
    \item \code{3}: Intercept only (default)
    \item \code{5}: Intercept and linear trend
  }}

\item{maxlag}{Integer. Maximum ARDL order for AIC/BIC selection. Default is
  10. Must be between 1 and 12.}

\item{ic}{Character. Information criterion for lag selection: \code{"aic"}
  (default) or \code{"bic"}.}

\item{fixlag}{Optional numeric vector of length 2, the fixed orders
  \code{c(p, q)} of the ARDL(p, q) model (see Details). Both must be at
  least 1. If provided, overrides automatic lag selection.}

\item{nboot}{Integer. Number of bootstrap replications. Default is 999.
  Minimum is 99.}

\item{level}{Numeric. Significance level of the two tests used for the case
  classification (0 to 1). Default is 0.05.}

\item{seed}{Optional integer. If supplied, the random number generator is
  seeded with it for the bootstrap and the caller's RNG state is restored on
  exit. Default \code{NULL} leaves the RNG untouched.}

\item{boot}{Logical. Whether to compute bootstrap critical values and
  p-values. Default is \code{TRUE}.}

\item{reps}{Deprecated name for \code{nboot}, kept for backward
  compatibility.}
}
\details{
The test regression is equation (1) of the paper,

\deqn{\Delta y_t = c_1 + c_2 t + b_1 y_{t-1} + \beta_2' x_{t-1}
  + \sum_{i=1}^{p-1} \phi_i \Delta y_{t-i}
  + \sum_{j=1}^{q-1} \Phi_j' \Delta x_{t-j} + \omega' \Delta x_t + u_t,}

estimated by OLS on \eqn{t = \max(p, q) + 1, \ldots, n}. The orders follow
equation (1) of the paper: an ARDL(p, q) model contains \eqn{p - 1} lagged
differences of y and \eqn{q - 1} lagged differences of x in addition to the
contemporaneous difference of x, so p and q are at least 1 and ARDL(1, 1)
contains no lagged differences at all. The paper is not consistent on this
point: equations (1) and (24) sum to \eqn{p - 1} and \eqn{q - 1}, whereas
Section 4.1 describes the ARDL(1, 1) data generating process as having one
lagged difference of each variable and Table 4 reports an ARDL(0, 2) model
with no lagged difference of y. This package follows equation (1); the
paper's ARDL(0, 2) corresponds to \code{fixlag = c(1, 3)} here, and
\code{fixlag = c(p, q)} in versions before 1.1.0 corresponds to
\code{fixlag = c(p + 1, q + 1)} in 1.1.0 (apart from the contemporaneous
difference of x, which earlier versions omitted).

Two statistics are computed (Section 3.3, equations (11) and (12) of the
paper):
\itemize{
  \item the t statistic on \eqn{y_{t-1}} for \eqn{H_0: b_1 = 0} against
    \eqn{H_1: b_1 < 0} (left tailed);
  \item the F (Wald) statistic for the joint hypothesis
    \eqn{H_0: \beta_2 = 0} on the lagged levels of all covariates (right
    tailed). With one covariate this equals the squared t statistic on
    \eqn{x_{t-1}}.
}

Their null distributions depend on nuisance parameters, so both are
bootstrapped separately with the respective null imposed: the restricted
regression (without \eqn{y_{t-1}} for the t test, without \eqn{x_{t-1}} for
the F test) is estimated, its recentred residuals are resampled with
replacement, a bootstrap series \eqn{y^*} is generated recursively from the
restricted estimated equation with the observed x held fixed, and the
unrestricted regression is re-estimated on \eqn{(y^*, x)}. The
contemporaneous difference of x, which equation (24) of the paper omits
although equation (1) has it, is included in the restricted regressions and
in the bootstrap, following equation (1). The initial values
\eqn{y^*_1, \ldots, y^*_{\max(p, q)}} are set to the observed y, so the
recursion starts at \eqn{t = \max(p, q) + 1} with the observed history.
Critical values
are the percentiles of the \code{nboot} bootstrap statistics, and the
bootstrap p-values are the proportions of bootstrap statistics at least as
extreme as the observed ones.

The decision at significance level \code{level} classifies the series
according to the four cases of Section 3.2 of the paper:
\itemize{
  \item \strong{Case I} (t not rejected, F not rejected): nonstationary
    process, no cointegration; y is I(1).
  \item \strong{Case II} (t rejected, F not rejected): stationary process;
    y is I(0).
  \item \strong{Case III} (t not rejected, F rejected): degenerate lagged
    dependent variable; y is I(2).
  \item \strong{Case IV} (both rejected): nonstationary process with
    cointegration between y and x; y is I(1).
}
The Monte Carlo experiments of the paper use a 5 percent level for both
tests, which is the default here. In the empirical application of the paper
the F statistic is reported as significant at the 10 percent level and
cointegration is concluded (Section 6); set \code{level = 0.10} to
reproduce that convention. Before version 1.1.0 \code{level} was a
confidence level (default 0.95); values of 0.5 or more are now rejected.

Lag orders are either fixed through \code{fixlag} or selected by AIC or BIC
over p = 1, ..., \code{maxlag} and q = 1, ..., \code{maxlag}. All candidate
models are estimated on the same observations (t = \code{maxlag} + 1, ...,
n) so that the criteria are comparable;
the selected model is then re-estimated on its own full sample.
}
\value{
An object of class \code{"mvardlurt"} containing:
\item{tstat}{t statistic for \eqn{H_0: b_1 = 0}}
\item{fstat}{F statistic for \eqn{H_0: \beta_2 = 0}}
\item{t_pval, f_pval}{Bootstrap p-values (NA when \code{boot = FALSE})}
\item{t_cv, f_cv}{Bootstrap critical values at the 10, 5, 2.5 and 1 percent
  levels}
\item{t_boot, f_boot}{The bootstrap distributions}
\item{b1, b1_se}{Estimate and standard error of \eqn{b_1}}
\item{beta2, beta2_se}{Estimates and standard errors of \eqn{\beta_2} (one
  per covariate)}
\item{lr_mult}{Long-run multipliers \eqn{-\beta_2 / b_1}}
\item{opt_p, opt_q}{Orders of the ARDL(p, q) model used}
\item{case, casename}{Deterministic case and its description}
\item{nboot}{Number of bootstrap replications}
\item{nobs}{Number of observations in the regression}
\item{aic, bic, r_squared, adj_r_squared}{Fit statistics of the model}
\item{ic_table}{Matrix of information criteria for all (p, q) candidates
  (NULL with \code{fixlag})}
\item{decision}{List with the decision at \code{level}: \code{reject_t},
  \code{reject_f}, \code{case_num}, \code{case_result}, \code{integration}
  (the implied order of integration of y) and the significance codes}
\item{model}{The fitted \code{lm} object, for inspection (its call refers
  to an internal data frame; to refit use
  \code{lm(formula(object$model), data = object$data)})}
\item{data}{The data frame of the test regression}
\item{y, x}{The data used (x as a matrix)}
\item{residuals}{Residuals from the fitted model}
}
\references{
Sam, C. Y., McNown, R., Goh, S. K. and Goh, K. L. (2025). A multivariate
autoregressive distributed lag unit root test. \emph{Studies in Economics
and Econometrics}, 49(1), 17-33.
\doi{10.1080/03796205.2024.2439101}
}
\examples{
# Cointegrated data: y is I(1) through x (Case IV)
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))
y <- 0.5 * x + rnorm(n, sd = 0.5)

result <- mvardlurt(y, x, case = 3, fixlag = c(2, 2), nboot = 199)
print(result)
summary(result)

# Automatic lag selection with two covariates
X <- cbind(x, z = cumsum(rnorm(n)))
result2 <- mvardlurt(y, X, maxlag = 3, nboot = 199)
result2$fstat
}
\author{Muhammad Alkhalaf}
