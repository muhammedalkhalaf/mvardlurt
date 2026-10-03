# mvardlurt: Multivariate ARDL Unit Root Test

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mvardlurt)](https://CRAN.R-project.org/package=mvardlurt)
<!-- badges: end -->

R implementation of the multivariate autoregressive distributed lag (ARDL) unit root test of Sam, McNown, Goh and Goh (2025).

## Overview

The test regression is equation (1) of the paper,

```
dy[t] = c1 + c2 t + b1 y[t-1] + beta2' x[t-1] + sum_{i=1}^{p-1} phi_i dy[t-i]
        + sum_{j=1}^{q-1} Phi_j' dx[t-j] + omega' dx[t] + u[t]
```

Two statistics are computed: the t statistic on `y[t-1]` (H0: b1 = 0) and the F statistic on the lagged levels of the covariates (H0: beta2 = 0). Their null distributions depend on nuisance parameters, so both are bootstrapped separately with the null imposed (residual bootstrap of Section 4.2 of the paper): the restricted regression is estimated, its recentred residuals are resampled, a bootstrap series y* is generated recursively from the restricted equation with the observed x held fixed, and the unrestricted regression is re-estimated on (y*, x). Bootstrap critical values at 10, 5, 2.5 and 1 percent and bootstrap p-values are reported. The lag orders selected (or fixed) on the observed data are held fixed in every bootstrap replication, as in Section 4.2 of the paper; the lag selection is not repeated on the bootstrap samples (see `?mvardlurt`).

An ARDL(p, q) model contains p - 1 lagged differences of y and q - 1 lagged differences of x in addition to the contemporaneous difference of x (p, q >= 1), following equation (1) of the paper. The paper is not consistent on this point (Section 4.1 and Table 4 count lagged differences directly); the paper's ARDL(0, 2) corresponds to `fixlag = c(1, 3)` here. The contemporaneous difference of x, omitted in equation (24) of the paper, is kept in the restricted regressions and in the bootstrap.

## Installation

Install from CRAN:

```r
install.packages("mvardlurt")
```

## Usage

```r
library(mvardlurt)

# Generate example data with cointegration
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))  # I(1) process
y <- 0.5 * x + rnorm(n, sd = 0.5)  # Cointegrated with x

# Run the test (lags selected by AIC, 999 bootstrap replications)
result <- mvardlurt(y, x, case = 3, nboot = 999, seed = 1)
print(result)
summary(result)

# Several covariates: the F test is joint on all lagged levels
X <- cbind(x, z = cumsum(rnorm(n)))
result2 <- mvardlurt(y, X, fixlag = c(2, 2), nboot = 999)

# Diagnostic plots
plot(result)
```

## The Four Cases (Section 3.2 of the paper)

The two tests are judged at the significance level `level` (default 0.05, as in the paper's simulations; in the paper's empirical application the F statistic is reported as significant at the 10 percent level and cointegration is concluded, Section 6). Before 1.1.0 `level` was a confidence level (default 0.95).

| Case | t-test (b1 = 0) | F-test (beta2 = 0) | y is | Interpretation |
|------|--------|--------|------|----------------|
| I | Not rejected | Not rejected | I(1) | Nonstationary, no cointegration |
| II | Rejected | Not rejected | I(0) | Stationary process |
| III | Not rejected | Rejected | I(2) | Degenerate lagged dependent variable |
| IV | Rejected | Rejected | I(1) | Nonstationary, cointegrated with x |

## Deterministic Cases

- **Case 1**: No deterministic terms
- **Case 3**: Intercept only (default)
- **Case 5**: Intercept and linear trend

## References

Sam, C. Y., McNown, R., Goh, S. K. and Goh, K. L. (2025). A multivariate autoregressive distributed lag unit root test. *Studies in Economics and Econometrics*, 49(1), 17-33. [doi:10.1080/03796205.2024.2439101](https://doi.org/10.1080/03796205.2024.2439101)

## License

GPL-3

## Author

Muhammad Alkhalaf
