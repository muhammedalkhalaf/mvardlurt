# mvardlurt: Multivariate ARDL Unit Root Test

<!-- badges: start -->
[![R-CMD-check](https://github.com/merwanroudane/mvardlurt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/merwanroudane/mvardlurt/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/mvardlurt)](https://CRAN.R-project.org/package=mvardlurt)
<!-- badges: end -->

R implementation of the multivariate autoregressive distributed lag (ARDL) unit root test proposed by Sam, McNown, Goh, and Goh (2024).

## Overview

The `mvardlurt` package implements a unit root test that augments the standard ADF regression with lagged levels of a covariate (independent variable) to improve power, especially when cointegration exists. Bootstrap critical values ensure correct size regardless of nuisance parameters.

## Installation

Install from CRAN:

```r
install.packages("mvardlurt")
```
Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("merwanroudane/mvardlurt")
```

## Usage

```r
library(mvardlurt)

# Generate example data with cointegration
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))  # I(1) process
y <- 0.5 * x + rnorm(n, sd = 0.5)  # Cointegrated with x

# Run the test
result <- mvardlurt(y, x, case = 3, reps = 1000)
print(result)
summary(result)

# Diagnostic plots
plot(result)
```

## The Four-Case Framework

Based on the test results, the package determines one of four cases:

| Case | t-test | F-test | Interpretation |
|------|--------|--------|----------------|
| I | Reject | Reject | Cointegration |
| II | Reject | Accept | Degenerate case 1 (y may be I(0)) |
| III | Accept | Reject | Degenerate case 2 (spurious) |
| IV | Accept | Accept | No cointegration |

## Deterministic Cases

- **Case 1**: No deterministic terms
- **Case 3**: Intercept only (default)
- **Case 5**: Intercept and linear trend

## References

Sam, C. Y., McNown, R., Goh, S. K., & Goh, K. L. (2024). A multivariate autoregressive distributed lag unit root test. *Studies in Economics and Econometrics*, 1-17. [doi:10.1080/03796205.2024.2439101](https://doi.org/10.1080/03796205.2024.2439101)

## License
GPL (>= 3)

## Author

Merwan Roudane (merwanroudane920@gmail.com)
