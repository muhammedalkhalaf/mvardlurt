## mvardlurt 1.1.1

This is a documentation-only update; the version on CRAN is 1.1.0. There is no numerical change: all statistics, critical values, p-values and decisions are identical to 1.1.0.

* Documentation only; no computation changed and all results are identical
  to 1.1.0. The help page of `mvardlurt()` (Details) and the README now state
  that the lag orders p and q of the bootstrap are those of the test
  regression (fixed through `fixlag` or selected by AIC or BIC on the
  observed data) and are held fixed in every bootstrap replication; the lag
  selection is not repeated on the bootstrap samples. This follows Section
  4.2 of Sam, McNown, Goh and Goh (2025), where the restricted regression
  (24) is estimated with the lag length of the data generating process and
  equation (1) is re-estimated on y* and x. The help page also reports a
  Monte Carlo experiment by the package author (n = 100, y and x
  independent random walks, 200 replications, 199 bootstrap replications,
  5 percent level): rejection frequencies 0.055 to 0.070 with fixed lag
  orders and 0.080 to 0.100 with AIC selection (Monte Carlo standard error
  about 0.02).

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
