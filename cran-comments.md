## mvardlurt 1.1.1

CRAN version: 1.1.0.

### Changes

* Documentation only; computations and results are identical to 1.1.0.
* Help page of mvardlurt() (Details) and README: the lag orders p and q used in the bootstrap are stated to be those of the test regression and to be held fixed in every bootstrap replication.
* Help page of mvardlurt(): rejection frequencies of the bootstrap tests from a Monte Carlo experiment (n = 100, 200 replications) are reported for fixed and AIC-selected lag orders.

### Test environments

* Local: R 4.3.3 on Linux (x86_64), R CMD check --as-cran.
* win-builder: R-devel (2026-10-30 r90805 ucrt), Status OK.

### R CMD check results

R CMD check --as-cran on R 4.3.3 (Linux, LANG=C.UTF-8, _R_CHECK_SYSTEM_CLOCK_=0): 0 errors, 0 warnings, 1 NOTE. The NOTE is from the network check: the build environment could not reach doi.org or github.com (proxy refused the connection), so the URLs were reported as possibly invalid; the DOI was verified separately through the CrossRef API.

win-builder, R-devel: Status OK for the same sources with Date 2026-11-01.
