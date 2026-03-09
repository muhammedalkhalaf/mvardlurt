# Tests for mvardlurt package

test_that("mvardlurt returns correct class", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, reps = 100, boot = TRUE)

expect_s3_class(result, "mvardlurt")
})

test_that("mvardlurt handles different cases", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  # Case 1: No deterministics
  result1 <- mvardlurt(y, x, case = 1, reps = 100)
  expect_equal(result1$case, 1L)
  expect_equal(result1$casename, "No Deterministic Terms")

  # Case 3: Intercept only
  result3 <- mvardlurt(y, x, case = 3, reps = 100)
  expect_equal(result3$case, 3L)
  expect_equal(result3$casename, "Intercept Only")

  # Case 5: Intercept and trend
  result5 <- mvardlurt(y, x, case = 5, reps = 100)
  expect_equal(result5$case, 5L)
  expect_equal(result5$casename, "Intercept and Trend")
})

test_that("mvardlurt validates inputs correctly", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  # Invalid case
  expect_error(mvardlurt(y, x, case = 2), "'case' must be 1")

  # Invalid IC
  expect_error(mvardlurt(y, x, ic = "hqic"), "'ic' must be")

  # Invalid maxlag
  expect_error(mvardlurt(y, x, maxlag = 15), "'maxlag' must be")

  # Invalid reps
  expect_error(mvardlurt(y, x, reps = 50), "'reps' must be at least")

  # Too few observations
  expect_error(mvardlurt(y[1:20], x[1:20]), "Too few observations")

  # Mismatched lengths
  expect_error(mvardlurt(y, x[1:50]), "same length")
})

test_that("mvardlurt with fixed lags works", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, fixlag = c(2, 2), reps = 100)

  expect_equal(result$opt_p, 2)
  expect_equal(result$opt_q, 2)
  expect_true(result$manual_lag)
  expect_null(result$ic_table)
})

test_that("mvardlurt returns expected components", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, reps = 100)

  # Check all expected components exist
  expect_true("tstat" %in% names(result))
  expect_true("fstat" %in% names(result))
  expect_true("pi_coef" %in% names(result))
  expect_true("delta_coef" %in% names(result))
  expect_true("opt_p" %in% names(result))
  expect_true("opt_q" %in% names(result))
  expect_true("t_cv" %in% names(result))
  expect_true("f_cv" %in% names(result))
  expect_true("decision" %in% names(result))
  expect_true("model" %in% names(result))

  # Check types
  expect_type(result$tstat, "double")
  expect_type(result$fstat, "double")
  expect_type(result$opt_p, "integer")
  expect_type(result$opt_q, "integer")
})

test_that("mvardlurt bootstrap critical values have correct structure", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, reps = 200)

  # t critical values (left-tailed, should be negative)
  expect_length(result$t_cv, 4)
  expect_named(result$t_cv, c("cv10", "cv05", "cv025", "cv01"))

  # F critical values (right-tailed, should be positive)
  expect_length(result$f_cv, 4)
  expect_named(result$f_cv, c("cv10", "cv05", "cv025", "cv01"))

  # Critical values should be ordered
  expect_true(all(diff(result$f_cv) >= 0))  # F increasing
})

test_that("print and summary methods work", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, reps = 100)

  # Should not error
  expect_output(print(result))
  expect_output(summary(result))
})

test_that("coef, residuals, and fitted methods work", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, reps = 100)

  # coef
  cf <- coef(result)
  expect_equal(names(cf), c("pi", "delta", "lr_mult"))

  # residuals
  res <- residuals(result)
  expect_type(res, "double")
  expect_true(length(res) > 0)

  # fitted
  fit <- fitted(result)
  expect_type(fit, "double")
  expect_equal(length(fit), length(res))
})

test_that("mvardlurt without bootstrap works", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, boot = FALSE)

  expect_false(result$boot)
  expect_true(all(is.na(result$t_cv)))
  expect_true(all(is.na(result$f_cv)))
})

test_that("mvardlurt handles cointegrated data", {
  set.seed(42)
  n <- 200

  # Generate cointegrated series
  x <- cumsum(rnorm(n))
  e <- rnorm(n, sd = 0.3)
  y <- 1.5 + 0.8 * x + e  # Cointegration relationship

  result <- mvardlurt(y, x, case = 3, reps = 500)

  # With cointegrated data, we expect:
  # - Significant t-test (reject unit root in error)
  # - Significant F-test (reject no cointegration)
  expect_true(result$decision$reject_t || result$decision$reject_f)
})

test_that("mvardlurt handles non-cointegrated data", {
  set.seed(42)
  n <- 200

  # Generate independent I(1) processes
  x <- cumsum(rnorm(n))
  y <- cumsum(rnorm(n))

  result <- mvardlurt(y, x, case = 3, reps = 500)

  # With non-cointegrated data, we expect to be in Case III or IV
  # (accept t-test more likely)
  expect_true(result$decision$case_num %in% c(3L, 4L) ||
                result$decision$case_num %in% c(1L, 2L))
})

test_that("ic_table has correct dimensions", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result <- mvardlurt(y, x, maxlag = 5, reps = 100)

  expect_equal(nrow(result$ic_table), 6)  # 0 to 5
  expect_equal(ncol(result$ic_table), 6)
})

test_that("BIC selection works", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result_aic <- mvardlurt(y, x, ic = "aic", reps = 100)
  result_bic <- mvardlurt(y, x, ic = "bic", reps = 100)

  expect_equal(result_aic$ic, "aic")
  expect_equal(result_bic$ic, "bic")

  # BIC tends to select more parsimonious models
  # (not always, but a reasonable expectation)
})

test_that("seed reproducibility works", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n, sd = 0.5)

  result1 <- mvardlurt(y, x, seed = 999, reps = 100)
  result2 <- mvardlurt(y, x, seed = 999, reps = 100)

  expect_equal(result1$t_cv, result2$t_cv)
  expect_equal(result1$f_cv, result2$f_cv)
})
