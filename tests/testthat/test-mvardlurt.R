# Tests for mvardlurt package

# Data from the paper's second DGP (equations (16) to (21)) with b1 = -0.2,
# b2 = 0.3, r = 0.3, c1 = 0.5: y is I(1) through cointegration with x.
gen_coint <- function(n = 120, seed = 1) {
  set.seed(seed)
  ex <- rnorm(n)
  ey <- 0.3 * ex + rnorm(n)
  x <- cumsum(ex)
  y <- numeric(n)
  for (t in 3:n) {
    dy <- 0.5 - 0.2 * y[t - 1] + 0.3 * x[t - 1] -
      0.3 * (y[t - 1] - y[t - 2]) - 0.3 * (x[t - 1] - x[t - 2]) + ey[t]
    y[t] <- y[t - 1] + dy
  }
  list(y = y, x = x)
}

test_that("statistics equal lm() and anova() computed by hand (ARDL(2, 2))", {
  g <- gen_coint()
  r <- mvardlurt(g$y, g$x, case = 3, fixlag = c(2, 2), boot = FALSE)

  dy <- diff(g$y); dx <- diff(g$x); N <- length(dy); s <- 2:N
  d <- data.frame(dy = dy[s], Ly = g$y[s], Lx = g$x[s], L1dy = dy[s - 1],
                  L1dx = dx[s - 1], dx0 = dx[s])
  m <- lm(dy ~ Ly + Lx + L1dy + L1dx + dx0, d)
  cs <- summary(m)$coefficients
  expect_equal(r$nobs, nrow(d))
  expect_equal(r$tstat, unname(cs["Ly", "t value"]), tolerance = 1e-10)
  expect_equal(r$b1, unname(cs["Ly", "Estimate"]), tolerance = 1e-10)
  expect_equal(unname(r$beta2), unname(cs["Lx", "Estimate"]), tolerance = 1e-10)
  expect_equal(r$fstat, unname(cs["Lx", "t value"])^2, tolerance = 1e-10)
  m0 <- lm(dy ~ Ly + L1dy + L1dx + dx0, d)
  expect_equal(r$fstat, anova(m0, m)$F[2], tolerance = 1e-10)
  expect_equal(unname(r$lr_mult), unname(-cs["Lx", 1] / cs["Ly", 1]))
  expect_true(is.na(r$t_pval) && is.na(r$f_pval))
  # the returned data frame refits the stored model
  refit <- lm(formula(r$model), data = r$data)
  expect_equal(coef(refit), coef(r$model))
  expect_equal(nrow(r$data), r$nobs)
})

test_that("ARDL(1, 1) with trend and no deterministics matches lm()", {
  g <- gen_coint(n = 80, seed = 2)
  dy <- diff(g$y); dx <- diff(g$x); N <- length(dy); s <- 1:N
  d <- data.frame(dy = dy[s], Ly = g$y[s], Lx = g$x[s], dx0 = dx[s],
                  trend = s + 1)
  r5 <- mvardlurt(g$y, g$x, case = 5, fixlag = c(1, 1), boot = FALSE)
  cs5 <- summary(lm(dy ~ Ly + Lx + dx0 + trend, d))$coefficients
  expect_equal(r5$nobs, nrow(d))
  expect_equal(r5$tstat, unname(cs5["Ly", 3]), tolerance = 1e-10)
  expect_equal(r5$fstat, unname(cs5["Lx", 3])^2, tolerance = 1e-10)
  r1 <- mvardlurt(g$y, g$x, case = 1, fixlag = c(1, 1), boot = FALSE)
  cs1 <- summary(lm(dy ~ Ly + Lx + dx0 - 1, d))$coefficients
  expect_equal(r1$tstat, unname(cs1["Ly", 3]), tolerance = 1e-10)
  expect_equal(r1$fstat, unname(cs1["Lx", 3])^2, tolerance = 1e-10)
})

test_that("joint F test with two covariates equals anova()", {
  g <- gen_coint(n = 100, seed = 3)
  set.seed(4)
  X <- cbind(a = g$x, b = cumsum(rnorm(100)))
  r <- mvardlurt(g$y, X, fixlag = c(2, 2), boot = FALSE)
  dy <- diff(g$y); N <- length(dy); s <- 2:N
  dA <- diff(X[, 1]); dB <- diff(X[, 2])
  d <- data.frame(dy = dy[s], Ly = g$y[s], La = X[s, 1], Lb = X[s, 2],
                  L1dy = dy[s - 1], L1da = dA[s - 1], L1db = dB[s - 1],
                  da = dA[s], db = dB[s])
  m <- lm(dy ~ ., d)
  m0 <- lm(dy ~ . - La - Lb, d)
  expect_equal(r$fstat, anova(m0, m)$F[2], tolerance = 1e-10)
  expect_equal(r$tstat, unname(summary(m)$coefficients["Ly", 3]),
               tolerance = 1e-10)
  expect_named(r$beta2, c("a", "b"))
  expect_named(coef(r), c("b1", "beta2.a", "beta2.b", "lr_mult.a",
                          "lr_mult.b"))
})

test_that("fast bootstrap statistics agree with lm() on the original data", {
  g <- gen_coint()
  X <- matrix(g$x, ncol = 1, dimnames = list(NULL, "x"))
  d <- mvardlurt:::.ardl_design(g$y, X, 3, 2, 3)
  s <- mvardlurt:::.ols_stats(d$Y, d$Z, d$ty, d$tx)
  m <- lm(mvardlurt:::.ardl_formula(d), data = d$data)
  cs <- summary(m)$coefficients
  expect_equal(unname(s["t"]), unname(cs["L.y", 3]), tolerance = 1e-10)
  expect_equal(unname(s["F"]), unname(cs["L.x", 3])^2, tolerance = 1e-10)
})

test_that("bootstrap recursion (filter) matches an explicit loop", {
  g <- gen_coint()
  y <- g$y
  X <- matrix(g$x, ncol = 1, dimnames = list(NULL, "x"))
  d <- mvardlurt:::.ardl_design(y, X, 3, 2, 3)
  keep <- setdiff(colnames(d$Z), "L.x")
  Zr <- d$Z[, keep]
  b <- lm.fit(Zr, d$Y)$coefficients
  set.seed(5)
  e <- rnorm(length(d$Y))
  xc <- c("(Intercept)", "d.x", "L1.dx")
  w <- as.numeric(Zr[, xc] %*% b[xc]) + e
  ys <- y
  for (i in seq_along(d$idx)) {
    t <- d$idx[i]
    ys[t] <- ys[t - 1] + b["L.y"] * ys[t - 1] +
      b["L1.dy"] * (ys[t - 1] - ys[t - 2]) +
      b["L2.dy"] * (ys[t - 2] - ys[t - 3]) + w[i]
  }
  phi <- b[c("L1.dy", "L2.dy")]
  a <- as.numeric(c(1 + b["L.y"] + phi[1], diff(phi), -phi[2]))
  st <- d$idx[1]
  ys2 <- y
  ys2[d$idx] <- as.numeric(stats::filter(w, a, method = "recursive",
                                         init = y[(st - 1):(st - 3)]))
  expect_equal(ys, ys2, tolerance = 1e-12)
})

test_that("case labelling follows Section 3.2 of the paper", {
  t_boot <- seq(-4, 1, length.out = 200)
  f_boot <- seq(0, 10, length.out = 200)
  t_cv <- quantile(t_boot, c(0.10, 0.05, 0.025, 0.01), names = FALSE)
  f_cv <- quantile(f_boot, c(0.90, 0.95, 0.975, 0.99), names = FALSE)
  names(t_cv) <- names(f_cv) <- c("cv10", "cv05", "cv025", "cv01")
  mk <- function(tstat, fstat) {
    mvardlurt:::.make_decision(tstat, fstat, t_boot, f_boot, t_cv, f_cv,
                               0.05, TRUE)
  }
  d <- mk(-1, 1)
  expect_equal(d$case_num, 1L); expect_equal(d$integration, "I(1)")
  expect_false(d$reject_t); expect_false(d$reject_f)
  d <- mk(-5, 1)
  expect_equal(d$case_num, 2L); expect_equal(d$integration, "I(0)")
  d <- mk(-1, 20)
  expect_equal(d$case_num, 3L); expect_equal(d$integration, "I(2)")
  d <- mk(-5, 20)
  expect_equal(d$case_num, 4L); expect_equal(d$integration, "I(1)")
  expect_match(d$case_result, "cointegration")
  expect_equal(d$t_pval, 0); expect_equal(d$f_pval, 0)
  # level is used: a t statistic between the 5% and 10% quantiles rejects
  # at 10% but not at 5%
  tmid <- mean(t_cv[c("cv10", "cv05")])
  expect_false(mk(tmid, 1)$reject_t)
  expect_true(mvardlurt:::.make_decision(tmid, 1, t_boot, f_boot, t_cv, f_cv,
                                         0.10, TRUE)$reject_t)
  expect_equal(mk(tmid, 1)$t_level, "10%")
})

test_that("bootstrap critical values depend on the data", {
  g <- gen_coint(n = 100, seed = 1)
  r1 <- mvardlurt(g$y, g$x, fixlag = c(2, 2), nboot = 199, seed = 7)
  set.seed(9)
  y2 <- cumsum(rnorm(100, sd = 3))
  x2 <- cumsum(rnorm(100))
  r2 <- mvardlurt(y2, x2, fixlag = c(2, 2), nboot = 199, seed = 7)
  expect_false(isTRUE(all.equal(r1$t_cv, r2$t_cv)))
  expect_false(isTRUE(all.equal(r1$f_cv, r2$f_cv)))
  expect_length(r1$t_boot, 199)
  expect_true(all(diff(r1$t_cv) <= 0))
  expect_true(all(diff(r1$f_cv) >= 0))
  expect_true(r1$t_pval >= 0 && r1$t_pval <= 1)
  expect_equal(r1$t_pval, mean(r1$t_boot <= r1$tstat))
  expect_equal(r1$f_pval, mean(r1$f_boot >= r1$fstat))
})

test_that("the cointegrated example is classified as Case IV", {
  g <- gen_coint(n = 120, seed = 1)
  r <- mvardlurt(g$y, g$x, case = 3, fixlag = c(2, 2), nboot = 199, seed = 7)
  expect_equal(r$decision$case_num, 4L)
  expect_equal(r$decision$integration, "I(1)")
  expect_true(r$decision$reject_t && r$decision$reject_f)
})

test_that("seed handling restores the caller's RNG state", {
  g <- gen_coint(n = 80, seed = 2)
  set.seed(100); a <- runif(3)
  set.seed(100)
  invisible(mvardlurt(g$y, g$x, fixlag = c(1, 1), nboot = 99, seed = 3))
  b <- runif(3)
  expect_identical(a, b)
  # seed = NULL does not reset the stream either
  set.seed(100)
  invisible(mvardlurt(g$y, g$x, fixlag = c(1, 1), nboot = 99))
  c1 <- runif(3)
  set.seed(100)
  invisible(mvardlurt(g$y, g$x, fixlag = c(1, 1), nboot = 99))
  c2 <- runif(3)
  expect_identical(c1, c2)
  expect_false(identical(a, c1))
  # reproducibility with a seed
  r1 <- mvardlurt(g$y, g$x, fixlag = c(1, 1), nboot = 99, seed = 11)
  r2 <- mvardlurt(g$y, g$x, fixlag = c(1, 1), nboot = 99, seed = 11)
  expect_identical(r1$t_cv, r2$t_cv)
  expect_identical(r1$f_cv, r2$f_cv)
})

test_that("lag selection uses a common sample and returns the IC table", {
  g <- gen_coint(n = 100, seed = 6)
  r <- mvardlurt(g$y, g$x, maxlag = 3, ic = "bic", boot = FALSE)
  expect_equal(dim(r$ic_table), c(3, 3))
  expect_false(r$manual_lag)
  X <- matrix(g$x, ncol = 1, dimnames = list(NULL, "x"))
  # every candidate is estimated on t = 4, ..., n
  for (p in 1:3) for (q in 1:3) {
    d <- mvardlurt:::.ardl_design(g$y, X, p, q, 3, start = 4)
    expect_equal(length(d$Y), 97)
    expect_equal(r$ic_table[p, q], BIC(lm(mvardlurt:::.ardl_formula(d),
                                          data = d$data)))
  }
  expect_equal(r$ic_table[r$opt_p, r$opt_q], min(r$ic_table))
  # the selected model is re-estimated on its own sample
  expect_equal(r$nobs, 100 - max(r$opt_p, r$opt_q))
})

test_that("input validation works", {
  g <- gen_coint(n = 100, seed = 1)
  expect_error(mvardlurt(g$y, g$x, case = 2), "'case' must be 1")
  expect_error(mvardlurt(g$y, g$x, ic = "hqic"), "'ic' must be")
  expect_error(mvardlurt(g$y, g$x, maxlag = 15), "'maxlag' must be")
  expect_error(mvardlurt(g$y, g$x, nboot = 50), "'nboot' must be at least")
  expect_error(mvardlurt(g$y, g$x, fixlag = c(0, 1)), "at least 1")
  expect_error(mvardlurt(g$y, g$x, level = 1.5), "'level' must be")
  expect_error(mvardlurt(g$y, g$x, level = 0.95), "did you mean 1 - level")
  expect_error(mvardlurt(g$y[1:20], g$x[1:20]), "Too few observations")
  expect_error(mvardlurt(g$y, g$x[1:50]), "same length")
  expect_warning(r <- mvardlurt(g$y, g$x, fixlag = c(1, 1), reps = 99),
                 "deprecated")
  expect_equal(r$nboot, 99L)
})

test_that("methods work", {
  g <- gen_coint(n = 80, seed = 2)
  r <- mvardlurt(g$y, g$x, fixlag = c(2, 2), nboot = 99, seed = 1)
  expect_s3_class(r, "mvardlurt")
  expect_equal(r$casename, "Intercept Only")
  expect_output(print(r), "CASE")
  expect_output(print(r), "y is I\\(")
  expect_output(summary(r), "Four Cases")
  expect_type(residuals(r), "double")
  expect_equal(length(fitted(r)), length(residuals(r)))
  r0 <- mvardlurt(g$y, g$x, fixlag = c(2, 2), boot = FALSE)
  expect_true(all(is.na(r0$t_cv)))
  expect_output(print(r0), "t-statistic")
})
