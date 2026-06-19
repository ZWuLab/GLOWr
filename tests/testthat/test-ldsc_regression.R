# ==============================================================================
# Tests for ldsc_regression()
# ==============================================================================
#
# Synthetic-data checks of the weighted two-pass LD Score regression: linear
# recovery of intercept/slope, the confounding-ratio extremes (pure confounding
# -> ~1, pure polygenicity -> ~0), positive jackknife SEs, and input validation.
# Numerical fidelity to the research script is checked separately by the
# Phase-2 reproduction of ldsc-results.csv (not a unit test; needs cohort data).

# Deterministic synthetic generator: E[chi2] = intercept + slope * ld, with
# small additive noise so the weighted fit recovers the line tightly.
.make_ldsc_data <- function(n = 8000L, intercept = 1.1, slope = 5e-4,
                            noise = 0.02, seed = 1L) {
  set.seed(seed)
  ld  <- stats::rgamma(n, shape = 2, scale = 50)          # mean ~100
  chr <- rep(1:22, length.out = n)
  chi2 <- intercept + slope * ld + stats::rnorm(n, 0, noise)
  list(chi2 = pmax(chi2, 1e-6), ld = ld, chr = chr)
}

test_that("ldsc_regression recovers the intercept and slope of the line", {
  d <- .make_ldsc_data(intercept = 1.12, slope = 6e-4)
  fit <- ldsc_regression(d$chi2, d$ld, d$chr)
  expect_s3_class(fit, "data.frame")
  expect_equal(nrow(fit), 1L)
  expect_equal(fit$intercept, 1.12, tolerance = 0.02)
  expect_equal(fit$slope,     6e-4, tolerance = 0.20)   # relative
  expect_equal(fit$n_variants, length(d$chi2))
})

test_that("confounding_ratio -> ~1 under pure confounding, ~0 under pure polygenicity", {
  # Pure confounding: a flat lift (slope 0) -> all inflation is the intercept.
  dc <- .make_ldsc_data(intercept = 1.4, slope = 0, seed = 2L)
  rc <- ldsc_regression(dc$chi2, dc$ld, dc$chr)
  expect_gt(rc$confounding_ratio, 0.85)
  expect_equal(rc$slope, 0, tolerance = 5e-4)

  # Pure polygenicity: intercept ~1, inflation rides the slope.
  dp <- .make_ldsc_data(intercept = 1.0, slope = 3e-3, seed = 3L)
  rp <- ldsc_regression(dp$chi2, dp$ld, dp$chr)
  expect_lt(rp$confounding_ratio, 0.20)
  expect_gt(rp$slope, 0)
})

test_that("jackknife standard errors are finite and positive", {
  d <- .make_ldsc_data()
  fit <- ldsc_regression(d$chi2, d$ld, d$chr)
  for (se in c(fit$intercept_se, fit$slope_se, fit$ratio_se)) {
    expect_true(is.finite(se) && se > 0)
  }
})

test_that("two_pass = FALSE differs from the two-pass fit but stays sensible", {
  d <- .make_ldsc_data(intercept = 1.2, slope = 8e-4)
  f2 <- ldsc_regression(d$chi2, d$ld, d$chr, two_pass = TRUE)
  f1 <- ldsc_regression(d$chi2, d$ld, d$chr, two_pass = FALSE)
  expect_equal(f1$intercept, 1.2, tolerance = 0.03)
  expect_false(isTRUE(all.equal(f1$intercept, f2$intercept)))  # weights differ
})

test_that("winsorize caps chi-square inside the fit", {
  d <- .make_ldsc_data(intercept = 1.1, slope = 5e-4)
  d$chi2[1:20] <- 5000                       # extreme outliers
  hi  <- ldsc_regression(d$chi2, d$ld, d$chr, winsorize = 80)
  expect_lt(hi$mean_chi2, 5)                 # mean uses winsorised chi2
})

test_that("input validation: too few pairs / too few chr blocks error", {
  d <- .make_ldsc_data(n = 50L)
  expect_error(ldsc_regression(d$chi2, d$ld, d$chr), "at least 100")
  d2 <- .make_ldsc_data(n = 500L)
  d2$chr <- rep(1L, length(d2$chr))          # single block
  expect_error(ldsc_regression(d2$chi2, d2$ld, d2$chr), "chr blocks")
})

test_that("non-finite (chi2, ld) pairs are dropped", {
  d <- .make_ldsc_data()
  d$chi2[c(1, 5)] <- NA; d$ld[c(2, 6)] <- Inf
  fit <- ldsc_regression(d$chi2, d$ld, d$chr)
  expect_equal(fit$n_variants, length(d$chi2) - 4L)
})
