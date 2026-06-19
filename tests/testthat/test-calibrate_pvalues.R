# ==============================================================================
# Tests for calibrate_pvalues() and compute_lambda_gc()
# ==============================================================================
#
# Genomic-control p-value calibration (GLOWr Phase 1; design M1). Covers the
# lambda_GC estimator, the general scalar rescaling (closed-form), the two
# methods (lambda_gc vs ldsc_intercept), and the NA / degenerate-factor
# pass-through contract.


# ---- compute_lambda_gc -------------------------------------------------------

test_that("compute_lambda_gc recovers ~1 on a uniform null and the true factor", {
  set.seed(20260530)
  p_null <- stats::runif(20000)
  expect_equal(compute_lambda_gc(p_null), 1, tolerance = 0.05)

  chi2 <- stats::rchisq(20000, df = 1) * 1.3              # inflated by 1.3
  p_inf <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)
  expect_equal(compute_lambda_gc(p_inf), 1.3, tolerance = 0.05)
})

test_that("compute_lambda_gc returns NA below 50 usable p-values and ignores junk", {
  expect_true(is.na(compute_lambda_gc(stats::runif(49))))
  expect_false(is.na(compute_lambda_gc(stats::runif(50))))
  # NA / out-of-(0,1] are dropped before counting and before the median.
  p <- c(stats::runif(60), NA, 1.5, -0.2, 0)
  expect_equal(compute_lambda_gc(p),
               compute_lambda_gc(p[is.finite(p) & p > 0 & p <= 1]))
})


# ---- calibrate_pvalues: the scalar rescaling (closed form) -------------------

test_that("calibrate_pvalues divides chi-square by the factor (closed form)", {
  p <- c(0.5, 0.1, 1e-4, 1e-8, 1e-20)
  f <- 1.5
  expected <- stats::pchisq(
    stats::qchisq(p, df = 1, lower.tail = FALSE) / f,
    df = 1, lower.tail = FALSE)
  res <- calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = f)
  expect_equal(res$p, expected)
  expect_identical(res$calibration_factor, f)
  expect_identical(res$method, "ldsc_intercept")
})

test_that("calibrate_pvalues with factor = 1 is a no-op on valid p-values", {
  p <- c(0.5, 0.1, 1e-4, 1e-8)
  res <- calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = 1)
  expect_equal(res$p, p, tolerance = 1e-12)
})

test_that("calibrate_pvalues is monotone in p (order preserved)", {
  p <- sort(c(0.9, 0.3, 0.05, 1e-3, 1e-6))
  res <- calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = 1.4)
  expect_false(is.unsorted(res$p))
  expect_length(res$p, length(p))
})


# ---- calibrate_pvalues: methods ----------------------------------------------

test_that("lambda_gc estimates the factor from p and maps lambda toward 1", {
  set.seed(11)
  chi2 <- stats::rchisq(20000, df = 1) * 1.3
  p <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)
  res <- calibrate_pvalues(p, method = "lambda_gc")
  expect_equal(res$calibration_factor, 1.3, tolerance = 0.05)   # = compute_lambda_gc(p)
  expect_equal(res$calibration_factor, compute_lambda_gc(p))
  expect_equal(compute_lambda_gc(res$p), 1, tolerance = 0.05)   # calibrated to ~1
})

test_that("ldsc_intercept requires an explicit calibration_factor", {
  expect_error(
    calibrate_pvalues(stats::runif(100), method = "ldsc_intercept"),
    "requires `calibration_factor`")
})


# ---- pass-through contract (NA / degenerate factor) --------------------------

test_that("NA and out-of-(0,1] entries flow through as NA; p = 0 maps to 0", {
  p <- c(0.5, NA, 1.5, 0, 0.01)
  res <- calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = 1.2)
  expect_true(is.na(res$p[2]))      # NA in -> NA out
  expect_true(is.na(res$p[3]))      # 1.5 (> 1) -> NA
  expect_false(is.na(res$p[1]))
  expect_equal(res$p[4], 0)         # p = 0 -> 0
  expect_length(res$p, length(p))
})

test_that("a non-positive or NA factor passes p through unchanged", {
  p <- c(0.5, 0.1, 1e-4)
  expect_identical(
    calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = 0)$p, p)
  expect_identical(
    calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = NA_real_)$p, p)
  # lambda_gc on < 50 p-values: compute_lambda_gc -> NA -> pass-through.
  expect_identical(calibrate_pvalues(p, method = "lambda_gc")$p, p)
})
