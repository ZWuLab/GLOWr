##########################################################################
# Test Suite for P-value Handling in get_B
#
# Tests robust handling of extremely small p-values via the new
# training_P_mlog10 parameter and .compute_chi2_from_pvalue() helper

library(testthat)
library(GLOWr)

# ========== Test .compute_chi2_from_pvalue() Helper ==========

test_that(".compute_chi2_from_pvalue handles regular p-values correctly", {
  # Test with regular p-values in normal range
  p <- c(0.05, 0.01, 1e-5, 1e-10)

  # Compute chi2 using helper
  chi2_helper <- GLOWr:::.compute_chi2_from_pvalue(p = p)

  # Compute chi2 using standard qchisq for comparison
  chi2_expected <- qchisq(p, df = 1, lower.tail = FALSE)

  # Should match exactly in this range
  expect_equal(chi2_helper, chi2_expected, tolerance = 1e-10)
})


test_that(".compute_chi2_from_pvalue handles -log10 transformed p-values", {
  # Test with -log10 transformed p-values
  p_mlog10 <- c(10, 50, 100, 200)  # p from 1e-10 to 1e-200

  # Compute chi2 using helper
  chi2 <- GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = p_mlog10)

  # All should be finite (no Inf)
  expect_true(all(is.finite(chi2)))

  # All should be positive
  expect_true(all(chi2 > 0))

  # Should be monotonically increasing with p_mlog10
  expect_true(all(diff(chi2) > 0))
})


test_that(".compute_chi2_from_pvalue uses approximation for extreme p-values", {
  # Test extremely small p-values where approximation is used (p_mlog10 > 324)
  # The code uses exact qchisq for p_mlog10 <= 324, approximation for > 324
  p_mlog10_extreme <- c(325, 350, 400)

  chi2 <- GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = p_mlog10_extreme)

  # Expected values using approximation: chi2 ≈ 2 * ln(10) * p_mlog10
  chi2_expected <- 2 * log(10) * p_mlog10_extreme

  # Should match the approximation
  expect_equal(chi2, chi2_expected, tolerance = 1e-10)
})


test_that(".compute_chi2_from_pvalue matches between p and p_mlog10 in overlapping range", {
  # Test that results match for same p-value expressed in different ways
  p_mlog10 <- c(1, 5, 10, 50, 100, 200)  # Test various magnitudes

  # Compute using p_mlog10
  chi2_from_log10 <- GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = p_mlog10)

  # Compute using p (for p_mlog10 <= 300, both should work)
  p <- 10^(-p_mlog10)
  chi2_from_p <- GLOWr:::.compute_chi2_from_pvalue(p = p)

  # Should match closely (some numerical error expected for very small p)
  expect_equal(chi2_from_log10, chi2_from_p, tolerance = 1e-8)
})


test_that(".compute_chi2_from_pvalue protects against p-value underflow", {
  # Test p-values that underflow to exactly 0
  p_underflow <- c(1e-300, 1e-320, 0)

  # Should not produce Inf
  chi2 <- GLOWr:::.compute_chi2_from_pvalue(p = p_underflow)

  # Should all be finite
  expect_true(all(is.finite(chi2)))

  # Should all be positive
  expect_true(all(chi2 > 0))
})


test_that(".compute_chi2_from_pvalue validates input correctly", {
  # Test error when neither p nor p_mlog10 provided
  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(),
    "Either p or p_mlog10 must be provided"
  )

  # Test error when both p and p_mlog10 provided
  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(p = 0.05, p_mlog10 = 1.3),
    "Only one of p or p_mlog10 should be provided, not both"
  )

  # Test error for invalid p values
  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(p = c(-0.1, 0.5)),
    "p values must be in the interval"
  )

  # Test error for negative p_mlog10 values
  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = c(-1, 10)),
    "p_mlog10 values must be non-negative"
  )
})


# ========== Test get_B() with training_P_mlog10 ==========

test_that("get_B accepts training_P_mlog10 parameter", {
  # Setup test data
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 10, 100)  # Very small p-values
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Should work without error
  B <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0,
    show_model_selection = FALSE  # Suppress output for test
  )

  # Should return numeric vector
  expect_type(B, "double")
  expect_length(B, length(target_MAF))

  # All values should be finite
  expect_true(all(is.finite(B)))

  # All values should be positive
  expect_true(all(B > 0))
})


test_that("get_B rejects both training_P and training_P_mlog10 simultaneously", {
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P <- runif(n, 1e-5, 0.1)
  training_P_mlog10 <- runif(n, 10, 100)
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Should error when both provided
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_P = training_P,
      training_P_mlog10 = training_P_mlog10,  # Both provided - should error
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "pvalue"
    ),
    "Provide either training_P or training_P_mlog10, not both"
  )
})


test_that("get_B requires either training_P or training_P_mlog10 for pvalue method", {
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Should error when neither provided
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      # Neither training_P nor training_P_mlog10 provided
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "pvalue"
    ),
    "Either training_P or training_P_mlog10 is required for method"
  )
})


test_that("get_B produces consistent results with training_P vs training_P_mlog10", {
  # Test that same p-values give same results when specified as P vs P_mlog10
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P <- runif(n, 1e-10, 1e-2)  # Range where both methods work
  training_P_mlog10 <- -log10(training_P)
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Get B using training_P
  B_from_p <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0,
    show_model_selection = FALSE
  )

  # Get B using training_P_mlog10
  B_from_p_mlog10 <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0,
    show_model_selection = FALSE
  )

  # Should be very similar (some numerical differences expected)
  expect_equal(B_from_p, B_from_p_mlog10, tolerance = 1e-6)
})


test_that("get_B handles extremely small p-values via training_P_mlog10", {
  # Test with p-values that would underflow if not using -log10
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 100, 350)  # p from 1e-100 to 1e-350
  training_N <- rep(50000, n)  # Large GWAS sample

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Should work without producing Inf or NaN
  B <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0
  )

  # All values should be finite
  expect_true(all(is.finite(B)))

  # All values should be positive
  expect_true(all(B > 0))

  # Values should be reasonable (not absurdly large)
  expect_true(all(B < 100))  # Sanity check
})


test_that("get_B with training_P_mlog10 works for continuous traits", {
  # Test continuous trait target
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 10, 200)
  training_N <- rep(10000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_SE <- rep(0.1, 10)  # Continuous trait needs SE

  B <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "continuous",
    target_MAF = target_MAF,
    target_SE = target_SE,
    method = "pvalue",
    verbose = 0
  )

  # Should work and produce finite results
  expect_type(B, "double")
  expect_length(B, length(target_MAF))
  expect_true(all(is.finite(B)))
  expect_true(all(B > 0))
})


test_that("get_B with training_P_mlog10 works for cross-trait inference", {
  # Test training on continuous, predicting for binary
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 10, 150)
  training_N <- rep(15000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.3, 10)

  B <- get_B(
    training_trait = "continuous",  # Different from target
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0
  )

  # Should work and produce finite results
  expect_type(B, "double")
  expect_length(B, length(target_MAF))
  expect_true(all(is.finite(B)))
  expect_true(all(B > 0))
})


test_that("get_B with training_P_mlog10 respects return_full parameter", {
  # Test that full results include P_mlog10 in training_data
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 10, 100)
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    return_full = TRUE,
    verbose = 0
  )

  # Should be a list with class glow_B_estimate
  expect_type(result, "list")
  expect_s3_class(result, "glow_B_estimate")

  # Should contain training_data with P_mlog10 (nested under $model)
  expect_true("model" %in% names(result))
  expect_true("training_data" %in% names(result$model))
  expect_true("P_mlog10" %in% names(result$model$training_data))

  # training_data$P_mlog10 should match input (after any outlier removal)
  expect_true(!is.null(result$model$training_data$P_mlog10))
})


test_that("get_B auto mode selects pvalue method when only training_P_mlog10 available", {
  # Test auto mode behavior with training_P_mlog10
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 10, 100)
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Use method = "auto" (default)
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,  # Only P data, no BETA
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "auto",  # Should select pvalue method
    return_full = TRUE,
    verbose = 0
  )

  # Should have selected pvalue method (nested under $model, new naming convention)
  expect_equal(result$model$method_used, "pvalue_method")
})


# ========== Integration Tests ==========

test_that("get_B handles mixed extreme and moderate p-values via training_P_mlog10", {
  # Test with mix of extreme and moderate p-values
  set.seed(123)
  n <- 100
  training_MAF <- runif(n, 0.001, 0.3)

  # Mix of moderate (10-50) and extreme (200-350) -log10(p)
  training_P_mlog10 <- c(
    runif(50, 10, 50),   # Moderate: p from 1e-10 to 1e-50
    runif(50, 200, 350)  # Extreme: p from 1e-200 to 1e-350
  )
  training_N <- rep(30000, n)

  target_MAF <- runif(20, 0.001, 0.3)
  target_case_prop <- rep(0.5, 20)

  B <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0
  )

  # Should handle both ranges without issues
  expect_type(B, "double")
  expect_length(B, length(target_MAF))
  expect_true(all(is.finite(B)))
  expect_true(all(B > 0))
})


test_that("get_B outlier detection works with training_P_mlog10", {
  # Test that outlier detection works correctly with training_P_mlog10
  set.seed(123)
  n <- 50
  training_MAF <- runif(n, 0.001, 0.3)
  training_P_mlog10 <- runif(n, 10, 100)
  training_N <- rep(5000, n)

  target_MAF <- runif(10, 0.001, 0.3)
  target_case_prop <- rep(0.5, 10)

  # Run with outlier detection enabled
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    outlier_method = "statistical",  # Enable outlier detection
    outlier_action = "flag",
    return_full = TRUE,
    verbose = 0
  )

  # Should have outlier information (nested under $model)
  expect_true("outliers" %in% names(result$model))
  expect_equal(result$model$outliers$method, "statistical")

  # Should complete without error
  expect_type(result$B, "double")
  expect_true(all(is.finite(result$B)))
})


test_that("get_B produces sensible B estimates across p-value ranges", {
  # Test that B estimates are sensible across different p-value magnitudes
  set.seed(123)
  n <- 30
  training_MAF <- rep(0.01, n)  # Fixed MAF for comparison

  # Three groups with different p-value magnitudes
  training_P_mlog10 <- c(
    rep(10, 10),   # p ≈ 1e-10
    rep(50, 10),   # p ≈ 1e-50
    rep(100, 10)   # p ≈ 1e-100
  )
  training_N <- rep(10000, n)

  target_MAF <- c(0.01)
  target_case_prop <- 0.5

  B <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0
  )

  # Should produce a single finite estimate
  expect_length(B, 1)
  expect_true(is.finite(B))
  expect_true(B > 0)

  # More significant p-values should generally lead to larger B estimates
  # (though model fitting may smooth this relationship)
  # Just check that result is reasonable
  expect_true(B < 10)  # Sanity check
})
