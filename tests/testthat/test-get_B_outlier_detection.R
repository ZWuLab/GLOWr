########## Tests for Outlier Detection Module ##########
#
# Comprehensive tests for detect_B_outliers() and related functions
#
# Test coverage:
# 1. Statistical outlier detection (Cook's distance)
# 2. Biological outlier detection (implausible effects)
# 3. Combined detection (both methods)
# 4. Edge cases and error handling
# 5. S3 methods (print, summary)
#

library(testthat)

test_that("detect_B_outliers detects statistical outliers correctly", {
  # Create synthetic data with clear high-influence point
  set.seed(123)
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45)
  # Most values follow a reasonable pattern
  values <- c(0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.05)
  # Add one high-influence outlier
  values[6] <- 2.0  # Much higher than expected

  result <- detect_B_outliers(MAF, values, method = "statistical", verbose = 0)

  # Should detect at least one statistical outlier
  expect_true(length(result$statistical) >= 1)
  expect_true(inherits(result, "glow_outliers"))
  expect_equal(result$method_used, "statistical")
  expect_true(6 %in% result$indices)  # The outlier we created
})

test_that("detect_B_outliers detects biological outliers correctly", {
  # Create data with common variant having unrealistic large effect
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45)
  # Common variant (MAF=0.4) has implausibly large effect
  values <- c(0.5, 0.3, 0.2, 0.15, 0.1, 15, 0.05)

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # Should detect the common variant with large effect
  expect_true(6 %in% result$biological)
  expect_true(6 %in% result$indices)
  expect_equal(result$method_used, "biological")
  expect_true(inherits(result, "glow_outliers"))
})

test_that("rare variants with large effects are NOT flagged as biological outliers", {
  # Rare variants can legitimately have large effects
  MAF <- c(0.001, 0.01, 0.02, 0.05, 0.1, 0.2)
  values <- c(20, 15, 10, 8, 0.3, 0.15)  # Large effects for rare variants

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # Should not flag rare variants (MAF <= 0.05) with large effects
  # Only common variants (MAF > 0.05 by default) should be flagged
  expect_false(1 %in% result$biological)  # MAF=0.001
  expect_false(2 %in% result$biological)  # MAF=0.01
  expect_false(3 %in% result$biological)  # MAF=0.02
  expect_false(4 %in% result$biological)  # MAF=0.05 (boundary)
})

test_that("combined detection method works correctly", {
  # Create data with both types of outliers
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45)
  values <- c(0.5, 0.3, 0.2, 0.15, 0.1, 12, 0.08, 0.05)
  # Index 6 should be biological outlier (common + large effect)
  # We'll also make it a statistical outlier by design

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # Should detect outliers
  expect_true(length(result$indices) >= 1)
  expect_equal(result$method_used, "both")

  # Check that both detection lists are present
  expect_true(!is.null(result$statistical))
  expect_true(!is.null(result$biological))
})

test_that("custom biological rules work correctly", {
  MAF <- c(0.01, 0.05, 0.12, 0.2, 0.3, 0.4)
  values <- c(0.5, 0.3, 6, 0.15, 0.1, 0.08)

  # Custom rule: MAF > 0.1 and |value| > 5
  custom_rules <- list(
    common_large = list(maf_min = 0.1, effect_max = 5)
  )

  result <- detect_B_outliers(MAF, values, method = "biological",
                              biological_rules = custom_rules, verbose = 0)

  # Should flag index 3 (MAF=0.12, value=6)
  expect_true(3 %in% result$biological)
  # Should not flag indices with MAF <= 0.1
  expect_false(1 %in% result$biological)
  expect_false(2 %in% result$biological)
})

test_that("cook_threshold parameter affects detection", {
  set.seed(456)
  MAF <- seq(0.01, 0.4, length.out = 20)
  values <- 0.5 / MAF + rnorm(20, 0, 0.1)
  # Add one outlier
  values[10] <- values[10] * 3

  # Stricter threshold (larger value) should detect fewer outliers
  result_strict <- detect_B_outliers(MAF, values, method = "statistical",
                                     cook_threshold = 8, verbose = 0)

  # Looser threshold (smaller value) should detect more outliers
  result_loose <- detect_B_outliers(MAF, values, method = "statistical",
                                    cook_threshold = 2, verbose = 0)

  # Loose threshold should detect at least as many as strict
  expect_true(length(result_loose$statistical) >= length(result_strict$statistical))
})

test_that("detect_B_outliers handles edge case: empty data", {
  expect_error(
    detect_B_outliers(numeric(0), numeric(0)),
    "cannot be empty"
  )
})

test_that("detect_B_outliers handles edge case: single observation", {
  # Single observation should work but likely detect no outliers
  result <- detect_B_outliers(MAF = 0.1, values = 0.5,
                              method = "statistical", verbose = 0)

  expect_true(inherits(result, "glow_outliers"))
  expect_equal(length(result$indices), 0)  # Can't compute Cook's D with n=1
})

test_that("detect_B_outliers handles case where all are outliers", {
  # All common variants with large effects
  MAF <- rep(0.3, 5)
  values <- rep(20, 5)

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # All should be flagged
  expect_equal(length(result$biological), 5)
  expect_equal(result$indices, 1:5)
})

test_that("detect_B_outliers handles case with no outliers", {
  # Well-behaved data
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3)
  values <- 0.3 / sqrt(MAF)  # Reasonable pattern

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # May or may not detect outliers depending on Cook's D, but shouldn't error
  expect_true(inherits(result, "glow_outliers"))
})

test_that("detect_B_outliers validates input correctly", {
  MAF <- c(0.01, 0.1, 0.2)
  values <- c(0.5, 0.3, 0.2)

  # Mismatched lengths
  expect_error(
    detect_B_outliers(MAF, c(0.5, 0.3)),
    "same length"
  )

  # NA values
  expect_error(
    detect_B_outliers(c(0.1, NA, 0.2), values),
    "cannot contain NA"
  )

  expect_error(
    detect_B_outliers(MAF, c(0.5, NA, 0.2)),
    "cannot contain NA"
  )

  # Invalid MAF range
  expect_error(
    detect_B_outliers(c(0, 0.1, 0.2), values),
    "must be in the range"
  )

  expect_error(
    detect_B_outliers(c(0.1, 0.2, 0.6), values),
    "must be in the range"
  )

  # Invalid method
  expect_error(
    detect_B_outliers(MAF, values, method = "invalid"),
    "must be one of"
  )

  # Invalid cook_threshold
  expect_error(
    detect_B_outliers(MAF, values, cook_threshold = -1),
    "must be positive"
  )

  # Inf values in MAF
  expect_error(
    detect_B_outliers(c(0.1, Inf, 0.2), values),
    "Inf"
  )

  expect_error(
    detect_B_outliers(c(0.1, -Inf, 0.2), values),
    "Inf"
  )

  # Inf values in values
  expect_error(
    detect_B_outliers(MAF, c(0.5, Inf, 0.2)),
    "Inf"
  )

  expect_error(
    detect_B_outliers(MAF, c(0.5, -Inf, 0.2)),
    "Inf"
  )

  # Non-numeric MAF
  expect_error(
    detect_B_outliers(c("0.1", "0.2", "0.3"), values),
    "numeric"
  )

  # Non-numeric values
  expect_error(
    detect_B_outliers(MAF, c("0.5", "0.3", "0.2")),
    "numeric"
  )
})

test_that("detect_B_outliers handles negative values correctly", {
  # BETA values can be negative
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  values <- c(-0.5, 0.3, 0.2, 0.15, -15, 0.08)  # Negative large effect at index 5

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # Should detect negative large effect for common variant
  expect_true(5 %in% result$biological)
})

test_that("outlier descriptions are generated correctly", {
  MAF <- c(0.01, 0.1, 0.3, 0.4)
  values <- c(0.5, 0.3, 12, 0.15)

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # Should have descriptions for detected outliers
  expect_true(length(result$reasons) > 0)

  # Check description format
  if (length(result$indices) > 0) {
    idx <- result$indices[1]
    desc <- result$reasons[[as.character(idx)]]
    expect_true(is.character(desc))
    expect_true(grepl("MAF=", desc))
    expect_true(grepl("value=", desc))
  }
})

test_that("print.glow_outliers produces output without errors", {
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  values <- c(0.5, 0.3, 0.2, 0.15, 15, 0.08)

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # Print should not error
  expect_output(print(result), "GLOW Outlier Detection")
  expect_output(print(result), "Total variants")
  expect_output(print(result), "Total outliers")
})

test_that("print.glow_outliers handles max_show parameter", {
  # Create data with many outliers
  MAF <- rep(0.3, 15)
  values <- rep(15, 15)

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # With max_show = 5, should show only 5 and indicate more exist
  output <- capture.output(print(result, max_show = 5))
  combined <- paste(output, collapse = " ")

  expect_true(any(grepl("and.*more", output, ignore.case = TRUE)))
})

test_that("summary.glow_outliers produces comprehensive output", {
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  values <- c(0.5, 0.3, 0.2, 0.15, 15, 0.08)

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # Summary should produce output and return data frame
  expect_output(summary_result <- summary(result), "Summary")
  expect_output(summary(result), "Detection method")

  # If outliers detected, should return data frame
  if (length(result$indices) > 0) {
    expect_true(is.data.frame(summary_result) || is.null(summary_result))
  }
})

test_that("summary.glow_outliers handles case with no outliers", {
  # Well-behaved data
  MAF <- c(0.01, 0.05, 0.1, 0.2)
  values <- c(0.5, 0.35, 0.25, 0.18)

  result <- detect_B_outliers(MAF, values, method = "biological", verbose = 0)

  # Summary should handle no outliers gracefully
  expect_output(summary(result), "Summary")
})

test_that("verbose parameter controls output correctly", {
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  values <- c(0.5, 0.3, 0.2, 0.15, 15, 0.08)

  # verbose = 0 should produce no messages
  expect_silent(detect_B_outliers(MAF, values, method = "both", verbose = 0))

  # verbose >= 1 should produce messages if outliers found
  # (we use expect_message or expect_output)
  expect_message(
    detect_B_outliers(MAF, values, method = "both", verbose = 1),
    "outlier"
  )

  # verbose >= 2 should produce more detailed messages
  expect_message(
    detect_B_outliers(MAF, values, method = "both", verbose = 2),
    "Computing|Detecting"
  )
})

test_that("result structure is complete and correct", {
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  values <- c(0.5, 0.3, 0.2, 0.15, 15, 0.08)

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # Check all required components exist
  expect_true(!is.null(result$indices))
  expect_true(!is.null(result$statistical))
  expect_true(!is.null(result$biological))
  expect_true(!is.null(result$reasons))
  expect_true(!is.null(result$outlier_MAF))
  expect_true(!is.null(result$outlier_values))
  expect_true(!is.null(result$n_total))
  expect_true(!is.null(result$method_used))

  # Check types
  expect_true(is.integer(result$indices) || is.numeric(result$indices))
  expect_true(is.list(result$reasons))
  expect_equal(length(result$outlier_MAF), length(result$indices))
  expect_equal(length(result$outlier_values), length(result$indices))
  expect_equal(result$n_total, length(MAF))

  # Check class
  expect_true(inherits(result, "glow_outliers"))
})

test_that("overlapping outliers are handled correctly", {
  # Create data where same point is both statistical and biological outlier
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.35, 0.4)
  values <- c(0.5, 0.3, 0.2, 0.15, 20, 0.08)  # Index 5 is extreme

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # If index 5 is detected by both methods
  if (5 %in% result$statistical && 5 %in% result$biological) {
    # Should appear only once in combined indices
    expect_equal(sum(result$indices == 5), 1)

    # Description should mention both
    desc <- result$reasons[["5"]]
    expect_true(grepl("Statistical.*Biological|Both", desc, ignore.case = TRUE))
  }
})

test_that("detect_B_outliers works with h-squared values", {
  # h-squared values (typically smaller than BETA^2)
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  h2_values <- c(0.001, 0.0005, 0.0003, 0.0002, 0.5, 0.00015)

  # Should work the same way
  result <- detect_B_outliers(MAF, h2_values, method = "biological",
                              biological_rules = list(
                                common_large = list(maf_min = 0.05, effect_max = 0.01)
                              ), verbose = 0)

  expect_true(inherits(result, "glow_outliers"))
  # Index 5 has large h2 for common variant
  expect_true(5 %in% result$biological)
})

test_that("Cook's distance calculation handles perfect fit edge case", {
  # Data that fits perfectly might have issues with Cook's distance
  MAF <- c(0.1, 0.2, 0.3)
  values <- c(1, 2, 3)  # Perfect linear relationship

  # Should handle gracefully without errors
  result <- detect_B_outliers(MAF, values, method = "statistical", verbose = 0)
  expect_true(inherits(result, "glow_outliers"))
})

test_that("result object stores outlier data efficiently", {
  MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  values <- c(0.5, 0.3, 0.2, 0.15, 15, 0.08)

  result <- detect_B_outliers(MAF, values, method = "both", verbose = 0)

  # Should store only outlier values (memory efficient)
  expect_equal(length(result$outlier_MAF), length(result$indices))
  expect_equal(length(result$outlier_values), length(result$indices))

  # Should store total number of variants
  expect_equal(result$n_total, length(MAF))

  # Outlier values should match the original MAF/values at outlier indices
  if (length(result$indices) > 0) {
    expect_equal(result$outlier_MAF, MAF[result$indices])
    expect_equal(result$outlier_values, values[result$indices])
  }
})
