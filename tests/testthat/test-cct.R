# This file tests the cct_test() function (Cauchy Combination Test).

test_that("cct_test handles regular p-values correctly", {
  # Test with typical p-values
  p_vals <- c(0.01, 0.05, 0.10, 0.20)
  result <- cct_test(p_vals)

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("cct", "pval_cct"))
  expect_type(result$cct, "double")
  expect_type(result$pval_cct, "double")

  # Check p-value is valid
  expect_true(result$pval_cct >= 0 && result$pval_cct <= 1)

  # Manually compute CCT statistic for verification
  expected_cct <- mean(tan(pi * (0.5 - p_vals)))
  expect_equal(result$cct, expected_cct, tolerance = 1e-12)

  # Verify p-value from Cauchy distribution
  expected_pval <- pcauchy(expected_cct, lower.tail = FALSE)
  expect_equal(result$pval_cct, expected_pval, tolerance = 1e-12)
})

test_that("cct_test handles very small p-values correctly", {
  # Test with very small p-values below threshold
  p_vals_small <- c(1e-20, 1e-18, 0.05)
  result <- cct_test(p_vals_small, thr.smallp = 1e-15)

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("cct", "pval_cct"))

  # Check p-value is valid
  expect_true(result$pval_cct >= 0 && result$pval_cct <= 1)

  # Manually compute with mixed formula
  is.small <- p_vals_small < 1e-15
  cct_stat <- p_vals_small
  cct_stat[!is.small] <- tan((0.5 - p_vals_small[!is.small]) * pi)
  cct_stat[is.small] <- 1 / p_vals_small[is.small] / pi
  expected_cct <- mean(cct_stat)

  expect_equal(result$cct, expected_cct, tolerance = 1e-10)

  # Very small p-values should give small combined p-value
  expect_true(result$pval_cct < 0.01)
})

test_that("cct_test handles edge case: p-value exactly at small threshold", {
  # Test p-value exactly at threshold
  p_vals <- c(1e-15, 0.05, 0.10)
  result <- cct_test(p_vals, thr.smallp = 1e-15)

  # At threshold, should be treated as regular p-value (not small)
  is.small <- p_vals < 1e-15
  expect_equal(sum(is.small), 0)

  # Should use standard formula for all
  expected_cct <- mean(tan(pi * (0.5 - p_vals)))
  expect_equal(result$cct, expected_cct, tolerance = 1e-12)
})

test_that("cct_test truncates large p-values correctly", {
  # Test with large p-values above threshold
  p_vals_large <- c(0.001, 0.95, 0.99)
  result <- cct_test(p_vals_large, thr.largp = 0.9)

  # Large p-values should be truncated to 0.9
  p_vals_truncated <- p_vals_large
  p_vals_truncated[p_vals_truncated > 0.9] <- 0.9

  # Compute expected CCT with truncated values
  expected_cct <- mean(tan(pi * (0.5 - p_vals_truncated)))
  expect_equal(result$cct, expected_cct, tolerance = 1e-12)

  # Check p-value is valid
  expect_true(result$pval_cct >= 0 && result$pval_cct <= 1)
})

test_that("cct_test handles edge case: p = 0 (treated as very small)", {
  # P-value of 0 should be handled as very small
  p_vals <- c(0, 0.05, 0.10)
  result <- cct_test(p_vals, thr.smallp = 1e-15)

  # p=0 is less than threshold, so uses special formula
  is.small <- p_vals < 1e-15
  expect_true(is.small[1])

  # Note: 1/(0*pi) = Inf, which should give very strong evidence
  # The function should handle this gracefully
  expect_type(result$cct, "double")
  expect_type(result$pval_cct, "double")

  # With p=0, CCT should be very large (or Inf)
  expect_true(is.infinite(result$cct) || result$cct > 1e10)

  # Combined p-value should be very small (or 0)
  expect_true(result$pval_cct < 1e-10 || result$pval_cct == 0)
})

test_that("cct_test handles edge case: p = 1 (truncated)", {
  # P-value of 1 should be truncated to thr.largp
  p_vals <- c(0.01, 0.05, 1.0)
  result <- cct_test(p_vals, thr.largp = 0.9)

  # p=1 should be truncated to 0.9
  p_vals_truncated <- c(0.01, 0.05, 0.9)
  expected_cct <- mean(tan(pi * (0.5 - p_vals_truncated)))
  expect_equal(result$cct, expected_cct, tolerance = 1e-12)

  # Check p-value is valid
  expect_true(result$pval_cct >= 0 && result$pval_cct <= 1)
})

test_that("cct_test handles single p-value", {
  # Single p-value
  p_val <- 0.05
  result <- cct_test(p_val)

  # CCT statistic should be the transformed single p-value
  expected_cct <- tan(pi * (0.5 - p_val))
  expect_equal(result$cct, expected_cct, tolerance = 1e-12)

  # P-value should match
  expected_pval <- pcauchy(expected_cct, lower.tail = FALSE)
  expect_equal(result$pval_cct, expected_pval, tolerance = 1e-12)
})

test_that("cct_test handles all identical p-values", {
  # All same p-values
  p_vals <- rep(0.05, 5)
  result <- cct_test(p_vals)

  # CCT should be same as single transformation
  expected_cct <- tan(pi * (0.5 - 0.05))
  expect_equal(result$cct, expected_cct, tolerance = 1e-12)

  # P-value should be valid
  expect_true(result$pval_cct >= 0 && result$pval_cct <= 1)
})

test_that("cct_test respects custom threshold parameters", {
  # Test with custom thresholds
  p_vals <- c(1e-10, 0.05, 0.85)
  result1 <- cct_test(p_vals, thr.largp = 0.8, thr.smallp = 1e-8)
  result2 <- cct_test(p_vals, thr.largp = 0.9, thr.smallp = 1e-12)

  # Results should differ due to different thresholds
  expect_false(isTRUE(all.equal(result1$cct, result2$cct)))

  # With thr.smallp = 1e-8, p=1e-10 is treated as very small
  # With thr.smallp = 1e-12, p=1e-10 is regular
  # So CCT statistics should differ

  # Both should produce valid p-values
  expect_true(result1$pval_cct >= 0 && result1$pval_cct <= 1)
  expect_true(result2$pval_cct >= 0 && result2$pval_cct <= 1)
})

test_that("cct_test produces smaller p-value with more significant inputs", {
  # More significant p-values should give smaller combined p-value
  p_vals_sig <- c(0.001, 0.005, 0.01)
  p_vals_nonsig <- c(0.3, 0.4, 0.5)

  result_sig <- cct_test(p_vals_sig)
  result_nonsig <- cct_test(p_vals_nonsig)

  # Significant inputs should give smaller combined p-value
  expect_true(result_sig$pval_cct < result_nonsig$pval_cct)
})

test_that("cct_test mixed small and regular p-values computes correctly", {
  # Mix of very small and regular p-values
  p_vals <- c(1e-20, 1e-16, 0.01, 0.05)
  result <- cct_test(p_vals, thr.smallp = 1e-15)

  # Manually compute with mixed formula
  is.small <- p_vals < 1e-15
  expect_equal(sum(is.small), 2)  # First two are small

  cct_stat <- p_vals
  cct_stat[!is.small] <- tan((0.5 - p_vals[!is.small]) * pi)
  cct_stat[is.small] <- 1 / p_vals[is.small] / pi
  expected_cct <- mean(cct_stat)

  expect_equal(result$cct, expected_cct, tolerance = 1e-10)

  # P-value should be valid
  expect_true(result$pval_cct >= 0 && result$pval_cct <= 1)
})

test_that("cct_test handles p-values just above and below small threshold", {
  # P-values straddling the threshold
  thr <- 1e-15
  p_vals <- c(thr * 0.9, thr * 1.1, 0.05)  # One below, one above threshold
  result <- cct_test(p_vals, thr.smallp = thr)

  # First p-value should be small, second should not
  is.small <- p_vals < thr
  expect_true(is.small[1])
  expect_false(is.small[2])

  # Should use mixed formula
  cct_stat <- p_vals
  cct_stat[!is.small] <- tan((0.5 - p_vals[!is.small]) * pi)
  cct_stat[is.small] <- 1 / p_vals[is.small] / pi
  expected_cct <- mean(cct_stat)

  expect_equal(result$cct, expected_cct, tolerance = 1e-10)
})

test_that("cct_test numerical stability with extreme values", {
  # Test numerical stability with very extreme values
  p_vals <- c(1e-100, 1e-50, 1e-10, 0.5)
  result <- cct_test(p_vals, thr.smallp = 1e-15)

  # First two values are very small (< 1e-15), third is 1e-10 which is > 1e-15
  is.small <- p_vals < 1e-15
  expect_equal(sum(is.small), 2)

  # Function should handle extreme values without error
  expect_type(result$cct, "double")
  expect_type(result$pval_cct, "double")

  # CCT statistic might be very large or Inf with such extreme p-values
  # But should not be NaN
  expect_false(is.nan(result$cct))
  expect_false(is.nan(result$pval_cct))

  # P-value should be very small or 0
  expect_true(result$pval_cct < 0.001)
})

test_that("cct_test validates against manual Cauchy calculation", {
  # Test that pcauchy is used correctly
  p_vals <- c(0.1, 0.2, 0.3)
  result <- cct_test(p_vals)

  # Manual calculation
  cct_stat <- mean(tan(pi * (0.5 - p_vals)))
  pval_manual <- 1 - pcauchy(cct_stat)  # upper.tail = 1 - lower.tail
  pval_function <- pcauchy(cct_stat, lower.tail = FALSE)

  # Both methods should give same result
  expect_equal(pval_manual, pval_function, tolerance = 1e-12)
  expect_equal(result$pval_cct, pval_function, tolerance = 1e-12)
})
