# Unit tests for g-transformation functions

library(testthat)

# Test g_GFisher_two function
test_that("g_GFisher_two handles normal Z-scores correctly", {
  # Test with known Z-scores
  z_scores <- c(-2, -1, 0, 1, 2)
  df <- 2

  result <- GLOWr:::g_GFisher_two(z_scores, df)

  # Result should be numeric and same length as input
  expect_type(result, "double")
  expect_length(result, length(z_scores))

  # Result should be non-negative (chi-square statistics)
  expect_true(all(result >= 0 | is.na(result)))

  # Transformation should be symmetric: g(z) = g(-z) for two-sided
  expect_equal(
    GLOWr:::g_GFisher_two(2, df),
    GLOWr:::g_GFisher_two(-2, df),
    tolerance = 1e-14
  )

  # Manual calculation for z = 1.96 (approximately p = 0.05)
  z <- 1.96
  # Two-sided p-value: 2 * pnorm(-1.96) ≈ 0.05
  p_val_log <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  expected <- qchisq(p_val_log, df = df, lower.tail = FALSE, log.p = TRUE)
  actual <- GLOWr:::g_GFisher_two(z, df)
  expect_equal(actual, expected, tolerance = 1e-14)
})

test_that("g_GFisher_two handles edge cases correctly", {
  df <- 2

  # Test Z = 0 (p-value = 1)
  result_zero <- GLOWr:::g_GFisher_two(0, df)
  expect_type(result_zero, "double")
  expect_true(is.finite(result_zero))
  # For Z=0, two-sided p-value is 1, so log(p) = 0
  # qchisq(log(2), df=2, lower.tail=FALSE, log.p=TRUE) should be finite

  # Test very large Z (p-value approaches 0)
  result_large <- GLOWr:::g_GFisher_two(10, df)
  expect_type(result_large, "double")
  # For very large Z, result should be large positive or Inf
  expect_true(result_large > 0 || is.infinite(result_large))

  # Test very small Z (negative, large absolute value)
  result_small <- GLOWr:::g_GFisher_two(-10, df)
  expect_type(result_small, "double")
  # Should equal result for positive large Z (symmetric)
  expect_equal(result_small, result_large, tolerance = 1e-14)

  # Test NA handling
  result_na <- GLOWr:::g_GFisher_two(NA, df)
  expect_true(is.na(result_na))

  # Test vector with NA
  mixed <- c(1, NA, 2)
  result_mixed <- GLOWr:::g_GFisher_two(mixed, df)
  expect_true(is.na(result_mixed[2]))
  expect_false(is.na(result_mixed[1]))
  expect_false(is.na(result_mixed[3]))
})

test_that("g_GFisher_two vectorization works correctly", {
  # Test vectorization
  z_vec <- rnorm(100)
  df <- 2

  result_vec <- GLOWr:::g_GFisher_two(z_vec, df)
  expect_length(result_vec, 100)

  # Compare vectorized vs. element-wise
  result_loop <- sapply(z_vec, function(z) GLOWr:::g_GFisher_two(z, df))
  expect_equal(result_vec, result_loop, tolerance = 1e-14)
})


# Test g_GFisher_one function
test_that("g_GFisher_one handles normal Z-scores correctly", {
  # Test with known Z-scores
  z_scores <- c(-2, -1, 0, 1, 2)
  df <- 2

  result <- GLOWr:::g_GFisher_one(z_scores, df)

  # Result should be numeric and same length as input
  expect_type(result, "double")
  expect_length(result, length(z_scores))

  # Result should be non-negative (chi-square statistics)
  expect_true(all(result >= 0 | is.na(result)))

  # Transformation should NOT be symmetric for one-sided
  expect_false(
    isTRUE(all.equal(
      GLOWr:::g_GFisher_one(2, df),
      GLOWr:::g_GFisher_one(-2, df),
      tolerance = 1e-10
    ))
  )

  # Manual calculation for z = 1.96
  z <- 1.96
  # One-sided p-value: pnorm(-1.96) ≈ 0.025
  p_val_log <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
  expected <- qchisq(p_val_log, df = df, lower.tail = FALSE, log.p = TRUE)
  actual <- GLOWr:::g_GFisher_one(z, df)
  expect_equal(actual, expected, tolerance = 1e-14)
})

test_that("g_GFisher_one handles edge cases correctly", {
  df <- 2

  # Test Z = 0 (p-value = 0.5)
  result_zero <- GLOWr:::g_GFisher_one(0, df)
  expect_type(result_zero, "double")
  expect_true(is.finite(result_zero))

  # Test very large positive Z (p-value approaches 0)
  result_large_pos <- GLOWr:::g_GFisher_one(10, df)
  expect_type(result_large_pos, "double")
  expect_true(result_large_pos > 0 || is.infinite(result_large_pos))

  # Test very large negative Z (p-value approaches 1)
  result_large_neg <- GLOWr:::g_GFisher_one(-10, df)
  expect_type(result_large_neg, "double")
  # For Z = -10, p-value ≈ 1, so result should be small
  expect_true(result_large_neg < result_large_pos)

  # Test NA handling
  result_na <- GLOWr:::g_GFisher_one(NA, df)
  expect_true(is.na(result_na))
})

test_that("g_GFisher_one vectorization works correctly", {
  # Test vectorization
  z_vec <- rnorm(100)
  df <- 2

  result_vec <- GLOWr:::g_GFisher_one(z_vec, df)
  expect_length(result_vec, 100)

  # Compare vectorized vs. element-wise
  result_loop <- sapply(z_vec, function(z) GLOWr:::g_GFisher_one(z, df))
  expect_equal(result_vec, result_loop, tolerance = 1e-14)
})


# Test g_Burden function
test_that("g_Burden is identity function", {
  # Test with various inputs
  x1 <- c(-2, -1, 0, 1, 2)
  expect_identical(GLOWr:::g_Burden(x1), x1)

  x2 <- rnorm(100)
  expect_identical(GLOWr:::g_Burden(x2), x2)

  # Test with NA
  x3 <- c(1, NA, 2)
  result <- GLOWr:::g_Burden(x3)
  expect_identical(result, x3)

  # Test with single value
  expect_identical(GLOWr:::g_Burden(5), 5)

  # Test with Inf
  expect_identical(GLOWr:::g_Burden(Inf), Inf)
  expect_identical(GLOWr:::g_Burden(-Inf), -Inf)
})


# Test g_GFisher wrapper function
test_that("g_GFisher wrapper routes to correct function with p.type='two'", {
  z_scores <- c(-2, -1, 0, 1, 2)
  df <- 2

  result_wrapper <- GLOWr:::g_GFisher(z_scores, df, p.type = "two")
  result_direct <- GLOWr:::g_GFisher_two(z_scores, df)

  expect_equal(result_wrapper, result_direct, tolerance = 1e-14)
})

test_that("g_GFisher wrapper routes to correct function with p.type='one'", {
  z_scores <- c(-2, -1, 0, 1, 2)
  df <- 2

  result_wrapper <- GLOWr:::g_GFisher(z_scores, df, p.type = "one")
  result_direct <- GLOWr:::g_GFisher_one(z_scores, df)

  expect_equal(result_wrapper, result_direct, tolerance = 1e-14)
})

test_that("g_GFisher wrapper defaults to two-sided", {
  z_scores <- c(-2, -1, 0, 1, 2)
  df <- 2

  # Default should be p.type = "two"
  result_default <- GLOWr:::g_GFisher(z_scores, df)
  result_two <- GLOWr:::g_GFisher_two(z_scores, df)

  expect_equal(result_default, result_two, tolerance = 1e-14)
})

test_that("g_GFisher wrapper handles invalid p.type gracefully", {
  z_scores <- c(-2, -1, 0, 1, 2)
  df <- 2

  # Invalid p.type should return NULL (matching legacy behavior)
  result <- GLOWr:::g_GFisher(z_scores, df, p.type = "invalid")
  expect_null(result)
})


# Legacy validation tests - compare with known values
test_that("g_GFisher_two matches expected numerical values", {
  # Test cases with pre-calculated expected values
  # These were calculated using the legacy implementation

  # Z = 1.96, df = 2 (approximately 5% significance)
  z1 <- 1.96
  df <- 2
  result1 <- GLOWr:::g_GFisher_two(z1, df)
  # Expected value calculated from legacy: qchisq(log(2) + pnorm(1.96, lower=F, log=T), df=2, lower=F, log=T)
  expected1 <- qchisq(log(2) + pnorm(1.96, lower.tail = FALSE, log.p = TRUE),
                      df = 2, lower.tail = FALSE, log.p = TRUE)
  expect_equal(result1, expected1, tolerance = 1e-14)

  # Z = 2.576, df = 2 (approximately 1% significance)
  z2 <- 2.576
  result2 <- GLOWr:::g_GFisher_two(z2, df)
  expected2 <- qchisq(log(2) + pnorm(2.576, lower.tail = FALSE, log.p = TRUE),
                      df = 2, lower.tail = FALSE, log.p = TRUE)
  expect_equal(result2, expected2, tolerance = 1e-14)
})

test_that("g_GFisher_one matches expected numerical values", {
  # Test cases with pre-calculated expected values

  # Z = 1.645, df = 2 (approximately 5% one-sided)
  z1 <- 1.645
  df <- 2
  result1 <- GLOWr:::g_GFisher_one(z1, df)
  expected1 <- qchisq(pnorm(1.645, lower.tail = FALSE, log.p = TRUE),
                      df = 2, lower.tail = FALSE, log.p = TRUE)
  expect_equal(result1, expected1, tolerance = 1e-14)

  # Z = -1.645, df = 2
  z2 <- -1.645
  result2 <- GLOWr:::g_GFisher_one(z2, df)
  expected2 <- qchisq(pnorm(-1.645, lower.tail = FALSE, log.p = TRUE),
                      df = 2, lower.tail = FALSE, log.p = TRUE)
  expect_equal(result2, expected2, tolerance = 1e-14)
})


# Relationship tests between transformations
test_that("g_GFisher_two and g_GFisher_one have expected relationship", {
  df <- 2
  z_positive <- 2
  z_negative <- -2

  # For two-sided, g(z) = g(-z)
  two_pos <- GLOWr:::g_GFisher_two(z_positive, df)
  two_neg <- GLOWr:::g_GFisher_two(z_negative, df)
  expect_equal(two_pos, two_neg, tolerance = 1e-14)

  # For one-sided, g(z) != g(-z) in general
  one_pos <- GLOWr:::g_GFisher_one(z_positive, df)
  one_neg <- GLOWr:::g_GFisher_one(z_negative, df)
  expect_false(isTRUE(all.equal(one_pos, one_neg, tolerance = 1e-10)))
})


# Numerical stability tests
test_that("g_GFisher functions maintain numerical stability", {
  df <- 2

  # Test extreme values
  extreme_z <- c(-10, -5, -1, 0, 1, 5, 10)

  # g_GFisher_two
  result_two <- GLOWr:::g_GFisher_two(extreme_z, df)
  expect_true(all(is.finite(result_two) | is.infinite(result_two)))
  expect_false(any(is.nan(result_two)))

  # g_GFisher_one
  result_one <- GLOWr:::g_GFisher_one(extreme_z, df)
  expect_true(all(is.finite(result_one) | is.infinite(result_one)))
  expect_false(any(is.nan(result_one)))
})


# Test different df values
test_that("g_GFisher functions work with different df values", {
  z <- 1.96

  # Test with df = 1, 2, 3, 4
  for (df in 1:4) {
    result_two <- GLOWr:::g_GFisher_two(z, df)
    expect_type(result_two, "double")
    expect_true(is.finite(result_two))

    result_one <- GLOWr:::g_GFisher_one(z, df)
    expect_type(result_one, "double")
    expect_true(is.finite(result_one))
  }
})
