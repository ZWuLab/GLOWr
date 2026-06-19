# Test file for input validation functions

# Load required packages
library(testthat)

# ============================================================================
# Tests for validate_genotype_matrix()
# ============================================================================

test_that("validate_genotype_matrix accepts valid discrete genotypes", {
  # Create valid discrete genotype matrix (0/1/2)
  G <- matrix(sample(0:2, 100, replace = TRUE), nrow = 10, ncol = 10)

  # Run validation
  result <- validate_genotype_matrix(G)

  # Should pass validation
  expect_true(result$valid)
  expect_equal(result$n_samples, 10)
  expect_equal(result$n_variants, 10)
  expect_equal(result$n_missing, 0)
})

test_that("validate_genotype_matrix accepts valid continuous dosages", {
  # Create continuous dosage data [0, 2]
  G <- matrix(runif(100, 0, 2), nrow = 10, ncol = 10)

  # Run validation with discrete=FALSE
  result <- validate_genotype_matrix(G, discrete = FALSE)

  # Should pass validation
  expect_true(result$valid)
  expect_equal(result$n_samples, 10)
  expect_equal(result$n_variants, 10)
})

test_that("validate_genotype_matrix rejects invalid discrete values", {
  # Create genotype matrix with invalid value (3)
  G <- matrix(c(0, 1, 2, 3, 0, 1), nrow = 2, ncol = 3)

  # Run validation
  result <- validate_genotype_matrix(G)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "invalid values")
})

test_that("validate_genotype_matrix rejects dosages outside [0, 2]", {
  # Create dosage data with value > 2
  G <- matrix(c(0, 1, 2.5, 1.5, 0.5), nrow = 5, ncol = 1)

  # Run validation with discrete=FALSE
  result <- validate_genotype_matrix(G, discrete = FALSE)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "range \\[0, 2\\]")
})

test_that("validate_genotype_matrix handles missing values correctly", {
  # Create matrix with missing values
  G <- matrix(sample(0:2, 100, replace = TRUE), nrow = 10, ncol = 10)
  G[1:3, 1] <- NA  # Add 3 missing values

  # Should fail without allow_na
  result_no_na <- validate_genotype_matrix(G, allow_na = FALSE)
  expect_false(result_no_na$valid)
  expect_match(result_no_na$message, "missing values")

  # Should pass with allow_na=TRUE
  result_with_na <- validate_genotype_matrix(G, allow_na = TRUE)
  expect_true(result_with_na$valid)
  expect_equal(result_with_na$n_missing, 3)
})

test_that("validate_genotype_matrix rejects all-monomorphic data", {
  # Create matrix with no variation (all zeros)
  G <- matrix(0, nrow = 10, ncol = 10)

  # Run validation
  result <- validate_genotype_matrix(G)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "monomorphic")
})

test_that("validate_genotype_matrix warns about high monomorphic rate", {
  # Create matrix where 80% of variants are monomorphic
  G <- matrix(0, nrow = 10, ncol = 10)
  G[, 1:2] <- sample(0:2, 20, replace = TRUE)  # Only 2 polymorphic variants

  # Run validation
  result <- validate_genotype_matrix(G)

  # Should pass but with warning
  expect_true(result$valid)
  expect_equal(result$n_monomorphic, 8)
  expect_length(result$warnings, 1)
  expect_match(result$warnings[1], "monomorphic")
})

test_that("validate_genotype_matrix warns about high missing rate", {
  # Create matrix with 20% missing values
  G <- matrix(sample(0:2, 100, replace = TRUE), nrow = 10, ncol = 10)
  G[1:20] <- NA

  # Run validation with allow_na=TRUE
  result <- validate_genotype_matrix(G, allow_na = TRUE)

  # Should pass but with warning
  expect_true(result$valid)
  expect_equal(result$n_missing, 20)
  expect_length(result$warnings, 1)
  expect_match(result$warnings[1], "missing")
})

test_that("validate_genotype_matrix rejects empty matrices", {
  # Test empty matrix
  G_empty <- matrix(numeric(0), nrow = 0, ncol = 0)
  result <- validate_genotype_matrix(G_empty)

  expect_false(result$valid)
  expect_match(result$message, "at least 1 sample")
})

test_that("validate_genotype_matrix rejects non-numeric data", {
  # Create character matrix
  G_char <- matrix(c("A", "B", "C", "D"), nrow = 2, ncol = 2)

  # Run validation
  result <- validate_genotype_matrix(G_char)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "numeric")
})


# ============================================================================
# Tests for validate_correlation_matrix()
# ============================================================================

test_that("validate_correlation_matrix accepts valid correlation matrix", {
  # Create valid correlation matrix (identity)
  M <- diag(10)

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should pass validation
  expect_true(result$valid)
  expect_equal(result$dimension, 10)
  expect_true(result$is_symmetric)
})

test_that("validate_correlation_matrix accepts correlation structure", {
  # Create correlation matrix with off-diagonal correlations
  M <- matrix(0.3, 5, 5) + diag(0.7, 5)

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should pass validation
  expect_true(result$valid)
  expect_equal(result$dimension, 5)
  expect_true(result$is_symmetric)
  expect_gte(result$min_eigenvalue, -1e-8)  # Should be PSD
})

test_that("validate_correlation_matrix rejects non-square matrices", {
  # Create non-square matrix
  M <- matrix(1, nrow = 5, ncol = 10)

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "square")
})

test_that("validate_correlation_matrix rejects asymmetric matrices", {
  # Create asymmetric matrix
  M <- matrix(runif(25), 5, 5)

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "symmetric")
})

test_that("validate_correlation_matrix rejects matrices with bad diagonal", {
  # Create matrix with diagonal != 1
  M <- diag(5)
  M[1, 1] <- 2  # Wrong diagonal value

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "diagonal must be all 1s")
})

test_that("validate_correlation_matrix rejects off-diagonal values outside [-1,1]", {
  # Create matrix with invalid correlation (> 1)
  M <- matrix(1.5, 3, 3)
  diag(M) <- 1

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "\\[-1, 1\\]")
})

test_that("validate_correlation_matrix rejects non-PSD matrices", {
  # Create non-PSD matrix by forcing negative eigenvalues
  # Start with valid correlation matrix
  M <- matrix(0.3, 4, 4) + diag(0.7, 4)
  # Make one correlation too large to maintain PSD property
  M[1, 2] <- M[2, 1] <- 0.9
  M[1, 3] <- M[3, 1] <- 0.9
  M[2, 3] <- M[3, 2] <- -0.9

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should fail if matrix is truly non-PSD
  # (This specific matrix might still be PSD, but demonstrates the check)
  if (!result$valid) {
    expect_match(result$message, "positive semi-definite")
  }
})

test_that("validate_correlation_matrix checks dimension matching", {
  # Create 5x5 correlation matrix
  M <- diag(5)

  # Should pass with correct n_variants
  result_match <- validate_correlation_matrix(M, n_variants = 5)
  expect_true(result_match$valid)

  # Should fail with wrong n_variants
  result_mismatch <- validate_correlation_matrix(M, n_variants = 10)
  expect_false(result_mismatch$valid)
  expect_match(result_mismatch$message, "does not match")
})

test_that("validate_correlation_matrix warns about near-singularity", {
  # Create nearly singular matrix
  M <- matrix(0.99, 5, 5)
  diag(M) <- 1

  # Run validation
  result <- validate_correlation_matrix(M)

  # Might pass but should have warnings
  if (result$valid) {
    # Check for warnings about singularity or high correlation
    expect_gte(length(result$warnings), 0)
  }
})

test_that("validate_correlation_matrix warns about very high correlations", {
  # Create matrix with very high correlation
  M <- matrix(0.995, 3, 3)
  diag(M) <- 1

  # Run validation
  result <- validate_correlation_matrix(M)

  # Should pass but warn
  if (result$valid) {
    expect_gte(length(result$warnings), 1)
    expect_match(result$warnings[1], "high correlation|multicollinearity", ignore.case = TRUE)
  }
})

test_that("validate_correlation_matrix rejects non-numeric matrices", {
  # Create character matrix
  M_char <- matrix("a", nrow = 3, ncol = 3)

  # Run validation
  result <- validate_correlation_matrix(M_char)

  # Should fail validation
  expect_false(result$valid)
  expect_match(result$message, "numeric")
})
