#' Unit Tests for getZ_marg_score_binary_SPA()
#'
#' This file contains comprehensive unit tests for the SPA-adjusted marginal
#' score statistics function. Tests cover:
#' - Basic functionality with binary traits
#' - Rare variant scenarios where SPA is critical
#' - Unbalanced case-control designs
#' - Edge cases and error handling
#' - Comparison with standard method
#' - Output structure validation
#'
#' File Log (reverse chronological order):
#' - 2026-06-19: Modified by Claude Code (Opus 4.8 / 1M ctx), prompted by ZWu.
#'   Added Test 7b: a deterministic regression test that a zero-variance
#'   (monomorphic) variant is reported as NA with an informative warning rather
#'   than aborting inside SPAtest's saddlepoint routine (the cause of a rare
#'   release-gate test flake). Pairs with the guard added to getZ_marg_score.R.
#' - 2025-10-18: Created by r-developer - Comprehensive tests for SPA function

library(testthat)
library(GLOWr)

# Load validation helper functions
source("helper-validation.R")

# Helper function: test data has X with intercept, but SPA function needs X without intercept
strip_intercept <- function(X) {
  if (ncol(X) > 0 && all(X[, 1] == 1)) {
    return(X[, -1, drop = FALSE])
  }
  return(X)
}

# ==============================================================================
# Test 1: Basic Functionality with Balanced Binary Trait
# ==============================================================================
test_that("getZ_marg_score_binary_SPA works with balanced binary trait", {
  # Load test data
  fixture_path <- system.file("validation/fixtures/test_data_binary.rds", package = "GLOWr")
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }
  test_data <- readRDS(fixture_path)

  # Extract components
  G <- test_data$G
  X <- strip_intercept(test_data$X)  # Remove intercept
  Y <- test_data$Y

  # Run SPA function
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Check output structure
  expect_type(result, "list")
  expect_named(result, c("Zscores", "scores", "M_Z", "M_s", "s0"))

  # Check dimensions
  p <- ncol(G)
  expect_length(result$Zscores, p)
  expect_equal(dim(result$M_Z), c(p, p))
  expect_equal(dim(result$M_s), c(p, p))
  expect_length(result$s0, 1)

  # Check properties
  expect_true(all(is.finite(result$Zscores)))
  expect_equal(result$s0, 1)  # Binary trait: s0 = 1

  # M_Z should be a correlation matrix (diagonal = 1)
  expect_equal(diag(result$M_Z), rep(1, p), tolerance = 1e-10)

  # M_Z should be symmetric
  expect_equal(result$M_Z, t(result$M_Z), tolerance = 1e-10)

  # M_s should be symmetric
  expect_equal(result$M_s, t(result$M_s), tolerance = 1e-10)
})

# ==============================================================================
# Test 2: Rare Variants with SPA Correction
# ==============================================================================
test_that("getZ_marg_score_binary_SPA handles rare variants correctly", {
  # Load rare variant test data
  fixture_path <- system.file("validation/fixtures/test_data_rare.rds", package = "GLOWr")
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }
  test_data <- readRDS(fixture_path)

  G <- test_data$G
  X <- strip_intercept(test_data$X)  # Remove intercept
  Y <- test_data$Y

  # Verify this is indeed rare variant data
  maf <- colMeans(G) / 2
  expect_true(all(maf < 0.05), info = "Test data should contain rare variants")

  # Run SPA function
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Check that function completes without errors
  expect_type(result, "list")
  expect_length(result$Zscores, ncol(G))

  # Z-scores should be finite (SPA prevents inflation)
  expect_true(all(is.finite(result$Zscores)))

  # No extreme outliers (SPA should control for this)
  # Allow for some large values but not absurdly large
  expect_true(all(abs(result$Zscores) < 100),
              info = "SPA should prevent extremely inflated Z-scores")
})

# ==============================================================================
# Test 3: Unbalanced Case-Control Design
# ==============================================================================
test_that("getZ_marg_score_binary_SPA handles unbalanced cases/controls", {
  # Create highly unbalanced data: 10% cases, 90% controls
  set.seed(98765)
  n <- 1000
  p <- 8
  n_case <- 100
  n_ctrl <- 900

  # Generate genotypes with varying MAF
  G <- matrix(0, n, p)
  for (j in 1:p) {
    maf <- runif(1, 0.01, 0.1)
    G[, j] <- rbinom(n, 2, maf)
  }

  # Unbalanced binary outcome
  Y <- c(rep(1, n_case), rep(0, n_ctrl))

  # Covariates
  X <- matrix(rnorm(n * 2), n, 2)

  # Run SPA function
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Check that function handles unbalanced design
  expect_type(result, "list")
  expect_length(result$Zscores, p)
  expect_true(all(is.finite(result$Zscores)))

  # Verify case-control ratio
  case_prop <- mean(Y)
  expect_equal(case_prop, 0.1, tolerance = 0.01)
})

# ==============================================================================
# Test 4: Comparison with Standard Method (Common Variants)
# ==============================================================================
test_that("SPA results similar to standard for common variants", {
  # For common variants with balanced design, SPA and standard should be close
  set.seed(11111)
  n <- 500
  p <- 10

  # Common variants (MAF > 0.1)
  G <- matrix(0, n, p)
  for (j in 1:p) {
    maf <- runif(1, 0.1, 0.4)  # Common MAF
    G[, j] <- rbinom(n, 2, maf)
  }

  # Balanced binary trait
  Y <- sample(c(0, 1), n, replace = TRUE)

  # Covariates
  X <- matrix(rnorm(n * 2), n, 2)

  # Run both methods
  result_spa <- getZ_marg_score_binary_SPA(G, X, Y)
  result_std <- getZ_marg_score(G, cbind(1, X), Y, trait = "binary")

  # For common variants, Z-scores should be reasonably correlated
  # (not identical due to SPA adjustment, but should be close)
  correlation <- cor(result_spa$Zscores, result_std$Zscores)
  expect_true(correlation > 0.95,
              info = paste("Correlation with standard method:", round(correlation, 4)))

  # Note: M_Z matrices will differ because SPA function uses X without intercept
  # while standard function uses X with intercept for projection.
  # Both are mathematically valid; they use different covariate specifications.
})

# ==============================================================================
# Test 5: Edge Case - Single Variant
# ==============================================================================
test_that("getZ_marg_score_binary_SPA works with single variant", {
  set.seed(22222)
  n <- 200
  p <- 1  # Single variant

  # Generate single variant
  G <- matrix(rbinom(n, 2, 0.2), n, 1)

  # Binary trait
  Y <- sample(c(0, 1), n, replace = TRUE)

  # Covariates
  X <- matrix(rnorm(n * 2), n, 2)

  # Run SPA function
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Check output
  expect_length(result$Zscores, 1)
  expect_equal(dim(result$M_Z), c(1, 1))
  expect_equal(result$M_Z[1, 1], 1)  # Correlation with itself is 1
  expect_true(is.finite(result$Zscores[1]))
})

# ==============================================================================
# Test 6: Edge Case - No Covariates (Intercept Only)
# ==============================================================================
test_that("getZ_marg_score_binary_SPA works with no covariates", {
  set.seed(33333)
  n <- 300
  p <- 5

  # Generate genotypes
  G <- matrix(0, n, p)
  for (j in 1:p) {
    G[, j] <- rbinom(n, 2, 0.15)
  }

  # Binary trait
  Y <- sample(c(0, 1), n, replace = TRUE)

  # No covariates except what function adds internally
  X <- matrix(0, n, 0)  # Empty matrix

  # Run SPA function
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Check output
  expect_type(result, "list")
  expect_length(result$Zscores, p)
  expect_true(all(is.finite(result$Zscores)))
})

# ==============================================================================
# Test 7: Edge Case - Monomorphic Variant
# ==============================================================================
test_that("getZ_marg_score_binary_SPA handles monomorphic variants", {
  set.seed(44444)
  n <- 200
  p <- 5

  # Generate genotypes with one monomorphic variant
  G <- matrix(0, n, p)
  for (j in 1:(p-1)) {
    G[, j] <- rbinom(n, 2, 0.2)
  }
  G[, p] <- 0  # Monomorphic (all 0s)

  # Binary trait
  Y <- sample(c(0, 1), n, replace = TRUE)

  # Covariates
  X <- matrix(rnorm(n * 2), n, 2)

  # Run SPA function - should handle gracefully
  # SPAtest should handle monomorphic variants
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Check that we get results (may include NA or 0 for monomorphic)
  expect_type(result, "list")
  expect_length(result$Zscores, p)
})

# ==============================================================================
# Test 7b: Degenerate (zero-variance) variants -> graceful NA, not a hard error
# ==============================================================================
# Regression test: a variant monomorphic in the sample has zero adjusted score
# variance (var1 == 0). SPAtest's saddlepoint routine divides by sqrt(var1) and
# fails with "missing value where TRUE/FALSE needed"; the wrapper must exclude
# such variants and report them as NA (not abort, not emit NaN correlations).
test_that("getZ_marg_score_binary_SPA returns NA for zero-variance variants", {
  set.seed(20260619)
  n <- 200
  maf <- c(0.30, 0.20, 0.25, 0.15, 0.35)
  G <- sapply(maf, function(f) rbinom(n, 2, f))
  G[, 3] <- 0  # monomorphic -> zero adjusted variance (the var1 == 0 case)
  X <- matrix(rnorm(n), n, 1)
  Y <- rbinom(n, 1, plogis(-0.2 + 0.3 * X[, 1]))

  # Must not abort; the wrapper warns about the dropped variant.
  expect_warning(
    result <- getZ_marg_score_binary_SPA(G, X, Y),
    "zero adjusted variance"
  )

  # Degenerate variant -> NA; well-conditioned variants -> finite.
  expect_length(result$Zscores, ncol(G))
  expect_true(is.na(result$Zscores[3]))
  expect_true(all(is.finite(result$Zscores[-3])))
  expect_true(all(is.na(result$M_Z[3, ])))
  expect_true(all(is.finite(result$M_Z[-3, -3])))
})

# ==============================================================================
# Test 8: Input Validation - Dimension Mismatch
# ==============================================================================
test_that("getZ_marg_score_binary_SPA catches dimension mismatches", {
  set.seed(55555)
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- sample(c(0, 1), n, replace = TRUE)

  # Mismatch in number of rows
  G_wrong <- matrix(rnorm((n+10) * p), n+10, p)
  expect_error(
    getZ_marg_score_binary_SPA(G_wrong, X, Y),
    regexp = "same number of observations"
  )

  X_wrong <- matrix(rnorm((n-5) * 2), n-5, 2)
  expect_error(
    getZ_marg_score_binary_SPA(G, X_wrong, Y),
    regexp = "same number of observations"
  )

  Y_wrong <- sample(c(0, 1), n+20, replace = TRUE)
  expect_error(
    getZ_marg_score_binary_SPA(G, X, Y_wrong),
    regexp = "same number of observations"
  )
})

# ==============================================================================
# Test 9: Input Validation - Non-Binary Y
# ==============================================================================
test_that("getZ_marg_score_binary_SPA rejects non-binary Y", {
  set.seed(66666)
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)

  # Non-binary Y
  Y_cont <- rnorm(n)
  expect_error(
    getZ_marg_score_binary_SPA(G, X, Y_cont),
    regexp = "must contain only 0 and 1"
  )

  Y_multi <- sample(c(0, 1, 2), n, replace = TRUE)
  expect_error(
    getZ_marg_score_binary_SPA(G, X, Y_multi),
    regexp = "must contain only 0 and 1"
  )
})

# ==============================================================================
# Test 10: Input Validation - Missing Values
# ==============================================================================
test_that("getZ_marg_score_binary_SPA rejects missing values", {
  set.seed(77777)
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- sample(c(0, 1), n, replace = TRUE)

  # Missing values in G
  G_na <- G
  G_na[1, 1] <- NA
  expect_error(
    getZ_marg_score_binary_SPA(G_na, X, Y),
    regexp = "contains missing values"
  )

  # Missing values in X
  X_na <- X
  X_na[2, 1] <- NA
  expect_error(
    getZ_marg_score_binary_SPA(G, X_na, Y),
    regexp = "contains missing values"
  )

  # Missing values in Y
  Y_na <- Y
  Y_na[3] <- NA
  expect_error(
    getZ_marg_score_binary_SPA(G, X, Y_na),
    regexp = "contains missing values"
  )
})

# ==============================================================================
# Test 11: Input Validation - Matrix Type
# ==============================================================================
test_that("getZ_marg_score_binary_SPA handles matrix type conversions", {
  set.seed(88888)
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- sample(c(0, 1), n, replace = TRUE)

  # X as data.frame (should be converted)
  X_df <- as.data.frame(X)
  expect_silent(result <- getZ_marg_score_binary_SPA(G, X_df, Y))
  expect_type(result, "list")

  # Y as matrix (should be converted)
  Y_mat <- matrix(Y, ncol = 1)
  expect_silent(result <- getZ_marg_score_binary_SPA(G, X, Y_mat))
  expect_type(result, "list")

  # G not as matrix (should error)
  G_df <- as.data.frame(G)
  expect_error(
    getZ_marg_score_binary_SPA(G_df, X, Y),
    regexp = "G must be a matrix"
  )
})

# ==============================================================================
# Test 12: Correlation Matrix Properties
# ==============================================================================
test_that("M_Z has correct correlation matrix properties", {
  # Load test data
  fixture_path <- system.file("validation/fixtures/test_data_binary.rds", package = "GLOWr")
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }
  test_data <- readRDS(fixture_path)

  result <- getZ_marg_score_binary_SPA(test_data$G, strip_intercept(test_data$X), test_data$Y)
  M_Z <- result$M_Z

  # Diagonal should be all 1s
  expect_equal(diag(M_Z), rep(1, nrow(M_Z)), tolerance = 1e-10)

  # Should be symmetric
  expect_equal(M_Z, t(M_Z), tolerance = 1e-10)

  # Off-diagonal elements should be in [-1, 1]
  off_diag <- M_Z[upper.tri(M_Z)]
  expect_true(all(off_diag >= -1 & off_diag <= 1))

  # Should be positive semi-definite (all eigenvalues >= 0)
  eigs <- eigen(M_Z, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10), info = "M_Z should be positive semi-definite")
})

# ==============================================================================
# Test 13: Covariance Matrix Properties
# ==============================================================================
test_that("M_s has correct covariance matrix properties", {
  # Load test data
  fixture_path <- system.file("validation/fixtures/test_data_binary.rds", package = "GLOWr")
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }
  test_data <- readRDS(fixture_path)

  result <- getZ_marg_score_binary_SPA(test_data$G, strip_intercept(test_data$X), test_data$Y)
  M_s <- result$M_s

  # Should be symmetric
  expect_equal(M_s, t(M_s), tolerance = 1e-10)

  # Diagonal elements should be positive (variances)
  expect_true(all(diag(M_s) > 0))

  # Should be positive semi-definite
  eigs <- eigen(M_s, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10), info = "M_s should be positive semi-definite")
})

# ==============================================================================
# Test 14: Consistency Across Runs (Reproducibility)
# ==============================================================================
test_that("getZ_marg_score_binary_SPA is reproducible", {
  set.seed(99999)
  n <- 200
  p <- 8

  G <- matrix(0, n, p)
  for (j in 1:p) {
    G[, j] <- rbinom(n, 2, 0.15)
  }

  Y <- sample(c(0, 1), n, replace = TRUE)
  X <- matrix(rnorm(n * 2), n, 2)

  # Run twice with same data
  result1 <- getZ_marg_score_binary_SPA(G, X, Y)
  result2 <- getZ_marg_score_binary_SPA(G, X, Y)

  # Results should be identical
  expect_equal(result1$Zscores, result2$Zscores)
  expect_equal(result1$M_Z, result2$M_Z)
  expect_equal(result1$M_s, result2$M_s)
  expect_equal(result1$s0, result2$s0)
})

# ==============================================================================
# Test 15: Very Rare Variants (MAC < 10)
# ==============================================================================
test_that("getZ_marg_score_binary_SPA handles very rare variants (MAC < 10)", {
  set.seed(12121)
  n <- 1000
  p <- 5

  # Create variants with very low minor allele counts
  G <- matrix(0, n, p)
  for (j in 1:p) {
    # MAC between 2 and 8
    mac <- sample(2:8, 1)
    carrier_idx <- sample(1:n, mac)
    G[carrier_idx, j] <- sample(c(1, 2), mac, replace = TRUE,
                                 prob = c(0.7, 0.3))
  }

  Y <- sample(c(0, 1), n, replace = TRUE)
  X <- matrix(rnorm(n * 2), n, 2)

  # Verify very rare
  mac <- colSums(G)
  expect_true(all(mac < 10), info = "All variants should have MAC < 10")

  # Run SPA (critical for very rare variants)
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # Should complete without errors
  expect_type(result, "list")
  expect_length(result$Zscores, p)
  expect_true(all(is.finite(result$Zscores)))
})

# ==============================================================================
# Test 16: Sign Preservation in Z-scores
# ==============================================================================
test_that("SPA Z-scores preserve direction of effect", {
  # Create data where we know the direction of association
  set.seed(13131)
  n <- 500
  p <- 3

  # Generate genotypes
  G <- matrix(0, n, p)
  for (j in 1:p) {
    G[, j] <- rbinom(n, 2, 0.2)
  }

  # Create Y with known positive association with G[,1]
  X <- matrix(rnorm(n * 2), n, 2)
  logit_p <- -0.5 + 0.8 * G[, 1] + 0.1 * X[, 1]
  p_case <- 1 / (1 + exp(-logit_p))
  Y <- rbinom(n, 1, p_case)

  # Run SPA
  result <- getZ_marg_score_binary_SPA(G, X, Y)

  # First variant should have positive Z-score (positive association)
  # (This is probabilistic, but with n=500 and strong effect should hold)
  expect_true(result$Zscores[1] > 0,
              info = "Known positive association should have positive Z-score")
})

# ==============================================================================
# Test 17: Output Consistency with Legacy Structure
# ==============================================================================
test_that("Output structure matches legacy implementation", {
  # Load test data
  fixture_path <- system.file("validation/fixtures/test_data_binary.rds", package = "GLOWr")
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }
  test_data <- readRDS(fixture_path)

  result <- getZ_marg_score_binary_SPA(test_data$G, strip_intercept(test_data$X), test_data$Y)

  # Check exact structure expected by downstream GLOW functions
  expect_true(is.list(result))
  expect_true(is.vector(result$Zscores))
  expect_true(is.matrix(result$M_Z))
  expect_true(is.matrix(result$M_s))
  expect_true(is.numeric(result$s0))
  expect_length(result$s0, 1)

  # No extra elements
  expect_setequal(names(result), c("Zscores", "scores", "M_Z", "M_s", "s0"))
})
