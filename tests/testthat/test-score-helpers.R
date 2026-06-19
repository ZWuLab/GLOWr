#' Unit Tests for Internal Score Calculation Helper Functions
#'
#' This file tests the internal helper functions used by getZ_marg_score() and
#' getZ_marg_score_binary_SPA():
#' 1. .validate_score_inputs() - Input validation
#' 2. .compute_projection_matrices() - Matrix computations
#' 3. .handle_zero_variance_snps() - Zero-variance handling
#'
#' File Log (reverse chronological order):
#' - 2025-10-29: Created by r-developer - Tests for refactored helper functions

library(testthat)
library(GLOWr)

# ============================================================================
# Test Suite 1: .validate_score_inputs()
# ============================================================================

test_that(".validate_score_inputs returns correct structure", {
  n <- 100
  p <- 10
  k <- 2

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * k), n, k)
  Y <- rnorm(n)

  # Access internal function using :::
  result <- GLOWr:::.validate_score_inputs(G, X, Y)

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("G", "X", "Y", "n", "p", "k"))

  # Check dimensions
  expect_equal(result$n, n)
  expect_equal(result$p, p)
  expect_equal(result$k, k)

  # Check matrices are returned
  expect_true(is.matrix(result$G))
  expect_true(is.matrix(result$X))
  expect_true(is.vector(result$Y))
})

test_that(".validate_score_inputs converts X to matrix", {
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- rnorm(n)  # Vector, not matrix
  Y <- rnorm(n)

  result <- GLOWr:::.validate_score_inputs(G, X, Y)

  # Should convert X to matrix
  expect_true(is.matrix(result$X))
  expect_equal(ncol(result$X), 1)
  expect_equal(nrow(result$X), n)
})

test_that(".validate_score_inputs converts Y from matrix to vector", {
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- matrix(rnorm(n), n, 1)  # Single-column matrix

  result <- GLOWr:::.validate_score_inputs(G, X, Y)

  # Should convert Y to vector
  expect_true(is.vector(result$Y))
  expect_length(result$Y, n)
})

test_that(".validate_score_inputs catches dimension mismatches", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  # Mismatched G
  G_wrong <- matrix(rnorm(50 * 10), 50, 10)
  expect_error(
    GLOWr:::.validate_score_inputs(G_wrong, X, Y),
    "same number of observations"
  )

  # Mismatched X
  X_wrong <- matrix(rnorm(80 * 2), 80, 2)
  expect_error(
    GLOWr:::.validate_score_inputs(G, X_wrong, Y),
    "same number of observations"
  )

  # Mismatched Y
  Y_wrong <- rnorm(120)
  expect_error(
    GLOWr:::.validate_score_inputs(G, X, Y_wrong),
    "same number of observations"
  )
})

test_that(".validate_score_inputs catches missing values", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  # Missing in G
  G_na <- G
  G_na[1, 1] <- NA
  expect_error(
    GLOWr:::.validate_score_inputs(G_na, X, Y),
    "G contains missing values"
  )

  # Missing in X
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(
    GLOWr:::.validate_score_inputs(G, X_na, Y),
    "X contains missing values"
  )

  # Missing in Y
  Y_na <- Y
  Y_na[1] <- NA
  expect_error(
    GLOWr:::.validate_score_inputs(G, X, Y_na),
    "Y contains missing values"
  )
})

test_that(".validate_score_inputs validates trait type when provided", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  # Invalid trait type
  expect_error(
    GLOWr:::.validate_score_inputs(G, X, Y, trait = "invalid"),
    "trait must be either 'binary' or 'continuous'"
  )

  # Valid trait types should work
  expect_silent(GLOWr:::.validate_score_inputs(G, X, Y, trait = "continuous"))
  expect_silent(GLOWr:::.validate_score_inputs(G, X, sample(c(0, 1), 100, replace = TRUE), trait = "binary"))
})

test_that(".validate_score_inputs validates binary Y values", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)

  # Non-binary Y with binary trait
  Y_cont <- rnorm(100)
  expect_error(
    GLOWr:::.validate_score_inputs(G, X, Y_cont, trait = "binary"),
    "For binary trait, Y must contain only 0 and 1"
  )

  # Binary Y should work
  Y_binary <- sample(c(0, 1), 100, replace = TRUE)
  expect_silent(GLOWr:::.validate_score_inputs(G, X, Y_binary, trait = "binary"))
})

test_that(".validate_score_inputs rejects non-matrix G", {
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  # Data frame instead of matrix
  G_df <- as.data.frame(matrix(rnorm(100 * 10), 100, 10))
  expect_error(
    GLOWr:::.validate_score_inputs(G_df, X, Y),
    "G must be a matrix"
  )
})


# ============================================================================
# Test Suite 2: .compute_projection_matrices()
# ============================================================================

test_that(".compute_projection_matrices returns correct structure", {
  n <- 100
  k <- 3
  p <- 10

  X <- matrix(rnorm(n * k), n, k)
  G <- matrix(rnorm(n * p), n, p)

  result <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("Hhalf", "Gtilde", "Xtilde"))

  # Check dimensions
  expect_equal(dim(result$Hhalf), c(n, k))
  expect_equal(dim(result$Gtilde), c(n, p))
  expect_equal(dim(result$Xtilde), c(n, k))
})

test_that(".compute_projection_matrices handles unweighted case", {
  n <- 100
  k <- 3
  p <- 10

  X <- matrix(rnorm(n * k), n, k)
  G <- matrix(rnorm(n * p), n, p)

  result <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)

  # In unweighted case, Gtilde and Xtilde should be same as inputs
  expect_equal(result$Gtilde, G)
  expect_equal(result$Xtilde, X)

  # Hhalf should satisfy: Hhalf %*% t(Hhalf) ≈ X %*% solve(t(X) %*% X) %*% t(X)
  # This is the projection matrix onto column space of X
  proj_expected <- X %*% solve(t(X) %*% X) %*% t(X)
  proj_computed <- result$Hhalf %*% t(result$Hhalf)
  expect_equal(proj_computed, proj_expected, tolerance = 1e-10)
})

test_that(".compute_projection_matrices handles weighted case", {
  set.seed(123)
  n <- 100
  k <- 3
  p <- 10

  X <- matrix(rnorm(n * k), n, k)
  G <- matrix(rnorm(n * p), n, p)
  weights <- sqrt(runif(n, 0.1, 0.9))  # Positive weights

  result <- GLOWr:::.compute_projection_matrices(X, G, weights = weights)

  # In weighted case, Gtilde and Xtilde should be weighted versions
  expect_equal(result$Gtilde, G * weights)
  expect_equal(result$Xtilde, X * weights)

  # Hhalf should have correct dimensions
  expect_equal(dim(result$Hhalf), c(n, k))
})

test_that(".compute_projection_matrices Hhalf is numerically stable", {
  set.seed(456)
  n <- 50
  k <- 2
  p <- 5

  X <- cbind(1, rnorm(n))
  G <- matrix(rnorm(n * p), n, p)

  result <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)

  # Check that Hhalf has no NaN or Inf values
  expect_false(any(is.nan(result$Hhalf)))
  expect_false(any(is.infinite(result$Hhalf)))
})

test_that(".compute_projection_matrices catches singular X", {
  n <- 100
  k <- 3
  p <- 10

  # Create singular X (duplicate column)
  X <- cbind(rnorm(n), rnorm(n), rnorm(n))
  X[, 3] <- X[, 2]  # Make singular
  G <- matrix(rnorm(n * p), n, p)

  # Should error with appropriate message
  expect_error(
    GLOWr:::.compute_projection_matrices(X, G, weights = NULL),
    "singular|collinear"
  )
})

test_that(".compute_projection_matrices catches singular weighted X", {
  n <- 100
  k <- 2
  p <- 10

  X <- cbind(rnorm(n), rnorm(n))
  X[, 2] <- X[, 1]  # Make singular
  G <- matrix(rnorm(n * p), n, p)
  weights <- sqrt(runif(n, 0.1, 0.9))

  # Should error with appropriate message for weighted case
  expect_error(
    GLOWr:::.compute_projection_matrices(X, G, weights = weights),
    "singular|separation|collinear"
  )
})

test_that(".compute_projection_matrices handles intercept-only model", {
  n <- 100
  p <- 5

  X <- matrix(1, n, 1)  # Intercept only
  G <- matrix(rnorm(n * p), n, p)

  result <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)

  # Should work without errors
  expect_type(result, "list")
  expect_equal(dim(result$Hhalf), c(n, 1))
})

test_that(".compute_projection_matrices projection properties", {
  set.seed(789)
  n <- 80
  k <- 3
  p <- 8

  X <- cbind(1, matrix(rnorm(n * (k-1)), n, k-1))
  G <- matrix(rnorm(n * p), n, p)

  result <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)

  # Projection matrix P = Hhalf %*% t(Hhalf) should be idempotent: P %*% P = P
  P <- result$Hhalf %*% t(result$Hhalf)
  P2 <- P %*% P
  expect_equal(P, P2, tolerance = 1e-10)

  # Projection should be symmetric
  expect_equal(P, t(P), tolerance = 1e-10)
})


# ============================================================================
# Test Suite 3: .handle_zero_variance_snps()
# ============================================================================

test_that(".handle_zero_variance_snps detects zero-variance SNPs", {
  # Create diagonal with some zero elements
  diag_GHG <- c(1.0, 2.5, 0.0, 3.2, 1e-20, 0.5)

  result <- GLOWr:::.handle_zero_variance_snps(diag_GHG)

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("zero_var_idx", "diag_GHG_safe", "has_zero_var"))

  # Should detect indices 3 and 5 (0.0 and 1e-20)
  expect_true(result$has_zero_var)
  expect_true(3 %in% result$zero_var_idx)
  expect_true(5 %in% result$zero_var_idx)
  expect_length(result$zero_var_idx, 2)
})

test_that(".handle_zero_variance_snps replaces zeros with 1", {
  diag_GHG <- c(1.0, 0.0, 2.5)

  result <- GLOWr:::.handle_zero_variance_snps(diag_GHG)

  # Safe version should have 1 instead of 0
  expect_equal(result$diag_GHG_safe[1], 1.0)
  expect_equal(result$diag_GHG_safe[2], 1.0)  # Replaced
  expect_equal(result$diag_GHG_safe[3], 2.5)
})

test_that(".handle_zero_variance_snps handles no zero-variance case", {
  diag_GHG <- c(1.0, 2.5, 3.2, 0.5, 1.8)

  result <- GLOWr:::.handle_zero_variance_snps(diag_GHG)

  # Should find no zero-variance SNPs
  expect_false(result$has_zero_var)
  expect_length(result$zero_var_idx, 0)

  # Safe version should be identical to input
  expect_equal(result$diag_GHG_safe, diag_GHG)
})

test_that(".handle_zero_variance_snps handles all zero-variance case", {
  diag_GHG <- rep(0.0, 10)

  result <- GLOWr:::.handle_zero_variance_snps(diag_GHG)

  # Should detect all as zero-variance
  expect_true(result$has_zero_var)
  expect_length(result$zero_var_idx, 10)

  # Safe version should be all 1s
  expect_equal(result$diag_GHG_safe, rep(1.0, 10))
})

test_that(".handle_zero_variance_snps uses machine epsilon threshold", {
  # Test values around machine epsilon
  diag_GHG <- c(
    .Machine$double.eps * 2,    # Above threshold - should NOT be flagged
    .Machine$double.eps * 0.5,  # Below threshold - should be flagged
    .Machine$double.eps * 0.1   # Well below threshold - should be flagged
  )

  result <- GLOWr:::.handle_zero_variance_snps(diag_GHG)

  # Should only detect 2nd and 3rd as zero-variance
  expect_true(result$has_zero_var)
  expect_equal(result$zero_var_idx, c(2, 3))
})

test_that(".handle_zero_variance_snps preserves non-zero values", {
  diag_GHG <- c(1.5, 0.0, 2.8, 1e-20, 0.3, 4.2)

  result <- GLOWr:::.handle_zero_variance_snps(diag_GHG)

  # Non-zero values should be unchanged
  expect_equal(result$diag_GHG_safe[1], 1.5)
  expect_equal(result$diag_GHG_safe[3], 2.8)
  expect_equal(result$diag_GHG_safe[5], 0.3)
  expect_equal(result$diag_GHG_safe[6], 4.2)

  # Zero values should be replaced with 1
  expect_equal(result$diag_GHG_safe[2], 1.0)
  expect_equal(result$diag_GHG_safe[4], 1.0)
})


# ============================================================================
# Test Suite 4: Integration Tests - Helpers Work Together
# ============================================================================

test_that("Helpers integrate correctly in binary trait workflow", {
  set.seed(999)
  n <- 100
  p <- 10
  k <- 2

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, matrix(rnorm(n * k), n, k))
  Y <- sample(c(0, 1), n, replace = TRUE)

  # Step 1: Validate inputs
  validated <- GLOWr:::.validate_score_inputs(G, X, Y, trait = "binary")

  # Step 2: Fit model and compute weights
  mod0 <- glm(validated$Y ~ validated$X, family = "binomial")
  Y0 <- mod0$fitted.values
  w <- sqrt(Y0 * (1 - Y0))

  # Step 3: Compute projection matrices
  proj_matrices <- GLOWr:::.compute_projection_matrices(validated$X, validated$G, weights = w)

  # Step 4: Compute GHG
  GHhalf <- t(proj_matrices$Gtilde) %*% proj_matrices$Hhalf
  GHG <- t(proj_matrices$Gtilde) %*% proj_matrices$Gtilde - GHhalf %*% t(GHhalf)

  # Step 5: Handle zero-variance
  zero_var <- GLOWr:::.handle_zero_variance_snps(diag(GHG))

  # All steps should complete without errors
  expect_type(validated, "list")
  expect_type(proj_matrices, "list")
  expect_true(is.matrix(GHG))
  expect_type(zero_var, "list")
})

test_that("Helpers integrate correctly in continuous trait workflow", {
  set.seed(888)
  n <- 100
  p <- 10
  k <- 2

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, matrix(rnorm(n * k), n, k))
  Y <- rnorm(n)

  # Step 1: Validate inputs
  validated <- GLOWr:::.validate_score_inputs(G, X, Y, trait = "continuous")

  # Step 2: Compute projection matrices (no weights)
  proj_matrices <- GLOWr:::.compute_projection_matrices(validated$X, validated$G, weights = NULL)

  # Step 3: Compute GHG
  GHhalf <- t(validated$G) %*% proj_matrices$Hhalf
  GHG <- t(validated$G) %*% validated$G - GHhalf %*% t(GHhalf)

  # Step 4: Handle zero-variance
  zero_var <- GLOWr:::.handle_zero_variance_snps(diag(GHG))

  # All steps should complete without errors
  expect_type(validated, "list")
  expect_type(proj_matrices, "list")
  expect_true(is.matrix(GHG))
  expect_type(zero_var, "list")
})


# ============================================================================
# Test Suite 5: Edge Cases
# ============================================================================

test_that("Helpers handle single SNP correctly", {
  n <- 100
  p <- 1  # Single SNP

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)

  validated <- GLOWr:::.validate_score_inputs(G, X, Y)
  expect_equal(validated$p, 1)

  proj_matrices <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)
  expect_equal(ncol(proj_matrices$Gtilde), 1)

  diag_GHG <- c(2.5)
  zero_var <- GLOWr:::.handle_zero_variance_snps(diag_GHG)
  expect_false(zero_var$has_zero_var)
})

test_that("Helpers handle empty covariate matrix correctly", {
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(0, n, 0)  # Empty matrix
  Y <- rnorm(n)

  validated <- GLOWr:::.validate_score_inputs(G, X, Y)
  expect_equal(validated$k, 0)
  expect_equal(ncol(validated$X), 0)
})

test_that("Helpers handle large number of SNPs", {
  n <- 50
  p <- 100  # More SNPs than samples

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)

  validated <- GLOWr:::.validate_score_inputs(G, X, Y)
  expect_equal(validated$p, 100)

  proj_matrices <- GLOWr:::.compute_projection_matrices(X, G, weights = NULL)
  expect_equal(ncol(proj_matrices$Gtilde), 100)
})
