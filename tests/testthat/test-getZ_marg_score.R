#' Unit Tests for getZ_marg_score()
#'
#' This file tests the getZ_marg_score() function for:
#' 1. Correctness with continuous traits
#' 2. Correctness with binary traits
#' 3. Edge cases (single SNP, collinearity, etc.)
#' 4. Input validation
#' 5. Validation against legacy implementation
#'
#' File Log (reverse chronological order):
#' - 2026-05-26: Cleaned up by Claude Code (Opus 4.7), prompted by ZWu. Fixed 9
#'   failures with two root causes, both test-file defects: (1) three legacy-
#'   comparison tests called source(legacy_code_path) into the GLOBAL env, which
#'   shadowed the package getZ_marg_score (legacy lacks the null_model arg and
#'   zero-variance handling) for every later test — now sourced into a local
#'   env and called as legacy_env$getZ_marg_score(); (2) the test file had each
#'   of [singular-X, constant-genotype, C++ GHG] triplicated (a copy-paste of
#'   the block after each legacy test) — removed the two redundant copies of
#'   each (36 -> 30 blocks). Strengthened the constant-genotype test to assert
#'   M_Z stays finite (paired with the getZ_marg_score.R cov2cor robustness fix).
#' - 2025-10-18: Created by r-developer - Initial test suite for getZ_marg_score()

library(testthat)
library(GLOWr)

# ============================================================================
# Test 1: Basic Functionality - Continuous Trait
# ============================================================================

test_that("getZ_marg_score works for continuous trait", {
  set.seed(12345)
  n <- 100
  p <- 10
  k <- 2

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, matrix(rnorm(n * k), n, k))
  Y <- rnorm(n)

  result <- getZ_marg_score(G, X, Y, trait = "continuous")

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("Zscores", "scores", "M_Z", "M_s", "s0"))

  # Check dimensions
  expect_length(result$Zscores, p)
  expect_length(result$scores, p)
  expect_equal(dim(result$M_Z), c(p, p))
  expect_equal(dim(result$M_s), c(p, p))
  expect_length(result$s0, 1)

  # Check types
  expect_type(result$Zscores, "double")
  expect_type(result$scores, "double")
  expect_true(is.matrix(result$M_Z))
  expect_true(is.matrix(result$M_s))
  expect_type(result$s0, "double")

  # Check that M_Z is a correlation matrix
  expect_true(all(diag(result$M_Z) - 1 < 1e-10))  # Diagonal should be 1
  expect_true(all(result$M_Z >= -1 - 1e-10 & result$M_Z <= 1 + 1e-10))  # Between -1 and 1
  expect_true(isSymmetric(result$M_Z))  # Should be symmetric

  # Check that M_s is symmetric
  expect_true(isSymmetric(result$M_s))

  # Check that s0 is positive for continuous trait
  expect_true(result$s0 > 0)

  # Check no NaN or Inf values
  expect_false(any(is.nan(result$Zscores)))
  expect_false(any(is.infinite(result$Zscores)))
})


# ============================================================================
# Test 2: Basic Functionality - Binary Trait
# ============================================================================

test_that("getZ_marg_score works for binary trait", {
  set.seed(23456)
  n <- 200
  p <- 15
  k <- 2

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, matrix(rnorm(n * k), n, k))
  Y <- sample(c(0, 1), n, replace = TRUE)

  result <- getZ_marg_score(G, X, Y, trait = "binary")

  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("Zscores", "scores", "M_Z", "M_s", "s0"))

  # Check dimensions
  expect_length(result$Zscores, p)
  expect_length(result$scores, p)
  expect_equal(dim(result$M_Z), c(p, p))
  expect_equal(dim(result$M_s), c(p, p))

  # Check that s0 is 1 for binary trait
  expect_equal(result$s0, 1)

  # Check that M_Z is a correlation matrix
  expect_true(all(abs(diag(result$M_Z) - 1) < 1e-10))
  expect_true(isSymmetric(result$M_Z))

  # Check no NaN or Inf values
  expect_false(any(is.nan(result$Zscores)))
  expect_false(any(is.infinite(result$Zscores)))
})


# ============================================================================
# Test 3: Continuous Trait with use_lm_t Option
# ============================================================================

test_that("getZ_marg_score use_lm_t option works for continuous trait", {
  set.seed(34567)
  n <- 100
  p <- 10
  k <- 2

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, matrix(rnorm(n * k), n, k))
  Y <- rnorm(n)

  result_default <- getZ_marg_score(G, X, Y, trait = "continuous", use_lm_t = FALSE)
  result_lm_t <- getZ_marg_score(G, X, Y, trait = "continuous", use_lm_t = TRUE)

  # Both should return valid results
  expect_type(result_default, "list")
  expect_type(result_lm_t, "list")

  # Z-scores should be similar but not identical
  # (They use different variance estimators)
  expect_false(identical(result_default$Zscores, result_lm_t$Zscores))

  # Correlation should be very high between the two methods
  cor_z <- cor(result_default$Zscores, result_lm_t$Zscores)
  expect_true(cor_z > 0.95)

  # M_Z and M_s should be identical (they don't depend on use_lm_t)
  expect_equal(result_default$M_Z, result_lm_t$M_Z)
  expect_equal(result_default$M_s, result_lm_t$M_s)
  expect_equal(result_default$s0, result_lm_t$s0)
})


# ============================================================================
# Test 4: Edge Case - Single SNP
# ============================================================================

test_that("getZ_marg_score works with single SNP", {
  set.seed(45678)
  n <- 100
  p <- 1  # Single SNP

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)

  result <- getZ_marg_score(G, X, Y, trait = "continuous")

  # Check dimensions
  expect_length(result$Zscores, 1)
  expect_equal(dim(result$M_Z), c(1, 1))
  expect_equal(dim(result$M_s), c(1, 1))

  # M_Z should be matrix with single element = 1
  expect_equal(result$M_Z[1, 1], 1)

  # No NaN or Inf
  expect_false(is.nan(result$Zscores[1]))
  expect_false(is.infinite(result$Zscores[1]))
})


# ============================================================================
# Test 5: Edge Case - Intercept Only (No Covariates)
# ============================================================================

test_that("getZ_marg_score works with intercept only", {
  set.seed(56789)
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(1, n, 1)  # Intercept only
  Y <- rnorm(n)

  result_cont <- getZ_marg_score(G, X, Y, trait = "continuous")
  expect_type(result_cont, "list")
  expect_length(result_cont$Zscores, p)

  Y_binary <- sample(c(0, 1), n, replace = TRUE)
  result_binary <- getZ_marg_score(G, X, Y_binary, trait = "binary")
  expect_type(result_binary, "list")
  expect_length(result_binary$Zscores, p)
})


# ============================================================================
# Test 6: Input Validation - Wrong Dimensions
# ============================================================================

test_that("getZ_marg_score validates input dimensions", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  # Mismatched G and X dimensions
  X_wrong <- matrix(rnorm(50 * 2), 50, 2)
  expect_error(
    getZ_marg_score(G, X_wrong, Y, trait = "continuous"),
    "must have the same number of observations"
  )

  # Mismatched Y dimension
  Y_wrong <- rnorm(50)
  expect_error(
    getZ_marg_score(G, X, Y_wrong, trait = "continuous"),
    "must have the same number of observations"
  )
})


# ============================================================================
# Test 7: Input Validation - Missing Values
# ============================================================================

test_that("getZ_marg_score rejects missing values", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  # Missing in G
  G_na <- G
  G_na[1, 1] <- NA
  expect_error(
    getZ_marg_score(G_na, X, Y, trait = "continuous"),
    "G contains missing values"
  )

  # Missing in X
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(
    getZ_marg_score(G, X_na, Y, trait = "continuous"),
    "X contains missing values"
  )

  # Missing in Y
  Y_na <- Y
  Y_na[1] <- NA
  expect_error(
    getZ_marg_score(G, X, Y_na, trait = "continuous"),
    "Y contains missing values"
  )
})


# ============================================================================
# Test 8: Input Validation - Invalid Trait Type
# ============================================================================

test_that("getZ_marg_score validates trait type", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)

  expect_error(
    getZ_marg_score(G, X, Y, trait = "invalid"),
    "trait must be either 'binary' or 'continuous'"
  )
})


# ============================================================================
# Test 9: Input Validation - Binary Trait with Non-0/1 Values
# ============================================================================

test_that("getZ_marg_score validates binary trait values", {
  G <- matrix(rnorm(100 * 10), 100, 10)
  X <- matrix(rnorm(100 * 2), 100, 2)
  Y <- rnorm(100)  # Continuous values

  expect_error(
    getZ_marg_score(G, X, Y, trait = "binary"),
    "For binary trait, Y must contain only 0 and 1"
  )
})


# ============================================================================
# Test 10: Input Validation - Matrix Types
# ============================================================================

test_that("getZ_marg_score handles different input types", {
  n <- 100
  p <- 5

  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rnorm(n)

  # X as vector (should be converted to matrix)
  X_vec <- rnorm(n)
  result <- getZ_marg_score(G, X_vec, Y, trait = "continuous")
  expect_type(result, "list")

  # Y as single-column matrix (should be converted to vector)
  Y_mat <- matrix(Y, n, 1)
  result <- getZ_marg_score(G, X, Y_mat, trait = "continuous")
  expect_type(result, "list")

  # G must be matrix
  G_df <- as.data.frame(G)
  expect_error(
    getZ_marg_score(G_df, X, Y, trait = "continuous"),
    "G must be a matrix"
  )

  # Y with multiple columns should fail
  Y_multi <- matrix(rnorm(n * 2), n, 2)
  expect_error(
    getZ_marg_score(G, X, Y_multi, trait = "continuous"),
    "Y must be a single column"
  )
})


# ============================================================================
# Test 11: Mathematical Properties - Correlation Matrix
# ============================================================================

test_that("M_Z is a valid correlation matrix", {
  set.seed(67890)
  n <- 150
  p <- 8

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)

  result <- getZ_marg_score(G, X, Y, trait = "continuous")

  M_Z <- result$M_Z

  # Test 1: Diagonal elements are 1
  expect_true(all(abs(diag(M_Z) - 1) < 1e-10))

  # Test 2: Off-diagonal elements are between -1 and 1
  M_Z_off <- M_Z
  diag(M_Z_off) <- 0
  expect_true(all(M_Z_off >= -1 - 1e-10 & M_Z_off <= 1 + 1e-10))

  # Test 3: Matrix is symmetric
  expect_true(isSymmetric(M_Z))

  # Test 4: Matrix is positive semi-definite
  eigenvalues <- eigen(M_Z, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))
})


# ============================================================================
# Test 12: Reproducibility
# ============================================================================

test_that("getZ_marg_score gives reproducible results", {
  set.seed(11111)
  n <- 100
  p <- 10

  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)

  result1 <- getZ_marg_score(G, X, Y, trait = "continuous")
  result2 <- getZ_marg_score(G, X, Y, trait = "continuous")

  # All components should be identical
  expect_identical(result1$Zscores, result2$Zscores)
  expect_identical(result1$scores, result2$scores)
  expect_identical(result1$M_Z, result2$M_Z)
  expect_identical(result1$M_s, result2$M_s)
  expect_identical(result1$s0, result2$s0)
})


# ============================================================================
# Test 13: Validation Against Legacy Implementation - Continuous Trait
# ============================================================================

test_that("getZ_marg_score matches legacy for continuous trait", {
  # Load test fixture
  fixture_path <- "../../inst/validation/fixtures/test_data_continuous.rds"
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }

  test_data <- readRDS(fixture_path)

  # Run new implementation
  result_new <- getZ_marg_score(
    G = test_data$G,
    X = test_data$X,
    Y = test_data$Y,
    trait = "continuous",
    use_lm_t = FALSE
  )

  # Run legacy implementation (source from legacy materials)
  legacy_code_path <- "../../../../legacy-materials/code/GLOW_R_pacakge/GLOW/R/getZ_marg_score.R"
  if (!file.exists(legacy_code_path)) {
    skip("Legacy code not found")
  }

  legacy_env <- new.env()
  source(legacy_code_path, local = legacy_env)
  result_legacy <- legacy_env$getZ_marg_score(
    G = test_data$G,
    X = test_data$X,
    Y = test_data$Y,
    trait = "continuous",
    use_lm_t = FALSE
  )

  # Validate: tolerance < 1e-10
  tol <- 1e-10

  expect_equal(result_new$Zscores, result_legacy$Zscores, tolerance = tol)
  expect_equal(result_new$scores, result_legacy$scores, tolerance = tol)
  expect_equal(result_new$M_Z, result_legacy$M_Z, tolerance = tol)
  expect_equal(result_new$M_s, result_legacy$M_s, tolerance = tol)
  expect_equal(result_new$s0, result_legacy$s0, tolerance = tol)
})


# ============================================================================
# Test 16: Robustness - Singular Covariate Matrix
# ============================================================================

test_that("getZ_marg_score handles singular X matrix", {
  set.seed(123)
  n <- 100
  p <- 5
  G <- matrix(rnorm(n * p), n, p)
  X <- cbind(1, rnorm(n), rnorm(n))
  X[, 3] <- X[, 2]  # Make X singular (duplicate column)
  Y <- rnorm(n)

  expect_error(
    getZ_marg_score(G, X, Y, trait = "continuous"),
    "singular|collinear"
  )
})


# ============================================================================
# Test 17: Robustness - Zero Variance SNP
# ============================================================================

test_that("getZ_marg_score handles constant genotype", {
  set.seed(456)
  n <- 100
  p <- 5
  G <- matrix(rnorm(n * p), n, p)
  G[, 1] <- 0  # Make first SNP constant
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)

  # Exactly the "zero variance" warning is expected — no secondary cov2cor
  # "non-positive diagonal" warning (the zero-variance column is sanitized
  # before cov2cor so M_Z stays finite).
  expect_warning(
    result <- getZ_marg_score(G, X, Y, trait = "continuous"),
    "zero variance"
  )
  expect_equal(result$Zscores[1], 0)  # Constant SNP should have Z=0
  expect_false(any(is.infinite(result$Zscores)))
  expect_false(any(is.nan(result$M_Z)))   # correlation matrix stays finite
})


# ============================================================================
# Test 18: C++ Implementation Matches R Implementation
# ============================================================================

test_that("C++ GHG matches R implementation", {
  set.seed(789)
  n <- 100
  p <- 10
  k <- 3
  G <- matrix(rnorm(n * p), n, p)
  Hhalf <- matrix(rnorm(n * k), n, k)

  # C++ version
  cpp_result <- compute_GHG_cpp(G, Hhalf)

  # R version
  GHhalf_r <- t(G) %*% Hhalf
  GHG_r <- t(G) %*% G - GHhalf_r %*% t(GHhalf_r)

  # Compare
  expect_equal(cpp_result$GHG, GHG_r, tolerance = 1e-12)
  expect_equal(cpp_result$GHhalf, GHhalf_r, tolerance = 1e-12)
})


# ============================================================================
# Test 14: Validation Against Legacy Implementation - Binary Trait
# ============================================================================

test_that("getZ_marg_score matches legacy for binary trait", {
  # Load test fixture
  fixture_path <- "../../inst/validation/fixtures/test_data_binary.rds"
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }

  test_data <- readRDS(fixture_path)

  # Run new implementation
  result_new <- getZ_marg_score(
    G = test_data$G,
    X = test_data$X,
    Y = test_data$Y,
    trait = "binary"
  )

  # Run legacy implementation
  legacy_code_path <- "../../../../legacy-materials/code/GLOW_R_pacakge/GLOW/R/getZ_marg_score.R"
  if (!file.exists(legacy_code_path)) {
    skip("Legacy code not found")
  }

  legacy_env <- new.env()
  source(legacy_code_path, local = legacy_env)
  result_legacy <- legacy_env$getZ_marg_score(
    G = test_data$G,
    X = test_data$X,
    Y = test_data$Y,
    trait = "binary"
  )

  # Validate: tolerance < 1e-10
  tol <- 1e-10

  expect_equal(result_new$Zscores, result_legacy$Zscores, tolerance = tol)
  expect_equal(result_new$scores, result_legacy$scores, tolerance = tol)
  expect_equal(result_new$M_Z, result_legacy$M_Z, tolerance = tol)
  expect_equal(result_new$M_s, result_legacy$M_s, tolerance = tol)
  expect_equal(result_new$s0, result_legacy$s0, tolerance = tol)
})


# ============================================================================
# Test 15: Validation Against Legacy Implementation - use_lm_t = TRUE
# ============================================================================

test_that("getZ_marg_score matches legacy for continuous trait with use_lm_t=TRUE", {
  # Load test fixture
  fixture_path <- "../../inst/validation/fixtures/test_data_simple.rds"
  if (!file.exists(fixture_path)) {
    skip("Test fixture not found")
  }

  test_data <- readRDS(fixture_path)

  # Run new implementation
  result_new <- getZ_marg_score(
    G = test_data$G,
    X = test_data$X,
    Y = test_data$Y,
    trait = "continuous",
    use_lm_t = TRUE
  )

  # Run legacy implementation
  legacy_code_path <- "../../../../legacy-materials/code/GLOW_R_pacakge/GLOW/R/getZ_marg_score.R"
  if (!file.exists(legacy_code_path)) {
    skip("Legacy code not found")
  }

  legacy_env <- new.env()
  source(legacy_code_path, local = legacy_env)
  result_legacy <- legacy_env$getZ_marg_score(
    G = test_data$G,
    X = test_data$X,
    Y = test_data$Y,
    trait = "continuous",
    use_lm_t = TRUE
  )

  # Validate: tolerance < 1e-10
  tol <- 1e-10

  expect_equal(result_new$Zscores, result_legacy$Zscores, tolerance = tol)
  expect_equal(result_new$scores, result_legacy$scores, tolerance = tol)
  expect_equal(result_new$M_Z, result_legacy$M_Z, tolerance = tol)
  expect_equal(result_new$M_s, result_legacy$M_s, tolerance = tol)
  expect_equal(result_new$s0, result_legacy$s0, tolerance = tol)
})


# ============================================================================
# Test 19: fit_null_model() Basic Functionality
# ============================================================================

test_that("fit_null_model() works for binary traits", {
  set.seed(100)
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  # Fit null model
  null_model <- fit_null_model(X, Y, trait = "binary")

  # Check class
  expect_s3_class(null_model, "glow_null_model")

  # Check components
  expect_equal(null_model$trait, "binary")
  expect_equal(null_model$n, n)
  expect_equal(nrow(null_model$X), n)
  expect_equal(ncol(null_model$X), 3)  # Should add intercept
  expect_equal(length(null_model$fitted_probs), n)
  expect_equal(length(null_model$Y), n)

  # Check fitted probabilities are valid
  expect_true(all(null_model$fitted_probs >= 0))
  expect_true(all(null_model$fitted_probs <= 1))
})

test_that("fit_null_model() works for continuous traits", {
  set.seed(101)
  n <- 100
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rnorm(n)

  # Fit null model
  null_model <- fit_null_model(X, Y, trait = "continuous")

  # Check class
  expect_s3_class(null_model, "glow_null_model")

  # Check components
  expect_equal(null_model$trait, "continuous")
  expect_equal(null_model$n, n)
  expect_equal(nrow(null_model$X), n)
  expect_equal(ncol(null_model$X), 3)  # Should add intercept
  expect_equal(length(null_model$residuals), n)
  expect_true(is.numeric(null_model$s0))
  expect_true(null_model$s0 > 0)
})


# ============================================================================
# Test 20: Equivalence Between Standard and Optimized Modes - Continuous Trait
# ============================================================================

test_that("getZ_marg_score() produces identical results with and without null_model (continuous)", {
  set.seed(200)
  n <- 150
  p <- 5
  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- rnorm(n)

  # Standard mode: provide X and Y
  Z_standard <- getZ_marg_score(G, X, Y, trait = "continuous")

  # Optimized mode: pre-compute null model
  null_model <- fit_null_model(X, Y, trait = "continuous")
  Z_optimized <- getZ_marg_score(G, null_model = null_model)

  # Should be identical
  expect_equal(Z_optimized, Z_standard, tolerance = 1e-12)
})

test_that("getZ_marg_score() with use_lm_t=TRUE produces identical results (continuous)", {
  set.seed(201)
  n <- 150
  p <- 5
  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- rnorm(n)

  # Standard mode with use_lm_t
  Z_standard <- getZ_marg_score(G, X, Y, trait = "continuous", use_lm_t = TRUE)

  # Optimized mode with use_lm_t
  null_model <- fit_null_model(X, Y, trait = "continuous")
  Z_optimized <- getZ_marg_score(G, null_model = null_model, use_lm_t = TRUE)

  # Should be identical
  expect_equal(Z_optimized, Z_standard, tolerance = 1e-12)
})


# ============================================================================
# Test 21: Equivalence Between Standard and Optimized Modes - Binary Trait
# ============================================================================

test_that("getZ_marg_score() produces identical results with and without null_model (binary)", {
  set.seed(300)
  n <- 150
  p <- 5
  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- rbinom(n, 1, 0.4)

  # Standard mode: provide X and Y
  Z_standard <- getZ_marg_score(G, X, Y, trait = "binary")

  # Optimized mode: pre-compute null model
  null_model <- fit_null_model(X, Y, trait = "binary")
  Z_optimized <- getZ_marg_score(G, null_model = null_model)

  # Should be identical
  expect_equal(Z_optimized, Z_standard, tolerance = 1e-12)
})

test_that("getZ_marg_score_binary_SPA() produces identical results with and without null_model", {
  set.seed(301)
  n <- 200
  p <- 8
  G <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  # Standard mode: provide X and Y
  Z_standard <- getZ_marg_score_binary_SPA(G, X, Y)

  # Optimized mode: pre-compute null model
  null_model <- fit_null_model(X, Y, trait = "binary")
  Z_optimized <- getZ_marg_score_binary_SPA(G, null_model = null_model)

  # Should be identical
  expect_equal(Z_optimized, Z_standard, tolerance = 1e-10)
})


# ============================================================================
# Test 22: Dimension Validation with null_model
# ============================================================================

test_that("getZ_marg_score() validates dimensions when using null_model", {
  set.seed(400)
  n <- 100
  G <- matrix(rnorm(n * 5), n, 5)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rnorm(n)

  # Fit null model with n=100
  null_model <- fit_null_model(X, Y, trait = "continuous")

  # Try to use with mismatched G (n=50)
  G_wrong <- matrix(rnorm(50 * 5), 50, 5)
  expect_error(
    getZ_marg_score(G_wrong, null_model = null_model),
    "Dimensions must match"
  )
})

test_that("getZ_marg_score_binary_SPA() validates dimensions when using null_model", {
  set.seed(401)
  n <- 100
  G <- matrix(rnorm(n * 5), n, 5)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.3)

  # Fit null model with n=100
  null_model <- fit_null_model(X, Y, trait = "binary")

  # Try to use with mismatched G (n=80)
  G_wrong <- matrix(rnorm(80 * 5), 80, 5)
  expect_error(
    getZ_marg_score_binary_SPA(G_wrong, null_model = null_model),
    "Dimensions must match"
  )
})


# ============================================================================
# Test 23: Error Handling for Invalid null_model
# ============================================================================

test_that("getZ_marg_score() rejects invalid null_model objects", {
  set.seed(500)
  n <- 100
  G <- matrix(rnorm(n * 5), n, 5)

  # Not a glow_null_model object
  fake_null <- list(trait = "binary", n = n)
  expect_error(
    getZ_marg_score(G, null_model = fake_null),
    "null_model must be created by fit_null_model"
  )
})

test_that("getZ_marg_score_binary_SPA() rejects continuous null_model", {
  set.seed(501)
  n <- 100
  G <- matrix(rnorm(n * 5), n, 5)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rnorm(n)

  # Fit continuous trait null model
  null_model <- fit_null_model(X, Y, trait = "continuous")

  # SPA function requires binary trait
  expect_error(
    getZ_marg_score_binary_SPA(G, null_model = null_model),
    "binary trait"
  )
})

test_that("getZ_marg_score() requires either (X,Y) or null_model", {
  set.seed(502)
  n <- 100
  G <- matrix(rnorm(n * 5), n, 5)

  # Neither mode provided
  expect_error(
    getZ_marg_score(G),
    "Either provide"
  )
})


# ============================================================================
# Test 24: Multiple SNP Sets with Same null_model (Realistic Scenario)
# ============================================================================

test_that("null_model can be reused across multiple SNP sets", {
  set.seed(600)
  n <- 200
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- rbinom(n, 1, 0.35)

  # Pre-compute null model once
  null_model <- fit_null_model(X, Y, trait = "binary")

  # Test 10 different SNP sets
  num_sets <- 10
  Z_standard_list <- vector("list", num_sets)
  Z_optimized_list <- vector("list", num_sets)

  for (i in 1:num_sets) {
    set.seed(600 + i)
    p <- sample(3:8, 1)
    G <- matrix(rnorm(n * p), n, p)

    # Standard mode
    Z_standard_list[[i]] <- getZ_marg_score(G, X, Y, trait = "binary")

    # Optimized mode (reuse null_model)
    Z_optimized_list[[i]] <- getZ_marg_score(G, null_model = null_model)
  }

  # All should match
  for (i in 1:num_sets) {
    expect_equal(Z_optimized_list[[i]], Z_standard_list[[i]], tolerance = 1e-12)
  }
})
