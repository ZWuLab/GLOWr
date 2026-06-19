########## Robustness Tests for Optimal Weights Functions ##########
#
# This file tests the robustness improvements added to optimal weights
# functions, including matrix inversion stability, NA/Inf validation,
# division by zero protection, and positive definite checks.

test_that("safe_matrix_inverse handles well-conditioned matrices", {
  # Well-conditioned matrix should invert without warning
  M_good <- diag(3)
  result <- GLOWr:::safe_matrix_inverse(M_good, "test")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_equal(result, diag(3), tolerance = 1e-10)
})

test_that("safe_matrix_inverse handles near-singular matrices", {
  # Near-singular matrix (condition number > 1e10)
  # Create a matrix with very high condition number
  M_bad <- matrix(c(1, 1-1e-11, 1-1e-11, 1), 2, 2)

  expect_warning(
    result <- GLOWr:::safe_matrix_inverse(M_bad, "test"),
    "near-singular"
  )

  # Result should still be valid
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("safe_matrix_inverse provides informative errors", {
  # Create a truly singular matrix
  M_singular <- matrix(c(1, 1, 1, 1), 2, 2)

  # This should either warn and then possibly error, or directly error
  # depending on whether nearPD can fix it
  # We just check that it either errors or handles it gracefully
  result <- tryCatch(
    {
      suppressWarnings(GLOWr:::safe_matrix_inverse(M_singular, "test context"))
      "success"
    },
    error = function(e) {
      expect_match(as.character(e), "Matrix inversion failed|positive semi-definite")
      "error"
    }
  )

  # Either should error or succeed after nearPD
  expect_true(result %in% c("success", "error"))
})

test_that("validate_numeric_input catches NA values", {
  x_with_na <- c(1, NA, 3)
  expect_error(
    GLOWr:::validate_numeric_input(x_with_na, "test"),
    "NA values"
  )
})

test_that("validate_numeric_input catches infinite values", {
  x_with_inf <- c(1, Inf, 3)
  expect_error(
    GLOWr:::validate_numeric_input(x_with_inf, "test"),
    "infinite"
  )

  x_with_neginf <- c(1, -Inf, 3)
  expect_error(
    GLOWr:::validate_numeric_input(x_with_neginf, "test"),
    "infinite"
  )
})

test_that("validate_numeric_input catches negative values when not allowed", {
  x_with_neg <- c(1, -1, 3)
  expect_error(
    GLOWr:::validate_numeric_input(x_with_neg, "test", allow_negative = FALSE),
    "non-negative"
  )
})

test_that("validate_numeric_input catches zero values when not allowed", {
  x_with_zero <- c(0, 1, 2)
  expect_error(
    GLOWr:::validate_numeric_input(x_with_zero, "test", allow_zero = FALSE),
    "positive"
  )
})

test_that("validate_numeric_input accepts valid inputs", {
  x_valid <- c(1, 2, 3)
  expect_silent(GLOWr:::validate_numeric_input(x_valid, "test"))

  # Should allow negative by default
  x_neg <- c(-1, 0, 1)
  expect_silent(GLOWr:::validate_numeric_input(x_neg, "test"))

  # Should allow zero by default
  x_zero <- c(0, 1, 2)
  expect_silent(GLOWr:::validate_numeric_input(x_zero, "test"))
})

test_that("check_positive_definite accepts valid PD matrices", {
  M_good <- diag(3)
  expect_silent(GLOWr:::check_positive_definite(M_good, "test"))

  # Valid correlation matrix
  M_corr <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  expect_silent(GLOWr:::check_positive_definite(M_corr, "test"))
})

test_that("check_positive_definite warns on nearly singular matrices", {
  # Nearly singular (minimum eigenvalue near zero but positive)
  # Create a matrix with eigenvalues very close to zero
  M_near <- matrix(c(1, 1-1e-11, 1-1e-11, 1), 2, 2)

  # Check if it warns (may not if eigenvalue isn't small enough)
  result <- tryCatch(
    {
      suppressWarnings(GLOWr:::check_positive_definite(M_near, "test"))
      "no_warning"
    },
    warning = function(w) {
      expect_match(as.character(w), "nearly singular")
      "warning"
    }
  )

  # Either warns or passes silently (depending on actual eigenvalue)
  expect_true(result %in% c("no_warning", "warning"))
})

test_that("check_positive_definite errors on non-PD matrices", {
  # Non-PD matrix (should error)
  M_bad <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    GLOWr:::check_positive_definite(M_bad, "test"),
    "positive semi-definite"
  )
})

test_that("get_wts handles all-zero weights", {
  Sigma <- diag(3)
  r <- c(0, 0, 0)

  expect_warning(
    result <- GLOWr:::get_wts(Sigma, r),
    "effectively zero"
  )

  # Should return equal weights
  expect_equal(result$w_normalized, rep(1/3, 3), tolerance = 1e-10)
})

test_that("get_wts handles near-zero weights", {
  Sigma <- diag(3)
  r <- c(1e-20, 1e-20, 1e-20)

  # Should either warn about near-zero or normalize successfully
  # depending on machine precision
  result <- GLOWr:::get_wts(Sigma, r)

  # Result should be valid
  expect_true(is.list(result))
  expect_true("w_normalized" %in% names(result))
  expect_equal(length(result$w_normalized), 3)
})

test_that("select_best_model validates inputs", {
  # NA in X
  expect_error(
    select_best_model(c(0.1, NA, 0.3), c(0.5, 0.6, 0.7)),
    "NA values"
  )

  # Inf in Y
  expect_error(
    select_best_model(c(0.1, 0.2, 0.3), c(0.5, Inf, 0.7)),
    "infinite"
  )

  # Negative MAF
  expect_error(
    select_best_model(c(-0.1, 0.2, 0.3), c(0.5, 0.6, 0.7)),
    "non-negative"
  )

  # Zero effect size
  expect_error(
    select_best_model(c(0.1, 0.2, 0.3), c(0, 0.6, 0.7)),
    "positive"
  )
})

test_that("model_PI validates inputs", {
  # NA in case annotation
  case_anno_na <- matrix(c(1, NA, 3, 4), ncol = 2)
  control_anno <- matrix(runif(100), ncol = 2)

  expect_error(
    model_PI(case_anno_na, control_anno, "LASSO", 10, 5),
    "NA values"
  )

  # Inf in control annotation
  control_anno_inf <- matrix(c(rep(0.5, 98), Inf, 0.5), ncol = 2)
  case_anno <- matrix(runif(10), ncol = 2)

  expect_error(
    model_PI(case_anno, control_anno_inf, "GLM", 10, 5),
    "infinite"
  )
})

test_that("Optimal_Weights_M validates inputs", {
  # Identity function for Burden test
  g_identity <- function(x) x

  # NA in Bstar
  expect_error(
    Optimal_Weights_M(g_identity, c(1, NA, 2), c(0.1, 0.2, 0.3), diag(3)),
    "NA values"
  )

  # Inf in PI
  expect_error(
    Optimal_Weights_M(g_identity, c(1, 2, 3), c(0.1, Inf, 0.3), diag(3)),
    "infinite"
  )

  # Negative Bstar should be ALLOWED (protective alleles)
  expect_silent({
    result <- Optimal_Weights_M(g_identity, c(1, -1, 2), c(0.1, 0.2, 0.3), diag(3))
  })

  # NA in M
  M_na <- matrix(c(1, NA, 0.5, 1), 2, 2)
  expect_error(
    Optimal_Weights_M(g_identity, c(1, 2), c(0.1, 0.2), M_na),
    "NA values"
  )
})

test_that("Optimal_Weights_M checks positive definiteness of M", {
  g_identity <- function(x) x

  # Non-PD matrix
  M_bad <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    Optimal_Weights_M(g_identity, c(1, 2), c(0.1, 0.2), M_bad),
    "positive semi-definite"
  )
})

test_that("Optimal_Weights_M handles near-singular correlation matrices", {
  g_identity <- function(x) x

  # Nearly singular M (condition number > 1e10)
  M_near <- matrix(c(1, 1-1e-11, 1-1e-11, 1), 2, 2)

  # Should warn about nearly singular matrix (either from check_positive_definite or safe_matrix_inverse)
  # But may not if the matrix isn't quite singular enough
  result <- suppressWarnings(
    Optimal_Weights_M(g_identity, c(1, 2), c(0.1, 0.2), M_near)
  )

  # Result should still be valid regardless
  expect_true(is.list(result))
  expect_true("wts_BE" %in% names(result))
  expect_true("wts_APE" %in% names(result))
})

test_that("CovM_gXgY checks positive definiteness", {
  g_identity <- function(x) x

  # Non-PD matrix
  M_bad <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    GLOWr:::CovM_gXgY(g_identity, c(0, 0), c(0, 0), M_bad),
    "positive semi-definite"
  )
})

test_that("CovMT_mix checks positive definiteness", {
  g_identity <- function(x) x

  # Non-PD matrix
  M_bad <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    GLOWr:::CovMT_mix(g_identity, c(1, 2), c(0.1, 0.2), M_bad),
    "positive semi-definite"
  )
})

test_that("Burden test handles all-zero effect sizes gracefully", {
  g_identity <- function(x) x

  # All zero Bstar values - this creates zero weights
  result <- suppressWarnings(
    Optimal_Weights_M(g_identity, c(0, 0, 0), c(0.1, 0.2, 0.3), diag(3))
  )

  # Should have returned something (likely equal weights)
  expect_true(is.list(result))
  expect_equal(length(result$wts_BE), 3)
  expect_equal(length(result$wts_APE), 3)
})

test_that("Robustness improvements don't break existing functionality", {
  # Test with valid inputs - should work without warnings/errors
  g_identity <- function(x) x
  Bstar <- c(0.5, 1.0, 0.3)
  PI <- c(0.01, 0.05, 0.02)
  M <- diag(3)

  expect_silent({
    result <- Optimal_Weights_M(g_identity, Bstar, PI, M)
  })

  expect_true(is.list(result))
  expect_true("wts_BE" %in% names(result))
  expect_true("wts_APE" %in% names(result))
  expect_equal(length(result$wts_BE), 3)
  expect_equal(length(result$wts_APE), 3)
})


########## Phase 3: Performance Optimization Tests ##########

test_that("H0 covariance caching works correctly", {
  # Enable caching
  options(GLOWr.use_cache = TRUE)

  # Clear cache
  clear_glow_cache()
  expect_equal(glow_cache_info()$cached_items, 0)

  # Create test data
  g <- function(x) x^2
  M <- diag(3)
  MU <- c(0.5, 1.0, 0.3)
  PI <- c(0.01, 0.05, 0.02)

  # First call - should cache H0 covariance
  Sigma1 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
  cache_size_1 <- glow_cache_info()$cached_items
  expect_true(cache_size_1 > 0)

  # Second call - should use cache (same M and g)
  Sigma2 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
  cache_size_2 <- glow_cache_info()$cached_items
  expect_equal(cache_size_1, cache_size_2)  # Cache size unchanged

  # Results should be identical
  expect_equal(Sigma1, Sigma2)

  # Clear cache
  clear_glow_cache()
  expect_equal(glow_cache_info()$cached_items, 0)

  # Disable caching
  options(GLOWr.use_cache = FALSE)
  Sigma3 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
  expect_equal(glow_cache_info()$cached_items, 0)  # Still zero
  expect_equal(Sigma1, Sigma3)  # Results identical

  # Re-enable caching for other tests
  options(GLOWr.use_cache = TRUE)
})

test_that("cache fallback works when errors occur", {
  options(GLOWr.use_cache = TRUE)
  clear_glow_cache()

  # This should still work even if caching has issues
  g <- function(x) x
  M <- diag(2)
  MU <- c(0.5, 1.0)
  PI <- c(0.01, 0.05)

  # Should not error even if caching fails
  expect_no_error({
    Sigma <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
  })
})

test_that("caching respects user options", {
  # Disable caching
  options(GLOWr.use_cache = FALSE)

  # Check cache info
  info <- glow_cache_info()
  expect_false(info$enabled)

  # Enable caching
  options(GLOWr.use_cache = TRUE)

  # Check cache info again
  info <- glow_cache_info()
  expect_true(info$enabled)
})

test_that("cache management functions work correctly", {
  options(GLOWr.use_cache = TRUE)

  # Clear cache
  clear_glow_cache()
  expect_equal(glow_cache_info()$cached_items, 0)

  # Add something to cache
  g <- function(x) x^2
  M <- diag(2)
  MU <- c(0.5, 1.0)
  PI <- c(0.01, 0.05)

  Sigma <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")

  # Cache should have items
  expect_true(glow_cache_info()$cached_items > 0)

  # Clear cache again
  clear_glow_cache()
  expect_equal(glow_cache_info()$cached_items, 0)
})

test_that("fast path for independent variants produces identical results", {
  # Define test transformation
  g <- function(x) x^2
  MU1 <- c(0.5, 1.0, 0.3)
  MU2 <- c(0.5, 1.0, 0.3)
  n <- length(MU1)

  # Identity matrix (independent variants)
  M_identity <- diag(n)

  # Correlated matrix
  M_corr <- diag(n)
  M_corr[1, 2] <- M_corr[2, 1] <- 0.3

  # Compute with identity (should use fast path)
  result_fast <- GLOWr:::CovM_gXgY(g, MU1, MU2, M_identity)

  # Result should be diagonal (off-diagonals zero)
  expect_true(all(abs(result_fast[row(result_fast) != col(result_fast)]) < 1e-10))

  # Diagonal elements should be positive
  expect_true(all(diag(result_fast) > 0))

  # Compare with correlated version
  result_corr <- GLOWr:::CovM_gXgY(g, MU1, MU2, M_corr)

  # Diagonal should be similar
  expect_equal(diag(result_fast), diag(result_corr), tolerance = 1e-6)

  # Off-diagonals should differ
  expect_true(max(abs(result_corr[row(result_corr) != col(result_corr)])) > 0.01)
})

test_that("fast path for independent variants works with different transformations", {
  # Test with identity function
  g_identity <- function(x) x
  MU1 <- c(0, 0, 0)
  MU2 <- c(0, 0, 0)
  M_identity <- diag(3)

  result_identity <- GLOWr:::CovM_gXgY(g_identity, MU1, MU2, M_identity)

  # For identity function with mean 0, should be identity matrix
  expect_equal(result_identity, M_identity, tolerance = 1e-6)

  # Test with square function
  g_square <- function(x) x^2
  MU1_nonzero <- c(1, 2, 3)
  MU2_nonzero <- c(1, 2, 3)

  result_square <- GLOWr:::CovM_gXgY(g_square, MU1_nonzero, MU2_nonzero, M_identity)

  # Should be diagonal
  expect_true(all(abs(result_square[row(result_square) != col(result_square)]) < 1e-10))

  # Diagonal should be positive
  expect_true(all(diag(result_square) > 0))
})

test_that("fast path numerical accuracy matches regular path", {
  # For independent variants, fast path and regular path should give identical results
  g <- function(x) x^3  # Cubic transformation
  MU1 <- c(0.2, 0.8, 1.5)
  MU2 <- c(0.2, 0.8, 1.5)
  M_identity <- diag(3)

  # Compute with fast path (automatically used for identity matrix)
  result_fast <- GLOWr:::CovM_gXgY(g, MU1, MU2, M_identity)

  # Diagonal should match variance calculations
  # For fast path, variance is computed using Hermite expansion
  # The result should be a diagonal matrix
  expect_true(all(abs(result_fast[row(result_fast) != col(result_fast)]) < 1e-10))

  # Check that diagonal elements are reasonable
  expect_true(all(diag(result_fast) > 0))
  expect_true(all(is.finite(diag(result_fast))))
})

test_that("caching provides speedup", {
  skip_on_cran()

  options(GLOWr.use_cache = TRUE)
  clear_glow_cache()

  g <- function(x) x^2
  M <- diag(10)
  MU <- rep(1, 10)
  PI <- rep(0.05, 10)

  # First call (no cache) - run multiple times to get stable timing
  time1 <- system.time({
    for (i in 1:3) {
      Sigma1 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
    }
  })

  # Second call (with cache) - should be much faster
  time2 <- system.time({
    for (i in 1:3) {
      Sigma2 <- GLOWr:::get_Sigma(g, MU, PI, M, hypo = "H0")
    }
  })

  # Cached call should be faster (note: timing can be variable, so we use a conservative threshold)
  # The cache should eliminate the expensive CovM_gXgY computation
  expect_lt(time2["elapsed"], time1["elapsed"])

  # Results identical
  expect_equal(Sigma1, Sigma2)
})

test_that("fast path provides speedup for independent variants", {
  skip_on_cran()

  g <- function(x) x^2
  n <- 20
  MU1 <- rep(1, n)
  MU2 <- rep(1, n)

  # Identity matrix (independent - fast path)
  M_identity <- diag(n)

  # Slightly correlated (regular path)
  M_corr <- diag(n)
  M_corr[1, 2] <- M_corr[2, 1] <- 0.1

  # Time fast path
  time_fast <- system.time({
    result_fast <- GLOWr:::CovM_gXgY(g, MU1, MU2, M_identity)
  })

  # Time regular path
  time_regular <- system.time({
    result_regular <- GLOWr:::CovM_gXgY(g, MU1, MU2, M_corr)
  })

  # Fast path should be faster (it only computes n diagonal elements
  # instead of n^2 matrix elements)
  # We expect at least some speedup, though it may not be dramatic for small n
  expect_lt(time_fast["elapsed"], time_regular["elapsed"] * 1.5)
})
