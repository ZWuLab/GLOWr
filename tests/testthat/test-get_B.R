########## Tests for get_B() Function ##########
#
# This file contains comprehensive tests for the get_B() function, including:
# - Tests with external summary statistics (same trait)
# - Tests with external summary statistics (different traits)
# - Model selection validation (via select_best_model)
# - Comparison with legacy implementation (tolerance < 1e-10)
# - Edge cases: missing data, single variant, extreme MAF values

library(testthat)

# ========== Helper Functions ==========

#' Source legacy get_B for comparison
source_legacy_get_B <- function() {
  legacy_file <- file.path(
    Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"),
    "legacy-materials/code/GLOW_R_pacakge/GLOW/R/get_B.R"
  )

  if (!file.exists(legacy_file)) {
    skip("Legacy get_B.R not found")
  }

  # Source the legacy file
  source(legacy_file, local = TRUE)

  # Also need legacy select_best_model
  legacy_helpers <- file.path(
    Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"),
    "legacy-materials/code/GLOW_R_pacakge/GLOW/R/helpers_optimalWeights.R"
  )

  if (!file.exists(legacy_helpers)) {
    skip("Legacy helpers_optimalWeights.R not found")
  }

  source(legacy_helpers, local = TRUE)

  # Return the function
  return(get_B)
}


# ========== Basic Functionality Tests ==========

test_that("get_B works with same binary trait", {
  set.seed(123)

  # Generate training data
  training_MAF <- runif(50, 0.0002, 0.3)
  training_BETA <- -log(training_MAF) / 10

  # Generate testing data
  target_MAF <- runif(10, 0.0002, 0.5)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  # Check output properties
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
  expect_true(all(B_estimates > 0))  # Effect sizes should be positive
})


test_that("get_B works with same continuous trait", {
  set.seed(456)

  # Generate training data
  training_MAF <- runif(100, 0.001, 0.4)
  training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF)) + rnorm(100, 0, 0.01)

  # Generate testing data
  target_MAF <- runif(20, 0.001, 0.5)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "continuous",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  # Check output properties
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
  expect_true(all(B_estimates > 0))
})


test_that("get_B works with different traits (continuous to binary)", {
  set.seed(789)

  # Generate training data (continuous trait)
  training_MAF <- runif(50, 0.0002, 0.3)
  training_N <- rep(1000, 50)
  training_Z_temp <- rnorm(50, mean = sqrt(2 * training_N * training_MAF *
                                       (1 - training_MAF) * 0.01))
  training_P <- pchisq(training_Z_temp^2, df = 1, lower.tail = FALSE)

  # Generate testing data (binary trait)
  target_MAF <- runif(10, 0.0002, 0.5)
  target_case_prop <- rep(0.1, 10)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop
  )

  # Check output properties
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
  expect_true(all(B_estimates > 0))
})


test_that("get_B works with different traits (binary to continuous)", {
  set.seed(321)

  # Generate training data (binary trait)
  training_MAF <- runif(50, 0.0002, 0.3)
  training_N <- rep(1000, 50)
  training_Z_temp <- rnorm(50, mean = sqrt(2 * training_N * training_MAF *
                                       (1 - training_MAF) * 0.01))
  training_P <- pchisq(training_Z_temp^2, df = 1, lower.tail = FALSE)

  # Generate testing data (continuous trait)
  target_MAF <- runif(10, 0.0002, 0.5)
  target_SE <- runif(10, 0.01, 0.1)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "continuous",
    target_MAF = target_MAF,
    target_SE = target_SE
  )

  # Check output properties
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
  expect_true(all(B_estimates > 0))
})


# ========== Comparison with Legacy Implementation ==========

test_that("get_B matches legacy implementation (same trait)", {
  skip_on_cran()

  # Source legacy function
  legacy_get_B <- source_legacy_get_B()

  set.seed(12345)

  # Generate training data
  training_MAF <- runif(50, 0.0002, 0.3)
  training_BETA <- -log(training_MAF) / 10

  # Generate testing data
  target_MAF <- runif(10, 0.0002, 0.5)

  # New implementation
  B_new <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  # Legacy implementation with error handling
  B_legacy <- tryCatch({
    legacy_get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
    show_model_selection = FALSE
    )
  }, error = function(e) {
    skip(paste("Legacy get_B failed:", e$message))
  })

  # Check for NaN values in legacy result
  if (any(is.nan(B_legacy))) {
    skip("Legacy implementation produced NaN values for this test case")
  }

  # Compare results (tolerance < 1e-10)
  expect_equal(B_new, B_legacy, tolerance = 1e-10)
})


test_that("get_B matches legacy implementation (different traits)", {
  skip_on_cran()

  # Source legacy function
  legacy_get_B <- source_legacy_get_B()

  set.seed(54321)

  # Generate training data (continuous)
  training_MAF <- runif(50, 0.0002, 0.3)
  training_N <- rep(1000, 50)
  training_Z_temp <- rnorm(50, mean = sqrt(2 * training_N * training_MAF *
                                       (1 - training_MAF) * 0.01))
  training_P <- pchisq(training_Z_temp^2, df = 1, lower.tail = FALSE)

  # Generate testing data (binary)
  target_MAF <- runif(10, 0.0002, 0.5)
  target_case_prop <- rep(0.1, 10)

  # New implementation
  B_new <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop
  )

  # Legacy implementation with error handling
  B_legacy <- tryCatch({
    legacy_get_B(
      training_trait = "continuous",
      training_MAF = training_MAF,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop
    )
  }, error = function(e) {
    skip(paste("Legacy get_B failed:", e$message))
  })

  # Check for NaN values in legacy result
  if (any(is.nan(B_legacy))) {
    skip("Legacy implementation produced NaN values for this test case")
  }

  # Compare results (tolerance < 1e-10)
  expect_equal(B_new, B_legacy, tolerance = 1e-10)
})


test_that("get_B matches legacy implementation (binary to continuous)", {
  skip_on_cran()

  # Source legacy function
  legacy_get_B <- source_legacy_get_B()

  set.seed(99999)

  # Generate training data (binary)
  training_MAF <- runif(50, 0.0002, 0.3)
  training_N <- rep(1000, 50)
  training_Z_temp <- rnorm(50, mean = sqrt(2 * training_N * training_MAF *
                                       (1 - training_MAF) * 0.01))
  training_P <- pchisq(training_Z_temp^2, df = 1, lower.tail = FALSE)

  # Generate testing data (continuous)
  target_MAF <- runif(10, 0.0002, 0.5)
  target_SE <- runif(10, 0.01, 0.1)

  # New implementation
  B_new <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "continuous",
    target_MAF = target_MAF,
    target_SE = target_SE
  )

  # Legacy implementation with error handling
  B_legacy <- tryCatch({
    legacy_get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_P = training_P,
      training_N = training_N,
      target_trait = "continuous",
      target_MAF = target_MAF,
      target_SE = target_SE
    )
  }, error = function(e) {
    skip(paste("Legacy get_B failed:", e$message))
  })

  # Check for NaN values in legacy result
  if (any(is.nan(B_legacy))) {
    skip("Legacy implementation produced NaN values for this test case")
  }

  # Compare results (tolerance < 1e-10)
  expect_equal(B_new, B_legacy, tolerance = 1e-10)
})


# ========== Input Validation Tests ==========

test_that("get_B validates trait types", {
  training_MAF <- runif(10, 0.01, 0.5)
  training_BETA <- runif(10, 0, 1)
  target_MAF <- runif(5, 0.01, 0.5)

  # Invalid training_trait
  expect_error(
    get_B(
      training_trait = "invalid",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
    show_model_selection = FALSE
    ),
    "training_trait must be 'binary', 'continuous', 'mixed', or NULL"
  )

  # Invalid target_trait
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "invalid",
      target_MAF = target_MAF,
    show_model_selection = FALSE
    ),
    "target_trait must be either 'binary', 'continuous', or NULL"
  )
})


test_that("get_B validates MAF ranges", {
  # training_MAF with 0 and 1: after folding 1->0, both 0s fail (0, 0.5] check
  training_MAF <- c(0, 0.1, 0.5, 1)  # Invalid: contains 0 (and 1 folds to 0)
  training_BETA <- runif(4, 0, 1)
  target_MAF <- runif(5, 0.01, 0.5)

  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
    show_model_selection = FALSE
    ),
    "training_MAF values must be in \\(0, 0\\.5\\]"
  )

  # target_MAF with -0.1 and 1.1: 1.1 folds to -0.1, both negative -> error
  training_MAF <- runif(10, 0.01, 0.5)
  target_MAF <- c(-0.1, 0.5, 1.1)  # Invalid: -0.1 negative, 1.1 folds to -0.1

  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
    show_model_selection = FALSE
    ),
    "target_MAF values must be in \\(0, 0\\.5\\]"
  )
})

test_that("get_B folds MAF > 0.5 with warning", {
  set.seed(42)
  # Training MAF with values > 0.5 should be folded
  training_MAF <- c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9, 0.15, 0.25, 0.35, 0.05)
  training_BETA <- sqrt(pmin(training_MAF, 1 - training_MAF)) * 0.1

  # target_MAF with value > 0.5 should also fold
  target_MAF <- c(0.1, 0.7, 0.3)

  # Should get warnings for both training and target MAF folding
  expect_warning(
    B <- get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      show_model_selection = FALSE,
      verbose = 0
    ),
    "value\\(s\\) > 0\\.5 detected"
  )
  expect_true(is.numeric(B))
  expect_equal(length(B), 3)
})


test_that("get_B requires training_BETA when traits are same", {
  training_MAF <- runif(10, 0.01, 0.5)
  target_MAF <- runif(5, 0.01, 0.5)

  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      # training_BETA missing
      target_trait = "binary",
      target_MAF = target_MAF
    ),
    "Insufficient data: need either training_BETA \\(for beta method\\) or training_P/training_P_mlog10 \\(for pvalue method\\)"
  )
})


test_that("get_B requires training_P and training_N when traits differ", {
  training_MAF <- runif(10, 0.01, 0.5)
  target_MAF <- runif(5, 0.01, 0.5)
  target_case_prop <- rep(0.1, 5)

  # Missing both P and N
  expect_error(
    get_B(
      training_trait = "continuous",
      training_MAF = training_MAF,
      # training_P missing
      # training_N missing
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop
    ),
    "Insufficient data: need either training_BETA \\(for beta method\\) or training_P/training_P_mlog10 \\(for pvalue method\\)"
  )

  # Missing N (has training_P but not training_N)
  expect_error(
    get_B(
      training_trait = "continuous",
      training_MAF = training_MAF,
      training_P = runif(10, 0.0001, 0.1),
      # training_N missing
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop
    ),
    "training_N is required for method 'pvalue'"
  )
})


test_that("get_B requires target_SE for continuous target trait (when traits differ)", {
  training_MAF <- runif(10, 0.01, 0.5)
  training_P <- runif(10, 0.0001, 0.1)
  training_N <- rep(1000, 10)
  target_MAF <- runif(5, 0.01, 0.5)

  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_P = training_P,
      training_N = training_N,
      target_trait = "continuous",
      target_MAF = target_MAF
      # target_SE missing
    ),
    "target_SE is required for continuous target trait with p-value method"
  )
})


test_that("get_B requires target_case_prop for binary target trait (when traits differ)", {
  training_MAF <- runif(10, 0.01, 0.5)
  training_P <- runif(10, 0.0001, 0.1)
  training_N <- rep(1000, 10)
  target_MAF <- runif(5, 0.01, 0.5)

  expect_error(
    get_B(
      training_trait = "continuous",
      training_MAF = training_MAF,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF
      # target_case_prop missing
    ),
    "target_case_prop is required for binary target trait with p-value method"
  )
})


test_that("get_B validates vector lengths", {
  training_MAF <- runif(10, 0.01, 0.5)
  training_BETA <- runif(5, 0, 1)  # Wrong length
  target_MAF <- runif(5, 0.01, 0.5)

  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
    show_model_selection = FALSE
    ),
    "training_BETA and training_MAF must have the same length"
  )
})


test_that("get_B validates prevalence range", {
  training_MAF <- runif(10, 0.01, 0.5)
  training_P <- runif(10, 0.0001, 0.1)
  training_N <- rep(1000, 10)
  target_MAF <- runif(5, 0.01, 0.5)
  target_case_prop <- c(0, 0.5, 1.1, 0.3, 0.2)  # Invalid: contains 0 and > 1, but same length

  expect_error(
    get_B(
      training_trait = "continuous",
      training_MAF = training_MAF,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop
    ),
    "target_case_prop values must be in the interval \\(0, 1\\)"
  )
})


# ========== Edge Cases ==========

test_that("get_B requires minimum training sample size", {
  set.seed(111)

  # Single training variant - too few for model fitting
  training_MAF <- 0.1
  training_BETA <- 0.5
  target_MAF <- 0.2

  # Should fail with informative error
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      show_model_selection = FALSE
    ),
    "Need at least 3 observations to fit models"
  )

  # With 3 training variants, should work
  training_MAF_min <- c(0.1, 0.2, 0.3)
  training_BETA_min <- c(0.5, 0.6, 0.4)

  B_estimate <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF_min,
    training_BETA = training_BETA_min,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  expect_type(B_estimate, "double")
  expect_length(B_estimate, 1)
  expect_true(is.finite(B_estimate))
})


test_that("get_B works with very small MAF values", {
  set.seed(222)

  # Very rare variants
  training_MAF <- runif(50, 1e-5, 1e-3)
  training_BETA <- -log(training_MAF) / 10

  target_MAF <- runif(10, 1e-5, 1e-3)

  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
})


test_that("get_B works with common variants", {
  set.seed(333)

  # Common variants
  training_MAF <- runif(50, 0.3, 0.49)
  training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF))

  target_MAF <- runif(10, 0.3, 0.49)

  B_estimates <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "continuous",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
})


test_that("get_B handles zero-valued BETA elements correctly", {
  set.seed(444)

  # Create training data with some zero-valued BETA
  n_train <- 50
  training_MAF <- runif(n_train, 0.01, 0.4)
  training_BETA <- -log(training_MAF) / 10

  # Set some BETA values to zero (simulating OR = 1.0 cases)
  zero_indices <- c(5, 15, 25)
  training_BETA[zero_indices] <- 0

  # Also need P and N for the p-value method
  training_P <- runif(n_train, 0.0001, 0.1)
  training_N <- rep(1000, n_train)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.5, 10)

  # Test with verbose output to verify zero-handling messages
  # Suppress messages but verify function works
  B_estimates <- suppressMessages(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "both",  # Test both methods
      verbose = 2,
      show_model_selection = FALSE
    )
  )

  # Verify the function handled zeros correctly
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
  expect_true(all(B_estimates > 0))

  # Test that it produces a message when verbose >= 1
  expect_message(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "both",
      show_model_selection = FALSE,
    verbose = 1,
    ),
    "zero-valued BETA"
  )
})


test_that("get_B errors when all BETA values are zero", {
  set.seed(445)

  # Create training data with all zero BETA values
  training_MAF <- runif(20, 0.01, 0.4)
  training_BETA <- rep(0, 20)

  target_MAF <- runif(5, 0.01, 0.4)

  # Should error because all BETA values are zero
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      method = "beta"
    ),
    "All training_BETA values are zero"
  )
})


# ========== Model Selection Integration ==========

test_that("get_B integrates with select_best_model correctly", {
  set.seed(444)

  # Create data where we know which model should be best
  # Use a simple linear relationship: Y ~ X
  training_MAF <- seq(0.01, 0.5, length.out = 50)
  training_BETA <- sqrt(2 * training_MAF)  # B^2 = 2*MAF, so Y = B^2 = 2*MAF

  target_MAF <- seq(0.05, 0.45, length.out = 10)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  # Expected B values (approximately)
  B_expected <- sqrt(2 * target_MAF)

  # The estimates should be close to expected (within 10%)
  expect_true(all(abs(B_estimates - B_expected) / B_expected < 0.1))
})


test_that("get_B handles log transformation correctly", {
  set.seed(555)

  # Create data where log transformation is better
  # Use relationship: log(Y) ~ log(X)
  training_MAF <- runif(100, 0.01, 0.5)
  # B^2 = 0.1 / MAF^0.5, so Y = B^2
  training_BETA <- sqrt(0.1 / sqrt(training_MAF))

  target_MAF <- runif(20, 0.01, 0.5)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "continuous",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  # Check that estimates are positive and finite
  expect_true(all(B_estimates > 0))
  expect_true(all(is.finite(B_estimates)))

  # The relationship should hold approximately
  B_expected <- sqrt(0.1 / sqrt(target_MAF))

  # Should be within 20% (log models can have more variability)
  expect_true(all(abs(B_estimates - B_expected) / B_expected < 0.2))
})


# ========== Integration with Optimal_Weights_M ==========

test_that("get_B output can be used in Optimal_Weights_M", {
  skip_on_cran()

  set.seed(666)

  # Generate training data
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- -log(training_MAF) / 10

  # Generate testing data
  target_MAF <- runif(10, 0.01, 0.3)

  # Estimate B
  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  # Create mock PI and M
  PI_estimates <- runif(10, 0.001, 0.1)
  M <- diag(10)  # Independent variants

  # Test with burden (identity function)
  g_identity <- function(x) x

  # This should work without error
  weights <- Optimal_Weights_M(g_identity, B_estimates, PI_estimates, M)

  expect_type(weights, "list")
  expect_true("wts_BE" %in% names(weights))
  expect_true("wts_APE" %in% names(weights))
  expect_length(weights$wts_BE, 10)
  expect_length(weights$wts_APE, 10)
})


# ========== Numerical Stability Tests ==========

test_that("get_B is numerically stable with extreme effect sizes", {
  set.seed(777)

  # Very large effect sizes
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- runif(50, 5, 10)  # Large effects

  target_MAF <- runif(10, 0.01, 0.3)

  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  expect_true(all(is.finite(B_estimates)))
  expect_true(all(B_estimates > 0))

  # Very small effect sizes
  training_BETA <- runif(50, 0.001, 0.01)  # Small effects

  B_estimates_small <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  expect_true(all(is.finite(B_estimates_small)))
  expect_true(all(B_estimates_small > 0))
})


# ========== Outlier Detection Integration Tests ==========

test_that("get_B integrates outlier detection correctly with flag action", {
  # Create training data with a clear outlier
  training_MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  training_BETA <- c(0.5, 0.3, 0.2, 0.15, 20, 0.08)  # index 5 is an outlier for MAF=0.3
  target_MAF <- c(0.15, 0.25)

  # Test with outlier detection (flag only)
  result_flagged <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    outlier_method = "biological",
    outlier_action = "flag",
    return_full = TRUE,
    verbose = 0,
  )

  expect_s3_class(result_flagged, "glow_B_estimate")
  expect_true(!is.null(result_flagged$model$outliers$beta_method))
  expect_true(length(result_flagged$model$outliers$beta_method$indices) > 0)
  expect_equal(result_flagged$model$training_summary$n_used, 6)  # All data used
  expect_equal(length(result_flagged$model$outliers$indices_removed), 0)  # None removed
})


test_that("get_B integrates outlier detection correctly with remove action", {
  # Create training data with a clear outlier
  training_MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  training_BETA <- c(0.5, 0.3, 0.2, 0.15, 20, 0.08)  # index 5 is an outlier
  target_MAF <- c(0.15, 0.25)

  # Test with outlier removal
  result_removed <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    outlier_method = "biological",
    outlier_action = "remove",
    return_full = TRUE,
    verbose = 0,
  )

  expect_equal(result_removed$model$training_summary$n_used, 5)  # One outlier removed
  expect_true(5 %in% result_removed$model$outliers$beta_method$indices)  # Index 5 is the outlier
  expect_equal(length(result_removed$model$outliers$indices_removed), 1)  # One removed
  expect_true(5 %in% result_removed$model$outliers$indices_removed)
})


test_that("get_B return_full includes comprehensive diagnostics", {
  training_MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3)
  training_BETA <- c(0.5, 0.3, 0.2, 0.15, 0.1)
  target_MAF <- c(0.15, 0.25)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0,
  )

  # Check structure (new glow_B_estimate nests model info under $model)
  expect_true(all(c("B", "B_beta_method", "B_pvalue_method",
                    "model", "target_summary") %in% names(result)))

  # Check training summary (nested under $model)
  expect_equal(result$model$training_summary$n_original, 5)
  expect_equal(result$model$training_summary$trait_type, "binary")

  # Check target summary
  expect_equal(result$target_summary$n_predictions, 2)
  expect_equal(result$target_summary$trait_type, "binary")

  # Check outlier info (nested under $model, should be present even if none detected)
  expect_true(!is.null(result$model$outliers))
  expect_equal(result$model$outliers$method, "none")  # Default
  expect_equal(result$model$outliers$action, "flag")  # Default
})


test_that("get_B handles no biological outliers", {
  # Data with no biological outliers (rare variants can have large effects)
  training_MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  training_BETA <- c(0.5, 0.3, 0.2, 0.15, 0.1, 0.08)
  target_MAF <- c(0.15, 0.25)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    outlier_method = "biological",  # Only biological, not statistical
    return_full = TRUE,
    verbose = 0,
  )

  # Should detect no biological outliers (rare variants with large effects are plausible)
  expect_equal(result$model$training_summary$n_outliers_detected, 0)
  expect_equal(result$model$training_summary$n_used, 6)
})


test_that("get_B outlier detection works with pvalue method", {
  set.seed(123)
  # Create training data with outlier
  # Need to create h² > 10 for outlier detection
  # h² = Z² / (2 * N * q * (1-q))
  # For MAF=0.3, N=10: h² = Z² / (2*10*0.3*0.7) = Z²/4.2
  # To get h² > 10, need Z² > 42
  training_MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  training_N <- c(1000, 1000, 1000, 1000, 10, 1000)  # Small N for index 5
  # P-value that gives Z² ~100: pchisq(100, 1, lower.tail=FALSE) ≈ 1.6e-23
  training_P <- c(0.01, 0.05, 0.1, 0.2, 1e-23, 0.3)  # index 5 is outlier
  target_MAF <- c(0.15, 0.25)
  target_case_prop <- rep(0.1, 2)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    outlier_method = "biological",
    outlier_action = "remove",
    return_full = TRUE,
    verbose = 0,
  )

  # Should detect outlier (common variant with unrealistically large h^2)
  expect_true(result$model$training_summary$n_outliers_detected > 0)
  expect_true(!is.null(result$model$outliers$pvalue_method))
})


test_that("get_B outlier detection works with both methods", {
  # Create training data with outlier
  training_MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
  training_BETA <- c(0.5, 0.3, 0.2, 0.15, 20, 0.08)  # index 5 is outlier
  training_N <- rep(1000, 6)
  training_P <- c(0.01, 0.05, 0.1, 0.2, 1e-10, 0.3)  # index 5 is outlier
  target_MAF <- c(0.15, 0.25)
  target_case_prop <- rep(0.1, 2)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "both",
    outlier_method = "biological",
    outlier_action = "flag",
    return_full = TRUE,
    verbose = 0,
    show_model_selection = FALSE
  )

  # Should have outlier info for both methods (nested under $model)
  expect_true(!is.null(result$model$outliers$beta_method))
  expect_true(!is.null(result$model$outliers$pvalue_method))
})


test_that("get_B validates outlier parameters correctly", {
  training_MAF <- runif(10, 0.01, 0.5)
  training_BETA <- runif(10, 0, 1)
  target_MAF <- runif(5, 0.01, 0.5)

  # Invalid outlier_method
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      outlier_method = "invalid"
    ),
    "outlier_method must be one of"
  )

  # Invalid outlier_action
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      outlier_action = "invalid"
    ),
    "outlier_action must be one of"
  )

  # Invalid cook_threshold
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      cook_threshold = -1
    ),
    "cook_threshold must be a positive"
  )
})


test_that("get_B backward compatibility - default outlier_method is none", {
  # Existing code should work without specifying outlier parameters
  set.seed(123)
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- -log(training_MAF) / 10
  target_MAF <- runif(10, 0.01, 0.5)

  # Should work with no outlier parameters specified
  B_estimates <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    show_model_selection = FALSE
  )

  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))

  # return_full should still work
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    return_full = TRUE,
    verbose = 0,
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_equal(result$model$outliers$method, "none")
})


# ========== Smart Primary Method Selection Tests (method="both") ==========

test_that("get_B selects primary method based on R2 criterion when method='both'", {
  set.seed(1001)

  # Create data where beta method should have higher R2
  # Linear relationship: BETA^2 ~ MAF
  training_MAF <- seq(0.01, 0.4, length.out = 50)
  training_BETA <- sqrt(2 * training_MAF)  # Strong linear relationship
  training_N <- rep(1000, 50)
  # Add noise to p-values so pvalue method has lower R2
  training_P <- runif(50, 0.001, 0.1)

  target_MAF <- seq(0.05, 0.35, length.out = 10)
  target_case_prop <- rep(0.5, 10)

  # Run with both methods and R2 criterion
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "both",
    selection_criterion = "R2",
    return_full = TRUE,
    verbose = 0,
      show_model_selection = FALSE,
  )

  # Check that comparison information is included (nested under $model)
  expect_true(!is.null(result$model$comparison))
  expect_true("method_selected" %in% names(result$model$comparison))
  expect_true("selection_criterion" %in% names(result$model$comparison))
  expect_true("criterion_beta_method" %in% names(result$model$comparison))
  expect_true("criterion_pvalue_method" %in% names(result$model$comparison))

  # Verify selection_criterion is recorded
  expect_equal(result$model$comparison$selection_criterion, "R2")

  # Verify method was selected (either beta_method or pvalue_method)
  expect_true(result$model$comparison$method_selected %in%
              c("beta_method", "pvalue_method"))

  # Primary B should match the selected method
  if (result$model$comparison$method_selected == "beta_method") {
    expect_equal(result$B, result$B_beta_method)
  } else {
    expect_equal(result$B, result$B_pvalue_method)
  }
})


test_that("get_B selects primary method based on adj_R2 criterion when method='both'", {
  set.seed(1002)

  training_MAF <- runif(60, 0.01, 0.4)
  training_BETA <- -log(training_MAF) / 10
  training_N <- rep(1000, 60)
  training_P <- runif(60, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.3, 10)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "both",
    selection_criterion = "adj_R2",
    return_full = TRUE,
    verbose = 0,
    show_model_selection = FALSE
  )

  # Check comparison structure (nested under $model)
  expect_equal(result$model$comparison$selection_criterion, "adj_R2")
  expect_true(result$model$comparison$method_selected %in%
              c("beta_method", "pvalue_method"))

  # For adj_R2, higher is better, so verify the selected method has higher adj_R2
  if (result$model$comparison$method_selected == "beta_method") {
    expect_true(result$model$comparison$criterion_beta_method >=
                result$model$comparison$criterion_pvalue_method)
    expect_equal(result$B, result$B_beta_method)
  } else {
    expect_true(result$model$comparison$criterion_pvalue_method >=
                result$model$comparison$criterion_beta_method)
    expect_equal(result$B, result$B_pvalue_method)
  }
})


test_that("get_B selects primary method based on CV_R2 criterion when method='both'", {
  set.seed(1003)

  training_MAF <- runif(40, 0.01, 0.4)
  training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF))
  training_N <- rep(2000, 40)
  training_P <- runif(40, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.2, 10)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "both",
    selection_criterion = "CV_R2",
    return_full = TRUE,
    verbose = 0,
    show_model_selection = FALSE
  )

  # Check comparison structure (nested under $model)
  expect_equal(result$model$comparison$selection_criterion, "CV_R2")
  expect_true(result$model$comparison$method_selected %in%
              c("beta_method", "pvalue_method"))

  # For CV_R2, higher is better
  if (result$model$comparison$method_selected == "beta_method") {
    expect_true(result$model$comparison$criterion_beta_method >=
                result$model$comparison$criterion_pvalue_method)
    expect_equal(result$B, result$B_beta_method)
  } else {
    expect_true(result$model$comparison$criterion_pvalue_method >=
                result$model$comparison$criterion_beta_method)
    expect_equal(result$B, result$B_pvalue_method)
  }
})


test_that("get_B selects primary method based on CV criterion when method='both'", {
  skip_on_cran()  # CV is slow

  set.seed(1004)

  training_MAF <- runif(30, 0.01, 0.4)  # Smaller sample for faster CV
  training_BETA <- -log(training_MAF) / 10
  training_N <- rep(1500, 30)
  training_P <- runif(30, 0.001, 0.1)

  target_MAF <- runif(5, 0.01, 0.4)
  target_case_prop <- rep(0.4, 5)

  result <- suppressMessages(  # CV may generate messages
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "both",
      selection_criterion = "CV",
      return_full = TRUE,
      verbose = 0,
        show_model_selection = FALSE,
    )
  )

  # Check comparison structure (nested under $model)
  # Note: "CV" is deprecated and converted to "CV_R2" internally
  expect_equal(result$model$comparison$selection_criterion, "CV_R2")
  expect_true(result$model$comparison$method_selected %in%
              c("beta_method", "pvalue_method"))

  # For CV_R2, higher is better
  if (result$model$comparison$method_selected == "beta_method") {
    expect_true(result$model$comparison$criterion_beta_method >=
                result$model$comparison$criterion_pvalue_method)
    expect_equal(result$B, result$B_beta_method)
  } else {
    expect_true(result$model$comparison$criterion_pvalue_method >=
                result$model$comparison$criterion_beta_method)
    expect_equal(result$B, result$B_pvalue_method)
  }
})


test_that("get_B smart selection provides informative message with verbose >= 1", {
  set.seed(1005)

  training_MAF <- runif(40, 0.01, 0.4)
  training_BETA <- sqrt(2 * training_MAF)
  training_N <- rep(1000, 40)
  training_P <- runif(40, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.5, 10)

  # Should get a message about which method was selected
  expect_message(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "both",
      selection_criterion = "R2",
      verbose = 1,
        show_model_selection = FALSE,
    ),
    "Selected .* as primary"
  )
})


test_that("get_B smart selection works with continuous traits", {
  set.seed(1006)

  training_MAF <- runif(50, 0.01, 0.4)
  training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF))
  training_N <- rep(2000, 50)
  training_P <- runif(50, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_SE <- runif(10, 0.01, 0.1)

  result <- get_B(
    training_trait = "continuous",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_N = training_N,
    target_trait = "continuous",
    target_MAF = target_MAF,
    target_SE = target_SE,
    method = "both",
    selection_criterion = "R2",
    return_full = TRUE,
    verbose = 0,
      show_model_selection = FALSE,
  )

  # Should have comparison information (nested under $model)
  expect_true(!is.null(result$model$comparison))
  expect_true("method_selected" %in% names(result$model$comparison))
  expect_true(result$model$comparison$method_selected %in%
              c("beta_method", "pvalue_method"))

  # Primary should match selected method
  if (result$model$comparison$method_selected == "beta_method") {
    expect_equal(result$B, result$B_beta_method)
  } else {
    expect_equal(result$B, result$B_pvalue_method)
  }
})


test_that("get_B smart selection backward compatible - return_full=FALSE works", {
  set.seed(1007)

  training_MAF <- runif(40, 0.01, 0.4)
  training_BETA <- -log(training_MAF) / 10
  training_N <- rep(1000, 40)
  training_P <- runif(40, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.5, 10)

  # With return_full=FALSE, should just return numeric vector
  B_estimates <- suppressMessages(  # Suppress selection messages
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "both",
      selection_criterion = "R2",
      return_full = FALSE,
      verbose = 1,  # Will generate message but shouldn't affect output type
      show_model_selection = FALSE
    )
  )

  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
})


test_that("get_B smart selection handles zero BETA values correctly", {
  set.seed(1008)

  # Create data with some zero BETAs
  training_MAF <- runif(50, 0.01, 0.4)
  training_BETA <- -log(training_MAF) / 10
  training_BETA[c(5, 15, 25)] <- 0  # Add some zeros
  training_N <- rep(1000, 50)
  training_P <- runif(50, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.5, 10)

  # Should work despite zero BETAs
  result <- suppressMessages(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "both",
      selection_criterion = "R2",
      return_full = TRUE,
      verbose = 2,
        show_model_selection = FALSE,
    )
  )

  # Should have selected a method (nested under $model)
  expect_true(!is.null(result$model$comparison$method_selected))
  expect_true(result$model$comparison$method_selected %in%
              c("beta_method", "pvalue_method"))

  # Primary B should be valid
  expect_type(result$B, "double")
  expect_true(all(is.finite(result$B)))
})


test_that("get_B comparison object has all expected fields when method='both'", {
  set.seed(1009)

  training_MAF <- runif(40, 0.01, 0.4)
  training_BETA <- sqrt(2 * training_MAF)
  training_N <- rep(1000, 40)
  training_P <- runif(40, 0.001, 0.1)

  target_MAF <- runif(10, 0.01, 0.4)
  target_case_prop <- rep(0.5, 10)

  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "both",
    selection_criterion = "CV_R2",
    return_full = TRUE,
    verbose = 0,
      show_model_selection = FALSE,
  )

  # Check training-time comparison fields (nested under $model$comparison)
  comp <- result$model$comparison
  expect_true(all(c("selection_criterion", "criterion_beta_method",
                    "criterion_pvalue_method", "method_selected") %in% names(comp)))

  # Check prediction-based comparison fields (nested under $model$comparison$prediction)
  pred_comp <- comp$prediction
  expect_true(all(c("correlation", "rmse", "mean_percent_diff",
                    "max_percent_diff", "summary_stats") %in% names(pred_comp)))

  # Check types
  expect_type(pred_comp$correlation, "double")
  expect_type(comp$criterion_beta_method, "double")
  expect_type(comp$criterion_pvalue_method, "double")
  expect_type(comp$method_selected, "character")
  expect_equal(comp$selection_criterion, "CV_R2")
})
