########## Tests for User Interface Wrappers ##########
#
# This file contains comprehensive tests for get_B_auto() and get_B_advanced()
# user interface wrappers.

# Load testthat
library(testthat)
library(GLOWr)

context("User Interface Wrappers for get_B")

# ========== Helper Function to Create Test Data ==========

.create_test_training_data <- function(n = 50, trait_type = "binary", include_beta = TRUE, include_pvalue = TRUE) {
  set.seed(12345)

  MAF <- runif(n, 0.01, 0.4)
  P <- runif(n, 1e-6, 0.1)
  N <- rep(1000, n)

  # Generate BETA based on MAF
  if (include_beta) {
    BETA <- sqrt(MAF * (1 - MAF)) * rnorm(n, 0.1, 0.02)
  } else {
    BETA <- NULL
  }

  # Create data.frame
  df <- data.frame(MAF = MAF)

  if (include_beta) {
    df$BETA <- BETA
  }

  if (include_pvalue) {
    df$P <- P
    df$N <- N
  }

  df$rsID <- paste0("rs", seq_len(n))

  return(df)
}


# ========== Tests for get_B_auto() ==========

test_that("get_B_auto works with data.frame input (BETA data)", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run get_B_auto
  result <- get_B_auto(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    verbose = 0
  )

  # Check output
  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_true(all(is.finite(result)))
  expect_true(all(result > 0))  # Effect sizes should be positive
})


test_that("get_B_auto works with data.frame input (P-value data)", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = FALSE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run get_B_auto - binary trait
  result <- get_B_auto(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    verbose = 0
  )

  # Check output
  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_true(all(is.finite(result)))
})


test_that("get_B_auto works with glow_training_data object", {

  # Create and prepare test data
  training_df <- .create_test_training_data(n = 40, include_beta = TRUE, include_pvalue = TRUE)

  prepared_data <- prepare_B_training_data(
    data = training_df,
    trait_type = "binary",
    verbose = 0
  )

  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run get_B_auto
  result <- get_B_auto(
    training_data = prepared_data,
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    verbose = 0
  )

  # Check output
  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_true(all(is.finite(result)))
})


test_that("get_B_auto auto-detects trait type from glow_training_data", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = TRUE)

  prepared_data <- prepare_B_training_data(
    data = training_df,
    trait_type = "continuous",
    verbose = 0
  )

  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run without specifying training_trait (should auto-detect)
  result <- get_B_auto(
    training_data = prepared_data,
    target_MAF = target_mafs,
    target_SE = 0.1,
    verbose = 0
  )

  expect_type(result, "double")
  expect_length(result, length(target_mafs))
})


test_that("get_B_auto errors when trait type not specified for data.frame", {

  # Create test data
  training_df <- .create_test_training_data(n = 30)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Should error without training_trait
  expect_error(
    get_B_auto(
      training_data = training_df,
      target_MAF = target_mafs,
      verbose = 0
    ),
    "training_trait must be specified"
  )
})


test_that("get_B_auto errors when required columns missing", {

  # Create data without MAF
  training_df <- data.frame(
    BETA = rnorm(30),
    P = runif(30, 0.001, 0.1),
    N = rep(1000, 30)
  )

  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Should error without MAF column
  expect_error(
    get_B_auto(
      training_data = training_df,
      training_trait = "binary",
      target_MAF = target_mafs,
      verbose = 0
    ),
    "MAF"
  )
})


test_that("get_B_auto errors when neither BETA nor P+N available", {

  # Create data with only MAF
  training_df <- data.frame(MAF = runif(30, 0.01, 0.4))
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Should error without estimation data
  expect_error(
    get_B_auto(
      training_data = training_df,
      training_trait = "binary",
      target_MAF = target_mafs,
      verbose = 0
    ),
    "BETA|P and N"
  )
})


test_that("get_B_auto handles cross-trait estimation", {

  # Create continuous trait training data
  training_df <- .create_test_training_data(n = 30, include_beta = FALSE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Estimate for binary target trait
  result <- get_B_auto(
    training_data = training_df,
    training_trait = "continuous",
    target_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    verbose = 0
  )

  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_true(all(is.finite(result)))
})


test_that("get_B_auto handles scalar target_SE and target_case_prop", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = FALSE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Test with scalar target_case_prop
  result1 <- get_B_auto(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.3,  # Single value
    verbose = 0
  )

  expect_type(result1, "double")
  expect_length(result1, length(target_mafs))

  # Test with scalar target_SE
  result2 <- get_B_auto(
    training_data = training_df,
    training_trait = "continuous",
    target_MAF = target_mafs,
    target_SE = 0.05,  # Single value
    verbose = 0
  )

  expect_type(result2, "double")
  expect_length(result2, length(target_mafs))
})


# ========== Tests for get_B_advanced() ==========

test_that("get_B_advanced works with all features", {

  # Create test data
  training_df <- .create_test_training_data(n = 40, include_beta = TRUE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with full features
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    method = "both",
    selection_criterion = "R2",
    outlier_method = "both",
    outlier_action = "flag",
    return_full = TRUE,
    verbose = 0
  )

  # Check output structure (new nested glow_B_estimate)
  expect_s3_class(result, "glow_B_estimate")
  expect_type(result$B, "double")
  expect_length(result$B, length(target_mafs))
  expect_equal(result$model$method_used, "both")
  expect_true(!is.null(result$B_beta_method))
  expect_true(!is.null(result$B_pvalue_method))
  expect_true(!is.null(result$model$comparison))
})


test_that("get_B_advanced cross-validation works", {

  # Create test data
  training_df <- .create_test_training_data(n = 50, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with cross-validation
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    method = "beta",
    cross_validation = TRUE,  # Uses LOOCV
    return_full = TRUE,
    verbose = 0
  )

  # Check that CV criterion was used (LOOCV), nested under $model
  # Note: "CV" is deprecated and converted to "CV_R2" internally
  expect_equal(result$model$selection_criterion, "CV_R2")
  expect_type(result$B, "double")
  expect_length(result$B, length(target_mafs))
})


test_that("get_B_advanced show_diagnostics runs without error", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with show_diagnostics = TRUE (should not error)
  # Note: show_diagnostics produces plot output, so we don't use expect_silent()
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    method = "both",
    show_diagnostics = TRUE,  # This triggers plot()
    return_full = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_estimate")
})


test_that("get_B_advanced show_outlier_details runs without error", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Capture output to avoid cluttering test results
  output <- capture.output({
    result <- get_B_advanced(
      training_data = training_df,
      training_trait = "binary",
      target_MAF = target_mafs,
      method = "beta",  # Use beta method since we don't have P values
      outlier_method = "statistical",
      show_outlier_details = TRUE,
      return_full = TRUE,
      verbose = 0
    )
  })

  expect_s3_class(result, "glow_B_estimate")
  expect_true(length(output) > 0)  # Should have printed something
})


test_that("get_B_advanced with method='beta' only", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with beta method only
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  expect_equal(result$model$method_used, "beta_method")
  expect_true(!is.null(result$B_beta_method))
  expect_null(result$B_pvalue_method)
  # When only one method is used, no prediction-based comparison
  expect_null(result$model$comparison$prediction)
})


test_that("get_B_advanced with method='pvalue' only", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with pvalue method only
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    method = "pvalue",
    return_full = TRUE,
    verbose = 0
  )

  expect_equal(result$model$method_used, "pvalue_method")
  expect_null(result$B_beta_method)
  expect_true(!is.null(result$B_pvalue_method))
  # When only one method is used, no prediction-based comparison
  expect_null(result$model$comparison$prediction)
})


test_that("get_B_advanced errors when method incompatible with data", {

  # Create data with BETA only
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Should error when requesting pvalue method without P+N
  expect_error(
    get_B_advanced(
      training_data = training_df,
      training_trait = "binary",
      target_MAF = target_mafs,
      target_case_prop = 0.5,
      method = "pvalue",
      verbose = 0
    ),
    "requires P/P_mlog10 and N"
  )
})


test_that("get_B_advanced custom_models works", {

  # Create test data
  training_df <- .create_test_training_data(n = 50, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Define custom models
  custom_formulas <- list(
    formula(Y ~ poly(X, 2)),
    formula(Y ~ poly(X, 3))
  )

  # Run with custom models
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    method = "beta",
    custom_models = custom_formulas,
    return_full = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_type(result$B, "double")
})


test_that("get_B_advanced outlier_action='remove' works", {

  # Create test data with some potential outliers
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)

  # Add an obvious outlier
  training_df$MAF[1] <- 0.4
  training_df$BETA[1] <- 5.0  # Very large effect for common variant

  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with outlier removal
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    method = "beta",  # Use beta method since we don't have P values
    outlier_method = "biological",
    outlier_action = "remove",
    return_full = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_type(result$B, "double")

  # Check that some outliers were detected (nested under $model in new structure)
  # (may or may not be removed depending on biological rules)
  expect_true(!is.null(result$model$outliers))
})


test_that("get_B_advanced return_full=FALSE works", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run with return_full = FALSE
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    method = "beta",  # Use beta method since we don't have P values
    return_full = FALSE,
    verbose = 0
  )

  # Should return simple vector
  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_false(inherits(result, "glow_B_estimate"))
})


test_that("get_B_advanced with glow_training_data object", {

  # Create and prepare test data
  training_df <- .create_test_training_data(n = 40, include_beta = TRUE, include_pvalue = TRUE)

  prepared_data <- prepare_B_training_data(
    data = training_df,
    trait_type = "binary",
    verbose = 0
  )

  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run get_B_advanced
  result <- get_B_advanced(
    training_data = prepared_data,
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    method = "both",
    return_full = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_equal(result$model$method_used, "both")
})


# ========== Tests Comparing Direct get_B with Wrappers ==========

test_that("get_B_auto produces same results as get_B", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run get_B_auto
  result_auto <- get_B_auto(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    verbose = 0
  )

  # Run get_B directly with same settings
  result_direct <- get_B(
    training_trait = "binary",
    training_MAF = training_df$MAF,
    training_BETA = training_df$BETA,
    target_trait = "binary",
    target_MAF = target_mafs,
    method = "auto",
    selection_criterion = "R2",
    outlier_method = "both",
    outlier_action = "flag",
    return_full = FALSE,
    verbose = 0
  )

  # Results should be identical (within numerical tolerance)
  expect_equal(result_auto, result_direct, tolerance = 1e-10)
})


test_that("get_B_advanced produces same results as get_B", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = TRUE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # Run get_B_advanced
  result_advanced <- get_B_advanced(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    method = "both",
    selection_criterion = "CV_R2",
    return_full = TRUE,
    verbose = 0
  )

  # Run get_B directly with same settings
  result_direct <- get_B(
    training_trait = "binary",
    training_MAF = training_df$MAF,
    training_BETA = training_df$BETA,
    training_P = training_df$P,
    training_N = training_df$N,
    target_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    method = "both",
    selection_criterion = "CV_R2",
    outlier_method = "both",
    outlier_action = "flag",
    return_full = TRUE,
    verbose = 0
  )

  # B estimates should be identical
  expect_equal(result_advanced$B, result_direct$B, tolerance = 1e-10)
  expect_equal(result_advanced$B_beta_method, result_direct$B_beta_method, tolerance = 1e-10)
  expect_equal(result_advanced$B_pvalue_method, result_direct$B_pvalue_method, tolerance = 1e-10)
})


# ========== Edge Cases and Boundary Conditions ==========

test_that("get_B_auto works with minimal training data", {

  # Create minimal data (n=10)
  training_df <- .create_test_training_data(n = 10, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- c(0.1, 0.2, 0.3)

  # Should work but may have warnings
  expect_warning(
    result <- get_B_auto(
      training_data = training_df,
      training_trait = "binary",
      target_MAF = target_mafs,
      verbose = 0
    ),
    NA  # Expect no warning or any warning is OK
  )

  expect_type(result, "double")
  expect_length(result, 3)
})


test_that("get_B_auto handles single target MAF", {

  # Create test data
  training_df <- .create_test_training_data(n = 30, include_beta = TRUE, include_pvalue = FALSE)

  # Single target MAF
  result <- get_B_auto(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = 0.2,
    verbose = 0
  )

  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.finite(result))
})


test_that("User interfaces handle verbose levels correctly", {

  # Create test data
  training_df <- .create_test_training_data(n = 20, include_beta = TRUE, include_pvalue = FALSE)
  target_mafs <- seq(0.05, 0.3, by = 0.05)

  # verbose = 0 should not produce messages (warnings from CV are OK)
  output <- capture.output({
    result <- suppressWarnings(
      get_B_auto(
        training_data = training_df,
        training_trait = "binary",
        target_MAF = target_mafs,
        verbose = 0
      )
    )
  }, type = "message")

  expect_equal(length(output), 0)  # No messages
  expect_type(result, "double")

  # verbose >= 1 should produce messages
  # Note: capture.output(type = "message") is unreliable inside testthat;
  # use expect_message() instead
  expect_message(
    suppressWarnings(
      get_B_auto(
        training_data = training_df,
        training_trait = "binary",
        target_MAF = target_mafs,
        verbose = 1
      )
    ),
    "Running B estimation"
  )
})


# ========== Integration Tests ==========

test_that("Full workflow: prepare_B_training_data -> get_B_auto", {

  # Create raw data
  raw_data <- data.frame(
    snp = paste0("rs", 1:50),
    maf = runif(50, 0.01, 0.4),
    effect = rnorm(50, 0, 0.1),
    pvalue = runif(50, 1e-8, 0.1),
    sample_size = rep(5000, 50)
  )

  # Prepare data
  prepared <- prepare_B_training_data(
    data = raw_data,
    column_mapping = list(
      rsID = "snp",
      MAF = "maf",
      BETA = "effect",
      P = "pvalue",
      N = "sample_size"
    ),
    trait_type = "continuous",
    verbose = 0
  )

  # Use in get_B_auto
  result <- get_B_auto(
    training_data = prepared,
    target_MAF = seq(0.05, 0.3, by = 0.05),
    target_SE = 0.1,
    verbose = 0
  )

  expect_type(result, "double")
  expect_true(all(is.finite(result)))
})


test_that("Full workflow: prepare_B_training_data -> get_B_advanced", {

  # Create raw data
  raw_data <- data.frame(
    MAF = runif(50, 0.01, 0.4),
    BETA = rnorm(50, 0, 0.1),
    P = runif(50, 1e-8, 0.1),
    N = rep(5000, 50)
  )

  # Prepare data
  prepared <- prepare_B_training_data(
    data = raw_data,
    trait_type = "binary",
    verbose = 0
  )

  # Use in get_B_advanced
  result <- get_B_advanced(
    training_data = prepared,
    target_MAF = seq(0.05, 0.3, by = 0.05),
    target_case_prop = 0.5,
    method = "both",
    cross_validation = TRUE,  # Uses LOOCV
    return_full = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_equal(result$model$method_used, "both")
  expect_equal(result$model$selection_criterion, "CV_R2")  # LOOCV (CV -> CV_R2)
})
