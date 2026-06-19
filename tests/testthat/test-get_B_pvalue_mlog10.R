########## Test get_B with P_mlog10 column support ##########
#
# This file tests the updated functionality to support P_mlog10 columns
# in training data for extreme p-values with better numerical stability.

test_that("get_B handles P_mlog10 column correctly", {

  # Generate test data with extreme p-values
  set.seed(123)
  n_snps <- 100

  # Create training data with both P and P_mlog10 columns
  training_data <- data.frame(
    MAF = runif(n_snps, 0.01, 0.5),
    BETA = rnorm(n_snps, 0, 0.1),
    P_mlog10 = runif(n_snps, 5, 350),  # -log10(p) from 1e-5 to 1e-350
    N = rep(10000, n_snps)
  )

  # Add regular P column for moderate values (where possible)
  moderate_idx <- training_data$P_mlog10 <= 300
  training_data$P <- NA
  training_data$P[moderate_idx] <- 10^(-training_data$P_mlog10[moderate_idx])

  # Test 1: get_B should use P_mlog10 when available
  target_mafs <- seq(0.05, 0.45, by = 0.1)

  result <- get_B(
    training_trait = "continuous",
    training_MAF = training_data$MAF,
    training_BETA = training_data$BETA,
    training_P_mlog10 = training_data$P_mlog10,
    training_N = training_data$N,
    target_trait = "continuous",
    target_MAF = target_mafs,
    target_SE = rep(0.05, length(target_mafs)),
    method = "pvalue",
    verbose = 0
  )

  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_true(all(is.finite(result)))

  # Test 2: Compare results using P vs P_mlog10 for moderate values
  training_moderate <- training_data[moderate_idx, ]

  result_p <- get_B(
    training_trait = "continuous",
    training_MAF = training_moderate$MAF,
    training_BETA = training_moderate$BETA,
    training_P = training_moderate$P,
    training_N = training_moderate$N,
    target_trait = "continuous",
    target_MAF = target_mafs,
    target_SE = rep(0.05, length(target_mafs)),
    method = "pvalue",
    verbose = 0
  )

  result_p_mlog10 <- get_B(
    training_trait = "continuous",
    training_MAF = training_moderate$MAF,
    training_BETA = training_moderate$BETA,
    training_P_mlog10 = training_moderate$P_mlog10,
    training_N = training_moderate$N,
    target_trait = "continuous",
    target_MAF = target_mafs,
    target_SE = rep(0.05, length(target_mafs)),
    method = "pvalue",
    verbose = 0
  )

  # Results should be nearly identical for moderate p-values
  expect_equal(result_p, result_p_mlog10, tolerance = 1e-10)

  # Test 3: Verify error when both P and P_mlog10 are provided
  expect_error(
    get_B(
      training_trait = "continuous",
      training_MAF = training_data$MAF,
      training_BETA = training_data$BETA,
      training_P = training_data$P[moderate_idx],
      training_P_mlog10 = training_data$P_mlog10[moderate_idx],
      training_N = training_data$N[moderate_idx],
      target_trait = "continuous",
      target_MAF = target_mafs,
      target_SE = rep(0.05, length(target_mafs)),
      method = "pvalue",
      verbose = 0
    ),
    "Provide either training_P or training_P_mlog10, not both"
  )
})

test_that("get_B_auto wrapper handles P_mlog10 column", {

  set.seed(456)
  n_snps <- 50

  # Create training data with P_mlog10 column
  training_df <- data.frame(
    MAF = runif(n_snps, 0.01, 0.5),
    P_mlog10 = runif(n_snps, 10, 100),  # -log10(p) from 1e-10 to 1e-100
    N = rep(5000, n_snps)
  )

  target_mafs <- seq(0.05, 0.45, by = 0.1)

  # Test that get_B_auto handles P_mlog10 column
  result <- get_B_auto(
    training_data = training_df,
    training_trait = "binary",
    target_MAF = target_mafs,
    target_case_prop = 0.5,
    verbose = 0
  )

  expect_type(result, "double")
  expect_length(result, length(target_mafs))
  expect_true(all(is.finite(result)))
})

test_that("get_B_advanced wrapper handles P_mlog10 column", {

  set.seed(789)
  n_snps <- 50

  # Create training data with both BETA and P_mlog10 columns
  training_df <- data.frame(
    MAF = runif(n_snps, 0.01, 0.5),
    BETA = rnorm(n_snps, 0, 0.1),
    P_mlog10 = runif(n_snps, 5, 50),  # -log10(p) from 1e-5 to 1e-50
    N = rep(8000, n_snps)
  )

  target_mafs <- seq(0.05, 0.45, by = 0.1)

  # Test that get_B_advanced handles P_mlog10 column with both methods
  result <- get_B_advanced(
    training_data = training_df,
    training_trait = "continuous",
    target_MAF = target_mafs,
    target_SE = 0.05,
    method = "both",
    return_full = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_length(result$B, length(target_mafs))
  expect_true(all(is.finite(result$B)))

  # Check that both methods were used (nested under $model in new structure)
  expect_true("model" %in% names(result))
  expect_true("method_used" %in% names(result$model))
  expect_true("comparison" %in% names(result$model))

  # Should have comparison data when method="both" (nested under $model)
  if (!is.null(result$model$comparison)) {
    expect_true("method_selected" %in% names(result$model$comparison))
    expect_true("criterion_beta_method" %in% names(result$model$comparison))
    expect_true("criterion_pvalue_method" %in% names(result$model$comparison))
  }
})

test_that(".compute_chi2_from_pvalue handles extreme p-values correctly", {

  # Test with regular p-values
  p_regular <- c(0.05, 0.01, 1e-10, 1e-50)
  chi2_regular <- GLOWr:::.compute_chi2_from_pvalue(p = p_regular)

  expect_type(chi2_regular, "double")
  expect_length(chi2_regular, length(p_regular))
  expect_true(all(chi2_regular > 0))

  # Test with -log10 p-values
  p_mlog10 <- c(5, 10, 50, 100, 350)  # p-values from 1e-5 to 1e-350
  chi2_mlog10 <- GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = p_mlog10)

  expect_type(chi2_mlog10, "double")
  expect_length(chi2_mlog10, length(p_mlog10))
  expect_true(all(chi2_mlog10 > 0))

  # Test asymptotic approximation for extreme values
  # For p_mlog10 > 300, chi2 ≈ 4.60517 * p_mlog10
  extreme_p_mlog10 <- c(350, 400, 500)
  chi2_extreme <- GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = extreme_p_mlog10)
  expected_chi2 <- 2 * log(10) * extreme_p_mlog10  # ≈ 4.60517 * p_mlog10

  # Should be very close to the approximation
  relative_diff <- abs(chi2_extreme - expected_chi2) / expected_chi2
  expect_true(all(relative_diff < 1e-10))

  # Test that results match for moderate values where both methods work
  p_moderate <- 1e-50
  p_mlog10_moderate <- 50

  chi2_from_p <- GLOWr:::.compute_chi2_from_pvalue(p = p_moderate)
  chi2_from_mlog10 <- GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = p_mlog10_moderate)

  expect_equal(chi2_from_p, chi2_from_mlog10, tolerance = 1e-10)

  # Test error handling
  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(),
    "Either p or p_mlog10 must be provided"
  )

  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(p = 0.05, p_mlog10 = 5),
    "Only one of p or p_mlog10 should be provided"
  )

  expect_error(
    GLOWr:::.compute_chi2_from_pvalue(p_mlog10 = -5),
    "p_mlog10 values must be non-negative"
  )
})
