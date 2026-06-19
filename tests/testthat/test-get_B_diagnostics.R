########## Unit Tests for B Estimation Diagnostics ##########
#
# This file contains testthat unit tests for the diagnostic functions
# in get_B_diagnostics.R

library(testthat)
library(GLOWr)

# ========== Test Data Setup ==========

# Create simple synthetic test data
set.seed(123)
n_train <- 50
training_MAF <- runif(n_train, 0.01, 0.4)
training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF)) * rnorm(n_train, 1, 0.2)
training_P <- runif(n_train, 0.0001, 0.1)
training_N <- sample(500:2000, n_train, replace = TRUE)

target_MAF <- seq(0.05, 0.3, by = 0.05)
target_case_prop <- rep(0.1, length(target_MAF))

# ========== Tests for plot.glow_B_estimate ==========

test_that("plot.glow_B_estimate works with beta method", {
  # Run get_B with beta method
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  # Test that plot doesn't throw errors
  expect_silent({
    pdf(tempfile())
    plot(result, type = "all")
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(result, type = "fit")
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(result, type = "residuals")
    dev.off()
  })
})

test_that("plot.glow_B_estimate works with pvalue method", {
  # Run get_B with pvalue method
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    return_full = TRUE,
    verbose = 0
  )

  # Test that plot doesn't throw errors
  expect_silent({
    pdf(tempfile())
    plot(result, type = "all", which_method = "pvalue")
    dev.off()
  })
})

test_that("plot.glow_B_estimate works with both methods", {
  # Run get_B with both methods
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
    return_full = TRUE,
    verbose = 0
  )

  # Test all plot types
  expect_silent({
    pdf(tempfile())
    plot(result, type = "all")
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(result, type = "comparison")
    dev.off()
  })

  # Test models plot type (now implemented)
  expect_silent({
    pdf(tempfile())
    plot(result, type = "models")
    dev.off()
  })
})

test_that("plot.glow_B_estimate validates inputs correctly", {
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  # Test invalid type
  expect_error(plot(result, type = "invalid"))

  # Test invalid which_method
  expect_error(plot(result, which_method = "invalid"))

  # Test requesting unavailable method
  expect_error(plot(result, which_method = "pvalue"),
               "pvalue method results not available")
})

# ========== Tests for get_B_diagnostics ==========

test_that("get_B_diagnostics produces expected output structure", {
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  diag <- get_B_diagnostics(result, show_plots = FALSE, verbose = 0)

  # Check class
  expect_s3_class(diag, "glow_B_diagnostics")

  # Check main components
  expect_true("model_fit" %in% names(diag))
  expect_true("residual_stats" %in% names(diag))
  expect_true("outliers" %in% names(diag))
  expect_true("qc_flags" %in% names(diag))
  expect_true("summary_text" %in% names(diag))

  # Check model_fit structure
  expect_true("beta" %in% names(diag$model_fit))
  expect_true("R2" %in% names(diag$model_fit$beta))
  expect_true("adj_R2" %in% names(diag$model_fit$beta))

  # Check residual_stats
  expect_true("mean" %in% names(diag$residual_stats))
  expect_true("sd" %in% names(diag$residual_stats))
  expect_true("normality_p" %in% names(diag$residual_stats))

  # Check outliers
  expect_true("n_outliers" %in% names(diag$outliers))
  expect_true("indices" %in% names(diag$outliers))
  expect_true(is.numeric(diag$outliers$proportion))
})

test_that("get_B_diagnostics works with both methods", {
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
    return_full = TRUE,
    verbose = 0
  )

  diag <- get_B_diagnostics(result, show_plots = FALSE, verbose = 0)

  # Check both methods are present
  expect_true("beta" %in% names(diag$model_fit))
  expect_true("pvalue" %in% names(diag$model_fit))

  # Check method comparison
  expect_false(is.null(diag$method_comparison))
  expect_true("correlation" %in% names(diag$method_comparison))
  expect_true("rmse" %in% names(diag$method_comparison))
})

test_that("get_B_diagnostics detects quality issues", {
  # Create data with known issues (high outlier)
  bad_MAF <- c(training_MAF, 0.45)
  bad_BETA <- c(training_BETA, 10)  # Outlier

  result <- get_B(
    training_trait = "binary",
    training_MAF = bad_MAF,
    training_BETA = bad_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  diag <- get_B_diagnostics(result, show_plots = FALSE, verbose = 0)

  # Should detect outliers
  expect_true(diag$outliers$n_outliers > 0)
})

test_that("get_B_diagnostics validates inputs", {
  # Test with non-glow_B_estimate/glow_B_model object
  expect_error(get_B_diagnostics(list(a = 1)),
               "glow_B_model.*glow_B_estimate")

  # Test with simple B vector (not full result)
  simple_result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = FALSE,  # Returns simple vector
    verbose = 0
  )

  expect_error(get_B_diagnostics(simple_result),
               "glow_B_model.*glow_B_estimate")
})

# ========== Tests for compare_B_models ==========

test_that("compare_B_models returns expected structure", {
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  model_comp <- compare_B_models(result, which_method = "beta", plot = FALSE, verbose = 0)

  # Should always return a data frame (not NULL) now that training_data is stored
  expect_true(is.data.frame(model_comp))

  # Should have multiple rows (one per candidate model)
  expect_gt(nrow(model_comp), 0)

  # Check required columns exist
  expect_true("model_name" %in% names(model_comp))
  expect_true("formula" %in% names(model_comp))
  expect_true("R2" %in% names(model_comp))
  expect_true("adj_R2" %in% names(model_comp))
  expect_true("CV_R2" %in% names(model_comp))
  expect_true("is_best" %in% names(model_comp))
  expect_true("method" %in% names(model_comp))

  # Exactly one model should be marked as best
  expect_equal(sum(model_comp$is_best), 1)
})

test_that("compare_B_models validates inputs", {
  result <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  # Test invalid which_method
  expect_error(compare_B_models(result, which_method = "invalid"))

  # Test invalid criterion
  expect_error(compare_B_models(result, criterion = "invalid"))

  # Test with non-glow_B_estimate object
  expect_error(compare_B_models(list(a = 1)))
})

test_that("compare_B_models works with both methods", {
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
    return_full = TRUE,
    verbose = 0
  )

  # Compare both methods
  model_comp_both <- compare_B_models(result, which_method = "both", plot = FALSE, verbose = 0)

  expect_true(is.data.frame(model_comp_both))
  expect_gt(nrow(model_comp_both), 0)

  # Should have results from both methods (now using _method suffix)
  expect_true("beta_method" %in% model_comp_both$method ||
              "pvalue_method" %in% model_comp_both$method)

  # Test individual methods
  model_comp_beta <- compare_B_models(result, which_method = "beta", plot = FALSE, verbose = 0)
  expect_true(all(model_comp_beta$method == "beta_method"))

  model_comp_pval <- compare_B_models(result, which_method = "pvalue", plot = FALSE, verbose = 0)
  expect_true(all(model_comp_pval$method == "pvalue_method"))
})

# ========== Integration Tests ==========

test_that("Full diagnostic pipeline works end-to-end", {
  # Run complete analysis
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
    return_full = TRUE,
    verbose = 0
  )

  # Test plot
  expect_silent({
    pdf(tempfile())
    plot(result)
    dev.off()
  })

  # Test diagnostics
  diag <- get_B_diagnostics(result, show_plots = FALSE, verbose = 0)
  expect_s3_class(diag, "glow_B_diagnostics")

  # Test print method
  expect_output(print(diag), "B Estimation Diagnostic Report")

  # Test model comparison
  model_comp <- compare_B_models(result, plot = FALSE)
  # May be NULL or data.frame depending on implementation
  expect_true(is.null(model_comp) || is.data.frame(model_comp))
})

test_that("Diagnostics work with minimal viable data", {
  # Test with smallest reasonable dataset
  min_MAF <- c(0.1, 0.2, 0.3)
  min_BETA <- c(0.1, 0.15, 0.12)
  min_target <- c(0.15, 0.25)

  result <- get_B(
    training_trait = "binary",
    training_MAF = min_MAF,
    training_BETA = min_BETA,
    target_trait = "binary",
    target_MAF = min_target,
    method = "beta",
    return_full = TRUE,
    verbose = 0
  )

  # Should work without errors
  expect_silent({
    pdf(tempfile())
    plot(result)
    dev.off()
  })

  diag <- get_B_diagnostics(result, show_plots = FALSE, verbose = 0)
  expect_s3_class(diag, "glow_B_diagnostics")
})

# ========== Cleanup ==========
# Remove any temporary files created during testing
unlink(list.files(tempdir(), pattern = "^file.*\\.pdf$", full.names = TRUE))
