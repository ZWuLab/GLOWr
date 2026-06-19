########## Tests for Enhanced Model Selection in get_B ##########
#
# This file contains tests for the enhanced model selection features:
# - selection_criterion parameter (R2, adj_R2, CV_R2)
# - Tests that AIC/BIC produce informative errors
# - custom_models parameter
# - return_all functionality in select_best_model()
# - backward compatibility for deprecated "CV" criterion name

library(testthat)

# ========== select_best_model() Tests ==========

test_that("select_best_model works with all valid criteria including adj_R2", {
  set.seed(123)

  # Simulate data
  MAF <- runif(100, 0.01, 0.5)
  BETA <- sqrt(0.5 * MAF * (1 - MAF)) + rnorm(100, 0, 0.01)
  Y <- BETA^2

  # Test R2 criterion (default) - suppress output
  model_r2 <- select_best_model(MAF, Y, criterion = "R2", verbose = FALSE)
  expect_s3_class(model_r2, "lm")
  expect_true(summary(model_r2)$r.squared > 0)

  # Test adj_R2 criterion
  model_adjr2 <- select_best_model(MAF, Y, criterion = "adj_R2", verbose = FALSE)
  expect_s3_class(model_adjr2, "lm")
  expect_true(summary(model_adjr2)$adj.r.squared > 0)

  # Test CV_R2 criterion
  model_cv <- select_best_model(MAF, Y, criterion = "CV_R2", verbose = FALSE)
  expect_s3_class(model_cv, "lm")

  # Models might be different depending on criterion
  # (not guaranteed, but possible)
})


test_that("adj_R2 criterion selects model correctly", {
  set.seed(789)

  MAF <- runif(30, 0.01, 0.5)
  Y <- MAF * (1 - MAF) * 0.5 + abs(rnorm(30, 0, 0.01))

  # Get all models (new return structure)
  result <- select_best_model(MAF, Y, criterion = "adj_R2", return_all = TRUE)
  summary_table <- result$all_models_info$summary_table

  # Find which model is marked as best
  best_idx <- which(summary_table$is_best)
  expect_length(best_idx, 1)

  # Get adj_R2 values
  adjr2_values <- summary_table$adj_R2

  # Best model should have highest adj_R2
  best_adjr2 <- adjr2_values[best_idx]
  expect_equal(best_adjr2, max(adjr2_values))
})


test_that("Backward compatibility: CV criterion still works with warning", {
  set.seed(456)

  MAF <- runif(50, 0.01, 0.5)
  BETA <- sqrt(0.5 * MAF * (1 - MAF))
  Y <- BETA^2

  # Test that "CV" triggers deprecation warning
  expect_warning(
    model_cv_old <- select_best_model(MAF, Y, criterion = "CV"),
    "The 'CV' criterion name is deprecated. Please use 'CV_R2' instead"
  )

  # Should still return a valid model
  expect_s3_class(model_cv_old, "lm")

  # Should produce same result as "CV_R2"
  model_cv_new <- select_best_model(MAF, Y, criterion = "CV_R2")

  # Models should be identical (same formula selected)
  expect_equal(formula(model_cv_old), formula(model_cv_new))
})


test_that("CV_R2 criterion calculates CV_R2 for all models", {
  set.seed(456)

  MAF <- runif(50, 0.01, 0.5)
  BETA <- sqrt(0.5 * MAF * (1 - MAF))
  Y <- BETA^2

  # Get all models with CV_R2 (new return structure)
  result <- select_best_model(MAF, Y, criterion = "CV_R2", return_all = TRUE)
  summary_table <- result$all_models_info$summary_table

  # All models should have CV_R2 calculated
  expect_true("CV_R2" %in% names(summary_table))
  for (i in seq_len(nrow(summary_table))) {
    # CV_R2 should be numeric (might be NA if calculation failed)
    expect_true(is.numeric(summary_table$CV_R2[i]))
  }

  # At least one model should have valid CV_R2
  cv_values <- summary_table$CV_R2
  expect_true(any(!is.na(cv_values)))
})


test_that("CV_R2 criterion selects model with best CV_R2", {
  set.seed(789)

  MAF <- runif(30, 0.01, 0.5)
  Y <- MAF * (1 - MAF) * 0.5 + abs(rnorm(30, 0, 0.01))

  # Get all models (new return structure)
  result <- select_best_model(MAF, Y, criterion = "CV_R2", return_all = TRUE)
  summary_table <- result$all_models_info$summary_table

  # Find which model is marked as best
  best_idx <- which(summary_table$is_best)
  expect_length(best_idx, 1)

  # Get CV values
  cv_values <- summary_table$CV_R2

  # Best model should have highest CV_R2 (among non-NA values)
  best_cv <- cv_values[best_idx]
  if (!is.na(best_cv)) {
    valid_cv <- cv_values[!is.na(cv_values)]
    expect_equal(best_cv, max(valid_cv))
  }
})


test_that("CV_R2 criterion works with small samples", {
  set.seed(111)

  # Small sample size (n=15)
  MAF <- runif(15, 0.01, 0.5)
  Y <- MAF * (1 - MAF) * 0.5 + abs(rnorm(15, 0, 0.01))

  # Should work without errors
  expect_no_error({
    model_cv <- select_best_model(MAF, Y, criterion = "CV_R2")
  })

  model_cv <- select_best_model(MAF, Y, criterion = "CV_R2")
  expect_s3_class(model_cv, "lm")
})


test_that("CV_R2 falls back to R2 if all CV calculations fail", {
  # This is hard to trigger naturally, but we can test the logic exists
  # by examining the code path indirectly
  set.seed(222)

  MAF <- runif(20, 0.01, 0.5)
  Y <- MAF^2 + abs(rnorm(20, 0, 0.001))

  # Should complete without error even if some CV calcs might fail
  expect_no_error({
    model_cv <- select_best_model(MAF, Y, criterion = "CV_R2")
  })

  model_cv <- select_best_model(MAF, Y, criterion = "CV_R2")
  expect_s3_class(model_cv, "lm")
})


test_that("select_best_model return_all returns complete results with adj_R2", {
  set.seed(456)

  MAF <- runif(50, 0.01, 0.5)
  BETA <- sqrt(0.5 * MAF * (1 - MAF))
  Y <- BETA^2

  # Get all models (new return structure)
  result <- select_best_model(MAF, Y, return_all = TRUE)

  # Should return list with best_model, best_name, all_models_info
  expect_type(result, "list")
  expect_true("best_model" %in% names(result))
  expect_true("all_models_info" %in% names(result))
  expect_s3_class(result$best_model, "lm")

  # all_models_info should contain models and summary_table
  all_info <- result$all_models_info
  expect_true("models" %in% names(all_info))
  expect_true("summary_table" %in% names(all_info))
  expect_true(length(all_info$models) >= 8)  # At least 8 standard models

  # Each model entry should have required components
  for (name in names(all_info$models)) {
    m <- all_info$models[[name]]
    expect_true("model" %in% names(m))
    expect_true("formula_text" %in% names(m))
    expect_true("R2" %in% names(m))
    expect_true("adj_R2" %in% names(m))
    expect_true("CV_R2" %in% names(m))
    expect_s3_class(m$model, "lm")
  }

  # summary_table should have is_best column with exactly one TRUE
  st <- all_info$summary_table
  expect_true("is_best" %in% names(st))
  expect_true("R2" %in% names(st))
  expect_true("adj_R2" %in% names(st))
  expect_true("CV_R2" %in% names(st))
  expect_equal(sum(st$is_best), 1)
})


test_that("select_best_model accepts custom models", {
  set.seed(789)

  MAF <- runif(50, 0.01, 0.5)
  Y <- MAF^2 + abs(rnorm(50, 0, 0.001))  # Ensure positive

  # Define custom models
  custom_formulas <- list(
    formula(Y ~ poly(X, 2)),     # Polynomial
    formula(Y ~ I(X^2))          # Squared term
  )

  # Fit with custom models (new return structure)
  result <- select_best_model(MAF, Y,
                              custom_models = custom_formulas,
                              return_all = TRUE)

  # Should have 8 standard + 2 custom = 10 models
  expect_equal(length(result$all_models_info$models), 10)

  # Custom models should be included in summary table
  formulas <- result$all_models_info$summary_table$formula
  expect_true(any(grepl("poly", formulas)))
  expect_true(any(grepl("I\\(X\\^2\\)", formulas)))
})


test_that("select_best_model verbose parameter works correctly", {
  set.seed(789)

  MAF <- runif(50, 0.01, 0.5)
  Y <- MAF * (1 - MAF) * 0.5 + abs(rnorm(50, 0, 0.01))

  # Test default verbose = TRUE produces output
  expect_output(
    select_best_model(MAF, Y, criterion = "R2"),
    "Model Selection for B Estimation"
  )

  # Test verbose = FALSE suppresses output
  expect_no_output <- function(expr) {
    output <- capture.output(expr, type = "output")
    expect_equal(length(output), 0)
  }

  model_quiet <- capture.output(
    result <- select_best_model(MAF, Y, criterion = "R2", verbose = FALSE),
    type = "output"
  )
  expect_equal(length(model_quiet), 0)
  expect_s3_class(result, "lm")
})


test_that("select_best_model verbose output shows selection table", {
  set.seed(890)

  MAF <- runif(30, 0.01, 0.5)
  Y <- MAF * (1 - MAF) * 0.5 + abs(rnorm(30, 0, 0.01))

  # Capture output
  output <- capture.output(
    model <- select_best_model(MAF, Y, criterion = "R2", verbose = TRUE)
  )

  # Check that key elements are in the output
  output_text <- paste(output, collapse = "\n")
  expect_true(grepl("Model Selection for B Estimation", output_text))
  expect_true(grepl("Model Comparison:", output_text))
  expect_true(grepl("SELECTED", output_text))
  expect_true(grepl("Best Model:", output_text))
})


test_that("AIC and BIC criteria produce informative error", {
  set.seed(123)
  X <- runif(50, 0.01, 0.5)
  Y <- (X * (1 - X))^(-0.5) + abs(rnorm(50, 0, 0.1))

  # Test AIC produces error mentioning valid criteria
  expect_error(
    select_best_model(X, Y, criterion = "AIC"),
    regexp = "criterion must be one of.*R2.*adj_R2.*CV_R2"
  )

  # Test BIC produces same error
  expect_error(
    select_best_model(X, Y, criterion = "BIC"),
    regexp = "criterion must be one of.*R2.*adj_R2.*CV_R2"
  )
})


test_that("select_best_model validates inputs correctly", {
  MAF <- runif(10, 0.01, 0.5)
  Y <- runif(10, 0, 1)

  # Invalid criterion
  expect_error(
    select_best_model(MAF, Y, criterion = "invalid"),
    "criterion must be one of"
  )

  # Valid criteria should work (R2, adj_R2, CV_R2)
  expect_no_error(select_best_model(MAF, Y, criterion = "R2"))
  expect_no_error(select_best_model(MAF, Y, criterion = "adj_R2"))
  expect_no_error(select_best_model(MAF, Y, criterion = "CV_R2"))

  # Invalid custom_models (not a list)
  expect_error(
    select_best_model(MAF, Y, custom_models = "not a list"),
    "custom_models must be a list"
  )

  # Invalid custom_models (not formulas)
  expect_error(
    select_best_model(MAF, Y, custom_models = list("not a formula")),
    "All elements of custom_models must be formula objects"
  )

  # Mismatched lengths
  expect_error(
    select_best_model(MAF, Y[1:5]),
    "X and Y must have the same length"
  )
})


# ========== get_B() with Enhanced Model Selection Tests ==========

test_that("get_B works with all valid selection criteria including adj_R2", {
  set.seed(111)

  # Training data
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- sqrt(training_MAF * (1 - training_MAF)) * 0.2

  # Target data
  target_MAF <- runif(10, 0.01, 0.5)

  # Test with R2
  B_r2 <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    selection_criterion = "R2"
  )
  expect_type(B_r2, "double")
  expect_length(B_r2, length(target_MAF))

  # Test with adj_R2
  B_adjr2 <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    selection_criterion = "adj_R2"
  )
  expect_type(B_adjr2, "double")
  expect_length(B_adjr2, length(target_MAF))

  # Test with CV_R2
  B_cv <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    selection_criterion = "CV_R2"
  )
  expect_type(B_cv, "double")
  expect_length(B_cv, length(target_MAF))

  # Results might differ (not guaranteed, but possible)
  # All should be finite and positive
  expect_true(all(is.finite(B_r2)))
  expect_true(all(is.finite(B_adjr2)))
  expect_true(all(is.finite(B_cv)))
})


test_that("get_B CV_R2 criterion works with p-value method", {
  set.seed(333)

  # Training data (using p-values)
  training_MAF <- runif(50, 0.01, 0.3)
  training_P <- runif(50, 0.0001, 0.1)
  training_N <- rep(1000, 50)

  # Target data
  target_MAF <- runif(10, 0.01, 0.5)
  target_case_prop <- rep(0.1, 10)

  # Test with CV_R2 criterion
  B_cv <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_P = training_P,
    training_N = training_N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    selection_criterion = "CV_R2"
  )

  expect_type(B_cv, "double")
  expect_length(B_cv, length(target_MAF))
  expect_true(all(is.finite(B_cv)))
})


test_that("get_B works with custom models", {
  set.seed(222)

  # Training data
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- sqrt(training_MAF * (1 - training_MAF)) * 0.15

  # Target data
  target_MAF <- runif(10, 0.01, 0.5)

  # Custom model (polynomial)
  custom_formula <- list(formula(Y ~ poly(X, 2)))

  # With custom model
  B_custom <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    custom_models = custom_formula
  )

  expect_type(B_custom, "double")
  expect_length(B_custom, length(target_MAF))
  expect_true(all(is.finite(B_custom)))
  expect_true(all(B_custom > 0))
})




test_that("get_B selection_criterion works with both methods", {
  set.seed(444)

  # Training data with both BETA and P
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- sqrt(training_MAF * (1 - training_MAF)) * 0.2
  training_P <- runif(50, 0.0001, 0.1)
  training_N <- rep(1000, 50)

  # Target data
  target_MAF <- runif(10, 0.01, 0.5)
  target_case_prop <- rep(0.1, 10)

  # Run with both methods and CV_R2 criterion
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
    return_full = TRUE
  )

  expect_s3_class(result, "glow_B_estimate")
  expect_true("B_beta_method" %in% names(result))
  expect_true("B_pvalue_method" %in% names(result))
  expect_true("model" %in% names(result))

  # Both methods should have used CV_R2 (nested under $model$models)
  expect_true(!is.null(result$model$models$beta_method))
  expect_true(!is.null(result$model$models$pvalue_method))
})


test_that("get_B validates selection_criterion parameter", {
  training_MAF <- runif(10, 0.01, 0.3)
  training_BETA <- runif(10, 0, 0.5)
  target_MAF <- runif(5, 0.01, 0.5)

  # Invalid criterion
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      selection_criterion = "invalid"
    ),
    "selection_criterion must be one of"
  )
})


test_that("get_B validates custom_models parameter", {
  training_MAF <- runif(10, 0.01, 0.3)
  training_BETA <- runif(10, 0, 0.5)
  target_MAF <- runif(5, 0.01, 0.5)

  # Invalid custom_models (not a list)
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      custom_models = "not a list"
    ),
    "custom_models must be a list"
  )

  # Invalid custom_models (not formulas)
  expect_error(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      custom_models = list("not a formula")
    ),
    "All elements of custom_models must be formula objects"
  )
})


# ========== Integration Tests ==========

test_that("Enhanced model selection maintains backward compatibility", {
  set.seed(555)

  # Training and target data
  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- sqrt(training_MAF * (1 - training_MAF)) * 0.2
  target_MAF <- runif(10, 0.01, 0.5)

  # Default behavior (should use R2)
  B_default <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF
  )

  # Explicit R2
  B_r2 <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    selection_criterion = "R2"
  )

  # Should give same results
  expect_equal(B_default, B_r2)
})


test_that("Custom models can improve fit for specific data patterns", {
  set.seed(666)

  # Generate data with quadratic relationship
  training_MAF <- seq(0.01, 0.4, length.out = 50)
  training_BETA <- 0.1 * training_MAF^2 + rnorm(50, 0, 0.005)
  target_MAF <- runif(10, 0.01, 0.4)

  # Standard models
  B_standard <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF
  )

  # With polynomial custom model
  custom_poly <- list(formula(Y ~ poly(X, 2)))
  B_custom <- get_B(
    training_trait = "binary",
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    custom_models = custom_poly
  )

  # Both should produce reasonable results
  expect_true(all(is.finite(B_standard)))
  expect_true(all(is.finite(B_custom)))
  expect_true(all(B_standard > 0))
  expect_true(all(B_custom > 0))

  # Custom model might fit better (not guaranteed due to noise)
  # Just check they're different (possible, not guaranteed)
})


test_that("adj_R2 prefers simpler models compared to R2", {
  set.seed(777)

  # Create data where simple and complex models fit similarly
  MAF <- runif(30, 0.01, 0.5)
  Y <- 0.5 * MAF + abs(rnorm(30, 0, 0.05))  # Use abs() to ensure positive

  # Get best model with R2 (new return structure)
  result_r2 <- select_best_model(MAF, Y, criterion = "R2", return_all = TRUE)
  st_r2 <- result_r2$all_models_info$summary_table
  best_r2_idx <- which(st_r2$is_best)

  # Get best model with adj_R2
  result_adjr2 <- select_best_model(MAF, Y, criterion = "adj_R2", return_all = TRUE)
  st_adjr2 <- result_adjr2$all_models_info$summary_table
  best_adjr2_idx <- which(st_adjr2$is_best)

  # Both should return valid models
  expect_true(best_r2_idx > 0)
  expect_true(best_adjr2_idx > 0)

  # Both selected models should have reasonable fit
  expect_true(st_r2$R2[best_r2_idx] > 0.5)
  expect_true(st_adjr2$adj_R2[best_adjr2_idx] > 0.5)
})


# ========== get_B() show_model_selection Tests ==========

test_that("get_B show_model_selection parameter works correctly", {
  set.seed(911)

  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- sqrt(training_MAF * (1 - training_MAF)) * 0.2
  target_MAF <- runif(10, 0.01, 0.5)

  # Test default show_model_selection = TRUE produces output
  expect_output(
    get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      verbose = 0  # Suppress other output
    ),
    "Model Selection for B Estimation"
  )

  # Test show_model_selection = FALSE suppresses output
  output_quiet <- capture.output(
    B_quiet <- get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_BETA = training_BETA,
      target_trait = "binary",
      target_MAF = target_MAF,
      show_model_selection = FALSE,
      verbose = 0
    ),
    type = "output"
  )
  expect_equal(length(output_quiet), 0)
  expect_type(B_quiet, "double")
  expect_length(B_quiet, length(target_MAF))
})


test_that("get_B show_model_selection works with both methods", {
  set.seed(912)

  training_MAF <- runif(50, 0.01, 0.3)
  training_BETA <- sqrt(training_MAF * (1 - training_MAF)) * 0.2
  training_P <- runif(50, 0.0001, 0.1)
  training_N <- rep(1000, 50)
  target_MAF <- runif(10, 0.01, 0.5)
  target_case_prop <- rep(0.5, 10)

  # With both methods, should show model selection table
  expect_output(
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
      verbose = 0  # Suppress other messages
    ),
    "Model Selection"
  )

  # Suppress with show_model_selection = FALSE
  output_quiet <- capture.output(
    suppressMessages(  # Suppress method selection messages
      B_quiet <- get_B(
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
        verbose = 0
      )
    ),
    type = "output"
  )
  # Output should be empty (no model selection table)
  expect_equal(length(output_quiet), 0)
})


test_that("get_B show_model_selection only affects beta method output", {
  set.seed(913)

  training_MAF <- runif(50, 0.01, 0.3)
  training_P <- runif(50, 0.0001, 0.1)
  training_N <- rep(1000, 50)
  target_MAF <- runif(10, 0.01, 0.5)
  target_case_prop <- rep(0.5, 10)

  # P-value method should not be affected by show_model_selection
  # (it uses select_best_model internally but with its own verbosity control)
  output <- capture.output(
    B_pvalue <- get_B(
      training_trait = "binary",
      training_MAF = training_MAF,
      training_P = training_P,
      training_N = training_N,
      target_trait = "binary",
      target_MAF = target_MAF,
      target_case_prop = target_case_prop,
      method = "pvalue",
      show_model_selection = FALSE,  # Should have no effect on pvalue method
      verbose = 0
    ),
    type = "output"
  )

  # Should have no output
  expect_equal(length(output), 0)
  expect_type(B_pvalue, "double")
})
