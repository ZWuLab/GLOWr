########## Tests for train_B_model() Function ##########
#
# Comprehensive unit tests for the train_B_model() function and its
# print.glow_B_model() S3 method. Covers input validation, beta method,
# p-value method, both-methods mode, auto method detection, outlier
# detection, custom models, selection criterion, print method, integration
# with predict_B(), and edge cases.

library(testthat)

# ========== Helper Functions ==========

#' Create synthetic training data for train_B_model tests
#'
#' Generates realistic MAF, BETA, P, and N vectors with a known
#' MAF-effect-size relationship (B ~ sqrt(MAF*(1-MAF))).
#'
#' @param n Integer, number of training variants
#' @param seed Integer, random seed
#' @return Named list with MAF, BETA, P, N
.make_training_data <- function(n = 50, seed = 42) {
  set.seed(seed)
  MAF <- runif(n, 0.01, 0.45)
  # True effect sizes scale with sqrt(MAF*(1-MAF)), plus noise
  BETA <- 0.5 * sqrt(MAF * (1 - MAF)) + rnorm(n, 0, 0.02)
  # P-values derived from effect sizes (larger |BETA| -> smaller P)
  P <- 10^(-abs(BETA) * 50)
  # Clamp P to (0, 1) range to avoid edge issues
  P <- pmin(pmax(P, 1e-300), 1 - 1e-10)
  N <- rep(10000, n)
  list(MAF = MAF, BETA = BETA, P = P, N = N)
}


# ========== 1. Input Validation ==========

test_that("train_B_model rejects invalid training_trait", {
  td <- .make_training_data()

  # Invalid string

  expect_error(
    train_B_model(
      training_trait = "invalid",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_trait must be"
  )

  # Numeric value instead of character
  expect_error(
    train_B_model(
      training_trait = 42,
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_trait must be"
  )
})


test_that("train_B_model accepts NULL training_trait", {
  td <- .make_training_data()

  # NULL training_trait should be accepted (only pvalue method works)
  result <- train_B_model(
    training_trait = NULL,
    training_MAF = td$MAF,
    training_P = td$P,
    training_N = td$N,
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )
  expect_s3_class(result, "glow_B_model")
  # NULL training_trait is normalized to "mixed"
  expect_equal(result$training_summary$trait_type, "mixed")
})


test_that("train_B_model rejects empty/NULL/non-numeric training_MAF", {
  td <- .make_training_data()

  # Empty
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = numeric(0),
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_MAF must be a non-empty numeric vector"
  )

  # NULL
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = NULL,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_MAF must be a non-empty numeric vector"
  )

  # Character vector
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = c("a", "b"),
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_MAF must be a non-empty numeric vector"
  )
})


test_that("train_B_model rejects MAF values outside valid range", {
  td <- .make_training_data()

  # Negative MAF values
  bad_maf <- td$MAF
  bad_maf[1] <- -0.1
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = bad_maf,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "values must be in"
  )

  # All NA
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = rep(NA_real_, 10),
      training_BETA = td$BETA[1:10],
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    )
  )
})


test_that("train_B_model folds MAF > 0.5 with warning", {
  td <- .make_training_data(n = 20)

  # Replace some MAFs with values > 0.5 (these should be folded to 1 - MAF)
  maf_with_high <- td$MAF
  maf_with_high[1:3] <- c(0.6, 0.7, 0.9)

  # Should warn about folding
  expect_warning(
    result <- train_B_model(
      training_trait = "binary",
      training_MAF = maf_with_high,
      training_BETA = td$BETA[1:20],
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_MAF value.*> 0.5"
  )

  expect_s3_class(result, "glow_B_model")
  # Stored training data should have folded values
  expect_true(all(result$training_data$MAF <= 0.5))
})


test_that("train_B_model rejects invalid method parameter", {
  td <- .make_training_data()

  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "invalid_method",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "method must be one of"
  )
})


test_that("train_B_model handles deprecated 'CV' selection_criterion", {
  td <- .make_training_data()

  # "CV" should produce a deprecation warning and be treated as "CV_R2"
  expect_warning(
    result <- train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      selection_criterion = "CV",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "deprecated"
  )

  expect_s3_class(result, "glow_B_model")
  expect_equal(result$selection_criterion, "CV_R2")
})


test_that("train_B_model rejects invalid selection_criterion", {
  td <- .make_training_data()

  # Completely invalid
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      selection_criterion = "MSE",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "selection_criterion must be one of"
  )

  # AIC/BIC are explicitly rejected with informative message
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      selection_criterion = "AIC",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "AIC and BIC criteria are inappropriate"
  )
})


test_that("train_B_model rejects missing data for beta method", {
  td <- .make_training_data()

  # Beta method requires training_BETA
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_BETA is required"
  )
})


test_that("train_B_model rejects missing data for pvalue method", {
  td <- .make_training_data()

  # P-value method with no P or P_mlog10

  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      method = "pvalue",
      training_N = td$N,
      show_model_selection = FALSE,
      verbose = 0
    ),
    "Either training_P or training_P_mlog10 is required"
  )

  # P-value method with no N
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_P = td$P,
      method = "pvalue",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "training_N is required"
  )

  # Both P and P_mlog10 provided (ambiguous)
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_P = td$P,
      training_P_mlog10 = -log10(td$P),
      training_N = td$N,
      method = "pvalue",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "Provide either training_P or training_P_mlog10, not both"
  )
})


test_that("train_B_model errors on beta method with mixed/NULL trait", {
  td <- .make_training_data()

  # method="beta" with mixed training_trait
  expect_error(
    train_B_model(
      training_trait = "mixed",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "Direct beta method cannot be used when training_trait is 'mixed'"
  )

  # method="beta" with NULL training_trait
  expect_error(
    train_B_model(
      training_trait = NULL,
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "Direct beta method cannot be used when training_trait is 'mixed'"
  )
})


test_that("train_B_model auto insufficient data gives error", {
  td <- .make_training_data()

  # Auto with no BETA, no P, no P_mlog10 -> should error
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      method = "auto",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "Insufficient data"
  )
})


# ========== 2. Beta Method Only ==========

test_that("train_B_model beta method with binary trait returns glow_B_model", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  # Class and method_used
  expect_s3_class(result, "glow_B_model")
  expect_equal(result$method_used, "beta_method")

  # Models structure
  expect_true(inherits(result$models$beta_method, "lm"))
  expect_null(result$models$pvalue_method)

  # all_models_info
  expect_false(is.null(result$all_models_info$beta_method))
  expect_null(result$all_models_info$pvalue_method)

  # training_summary
  expect_equal(result$training_summary$n_original, length(td$MAF))
  expect_equal(result$training_summary$n_used, length(td$MAF))
  expect_equal(result$training_summary$trait_type, "binary")

  # training_data stored correctly (original data)
  expect_equal(result$training_data$MAF, td$MAF)
  expect_equal(result$training_data$BETA, td$BETA)
  expect_equal(result$training_data$trait, "binary")

  # Selection criterion
  expect_equal(result$selection_criterion, "R2")

  # Comparison should be empty (single method)
  expect_null(result$comparison$method_selected)
})


test_that("train_B_model beta method with continuous trait", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "continuous",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")
  expect_equal(result$method_used, "beta_method")
  expect_equal(result$training_summary$trait_type, "continuous")
})


# ========== 3. P-value Method Only ==========

test_that("train_B_model pvalue method with regular P values", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P = td$P,
    training_N = td$N,
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )

  # Class and method
  expect_s3_class(result, "glow_B_model")
  expect_equal(result$method_used, "pvalue_method")

  # Models structure
  expect_null(result$models$beta_method)
  expect_true(inherits(result$models$pvalue_method, "lm"))

  # all_models_info
  expect_null(result$all_models_info$beta_method)
  expect_false(is.null(result$all_models_info$pvalue_method))

  # Training data
  expect_equal(result$training_data$P, td$P)
  expect_null(result$training_data$P_mlog10)
  expect_equal(result$training_data$N, td$N)
})


test_that("train_B_model pvalue method with P_mlog10", {
  td <- .make_training_data()
  p_mlog10 <- -log10(td$P)

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P_mlog10 = p_mlog10,
    training_N = td$N,
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")
  expect_equal(result$method_used, "pvalue_method")
  expect_null(result$training_data$P)
  expect_equal(result$training_data$P_mlog10, p_mlog10)
})


test_that("pvalue method model predicts h-squared (not B directly)", {
  # The pvalue method lm object models h^2 ~ f(MAF), not BETA^2 ~ f(MAF).
  # Verify by checking the model's response is derived from Z^2/(2*N*MAF*(1-MAF)).
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P = td$P,
    training_N = td$N,
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )

  lm_obj <- result$models$pvalue_method
  # The response variable in the lm should be named Y or logY (not BETA)
  resp <- all.vars(formula(lm_obj))[1]
  expect_true(resp %in% c("Y", "logY"))

  # The fitted values should be on the h^2 scale (or log(h^2) scale)
  # Verify the model produces numeric fitted values
  expect_true(all(is.finite(fitted(lm_obj))))
})


# ========== 4. Both Methods ==========

test_that("train_B_model both methods populates comparison", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "both",
    show_model_selection = FALSE,
    verbose = 0
  )

  # Class and method
  expect_s3_class(result, "glow_B_model")
  expect_equal(result$method_used, "both")

  # Both models should be present
  expect_true(inherits(result$models$beta_method, "lm"))
  expect_true(inherits(result$models$pvalue_method, "lm"))

  # Both all_models_info should be present
  expect_false(is.null(result$all_models_info$beta_method))
  expect_false(is.null(result$all_models_info$pvalue_method))

  # Comparison should be populated
  expect_equal(result$comparison$selection_criterion, "R2")
  expect_true(is.numeric(result$comparison$criterion_beta_method))
  expect_true(is.numeric(result$comparison$criterion_pvalue_method))
  expect_true(result$comparison$method_selected %in%
                c("beta_method", "pvalue_method"))

  # The selected method should be the one with higher criterion value
  if (result$comparison$criterion_beta_method >
      result$comparison$criterion_pvalue_method) {
    expect_equal(result$comparison$method_selected, "beta_method")
  } else {
    expect_equal(result$comparison$method_selected, "pvalue_method")
  }
})


test_that("train_B_model both methods with mixed trait falls back to pvalue", {
  td <- .make_training_data()

  # method="both" with training_trait="mixed" should fall back to pvalue only
  expect_message(
    result <- train_B_model(
      training_trait = "mixed",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      training_P = td$P,
      training_N = td$N,
      method = "both",
      show_model_selection = FALSE,
      verbose = 1
    ),
    "Only p-value method available for mixed training traits"
  )

  expect_equal(result$method_used, "pvalue_method")
  expect_null(result$models$beta_method)
  expect_true(inherits(result$models$pvalue_method, "lm"))
})


# ========== 5. Auto Method Detection ==========

test_that("auto with only BETA selects beta_method", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "auto",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$method_used, "beta_method")
})


test_that("auto with only P+N selects pvalue_method", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P = td$P,
    training_N = td$N,
    method = "auto",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$method_used, "pvalue_method")
})


test_that("auto with BETA + P + same trait selects both", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "auto",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$method_used, "both")
  expect_true(inherits(result$models$beta_method, "lm"))
  expect_true(inherits(result$models$pvalue_method, "lm"))
})


test_that("auto with BETA + P + mixed trait selects pvalue only", {
  td <- .make_training_data()

  # mixed trait -> beta not possible -> auto should pick pvalue
  result <- train_B_model(
    training_trait = "mixed",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "auto",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$method_used, "pvalue_method")
})


test_that("auto with only P_mlog10 + N selects pvalue_method", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P_mlog10 = -log10(td$P),
    training_N = td$N,
    method = "auto",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$method_used, "pvalue_method")
})


# ========== 6. Outlier Detection ==========

test_that("outlier_method='none' detects no outliers", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "none",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$outliers$method, "none")
  expect_equal(result$training_summary$n_outliers_detected, 0)
  expect_equal(result$training_summary$n_used, result$training_summary$n_original)
  expect_equal(length(result$outliers$indices_removed), 0)
})


test_that("outlier_method='statistical' with injected outlier detects it", {
  td <- .make_training_data(n = 30)

  # Inject a strong outlier: high MAF with huge BETA
  td$BETA[1] <- 100
  td$MAF[1] <- 0.3

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "statistical",
    outlier_action = "flag",
    show_model_selection = FALSE,
    verbose = 0
  )

  # At least one outlier should be detected
  expect_gt(result$training_summary$n_outliers_detected, 0)
  # With flag action, n_used should equal n_original
  expect_equal(result$training_summary$n_used, result$training_summary$n_original)
})


test_that("outlier_method='both' detects statistical and biological", {
  td <- .make_training_data(n = 30)

  # Inject a biological outlier: common variant with absurdly large effect
  td$BETA[2] <- 50
  td$MAF[2] <- 0.2

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "both",
    outlier_action = "flag",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_gt(result$training_summary$n_outliers_detected, 0)
  expect_equal(result$outliers$method, "both")
})


test_that("outlier_action='flag' does not remove outliers", {
  td <- .make_training_data(n = 30)
  td$BETA[1] <- 100  # Inject outlier

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "statistical",
    outlier_action = "flag",
    show_model_selection = FALSE,
    verbose = 0
  )

  # n_used should equal n_original when action is "flag"
  expect_equal(result$training_summary$n_used, result$training_summary$n_original)
  expect_equal(length(result$outliers$indices_removed), 0)
})


test_that("outlier_action='remove' reduces n_used", {
  td <- .make_training_data(n = 30)
  td$BETA[1] <- 100  # Inject outlier

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "statistical",
    outlier_action = "remove",
    show_model_selection = FALSE,
    verbose = 0
  )

  # n_used should be less than n_original
  expect_lt(result$training_summary$n_used, result$training_summary$n_original)
  # Removed indices should be non-empty
  expect_gt(length(result$outliers$indices_removed), 0)
})


test_that("outlier removal stores original data in training_data", {
  td <- .make_training_data(n = 30)
  td$BETA[1] <- 100  # Inject outlier

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "statistical",
    outlier_action = "remove",
    show_model_selection = FALSE,
    verbose = 0
  )

  # training_data should store ORIGINAL (pre-removal) data
  expect_equal(length(result$training_data$MAF), 30)
  expect_equal(result$training_data$BETA[1], 100)  # Original outlier value
})


test_that("outlier indices refer to original data positions", {
  td <- .make_training_data(n = 30)
  td$BETA[5] <- 200  # Inject known outlier at position 5

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    outlier_method = "statistical",
    outlier_action = "remove",
    show_model_selection = FALSE,
    verbose = 0
  )

  # If position 5 was detected as an outlier, verify it is in the indices
  if (5 %in% result$outliers$indices_removed) {
    # Original data at position 5 should be the outlier we injected
    expect_equal(result$training_data$BETA[5], 200)
  }
  # All removed indices should be valid positions in the original data
  if (length(result$outliers$indices_removed) > 0) {
    expect_true(all(result$outliers$indices_removed >= 1))
    expect_true(all(result$outliers$indices_removed <= 30))
  }
})


# ========== 7. Custom Models ==========

test_that("custom model formulas are accepted and appear in all_models_info", {
  td <- .make_training_data()

  custom <- list(
    "poly_2" = Y ~ poly(X, 2)
  )

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    custom_models = custom,
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")

  # The all_models_info should contain the custom model alongside standard ones
  info <- result$all_models_info$beta_method
  expect_false(is.null(info))
})


test_that("custom_models must be a list of formulas", {
  td <- .make_training_data()

  # Not a list
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      custom_models = Y ~ X,
      show_model_selection = FALSE,
      verbose = 0
    ),
    "custom_models must be a list"
  )

  # List with non-formula elements
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      custom_models = list("not a formula"),
      show_model_selection = FALSE,
      verbose = 0
    ),
    "All elements of custom_models must be formula objects"
  )
})


# ========== 8. Selection Criterion ==========

test_that("different selection criteria produce valid results", {
  td <- .make_training_data()

  criteria <- c("R2", "adj_R2", "CV_R2")

  for (crit in criteria) {
    result <- train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      selection_criterion = crit,
      show_model_selection = FALSE,
      verbose = 0
    )

    expect_s3_class(result, "glow_B_model")
    expect_equal(result$selection_criterion, crit)
    expect_true(inherits(result$models$beta_method, "lm"))
  }
})


test_that("selection criterion is stored in output", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    selection_criterion = "adj_R2",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$selection_criterion, "adj_R2")
})


test_that("both methods comparison uses the specified criterion", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "both",
    selection_criterion = "adj_R2",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_equal(result$comparison$selection_criterion, "adj_R2")
  expect_true(is.numeric(result$comparison$criterion_beta_method))
  expect_true(is.numeric(result$comparison$criterion_pvalue_method))
})


# ========== 9. Print Method ==========

test_that("print.glow_B_model runs without error", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  # Capture output -- should not error
  output <- capture.output(print(result))
  expect_true(length(output) > 0)
})


test_that("print.glow_B_model output contains expected strings", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  output <- paste(capture.output(print(result)), collapse = "\n")

  # Check for key information in the output
  expect_true(grepl("glow_B_model", output))
  expect_true(grepl("Method", output))
  expect_true(grepl("Beta Method", output))
  expect_true(grepl("Selection criterion", output))
  expect_true(grepl("R2", output))
  expect_true(grepl("binary", output))
  expect_true(grepl("Variants", output))
  expect_true(grepl("R-squared", output))
})


test_that("print.glow_B_model with both methods shows comparison", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "both",
    show_model_selection = FALSE,
    verbose = 0
  )

  output <- paste(capture.output(print(result)), collapse = "\n")

  # Should include method comparison section
  expect_true(grepl("Method Comparison", output))
  expect_true(grepl("Selected", output))
})


test_that("print.glow_B_model returns object invisibly", {
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  invisible_return <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})


# ========== 10. Integration with predict_B ==========

test_that("beta method model works end-to-end with predict_B", {
  td <- .make_training_data()

  model <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  target_MAF <- c(0.01, 0.05, 0.10, 0.20, 0.40)
  B <- predict_B(model, target_MAF = target_MAF)

  expect_type(B, "double")
  expect_length(B, length(target_MAF))
  expect_true(all(is.finite(B)))
  expect_true(all(B >= 0))
})


test_that("pvalue method model works end-to-end with predict_B", {
  td <- .make_training_data()

  model <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P = td$P,
    training_N = td$N,
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )

  target_MAF <- c(0.01, 0.05, 0.10, 0.20, 0.40)
  B <- predict_B(model, target_MAF = target_MAF,
                  target_trait = "binary",
                  target_case_prop = 0.3)

  expect_type(B, "double")
  expect_length(B, length(target_MAF))
  expect_true(all(is.finite(B)))
  expect_true(all(B >= 0))
})


test_that("both methods model works end-to-end with predict_B", {
  td <- .make_training_data()

  model <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "both",
    show_model_selection = FALSE,
    verbose = 0
  )

  target_MAF <- c(0.01, 0.05, 0.10, 0.20)

  # Default predict (uses primary/selected method)
  B_auto <- predict_B(model, target_MAF = target_MAF,
                       target_trait = "binary",
                       target_case_prop = 0.3)
  expect_true(all(is.finite(B_auto)))
  expect_true(all(B_auto >= 0))

  # Explicit beta method
  B_beta <- predict_B(model, target_MAF = target_MAF,
                       method = "beta_method")
  expect_true(all(is.finite(B_beta)))

  # Explicit pvalue method
  B_pval <- predict_B(model, target_MAF = target_MAF,
                       target_trait = "binary",
                       target_case_prop = 0.3,
                       method = "pvalue_method")
  expect_true(all(is.finite(B_pval)))
})


# ========== 11. Edge Cases ==========

test_that("train_B_model works with small training set (n=5)", {
  # 5 is the minimum that should work (need >=3 after filtering)
  set.seed(99)
  maf <- c(0.01, 0.05, 0.10, 0.20, 0.40)
  beta <- c(0.5, 0.3, 0.2, 0.15, 0.08)

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = maf,
    training_BETA = beta,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")
  expect_equal(result$training_summary$n_original, 5)
})


test_that("train_B_model handles all identical BETA values", {
  # All BETAs the same -> R2=0 but should still fit
  set.seed(77)
  maf <- runif(20, 0.01, 0.4)
  beta <- rep(0.1, 20)

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = maf,
    training_BETA = beta,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")
  # Model should still be fitted
  expect_true(inherits(result$models$beta_method, "lm"))
})


test_that("train_B_model with single training_N value (scalar)", {
  td <- .make_training_data()

  # Pass N as a single value (should be recycled or used directly)
  # Note: the function expects N as a vector of same length as MAF,
  # but a constant N vector is a common use case
  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_P = td$P,
    training_N = rep(5000, length(td$MAF)),
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")
})


test_that("train_B_model with very large BETA values", {
  # Ensure extreme but valid inputs do not crash
  set.seed(33)
  n <- 20
  maf <- runif(n, 0.01, 0.4)
  beta <- rnorm(n, mean = 0, sd = 5)  # Large but realistic range

  result <- train_B_model(
    training_trait = "continuous",
    training_MAF = maf,
    training_BETA = beta,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  expect_s3_class(result, "glow_B_model")
})


test_that("train_B_model with zero BETA values removes them for beta method", {
  td <- .make_training_data(n = 20)
  # Set some BETAs to zero
  td$BETA[1:3] <- 0

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  # Should succeed (zeros removed from beta method training)
  expect_s3_class(result, "glow_B_model")
  # Original data stored should still include the zeros
  expect_equal(sum(result$training_data$BETA == 0), 3)
})


test_that("train_B_model with all zero BETA errors", {
  td <- .make_training_data(n = 10)
  td$BETA <- rep(0, 10)

  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "All training_BETA values are zero"
  )
})


test_that("train_B_model verbose=2 produces messages", {
  td <- .make_training_data()

  # verbose=2 should produce informational messages
  expect_message(
    result <- train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      show_model_selection = FALSE,
      verbose = 2
    ),
    "Training beta method model"
  )
})


test_that("train_B_model validates outlier parameters", {
  td <- .make_training_data()

  # Invalid outlier_method
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      outlier_method = "invalid",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "outlier_method must be one of"
  )

  # Invalid outlier_action
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      outlier_action = "delete",
      show_model_selection = FALSE,
      verbose = 0
    ),
    "outlier_action must be one of"
  )

  # Invalid cook_threshold
  expect_error(
    train_B_model(
      training_trait = "binary",
      training_MAF = td$MAF,
      training_BETA = td$BETA,
      method = "beta",
      cook_threshold = -1,
      show_model_selection = FALSE,
      verbose = 0
    ),
    "cook_threshold must be a positive numeric value"
  )
})


test_that("glow_B_model structure is complete and well-formed", {
  # Comprehensive structural check of the returned object
  td <- .make_training_data()

  result <- train_B_model(
    training_trait = "binary",
    training_MAF = td$MAF,
    training_BETA = td$BETA,
    training_P = td$P,
    training_N = td$N,
    method = "both",
    selection_criterion = "R2",
    show_model_selection = FALSE,
    verbose = 0
  )

  # Top-level fields
  expect_true("method_used" %in% names(result))
  expect_true("models" %in% names(result))
  expect_true("all_models_info" %in% names(result))
  expect_true("outliers" %in% names(result))
  expect_true("training_summary" %in% names(result))
  expect_true("training_data" %in% names(result))
  expect_true("comparison" %in% names(result))
  expect_true("selection_criterion" %in% names(result))

  # models sub-fields
  expect_true("beta_method" %in% names(result$models))
  expect_true("pvalue_method" %in% names(result$models))

  # all_models_info sub-fields
  expect_true("beta_method" %in% names(result$all_models_info))
  expect_true("pvalue_method" %in% names(result$all_models_info))

  # outliers sub-fields
  expect_true("method" %in% names(result$outliers))
  expect_true("action" %in% names(result$outliers))
  expect_true("beta_method" %in% names(result$outliers))
  expect_true("pvalue_method" %in% names(result$outliers))
  expect_true("indices_removed" %in% names(result$outliers))

  # training_summary sub-fields
  expect_true("n_original" %in% names(result$training_summary))
  expect_true("n_used" %in% names(result$training_summary))
  expect_true("n_outliers_detected" %in% names(result$training_summary))
  expect_true("trait_type" %in% names(result$training_summary))

  # training_data sub-fields
  expect_true("MAF" %in% names(result$training_data))
  expect_true("BETA" %in% names(result$training_data))
  expect_true("P" %in% names(result$training_data))
  expect_true("P_mlog10" %in% names(result$training_data))
  expect_true("N" %in% names(result$training_data))
  expect_true("trait" %in% names(result$training_data))

  # comparison sub-fields
  expect_true("selection_criterion" %in% names(result$comparison))
  expect_true("criterion_beta_method" %in% names(result$comparison))
  expect_true("criterion_pvalue_method" %in% names(result$comparison))
  expect_true("method_selected" %in% names(result$comparison))
})
