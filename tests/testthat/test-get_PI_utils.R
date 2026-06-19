# Tests for PI Model Utility Functions

# ==============================================================================
# Test Data Setup
# ==============================================================================

# Create mock LASSO model (mimics glmnet output)
.create_mock_lasso_model <- function(coefs, lambda = 0.01) {
  # coefs should be a named numeric vector
  beta <- Matrix::Matrix(coefs, ncol = 1, sparse = TRUE)
  rownames(beta) <- names(coefs)

  model <- list(
    beta = beta,
    lambda = lambda
  )
  class(model) <- c("elnet", "glmnet")
  return(model)
}

# Create real GLM model with named features for testing
# Uses real glm() so summary.glm works properly in package namespace
.create_real_glm_model <- function(seed = 42, n = 200) {
  set.seed(seed)
  synthetic_data <- data.frame(
    y = rbinom(n, 1, 0.5),
    feat1 = rnorm(n),
    feat2 = rnorm(n),
    feat3 = rnorm(n)
  )
  glm(y ~ feat1 + feat2 + feat3, data = synthetic_data, family = binomial())
}


# ==============================================================================
# Tests: .detect_PI_model_type
# ==============================================================================

test_that(".detect_PI_model_type correctly identifies LASSO models", {
  coefs <- c(a = 0.1, b = 0.2, c = 0)
  mock_model <- .create_mock_lasso_model(coefs)

  result <- GLOWr:::.detect_PI_model_type(mock_model)
  expect_equal(result, "LASSO")
})

test_that(".detect_PI_model_type correctly identifies GLM models", {
  real_model <- .create_real_glm_model(seed = 42)

  result <- GLOWr:::.detect_PI_model_type(real_model)
  expect_equal(result, "GLM")
})


# ==============================================================================
# Tests: .compute_auc
# ==============================================================================

test_that(".compute_auc returns 1.0 for perfect classifier", {
  labels <- c(1, 1, 1, 0, 0, 0)
  predictions <- c(0.9, 0.8, 0.7, 0.3, 0.2, 0.1)

  auc <- GLOWr:::.compute_auc(labels, predictions)
  expect_equal(auc, 1.0)
})

test_that(".compute_auc returns 0.0 for perfectly wrong classifier", {
  labels <- c(1, 1, 1, 0, 0, 0)
  predictions <- c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9)

  auc <- GLOWr:::.compute_auc(labels, predictions)
  expect_equal(auc, 0.0)
})

test_that(".compute_auc returns ~0.5 for random predictions", {
  set.seed(123)
  labels <- c(rep(1, 100), rep(0, 100))
  predictions <- runif(200)

  auc <- GLOWr:::.compute_auc(labels, predictions)
  # Should be approximately 0.5 with some variance
  expect_true(auc > 0.3 && auc < 0.7)
})

test_that(".compute_auc handles NA values", {
  labels <- c(1, 1, NA, 0, 0)
  predictions <- c(0.9, 0.8, 0.5, 0.2, 0.1)

  auc <- GLOWr:::.compute_auc(labels, predictions)
  expect_true(!is.na(auc))
  expect_true(auc >= 0 && auc <= 1)
})

test_that(".compute_auc warns when only one class present", {
  labels <- c(1, 1, 1)
  predictions <- c(0.9, 0.8, 0.7)

  expect_warning(
    auc <- GLOWr:::.compute_auc(labels, predictions),
    "need both positive and negative"
  )
  expect_true(is.na(auc))
})


# ==============================================================================
# Tests: .compute_roc_coords
# ==============================================================================

test_that(".compute_roc_coords returns correct structure", {
  labels <- c(1, 1, 0, 0)
  predictions <- c(0.9, 0.7, 0.4, 0.2)

  roc <- GLOWr:::.compute_roc_coords(labels, predictions)

  expect_true(is.list(roc))
  expect_true("fpr" %in% names(roc))
  expect_true("tpr" %in% names(roc))
  expect_equal(length(roc$fpr), length(roc$tpr))
})

test_that(".compute_roc_coords starts at (0,0) and ends at (1,1)", {
  labels <- c(1, 1, 0, 0)
  predictions <- c(0.9, 0.7, 0.4, 0.2)

  roc <- GLOWr:::.compute_roc_coords(labels, predictions)

  expect_equal(roc$fpr[1], 0)
  expect_equal(roc$tpr[1], 0)
  expect_equal(roc$fpr[length(roc$fpr)], 1)
  expect_equal(roc$tpr[length(roc$tpr)], 1)
})


# ==============================================================================
# Tests: .extract_coef_lasso
# ==============================================================================

test_that(".extract_coef_lasso extracts coefficients correctly", {
  coefs <- c(feat1 = 0.1, feat2 = 0, feat3 = -0.2)
  mock_model <- .create_mock_lasso_model(coefs)

  result <- GLOWr:::.extract_coef_lasso(mock_model)

  expect_equal(length(result), 3)
  expect_equal(names(result), c("feat1", "feat2", "feat3"))
  expect_equal(result["feat1"], c(feat1 = 0.1))
  expect_equal(result["feat2"], c(feat2 = 0))
})


# ==============================================================================
# Tests: load_PI_models
# ==============================================================================

test_that("load_PI_models errors on non-existent directory", {
  expect_error(
    load_PI_models("/nonexistent/path"),
    "does not exist"
  )
})

test_that("load_PI_models errors on empty directory", {
  # Create temp empty directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  expect_error(
    load_PI_models(temp_dir),
    "No model files"
  )
})

test_that("load_PI_models loads models from directory", {
  pi_dir <- Sys.getenv("GLOW_PI_MODELS_DIR", unset = "")
  skip_if_not(nzchar(pi_dir) && dir.exists(pi_dir))

  result <- load_PI_models(pi_dir)

  expect_true(is.list(result))
  expect_true("models" %in% names(result))
  expect_true("model_type" %in% names(result))
  expect_true("n_models" %in% names(result))
  expect_equal(result$n_models, length(result$models))
  expect_true(result$model_type %in% c("LASSO", "GLM"))
})


# ==============================================================================
# Tests: summarize_PI_coefficients
# ==============================================================================

test_that("summarize_PI_coefficients returns correct structure", {
  # Create mock models
  coefs1 <- c(feat1 = 0.1, feat2 = 0.2, feat3 = 0)
  coefs2 <- c(feat1 = 0.15, feat2 = 0, feat3 = -0.1)
  models <- list(
    .create_mock_lasso_model(coefs1),
    .create_mock_lasso_model(coefs2)
  )

  result <- summarize_PI_coefficients(models, model_type = "LASSO")

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("coef_matrix" %in% names(result))
  expect_true("model_type" %in% names(result))
  expect_true("n_models" %in% names(result))

  expect_true(is.data.frame(result$summary))
  expect_true(is.matrix(result$coef_matrix))
  expect_equal(result$n_models, 2)
})

test_that("summarize_PI_coefficients computes selection frequency correctly", {
  # Create 4 models where feat1 is always selected, feat2 is selected 2/4 times
  coefs1 <- c(feat1 = 0.1, feat2 = 0.2, feat3 = 0)
  coefs2 <- c(feat1 = 0.15, feat2 = 0, feat3 = 0)
  coefs3 <- c(feat1 = 0.12, feat2 = 0.1, feat3 = 0)
  coefs4 <- c(feat1 = 0.08, feat2 = 0, feat3 = 0)

  models <- list(
    .create_mock_lasso_model(coefs1),
    .create_mock_lasso_model(coefs2),
    .create_mock_lasso_model(coefs3),
    .create_mock_lasso_model(coefs4)
  )

  result <- summarize_PI_coefficients(models, model_type = "LASSO")

  # feat1 selected 4/4 = 100%
  feat1_row <- result$summary[result$summary$feature == "feat1", ]
  expect_equal(feat1_row$n_selected, 4)
  expect_equal(feat1_row$pct_selected, 100)

  # feat2 selected 2/4 = 50%
  feat2_row <- result$summary[result$summary$feature == "feat2", ]
  expect_equal(feat2_row$n_selected, 2)
  expect_equal(feat2_row$pct_selected, 50)

  # feat3 selected 0/4 = 0%
  feat3_row <- result$summary[result$summary$feature == "feat3", ]
  expect_equal(feat3_row$n_selected, 0)
  expect_equal(feat3_row$pct_selected, 0)
})


# ==============================================================================
# Tests: .extract_coef_pval_glm
# ==============================================================================

test_that(".extract_coef_pval_glm extracts coefficients and p-values correctly", {
  real_model <- .create_real_glm_model(seed = 42)

  result <- GLOWr:::.extract_coef_pval_glm(real_model)

  # Check structure
  expect_true(is.list(result))
  expect_true("coef" %in% names(result))
  expect_true("pval" %in% names(result))
  expect_equal(length(result$coef), 3)
  expect_equal(length(result$pval), 3)
  expect_equal(names(result$coef), c("feat1", "feat2", "feat3"))
  expect_equal(names(result$pval), c("feat1", "feat2", "feat3"))

  # Verify values match what summary() returns (excluding intercept)
  coef_table <- summary(real_model)$coefficients
  expected_coefs <- coef_table[-1, "Estimate"]
  expected_pvals <- coef_table[-1, "Pr(>|z|)"]
  expect_equal(as.numeric(result$coef), as.numeric(expected_coefs))
  expect_equal(as.numeric(result$pval), as.numeric(expected_pvals))
})


# ==============================================================================
# Tests: summarize_PI_coefficients for GLM (significance frequency)
# ==============================================================================

test_that("summarize_PI_coefficients computes significance frequency for GLM", {
  # Create 4 real GLM models with different seeds
  models <- list(
    .create_real_glm_model(seed = 101),
    .create_real_glm_model(seed = 202),
    .create_real_glm_model(seed = 303),
    .create_real_glm_model(seed = 404)
  )

  result <- summarize_PI_coefficients(models, model_type = "GLM")

  # Check structure
  expect_true("pval_matrix" %in% names(result))
  expect_true(!is.null(result$pval_matrix))
  expect_true("n_significant" %in% names(result$summary))
  expect_true("pct_significant" %in% names(result$summary))
  expect_true("mean_neg_log10_pval" %in% names(result$summary))

  # Check summary has correct features
  expect_equal(sort(result$summary$feature), c("feat1", "feat2", "feat3"))

  # Check pval_matrix dimensions (4 models x 3 features)
  expect_equal(nrow(result$pval_matrix), 4)
  expect_equal(ncol(result$pval_matrix), 3)

  # Check significance counts are within valid range
  expect_true(all(result$summary$n_significant >= 0))
  expect_true(all(result$summary$n_significant <= 4))
  expect_true(all(result$summary$pct_significant >= 0))
  expect_true(all(result$summary$pct_significant <= 100))
})

test_that("summarize_PI_coefficients returns NULL pval_matrix for LASSO", {
  coefs1 <- c(feat1 = 0.1, feat2 = 0.2, feat3 = 0)
  coefs2 <- c(feat1 = 0.15, feat2 = 0, feat3 = -0.1)
  models <- list(
    .create_mock_lasso_model(coefs1),
    .create_mock_lasso_model(coefs2)
  )

  result <- summarize_PI_coefficients(models, model_type = "LASSO")

  expect_true(is.null(result$pval_matrix))
  expect_true("n_selected" %in% names(result$summary))
  expect_true("pct_selected" %in% names(result$summary))
  expect_false("n_significant" %in% names(result$summary))
})


# ==============================================================================
# Tests: predict_PI_ensemble
# ==============================================================================

test_that("predict_PI_ensemble returns correct structure", {
  skip("Requires real models for prediction - use integration tests")
})


# ==============================================================================
# Tests: evaluate_PI_models
# ==============================================================================

test_that("evaluate_PI_models returns correct structure", {
  skip("Requires real models and data - use integration tests")
})


# ==============================================================================
# Tests: plot functions (basic execution)
# ==============================================================================

test_that("plot_PI_coefficient_summary runs without error for LASSO", {
  coefs1 <- c(feat1 = 0.1, feat2 = 0.2, feat3 = 0)
  coefs2 <- c(feat1 = 0.15, feat2 = 0, feat3 = -0.1)
  models <- list(
    .create_mock_lasso_model(coefs1),
    .create_mock_lasso_model(coefs2)
  )

  coef_summary <- summarize_PI_coefficients(models, model_type = "LASSO")

  # Should not error (simplified interface - bar plot only, no type parameter)
  expect_silent({
    pdf(NULL)  # Null device to avoid actual plotting
    plot_PI_coefficient_summary(coef_summary)
    dev.off()
  })
})

test_that("plot_PI_coefficient_summary runs without error for GLM", {
  models <- list(
    .create_real_glm_model(seed = 501),
    .create_real_glm_model(seed = 502)
  )

  coef_summary <- summarize_PI_coefficients(models, model_type = "GLM")

  # Should not error
  expect_silent({
    pdf(NULL)
    plot_PI_coefficient_summary(coef_summary)
    dev.off()
  })
})

test_that("plot_PI_roc runs without error on mock data", {
  # Create mock evaluation result
  mock_eval <- list(
    roc_data = list(
      labels = c(1, 1, 1, 0, 0, 0),
      predictions = matrix(c(0.9, 0.8, 0.7, 0.3, 0.2, 0.1,
                             0.85, 0.75, 0.65, 0.35, 0.25, 0.15),
                           nrow = 6, ncol = 2),
      ensemble_pi = c(0.875, 0.775, 0.675, 0.325, 0.225, 0.125)
    ),
    ensemble = list(auc = 1.0),
    summary = list(mean_auc = 1.0)
  )

  expect_silent({
    pdf(NULL)
    plot_PI_roc(mock_eval)
    dev.off()
  })
})
