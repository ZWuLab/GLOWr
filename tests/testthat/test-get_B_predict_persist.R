########## Tests for predict_B(), save_B_model(), load_B_model() ##########
#
# This file contains unit tests for the B estimation predict and persist
# functions that support the "train once -> save -> predict many times"
# workflow.

library(testthat)

# ========== Helper Functions ==========

#' Create a minimal glow_B_model object for testing (beta method)
#' @param model_id Integer 1-8, which candidate model to use
#' @return A glow_B_model object with the specified model
.make_test_glow_B_model <- function(model_id = 4L, method = "beta_method") {
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)

  if (model_id <= 4L) {
    # Y-response models
    x <- switch(as.character(model_id),
      "1" = maf,
      "2" = maf * (1 - maf),
      "3" = log(maf),
      "4" = log(maf * (1 - maf))
    )
    var_name <- c("X", "fX", "logX", "logfX")[model_id]
    Y <- -0.14 * x + 0.5 + rnorm(50, 0, 0.02)
    df <- data.frame(Y = Y)
    df[[var_name]] <- x
    fmla <- as.formula(paste("Y ~", var_name))
    lm_obj <- lm(fmla, data = df)
  } else {
    # log(Y)-response models
    x_type <- model_id - 4L
    x <- switch(as.character(x_type),
      "1" = maf,
      "2" = maf * (1 - maf),
      "3" = log(maf),
      "4" = log(maf * (1 - maf))
    )
    var_name <- c("X", "fX", "logX", "logfX")[x_type]
    logY <- -0.14 * x - 0.20 + rnorm(50, 0, 0.1)
    df <- data.frame(logY = logY)
    df[[var_name]] <- x
    fmla <- as.formula(paste("logY ~", var_name))
    lm_obj <- lm(fmla, data = df)
  }

  attr(lm_obj, "model_id") <- model_id

  structure(list(
    method_used = method,
    models = list(
      beta_method = if (method %in% c("beta_method", "both")) lm_obj else NULL,
      pvalue_method = if (method %in% c("pvalue_method", "both")) lm_obj else NULL
    ),
    all_models_info = list(
      beta_method = NULL,
      pvalue_method = NULL
    ),
    outliers = list(
      method = "none", action = "flag",
      beta_method = NULL, pvalue_method = NULL,
      indices_removed = integer(0)
    ),
    training_summary = list(
      n_original = 50L, n_used = 50L,
      n_outliers_detected = 0L,
      trait_type = "binary"
    ),
    training_data = list(
      MAF = maf, BETA = NULL, P = NULL,
      P_mlog10 = NULL, N = NULL, trait = "binary"
    ),
    comparison = list(
      selection_criterion = NULL,
      criterion_beta_method = NULL,
      criterion_pvalue_method = NULL,
      method_selected = NULL
    ),
    selection_criterion = "R2"
  ), class = "glow_B_model")
}


#' Create a glow_B_model with both methods for testing
#' @return A glow_B_model with method_used="both"
.make_test_glow_B_model_both <- function() {
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)

  # Beta method model
  logX <- log(maf)
  Y <- -0.14 * logX + 0.5 + rnorm(50, 0, 0.02)
  lm_beta <- lm(Y ~ logX, data = data.frame(Y = Y, logX = logX))
  attr(lm_beta, "model_id") <- 3L

  # P-value method model (different coefficients)
  Y2 <- -0.2 * logX + 0.8 + rnorm(50, 0, 0.03)
  lm_pvalue <- lm(Y ~ logX, data = data.frame(Y = Y2, logX = logX))
  attr(lm_pvalue, "model_id") <- 3L

  structure(list(
    method_used = "both",
    models = list(
      beta_method = lm_beta,
      pvalue_method = lm_pvalue
    ),
    all_models_info = list(
      beta_method = NULL,
      pvalue_method = NULL
    ),
    outliers = list(
      method = "none", action = "flag",
      beta_method = NULL, pvalue_method = NULL,
      indices_removed = integer(0)
    ),
    training_summary = list(
      n_original = 50L, n_used = 50L,
      n_outliers_detected = 0L,
      trait_type = "binary"
    ),
    training_data = list(
      MAF = maf, BETA = NULL, P = NULL,
      P_mlog10 = NULL, N = NULL, trait = "binary"
    ),
    comparison = list(
      selection_criterion = "R2",
      criterion_beta_method = 0.85,
      criterion_pvalue_method = 0.80,
      method_selected = "beta_method"
    ),
    selection_criterion = "R2"
  ), class = "glow_B_model")
}


# ========== predict_B: glow_B_model tests (beta method) ==========

test_that("predict_B works with glow_B_model (beta method)", {
  model <- .make_test_glow_B_model()
  target_maf <- c(0.01, 0.05, 0.1, 0.2, 0.3)

  B <- predict_B(model, target_maf)
  expect_true(is.numeric(B))
  expect_equal(length(B), 5)
  expect_true(all(B >= 0))
})

test_that("predict_B works with raw lm object", {
  set.seed(42)
  x <- log(runif(30, 0.01, 0.3))
  Y <- -0.5 * x + 1 + rnorm(30, 0, 0.1)
  lm_obj <- lm(Y ~ logX, data = data.frame(Y = Y, logX = x))
  attr(lm_obj, "model_id") <- 3L  # Model 3: Y ~ log(X)

  B <- predict_B(lm_obj, target_MAF = c(0.01, 0.05, 0.1))
  expect_true(is.numeric(B))
  expect_equal(length(B), 3)
  expect_true(all(B >= 0))
})

test_that("predict_B handles MAF = 0", {
  model <- .make_test_glow_B_model()

  B <- predict_B(model, target_MAF = c(0, 0.05, 0.5))
  expect_equal(B[1], 0)
  expect_true(all(B >= 0))
})

test_that("predict_B handles all-zero MAF input", {
  model <- .make_test_glow_B_model()
  B <- predict_B(model, target_MAF = c(0, 0, 0))
  expect_equal(B, c(0, 0, 0))
})

test_that("predict_B handles empty input", {
  model <- .make_test_glow_B_model()
  B <- predict_B(model, target_MAF = numeric(0))
  expect_equal(length(B), 0)
  expect_true(is.numeric(B))
})

test_that("predict_B errors on NA in target_MAF", {
  model <- .make_test_glow_B_model()
  expect_error(predict_B(model, c(0.1, NA, 0.2)), "NA")
})

test_that("predict_B folds MAF > 0.5 with warning", {
  model <- .make_test_glow_B_model()

  # MAF > 0.5 should be folded to 1-MAF with a warning, not an error

  expect_warning(
    B <- predict_B(model, c(0.1, 0.6)),
    "target_MAF value\\(s\\) > 0.5 detected"
  )
  expect_true(is.numeric(B))
  expect_equal(length(B), 2)
  expect_true(all(B >= 0))

  # The folded result for 0.6 should equal the direct result for 0.4
  B_direct <- predict_B(model, c(0.1, 0.4))
  expect_equal(B, B_direct)
})

test_that("predict_B errors on negative MAF", {
  model <- .make_test_glow_B_model()
  expect_error(predict_B(model, c(-0.1, 0.1)), "0")
})

test_that("predict_B errors on wrong object type", {
  expect_error(predict_B("not_a_model", c(0.1)), "glow_B_model|lm")
})

test_that("predict_B errors on lm without model_id", {
  lm_obj <- lm(mpg ~ wt, data = mtcars)
  expect_error(predict_B(lm_obj, c(0.1)), "model_id")
})

test_that("predict_B errors on invalid method", {
  model <- .make_test_glow_B_model()
  expect_error(predict_B(model, c(0.1), method = "invalid"), "method must be")
})

test_that("predict_B back-transforms correctly for Y-response models (1-4)", {
  for (mid in 1:4) {
    model <- .make_test_glow_B_model(model_id = mid)
    maf <- c(0.05, 0.1, 0.2)
    B <- predict_B(model, maf)
    expect_true(all(is.finite(B)), info = paste("model_id =", mid))
    expect_true(all(B >= 0), info = paste("model_id =", mid))
  }
})

test_that("predict_B back-transforms correctly for log-Y models (5-8)", {
  for (mid in 5:8) {
    model <- .make_test_glow_B_model(model_id = mid)
    maf <- c(0.05, 0.1, 0.2)
    B <- predict_B(model, maf)
    expect_true(all(is.finite(B)), info = paste("model_id =", mid))
    expect_true(all(B > 0), info = paste("model_id =", mid))
  }
})


# ========== predict_B: method parameter selection ==========

test_that("predict_B respects method parameter for glow_B_model", {
  model <- .make_test_glow_B_model()
  # Default (auto) should work
  B_auto <- predict_B(model, c(0.1, 0.2))
  B_beta <- predict_B(model, c(0.1, 0.2), method = "beta_method")
  expect_equal(B_auto, B_beta)

  # Requesting unavailable method should error
  expect_error(predict_B(model, c(0.1), method = "pvalue_method"),
               "not available")
})

test_that("predict_B handles method_used='both' in glow_B_model", {
  model <- .make_test_glow_B_model_both()

  # Auto should use beta_method (the selected method)
  B <- predict_B(model, c(0.1, 0.2))
  expect_true(is.numeric(B))
  expect_equal(length(B), 2)

  # Can explicitly request either method
  B_beta <- predict_B(model, c(0.1, 0.2), method = "beta_method")
  expect_true(is.numeric(B_beta))

  # P-value method without target_trait issues a warning (returns sqrt(h2))
  expect_warning(
    B_pval <- predict_B(model, c(0.1, 0.2), method = "pvalue_method"),
    "target_trait"
  )
  expect_true(is.numeric(B_pval))

  # P-value method with proper params works without warning
  B_pval2 <- predict_B(model, c(0.1, 0.2), method = "pvalue_method",
                        target_trait = "binary", target_case_prop = 0.3)
  expect_true(is.numeric(B_pval2))
})


# ========== predict_B: p-value method h2-to-B conversion ==========

test_that("predict_B with pvalue model and binary trait conversion", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  B <- predict_B(model, c(0.05, 0.1, 0.2),
                 target_trait = "binary",
                 target_case_prop = 0.3)
  expect_true(is.numeric(B))
  expect_equal(length(B), 3)
  expect_true(all(B >= 0))
})

test_that("predict_B with pvalue model and continuous trait conversion", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  B <- predict_B(model, c(0.05, 0.1, 0.2),
                 target_trait = "continuous",
                 target_SE = 0.5)
  expect_true(is.numeric(B))
  expect_equal(length(B), 3)
  expect_true(all(B >= 0))
})

test_that("predict_B with pvalue model warns when target_trait is NULL", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  expect_warning(
    B <- predict_B(model, c(0.05, 0.1, 0.2)),
    "target_trait"
  )
  expect_true(is.numeric(B))
  expect_equal(length(B), 3)
})

test_that("predict_B pvalue method errors without required binary params", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  expect_error(
    predict_B(model, c(0.1), target_trait = "binary"),
    "target_case_prop"
  )
})

test_that("predict_B pvalue method errors without required continuous params", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  expect_error(
    predict_B(model, c(0.1), target_trait = "continuous"),
    "target_SE"
  )
})

test_that("predict_B pvalue method validates target_case_prop range", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  expect_error(
    predict_B(model, c(0.1), target_trait = "binary",
              target_case_prop = 0),
    "target_case_prop"
  )

  expect_error(
    predict_B(model, c(0.1), target_trait = "binary",
              target_case_prop = 1),
    "target_case_prop"
  )
})

test_that("predict_B pvalue method accepts scalar target_SE", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  # Scalar should be recycled
  B <- predict_B(model, c(0.05, 0.1, 0.2),
                 target_trait = "continuous",
                 target_SE = 0.5)
  expect_equal(length(B), 3)
})

test_that("predict_B pvalue method accepts vector target_SE", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  # Vector must match length
  B <- predict_B(model, c(0.05, 0.1, 0.2),
                 target_trait = "continuous",
                 target_SE = c(0.5, 0.6, 0.7))
  expect_equal(length(B), 3)

  # Wrong length should error
  expect_error(
    predict_B(model, c(0.05, 0.1, 0.2),
              target_trait = "continuous",
              target_SE = c(0.5, 0.6)),
    "same length"
  )
})

test_that("predict_B pvalue method accepts scalar target_case_prop", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  B <- predict_B(model, c(0.05, 0.1, 0.2),
                 target_trait = "binary",
                 target_case_prop = 0.3)
  expect_equal(length(B), 3)
})

test_that("predict_B pvalue method accepts vector target_case_prop", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  B <- predict_B(model, c(0.05, 0.1, 0.2),
                 target_trait = "binary",
                 target_case_prop = c(0.3, 0.3, 0.3))
  expect_equal(length(B), 3)

  # Wrong length should error
  expect_error(
    predict_B(model, c(0.05, 0.1, 0.2),
              target_trait = "binary",
              target_case_prop = c(0.3, 0.3)),
    "same length"
  )
})

test_that("predict_B pvalue method errors on invalid target_trait", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  expect_error(
    predict_B(model, c(0.1), target_trait = "invalid"),
    "target_trait"
  )
})

test_that("predict_B pvalue method handles zero-MAF entries with conversion", {
  model <- .make_test_glow_B_model(method = "pvalue_method")

  B <- predict_B(model, c(0, 0.05, 0.1, 0),
                 target_trait = "binary",
                 target_case_prop = 0.3)
  expect_equal(B[1], 0)
  expect_equal(B[4], 0)
  expect_true(all(B >= 0))
})

test_that("predict_B binary conversion scales with case proportion", {
  model <- .make_test_glow_B_model(method = "pvalue_method")
  maf <- c(0.05, 0.1, 0.2)

  # Lower case prop -> larger B (more extreme imbalance)
  B_low <- predict_B(model, maf, target_trait = "binary",
                      target_case_prop = 0.1)
  B_balanced <- predict_B(model, maf, target_trait = "binary",
                           target_case_prop = 0.5)

  # B = h / sqrt(p * (1-p)), so balanced (p=0.5) gives smaller B
  expect_true(all(B_low > B_balanced))
})

test_that("predict_B continuous conversion scales with SE", {
  model <- .make_test_glow_B_model(method = "pvalue_method")
  maf <- c(0.05, 0.1, 0.2)

  # B = h * SE, so larger SE gives larger B
  B_small_se <- predict_B(model, maf, target_trait = "continuous",
                           target_SE = 0.1)
  B_large_se <- predict_B(model, maf, target_trait = "continuous",
                           target_SE = 1.0)

  expect_true(all(B_large_se > B_small_se))
})


# ========== predict_B: beta method ignores target params ==========

test_that("predict_B beta method ignores target_trait and related params", {
  model <- .make_test_glow_B_model(method = "beta_method")
  maf <- c(0.05, 0.1, 0.2)

  # Beta method: target params should not affect result
  B_no_params <- predict_B(model, maf)
  B_with_params <- predict_B(model, maf, target_trait = "binary",
                              target_case_prop = 0.3)
  expect_equal(B_no_params, B_with_params)
})


# ========== save_B_model / load_B_model tests ==========

test_that("save_B_model accepts glow_B_model directly and round-trips", {
  model <- .make_test_glow_B_model()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  # save_B_model now accepts glow_B_model as primary input
  save_B_model(model, tmp)
  loaded <- load_B_model(tmp)

  expect_s3_class(loaded, "glow_B_model")
  expect_equal(loaded$method_used, model$method_used)
})

test_that("save_B_model accepts glow_B_estimate and extracts $model", {
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)
  logX <- log(maf)
  Y <- -0.14 * logX + 0.5 + rnorm(50, 0, 0.02)
  lm_obj <- lm(Y ~ logX, data = data.frame(Y = Y, logX = logX))
  attr(lm_obj, "model_id") <- 3L

  inner_model <- structure(list(
    method_used = "beta_method",
    models = list(beta_method = lm_obj, pvalue_method = NULL),
    all_models_info = list(beta_method = NULL, pvalue_method = NULL),
    outliers = list(method = "none", action = "flag",
                    beta_method = NULL, pvalue_method = NULL,
                    indices_removed = integer(0)),
    training_summary = list(n_original = 50L, n_used = 50L,
                             n_outliers_detected = 0L,
                             trait_type = "binary"),
    training_data = list(MAF = maf, BETA = NULL, P = NULL,
                         P_mlog10 = NULL, N = NULL, trait = "binary"),
    comparison = list(criterion_beta_method = 0.85, method_selected = "beta_method"),
    selection_criterion = "R2"
  ), class = "glow_B_model")

  estimate <- structure(list(
    B = sqrt(abs(predict(lm_obj))),
    B_beta_method = sqrt(abs(predict(lm_obj))),
    B_pvalue_method = NULL,
    model = inner_model,
    target_summary = list(n_predictions = 50L, trait_type = "binary",
                           MAF_range = c(0.01, 0.3))
  ), class = "glow_B_estimate")

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  # Should save the inner glow_B_model, not the full estimate
  save_B_model(estimate, tmp)
  loaded <- load_B_model(tmp)

  expect_s3_class(loaded, "glow_B_model")
  expect_equal(loaded$method_used, "beta_method")
})

test_that("save -> load -> predict gives identical results", {
  model <- .make_test_glow_B_model()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  target_maf <- c(0.01, 0.05, 0.1, 0.2, 0.3)
  B_before <- predict_B(model, target_maf)

  save_B_model(model, tmp)
  loaded <- load_B_model(tmp)
  B_after <- predict_B(loaded, target_maf)

  expect_equal(B_before, B_after)
})

test_that("save_B_model rejects invalid input", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  expect_error(save_B_model("not_a_model", tmp), "glow_B_model.*glow_B_estimate")
})

test_that("load_B_model validates structure", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  # Save a non-model object
  saveRDS(list(foo = "bar"), tmp)
  expect_error(load_B_model(tmp), "model|metadata|format")
})

test_that("load_B_model errors on nonexistent file", {
  expect_error(load_B_model("/nonexistent/path.rds"), "not found")
})

test_that("load_B_model returns glow_B_model (not glow_B_estimate)", {
  model <- .make_test_glow_B_model()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  save_B_model(model, tmp)
  loaded <- load_B_model(tmp)

  # Must return glow_B_model, not glow_B_estimate
  expect_s3_class(loaded, "glow_B_model")
  expect_false(inherits(loaded, "glow_B_estimate"))
})

test_that("load_B_model handles legacy glow_B_estimate envelope", {
  # Simulate old save format where envelope$model was a glow_B_estimate
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)
  logX <- log(maf)
  Y <- -0.14 * logX + 0.5 + rnorm(50, 0, 0.02)
  lm_obj <- lm(Y ~ logX, data = data.frame(Y = Y, logX = logX))
  attr(lm_obj, "model_id") <- 3L

  inner_model <- structure(list(
    method_used = "beta_method",
    models = list(beta_method = lm_obj, pvalue_method = NULL),
    all_models_info = list(beta_method = NULL, pvalue_method = NULL),
    outliers = list(method = "none", action = "flag",
                    beta_method = NULL, pvalue_method = NULL,
                    indices_removed = integer(0)),
    training_summary = list(n_original = 50L, n_used = 50L,
                             n_outliers_detected = 0L,
                             trait_type = "binary"),
    training_data = list(MAF = maf, BETA = NULL, P = NULL,
                         P_mlog10 = NULL, N = NULL, trait = "binary"),
    comparison = list(criterion_beta_method = 0.85, method_selected = "beta_method"),
    selection_criterion = "R2"
  ), class = "glow_B_model")

  old_estimate <- structure(list(
    B = sqrt(abs(predict(lm_obj))),
    B_beta_method = sqrt(abs(predict(lm_obj))),
    B_pvalue_method = NULL,
    model = inner_model,
    target_summary = list(n_predictions = 50L, trait_type = "binary",
                           MAF_range = c(0.01, 0.3))
  ), class = "glow_B_estimate")

  # Old envelope format: model field is glow_B_estimate
  envelope <- list(
    model = old_estimate,
    metadata = list(
      glowr_version = as.character(utils::packageVersion("GLOWr")),
      save_date = Sys.Date(),
      training_n = 50L,
      best_method = "beta_method",
      best_model_id = 3L,
      R2 = 0.85
    ),
    format_version = 1L
  )

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(envelope, tmp)

  # load_B_model should extract the inner glow_B_model from legacy format
  loaded <- load_B_model(tmp)
  expect_s3_class(loaded, "glow_B_model")
  expect_equal(loaded$method_used, "beta_method")
})

test_that("load_B_model warns on version mismatch", {
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)
  logX <- log(maf)
  Y <- -0.14 * logX + 0.5 + rnorm(50, 0, 0.02)
  lm_obj <- lm(Y ~ logX, data = data.frame(Y = Y, logX = logX))
  attr(lm_obj, "model_id") <- 3L

  b_model <- structure(list(
    method_used = "beta_method",
    models = list(beta_method = lm_obj, pvalue_method = NULL),
    all_models_info = list(beta_method = NULL, pvalue_method = NULL),
    outliers = list(method = "none", action = "flag",
                    beta_method = NULL, pvalue_method = NULL,
                    indices_removed = integer(0)),
    training_summary = list(n_original = 50L, n_used = 50L,
                             n_outliers_detected = 0L,
                             trait_type = "binary"),
    training_data = list(MAF = maf, BETA = NULL, P = NULL,
                         P_mlog10 = NULL, N = NULL, trait = "binary"),
    comparison = list(criterion_beta_method = 0.85, method_selected = "beta_method"),
    selection_criterion = "R2"
  ), class = "glow_B_model")

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  # Save with fake old version
  envelope <- list(
    model = b_model,
    metadata = list(
      glowr_version = "0.0.1",
      save_date = Sys.Date(),
      training_n = 50L,
      best_method = "beta_method",
      best_model_id = 3L,
      R2 = 0.85
    ),
    format_version = 1L
  )
  saveRDS(envelope, tmp)

  expect_warning(load_B_model(tmp), "version")
})

test_that("save_B_model includes correct metadata", {
  model <- .make_test_glow_B_model()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))

  save_B_model(model, tmp)
  envelope <- readRDS(tmp)

  expect_true("metadata" %in% names(envelope))
  expect_true("format_version" %in% names(envelope))
  expect_equal(envelope$format_version, 1L)

  # Envelope should contain glow_B_model (not glow_B_estimate)
  expect_s3_class(envelope$model, "glow_B_model")

  meta <- envelope$metadata
  expect_equal(meta$best_method, "beta_method")
  expect_equal(meta$training_n, 50L)
  expect_true(!is.null(meta$glowr_version))
  expect_true(!is.null(meta$save_date))
  expect_true(!is.null(meta$R2))
})


# ========== Integration: model_id attribute in .fit_all_candidate_models ==========

test_that("model_id attribute is attached by .fit_all_candidate_models", {
  set.seed(42)
  X <- runif(50, 0.01, 0.3)
  Y <- X^2 + 0.1 * rnorm(50)
  Y <- abs(Y) + 0.01  # Ensure positive

  all_models <- GLOWr:::.fit_all_candidate_models(X, Y, verbose = 0)

  # Check that each model has a model_id attribute
  for (i in seq_along(all_models)) {
    mid <- attr(all_models[[i]]$model, "model_id")
    expect_true(!is.null(mid),
                info = paste("Model", names(all_models)[i], "missing model_id"))
    expect_true(mid >= 1L && mid <= 8L,
                info = paste("Model", names(all_models)[i], "has invalid model_id:", mid))
  }
})

test_that("predict_B works end-to-end with select_best_model", {
  set.seed(42)
  X <- runif(50, 0.01, 0.3)
  Y <- 0.5 * log(X * (1 - X)) + 2 + rnorm(50, 0, 0.1)
  Y <- abs(Y) + 0.01

  best <- select_best_model(X, Y, verbose = FALSE)

  # select_best_model should now return lm with model_id
  mid <- attr(best, "model_id")
  expect_true(!is.null(mid), info = "select_best_model lm should have model_id")

  # predict_B should work with this lm object
  B <- predict_B(best, c(0.05, 0.1, 0.2))
  expect_true(is.numeric(B))
  expect_equal(length(B), 3)
  expect_true(all(B >= 0))
})


# ========== Integration: predict_B with train_B_model ==========

test_that("predict_B works with train_B_model beta method (integration)", {
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)
  beta <- sqrt(0.5 * maf * (1 - maf)) * 0.1 + rnorm(50, 0, 0.001)

  b_model <- train_B_model(
    training_trait = "binary",
    training_MAF = maf,
    training_BETA = beta,
    method = "beta",
    show_model_selection = FALSE,
    verbose = 0
  )

  B <- predict_B(b_model, target_MAF = c(0.01, 0.05, 0.1, 0.2))
  expect_true(is.numeric(B))
  expect_equal(length(B), 4)
  expect_true(all(B >= 0))
  expect_true(all(is.finite(B)))
})

test_that("predict_B works with train_B_model pvalue method (integration)", {
  set.seed(42)
  maf <- runif(50, 0.01, 0.3)
  p_vals <- runif(50, 1e-8, 0.01)
  n_vals <- rep(5000, 50)

  b_model <- train_B_model(
    training_trait = "binary",
    training_MAF = maf,
    training_P = p_vals,
    training_N = n_vals,
    method = "pvalue",
    show_model_selection = FALSE,
    verbose = 0
  )

  # P-value method with binary trait conversion
  B <- predict_B(b_model, target_MAF = c(0.01, 0.05, 0.1, 0.2),
                 target_trait = "binary",
                 target_case_prop = 0.3)
  expect_true(is.numeric(B))
  expect_equal(length(B), 4)
  expect_true(all(B >= 0))
  expect_true(all(is.finite(B)))
})
