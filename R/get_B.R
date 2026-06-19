########## Effect Size (B) Estimation for Optimal Weights ##########
#
# This file contains get_B(), the primary user-facing function for estimating
# effect sizes based on minor allele frequency (MAF). It is a convenience
# wrapper around train_B_model() + predict_B().
#
# EXPORTED FUNCTIONS:
#   get_B()                  - Train model and predict B in one call
#   print.glow_B_estimate()  - S3 print method for glow_B_estimate objects

#################### EXPORTED MAIN FUNCTIONS ####################

#' Estimate Effect Sizes (B) for Variants Based on MAF
#'
#' @description
#' Trains the best regression model from training data and predicts allelic
#' effect sizes (B) for target MAF values in one call. This is a convenience
#' wrapper around \code{\link{train_B_model}} (training) and
#' \code{\link{predict_B}} (prediction).
#'
#' Two estimation methods are supported: (1) direct beta method using effect
#' sizes, and (2) p-value/Z-score method using trait-independent h-squared.
#' When both are applicable, the function compares them and selects the
#' better-fitting model as the primary result.
#'
#' For workflows that train once and predict many times (e.g., per-chromosome
#' analysis), use \code{\link{train_B_model}} and \code{\link{predict_B}}
#' directly instead.
#'
#' @param training_trait Character string or NULL, trait type of training data.
#'   Can be "binary", "continuous", "mixed" (when training data contains both
#'   trait types), or NULL (unspecified). When "mixed" or NULL, only the
#'   p-value/Z-score method is available since the direct beta method requires
#'   a single trait type
#' @param training_MAF Numeric vector of minor allele frequencies from training
#'   data, used to fit the B prediction model. Values should be in (0, 0.5].
#'   Values > 0.5 are automatically folded to 1 - MAF with a warning
#' @param training_BETA Numeric vector of effect sizes (regression coefficients)
#'   from training data. Required for direct beta method. Should have the same
#'   length as training_MAF
#' @param training_P Numeric vector of p-values from training data. Used for
#'   p-value/Z-score method. Values should be in (0, 1). Either training_P OR
#'   training_P_mlog10 should be provided (not both) for p-value method
#' @param training_P_mlog10 Numeric vector of -log10 transformed p-values from
#'   training data. Use this instead of training_P when dealing with extremely
#'   small p-values (e.g., from GWAS) to avoid numerical underflow. For example,
#'   if p = 1e-300, use training_P_mlog10 = 300. Either training_P OR
#'   training_P_mlog10 should be provided (not both) for p-value method
#' @param training_N Numeric vector of sample sizes from training data. Required
#'   for p-value/Z-score method
#' @param target_trait Character string or NULL, trait type of target data.
#'   Must be either "binary", "continuous", or NULL (default). For the beta
#'   method, target_trait can be omitted if you are confident that training
#'   and target data share the same trait type. For the p-value method,
#'   target_trait is required. When NULL, appropriate warnings will be issued
#'   to ensure scientific validity
#' @param target_MAF Numeric vector of minor allele frequencies from target
#'   data for which to predict effect sizes. Values should be in (0, 0.5].
#'   Values > 0.5 are automatically folded to 1 - MAF with a warning
#' @param target_SE Numeric vector of standard errors from target data.
#'   Required when target_trait == "continuous" for p-value/Z-score method
#' @param target_case_prop Numeric value or vector specifying the proportion of
#'   cases in the target case-control study. This is calculated as:
#'   number_of_cases / (number_of_cases + number_of_controls).
#'   For example, in a balanced study with 500 cases and 500 controls,
#'   target_case_prop = 500/1000 = 0.5.
#'   IMPORTANT: This is NOT the population disease prevalence. For a rare disease
#'   with 0.1\% population prevalence, a typical balanced case-control study would
#'   have target_case_prop = 0.5, not 0.001.
#'   Required when target_trait == "binary" for p-value/Z-score method.
#'   Values must be in (0, 1)
#' @param method Character string specifying which method to use: "auto" (default),
#'   "beta", "pvalue", or "both". "auto" uses both methods when same trait and
#'   beta available, otherwise uses pvalue method. "both" runs both methods for
#'   comparison when both are applicable
#' @param selection_criterion Character string specifying model selection
#'   criterion: "R2" (default), "adj_R2", or "CV_R2". See
#'   \code{\link{select_best_model}} for detailed guidance. Note: "CV" is
#'   deprecated, use "CV_R2" instead
#' @param custom_models Optional list of custom model formulas to include in
#'   candidate model set. Each element should be a formula object (e.g.,
#'   \code{formula(Y ~ poly(X, 2))}). Custom models are fit alongside the 8
#'   standard models
#' @param outlier_method Character string specifying outlier detection method:
#'   "none" (default, no detection), "statistical" (Cook's distance),
#'   "biological" (implausible MAF-effect combinations), or "both"
#' @param outlier_action Character string specifying action for detected outliers:
#'   "flag" (default, report but don't remove) or "remove" (exclude from analysis)
#' @param cook_threshold Numeric threshold multiplier for Cook's distance outlier
#'   detection (default: 4). Outliers are defined as points with
#'   Cook's D > cook_threshold/n
#' @param biological_rules Optional list of custom rules for biological outlier
#'   detection. Default uses common variants (MAF > 0.05) with large effects
#'   (|effect| > 10)
#' @param show_model_selection Logical. If TRUE (default), prints the model
#'   selection comparison table when estimating B. Set to FALSE to suppress
#'   output in loops or simulations
#' @param return_full Logical. If FALSE (default), returns only the primary B
#'   estimates as a numeric vector. If TRUE, returns a comprehensive result
#'   object of class \code{glow_B_estimate}
#' @param verbose Integer controlling verbosity: 0=silent, 1=warnings (default),
#'   2=info messages, 3=debug output
#'
#' @return If \code{return_full = FALSE} (default): Numeric vector of predicted
#'   effect sizes (B) for each variant in \code{target_MAF}.
#'
#'   If \code{return_full = TRUE}: An object of class \code{glow_B_estimate}
#'   containing:
#'   \describe{
#'     \item{B}{Primary B predictions (numeric vector)}
#'     \item{B_beta_method}{B predictions from beta method (NULL if not used)}
#'     \item{B_pvalue_method}{B predictions from p-value method (NULL if not
#'       used)}
#'     \item{model}{A \code{glow_B_model} object (from
#'       \code{\link{train_B_model}}) containing all training results: fitted
#'       models, outlier info, training data, and training-time comparison}
#'     \item{target_summary}{List with n_predictions, trait_type, MAF_range}
#'   }
#'
#'   When both methods were used, the \code{model$comparison} field is
#'   augmented with prediction-based comparison statistics (correlation,
#'   RMSE, percent difference).
#'
#' @details
#' This function is a convenience wrapper equivalent to:
#' \preformatted{
#'   b_model <- train_B_model(...)
#'   B <- predict_B(b_model, target_MAF, ...)
#' }
#'
#' For "train once, predict many times" workflows, call
#' \code{\link{train_B_model}} and \code{\link{predict_B}} separately.
#'
#' \strong{Two Estimation Methods:}
#'
#' \strong{Method 1: Direct Beta Method}
#' \itemize{
#'   \item Uses training_BETA values directly
#'   \item Fits regression model: \eqn{BETA^2 \sim f(MAF)}
#'   \item Best when training and target have the same trait type
#'   \item Provides most direct relationship between MAF and effect size
#' }
#'
#' \strong{Method 2: P-value/Z-score Method}
#' \itemize{
#'   \item Converts p-values to Z-squared using inverse chi-squared CDF
#'   \item Calculates trait-independent function:
#'     \eqn{h^2(q) = Z^2 / (2 * N * q * (1-q))}
#'   \item Fits regression model: \eqn{h^2 \sim f(MAF)}
#'   \item Converts back to trait-specific B:
#'     \itemize{
#'       \item Continuous: \eqn{B = h(q) * SE}
#'       \item Binary: \eqn{B = h(q) / \sqrt{p_0 * (1-p_0)}}
#'       where \eqn{p_0} is the case proportion in the case-control study
#'     }
#'   \item Can transfer information between different trait types
#' }
#'
#' \strong{Method Selection Logic:}
#' \itemize{
#'   \item method="auto" (default): Uses "both" if same trait + BETA available,
#'     else "pvalue"
#'   \item method="beta": Only direct beta method (requires training_BETA)
#'   \item method="pvalue": Only p-value/Z-score method (requires
#'     training_P/training_N)
#'   \item method="both": Runs both methods for comparison
#' }
#'
#' \strong{Primary Method Selection (when method="both"):}
#'
#' When both methods are run, the function selects the better method
#' as primary based on model quality. The method with the higher
#' selection_criterion value (R2/adj_R2/CV_R2) is chosen. The comparison
#' is done on the training data models (BETA^2 vs h^2).
#'
#' \strong{Outlier Detection:}
#'
#' When \code{outlier_method} is not "none", outliers are detected using
#' statistical (Cook's distance), biological (implausible combinations),
#' or both methods. The \code{outlier_action} parameter controls whether
#' outliers are flagged or removed.
#'
#' @section Computational Complexity:
#' O(8 * n_train + n_test) where:
#' \itemize{
#'   \item n_train = length(training_MAF) (fitting 8 models)
#'   \item n_test = length(target_MAF) (prediction step)
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Direct beta method (same trait type)
#' set.seed(123)
#' training_MAF <- runif(50, 0.001, 0.3)
#' training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF)) * 0.1
#' target_MAF <- runif(10, 0.001, 0.5)
#'
#' B_estimates <- get_B(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   method = "beta"
#' )
#'
#' # Example 2: P-value/Z-score method
#' training_P <- runif(50, 0.0001, 0.1)
#' training_N <- rep(1000, 50)
#'
#' B_pval <- get_B(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_P = training_P,
#'   training_N = training_N,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   target_case_prop = 0.1,
#'   method = "pvalue"
#' )
#'
#' # Example 3: Full result with both methods
#' result_full <- get_B(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   training_P = training_P,
#'   training_N = training_N,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   target_case_prop = 0.1,
#'   method = "both",
#'   return_full = TRUE,
#'   verbose = 2
#' )
#' print(result_full)
#' print(result_full$model)  # nested glow_B_model
#'
#' # Example 4: Transfer across trait types (continuous to binary)
#' training_P_cont <- runif(50, 0.0001, 0.1)
#' training_N_cont <- rep(2000, 50)
#'
#' B_cross_trait <- get_B(
#'   training_trait = "continuous",
#'   training_MAF = training_MAF,
#'   training_P = training_P_cont,
#'   training_N = training_N_cont,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   target_case_prop = 0.5,
#'   method = "pvalue"  # Only pvalue method works for different traits
#' )
#'
#' # Example 5: Suppress model selection output (useful in loops)
#' for (i in 1:5) {
#'   B <- get_B(
#'     training_trait = "binary",
#'     training_MAF = training_MAF,
#'     training_BETA = training_BETA,
#'     target_trait = "binary",
#'     target_MAF = target_MAF,
#'     show_model_selection = FALSE  # Suppress table in loop
#'   )
#' }
#'
#' # Example 6: Extremely small p-values using -log10 transformation
#' # This is useful for GWAS data with very significant hits
#' training_P_mlog10 <- runif(50, 10, 350)  # -log10(p) from 10 to 350
#' # Equivalent to p-values from 1e-10 to 1e-350
#' training_N_gwas <- rep(50000, 50)  # Large GWAS sample size
#'
#' B_extreme_pval <- get_B(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_P_mlog10 = training_P_mlog10,  # Use -log10(p) for numerical stability
#'   training_N = training_N_gwas,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   target_case_prop = 0.5,
#'   method = "pvalue"
#' )
#'
#' # Example 7: Use in optimal weight calculation
#' PI_estimates <- runif(10, 0.001, 0.1)
#' M <- diag(10)  # Independent variants
#' g_identity <- function(x) x
#'
#' weights <- Optimal_Weights_M(g_identity, B_estimates, PI_estimates, M)
#' }
#'
#' @section References:
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @seealso
#' \code{\link{train_B_model}} for training models separately
#' \code{\link{predict_B}} for applying trained models to new data
#' \code{\link{select_best_model}} for the model selection algorithm
#' \code{\link{Optimal_Weights_M}} for calculating optimal weights using B
#'
#' @export
get_B <- function(training_trait,
                  training_MAF,
                  training_BETA = NULL,
                  training_P = NULL,
                  training_P_mlog10 = NULL,
                  training_N = NULL,
                  target_trait = NULL,
                  target_MAF,
                  target_SE = NULL,
                  target_case_prop = NULL,
                  method = "auto",
                  selection_criterion = "R2",
                  custom_models = NULL,
                  outlier_method = "none",
                  outlier_action = "flag",
                  cook_threshold = 4,
                  biological_rules = NULL,
                  show_model_selection = TRUE,
                  return_full = FALSE,
                  verbose = 1) {

  # ========== Validate Target-Specific Parameters ==========
  # (target_MAF, target_trait, target_SE, target_case_prop are not passed to
  #  train_B_model, so we validate them here before proceeding)

  if (is.null(target_MAF)) {
    stop("target_MAF is required")
  }
  if (!is.numeric(target_MAF) || length(target_MAF) == 0) {
    stop("target_MAF must be a non-empty numeric vector")
  }

  # Fold MAF values > 0.5 to minor allele frequency convention
  # get_B() requires MAF > 0 (unlike predict_B which allows MAF = 0),
  # so use allow_zero = FALSE for consistency with training data requirements
  target_MAF <- .fold_maf(target_MAF, allow_zero = FALSE,
                          context = "target_MAF")

  # Validate target_trait
 valid_target_traits <- c("binary", "continuous")
  if (!is.null(target_trait) && !target_trait %in% valid_target_traits) {
    stop("target_trait must be either 'binary', 'continuous', or NULL, got: ",
         target_trait)
  }

  # ========== Validate Trait Compatibility ==========
  # This check is specific to get_B() -- train_B_model() only knows about
  # training traits. We need to check training vs. target compatibility here.

  # Normalize training_trait for compatibility check
  training_trait_normalized <- if (is.null(training_trait)) "mixed" else training_trait

  # Determine which method will be used (approximate -- train_B_model will

  # do the definitive determination, but we need to know for trait validation)
  beta_method_possible <- training_trait_normalized %in% c("binary", "continuous") &&
                          !is.null(training_BETA)
  has_pvalue_data <- !is.null(training_P) || !is.null(training_P_mlog10)

  if (method == "auto") {
    validation_method <- if (beta_method_possible && has_pvalue_data) {
      "beta"  # "both" mode, validate beta compatibility
    } else if (beta_method_possible) {
      "beta"
    } else {
      "pvalue"
    }
  } else if (method == "both") {
    validation_method <- "beta"
  } else {
    validation_method <- method
  }

  # Only validate if the beta method will be used (pvalue method is flexible)
  if (validation_method == "beta") {
    validation_result <- .validate_trait_compatibility(
      training_trait = training_trait,
      target_trait = target_trait,
      method = "beta"
    )

    if (!validation_result$valid) {
      stop(validation_result$message)
    } else if (validation_result$severity == "warning" && verbose >= 1) {
      warning(validation_result$message, call. = FALSE)
    }
  }

  # ========== Step 1: Train the Model ==========

  b_model <- train_B_model(
    training_trait = training_trait,
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    method = method,
    selection_criterion = selection_criterion,
    custom_models = custom_models,
    outlier_method = outlier_method,
    outlier_action = outlier_action,
    cook_threshold = cook_threshold,
    biological_rules = biological_rules,
    show_model_selection = show_model_selection,
    verbose = verbose
  )

  # ========== Step 2: Predict B for Each Active Method ==========

  B_beta_method <- NULL
  B_pvalue_method <- NULL

  if (b_model$method_used %in% c("beta_method", "both")) {
    B_beta_method <- predict_B(
      model = b_model,
      target_MAF = target_MAF,
      target_trait = target_trait,
      target_SE = target_SE,
      target_case_prop = target_case_prop,
      method = "beta_method"
    )
  }

  if (b_model$method_used %in% c("pvalue_method", "both")) {
    B_pvalue_method <- predict_B(
      model = b_model,
      target_MAF = target_MAF,
      target_trait = target_trait,
      target_SE = target_SE,
      target_case_prop = target_case_prop,
      method = "pvalue_method"
    )
  }

  # ========== Step 3: Determine Primary B Values ==========

  if (b_model$method_used == "both") {
    # Use the method selected by training-time criterion comparison
    method_selected <- b_model$comparison$method_selected
    B_primary <- if (method_selected == "beta_method") {
      B_beta_method
    } else {
      B_pvalue_method
    }
  } else if (b_model$method_used == "beta_method") {
    B_primary <- B_beta_method
  } else {
    # pvalue_method
    B_primary <- B_pvalue_method
  }

  # ========== Step 4: Return Results ==========

  if (return_full) {

    # Build glow_B_estimate with nested glow_B_model
    result <- structure(list(
      B = B_primary,
      B_beta_method = B_beta_method,
      B_pvalue_method = B_pvalue_method,
      model = b_model,
      target_summary = list(
        n_predictions = length(target_MAF),
        trait_type = target_trait,
        MAF_range = range(target_MAF)
      )
    ), class = "glow_B_estimate")

    # Add full prediction-based comparison when both methods were used
    if (b_model$method_used == "both" &&
        !is.null(B_beta_method) && !is.null(B_pvalue_method)) {
      prediction_comparison <- .compare_B_methods(
        B_beta_method, B_pvalue_method, target_MAF
      )

      # Augment the model's training-time comparison with prediction stats
      result$model$comparison$prediction <- prediction_comparison

      if (verbose >= 2) {
        message(sprintf(
          "Prediction comparison: correlation = %.3f, mean %% diff = %.1f%%",
          prediction_comparison$correlation,
          prediction_comparison$mean_percent_diff
        ))
      }
    }

    return(result)

  } else {
    # Simple vector return
    return(B_primary)
  }
}


#' Print Method for glow_B_estimate Objects
#'
#' @description
#' Displays a concise summary of B estimation results, including prediction
#' summary and nested model information. Delegates model details to the
#' \code{\link{print.glow_B_model}} method.
#'
#' @param x A \code{glow_B_estimate} object from \code{\link{get_B}} with
#'   \code{return_full = TRUE}
#' @param ... Additional arguments (ignored)
#'
#' @return The input object, invisibly
#'
#' @examples
#' \dontrun{
#' result <- get_B(
#'   training_trait = "continuous",
#'   training_MAF = runif(50, 0.001, 0.3),
#'   training_BETA = rnorm(50, 0, 0.1),
#'   target_MAF = runif(20, 0.001, 0.5),
#'   return_full = TRUE
#' )
#' print(result)
#' }
#'
#' @export
print.glow_B_estimate <- function(x, ...) {

  cat("\n")
  cat("========================================\n")
  cat("GLOW B Estimation Result (glow_B_estimate)\n")
  cat("========================================\n\n")

  # Prediction summary
  cat("--- Prediction Summary ---\n")
  cat("  N predictions:     ", x$target_summary$n_predictions, "\n")
  cat("  Target trait:      ",
      if (is.null(x$target_summary$trait_type)) "NULL" else
        x$target_summary$trait_type, "\n")
  cat("  Target MAF range:  ",
      sprintf("[%.4f, %.4f]",
              x$target_summary$MAF_range[1],
              x$target_summary$MAF_range[2]), "\n")

  # B value summary
  cat("\n--- B Values ---\n")
  cat("  Primary B range:   ",
      sprintf("[%.4g, %.4g]", min(x$B), max(x$B)), "\n")
  cat("  Primary B mean:    ", sprintf("%.4g", mean(x$B)), "\n")

  if (!is.null(x$B_beta_method)) {
    cat("  Beta method range:  ",
        sprintf("[%.4g, %.4g]",
                min(x$B_beta_method), max(x$B_beta_method)), "\n")
  }
  if (!is.null(x$B_pvalue_method)) {
    cat("  P-value method range: ",
        sprintf("[%.4g, %.4g]",
                min(x$B_pvalue_method), max(x$B_pvalue_method)), "\n")
  }

  # Prediction-based comparison (if available)
  if (!is.null(x$model$comparison$prediction)) {
    pred_comp <- x$model$comparison$prediction
    cat("\n--- Prediction Comparison ---\n")
    cat("  Correlation:       ",
        sprintf("%.4f", pred_comp$correlation), "\n")
    cat("  RMSE:              ",
        sprintf("%.4g", pred_comp$rmse), "\n")
    cat("  Mean %% diff:       ",
        sprintf("%.1f%%", pred_comp$mean_percent_diff), "\n")
  }

  # Delegate model details
  cat("\n--- Nested Model ---\n")
  print(x$model)

  invisible(x)
}
