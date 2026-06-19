########## Diagnostic Functions for B Estimation ##########
#
# This file contains diagnostic and visualization functions for evaluating
# B estimation quality from get_B() and train_B_model(). These tools help
# users assess model fit, detect outliers, and compare estimation methods.
#
# EXPORTED FUNCTIONS:
#   plot.glow_B_estimate()  - S3 plot method for glow_B_estimate objects
#   plot.glow_B_model()     - S3 plot method for glow_B_model objects
#   get_B_diagnostics()     - Comprehensive diagnostic report
#   compare_B_models()      - Compare all 8 candidate models
#   print.glow_B_diagnostics() - S3 print method for diagnostics
#
# INTERNAL HELPER FUNCTIONS:
#   .plot_dashboard()               - 2x2 dashboard of diagnostic plots
#   .plot_model_fit()               - Model fit plot (training data vs fitted)
#   .plot_residuals()               - Residual diagnostic plots
#   .plot_method_comparison()       - Detailed method comparison (4 panels)
#   .plot_method_comparison_simple() - Simple method comparison (1 panel)
#   .plot_all_models()              - All 8 candidate models visualization
#   .generate_diagnostic_summary()  - Text summary generator
#   .plot_model_comparison_table()  - Bar chart of model R2 values
#   .extract_b_model()              - Extract glow_B_model from either class

#################### QUALITY CONTROL THRESHOLDS ####################

# Quality control thresholds for diagnostic flags
QC_OUTLIER_PROPORTION_THRESHOLD <- 0.1     # Flag if >10% of observations are outliers
QC_CORRELATION_THRESHOLD <- 0.7            # Flag if method correlation < 0.7
QC_PERCENT_DIFF_THRESHOLD <- 20            # Flag if mean % difference > 20%
QC_R2_LOW_THRESHOLD <- 0.3                 # Flag if R^2 < 0.3 (poor fit)
QC_R2_MODERATE_THRESHOLD <- 0.5            # Note if R^2 < 0.5 (moderate fit)

#################### EXPORTED MAIN FUNCTIONS ####################

#' Plot Method for B Estimation Results
#'
#' @description
#' S3 plot method for objects of class "glow_B_estimate" returned by
#' \code{\link{get_B}} when \code{return_full = TRUE}. Provides visual
#' diagnostics for model fit quality, residuals, and method comparison.
#'
#' @param x Object of class "glow_B_estimate" from \code{\link{get_B}}
#' @param type Character string specifying plot type:
#'   \itemize{
#'     \item "all" (default): 2x2 dashboard with multiple diagnostic plots
#'     \item "fit": Model fit plot showing data vs fitted curve
#'     \item "residuals": Residual diagnostic plots (Q-Q, residuals vs fitted)
#'     \item "comparison": Compare beta and pvalue methods (if both available)
#'     \item "models": Compare all 8 candidate model fits
#'   }
#' @param which_method Character string, which method to plot: "beta", "pvalue",
#'   or "primary" (default, uses the primary B estimate)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return NULL (invisibly). Plots are displayed as a side effect.
#'
#' @details
#' The function provides several visualization options:
#'
#' \strong{Dashboard View (type = "all"):}
#' \itemize{
#'   \item Top-left: Fitted model vs training data
#'   \item Top-right: Residual Q-Q plot
#'   \item Bottom-left: MAF vs predicted B
#'   \item Bottom-right: Method comparison (if both methods used) or Cook's distance
#' }
#'
#' \strong{Model Fit (type = "fit"):}
#' Shows training data points (MAF vs B^2 or h^2) with fitted curve overlaid.
#'
#' \strong{Residuals (type = "residuals"):}
#' Standard diagnostic plots including Q-Q plot and residuals vs fitted values.
#'
#' \strong{Comparison (type = "comparison"):}
#' Scatterplot comparing beta method vs pvalue method estimates with
#' correlation and agreement statistics.
#'
#' \strong{Computational Complexity:}
#' O(n log n) where n is the number of training observations, due to sorting
#' operations in plotting. For type = "models", complexity is O(8n) as it fits
#' all 8 candidate models.
#'
#' @examples
#' \dontrun{
#' # Load ALS data and run get_B with both methods
#' data(ALS_snvs_B_training)
#'
#' result <- get_B(
#'   training_trait = "binary",
#'   training_MAF = ALS_snvs_B_training$MAF,
#'   training_BETA = ALS_snvs_B_training$BETA,
#'   training_P = ALS_snvs_B_training$P,
#'   training_N = ALS_snvs_B_training$N,
#'   target_trait = "binary",
#'   target_MAF = seq(0.01, 0.4, by = 0.05),
#'   target_case_prop = 0.001,
#'   method = "both",
#'   return_full = TRUE
#' )
#'
#' # Dashboard view
#' plot(result)
#'
#' # Individual plots
#' plot(result, type = "fit")
#' plot(result, type = "residuals")
#' plot(result, type = "comparison")
#' }
#'
#' @seealso
#' \code{\link{get_B}} for B estimation
#' \code{\link{plot.glow_B_model}} for model-level diagnostics
#' \code{\link{get_B_diagnostics}} for detailed diagnostic reports
#' \code{\link{compare_B_models}} for comparing candidate models
#'
#' @export
plot.glow_B_estimate <- function(x, type = "all", which_method = "primary", ...) {

  # Input validation
  if (!inherits(x, "glow_B_estimate")) {
    stop("x must be an object of class 'glow_B_estimate'")
  }

  valid_types <- c("all", "fit", "residuals", "comparison", "models")
  if (!type %in% valid_types) {
    stop("type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Determine which model to use (using nested glow_B_model)
  if (which_method == "primary") {
    method_used <- x$model$method_used
    if (method_used == "both" && !is.null(x$model$comparison$method_selected)) {
      which_method <- sub("_method$", "", x$model$comparison$method_selected)
    } else {
      which_method <- if (method_used == "pvalue_method") "pvalue" else "beta"
    }
  }

  if (!which_method %in% c("beta", "pvalue")) {
    stop("which_method must be 'beta', 'pvalue', or 'primary'")
  }

  # Map short names to model field names
  model_key <- paste0(which_method, "_method")

  # Check if requested method is available (in nested model)
  if (is.null(x$model$models[[model_key]])) {
    stop(which_method, " method results not available. ",
         "Use which_method = '",
         if (which_method == "beta") "pvalue" else "beta",
         "' or method = '", which_method, "' in get_B()")
  }

  # Route to appropriate plotting function
  if (type == "all") {
    .plot_dashboard(x, which_method = which_method, ...)
  } else if (type == "fit") {
    .plot_model_fit(x, which_method = which_method, ...)
  } else if (type == "residuals") {
    .plot_residuals(x, which_method = which_method, ...)
  } else if (type == "comparison") {
    .plot_method_comparison(x, ...)
  } else if (type == "models") {
    .plot_all_models(x, which_method = which_method, ...)
  }

  invisible(NULL)
}


#' Plot Method for B Estimation Models
#'
#' @description
#' S3 plot method for objects of class "glow_B_model" returned by
#' \code{\link{train_B_model}}. Provides model-level diagnostic plots
#' including model fit, residuals, Q-Q plots, and model comparison.
#'
#' Unlike \code{\link{plot.glow_B_estimate}}, this method works without
#' target/prediction data. It focuses on training-time model diagnostics.
#'
#' @param x Object of class "glow_B_model" from \code{\link{train_B_model}}
#' @param type Character string specifying plot type:
#'   \itemize{
#'     \item "all" (default): 2x2 dashboard with model diagnostics
#'     \item "fit": Model fit plot showing training data vs fitted curve
#'     \item "residuals": Standard residual diagnostic plots (Q-Q, residuals vs fitted)
#'     \item "models": Compare all 8 candidate model fits (bar chart)
#'   }
#' @param which_method Character string, which method to plot: "beta", "pvalue",
#'   or "primary" (default, uses the primary/selected method)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return NULL (invisibly). Plots are displayed as a side effect.
#'
#' @details
#' This method produces model-level diagnostics only (no prediction-based
#' plots, since a \code{glow_B_model} does not contain B predictions).
#'
#' \strong{Dashboard View (type = "all"):}
#' \itemize{
#'   \item Top-left: Fitted model vs training data
#'   \item Top-right: Normal Q-Q plot of residuals
#'   \item Bottom-left: Residuals vs fitted values
#'   \item Bottom-right: Cook's distance
#' }
#'
#' \strong{Comparison Note:}
#' The "comparison" type from \code{\link{plot.glow_B_estimate}} is not
#' available here because it requires predicted B vectors. Use
#' \code{type = "models"} to compare all candidate model fits instead.
#'
#' @examples
#' \dontrun{
#' b_model <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   method = "beta"
#' )
#'
#' # Dashboard view
#' plot(b_model)
#'
#' # Individual plots
#' plot(b_model, type = "fit")
#' plot(b_model, type = "residuals")
#' plot(b_model, type = "models")
#' }
#'
#' @seealso
#' \code{\link{train_B_model}} for training B estimation models
#' \code{\link{plot.glow_B_estimate}} for prediction-level diagnostics
#'
#' @export
plot.glow_B_model <- function(x, type = "all", which_method = "primary", ...) {

  # Input validation
  if (!inherits(x, "glow_B_model")) {
    stop("x must be an object of class 'glow_B_model'")
  }

  valid_types <- c("all", "fit", "residuals", "models")
  if (!type %in% valid_types) {
    stop("type must be one of: ", paste(valid_types, collapse = ", "),
         "\nNote: 'comparison' requires a glow_B_estimate (predicted B vectors)")
  }

  # Determine which model to use
  if (which_method == "primary") {
    method_used <- x$method_used
    which_method <- if (method_used == "pvalue_method") "pvalue" else "beta"
    # If "both", use the selected method
    if (method_used == "both" && !is.null(x$comparison$method_selected)) {
      which_method <- sub("_method$", "", x$comparison$method_selected)
    }
  }

  if (!which_method %in% c("beta", "pvalue")) {
    stop("which_method must be 'beta', 'pvalue', or 'primary'")
  }

  # Map short names to model field names
  model_key <- paste0(which_method, "_method")

  # Check if requested method is available
  if (is.null(x$models[[model_key]])) {
    stop(which_method, " method model not available. ",
         "Use which_method = '",
         if (which_method == "beta") "pvalue" else "beta", "'")
  }

  # Create a lightweight wrapper for the internal plotting functions.
  # These functions expect an object with $model$models structure.
  # We wrap the glow_B_model to match.
  wrapper <- list(model = x)

  # Route to appropriate plotting function
  if (type == "all") {
    .plot_dashboard_model(wrapper, which_method = which_method, ...)
  } else if (type == "fit") {
    .plot_model_fit(wrapper, which_method = which_method, ...)
  } else if (type == "residuals") {
    .plot_residuals(wrapper, which_method = which_method, ...)
  } else if (type == "models") {
    .plot_all_models(wrapper, which_method = which_method, ...)
  }

  invisible(NULL)
}


#' Generate Comprehensive Diagnostics for B Estimation
#'
#' @description
#' Generates a comprehensive diagnostic report for B estimation results,
#' including model fit statistics, outlier detection, method comparison
#' (if applicable), and quality control flags.
#'
#' @param B_result Object of class \code{"glow_B_estimate"} (from
#'   \code{\link{get_B}} with \code{return_full = TRUE}) or
#'   \code{"glow_B_model"} (from \code{\link{train_B_model}}).
#'   When a \code{glow_B_model} is provided, prediction-based diagnostics
#'   (method comparison) are skipped.
#' @param show_plots Logical. If TRUE (default), displays diagnostic plots
#'   using \code{\link{plot.glow_B_estimate}} or \code{\link{plot.glow_B_model}}
#' @param outlier_threshold Numeric. Cook's distance threshold multiplier for
#'   outlier detection. Outliers defined as Cook's D > outlier_threshold/n where
#'   n is sample size (default = 4, following standard statistical convention)
#' @param verbose Integer. Verbosity level: 0=silent, 1=summary (default),
#'   2=detailed output
#'
#' @return A list of class "glow_B_diagnostics" containing:
#' \itemize{
#'   \item model_fit: List with R^2, adj_R2, formula for each method
#'   \item residual_stats: Summary of residuals (mean, SD, range, normality test)
#'   \item outliers: Indices and details of potential outliers
#'   \item method_comparison: Comparison metrics if both methods used
#'   \item qc_flags: Quality control warnings and recommendations
#'   \item summary_text: Human-readable summary
#' }
#'
#' @details
#' This function provides a comprehensive assessment of B estimation quality:
#'
#' \strong{Model Fit Statistics:}
#' \itemize{
#'   \item R^2: Proportion of variance explained (higher is better, >0.5 preferred)
#'   \item Adjusted R^2: R^2 adjusted for number of predictors
#'   \item Formula: The selected model formula
#' }
#'
#' \strong{Outlier Detection:}
#' Uses Cook's distance to identify influential observations that may affect
#' model fit. Threshold is typically 4/n where n is sample size.
#'
#' \strong{Method Comparison (if both methods used):}
#' \itemize{
#'   \item Correlation between methods
#'   \item RMSE (root mean squared error)
#'   \item Percent difference statistics
#' }
#'
#' \strong{QC Flags:}
#' Automatic warnings for:
#' \itemize{
#'   \item Low R^2 (<0.3)
#'   \item Non-normal residuals
#'   \item High number of outliers (>10\%)
#'   \item Large method disagreement (if both methods used)
#' }
#'
#' @section Computational Complexity:
#' O(n) where n is the number of training variants
#'
#' @examples
#' \dontrun{
#' # Run get_B with full output
#' result <- get_B(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   training_P = training_P,
#'   training_N = training_N,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   target_case_prop = 0.1,
#'   method = "both",
#'   return_full = TRUE
#' )
#'
#' # Generate diagnostics
#' diag <- get_B_diagnostics(result, show_plots = TRUE)
#'
#' # View summary
#' print(diag)
#'
#' # Access specific components
#' diag$model_fit$beta$R2
#' diag$outliers$indices
#' diag$qc_flags
#' }
#'
#' @seealso
#' \code{\link{get_B}} for B estimation
#' \code{\link{train_B_model}} for training B estimation models
#' \code{\link{plot.glow_B_estimate}} for visualization
#' \code{\link{compare_B_models}} for comparing all candidate models
#'
#' @export
get_B_diagnostics <- function(B_result,
                               show_plots = TRUE,
                               outlier_threshold = 4,
                               verbose = 1) {

  # Extract glow_B_model, keeping track of whether we have predictions
  b_model <- .extract_b_model(B_result)
  has_predictions <- inherits(B_result, "glow_B_estimate")

  # Initialize diagnostics list
  diagnostics <- list(
    model_fit = list(),
    residual_stats = list(),
    outliers = list(),
    method_comparison = NULL,
    qc_flags = character(),
    summary_text = ""
  )

  # ========== Model Fit Statistics ==========

  # Beta method
  if (!is.null(b_model$models$beta_method)) {
    model_beta <- b_model$models$beta_method
    summ_beta <- summary(model_beta)  # Cache summary
    diagnostics$model_fit$beta <- list(
      R2 = summ_beta$r.squared,
      adj_R2 = summ_beta$adj.r.squared,
      formula = format(formula(model_beta)),
      n_obs = length(residuals(model_beta))
    )
  }

  # P-value method
  if (!is.null(b_model$models$pvalue_method)) {
    model_pvalue <- b_model$models$pvalue_method
    summ_pvalue <- summary(model_pvalue)  # Cache summary
    diagnostics$model_fit$pvalue <- list(
      R2 = summ_pvalue$r.squared,
      adj_R2 = summ_pvalue$adj.r.squared,
      formula = format(formula(model_pvalue)),
      n_obs = length(residuals(model_pvalue))
    )
  }

  # ========== Residual Statistics ==========

  # Choose primary model for residual analysis
  primary_model <- if (!is.null(b_model$models$beta_method)) {
    b_model$models$beta_method
  } else {
    b_model$models$pvalue_method
  }

  if (!is.null(primary_model)) {
    resids <- residuals(primary_model)

    # Shapiro-Wilk test for normality (if n <= 5000)
    shapiro_test_result <- if (length(resids) <= 5000) {
      shapiro.test(resids)
    } else {
      list(p.value = NA, method = "Sample too large for Shapiro-Wilk test")
    }

    diagnostics$residual_stats <- list(
      mean = mean(resids),
      sd = sd(resids),
      min = min(resids),
      max = max(resids),
      q1 = quantile(resids, 0.25),
      median = median(resids),
      q3 = quantile(resids, 0.75),
      normality_p = shapiro_test_result$p.value,
      normality_test = shapiro_test_result$method
    )

    # Flag non-normal residuals
    if (!is.na(shapiro_test_result$p.value) && shapiro_test_result$p.value < 0.05) {
      diagnostics$qc_flags <- c(diagnostics$qc_flags,
                                "Residuals deviate from normality (Shapiro-Wilk p < 0.05)")
    }
  }

  # ========== Outlier Detection ==========

  if (!is.null(primary_model)) {
    cooks_d <- cooks.distance(primary_model)
    n <- length(cooks_d)
    threshold <- outlier_threshold / n

    outlier_indices <- which(cooks_d > threshold)

    diagnostics$outliers <- list(
      n_outliers = length(outlier_indices),
      proportion = length(outlier_indices) / n,
      indices = outlier_indices,
      cooks_d = if (length(outlier_indices) > 0) cooks_d[outlier_indices] else numeric(0),
      threshold = threshold
    )

    # Flag high proportion of outliers
    if (diagnostics$outliers$proportion > QC_OUTLIER_PROPORTION_THRESHOLD) {
      diagnostics$qc_flags <- c(diagnostics$qc_flags,
                                sprintf("High proportion of outliers detected: %.1f%%",
                                       diagnostics$outliers$proportion * 100))
    }
  }

  # ========== Method Comparison ==========
  # Only available when we have prediction-based comparison (glow_B_estimate)

  if (has_predictions &&
      b_model$method_used == "both" &&
      !is.null(b_model$comparison$prediction)) {
    # Use prediction-based comparison stats
    pred_comp <- b_model$comparison$prediction
    diagnostics$method_comparison <- pred_comp

    # Flag large method disagreement
    if (pred_comp$correlation < QC_CORRELATION_THRESHOLD) {
      diagnostics$qc_flags <- c(diagnostics$qc_flags,
                                sprintf("Low correlation between methods: %.3f",
                                       pred_comp$correlation))
    }
    if (pred_comp$mean_percent_diff > QC_PERCENT_DIFF_THRESHOLD) {
      diagnostics$qc_flags <- c(diagnostics$qc_flags,
                                sprintf("Large mean difference between methods: %.1f%%",
                                       pred_comp$mean_percent_diff))
    }
  }

  # ========== QC Flags ==========

  # Check R^2 for primary model
  primary_r2 <- if (!is.null(diagnostics$model_fit$beta)) {
    diagnostics$model_fit$beta$R2
  } else if (!is.null(diagnostics$model_fit$pvalue)) {
    diagnostics$model_fit$pvalue$R2
  } else {
    NA
  }

  if (!is.na(primary_r2)) {
    if (primary_r2 < QC_R2_LOW_THRESHOLD) {
      diagnostics$qc_flags <- c(diagnostics$qc_flags,
                                sprintf("Low R^2 (%.3f): Model explains <%.0f%% of variance",
                                       primary_r2, QC_R2_LOW_THRESHOLD * 100))
    } else if (primary_r2 < QC_R2_MODERATE_THRESHOLD) {
      diagnostics$qc_flags <- c(diagnostics$qc_flags,
                                sprintf("Moderate R^2 (%.3f): Consider if additional data might improve fit",
                                       primary_r2))
    }
  }

  # ========== Generate Summary Text ==========

  diagnostics$summary_text <- .generate_diagnostic_summary(diagnostics, B_result)

  # ========== Display Plots if Requested ==========

  if (show_plots) {
    plot(B_result, type = "all")
  }

  # ========== Print Summary if Verbose ==========

  if (verbose >= 1) {
    cat(diagnostics$summary_text)
  }

  # Set class and return
  class(diagnostics) <- "glow_B_diagnostics"
  return(diagnostics)
}


#' Compare All Candidate Models for B Estimation
#'
#' @description
#' Fits and compares all 8 candidate regression models for B estimation,
#' returning a summary table with fit statistics. This helps users understand
#' which transformations best capture the MAF-effect size relationship.
#'
#' @param B_result Object of class \code{"glow_B_estimate"} (from
#'   \code{\link{get_B}} with \code{return_full = TRUE}) or
#'   \code{"glow_B_model"} (from \code{\link{train_B_model}}).
#' @param which_method Character string specifying which method to analyze:
#'   "beta", "pvalue", or "both" (default)
#' @param criterion Character string specifying selection criterion to highlight:
#'   "R2" (default), "adj_R2", or "CV_R2"
#' @param plot Logical. If TRUE (default), creates a visualization comparing models
#' @param verbose Integer. Verbosity level: 0=silent, 1=warnings (default),
#'   2=info messages
#'
#' @return A data frame with one row per model containing:
#' \itemize{
#'   \item model_name: Descriptive name (e.g., "Y ~ X", "log(Y) ~ log(X)")
#'   \item formula: Model formula
#'   \item R2: R-squared value
#'   \item adj_R2: Adjusted R-squared
#'   \item CV_R2: Cross-validated R-squared (if requested)
#'   \item is_best: Logical indicating if this is the selected model (by criterion)
#'   \item method: Which estimation method ("beta_method" or "pvalue_method")
#' }
#'
#' @details
#' The function fits all 8 candidate models:
#' \enumerate{
#'   \item Y ~ X
#'   \item Y ~ f(X) where f(X) = X(1-X)
#'   \item Y ~ log(X)
#'   \item Y ~ log(f(X))
#'   \item log(Y) ~ X
#'   \item log(Y) ~ f(X)
#'   \item log(Y) ~ log(X)
#'   \item log(Y) ~ log(f(X))
#' }
#'
#' For the beta method, Y = BETA^2. For the p-value method, Y = h^2.
#'
#' \strong{Model Selection Criteria:}
#' \itemize{
#'   \item \strong{R^2:} Proportion of variance explained. Higher is better.
#'     Preferred when sample size is large (>100).
#'   \item \strong{Adjusted R^2:} R^2 adjusted for number of predictors. Higher is better.
#'     Balances fit and parsimony. Good for small-medium samples.
#'   \item \strong{CV_R2:} Cross-validated R^2. Higher is better.
#'     Most robust against overfitting. Recommended for small samples.
#' }
#'
#' @section Computational Complexity:
#' O(8 * n) where n is the number of training variants
#'
#' @examples
#' \dontrun{
#' # Run get_B with both methods
#' result <- get_B(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   training_P = training_P,
#'   training_N = training_N,
#'   target_trait = "binary",
#'   target_MAF = target_MAF,
#'   target_case_prop = 0.1,
#'   method = "both",
#'   return_full = TRUE
#' )
#'
#' # Compare all models
#' model_comparison <- compare_B_models(result, criterion = "R2", plot = TRUE)
#'
#' # View results sorted by R^2
#' model_comparison[order(-model_comparison$R2), ]
#'
#' # Also works with glow_B_model directly
#' b_model <- train_B_model(...)
#' compare_B_models(b_model, which_method = "beta", plot = TRUE)
#' }
#'
#' @seealso
#' \code{\link{get_B}} for B estimation
#' \code{\link{train_B_model}} for training B estimation models
#' \code{\link{select_best_model}} for the internal model selection function
#' \code{\link{get_B_diagnostics}} for comprehensive diagnostics
#'
#' @export
compare_B_models <- function(B_result,
                              which_method = "both",
                              criterion = "R2",
                              plot = TRUE,
                              verbose = 1) {

  # Extract glow_B_model from either class
  b_model <- .extract_b_model(B_result)

  valid_methods <- c("beta", "pvalue", "both")
  if (!which_method %in% valid_methods) {
    stop("which_method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  valid_criteria <- c("R2", "adj_R2", "CV_R2")
  if (!criterion %in% valid_criteria) {
    stop("criterion must be one of: ", paste(valid_criteria, collapse = ", "))
  }

  # Determine which methods to compare
  methods_to_compare <- character()
  if (which_method == "both" || which_method == "beta") {
    if (!is.null(b_model$models$beta_method)) {
      methods_to_compare <- c(methods_to_compare, "beta_method")
    }
  }
  if (which_method == "both" || which_method == "pvalue") {
    if (!is.null(b_model$models$pvalue_method)) {
      methods_to_compare <- c(methods_to_compare, "pvalue_method")
    }
  }

  if (length(methods_to_compare) == 0) {
    stop("No models available for comparison")
  }

  # Check if training data is available
  if (is.null(b_model$training_data)) {
    stop("Training data not available in B_result object. ",
         "Ensure get_B() was run with return_full = TRUE and ",
         "is from a recent version that stores training data.")
  }

  # Extract training data
  training_MAF <- b_model$training_data$MAF
  training_trait <- b_model$training_data$trait

  # Initialize results data frame
  results <- data.frame(
    model_name = character(),
    formula = character(),
    R2 = numeric(),
    adj_R2 = numeric(),
    CV_R2 = numeric(),
    is_best = logical(),
    method = character(),
    stringsAsFactors = FALSE
  )

  # ========== Extract pre-computed results from all_models_info ==========

  if (is.null(b_model$all_models_info)) {
    stop("Model comparison requires get_B(..., return_full = TRUE)\n",
         "Please re-run get_B() with return_full = TRUE to store complete model information")
  }

  # Extract pre-computed results for each method
  for (method_name in methods_to_compare) {
    method_info <- b_model$all_models_info[[method_name]]

    if (is.null(method_info)) {
      warning("Model information not available for ", method_name, " method. Skipping.")
      next
    }

    # Extract the summary table which already contains all model statistics
    method_table <- method_info$summary_table

    # Add method column
    method_table$method <- method_name

    # Append to results
    results <- rbind(results, method_table)
  }

  # Sort by criterion (all criteria: higher is better)
  if (criterion == "R2") {
    results <- results[order(-results$R2), ]
  } else if (criterion == "adj_R2") {
    results <- results[order(-results$adj_R2), ]
  } else if (criterion == "CV_R2") {
    # Handle NAs by putting them last
    results <- results[order(-ifelse(is.na(results$CV_R2), -Inf, results$CV_R2)), ]
  }

  # Plot if requested
  if (plot && nrow(results) > 0) {
    .plot_model_comparison_table(results, criterion)
  }

  return(results)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Extract glow_B_model from Either Input Class
#'
#' @description
#' Internal helper that accepts either a glow_B_model or glow_B_estimate
#' and returns the glow_B_model. Used by diagnostic functions that accept
#' both classes.
#'
#' @param x Object of class "glow_B_model" or "glow_B_estimate"
#'
#' @return A glow_B_model object
#'
#' @keywords internal
#' @noRd
.extract_b_model <- function(x) {
  if (inherits(x, "glow_B_model")) {
    return(x)
  } else if (inherits(x, "glow_B_estimate")) {
    return(x$model)
  } else {
    stop("B_result must be a 'glow_B_model' or 'glow_B_estimate' object")
  }
}


#' Create Dashboard Plot for B Estimation Diagnostics
#'
#' @description
#' Internal function to create a 2x2 dashboard of diagnostic plots
#' for glow_B_estimate objects. Includes method comparison panel when
#' both methods are available.
#'
#' @param x Object of class "glow_B_estimate"
#' @param which_method Character string: "beta" or "pvalue"
#' @param ... Additional arguments (unused)
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_dashboard <- function(x, which_method = "beta", ...) {
  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Set up 2x2 layout
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  # Plot 1: Model fit
  .plot_model_fit(x, which_method = which_method, add_title = FALSE)

  # Plot 2: Q-Q plot
  model_key <- paste0(which_method, "_method")
  model <- x$model$models[[model_key]]
  qqnorm(residuals(model), main = "Normal Q-Q Plot")
  qqline(residuals(model), col = "red", lwd = 2)

  # Plot 3: Residuals vs Fitted
  plot(fitted(model), residuals(model),
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted",
       pch = 16, col = rgb(0, 0, 0, 0.5))
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(lowess(fitted(model), residuals(model)), col = "blue", lwd = 2)

  # Plot 4: Method comparison or Cook's distance
  if (x$model$method_used == "both") {
    # Method comparison plot (uses B predictions from glow_B_estimate)
    .plot_method_comparison_simple(x)
  } else {
    # Cook's distance
    cooks_d <- cooks.distance(model)
    n <- length(cooks_d)
    plot(cooks_d, type = "h", main = "Cook's Distance",
         xlab = "Observation", ylab = "Cook's distance",
         col = ifelse(cooks_d > 4/n, "red", "black"))
    abline(h = 4/n, col = "red", lty = 2)
  }

  # Overall title
  method_text <- if (which_method == "beta") "Direct Beta Method" else "P-value/Z-score Method"
  mtext(paste("B Estimation Diagnostics -", method_text),
        outer = TRUE, cex = 1.2, font = 2)

  invisible(NULL)
}


#' Create Dashboard Plot for glow_B_model (No Predictions)
#'
#' @description
#' Internal function to create a 2x2 dashboard of model-level diagnostic
#' plots for glow_B_model objects. Since no B predictions are available,
#' the fourth panel always shows Cook's distance (no method comparison).
#'
#' @param x List with $model containing a glow_B_model
#' @param which_method Character string: "beta" or "pvalue"
#' @param ... Additional arguments (unused)
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_dashboard_model <- function(x, which_method = "beta", ...) {
  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Set up 2x2 layout
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  # Plot 1: Model fit
  .plot_model_fit(x, which_method = which_method, add_title = FALSE)

  # Plot 2: Q-Q plot
  model_key <- paste0(which_method, "_method")
  model <- x$model$models[[model_key]]
  qqnorm(residuals(model), main = "Normal Q-Q Plot")
  qqline(residuals(model), col = "red", lwd = 2)

  # Plot 3: Residuals vs Fitted
  plot(fitted(model), residuals(model),
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted",
       pch = 16, col = rgb(0, 0, 0, 0.5))
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  lines(lowess(fitted(model), residuals(model)), col = "blue", lwd = 2)

  # Plot 4: Cook's distance (always, since no predictions available)
  cooks_d <- cooks.distance(model)
  n <- length(cooks_d)
  plot(cooks_d, type = "h", main = "Cook's Distance",
       xlab = "Observation", ylab = "Cook's distance",
       col = ifelse(cooks_d > 4/n, "red", "black"))
  abline(h = 4/n, col = "red", lty = 2)

  # Overall title
  method_text <- if (which_method == "beta") "Direct Beta Method" else "P-value/Z-score Method"
  mtext(paste("B Model Diagnostics -", method_text),
        outer = TRUE, cex = 1.2, font = 2)

  invisible(NULL)
}


#' Plot Model Fit
#'
#' @description
#' Internal function to plot training data with fitted curve.
#'
#' @param x Object with $model containing a glow_B_model
#' @param which_method Character string: "beta" or "pvalue"
#' @param add_title Logical, whether to add title
#' @param ... Additional arguments
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_model_fit <- function(x, which_method = "beta", add_title = TRUE, ...) {
  # Get model from nested glow_B_model
  model_key <- paste0(which_method, "_method")
  model <- x$model$models[[model_key]]

  # Extract data from model
  model_data <- model$model
  Y <- model_data[[1]]
  X_predictor <- model_data[[2]]

  # Try to extract or reconstruct original X (MAF)
  # Check if X is in the model data
  if ("X" %in% names(model_data)) {
    X_maf <- model_data$X
  } else if ("fX" %in% names(model_data)) {
    # f(X) = X(1-X), solve for X is complex
    # Use fX values as-is for plotting
    X_maf <- model_data$fX
  } else if ("logX" %in% names(model_data)) {
    X_maf <- exp(model_data$logX)
  } else if ("logfX" %in% names(model_data)) {
    X_maf <- exp(model_data$logfX)
  } else {
    # Fallback: use predictor as X
    X_maf <- X_predictor
  }

  # Check if Y was log-transformed
  response_name <- names(model_data)[1]
  y_is_log <- grepl("^logY", response_name)

  # Get fitted values on transformed scale
  fitted_values_on_model_scale <- fitted(model)

  # Plot
  if (y_is_log) {
    # If Y was log-transformed, plot on log scale or back-transform
    Y_plot <- Y  # Already log(Y)
    fitted_plot <- fitted_values_on_model_scale
    ylab <- expression(log(Y))
  } else {
    Y_plot <- Y
    fitted_plot <- fitted_values_on_model_scale
    ylab <- "Y (effect size measure)"
  }

  # Order for plotting fitted line
  maf_order <- order(X_maf)

  # Plot data points
  plot(X_maf, Y_plot,
       xlab = "MAF (training data)",
       ylab = ylab,
       main = if (add_title) "Model Fit: Training Data vs Fitted Curve" else "",
       pch = 16, col = rgb(0, 0, 1, 0.5), cex = 0.8)

  # Add fitted line
  lines(X_maf[maf_order], fitted_plot[maf_order], col = "red", lwd = 2)

  # Add R^2
  r2 <- summary(model)$r.squared
  legend("topright", legend = sprintf("R^2 = %.3f", r2),
         bty = "n", text.col = "red")

  invisible(NULL)
}


#' Plot Residual Diagnostics
#'
#' @description
#' Internal function to create residual diagnostic plots.
#'
#' @param x Object with $model containing a glow_B_model
#' @param which_method Character string: "beta" or "pvalue"
#' @param ... Additional arguments
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_residuals <- function(x, which_method = "beta", ...) {
  model_key <- paste0(which_method, "_method")
  model <- x$model$models[[model_key]]

  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Create 2x2 layout for standard diagnostic plots
  par(mfrow = c(2, 2))
  plot(model)

  invisible(NULL)
}


#' Plot Method Comparison
#'
#' @description
#' Internal function to create detailed method comparison plot.
#' Requires a glow_B_estimate with B predictions from both methods.
#'
#' @param x Object of class "glow_B_estimate"
#' @param ... Additional arguments
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_method_comparison <- function(x, ...) {
  if (x$model$method_used != "both") {
    stop("Method comparison requires both methods to be run. Use method = 'both' in get_B()")
  }

  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Create 2x2 layout
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  # Plot 1: Scatterplot with diagonal
  plot(x$B_beta_method, x$B_pvalue_method,
       xlab = "B (Direct Beta Method)",
       ylab = "B (P-value/Z-score Method)",
       main = "Method Comparison",
       pch = 16, col = rgb(0, 0, 1, 0.5))
  abline(0, 1, col = "red", lwd = 2, lty = 2)  # y = x line

  # Add correlation
  pred_comp <- x$model$comparison$prediction
  cor_val <- if (!is.null(pred_comp)) pred_comp$correlation else
    cor(x$B_beta_method, x$B_pvalue_method, use = "complete.obs")
  legend("topleft", legend = sprintf("r = %.3f", cor_val),
         bty = "n", text.col = "blue")

  # Plot 2: Bland-Altman plot (difference vs average)
  avg <- (x$B_beta_method + x$B_pvalue_method) / 2
  diff <- x$B_beta_method - x$B_pvalue_method
  plot(avg, diff,
       xlab = "Average B",
       ylab = "Difference (Beta - P-value)",
       main = "Bland-Altman Plot",
       pch = 16, col = rgb(0, 0, 1, 0.5))
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  abline(h = mean(diff) + 1.96 * sd(diff), col = "red", lwd = 1, lty = 3)
  abline(h = mean(diff) - 1.96 * sd(diff), col = "red", lwd = 1, lty = 3)

  # Plot 3: Percent difference
  percent_diff <- 100 * abs(x$B_beta_method - x$B_pvalue_method) / avg
  plot(avg, percent_diff,
       xlab = "Average B",
       ylab = "Percent Difference (%)",
       main = "Percent Difference",
       pch = 16, col = rgb(0, 0, 1, 0.5))
  abline(h = 10, col = "orange", lwd = 1, lty = 2)
  abline(h = 20, col = "red", lwd = 1, lty = 2)

  # Plot 4: Distribution of differences
  hist(diff, breaks = 20,
       main = "Distribution of Differences",
       xlab = "Difference (Beta - P-value)",
       col = rgb(0, 0, 1, 0.5), border = "white")
  abline(v = 0, col = "red", lwd = 2, lty = 2)

  par(old_par)
  invisible(NULL)
}


#' Simple Method Comparison Plot
#'
#' @description
#' Internal function for simplified method comparison (for dashboard).
#' Requires a glow_B_estimate with B predictions from both methods.
#'
#' @param x Object of class "glow_B_estimate"
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_method_comparison_simple <- function(x) {
  plot(x$B_beta_method, x$B_pvalue_method,
       xlab = "Beta Method", ylab = "P-value Method",
       main = "Method Comparison",
       pch = 16, col = rgb(0, 0, 1, 0.5))
  abline(0, 1, col = "red", lwd = 2, lty = 2)

  pred_comp <- x$model$comparison$prediction
  cor_val <- if (!is.null(pred_comp)) pred_comp$correlation else
    cor(x$B_beta_method, x$B_pvalue_method, use = "complete.obs")
  legend("topleft", legend = sprintf("r = %.3f", cor_val),
         bty = "n", text.col = "blue", cex = 0.9)

  invisible(NULL)
}


#' Plot All Candidate Models
#'
#' @description
#' Internal function to visualize all 8 candidate models by calling
#' compare_B_models() with plotting enabled.
#'
#' @param x Object with $model containing a glow_B_model
#' @param which_method Character string: "beta" or "pvalue"
#' @param ... Additional arguments
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_all_models <- function(x, which_method = "beta", ...) {
  # Check if training data is available
  if (is.null(x$model$training_data)) {
    stop("Cannot create model comparison plot: training data not available.\n",
         "Possible causes:\n",
         "  1. get_B() was not run with return_full = TRUE\n",
         "  2. Using an older version of get_B() that doesn't store training data\n",
         "Solution: Re-run get_B() with return_full = TRUE")
  }

  # Determine which method to use
  if (which_method == "primary") {
    which_method <- if (x$model$method_used == "pvalue_method") "pvalue" else "beta"
  }

  # Check if the requested method has required data
  if (which_method == "beta" && is.null(x$model$training_data$BETA)) {
    stop("Cannot create model comparison plot for beta method: BETA values not available.\n",
         "The get_B() function was run with p-value method only.\n",
         "Solution: Either use which_method = 'pvalue' or re-run get_B() with method = 'beta' or 'both'")
  }

  if (which_method == "pvalue" &&
      ((is.null(x$model$training_data$P) && is.null(x$model$training_data$P_mlog10)) ||
       is.null(x$model$training_data$N))) {
    stop("Cannot create model comparison plot for p-value method: P-values or sample sizes not available.\n",
         "The get_B() function was run with beta method only.\n",
         "Solution: Either use which_method = 'beta' or re-run get_B() with method = 'pvalue' or 'both'")
  }

  # Determine criterion to use
  criterion <- if (!is.null(x$model$selection_criterion)) {
    x$model$selection_criterion
  } else {
    "R2"
  }

  # Check if a plotting device is available
  if (dev.cur() == 1) {
    # No device is open, open a default one
    tryCatch({
      dev.new()
    }, error = function(e) {
      # If dev.new() fails (e.g., in non-interactive mode), try pdf
      pdf(file = tempfile(fileext = ".pdf"))
      on.exit(dev.off(), add = TRUE)
    })
  }

  # Construct input for compare_B_models (accepts either class)
  # If x has a $model field (our wrapper), unwrap to the original input
  compare_input <- if (inherits(x, "glow_B_estimate") || inherits(x, "glow_B_model")) {
    x
  } else {
    # x is a wrapper list with $model being the glow_B_model
    x$model
  }

  # Use compare_B_models to fit all models and create plot
  tryCatch({
    model_comparison <- compare_B_models(
      B_result = compare_input,
      which_method = which_method,
      criterion = criterion,
      plot = TRUE,
      verbose = 0  # Suppress warnings during plotting
    )

    # Only add title if a plot was actually created
    tryCatch({
      par("usr")  # This will error if no plot exists
      title(sub = sprintf("All 8 candidate models compared using %s", criterion),
            cex.sub = 0.8, line = 4)
    }, error = function(e) {
      # No plot exists, skip title
    })

  }, error = function(e) {
    stop("Failed to create model comparison plot.\n",
         "Error: ", e$message, "\n",
         "Debugging suggestions:\n",
         "  1. Check if training data has valid values (no NAs, positive values for log transforms)\n",
         "  2. Try calling compare_B_models() directly:\n",
         "     compare_B_models(your_result, which_method = '", which_method,
         "', verbose = 1)\n",
         "  3. Check if a plotting device is open (try pdf() or dev.new() first)")
  })

  invisible(NULL)
}


#' Generate Diagnostic Summary Text
#'
#' @description
#' Internal function to create human-readable diagnostic summary.
#'
#' @param diagnostics List of diagnostic results
#' @param B_result Original B estimation result object (glow_B_estimate or glow_B_model)
#'
#' @return Character string with formatted summary
#'
#' @keywords internal
#' @noRd
.generate_diagnostic_summary <- function(diagnostics, B_result) {
  # Extract glow_B_model for field access
  b_model <- .extract_b_model(B_result)

  lines <- character()

  lines <- c(lines, "========================================")
  lines <- c(lines, "B Estimation Diagnostic Report")
  lines <- c(lines, "========================================\n")

  # Method used
  lines <- c(lines, sprintf("Method(s) used: %s\n", b_model$method_used))

  # Model fit statistics
  lines <- c(lines, "Model Fit Statistics:")
  lines <- c(lines, "---------------------")

  if (!is.null(diagnostics$model_fit$beta)) {
    lines <- c(lines, "\nDirect Beta Method:")
    lines <- c(lines, sprintf("  R^2 = %.3f (Adjusted R^2 = %.3f)",
                             diagnostics$model_fit$beta$R2,
                             diagnostics$model_fit$beta$adj_R2))
    lines <- c(lines, sprintf("  Formula: %s", diagnostics$model_fit$beta$formula))
    lines <- c(lines, sprintf("  N observations: %d", diagnostics$model_fit$beta$n_obs))
  }

  if (!is.null(diagnostics$model_fit$pvalue)) {
    lines <- c(lines, "\nP-value/Z-score Method:")
    lines <- c(lines, sprintf("  R^2 = %.3f (Adjusted R^2 = %.3f)",
                             diagnostics$model_fit$pvalue$R2,
                             diagnostics$model_fit$pvalue$adj_R2))
    lines <- c(lines, sprintf("  Formula: %s", diagnostics$model_fit$pvalue$formula))
    lines <- c(lines, sprintf("  N observations: %d", diagnostics$model_fit$pvalue$n_obs))
  }

  # Residual statistics
  if (length(diagnostics$residual_stats) > 0) {
    lines <- c(lines, "\n\nResidual Statistics:")
    lines <- c(lines, "--------------------")
    lines <- c(lines, sprintf("  Mean: %.4f, SD: %.4f",
                             diagnostics$residual_stats$mean,
                             diagnostics$residual_stats$sd))
    lines <- c(lines, sprintf("  Range: [%.4f, %.4f]",
                             diagnostics$residual_stats$min,
                             diagnostics$residual_stats$max))
    lines <- c(lines, sprintf("  Median: %.4f (Q1=%.4f, Q3=%.4f)",
                             diagnostics$residual_stats$median,
                             diagnostics$residual_stats$q1,
                             diagnostics$residual_stats$q3))

    if (!is.na(diagnostics$residual_stats$normality_p)) {
      lines <- c(lines, sprintf("  Normality test p-value: %.4f",
                               diagnostics$residual_stats$normality_p))
    }
  }

  # Outlier detection
  if (length(diagnostics$outliers) > 0) {
    lines <- c(lines, "\n\nOutlier Detection:")
    lines <- c(lines, "------------------")
    lines <- c(lines, sprintf("  Number of outliers: %d (%.1f%%)",
                             diagnostics$outliers$n_outliers,
                             diagnostics$outliers$proportion * 100))
    if (diagnostics$outliers$n_outliers > 0) {
      lines <- c(lines, sprintf("  Outlier indices: %s",
                               paste(utils::head(diagnostics$outliers$indices, 10), collapse = ", ")))
      if (diagnostics$outliers$n_outliers > 10) {
        lines <- c(lines, sprintf("  ... and %d more",
                                 diagnostics$outliers$n_outliers - 10))
      }
    }
  }

  # Method comparison
  if (!is.null(diagnostics$method_comparison)) {
    lines <- c(lines, "\n\nMethod Comparison:")
    lines <- c(lines, "------------------")
    lines <- c(lines, sprintf("  Correlation: %.4f",
                             diagnostics$method_comparison$correlation))
    lines <- c(lines, sprintf("  RMSE: %.4f",
                             diagnostics$method_comparison$rmse))
    lines <- c(lines, sprintf("  Mean percent difference: %.2f%%",
                             diagnostics$method_comparison$mean_percent_diff))
    lines <- c(lines, sprintf("  Max percent difference: %.2f%%",
                             diagnostics$method_comparison$max_percent_diff))
  }

  # QC flags
  if (length(diagnostics$qc_flags) > 0) {
    lines <- c(lines, "\n\nQuality Control Flags:")
    lines <- c(lines, "----------------------")
    for (flag in diagnostics$qc_flags) {
      lines <- c(lines, sprintf("  [!] %s", flag))
    }
  } else {
    lines <- c(lines, "\n\n[OK] No quality control issues detected")
  }

  lines <- c(lines, "\n========================================\n")

  return(paste(lines, collapse = "\n"))
}


#' Plot Model Comparison Table
#'
#' @description
#' Internal function to visualize model comparison results.
#'
#' @param results Data frame from compare_B_models
#' @param criterion Character string: selection criterion
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.plot_model_comparison_table <- function(results, criterion) {
  # Save current par settings
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Increase right margin to accommodate value labels
  par(mar = c(5, 8, 4, 4))

  # Bar plot of criterion values
  if (criterion == "R2") {
    values <- results$R2
    xlab <- "R^2 (higher is better)"
    colors <- ifelse(results$is_best, "darkgreen", "lightblue")
  } else {
    values <- results[[criterion]]
    xlab <- paste(criterion, "(higher is better)")
    colors <- ifelse(results$is_best, "darkgreen", "lightblue")
  }

  # Extend xlim to provide space for labels
  max_val <- max(values, na.rm = TRUE)
  xlim_range <- c(0, max_val * 1.15)  # Add 15% extra space for labels

  # Create bar plot with extended x-axis
  bp <- barplot(values, names.arg = results$model_name,
                horiz = TRUE, las = 1, col = colors,
                xlab = xlab, main = "Comparison of All Candidate Models",
                cex.names = 0.7, xlim = xlim_range)

  # Add value labels on bars (slightly offset from bar ends)
  label_x <- values + max_val * 0.01  # Small offset from bar end
  text(x = label_x, y = bp, labels = sprintf("%.3f", values),
       pos = 4, cex = 0.7, xpd = TRUE)

  # Place legend in upper-right corner to avoid overlapping with bars
  legend("topright", legend = c("Selected Model", "Other Models"),
         fill = c("darkgreen", "lightblue"), bty = "n", cex = 0.9)

  invisible(NULL)
}


#' Print Method for B Diagnostics
#'
#' @description
#' S3 print method for glow_B_diagnostics objects.
#'
#' @param x Object of class "glow_B_diagnostics"
#' @param ... Additional arguments (unused)
#'
#' @return The object x (invisibly)
#'
#' @export
print.glow_B_diagnostics <- function(x, ...) {
  cat(x$summary_text)
  invisible(x)
}
