########## Model Selection Functions for B Estimation ##########
#
# This file contains functions for selecting the best regression model
# when estimating effect sizes (B) from minor allele frequency (MAF).
# It supports multiple selection criteria and custom models.

#################### EXPORTED FUNCTIONS ####################

#' Select Best Regression Model for B Estimation
#'
#' @description
#' Selects the best linear regression model for predicting effect size B based
#' on minor allele frequency (MAF). This function fits 8 different models with
#' various transformations of X (MAF) and Y (effect size measure) and selects
#' the best model based on the specified criterion.
#'
#' @param X Numeric vector of minor allele frequencies (MAF values)
#' @param Y Numeric vector of effect size measures (typically BETA^2 or related)
#' @param criterion Character string specifying selection criterion:
#'   \itemize{
#'     \item "R2" (default): R-squared, proportion of variance explained
#'     \item "adj_R2": Adjusted R-squared, penalizes model complexity
#'     \item "CV_R2": Cross-validated R-squared (leave-one-out)
#'   }
#'   See Details section for guidance on choosing a criterion.
#'   Note: "CV" is deprecated but still accepted (use "CV_R2" instead).
#'   AIC and BIC have been removed - see Details for explanation
#' @param custom_models Optional list of custom model formulas to include in
#'   candidate set. Each element should be a formula object
#' @param return_all Logical. If TRUE, returns list of all fitted models with
#'   their fit statistics. If FALSE (default), returns only the best model
#' @param verbose Logical. If TRUE (default), prints detailed model comparison
#'   table and selection summary. Set to FALSE to suppress output in loops or
#'   simulations. Default is TRUE for transparency in model selection.
#'
#' @return If return_all=FALSE: A fitted linear model object (class "lm")
#'   representing the best model based on the specified criterion.
#'   If return_all=TRUE: A list with elements:
#'   \itemize{
#'     \item best_model: The selected fitted lm object
#'     \item best_name: Character name of the selected model
#'     \item all_models_info: A list with
#'       \itemize{
#'         \item models: A named list with one element per candidate model, each
#'           a list containing model (the fitted lm), formula, formula_text, R2,
#'           adj_R2, CV_R2, sigma, and coefficients
#'         \item summary_table: A data frame with columns model_name, formula,
#'           R2, adj_R2, CV_R2, and is_best (one row per model)
#'         \item criterion_used: The criterion used for selection
#'         \item training_data: A list with the original X and Y
#'       }
#'   }
#'
#' @details
#' The function considers the following transformations:
#' \itemize{
#'   \item X transformations: X, f(X) = X(1-X), log(X), log(f(X))
#'   \item Y transformations: Y, log(Y)
#' }
#'
#' This results in 8 possible models:
#' \enumerate{
#'   \item Y ~ X
#'   \item Y ~ f(X)
#'   \item Y ~ log(X)
#'   \item Y ~ log(f(X))
#'   \item log(Y) ~ X
#'   \item log(Y) ~ f(X)
#'   \item log(Y) ~ log(X)
#'   \item log(Y) ~ log(f(X))
#' }
#'
#' The model with the best criterion value is returned. This approach captures
#' the relationship between MAF and effect size, which can vary across different
#' genetic architectures.
#'
#' \strong{Selection Criteria - When to Use Which:}
#'
#' \describe{
#'   \item{\strong{R^2 (R-squared)}}{
#'     \itemize{
#'       \item \strong{What it measures:} Proportion of variance explained in training data
#'       \item \strong{Higher is better}
#'       \item \strong{Best for:} Large samples (n > 100), exploratory analysis
#'       \item \strong{Advantages:} Fast, intuitive, widely used
#'       \item \strong{Disadvantages:} Can overfit with small samples, always favors complex models
#'       \item \strong{Default for:} Historical compatibility
#'     }
#'   }
#'   \item{\strong{Adjusted R^2 (adj_R2)}}{
#'     \itemize{
#'       \item \strong{What it measures:} R^2 adjusted for number of predictors, penalizes complexity
#'       \item \strong{Higher is better}
#'       \item \strong{Best for:} Small-medium samples (n < 100), comparing models with different complexity
#'       \item \strong{Advantages:} Fast, balances fit and parsimony, more conservative than R^2
#'       \item \strong{Disadvantages:} Still uses in-sample fit (not true out-of-sample validation)
#'       \item \strong{Good choice:} When you want simplicity preference without CV computational cost
#'       \item \strong{Comparison with R^2:} Will prefer simpler models when fit is similar
#'     }
#'   }
#'   \item{\strong{CV_R2 (Cross-Validation R^2)}}{
#'     \itemize{
#'       \item \strong{What it measures:} Out-of-sample prediction accuracy (LOOCV)
#'       \item \strong{Higher is better}
#'       \item \strong{Best for:} Small-medium samples (10 < n < 200), when generalization matters
#'       \item \strong{Advantages:} Most robust against overfitting, best generalization estimate
#'       \item \strong{Disadvantages:} Slower (n model fits per candidate), can be noisy with very small n
#'       \item \strong{Recommended for:} B estimation from literature data (typical n=20-100)
#'     }
#'   }
#' }
#'
#' \strong{Quick Guide:}
#' \itemize{
#'   \item \strong{Small sample (n < 50):} Use \strong{CV_R2} (most robust)
#'   \item \strong{Medium sample (50 < n < 200):} Use \strong{CV_R2} or \strong{adj_R2}
#'   \item \strong{Large sample (n > 200):} Any criterion works, \strong{R^2} is fastest
#'   \item \strong{When in doubt:} Use \strong{CV_R2} (most robust)
#'   \item \strong{For speed:} Use \strong{adj_R2} (no iteration)
#'   \item \strong{Balance of speed and robustness:} Use \strong{adj_R2}
#' }
#'
#' @section Computational Complexity:
#' \itemize{
#'   \item \strong{R^2, adj_R2:} O(k) where k = number of candidate models (typically 8)
#'   \item \strong{CV_R2:} O(k * n) where n = sample size (LOOCV requires n fits per model)
#' }
#' For typical use: n=50, k=8 -> ~400 model fits for CV_R2. Fast with lm() (< 1 second)
#'
#' @section Assumptions:
#' \itemize{
#'   \item X and Y have the same length
#'   \item X values are in (0, 1) for MAF
#'   \item Y values are positive for log transformation
#'   \item X values are not exactly 0 or 1 for log transformations
#' }
#'
#' @note AIC and BIC are not available as selection criteria because they are
#' not appropriate for comparing models with different response transformations
#' (e.g., Y vs log(Y) vs sqrt(Y)). Use R^2, adj_R^2, or CV_R^2 instead, which
#' properly account for scale differences between transformations.
#'
#' @examples
#' # Simulate MAF and effect sizes
#' set.seed(123)
#' MAF <- runif(100, 0.01, 0.5)
#' BETA <- sqrt(0.5 * MAF * (1 - MAF)) + rnorm(100, 0, 0.1)
#' Y <- BETA^2
#'
#' # Select best model using R^2 (default, fastest)
#' best_model_r2 <- select_best_model(MAF, Y)
#' summary(best_model_r2)
#'
#' # Select using CV_R2 (recommended for small samples)
#' best_model_cv <- select_best_model(MAF, Y, criterion = "CV_R2")
#'
#' # Select using adjusted R^2 (balances fit and simplicity)
#' best_model_adjr2 <- select_best_model(MAF, Y, criterion = "adj_R2")
#'
#' # Suppress verbose output (default is verbose = TRUE)
#' best_model_quiet <- select_best_model(MAF, Y, criterion = "CV_R2", verbose = FALSE)
#'
#' # Get all models for comparison (returns a richer structure, see Value)
#' all_models <- select_best_model(MAF, Y, return_all = TRUE)
#'
#' # The ready-made summary table has R2, adj_R2, and CV_R2 for every candidate
#' print(all_models$all_models_info$summary_table)
#'
#' # Equivalently, build the comparison from the per-model list
#' model_comparison <- sapply(all_models$all_models_info$models, function(m) {
#'   c(R2 = m$R2, adj_R2 = m$adj_R2, CV_R2 = m$CV_R2)
#' })
#' print(model_comparison)
#'
#' # Add custom model
#' custom_formula <- formula(Y ~ poly(X, 2))  # Polynomial model
#' best_custom <- select_best_model(MAF, Y,
#'                                  custom_models = list(custom_formula))
#'
#' @seealso
#' \code{\link{get_B}} for B estimation using the selected model
#' \code{\link{compare_B_models}} for comparing fitted B models
#'
#' @export
select_best_model <- function(X, Y,
                              criterion = "R2",
                              custom_models = NULL,
                              return_all = FALSE,
                              verbose = TRUE) {
  # Input validation for NA/Inf
  validate_numeric_input(X, "X (MAF)", allow_negative = FALSE, allow_zero = FALSE)
  validate_numeric_input(Y, "Y (effect size)", allow_negative = FALSE, allow_zero = FALSE)

  # Additional validation
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length")
  }
  if (any(X <= 0 | X >= 1)) {
    stop("X (MAF) values must be in the interval (0, 1)")
  }
  if (any(Y <= 0)) {
    stop("Y values must be positive for log transformation")
  }

  # Handle backward compatibility: "CV" -> "CV_R2" with deprecation warning
  if (criterion == "CV") {
    warning("The 'CV' criterion name is deprecated. Please use 'CV_R2' instead. ",
            "The old name will be removed in a future version.")
    criterion <- "CV_R2"
  }

  # Validate criterion input
  valid_criteria <- c("R2", "adj_R2", "CV_R2")
  if (!criterion %in% valid_criteria) {
    stop("criterion must be one of: ", paste(valid_criteria, collapse = ", "))
  }

  # Validate custom_models if provided
  if (!is.null(custom_models)) {
    if (!is.list(custom_models)) {
      stop("custom_models must be a list of formula objects")
    }
    if (!all(sapply(custom_models, inherits, "formula"))) {
      stop("All elements of custom_models must be formula objects")
    }
  }

  # ==========  Use shared utility to fit all models ==========

  # Determine if we need to compute CV
  compute_cv <- (criterion == "CV_R2" || return_all)

  # Fit all 8 candidate models plus any custom models using shared utility
  all_models <- .fit_all_candidate_models(
    X = X,
    Y = Y,
    custom_models = custom_models,
    compute_cv = compute_cv,
    verbose = if (verbose) 1 else 0
  )

  # Check if any models succeeded
  if (length(all_models) == 0) {
    stop("All models failed to fit")
  }

  # Create summary table from model results
  summary_table <- do.call(rbind, lapply(names(all_models), function(name) {
    m <- all_models[[name]]
    data.frame(
      model_name = name,
      formula = m$formula_text,
      R2 = m$R2,
      adj_R2 = m$adj_R2,
      CV_R2 = m$CV_R2,
      stringsAsFactors = FALSE
    )
  }))

  # Select best model based on criterion
  best_idx <- switch(criterion,
    R2 = which.max(summary_table$R2),
    adj_R2 = which.max(summary_table$adj_R2),
    CV_R2 = {
      # Get CV R^2 values, handling NAs
      cv_values <- summary_table$CV_R2
      if (all(is.na(cv_values))) {
        warning("All CV calculations failed. Falling back to R^2 criterion.")
        which.max(summary_table$R2)
      } else {
        # Select model with highest CV R^2 (ignoring NAs)
        which.max(ifelse(is.na(cv_values), -Inf, cv_values))
      }
    }
  )

  best_name <- names(all_models)[best_idx]
  best_model <- all_models[[best_name]]$model

  # Mark best model in summary table
  summary_table$is_best <- FALSE
  summary_table$is_best[best_idx] <- TRUE

  # Print verbose output if requested
  if (verbose) {
    # Convert summary table to old format for printing
    results_for_print <- lapply(seq_len(nrow(summary_table)), function(i) {
      list(
        model = all_models[[i]]$model,
        formula = summary_table$formula[i],
        R2 = summary_table$R2[i],
        adj_R2 = summary_table$adj_R2[i],
        CV_R2 = summary_table$CV_R2[i],
        is_best = summary_table$is_best[i]
      )
    })
    .print_model_comparison_table(results_for_print, criterion)
  }

  # Return results
  if (return_all) {
    # Return comprehensive structure with all model information
    return(list(
      best_model = best_model,
      best_name = best_name,
      all_models_info = list(
        models = all_models,           # Full lm objects + metadata for all 8 models
        summary_table = summary_table, # Data frame with all statistics
        criterion_used = criterion,    # Which criterion was used for selection
        training_data = list(X = X, Y = Y)  # Original training data
      )
    ))
  } else {
    # Current behavior - just return best model
    return(best_model)
  }
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Print Model Comparison Table (Internal)
#'
#' @description
#' Internal helper function to print a formatted table comparing all candidate
#' models with their fit statistics. Called when verbose = TRUE in
#' select_best_model().
#'
#' @param results List of model results from select_best_model(), where each
#'   element contains: model, formula, R2, adj_R2, CV_R2, is_best
#' @param criterion Character string indicating which criterion is being used
#'   for selection: "R2", "adj_R2", or "CV_R2"
#'
#' @return NULL (called for side effect of printing to console)
#'
#' @details
#' This function formats and prints:
#' \itemize{
#'   \item Header with title
#'   \item Progress message about fitting models
#'   \item Formatted table with all model statistics
#'   \item Visual indicator of which model was selected
#'   \item Summary line showing the selection criterion and best model
#' }
#'
#' The table uses fixed-width formatting with sprintf() for alignment.
#' CV_R2 column shows "Not computed" if CV was not calculated for a model.
#'
#' @keywords internal
#' @noRd
.print_model_comparison_table <- function(results, criterion) {
  # Extract number of models
  n_models <- length(results)

  # Print header
  cat("\n")
  cat("========================================\n")
  cat("Model Selection for B Estimation\n")
  cat("========================================\n")
  cat(sprintf("Fitting %d candidate models...\n\n", n_models))

  # Determine which criterion column to highlight
  # For display purposes: higher is better for R2/adj_R2/CV_R2
  criterion_display <- switch(criterion,
    R2 = "R^2",
    adj_R2 = "Adj. R^2",
    CV_R2 = "CV R^2"
  )

  # Find the best model index
  best_idx <- which(sapply(results, function(x) x$is_best))

  # Print table header
  cat("Model Comparison:\n")
  cat(sprintf("%-4s %-20s %8s %8s %10s\n",
              "ID", "Formula", "R^2", "Adj.R^2", "CV R^2"))
  cat(sprintf("%-4s %-20s %8s %8s %10s\n",
              "----", "--------------------", "--------", "--------",
              "----------"))

  # Print each model's statistics
  for (i in seq_along(results)) {
    m <- results[[i]]

    # Format formula (truncate if too long)
    fmla_str <- m$formula
    if (nchar(fmla_str) > 20) {
      fmla_str <- paste0(substr(fmla_str, 1, 17), "...")
    }

    # Format CV_R2 (show "Not computed" if NA)
    cv_str <- if (is.na(m$CV_R2)) {
      "---"
    } else {
      sprintf("%.4f", m$CV_R2)
    }

    # Add selection marker for best model
    marker <- if (m$is_best) " *" else "  "

    # Print row
    cat(sprintf("%s%-2d %-20s %8.4f %8.4f %10s",
                marker, i, fmla_str,
                m$R2, m$adj_R2, cv_str))

    # Add arrow indicator for selected model
    if (m$is_best) {
      cat("  <- SELECTED")
    }
    cat("\n")
  }

  cat("\n")

  # Print selection summary
  best_model <- results[[best_idx]]
  cat(sprintf("Selection Criterion: %s\n", criterion_display))
  cat(sprintf("Best Model: %s", best_model$formula))

  # Add the criterion value for the best model
  criterion_value <- switch(criterion,
    R2 = sprintf(" (R^2 = %.4f)", best_model$R2),
    adj_R2 = sprintf(" (Adj. R^2 = %.4f)", best_model$adj_R2),
    CV_R2 = if (!is.na(best_model$CV_R2)) {
      sprintf(" (CV R^2 = %.4f)", best_model$CV_R2)
    } else {
      " (CV R^2 not available)"
    }
  )
  cat(criterion_value)
  cat("\n")

  cat("========================================\n")
  cat("\n")

  invisible(NULL)
}
