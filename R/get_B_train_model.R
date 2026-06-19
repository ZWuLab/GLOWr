########## Train B Estimation Model ##########
#
# This file contains train_B_model() for training regression models that
# predict effect sizes (B) from minor allele frequency (MAF). It also
# defines the glow_B_model S3 class returned by train_B_model().
#
# This function extracts the training logic from get_B(), separating model
# training from prediction. The companion predict_B() function applies the
# trained model to new MAF values.
#
# EXPORTED FUNCTIONS:
#   train_B_model()       - Train B estimation model(s) from training data
#   print.glow_B_model()  - S3 print method for glow_B_model objects

#################### EXPORTED MAIN FUNCTIONS ####################

#' Train a B Estimation Model from Training Data
#'
#' @description
#' Trains a regression model to predict allelic effect sizes (B) from minor
#' allele frequency (MAF), without requiring target data. This function
#' supports two estimation methods: (1) direct beta method using effect
#' sizes, and (2) p-value/Z-score method using a trait-independent h-squared
#' transformation. The trained model can be applied to new MAF values via
#' \code{\link{predict_B}}.
#'
#' This function performs the training step only. It does not produce B
#' predictions -- use \code{\link{predict_B}} for that. This separation
#' enables "train once, predict many times" workflows.
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
#'   training_P_mlog10 should be provided (not both)
#' @param training_P_mlog10 Numeric vector of -log10 transformed p-values from
#'   training data. Use this instead of training_P when dealing with extremely
#'   small p-values to avoid numerical underflow
#' @param training_N Numeric vector of sample sizes from training data. Required
#'   for p-value/Z-score method
#' @param method Character string specifying which method to use: "auto"
#'   (default), "beta", "pvalue", or "both". "auto" uses both methods when
#'   the beta method is possible and p-value data is available, otherwise uses
#'   the available method. "both" runs both methods for comparison
#' @param selection_criterion Character string specifying model selection
#'   criterion: "R2" (default), "adj_R2", or "CV_R2". See
#'   \code{\link{select_best_model}} for guidance on choosing a criterion
#' @param custom_models Optional list of custom model formulas to include in
#'   candidate model set. Each element should be a formula object
#' @param outlier_method Character string specifying outlier detection method:
#'   "none" (default), "statistical" (Cook's distance), "biological"
#'   (implausible MAF-effect combinations), or "both"
#' @param outlier_action Character string specifying action for detected
#'   outliers: "flag" (default, report but don't remove) or "remove"
#'   (exclude from analysis)
#' @param cook_threshold Numeric threshold multiplier for Cook's distance
#'   outlier detection (default: 4). Outliers are defined as points with
#'   Cook's D > cook_threshold/n
#' @param biological_rules Optional list of custom rules for biological
#'   outlier detection
#' @param show_model_selection Logical. If TRUE (default), prints model
#'   selection comparison table. Set to FALSE to suppress in loops
#' @param verbose Integer controlling verbosity: 0=silent, 1=warnings
#'   (default), 2=info messages, 3=debug output
#'
#' @return An object of class \code{glow_B_model} containing:
#'   \describe{
#'     \item{method_used}{Character: "beta_method", "pvalue_method", or "both"}
#'     \item{models}{List with \code{beta_method} and/or \code{pvalue_method}
#'       containing the best fitted lm objects}
#'     \item{all_models_info}{List with \code{beta_method} and/or
#'       \code{pvalue_method} containing complete model comparison info}
#'     \item{outliers}{List with outlier detection results}
#'     \item{training_summary}{List with n_original, n_used,
#'       n_outliers_detected, trait_type}
#'     \item{training_data}{List with MAF, BETA, P, P_mlog10, N, trait
#'       (retained for diagnostics)}
#'     \item{comparison}{List with training-time model comparison when
#'       method="both": selection_criterion, criterion values, and
#'       method_selected}
#'     \item{selection_criterion}{The model selection criterion used}
#'   }
#'
#' @details
#' \strong{Two Estimation Methods:}
#'
#' \strong{Method 1: Beta Method} (\code{beta_method})
#' \itemize{
#'   \item Uses training_BETA values directly
#'   \item Fits regression model: \eqn{BETA^2 \sim f(MAF)}
#'   \item Best when training and target have the same trait type
#' }
#'
#' \strong{Method 2: P-value Method} (\code{pvalue_method})
#' \itemize{
#'   \item Converts p-values to trait-independent h-squared
#'   \item Fits regression model: \eqn{h^2 \sim f(MAF)}
#'   \item Can transfer information between different trait types
#' }
#'
#' \strong{Method Selection Logic:}
#' \itemize{
#'   \item method="auto": Uses "both" if beta method possible + p-value data
#'     available; otherwise uses the available method
#'   \item method="beta": Only beta method (requires training_BETA)
#'   \item method="pvalue": Only p-value method (requires training_P/training_N)
#'   \item method="both": Runs both, compares at training time using
#'     selection_criterion
#' }
#'
#' \strong{Training-Time Comparison (method="both"):}
#'
#' When both methods are run, the function compares model quality using the
#' selection_criterion (R2/adj_R2/CV_R2) applied to each method's best model.
#' The method with the better criterion value is recorded as
#' \code{comparison$method_selected}. Note that the two methods model different
#' response variables (BETA^2 vs h^2), so this comparison reflects which
#' transformation better captures the MAF-effect relationship.
#'
#' Full prediction-based comparison (correlation, RMSE between B vectors)
#' requires target data and is performed by \code{\link{get_B}} after
#' calling \code{\link{predict_B}}.
#'
#' @section Computational Complexity:
#' O(8 * n_train) where n_train = length(training_MAF) (fitting 8 models).
#' With CV_R2 criterion: O(8 * n_train^2) due to leave-one-out CV.
#'
#' @examples
#' \dontrun{
#' # Example 1: Train with beta method
#' set.seed(123)
#' training_MAF <- runif(50, 0.001, 0.3)
#' training_BETA <- sqrt(0.5 * training_MAF * (1 - training_MAF)) * 0.1
#'
#' b_model <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   method = "beta"
#' )
#' print(b_model)
#'
#' # Predict B for new MAF values
#' B_values <- predict_B(b_model, target_MAF = c(0.01, 0.05, 0.1))
#'
#' # Example 2: Train with p-value method
#' training_P <- runif(50, 1e-8, 0.01)
#' training_N <- rep(5000, 50)
#'
#' b_model_pval <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_P = training_P,
#'   training_N = training_N,
#'   method = "pvalue"
#' )
#'
#' # Example 3: Train with both methods
#' b_model_both <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = training_MAF,
#'   training_BETA = training_BETA,
#'   training_P = training_P,
#'   training_N = training_N,
#'   method = "both",
#'   verbose = 2
#' )
#'
#' # Example 4: Save and reload
#' save_B_model(b_model, "my_model.rds")
#' loaded <- load_B_model("my_model.rds")
#' }
#'
#' @seealso
#' \code{\link{predict_B}} for applying trained models to new data
#' \code{\link{get_B}} for combined train + predict in one call
#' \code{\link{select_best_model}} for model selection details
#' \code{\link{save_B_model}} for persisting trained models
#'
#' @export
train_B_model <- function(training_trait,
                          training_MAF,
                          training_BETA = NULL,
                          training_P = NULL,
                          training_P_mlog10 = NULL,
                          training_N = NULL,
                          method = "auto",
                          selection_criterion = "R2",
                          custom_models = NULL,
                          outlier_method = "none",
                          outlier_action = "flag",
                          cook_threshold = 4,
                          biological_rules = NULL,
                          show_model_selection = TRUE,
                          verbose = 1) {

  # ========== Input Validation ==========

  # Validate trait types
  valid_training_traits <- c("binary", "continuous", "mixed")

  # training_trait can be NULL, "binary", "continuous", or "mixed"
  if (!is.null(training_trait) && !training_trait %in% valid_training_traits) {
    stop("training_trait must be 'binary', 'continuous', 'mixed', or NULL, got: ",
         training_trait)
  }

  # Normalize NULL training_trait to "mixed" for downstream processing
  training_trait_normalized <- if (is.null(training_trait)) "mixed" else training_trait

  # Validate MAF values: must be numeric, non-empty, no NAs in valid range
  if (!is.numeric(training_MAF) || length(training_MAF) == 0) {
    stop("training_MAF must be a non-empty numeric vector")
  }

  # Fold MAF values > 0.5 to minor allele frequency convention
  # After folding, validates that all values are in (0, 0.5]
  training_MAF <- .fold_maf(training_MAF, allow_zero = FALSE,
                            context = "training_MAF")

  # Validate method parameter
  valid_methods <- c("auto", "beta", "pvalue", "both")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Validate selection_criterion parameter
  # Handle backward compatibility for CV -> CV_R2
  if (selection_criterion == "CV") {
    warning("The 'CV' criterion name is deprecated. Please use 'CV_R2' instead. ",
            "The old name will be removed in a future version.")
    selection_criterion <- "CV_R2"
  }

  valid_criteria <- c("R2", "adj_R2", "CV_R2")
  if (!selection_criterion %in% valid_criteria) {
    if (selection_criterion %in% c("AIC", "BIC")) {
      stop(
        "AIC and BIC criteria are inappropriate for comparing models with ",
        "different response transformations (Y vs log(Y)).\n\n",
        "Valid criteria: 'R2', 'adj_R2', 'CV_R2'."
      )
    } else {
      stop("selection_criterion must be one of: ",
           paste(valid_criteria, collapse = ", "))
    }
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

  # Validate outlier_method parameter
  valid_outlier_methods <- c("none", "statistical", "biological", "both")
  if (!outlier_method %in% valid_outlier_methods) {
    stop("outlier_method must be one of: ",
         paste(valid_outlier_methods, collapse = ", "))
  }

  # Validate outlier_action parameter
  valid_outlier_actions <- c("flag", "remove")
  if (!outlier_action %in% valid_outlier_actions) {
    stop("outlier_action must be one of: ",
         paste(valid_outlier_actions, collapse = ", "))
  }

  # Validate cook_threshold
  if (!is.numeric(cook_threshold) || cook_threshold <= 0) {
    stop("cook_threshold must be a positive numeric value")
  }

  # ========== Determine Method to Use ==========

  # Check if beta method is possible (requires BETA and single training trait)
  beta_method_possible <- training_trait_normalized %in% c("binary", "continuous") &&
                          !is.null(training_BETA)

  if (method == "auto") {
    # Auto mode: Use both if BETA + P available, else whichever is available
    has_pvalue_data <- !is.null(training_P) || !is.null(training_P_mlog10)

    if (beta_method_possible && has_pvalue_data) {
      method_to_use <- "both"
      if (verbose >= 2) {
        message("Auto mode: Using both methods for comparison ",
                "(same trait, both BETA and P available)")
      }
    } else if (beta_method_possible) {
      method_to_use <- "beta"
      if (verbose >= 2) {
        message("Auto mode: Using direct beta method (BETA available)")
      }
    } else if (has_pvalue_data) {
      method_to_use <- "pvalue"
      if (verbose >= 2) {
        if (training_trait_normalized == "mixed") {
          message("Auto mode: Using p-value/Z-score method (mixed training traits)")
        } else {
          message("Auto mode: Using p-value/Z-score method")
        }
      }
    } else {
      stop("Insufficient data: need either training_BETA (for beta method) ",
           "or training_P/training_P_mlog10 (for pvalue method)")
    }
  } else {
    method_to_use <- method
  }

  # Validate that beta method is not requested when training trait is mixed
  if (method_to_use %in% c("beta", "both") &&
      training_trait_normalized == "mixed") {
    if (method == "beta") {
      stop("Direct beta method cannot be used when training_trait is 'mixed' ",
           "or NULL. Use method = 'pvalue' instead.")
    } else if (method == "both") {
      if (verbose >= 1) {
        message("Note: Only p-value method available for mixed training traits, ",
                "ignoring 'both' request")
      }
      method_to_use <- "pvalue"
    }
  }

  # Validate required data for chosen method
  if (method_to_use %in% c("beta", "both") && is.null(training_BETA)) {
    stop("training_BETA is required for method '", method_to_use, "'")
  }

  # For p-value method, need either training_P or training_P_mlog10 (not both)
  if (method_to_use %in% c("pvalue", "both")) {
    if (is.null(training_P) && is.null(training_P_mlog10)) {
      stop("Either training_P or training_P_mlog10 is required for method '",
           method_to_use, "'")
    }
    if (!is.null(training_P) && !is.null(training_P_mlog10)) {
      stop("Provide either training_P or training_P_mlog10, not both")
    }
    if (is.null(training_N)) {
      stop("training_N is required for method '", method_to_use, "'")
    }
  }

  # ========== Save Original Training Data Before Outlier Removal ==========
  # These copies are used for the output training_data field so that outlier
  # indices remain consistent with the original data (Issue: outlier removal
  # was overwriting the local vectors, making stored indices invalid).
  original_training_MAF <- training_MAF
  original_training_BETA <- training_BETA
  original_training_P <- training_P
  original_training_P_mlog10 <- training_P_mlog10
  original_training_N <- training_N

  # ========== Outlier Detection ==========

  outliers_beta_method <- NULL
  outliers_pvalue_method <- NULL
  training_indices_kept <- seq_along(training_MAF)

  if (outlier_method != "none") {

    # Detect outliers for beta method if applicable
    if (method_to_use %in% c("beta", "both") && !is.null(training_BETA)) {
      if (verbose >= 2) {
        message("Detecting outliers for beta method...")
      }

      outliers_beta_method <- detect_B_outliers(
        MAF = training_MAF,
        values = training_BETA^2,  # Use BETA^2 as in the model
        method = outlier_method,
        cook_threshold = cook_threshold,
        biological_rules = biological_rules,
        verbose = verbose
      )

      if (length(outliers_beta_method$indices) > 0 && verbose >= 1) {
        message(sprintf(
          "Beta method: Found %d outliers (%.1f%% of training data)",
          length(outliers_beta_method$indices),
          100 * length(outliers_beta_method$indices) / length(training_MAF)))
      }
    }

    # Detect outliers for pvalue method if applicable
    if (method_to_use %in% c("pvalue", "both") &&
        (!is.null(training_P) || !is.null(training_P_mlog10))) {
      if (verbose >= 2) {
        message("Detecting outliers for p-value method...")
      }

      # Calculate h-squared for outlier detection
      Z2 <- .compute_chi2_from_pvalue(
        p = training_P, p_mlog10 = training_P_mlog10)
      h2 <- Z2 / (2 * training_N * training_MAF * (1 - training_MAF))

      outliers_pvalue_method <- detect_B_outliers(
        MAF = training_MAF,
        values = h2,
        method = outlier_method,
        cook_threshold = cook_threshold,
        biological_rules = biological_rules,
        verbose = verbose
      )

      if (length(outliers_pvalue_method$indices) > 0 && verbose >= 1) {
        message(sprintf(
          "P-value method: Found %d outliers (%.1f%% of training data)",
          length(outliers_pvalue_method$indices),
          100 * length(outliers_pvalue_method$indices) / length(training_MAF)))
      }
    }

    # Combine outliers if both methods used
    if (method_to_use == "both") {
      combined_outliers <- unique(c(
        outliers_beta_method$indices,
        outliers_pvalue_method$indices))
    } else if (method_to_use == "beta") {
      combined_outliers <- outliers_beta_method$indices
    } else {
      combined_outliers <- outliers_pvalue_method$indices
    }

    # Handle outliers based on action parameter
    if (outlier_action == "remove" && length(combined_outliers) > 0) {
      if (verbose >= 1) {
        message(sprintf("Removing %d outliers from training data",
                        length(combined_outliers)))
      }

      # Remove outliers from training data
      training_indices_kept <- setdiff(
        seq_along(training_MAF), combined_outliers)
      training_MAF <- training_MAF[training_indices_kept]
      if (!is.null(training_BETA)) {
        training_BETA <- training_BETA[training_indices_kept]
      }
      if (!is.null(training_P)) {
        training_P <- training_P[training_indices_kept]
      }
      if (!is.null(training_P_mlog10)) {
        training_P_mlog10 <- training_P_mlog10[training_indices_kept]
      }
      if (!is.null(training_N)) {
        training_N <- training_N[training_indices_kept]
      }

      # Check if we have enough data left
      if (length(training_MAF) < 3) {
        warning("After removing outliers, only ", length(training_MAF),
                " training points remain. Results may be unreliable.")
      }
    } else if (outlier_action == "flag" && length(combined_outliers) > 0 &&
               verbose >= 1) {
      message("Outliers flagged but not removed ",
              "(use outlier_action='remove' to remove them)")
    }
  }

  # ========== Train Models ==========

  model_beta_method <- NULL
  model_pvalue_method <- NULL
  all_models_beta_method <- NULL
  all_models_pvalue_method <- NULL
  training_MAF_beta <- NULL      # Track beta-specific MAF (after zero removal)
  training_BETA_beta <- NULL     # Track beta-specific BETA (after zero removal)

  # Train beta method model if requested
  if (method_to_use %in% c("beta", "both")) {
    if (verbose >= 2) {
      message("Training beta method model...")
    }

    # Detect and handle zero-valued BETA elements
    # The model selection requires Y = BETA^2 > 0 for log transformation
    zero_indices <- which(training_BETA == 0)

    if (length(zero_indices) > 0) {
      if (verbose >= 1) {
        message(sprintf(
          "Detected %d zero-valued BETA elements (%.1f%% of training data)",
          length(zero_indices),
          100 * length(zero_indices) / length(training_BETA)))
        message("Removing these SNVs from beta method training data")
      }

      # Check if all BETA values are zero
      if (length(zero_indices) == length(training_BETA)) {
        stop("All training_BETA values are zero. ",
             "Cannot estimate effect sizes with beta method.")
      }

      # Remove zero-valued BETA and corresponding MAF values
      training_MAF_beta <- training_MAF[-zero_indices]
      training_BETA_beta <- training_BETA[-zero_indices]

      if (verbose >= 2) {
        message(sprintf(
          "Using %d non-zero BETA values for beta method (%.1f%% of original)",
          length(training_BETA_beta),
          100 * length(training_BETA_beta) / length(training_BETA)))
      }
    } else {
      # No zero values, use all data
      training_MAF_beta <- training_MAF
      training_BETA_beta <- training_BETA
    }

    # Call the renamed training helper
    result_beta <- .train_B_beta_method(
      training_MAF = training_MAF_beta,
      training_BETA = training_BETA_beta,
      selection_criterion = selection_criterion,
      custom_models = custom_models,
      show_model_selection = show_model_selection,
      verbose = verbose
    )
    model_beta_method <- result_beta$model
    all_models_beta_method <- result_beta$all_models_info
  }

  # Train p-value method model if requested
  if (method_to_use %in% c("pvalue", "both")) {
    if (verbose >= 2) {
      message("Training p-value method model...")
    }

    result_pvalue <- .train_B_pvalue_method(
      training_MAF = training_MAF,
      training_P = training_P,
      training_P_mlog10 = training_P_mlog10,
      training_N = training_N,
      selection_criterion = selection_criterion,
      custom_models = custom_models,
      show_model_selection = show_model_selection,
      verbose = verbose
    )
    model_pvalue_method <- result_pvalue$model
    all_models_pvalue_method <- result_pvalue$all_models_info
  }

  # ========== Training-Time Method Comparison ==========

  # When both methods are trained, compare model quality using the
  # selection criterion. This is training-time only: we compare model
  # fit metrics, NOT B predictions (which require target data).
  comparison <- list(
    selection_criterion = NULL,
    criterion_beta_method = NULL,
    criterion_pvalue_method = NULL,
    method_selected = NULL
  )

  if (method_to_use == "both") {
    if (verbose >= 2) {
      message("Comparing model quality to select primary method...")
    }

    # Extract criterion values for both models
    # Beta method: Y = BETA^2 ~ f(MAF)
    Y_beta <- training_BETA_beta^2
    criterion_beta <- .extract_model_criterion(
      model = model_beta_method,
      X = training_MAF_beta,
      Y = Y_beta,
      criterion = selection_criterion,
      verbose = verbose
    )

    # P-value method: Y = h^2 ~ f(MAF)
    Z2 <- .compute_chi2_from_pvalue(
      p = training_P, p_mlog10 = training_P_mlog10)
    Y_pvalue <- Z2 / (2 * training_N * training_MAF * (1 - training_MAF))
    criterion_pvalue <- .extract_model_criterion(
      model = model_pvalue_method,
      X = training_MAF,
      Y = Y_pvalue,
      criterion = selection_criterion,
      verbose = verbose
    )

    # Determine which model is better (higher is better for all criteria)
    beta_is_better <- criterion_beta > criterion_pvalue

    method_selected <- if (beta_is_better) "beta_method" else "pvalue_method"

    if (verbose >= 1) {
      message(sprintf(
        "Selected %s as primary (%s: beta_method=%.4f, pvalue_method=%.4f)",
        method_selected,
        selection_criterion,
        criterion_beta,
        criterion_pvalue
      ))
    }

    comparison <- list(
      selection_criterion = selection_criterion,
      criterion_beta_method = criterion_beta,
      criterion_pvalue_method = criterion_pvalue,
      method_selected = method_selected
    )
  }

  # ========== Compile Training Summary ==========

  # Calculate original number of training samples
  # Use original_training_MAF which was saved before any outlier removal
  n_original_training <- length(original_training_MAF)

  # Count total outliers detected (before any removal)
  if (method_to_use == "both") {
    n_outliers_detected <- length(unique(c(
      if (!is.null(outliers_beta_method)) {
        outliers_beta_method$indices
      } else {
        integer(0)
      },
      if (!is.null(outliers_pvalue_method)) {
        outliers_pvalue_method$indices
      } else {
        integer(0)
      }
    )))
  } else if (method_to_use == "beta" && !is.null(outliers_beta_method)) {
    n_outliers_detected <- length(outliers_beta_method$indices)
  } else if (method_to_use == "pvalue" && !is.null(outliers_pvalue_method)) {
    n_outliers_detected <- length(outliers_pvalue_method$indices)
  } else {
    n_outliers_detected <- 0
  }

  # Translate method_to_use to the _method naming convention
  method_used_name <- switch(method_to_use,
    beta = "beta_method",
    pvalue = "pvalue_method",
    both = "both"
  )

  # ========== Build glow_B_model Object ==========

  result <- structure(list(
    # Method identifier
    method_used = method_used_name,

    # Best fitted lm objects
    models = list(
      beta_method = model_beta_method,
      pvalue_method = model_pvalue_method
    ),

    # Complete model info for model comparison and diagnostics
    all_models_info = list(
      beta_method = all_models_beta_method,
      pvalue_method = all_models_pvalue_method
    ),

    # Outlier information
    outliers = list(
      method = outlier_method,
      action = outlier_action,
      beta_method = outliers_beta_method,
      pvalue_method = outliers_pvalue_method,
      indices_removed = if (outlier_action == "remove" &&
                            outlier_method != "none") {
        all_original <- seq_len(n_original_training)
        setdiff(all_original, training_indices_kept)
      } else {
        integer(0)
      }
    ),

    # Training data summary
    training_summary = list(
      n_original = n_original_training,
      n_used = length(training_MAF),
      n_outliers_detected = n_outliers_detected,
      trait_type = training_trait_normalized
    ),

    # Training data (retained for diagnostics and prediction)
    # Uses original (pre-outlier-removal) data so that outlier indices
    # stored in $outliers$beta_method$indices are consistent with these vectors
    training_data = list(
      MAF = original_training_MAF,
      BETA = original_training_BETA,
      P = original_training_P,
      P_mlog10 = original_training_P_mlog10,
      N = original_training_N,
      trait = training_trait_normalized
    ),

    # Method comparison (populated when method = "both")
    comparison = comparison,

    # Selection criterion used
    selection_criterion = selection_criterion

  ), class = "glow_B_model")

  return(result)
}


#' Print Method for glow_B_model Objects
#'
#' @description
#' Displays a concise summary of a trained B estimation model, including
#' method used, model formula, fit statistics, and training data summary.
#'
#' @param x A \code{glow_B_model} object from \code{\link{train_B_model}}
#' @param ... Additional arguments (ignored)
#'
#' @return The input object, invisibly
#'
#' @examples
#' \dontrun{
#' b_model <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = runif(50, 0.001, 0.3),
#'   training_BETA = rnorm(50, 0, 0.1),
#'   method = "beta"
#' )
#' print(b_model)
#' }
#'
#' @export
print.glow_B_model <- function(x, ...) {

  cat("\n")
  cat("========================================\n")
  cat("GLOW B Estimation Model (glow_B_model)\n")
  cat("========================================\n\n")

  # Method used
  method_display <- switch(x$method_used,
    beta_method = "Beta Method (BETA^2 ~ f(MAF))",
    pvalue_method = "P-value Method (h^2 ~ f(MAF))",
    both = "Both Methods"
  )
  cat("Method:              ", method_display, "\n")
  cat("Selection criterion: ", x$selection_criterion, "\n")

  # Training summary
  cat("\n--- Training Data ---\n")
  cat("  Trait type:         ", x$training_summary$trait_type, "\n")
  cat("  Variants (original):", x$training_summary$n_original, "\n")
  cat("  Variants (used):    ", x$training_summary$n_used, "\n")
  if (x$training_summary$n_outliers_detected > 0) {
    cat("  Outliers detected:  ", x$training_summary$n_outliers_detected, "\n")
    cat("  Outlier action:     ", x$outliers$action, "\n")
  }

  # Model details for each method
  .print_method_model_info <- function(method_name, lm_obj) {
    if (is.null(lm_obj)) return()

    display_name <- switch(method_name,
      beta_method = "Beta Method",
      pvalue_method = "P-value Method"
    )

    cat("\n--- ", display_name, " Model ---\n")

    # Formula
    model_terms <- all.vars(formula(lm_obj))
    formula_str <- paste(deparse(formula(lm_obj)), collapse = " ")
    cat("  Formula:  ", formula_str, "\n")

    # Model ID
    model_id <- attr(lm_obj, "model_id")
    if (!is.null(model_id)) {
      cat("  Model ID: ", model_id, "\n")
    }

    # Fit statistics
    summ <- summary(lm_obj)
    cat("  R-squared:", sprintf("%.4f", summ$r.squared), "\n")
    cat("  Adj R-sq: ", sprintf("%.4f", summ$adj.r.squared), "\n")
    cat("  Sigma:    ", sprintf("%.6g", summ$sigma), "\n")
  }

  # Print info for available methods
  if (!is.null(x$models$beta_method)) {
    .print_method_model_info("beta_method", x$models$beta_method)
  }
  if (!is.null(x$models$pvalue_method)) {
    .print_method_model_info("pvalue_method", x$models$pvalue_method)
  }

  # Comparison info (when both methods)
  if (x$method_used == "both" && !is.null(x$comparison$method_selected)) {
    cat("\n--- Method Comparison ---\n")
    cat("  Criterion:        ", x$comparison$selection_criterion, "\n")
    cat("  Beta method:      ",
        sprintf("%.4f", x$comparison$criterion_beta_method), "\n")
    cat("  P-value method:   ",
        sprintf("%.4f", x$comparison$criterion_pvalue_method), "\n")
    cat("  Selected:         ", x$comparison$method_selected, "\n")
  }

  cat("\n========================================\n\n")

  invisible(x)
}
