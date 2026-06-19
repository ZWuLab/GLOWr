########## User Interface Wrappers for get_B ##########
#
# This file contains simplified and advanced user interface wrappers for the
# get_B() function. These wrappers provide intuitive interfaces with sensible
# defaults for beginners (get_B_auto) and full control for advanced users
# (get_B_advanced).
#
# EXPORTED FUNCTIONS:
#   get_B_auto()     - Simplified interface with auto-detection
#   get_B_advanced() - Advanced interface with full control
#
# INTERNAL HELPER FUNCTIONS:
#   .print_outlier_details() - Print detailed outlier information

#################### EXPORTED MAIN FUNCTIONS ####################

#' Simplified Interface for B Estimation with Auto-Detection
#'
#' @description
#' A beginner-friendly wrapper for \code{\link{get_B}} that automatically
#' detects data format, selects appropriate methods, and uses sensible defaults.
#' This function minimizes the number of required parameters and handles most
#' common use cases automatically.
#'
#' @param training_data Either:
#'   \itemize{
#'     \item A data.frame with columns: MAF, BETA or (P/P_mlog10, N), and optionally rsID, CHR, POS.
#'           P_mlog10 represents -log10(p-value) and is preferred for extreme p-values.
#'     \item An object of class \code{glow_training_data} from \code{\link{prepare_B_training_data}}
#'   }
#' @param target_MAF Numeric vector of target MAF values for which to predict
#'   effect sizes. Values should be in (0, 1)
#' @param training_trait Character string specifying training trait type:
#'   "binary" or "continuous". If NULL and \code{training_data} is a
#'   \code{glow_training_data} object, automatically extracted. For data.frame
#'   input, this parameter is required
#' @param target_trait Character string specifying target trait type:
#'   "binary" or "continuous". If NULL (default), uses \code{training_trait}
#' @param target_SE Numeric value or vector of standard errors for target data.
#'   Required when target_trait == "continuous" and p-value method is used.
#'   Can be a single value (applied to all variants) or a vector matching
#'   length of \code{target_MAF}
#' @param target_case_prop Numeric value or vector specifying the proportion of
#'   cases in target case-control study. Required when target_trait == "binary"
#'   and p-value method is used. This is the study case proportion
#'   (cases / total sample size), not population prevalence. Can be a single
#'   value or a vector matching length of \code{target_MAF}
#' @param verbose Integer controlling verbosity: 0=silent, 1=warnings and
#'   summary (default), 2=detailed progress messages
#'
#' @return Numeric vector of predicted effect sizes (B) for each variant in
#'   \code{target_MAF}
#'
#' @details
#' \strong{Automatic Features:}
#'
#' This function automatically:
#' \itemize{
#'   \item Detects available data (BETA vs P-values) and chooses appropriate method
#'   \item Extracts trait type from \code{glow_training_data} objects
#'   \item Uses target_trait = training_trait if not specified
#'   \item Selects best model using R-squared criterion
#'   \item Detects both statistical and biological outliers (but keeps them)
#' }
#'
#' \strong{Default Settings:}
#' \itemize{
#'   \item method = "auto" (use both methods if same trait + both data available, else pvalue)
#'   \item selection_criterion = "R2" (maximize explained variance)
#'   \item outlier_method = "both" (detect statistical + biological outliers)
#'   \item outlier_action = "flag" (report but don't remove outliers)
#'   \item return_full = FALSE (return simple vector)
#' }
#'
#' \strong{When to Use:}
#'
#' Use this function when:
#' \itemize{
#'   \item You're new to the GLOW framework
#'   \item You want quick results with minimal setup
#'   \item You trust automatic method selection
#'   \item You don't need detailed diagnostics
#' }
#'
#' For more control, use \code{\link{get_B_advanced}} or \code{\link{get_B}}
#' directly.
#'
#' @section Computational Complexity:
#' O(8 * n_train + n_test) where n_train = number of training variants,
#' n_test = number of target variants
#'
#' @examples
#' \dontrun{
#' # Example 1: Using glow_training_data object
#' training_data <- prepare_B_training_data(
#'   data = "gwas_results.csv",
#'   trait_type = "binary"
#' )
#'
#' target_mafs <- seq(0.01, 0.4, by = 0.01)
#'
#' B_estimates <- get_B_auto(
#'   training_data = training_data,
#'   target_MAF = target_mafs,
#'   target_case_prop = 0.5
#' )
#'
#' # Example 2: Using data.frame directly
#' my_gwas <- data.frame(
#'   MAF = runif(100, 0.01, 0.5),
#'   BETA = rnorm(100, 0, 0.1),
#'   P = runif(100, 1e-8, 0.1),
#'   N = rep(5000, 100)
#' )
#'
#' B_estimates <- get_B_auto(
#'   training_data = my_gwas,
#'   training_trait = "continuous",
#'   target_MAF = seq(0.01, 0.4, by = 0.01),
#'   target_SE = 0.05  # Single value for all variants
#' )
#'
#' # Example 3: Cross-trait estimation (continuous to binary)
#' B_estimates <- get_B_auto(
#'   training_data = my_gwas,
#'   training_trait = "continuous",
#'   target_trait = "binary",
#'   target_MAF = seq(0.01, 0.4, by = 0.01),
#'   target_case_prop = 0.1
#' )
#' }
#'
#' @seealso
#' \code{\link{get_B}} for the main B estimation function
#' \code{\link{get_B_advanced}} for advanced interface with full control
#' \code{\link{prepare_B_training_data}} for preparing training data
#'
#' @export
get_B_auto <- function(
  training_data,
  target_MAF,
  training_trait = NULL,
  target_trait = NULL,
  target_SE = NULL,
  target_case_prop = NULL,
  verbose = 1
) {

  # ========== Input Validation ==========

  if (is.null(training_data)) {
    stop("training_data is required")
  }
  if (is.null(target_MAF)) {
    stop("target_MAF is required")
  }

  # ========== Extract Data from training_data ==========

  if (inherits(training_data, "glow_training_data")) {
    # Extract from glow_training_data object
    if (verbose >= 2) {
      message("Detected glow_training_data object")
    }

    # Extract trait type if not provided
    if (is.null(training_trait)) {
      training_trait <- training_data$metadata$trait_type
      if (verbose >= 2) {
        message("Auto-detected training trait type: ", training_trait)
      }
    }

    # Extract data columns
    data_df <- training_data$data

  } else if (is.data.frame(training_data)) {
    # Use data.frame directly
    if (verbose >= 2) {
      message("Using data.frame input")
    }

    if (is.null(training_trait)) {
      stop("training_trait must be specified when using data.frame input. ",
           "Use prepare_B_training_data() for automatic trait detection.")
    }

    data_df <- training_data

  } else {
    stop("training_data must be either a data.frame or glow_training_data object")
  }

  # ========== Extract Required Columns ==========

  # Check for required columns
  if (!"MAF" %in% names(data_df)) {
    stop("training_data must contain a 'MAF' column. ",
         "Use prepare_B_training_data() or standardize_column_names() to prepare your data.")
  }

  training_MAF <- data_df$MAF

  # Check for BETA column
  has_beta <- "BETA" %in% names(data_df)
  training_BETA <- if (has_beta) data_df$BETA else NULL

  # Check for P-value columns (either P or P_mlog10) and N columns
  has_p_regular <- "P" %in% names(data_df)
  has_p_mlog10 <- "P_mlog10" %in% names(data_df)
  has_n <- "N" %in% names(data_df)

  # Check if we have p-value data in either format
  has_pvalue <- (has_p_regular || has_p_mlog10) && has_n

  # Extract the appropriate p-value column
  training_P <- NULL
  training_P_mlog10 <- NULL
  if (has_p_mlog10) {
    # Prefer P_mlog10 for better numerical stability with extreme p-values
    training_P_mlog10 <- data_df$P_mlog10
    if (verbose >= 2) {
      message("Using P_mlog10 column for p-value method (better numerical stability)")
    }
  } else if (has_p_regular) {
    training_P <- data_df$P
  }

  training_N <- if (has_n) data_df$N else NULL

  # Validate that we have at least one method's data
  if (!has_beta && !has_pvalue) {
    stop("training_data must contain either BETA column or P/P_mlog10 and N columns. ",
         "Use prepare_B_training_data() to prepare your data.")
  }

  # ========== Determine Target Trait ==========

  if (is.null(target_trait)) {
    target_trait <- training_trait
    if (verbose >= 2) {
      message("Using target_trait = training_trait: ", target_trait)
    }
  }

  # ========== Validate Target-Specific Parameters ==========

  # For p-value method, check if we need target_SE or target_case_prop
  if (has_pvalue) {
    if (target_trait == "continuous" && is.null(target_SE)) {
      if (training_trait == target_trait && has_beta) {
        # We have beta method available, so pvalue method is optional
        if (verbose >= 1) {
          message("Note: target_SE not provided. Will use direct beta method only.")
        }
      } else {
        stop("target_SE is required when target_trait is 'continuous' and using p-value method. ",
             "Provide a single value or a vector matching length of target_MAF.")
      }
    }

    if (target_trait == "binary" && is.null(target_case_prop)) {
      if (training_trait == target_trait && has_beta) {
        # We have beta method available, so pvalue method is optional
        if (verbose >= 1) {
          message("Note: target_case_prop not provided. Will use direct beta method only.")
        }
      } else {
        stop("target_case_prop is required when target_trait is 'binary' and using p-value method. ",
             "Provide the proportion of cases in your case-control study (e.g., 0.5 for balanced design).")
      }
    }
  }

  # ========== Call get_B with Auto Settings ==========

  if (verbose >= 1) {
    message("Running B estimation with automatic settings...")
  }

  B_result <- get_B(
    training_trait = training_trait,
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = target_trait,
    target_MAF = target_MAF,
    target_SE = target_SE,
    target_case_prop = target_case_prop,
    method = "auto",
    selection_criterion = "R2",
    custom_models = NULL,
    outlier_method = "both",
    outlier_action = "flag",
    cook_threshold = 4,
    biological_rules = NULL,
    return_full = FALSE,
    verbose = verbose
  )

  return(B_result)
}


#' Advanced Interface for B Estimation with Full Control
#'
#' @description
#' An advanced wrapper for \code{\link{get_B}} that provides full control over
#' all parameters including cross-validation, outlier handling, and diagnostic
#' integration. This function is designed for users who need fine-grained
#' control over the estimation process.
#'
#' @param training_data Either:
#'   \itemize{
#'     \item A data.frame with columns: MAF, BETA or (P/P_mlog10, N), and optionally rsID, CHR, POS.
#'           P_mlog10 represents -log10(p-value) and is preferred for extreme p-values.
#'     \item An object of class \code{glow_training_data} from \code{\link{prepare_B_training_data}}
#'   }
#' @param target_MAF Numeric vector of target MAF values for which to predict
#'   effect sizes. Values should be in (0, 1)
#' @param training_trait Character string specifying training trait type:
#'   "binary" or "continuous". If NULL and \code{training_data} is a
#'   \code{glow_training_data} object, automatically extracted
#' @param target_trait Character string specifying target trait type:
#'   "binary" or "continuous". If NULL (default), uses \code{training_trait}
#' @param target_SE Numeric value or vector of standard errors for target data.
#'   Required when target_trait == "continuous" and p-value method is used
#' @param target_case_prop Numeric value or vector specifying the proportion of
#'   cases in target case-control study. Required when target_trait == "binary"
#'   and p-value method is used
#' @param method Character string specifying which method(s) to use:
#'   \itemize{
#'     \item "auto" (default): Automatic selection based on available data
#'     \item "beta": Direct beta method only
#'     \item "pvalue": P-value/Z-score method only
#'     \item "both": Run both methods for comparison
#'   }
#' @param selection_criterion Character string specifying model selection
#'   criterion: "R2" (default), "adj_R2", or "CV_R2" (cross-validation).
#'   When \code{cross_validation = TRUE}, this is automatically set to "CV_R2"
#' @param custom_models Optional list of custom model formulas to include in
#'   candidate model set. Each element should be a formula object
#' @param outlier_method Character string specifying outlier detection method:
#'   "none", "statistical" (Cook's distance), "biological" (implausible
#'   combinations), or "both" (default)
#' @param outlier_action Character string specifying action for detected outliers:
#'   "flag" (default, report but don't remove) or "remove" (exclude from analysis)
#' @param cook_threshold Numeric threshold multiplier for Cook's distance
#'   outlier detection (default: 4). Outliers are defined as points with
#'   Cook's D > cook_threshold/n
#' @param biological_rules Optional list of custom rules for biological outlier
#'   detection. Default uses common variants (MAF > 0.05) with large effects (|effect| > 10)
#' @param cross_validation Logical. If TRUE, enables leave-one-out cross-validation
#'   (LOOCV) for model selection and automatically sets \code{selection_criterion = "CV_R2"}.
#'   Note: Currently implements LOOCV only; k-fold CV support is planned for a future release
#' @param show_diagnostics Logical. If TRUE, automatically displays diagnostic
#'   plots after estimation (default: FALSE)
#' @param show_outlier_details Logical. If TRUE, prints detailed outlier
#'   information (default: FALSE)
#' @param return_full Logical. If FALSE, returns only B estimates. If TRUE
#'   (default), returns comprehensive result object with diagnostics
#' @param verbose Integer controlling verbosity: 0=silent, 1=warnings,
#'   2=detailed progress (default), 3=debug output
#'
#' @return If \code{return_full = FALSE}: Numeric vector of predicted effect
#'   sizes (B) for each variant in \code{target_MAF}.
#'
#'   If \code{return_full = TRUE} (default): A list with class "glow_B_estimate"
#'   containing comprehensive results including B estimates, fitted models,
#'   outlier detection results, method comparison (if applicable), and
#'   diagnostics. See \code{\link{get_B}} for full details.
#'
#' @details
#' \strong{Advanced Features:}
#'
#' This function provides:
#' \itemize{
#'   \item Full control over estimation method selection
#'   \item Customizable model selection criteria (R2, adj_R2, CV_R2)
#'   \item Cross-validation for robust model selection
#'   \item Flexible outlier detection and handling
#'   \item Automatic diagnostic display
#'   \item Method comparison when both methods applicable
#'   \item Custom model specification
#' }
#'
#' \strong{Cross-Validation:}
#'
#' When \code{cross_validation = TRUE}:
#' \itemize{
#'   \item Performs leave-one-out cross-validation (LOOCV)
#'   \item Automatically sets \code{selection_criterion = "CV_R2"}
#'   \item Provides more robust model selection for small samples
#'   \item Increases computation time by factor of n_train
#'   \item Note: K-fold CV support is planned for a future release
#' }
#'
#' \strong{Diagnostic Integration:}
#'
#' When \code{show_diagnostics = TRUE}:
#' \itemize{
#'   \item Displays 2x2 dashboard of diagnostic plots
#'   \item Shows model fit, residuals, and method comparison
#'   \item Automatically called after estimation completes
#' }
#'
#' When \code{show_outlier_details = TRUE}:
#' \itemize{
#'   \item Prints detailed information about detected outliers
#'   \item Shows Cook's distance values
#'   \item Reports biological rule violations
#' }
#'
#' \strong{When to Use:}
#'
#' Use this function when:
#' \itemize{
#'   \item You need full control over estimation parameters
#'   \item You want to compare different methods
#'   \item You need cross-validation for small samples
#'   \item You want comprehensive diagnostics
#'   \item You're developing or validating methods
#' }
#'
#' For simpler use cases, use \code{\link{get_B_auto}}.
#'
#' @section Computational Complexity:
#' O(k * 8 * n_train + n_test) where k = n_train for LOOCV (1 if no CV),
#' n_train = number of training variants, n_test = number of target variants.
#' LOOCV evaluates all 8 candidate models n_train times.
#'
#' @examples
#' \dontrun{
#' # Example 1: Advanced estimation with cross-validation
#' result <- get_B_advanced(
#'   training_data = prepared_data,
#'   target_MAF = seq(0.01, 0.4, by = 0.01),
#'   target_case_prop = 0.5,
#'   method = "both",
#'   cross_validation = TRUE,  # Uses LOOCV
#'   show_diagnostics = TRUE,
#'   verbose = 2
#' )
#'
#' # Access results
#' B_primary <- result$B
#' comparison <- result$model$comparison
#' model_stats <- summary(result$model$models$beta_method)
#'
#' # Example 2: Custom outlier handling
#' result <- get_B_advanced(
#'   training_data = my_data,
#'   training_trait = "binary",
#'   target_MAF = target_mafs,
#'   target_case_prop = 0.5,
#'   outlier_method = "both",
#'   outlier_action = "remove",  # Actually remove outliers
#'   cook_threshold = 3,  # More stringent threshold
#'   show_outlier_details = TRUE
#' )
#'
#' # Example 3: Custom model selection
#' custom_formulas <- list(
#'   formula(Y ~ poly(X, 3)),
#'   formula(Y ~ poly(X, 4))
#' )
#'
#' result <- get_B_advanced(
#'   training_data = my_data,
#'   training_trait = "continuous",
#'   target_MAF = target_mafs,
#'   target_SE = 0.1,
#'   custom_models = custom_formulas,
#'   selection_criterion = "adj_R2",  # Prefer simpler models
#'   return_full = TRUE
#' )
#'
#' # Example 4: Compare all selection criteria
#' criteria <- c("R2", "adj_R2", "CV_R2")
#' results <- lapply(criteria, function(crit) {
#'   get_B_advanced(
#'     training_data = my_data,
#'     training_trait = "binary",
#'     target_MAF = target_mafs,
#'     target_case_prop = 0.5,
#'     selection_criterion = crit,
#'     return_full = TRUE,
#'     verbose = 0
#'   )
#' })
#' names(results) <- criteria
#' }
#'
#' @seealso
#' \code{\link{get_B}} for the main B estimation function
#' \code{\link{get_B_auto}} for simplified automatic interface
#' \code{\link{get_B_diagnostics}} for detailed diagnostic analysis
#' \code{\link{plot.glow_B_estimate}} for visualization
#'
#' @export
get_B_advanced <- function(
  training_data,
  target_MAF,
  training_trait = NULL,
  target_trait = NULL,
  target_SE = NULL,
  target_case_prop = NULL,
  method = "both",
  selection_criterion = "R2",
  custom_models = NULL,
  outlier_method = "both",
  outlier_action = "flag",
  cook_threshold = 4,
  biological_rules = NULL,
  cross_validation = FALSE,
  show_diagnostics = FALSE,
  show_outlier_details = FALSE,
  return_full = TRUE,
  verbose = 2
) {

  # ========== Input Validation ==========

  if (is.null(training_data)) {
    stop("training_data is required")
  }
  if (is.null(target_MAF)) {
    stop("target_MAF is required")
  }

  # ========== Handle Cross-Validation Setting ==========

  if (cross_validation) {
    if (verbose >= 2) {
      message("Cross-validation enabled (LOOCV): setting selection_criterion = 'CV_R2'")
    }
    selection_criterion <- "CV_R2"
  }

  # ========== Extract Data from training_data ==========

  if (inherits(training_data, "glow_training_data")) {
    # Extract from glow_training_data object
    if (verbose >= 2) {
      message("Detected glow_training_data object")
    }

    # Extract trait type if not provided
    if (is.null(training_trait)) {
      training_trait <- training_data$metadata$trait_type
      if (verbose >= 2) {
        message("Auto-detected training trait type: ", training_trait)
      }
    }

    # Extract data columns
    data_df <- training_data$data

  } else if (is.data.frame(training_data)) {
    # Use data.frame directly
    if (verbose >= 2) {
      message("Using data.frame input")
    }

    if (is.null(training_trait)) {
      stop("training_trait must be specified when using data.frame input")
    }

    data_df <- training_data

  } else {
    stop("training_data must be either a data.frame or glow_training_data object")
  }

  # ========== Extract Required Columns ==========

  # Check for required columns
  if (!"MAF" %in% names(data_df)) {
    stop("training_data must contain a 'MAF' column")
  }

  training_MAF <- data_df$MAF

  # Check for BETA column
  has_beta <- "BETA" %in% names(data_df)
  training_BETA <- if (has_beta) data_df$BETA else NULL

  # Check for P-value columns (either P or P_mlog10) and N columns
  has_p_regular <- "P" %in% names(data_df)
  has_p_mlog10 <- "P_mlog10" %in% names(data_df)
  has_n <- "N" %in% names(data_df)

  # Check if we have p-value data in either format
  has_pvalue <- (has_p_regular || has_p_mlog10) && has_n

  # Extract the appropriate p-value column
  training_P <- NULL
  training_P_mlog10 <- NULL
  if (has_p_mlog10) {
    # Prefer P_mlog10 for better numerical stability with extreme p-values
    training_P_mlog10 <- data_df$P_mlog10
    if (verbose >= 2) {
      message("Using P_mlog10 column for p-value method (better numerical stability)")
    }
  } else if (has_p_regular) {
    training_P <- data_df$P
  }

  training_N <- if (has_n) data_df$N else NULL

  # Validate method compatibility with available data
  if (method == "beta" && !has_beta) {
    stop("method = 'beta' requires BETA column in training_data")
  }
  if (method == "pvalue" && !has_pvalue) {
    stop("method = 'pvalue' requires P/P_mlog10 and N columns in training_data")
  }
  if (method == "both" && (!has_beta || !has_pvalue)) {
    stop("method = 'both' requires BETA, P/P_mlog10, and N columns in training_data")
  }

  # ========== Determine Target Trait ==========

  if (is.null(target_trait)) {
    target_trait <- training_trait
    if (verbose >= 2) {
      message("Using target_trait = training_trait: ", target_trait)
    }
  }

  # ========== Call get_B ==========

  if (verbose >= 1) {
    message("Running B estimation with advanced settings...")
    if (method == "both") {
      message("  Running both methods for comparison")
    }
    if (cross_validation) {
      message("  Using leave-one-out cross-validation (LOOCV)")
    }
  }

  B_result <- get_B(
    training_trait = training_trait,
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    training_P = training_P,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    target_trait = target_trait,
    target_MAF = target_MAF,
    target_SE = target_SE,
    target_case_prop = target_case_prop,
    method = method,
    selection_criterion = selection_criterion,
    custom_models = custom_models,
    outlier_method = outlier_method,
    outlier_action = outlier_action,
    cook_threshold = cook_threshold,
    biological_rules = biological_rules,
    return_full = return_full,
    verbose = verbose
  )

  # ========== Show Outlier Details if Requested ==========

  if (show_outlier_details && return_full && !is.null(B_result$model$outliers)) {
    .print_outlier_details(B_result$model$outliers, verbose = verbose)
  }

  # ========== Show Diagnostics if Requested ==========

  if (show_diagnostics && return_full) {
    if (verbose >= 1) {
      message("\nDisplaying diagnostic plots...")
    }
    plot(B_result, type = "all")
  }

  return(B_result)
}


#################### INTERNAL HELPER FUNCTIONS ####################


#' Print Detailed Outlier Information
#'
#' @description
#' Internal helper function to print detailed information about detected outliers.
#'
#' @param outliers List of outlier detection results from get_B
#' @param verbose Verbosity level
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#' @noRd
.print_outlier_details <- function(outliers, verbose = 1) {

  cat("\n")
  cat("========================================\n")
  cat("Outlier Detection Details\n")
  cat("========================================\n\n")

  cat("Detection method: ", outliers$method, "\n")
  cat("Action taken:     ", outliers$action, "\n\n")

  # Beta method outliers
  if (!is.null(outliers$beta_method)) {
    cat("Direct Beta Method:\n")
    cat("-------------------\n")

    if (length(outliers$beta_method$indices) == 0) {
      cat("  No outliers detected\n\n")
    } else {
      cat(sprintf("  Found %d outliers\n", length(outliers$beta_method$indices)))
      cat("  Indices: ", paste(head(outliers$beta_method$indices, 20), collapse = ", "))
      if (length(outliers$beta_method$indices) > 20) {
        cat(", ... (", length(outliers$beta_method$indices) - 20, " more)")
      }
      cat("\n\n")

      # Show details for first few outliers
      if (!is.null(outliers$beta_method$details)) {
        cat("  Details for first 5 outliers:\n")
        print(head(outliers$beta_method$details, 5))
        cat("\n")
      }
    }
  }

  # P-value method outliers
  if (!is.null(outliers$pvalue_method)) {
    cat("P-value/Z-score Method:\n")
    cat("-----------------------\n")

    if (length(outliers$pvalue_method$indices) == 0) {
      cat("  No outliers detected\n\n")
    } else {
      cat(sprintf("  Found %d outliers\n", length(outliers$pvalue_method$indices)))
      cat("  Indices: ", paste(head(outliers$pvalue_method$indices, 20), collapse = ", "))
      if (length(outliers$pvalue_method$indices) > 20) {
        cat(", ... (", length(outliers$pvalue_method$indices) - 20, " more)")
      }
      cat("\n\n")

      # Show details for first few outliers
      if (!is.null(outliers$pvalue_method$details)) {
        cat("  Details for first 5 outliers:\n")
        print(head(outliers$pvalue_method$details, 5))
        cat("\n")
      }
    }
  }

  # Indices removed
  if (outliers$action == "remove" && length(outliers$indices_removed) > 0) {
    cat("Removed from analysis:\n")
    cat("----------------------\n")
    cat(sprintf("  %d training variants removed\n", length(outliers$indices_removed)))
    cat("  Indices: ", paste(head(outliers$indices_removed, 20), collapse = ", "))
    if (length(outliers$indices_removed) > 20) {
      cat(", ... (", length(outliers$indices_removed) - 20, " more)")
    }
    cat("\n")
  }

  cat("\n========================================\n\n")

  invisible(NULL)
}
