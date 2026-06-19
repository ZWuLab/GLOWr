########## Outlier Detection for B Estimation ##########
#
# This file contains functions for detecting outliers in GWAS training data
# used for B estimation. Supports both statistical (Cook's distance) and
# biological (implausible MAF-effect combinations) outlier detection methods.

#################### EXPORTED MAIN FUNCTIONS ####################

#' Detect Outliers in Training Data for B Estimation
#'
#' @description
#' Identifies outliers in GWAS training data using statistical (Cook's distance)
#' and/or biological (implausible MAF-effect combinations) criteria. This function
#' helps identify potentially problematic data points that may unduly influence
#' the B estimation model.
#'
#' @param MAF Numeric vector of minor allele frequencies (0 < MAF <= 0.5)
#' @param values Numeric vector of effect measures (BETA or h-squared values)
#' @param method Character string specifying detection method: "statistical"
#'   (Cook's distance only), "biological" (biological rules only), or "both"
#'   (default, applies both methods)
#' @param cook_threshold Numeric value for Cook's distance threshold multiplier.
#'   Outliers are defined as points with Cook's D > cook_threshold / n, where n
#'   is the number of observations. Default is 4.
#' @param biological_rules List of biological outlier rules. If NULL (default),
#'   uses default rule: common variants (MAF > 0.05) with large effects
#'   (|value| > 10). Custom rules can be provided as a list with elements
#'   containing 'maf_min' and 'effect_max' thresholds.
#' @param verbose Integer controlling verbosity level:
#'   0 = silent, 1 = warnings only, 2 = informational messages, 3 = detailed debug
#'
#' @return An S3 object of class "glow_outliers" containing:
#'   \describe{
#'     \item{indices}{Integer vector of all outlier indices (combined from both methods)}
#'     \item{statistical}{Integer vector of statistical outlier indices (Cook's distance)}
#'     \item{biological}{Integer vector of biological outlier indices}
#'     \item{reasons}{Named list of descriptions for each outlier}
#'     \item{outlier_MAF}{MAF values for detected outliers only (memory efficient)}
#'     \item{outlier_values}{Effect values for detected outliers only (memory efficient)}
#'     \item{n_total}{Total number of variants analyzed}
#'     \item{method_used}{Character string indicating which detection method was used}
#'   }
#'
#' @details
#' **Statistical Outliers:**
#' Uses Cook's distance from a simple linear regression (values ~ MAF) to identify
#' high-influence points. Cook's distance measures how much the fitted values
#' change when an observation is removed. The threshold is typically 4/n, where
#' n is the sample size.
#'
#' **Biological Outliers:**
#' Identifies variants with effect sizes that are biologically implausible given
#' their allele frequency. The default rule flags common variants (MAF > 0.05)
#' with very large effects (|value| > 10), as such combinations are rare in
#' complex traits and may indicate genotyping errors, population stratification,
#' or other artifacts.
#'
#' **Computational Complexity:**
#' O(n) for biological detection, O(n) for Cook's distance calculation via
#' leverages and residuals. Overall complexity is O(n) where n is the number
#' of variants.
#'
#' @examples
#' # Simulate training data with one clear outlier
#' set.seed(123)
#' MAF <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45)
#' values <- c(0.5, 0.3, 0.2, 0.15, 0.1, 15, 0.05)  # value[6] is outlier
#'
#' # Detect outliers using both methods
#' outliers <- detect_B_outliers(MAF, values, method = "both", verbose = 2)
#' print(outliers)
#'
#' # Use only statistical outlier detection
#' outliers_stat <- detect_B_outliers(MAF, values, method = "statistical")
#'
#' # Use custom biological rules
#' custom_rules <- list(
#'   common_large = list(maf_min = 0.1, effect_max = 5)
#' )
#' outliers_custom <- detect_B_outliers(MAF, values,
#'                                      method = "biological",
#'                                      biological_rules = custom_rules)
#'
#' @references
#' Cook, R. D. (1977). Detection of Influential Observations in Linear Regression.
#' Technometrics, 19(1), 15-22.
#'
#' @export
detect_B_outliers <- function(
  MAF,
  values,
  method = "both",
  cook_threshold = 4,
  biological_rules = NULL,
  verbose = 1
) {

  # Input validation
  # Check that inputs are numeric vectors
  if (!is.numeric(MAF) || !is.numeric(values)) {
    stop("MAF and values must be numeric vectors")
  }

  # Check for equal length
  if (length(MAF) != length(values)) {
    stop("MAF and values must have the same length")
  }

  # Check for non-empty inputs
  if (length(MAF) == 0) {
    stop("MAF and values cannot be empty")
  }

  # Check for NA values
  if (any(is.na(MAF)) || any(is.na(values))) {
    stop("MAF and values cannot contain NA values")
  }

  # Check for Inf/-Inf values which would cause calculation errors
  if (any(is.infinite(MAF)) || any(is.infinite(values))) {
    stop("MAF and values cannot contain Inf or -Inf values")
  }

  # Check MAF range (must be in (0, 0.5] for minor allele frequency)
  if (any(MAF <= 0) || any(MAF > 0.5)) {
    stop("MAF values must be in the range (0, 0.5]")
  }

  # Validate method parameter (must be one of the supported detection methods)
  if (!method %in% c("statistical", "biological", "both")) {
    stop("method must be one of: 'statistical', 'biological', 'both'")
  }

  # Validate Cook's distance threshold (must be positive for statistical detection)
  if (cook_threshold <= 0) {
    stop("cook_threshold must be positive")
  }

  # Initialize outlier lists to store indices of detected outliers
  outliers <- list()
  outliers$statistical <- integer(0)  # Will hold indices from Cook's distance method
  outliers$biological <- integer(0)   # Will hold indices from biological rules

  # Statistical outliers (Cook's distance)
  # Cook's distance measures influence: how much fitted values change when a point is removed
  if (method %in% c("statistical", "both")) {
    if (verbose >= 2) {
      message("Computing statistical outliers using Cook's distance...")
    }

    # Fit simple linear model for Cook's distance calculation
    # Use tryCatch in case model fitting fails (e.g., perfect collinearity, insufficient data)
    tryCatch({
      # Fit linear regression: effect size ~ MAF
      # This establishes the expected relationship between allele frequency and effect size
      simple_model <- lm(values ~ MAF)

      # Calculate Cook's distance for each observation
      # Cook's D combines leverage (extremeness of predictor) and residual (fit quality)
      cooks_d <- cooks.distance(simple_model)

      # Calculate threshold: common rule is 4/n for identifying influential points
      n <- length(values)
      threshold <- cook_threshold / n

      # Identify observations exceeding the Cook's distance threshold
      # These are high-influence points that strongly affect the regression line
      outliers$statistical <- which(cooks_d > threshold)

      # Report findings if verbose output is enabled
      if (length(outliers$statistical) > 0 && verbose >= 2) {
        message(sprintf("  Found %d statistical outlier(s) (Cook's D > %.4f)",
                       length(outliers$statistical), threshold))
      } else if (verbose >= 2) {
        message("  No statistical outliers detected")
      }
    }, error = function(e) {
      # If model fitting fails, warn user and continue without statistical outliers
      warning(sprintf("Failed to compute Cook's distance: %s. Skipping statistical outliers.", e$message))
      # outliers$statistical remains empty (initialized above)
    })
  }

  # Biological outliers
  # Flags variants with implausible MAF-effect combinations that suggest data artifacts
  if (method %in% c("biological", "both")) {
    if (verbose >= 2) {
      message("Detecting biological outliers...")
    }

    # Set default biological rules if none provided
    # Default rule: common variants (MAF > 0.05) should not have very large effects (> 10)
    # Such combinations are rare in complex traits and often indicate genotyping errors
    if (is.null(biological_rules)) {
      biological_rules <- list(
        common_large = list(maf_min = 0.05, effect_max = 10)
      )
    }

    # Pre-allocate list to collect outliers from each rule (avoids O(k²) vector growth)
    biological_outliers_list <- vector("list", length(biological_rules))

    # Apply each biological rule to identify implausible variants
    for (i in seq_along(biological_rules)) {
      rule_name <- names(biological_rules)[i]
      rule <- biological_rules[[rule_name]]

      # Check for common variant with large effect (if rule specifies both thresholds)
      if (!is.null(rule$maf_min) && !is.null(rule$effect_max)) {
        # Identify variants where MAF exceeds minimum AND absolute effect exceeds maximum
        # This flags common alleles with unexpectedly large effects
        idx <- which(MAF > rule$maf_min & abs(values) > rule$effect_max)

        # Store indices in pre-allocated list (efficient: no repeated concatenation)
        biological_outliers_list[[i]] <- idx

        # Report findings for this rule
        if (length(idx) > 0 && verbose >= 2) {
          message(sprintf("  Found %d biological outlier(s) from rule '%s' (MAF > %.2f, |effect| > %.1f)",
                         length(idx), rule_name, rule$maf_min, rule$effect_max))
        }
      }
    }

    # Combine all outliers from all rules and remove duplicates (single efficient operation)
    outliers$biological <- unique(unlist(biological_outliers_list))

    # Report overall biological outlier count
    if (length(outliers$biological) == 0 && verbose >= 2) {
      message("  No biological outliers detected")
    }
  }

  # Combine all outliers from both methods (remove duplicates)
  # A variant may be flagged by both statistical and biological methods
  all_outliers <- unique(c(outliers$statistical, outliers$biological))

  # Report total outlier count to user
  if (verbose >= 1 && length(all_outliers) > 0) {
    message(sprintf("Total of %d outlier(s) detected", length(all_outliers)))
  }

  # Generate human-readable descriptions explaining why each outlier was flagged
  descriptions <- .describe_outliers(MAF, values, outliers)

  # Create result object - store only outlier values for memory efficiency
  # For large datasets (e.g., 1M variants), this reduces memory from ~16MB to ~1-10KB
  # This is critical for genome-wide training data in get_B() pipeline
  result <- structure(
    list(
      indices = all_outliers,                # Combined outlier indices
      statistical = outliers$statistical,    # Subset: statistical outliers only
      biological = outliers$biological,      # Subset: biological outliers only
      reasons = descriptions,                # Explanation for each outlier
      outlier_MAF = MAF[all_outliers],      # Only store outlier MAF values (memory efficient)
      outlier_values = values[all_outliers], # Only store outlier effect values (memory efficient)
      n_total = length(MAF),                 # Total variants analyzed (for percentage calculations)
      method_used = method                   # Which detection method(s) were applied
    ),
    class = "glow_outliers"  # S3 class for print/summary methods
  )

  return(result)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' Generate Descriptions for Detected Outliers
#'
#' @description
#' Internal helper function to generate human-readable descriptions for each
#' detected outlier, indicating why it was flagged.
#'
#' @param MAF Numeric vector of minor allele frequencies
#' @param values Numeric vector of effect measures
#' @param outliers List containing statistical and biological outlier indices
#'
#' @return Named list of character strings describing each outlier
#'
#' @keywords internal
#' @noRd
.describe_outliers <- function(MAF, values, outliers) {

  descriptions <- list()

  # Describe statistical outliers
  if (length(outliers$statistical) > 0) {
    for (idx in outliers$statistical) {
      key <- as.character(idx)
      desc <- sprintf(
        "Statistical outlier: MAF=%.4f, value=%.4f (high influence)",
        MAF[idx], values[idx]
      )

      # If not already described (avoid overwriting biological description)
      if (is.null(descriptions[[key]])) {
        descriptions[[key]] <- desc
      } else {
        # Append to existing description
        descriptions[[key]] <- paste(descriptions[[key]], "; ", desc, sep = "")
      }
    }
  }

  # Describe biological outliers
  if (length(outliers$biological) > 0) {
    for (idx in outliers$biological) {
      key <- as.character(idx)
      desc <- sprintf(
        "Biological outlier: MAF=%.4f, value=%.4f (implausible effect for MAF)",
        MAF[idx], values[idx]
      )

      # If not already described
      if (is.null(descriptions[[key]])) {
        descriptions[[key]] <- desc
      } else {
        # Check if it's a statistical outlier too
        if (idx %in% outliers$statistical) {
          # Replace to show both
          descriptions[[key]] <- sprintf(
            "Statistical & Biological outlier: MAF=%.4f, value=%.4f (high influence and implausible effect)",
            MAF[idx], values[idx]
          )
        }
      }
    }
  }

  return(descriptions)
}

#################### S3 METHODS ####################

#' Print Method for glow_outliers Objects
#'
#' @description
#' Prints a summary of detected outliers in a human-readable format.
#'
#' @param x An object of class "glow_outliers"
#' @param max_show Maximum number of outliers to display in detail. Default is 10.
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.glow_outliers <- function(x, max_show = 10, ...) {

  cat("GLOW Outlier Detection Results\n")
  cat("==============================\n\n")

  # Summary counts
  cat(sprintf("Method used: %s\n", x$method_used))
  cat(sprintf("Total variants: %d\n", x$n_total))
  cat(sprintf("Total outliers: %d (%.1f%%)\n",
             length(x$indices),
             100 * length(x$indices) / x$n_total))

  if (x$method_used %in% c("statistical", "both")) {
    cat(sprintf("  Statistical outliers: %d\n", length(x$statistical)))
  }

  if (x$method_used %in% c("biological", "both")) {
    cat(sprintf("  Biological outliers: %d\n", length(x$biological)))
  }

  # Show detailed outlier information
  if (length(x$indices) > 0) {
    cat("\nDetailed Outlier Information:\n")
    cat("-----------------------------\n")

    # Determine how many to show
    n_show <- min(length(x$indices), max_show)

    for (i in 1:n_show) {
      idx <- x$indices[i]
      cat(sprintf("[%d] ", idx))
      cat(x$reasons[[as.character(idx)]], "\n")
    }

    # If there are more outliers than shown
    if (length(x$indices) > max_show) {
      cat(sprintf("\n... and %d more outlier(s). Use summary() for full list.\n",
                 length(x$indices) - max_show))
    }
  } else {
    cat("\nNo outliers detected.\n")
  }

  invisible(x)
}

#' Summary Method for glow_outliers Objects
#'
#' @description
#' Provides a comprehensive summary of detected outliers including detailed
#' statistics and all outlier descriptions.
#'
#' @param object An object of class "glow_outliers"
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns a summary data frame
#'
#' @export
summary.glow_outliers <- function(object, ...) {

  cat("GLOW Outlier Detection Summary\n")
  cat("==============================\n\n")

  # Basic information
  cat(sprintf("Detection method: %s\n", object$method_used))
  cat(sprintf("Total variants analyzed: %d\n", object$n_total))
  cat(sprintf("Total outliers detected: %d (%.1f%%)\n\n",
             length(object$indices),
             100 * length(object$indices) / object$n_total))

  # Method-specific counts
  if (object$method_used %in% c("statistical", "both")) {
    cat(sprintf("Statistical outliers (Cook's distance): %d\n",
               length(object$statistical)))
  }

  if (object$method_used %in% c("biological", "both")) {
    cat(sprintf("Biological outliers (implausible effects): %d\n",
               length(object$biological)))
  }

  # Overlap between methods (if both were used)
  if (object$method_used == "both") {
    overlap <- length(intersect(object$statistical, object$biological))
    cat(sprintf("Outliers detected by both methods: %d\n", overlap))
  }

  # Summary statistics for outliers
  if (length(object$indices) > 0) {
    cat("\nOutlier Statistics:\n")
    cat("-------------------\n")

    # Use stored outlier values (already subset)
    outlier_mafs <- object$outlier_MAF
    outlier_values <- object$outlier_values

    cat(sprintf("MAF range: %.4f - %.4f\n", min(outlier_mafs), max(outlier_mafs)))
    cat(sprintf("MAF mean: %.4f\n", mean(outlier_mafs)))
    cat(sprintf("Value range: %.4f - %.4f\n", min(outlier_values), max(outlier_values)))
    cat(sprintf("Value mean: %.4f\n", mean(outlier_values)))

    cat("\nComplete Outlier List:\n")
    cat("----------------------\n")

    # Create summary data frame using stored outlier values
    outlier_df <- data.frame(
      Index = object$indices,
      MAF = object$outlier_MAF,
      Value = object$outlier_values,
      Type = sapply(object$indices, function(idx) {
        is_stat <- idx %in% object$statistical
        is_bio <- idx %in% object$biological
        if (is_stat && is_bio) "Both"
        else if (is_stat) "Statistical"
        else "Biological"
      }),
      stringsAsFactors = FALSE
    )

    print(outlier_df, row.names = FALSE)

    cat("\nDetailed Reasons:\n")
    cat("-----------------\n")
    for (idx in object$indices) {
      cat(sprintf("[%d] %s\n", idx, object$reasons[[as.character(idx)]]))
    }

    # Return data frame invisibly
    invisible(outlier_df)
  } else {
    cat("\nNo outliers detected.\n")
    invisible(data.frame())
  }
}
