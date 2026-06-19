#' Helper Functions for GLOWr Validation Framework
#'
#' This file provides helper functions for validating GLOWr implementations
#' against the legacy GLOW package. These functions are used in test files
#' to ensure numerical accuracy and consistency with the original implementation.
#'
#' File Log (reverse chronological order):
#' - 2025-10-18: Created by r-developer - Initial validation framework

# Load required packages for validation
# suppressPackageStartupMessages({
#   if (!requireNamespace("glmnet", quietly = TRUE)) stop("glmnet required for validation")
#   if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("mvtnorm required for validation")
#   if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix required for validation")
# })

# Global variable to store loaded legacy functions
.legacy_env <- new.env()

#' Load Legacy GLOW Package Functions
#'
#' Sources all R files from the legacy GLOW package directory into a
#' separate environment to avoid namespace conflicts. This allows
#' comparison testing against the original implementation.
#'
#' @param legacy_dir Character string. Path to legacy GLOW/R directory.
#'   Defaults to \code{$GLOW_LEGACY_ROOT/legacy-materials/code/GLOW_R_pacakge/GLOW/R};
#'   set the \code{GLOW_LEGACY_ROOT} environment variable to enable the comparison.
#' @param verbose Logical. If TRUE, prints messages about loaded files.
#'
#' @return An environment containing all loaded legacy functions.
#'   The environment is also stored in `.legacy_env` for reuse.
#'
#' @details
#' The function sources all .R files in alphabetical order to handle
#' dependencies correctly. Functions are loaded into a separate environment
#' to prevent conflicts with GLOWr implementations.
#'
#' Dependencies loaded:
#' - get_B.R: Effect size estimation functions
#' - get_PI.R: Variant-importance score estimation functions
#' - getZ_marg_score.R: Marginal score calculation (continuous)
#' - getZ_marg_score_binary_SPA.R: Marginal score (binary, SPA-adjusted)
#' - GLOW_Burden.R: Burden test implementation
#' - GLOW_SKAT.R: SKAT test implementation
#' - GLOW_Fisher.R: Fisher test implementation
#' - GLOW_Omni.R: Omnibus test implementation
#' - GLOW_Omni_byP.R: Omnibus test by p-values
#' - helpers_GFisher.R: Helper functions for generalized Fisher test
#' - helpers_GLOWtests.R: Helper functions for GLOW tests
#' - helpers_optimalWeights.R: Helper functions for optimal weight calculation
#'
#' @examples
#' \dontrun{
#' legacy_env <- load_legacy_glow()
#' # Access legacy function
#' result <- legacy_env$getZ_marg_score(G, X, Y)
#' }
#'
#' @export
load_legacy_glow <- function(
    legacy_dir = file.path(Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"),
                           "legacy-materials/code/GLOW_R_pacakge/GLOW/R"),
    verbose = FALSE
) {
  # Verify directory exists
  if (!dir.exists(legacy_dir)) {
    stop("Legacy directory not found: ", legacy_dir)
  }

  # Get all R files in alphabetical order
  r_files <- list.files(
    path = legacy_dir,
    pattern = "\\.R$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(r_files) == 0) {
    stop("No R files found in legacy directory: ", legacy_dir)
  }

  # Sort files to ensure consistent loading order
  # Load helpers first to resolve dependencies
  priority_files <- c(
    "helpers_optimalWeights.R",
    "helpers_GFisher.R",
    "helpers_GLOWtests.R",
    "get_B.R",
    "get_PI.R",
    "getZ_marg_score.R",
    "getZ_marg_score_binary_SPA.R"
  )

  # Reorder files: priority files first, then others
  ordered_files <- character(0)
  for (pf in priority_files) {
    matching <- r_files[grepl(pf, r_files, fixed = TRUE)]
    if (length(matching) > 0) {
      ordered_files <- c(ordered_files, matching)
      r_files <- setdiff(r_files, matching)
    }
  }
  ordered_files <- c(ordered_files, r_files)

  # Source each file into the legacy environment
  for (file in ordered_files) {
    if (verbose) {
      message("Loading legacy file: ", basename(file))
    }
    tryCatch({
      source(file, local = .legacy_env)
    }, error = function(e) {
      warning("Error loading ", basename(file), ": ", e$message)
    })
  }

  if (verbose) {
    message("Loaded ", length(ordered_files), " legacy R files")
    message("Available functions: ",
            paste(head(ls(.legacy_env), 10), collapse = ", "),
            if (length(ls(.legacy_env)) > 10) "..." else "")
  }

  return(.legacy_env)
}


#' Compare GLOWr Results with Legacy Results
#'
#' Compares outputs from GLOWr and legacy GLOW implementations,
#' checking for numerical equivalence within a specified tolerance.
#'
#' @param glowr_result Result object from GLOWr function.
#' @param legacy_result Result object from legacy GLOW function.
#' @param tolerance Numeric. Maximum allowed difference for numerical comparison.
#'   Default is 1e-10.
#' @param type Character. Type of comparison: "numeric", "list", "matrix", "vector".
#'   If NULL (default), automatically detected.
#' @param test_name Character. Optional name for the test (used in reports).
#'
#' @return A list with components:
#'   \item{passed}{Logical. TRUE if comparison passed within tolerance.}
#'   \item{max_diff}{Numeric. Maximum absolute difference found.}
#'   \item{mean_diff}{Numeric. Mean absolute difference.}
#'   \item{details}{Character. Detailed comparison information.}
#'   \item{mismatches}{List. Elements that differ beyond tolerance (if any).}
#'
#' @details
#' The function handles different data types:
#' - Numeric vectors/matrices: Element-wise comparison
#' - Lists: Recursive comparison of matching elements
#' - NULL values: Checked for presence in both
#' - Mixed types: Reports type mismatch
#'
#' For list comparisons, the function recursively compares matching
#' elements by name. Unnamed elements are compared by position.
#'
#' @examples
#' \dontrun{
#' glowr_res <- GLOWr::getZ_marg_score(G, X, Y)
#' legacy_res <- .legacy_env$getZ_marg_score(G, X, Y)
#' comparison <- compare_with_legacy(glowr_res, legacy_res)
#' if (!comparison$passed) {
#'   print(comparison$details)
#' }
#' }
#'
#' @export
compare_with_legacy <- function(
    glowr_result,
    legacy_result,
    tolerance = 1e-10,
    type = NULL,
    test_name = NULL
) {
  # Initialize result list
  result <- list(
    passed = TRUE,
    max_diff = 0,
    mean_diff = 0,
    details = "",
    mismatches = list(),
    test_name = test_name
  )

  # Auto-detect type if not specified
  if (is.null(type)) {
    if (is.list(glowr_result) && !is.data.frame(glowr_result)) {
      type <- "list"
    } else if (is.matrix(glowr_result)) {
      type <- "matrix"
    } else if (is.numeric(glowr_result)) {
      type <- "vector"
    } else {
      type <- "other"
    }
  }

  # Check for NULL
  if (is.null(glowr_result) && is.null(legacy_result)) {
    result$details <- "Both results are NULL"
    return(result)
  }

  if (is.null(glowr_result) || is.null(legacy_result)) {
    result$passed <- FALSE
    result$details <- sprintf(
      "One result is NULL: GLOWr=%s, Legacy=%s",
      is.null(glowr_result), is.null(legacy_result)
    )
    return(result)
  }

  # Type checking
  if (class(glowr_result)[1] != class(legacy_result)[1]) {
    result$passed <- FALSE
    result$details <- sprintf(
      "Type mismatch: GLOWr=%s, Legacy=%s",
      class(glowr_result)[1], class(legacy_result)[1]
    )
    return(result)
  }

  # Comparison based on type
  if (type == "list") {
    # List comparison
    names_glowr <- names(glowr_result)
    names_legacy <- names(legacy_result)

    # Check for matching names
    if (!identical(sort(names_glowr), sort(names_legacy))) {
      result$passed <- FALSE
      result$details <- sprintf(
        "List names differ. GLOWr: %s, Legacy: %s",
        paste(names_glowr, collapse = ", "),
        paste(names_legacy, collapse = ", ")
      )
      return(result)
    }

    # Compare each element recursively
    all_diffs <- numeric(0)
    for (nm in names_glowr) {
      elem_comp <- compare_with_legacy(
        glowr_result[[nm]],
        legacy_result[[nm]],
        tolerance = tolerance,
        test_name = paste0(test_name, "$", nm)
      )

      if (!elem_comp$passed) {
        result$passed <- FALSE
        result$mismatches[[nm]] <- elem_comp
      }

      all_diffs <- c(all_diffs, elem_comp$max_diff)
    }

    result$max_diff <- max(all_diffs, na.rm = TRUE)
    result$mean_diff <- mean(all_diffs, na.rm = TRUE)
    result$details <- sprintf(
      "List with %d elements: max_diff=%.2e, mean_diff=%.2e",
      length(names_glowr), result$max_diff, result$mean_diff
    )

  } else if (type %in% c("vector", "matrix")) {
    # Numeric comparison
    if (!is.numeric(glowr_result) || !is.numeric(legacy_result)) {
      result$passed <- FALSE
      result$details <- "Non-numeric values in numeric comparison"
      return(result)
    }

    # Dimension checking
    if (!identical(dim(glowr_result), dim(legacy_result))) {
      result$passed <- FALSE
      result$details <- sprintf(
        "Dimension mismatch: GLOWr=%s, Legacy=%s",
        paste(dim(glowr_result), collapse = "x"),
        paste(dim(legacy_result), collapse = "x")
      )
      return(result)
    }

    # Calculate differences
    abs_diff <- abs(glowr_result - legacy_result)
    result$max_diff <- max(abs_diff, na.rm = TRUE)
    result$mean_diff <- mean(abs_diff, na.rm = TRUE)

    # Check tolerance
    if (result$max_diff > tolerance) {
      result$passed <- FALSE

      # Find where max difference occurs
      if (is.matrix(glowr_result)) {
        max_idx <- which(abs_diff == result$max_diff, arr.ind = TRUE)[1, ]
        result$details <- sprintf(
          "Max difference %.2e at [%d,%d] exceeds tolerance %.2e",
          result$max_diff, max_idx[1], max_idx[2], tolerance
        )
      } else {
        max_idx <- which.max(abs_diff)
        result$details <- sprintf(
          "Max difference %.2e at position %d exceeds tolerance %.2e",
          result$max_diff, max_idx, tolerance
        )
      }

      # Store mismatch details
      result$mismatches$max_location <- max_idx
      result$mismatches$glowr_value <- glowr_result[max_idx]
      result$mismatches$legacy_value <- legacy_result[max_idx]
    } else {
      result$details <- sprintf(
        "Numerical comparison passed: max_diff=%.2e, mean_diff=%.2e",
        result$max_diff, result$mean_diff
      )
    }

  } else {
    # Other types - use identical check
    if (!identical(glowr_result, legacy_result)) {
      result$passed <- FALSE
      result$details <- "Results not identical (non-numeric comparison)"
    } else {
      result$details <- "Results identical"
    }
  }

  return(result)
}


#' Generate Test Data for Validation
#'
#' Creates simulated genetic data with known properties for testing
#' GLOW methods. Includes options for continuous/binary outcomes,
#' correlated variants, and rare variants.
#'
#' @param n Integer. Sample size.
#' @param p Integer. Number of genetic variants.
#' @param k Integer. Number of covariates (excluding intercept). Default is 2.
#' @param binary Logical. If TRUE, generates binary outcome. Default is FALSE.
#' @param maf Numeric vector of length p. Minor allele frequencies.
#'   If NULL, generated uniformly between 0.05 and 0.5.
#' @param rare Logical. If TRUE, uses rare variant MAF (0.001 to 0.01).
#' @param ld_structure Logical. If TRUE, introduces LD structure among variants.
#' @param ld_strength Numeric. Correlation strength for LD (0 to 1). Default 0.3.
#' @param n_causal Integer. Number of causal variants. Default is p/5.
#' @param effect_size Numeric. Effect size for causal variants. Default is 0.5.
#' @param seed Integer. Random seed for reproducibility. Default is NULL.
#' @param prevalence Numeric. Disease prevalence for binary outcome (0 to 1).
#'   Default is 0.3.
#'
#' @return A list with components:
#'   \item{G}{Matrix (n x p). Genotype matrix with values 0, 1, 2.}
#'   \item{X}{Matrix (n x (k+1)). Covariate matrix including intercept.}
#'   \item{Y}{Vector (length n). Outcome variable (continuous or binary).}
#'   \item{true_B}{Vector (length p). True effect sizes for each variant.}
#'   \item{true_PI}{Vector (length p). Causal probabilities (1 for causal, 0 otherwise).}
#'   \item{causal_idx}{Integer vector. Indices of causal variants.}
#'   \item{maf}{Vector (length p). Minor allele frequencies.}
#'   \item{description}{Character. Description of the simulated data.}
#'
#' @details
#' The function simulates genetic data under different scenarios:
#'
#' **Genotype Generation:**
#' - Variants generated from binomial(2, maf) distribution
#' - If ld_structure=TRUE, uses multivariate normal with AR(1) correlation
#'
#' **Outcome Generation:**
#' - Continuous: Y = X %*% beta_X + G %*% beta_G + epsilon
#' - Binary: logit(P(Y=1)) = X %*% beta_X + G %*% beta_G
#'
#' **Causal Variants:**
#' - Randomly selected n_causal variants have non-zero effects
#' - Effect sizes follow N(effect_size, effect_size/2)
#' - true_PI is indicator vector for causal variants
#'
#' @examples
#' \dontrun{
#' # Simple continuous outcome
#' data <- generate_test_data(n = 100, p = 10)
#'
#' # Binary outcome with rare variants
#' data <- generate_test_data(n = 500, p = 15, binary = TRUE, rare = TRUE)
#'
#' # Correlated variants
#' data <- generate_test_data(n = 200, p = 20, ld_structure = TRUE, ld_strength = 0.5)
#' }
#'
#' @export
generate_test_data <- function(
    n,
    p,
    k = 2,
    binary = FALSE,
    maf = NULL,
    rare = FALSE,
    ld_structure = FALSE,
    ld_strength = 0.3,
    n_causal = max(1, floor(p / 5)),
    effect_size = 0.5,
    seed = NULL,
    prevalence = 0.3
) {
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate MAF if not provided
  if (is.null(maf)) {
    if (rare) {
      maf <- runif(p, 0.001, 0.01)
    } else {
      maf <- runif(p, 0.05, 0.5)
    }
  }

  # Generate genotype matrix
  if (ld_structure) {
    # Create correlation matrix with AR(1) structure
    rho_mat <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        rho_mat[i, j] <- ld_strength^abs(i - j)
      }
    }

    # Generate correlated normal variables
    Z <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = rho_mat)

    # Convert to genotypes using MAF
    G <- matrix(0, n, p)
    for (j in 1:p) {
      # Transform to uniform then to genotype
      u <- pnorm(Z[, j])
      G[, j] <- as.numeric(u < maf[j]^2) * 2 +
                as.numeric(u >= maf[j]^2 & u < (maf[j]^2 + 2*maf[j]*(1-maf[j]))) * 1
    }
  } else {
    # Independent genotypes
    G <- matrix(0, n, p)
    for (j in 1:p) {
      G[, j] <- rbinom(n, 2, maf[j])
    }
  }

  # Generate covariates (intercept + k covariates)
  X <- cbind(1, matrix(rnorm(n * k), n, k))
  colnames(X) <- c("Intercept", paste0("Cov", 1:k))

  # Determine causal variants
  causal_idx <- sample(1:p, n_causal)
  true_B <- numeric(p)
  true_B[causal_idx] <- rnorm(n_causal, mean = effect_size, sd = effect_size / 2)
  true_PI <- numeric(p)
  true_PI[causal_idx] <- 1

  # Generate outcome
  beta_X <- rnorm(k + 1, 0, 0.3)
  linear_pred <- X %*% beta_X + G %*% true_B

  if (binary) {
    # Binary outcome
    # Adjust intercept to match prevalence
    prob <- 1 / (1 + exp(-linear_pred))
    # Calibrate intercept
    beta_X[1] <- beta_X[1] + log(prevalence / (1 - prevalence)) - mean(log(prob / (1 - prob)))
    linear_pred <- X %*% beta_X + G %*% true_B
    prob <- 1 / (1 + exp(-linear_pred))
    Y <- rbinom(n, 1, prob)
  } else {
    # Continuous outcome
    Y <- as.vector(linear_pred + rnorm(n, 0, 1))
  }

  # Create description
  desc <- sprintf(
    "n=%d, p=%d, k=%d, %s outcome, %s variants (MAF: %.3f-%.3f), %s, %d causal",
    n, p, k,
    ifelse(binary, "binary", "continuous"),
    ifelse(rare, "rare", "common"),
    min(maf), max(maf),
    ifelse(ld_structure, sprintf("LD r=%.2f", ld_strength), "independent"),
    n_causal
  )

  # Return list
  return(list(
    G = G,
    X = X,
    Y = Y,
    true_B = true_B,
    true_PI = true_PI,
    causal_idx = causal_idx,
    maf = maf,
    description = desc
  ))
}


#' Generate Validation Report
#'
#' Creates a formatted report summarizing validation test results
#' and writes it to the inst/validation/reports directory.
#'
#' @param test_results List of comparison results from compare_with_legacy().
#' @param report_name Character. Name for the report file (without extension).
#' @param report_dir Character. Directory for saving reports.
#'   Defaults to inst/validation/reports.
#' @param format Character. Report format: "txt" or "html". Default is "txt".
#'
#' @return Character. Path to the generated report file.
#'
#' @details
#' The report includes:
#' - Summary statistics: total tests, passed, failed
#' - Details for each test result
#' - Maximum and mean differences across all tests
#' - List of failed tests with details
#' - Timestamp and R session information
#'
#' The report is saved as a text file with timestamp in the filename.
#'
#' @examples
#' \dontrun{
#' results <- list(
#'   test1 = compare_with_legacy(glowr_res1, legacy_res1),
#'   test2 = compare_with_legacy(glowr_res2, legacy_res2)
#' )
#' report_path <- validation_report(results, "score_calculation_validation")
#' }
#'
#' @export
validation_report <- function(
    test_results,
    report_name = "validation_report",
    report_dir = NULL,
    format = "txt"
) {
  # Determine report directory
  if (is.null(report_dir)) {
    # Try to find package root
    pkg_root <- tryCatch({
      rprojroot::find_package_root_file()
    }, error = function(e) {
      # Fallback to relative path from tests/testthat
      "../../inst/validation/reports"
    })

    if (grepl("testthat$", pkg_root)) {
      report_dir <- file.path(pkg_root, "..", "..", "inst", "validation", "reports")
    } else {
      report_dir <- file.path(pkg_root, "inst", "validation", "reports")
    }
  }

  # Create directory if it doesn't exist
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }

  # Generate filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  report_file <- file.path(
    report_dir,
    paste0(report_name, "_", timestamp, ".", format)
  )

  # Calculate summary statistics
  n_tests <- length(test_results)
  n_passed <- sum(sapply(test_results, function(x) x$passed))
  n_failed <- n_tests - n_passed

  all_max_diffs <- sapply(test_results, function(x) x$max_diff)
  all_mean_diffs <- sapply(test_results, function(x) x$mean_diff)

  # Generate report content
  report_lines <- c(
    paste(rep("=", 70), collapse = ""),
    "GLOWr Validation Report",
    paste(rep("=", 70), collapse = ""),
    "",
    sprintf("Report: %s", report_name),
    sprintf("Generated: %s", Sys.time()),
    sprintf("R version: %s", R.version.string),
    "",
    paste(rep("-", 70), collapse = ""),
    "SUMMARY",
    paste(rep("-", 70), collapse = ""),
    sprintf("Total tests: %d", n_tests),
    sprintf("Passed: %d (%.1f%%)", n_passed, 100 * n_passed / n_tests),
    sprintf("Failed: %d (%.1f%%)", n_failed, 100 * n_failed / n_tests),
    "",
    sprintf("Overall max difference: %.2e", max(all_max_diffs, na.rm = TRUE)),
    sprintf("Overall mean difference: %.2e", mean(all_mean_diffs, na.rm = TRUE)),
    "",
    paste(rep("-", 70), collapse = ""),
    "DETAILED RESULTS",
    paste(rep("-", 70), collapse = "")
  )

  # Add details for each test
  for (i in seq_along(test_results)) {
    test_name <- names(test_results)[i]
    if (is.null(test_name)) test_name <- paste("Test", i)

    res <- test_results[[i]]

    report_lines <- c(
      report_lines,
      "",
      sprintf("[%d] %s", i, test_name),
      sprintf("  Status: %s", ifelse(res$passed, "PASSED", "FAILED")),
      sprintf("  Max diff: %.2e", res$max_diff),
      sprintf("  Mean diff: %.2e", res$mean_diff),
      sprintf("  Details: %s", res$details)
    )

    # Add mismatch details if failed
    if (!res$passed && length(res$mismatches) > 0) {
      report_lines <- c(
        report_lines,
        "  Mismatches:"
      )
      for (mm_name in names(res$mismatches)) {
        report_lines <- c(
          report_lines,
          sprintf("    - %s: %s", mm_name,
                  substr(as.character(res$mismatches[[mm_name]]), 1, 60))
        )
      }
    }
  }

  # Add failed tests summary
  if (n_failed > 0) {
    report_lines <- c(
      report_lines,
      "",
      paste(rep("-", 70), collapse = ""),
      "FAILED TESTS SUMMARY",
      paste(rep("-", 70), collapse = "")
    )

    failed_idx <- which(sapply(test_results, function(x) !x$passed))
    for (idx in failed_idx) {
      test_name <- names(test_results)[idx]
      if (is.null(test_name)) test_name <- paste("Test", idx)
      res <- test_results[[idx]]

      report_lines <- c(
        report_lines,
        sprintf("- %s: %s", test_name, res$details)
      )
    }
  }

  # Footer
  report_lines <- c(
    report_lines,
    "",
    paste(rep("=", 70), collapse = ""),
    "End of Report",
    paste(rep("=", 70), collapse = "")
  )

  # Write report
  writeLines(report_lines, report_file)

  message("Validation report written to: ", report_file)
  return(report_file)
}
