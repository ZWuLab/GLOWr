# This file contains input validation helper functions for GLOW tests.

#' Validate genotype matrix for GLOW tests
#'
#' Checks that a genotype matrix meets basic requirements for variant-set
#' association testing. This includes checking dimensions, missing values,
#' valid genotype encoding, and sufficient variation.
#'
#' @description
#' This function performs comprehensive validation of genotype data to ensure
#' it meets the requirements for GLOW variant-set tests. It checks:
#' - Matrix structure (non-empty, numeric)
#' - Missing values (optionally allowed or disallowed)
#' - Genotype encoding (0/1/2 for additive, or continuous dosages)
#' - Variant variation (at least some variants must be polymorphic)
#'
#' @details
#' **Validation checks performed**:
#'
#' 1. **Matrix structure**:
#'    - Must be a matrix or coercible to matrix
#'    - Must have at least 1 sample and 1 variant
#'    - Must contain numeric values
#'
#' 2. **Missing values**:
#'    - If allow_na=FALSE: No missing values allowed
#'    - If allow_na=TRUE: Missing values are permitted but counted
#'
#' 3. **Genotype encoding**:
#'    - If discrete=TRUE: All non-NA values must be in {0, 1, 2}
#'    - If discrete=FALSE: Any numeric values in \eqn{[0, 2]} are allowed (dosages)
#'
#' 4. **Variation check**:
#'    - At least one variant must have standard deviation > min_sd
#'    - This ensures at least some genetic variation exists
#'    - Monomorphic variants (all 0s or all same value) are flagged
#'
#' **Computational Complexity**: O(n * p) where n=samples, p=variants
#'
#' **Returns**: List with validation status and diagnostic information
#'
#' @param G Genotype matrix (samples x variants). Rows are individuals,
#'   columns are genetic variants. Expected encoding: 0/1/2 for discrete
#'   genotypes, or \eqn{[0,2]} for dosages.
#' @param allow_na Logical; if TRUE, allow missing values (NA). Default: FALSE
#' @param discrete Logical; if TRUE, require discrete 0/1/2 encoding.
#'   If FALSE, allow continuous dosages in \eqn{[0,2]}. Default: TRUE
#' @param min_sd Minimum standard deviation required for at least one variant.
#'   Ensures sufficient genetic variation. Default: 1e-6
#'
#' @return A list with elements:
#' \describe{
#'   \item{valid}{Logical; TRUE if all checks pass}
#'   \item{message}{Character string with validation result or error message}
#'   \item{n_samples}{Number of samples (rows)}
#'   \item{n_variants}{Number of variants (columns)}
#'   \item{n_missing}{Number of missing values (if allow_na=TRUE)}
#'   \item{n_monomorphic}{Number of monomorphic variants (SD < min_sd)}
#'   \item{warnings}{Character vector of non-fatal warnings}
#' }
#'
#' @examples
#' # Valid genotype matrix (discrete, no missing)
#' G <- matrix(sample(0:2, 100, replace=TRUE), nrow=10, ncol=10)
#' result <- validate_genotype_matrix(G)
#' result$valid  # TRUE
#'
#' # Genotype matrix with missing values
#' G_na <- G
#' G_na[1:3, 1] <- NA
#' result_na <- validate_genotype_matrix(G_na, allow_na=TRUE)
#' result_na$n_missing  # 3
#'
#' # Continuous dosages
#' G_dosage <- matrix(runif(100, 0, 2), nrow=10, ncol=10)
#' result_dosage <- validate_genotype_matrix(G_dosage, discrete=FALSE)
#'
#' # Invalid: all monomorphic
#' G_mono <- matrix(0, nrow=10, ncol=10)
#' result_mono <- validate_genotype_matrix(G_mono)
#' result_mono$valid  # FALSE
#'
#' @export
validate_genotype_matrix <- function(G, allow_na = FALSE, discrete = TRUE, min_sd = 1e-6) {

  # Initialize result list
  result <- list(
    valid = FALSE,
    message = "",
    n_samples = 0,
    n_variants = 0,
    n_missing = 0,
    n_monomorphic = 0,
    warnings = character(0)
  )

  # Check 1: Is it a matrix?
  if (!is.matrix(G)) {
    # Try to coerce to matrix
    G <- tryCatch(
      as.matrix(G),
      error = function(e) {
        result$message <<- "Input G cannot be coerced to a matrix"
        return(NULL)
      }
    )
    if (is.null(G)) {
      return(result)
    }
  }

  # Check 2: Dimensions
  result$n_samples <- nrow(G)
  result$n_variants <- ncol(G)

  if (result$n_samples == 0 || result$n_variants == 0) {
    result$message <- sprintf(
      "Genotype matrix must have at least 1 sample and 1 variant (found %d x %d)",
      result$n_samples, result$n_variants
    )
    return(result)
  }

  # Check 3: Numeric values
  if (!is.numeric(G)) {
    result$message <- "Genotype matrix must contain numeric values"
    return(result)
  }

  # Check 4: Missing values
  result$n_missing <- sum(is.na(G))
  if (result$n_missing > 0 && !allow_na) {
    result$message <- sprintf(
      "Genotype matrix contains %d missing values. Set allow_na=TRUE to permit NAs",
      result$n_missing
    )
    return(result)
  }

  # Check 5: Genotype encoding
  non_na_values <- G[!is.na(G)]

  if (discrete) {
    # For discrete genotypes: must be in {0, 1, 2}
    if (!all(non_na_values %in% c(0, 1, 2))) {
      invalid_vals <- unique(non_na_values[!(non_na_values %in% c(0, 1, 2))])
      result$message <- sprintf(
        "Discrete genotypes must be 0, 1, or 2. Found invalid values: %s",
        paste(head(invalid_vals, 5), collapse=", ")
      )
      return(result)
    }
  } else {
    # For continuous dosages: must be in [0, 2]
    if (any(non_na_values < 0) || any(non_na_values > 2)) {
      result$message <- "Genotype dosages must be in the range [0, 2]"
      return(result)
    }
  }

  # Check 6: Variation (at least one variant should be polymorphic)
  # Calculate standard deviation for each variant (column)
  variant_sds <- apply(G, 2, function(x) {
    # Remove NAs for SD calculation
    x_noNA <- x[!is.na(x)]
    if (length(x_noNA) < 2) return(0)
    return(sd(x_noNA))
  })

  result$n_monomorphic <- sum(variant_sds < min_sd)

  if (all(variant_sds < min_sd)) {
    result$message <- sprintf(
      "All %d variants are monomorphic (SD < %g). No genetic variation detected",
      result$n_variants, min_sd
    )
    return(result)
  }

  # Warn if many monomorphic variants
  if (result$n_monomorphic > 0) {
    pct_mono <- 100 * result$n_monomorphic / result$n_variants
    if (pct_mono > 50) {
      result$warnings <- c(
        result$warnings,
        sprintf("%.1f%% of variants are monomorphic (%d/%d)",
                pct_mono, result$n_monomorphic, result$n_variants)
      )
    }
  }

  # Warn if many missing values
  if (result$n_missing > 0) {
    pct_missing <- 100 * result$n_missing / (result$n_samples * result$n_variants)
    if (pct_missing > 10) {
      result$warnings <- c(
        result$warnings,
        sprintf("%.1f%% of genotypes are missing (%d total)",
                pct_missing, result$n_missing)
      )
    }
  }

  # All checks passed
  result$valid <- TRUE
  result$message <- sprintf(
    "Valid genotype matrix: %d samples x %d variants",
    result$n_samples, result$n_variants
  )

  return(result)
}


#' Validate correlation matrix
#'
#' Checks that a matrix is a valid correlation matrix suitable for use in
#' GLOW tests. This includes checking symmetry, diagonal values, positive
#' semi-definiteness, and appropriate dimensions.
#'
#' @description
#' Validates that a matrix can serve as a correlation matrix for variant-set
#' tests. Correlation matrices are used to account for linkage disequilibrium
#' among genetic variants.
#'
#' @details
#' **Validation checks performed**:
#'
#' 1. **Basic structure**:
#'    - Must be a square matrix
#'    - Must be numeric
#'
#' 2. **Correlation properties**:
#'    - Must be symmetric (M = t(M) within tolerance)
#'    - Diagonal must be all 1s (within tolerance)
#'    - Off-diagonal values must be in \eqn{[-1, 1]}
#'
#' 3. **Positive semi-definiteness**:
#'    - All eigenvalues must be >= -tol
#'    - Required for valid correlation matrices
#'    - Small negative eigenvalues (< tol) are acceptable due to numerical error
#'
#' 4. **Dimension matching** (if n_variants provided):
#'    - Matrix dimension must match expected number of variants
#'
#' **Computational Complexity**: O(n^3) for eigenvalue check where n=dim(M)
#'
#' **Returns**: List with validation status and diagnostic information
#'
#' @param M Matrix to validate as correlation matrix
#' @param n_variants Expected dimension (optional). If provided, checks that
#'   nrow(M) == ncol(M) == n_variants
#' @param tol Numerical tolerance for symmetry and positive definiteness checks.
#'   Default: 1e-8
#'
#' @return A list with elements:
#' \describe{
#'   \item{valid}{Logical; TRUE if all checks pass}
#'   \item{message}{Character string with validation result or error message}
#'   \item{dimension}{Dimension of the matrix}
#'   \item{is_symmetric}{Logical; TRUE if symmetric within tolerance}
#'   \item{min_eigenvalue}{Smallest eigenvalue (for PSD check)}
#'   \item{warnings}{Character vector of non-fatal warnings}
#' }
#'
#' @examples
#' # Valid correlation matrix
#' M <- diag(10)
#' result <- validate_correlation_matrix(M)
#' result$valid  # TRUE
#'
#' # Valid with correlation structure
#' M_cor <- matrix(0.3, 5, 5) + diag(0.7, 5)
#' result_cor <- validate_correlation_matrix(M_cor)
#'
#' # Invalid: not symmetric
#' M_bad <- matrix(runif(25), 5, 5)
#' result_bad <- validate_correlation_matrix(M_bad)
#' result_bad$valid  # FALSE
#'
#' # Check dimension matching
#' result_dim <- validate_correlation_matrix(M_cor, n_variants=5)
#' result_dim$valid  # TRUE
#'
#' @export
validate_correlation_matrix <- function(M, n_variants = NULL, tol = 1e-8) {

  # Initialize result list
  result <- list(
    valid = FALSE,
    message = "",
    dimension = 0,
    is_symmetric = FALSE,
    min_eigenvalue = NA,
    warnings = character(0)
  )

  # Check 1: Is it a matrix?
  if (!is.matrix(M)) {
    result$message <- "Input M must be a matrix"
    return(result)
  }

  # Check 2: Is it numeric?
  if (!is.numeric(M)) {
    result$message <- "Correlation matrix M must contain numeric values"
    return(result)
  }

  # Check 3: Is it square?
  if (nrow(M) != ncol(M)) {
    result$message <- sprintf(
      "Correlation matrix must be square (found %d x %d)",
      nrow(M), ncol(M)
    )
    return(result)
  }

  result$dimension <- nrow(M)

  # Check 4: Dimension matching (if specified)
  if (!is.null(n_variants)) {
    if (result$dimension != n_variants) {
      result$message <- sprintf(
        "Correlation matrix dimension (%d) does not match n_variants (%d)",
        result$dimension, n_variants
      )
      return(result)
    }
  }

  # Check 5: Is it symmetric?
  max_asymmetry <- max(abs(M - t(M)))
  result$is_symmetric <- (max_asymmetry < tol)

  if (!result$is_symmetric) {
    result$message <- sprintf(
      "Correlation matrix must be symmetric (max asymmetry = %.2e, tol = %.2e)",
      max_asymmetry, tol
    )
    return(result)
  }

  # Check 6: Are diagonal elements all 1?
  diag_vals <- diag(M)
  if (!all(abs(diag_vals - 1) < tol)) {
    bad_diag <- which(abs(diag_vals - 1) >= tol)
    result$message <- sprintf(
      "Correlation matrix diagonal must be all 1s. Found non-1 values at positions: %s",
      paste(head(bad_diag, 10), collapse=", ")
    )
    return(result)
  }

  # Check 7: Are off-diagonal elements in [-1, 1]?
  off_diag <- M[row(M) != col(M)]
  if (any(off_diag < -1 - tol) || any(off_diag > 1 + tol)) {
    result$message <- "Correlation matrix off-diagonal values must be in [-1, 1]"
    return(result)
  }

  # Check 8: Is it positive semi-definite?
  # Calculate eigenvalues to check PSD property
  eigenvalues <- tryCatch(
    eigen(M, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) {
      result$message <<- sprintf("Failed to compute eigenvalues: %s", e$message)
      return(NULL)
    }
  )

  if (is.null(eigenvalues)) {
    return(result)
  }

  result$min_eigenvalue <- min(eigenvalues)

  # Allow small negative eigenvalues due to numerical error
  if (result$min_eigenvalue < -tol) {
    result$message <- sprintf(
      "Correlation matrix is not positive semi-definite (min eigenvalue = %.2e < -%.2e)",
      result$min_eigenvalue, tol
    )
    return(result)
  }

  # Warn if matrix is nearly singular
  max_eigenvalue <- max(eigenvalues)
  condition_number <- max_eigenvalue / max(abs(result$min_eigenvalue), 1e-16)

  if (condition_number > 1e10) {
    result$warnings <- c(
      result$warnings,
      sprintf("Correlation matrix is nearly singular (condition number = %.2e)",
              condition_number)
    )
  }

  # Warn if any very high correlations (may indicate multicollinearity)
  max_off_diag <- max(abs(off_diag))
  if (max_off_diag > 0.99) {
    result$warnings <- c(
      result$warnings,
      sprintf("Very high correlation detected (max = %.4f), may indicate multicollinearity",
              max_off_diag)
    )
  }

  # All checks passed
  result$valid <- TRUE
  result$message <- sprintf(
    "Valid correlation matrix: %d x %d",
    result$dimension, result$dimension
  )

  return(result)
}
