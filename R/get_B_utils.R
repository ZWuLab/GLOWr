########## Shared Utilities for B Estimation ##########
#
# This file contains internal helper functions shared across the B estimation
# framework. These utilities handle data preparation and model fitting to ensure
# consistency and avoid code duplication.
#
# INTERNAL HELPERS (all @noRd):
#   .fold_maf()                - Fold MAF values > 0.5 to minor allele frequency
#   .compute_chi2_from_pvalue() - Convert p-values to chi-squared statistics
#   .prepare_training_data_for_models() - Prepare training data for model fitting
#   .fit_all_candidate_models() - Fit all 8 candidate regression models
#   .validate_trait_compatibility() - Validate trait compatibility between datasets
#   .extract_model_criterion() - Extract model quality metric for comparison
#   .train_B_beta_method() - Train B model using direct beta method (model only)
#   .train_B_pvalue_method() - Train B model using p-value/Z-score method (model only)
#   .compare_B_methods() - Compare prediction results from both methods

#################### INTERNAL HELPER FUNCTIONS ####################

#' Fold MAF Values Greater Than 0.5 to Minor Allele Frequency
#'
#' @description
#' Internal helper that converts allele frequency values greater than 0.5
#' to their complement (1 - MAF), ensuring the minor allele frequency
#' convention where MAF is in (0, 0.5] or [0, 0.5]. Issues a warning
#' when folding is performed.
#'
#' @param maf Numeric vector of allele frequencies
#' @param allow_zero Logical. If TRUE, allow MAF = 0 (for prediction where
#'   MAF = 0 returns B = 0). Default FALSE (for training, where MAF = 0
#'   is not informative).
#' @param context Character string describing the parameter being validated,
#'   used in warning/error messages (e.g., "training_MAF", "target_MAF").
#'
#' @return Numeric vector with all values in (0, 0.5] (if allow_zero = FALSE)
#'   or [0, 0.5] (if allow_zero = TRUE).
#'
#' @details
#' The function performs two steps:
#' \enumerate{
#'   \item Fold: Any non-NA value > 0.5 is replaced with 1 - value, and a
#'     warning is issued indicating how many values were folded.
#'   \item Validate: After folding, checks that all non-NA values fall within
#'     the expected range. Stops with an error if any values are out of range.
#' }
#'
#' This is useful when input data may contain allele frequencies coded as
#' major allele frequency (e.g., 0.7 instead of 0.3). The folding converts
#' these to the minor allele frequency convention.
#'
#' @keywords internal
#' @noRd
.fold_maf <- function(maf, allow_zero = FALSE, context = "MAF") {
  # Fold values > 0.5 to their complement
  needs_fold <- !is.na(maf) & maf > 0.5
  if (any(needs_fold)) {
    n_folded <- sum(needs_fold)
    warning(n_folded, " ", context, " value(s) > 0.5 detected. ",
            "Converting to 1 - MAF (minor allele frequency convention).",
            call. = FALSE)
    maf[needs_fold] <- 1 - maf[needs_fold]
  }

  # Validate range after folding
  lower <- if (allow_zero) 0 else .Machine$double.eps
  if (any(!is.na(maf) & (maf < lower | maf > 0.5))) {
    stop(context, " values must be in ",
         if (allow_zero) "[0, 0.5]" else "(0, 0.5]",
         " after folding. Found values outside this range.",
         call. = FALSE)
  }

  maf
}


#' Compute Chi-Squared Statistic from P-value
#'
#' @description
#' Internal helper function that converts p-values to chi-squared statistics
#' with robust handling of extremely small p-values that cause numerical
#' underflow. Supports both regular p-values and -log10 transformed p-values
#' for numerical stability.
#'
#' @param p Numeric vector of p-values in (0, 1). Either \code{p} or
#'   \code{p_mlog10} must be provided, but not both
#' @param p_mlog10 Numeric vector of -log10 transformed p-values. When provided,
#'   this is used preferentially over \code{p} for numerical stability
#'
#' @details
#' The function uses different computational approaches based on the magnitude
#' of the p-value:
#'
#' \itemize{
#'   \item For -log10(p) <= 300: Uses the exact inverse CDF via
#'     \code{qchisq(10^(-p_mlog10), df = 1, lower.tail = FALSE)}
#'   \item For -log10(p) > 300: Uses the asymptotic approximation
#'     \eqn{\chi^2 \approx 2 \times \ln(10) \times (-\log_{10}(p))}
#' }
#'
#' The threshold of 300 is chosen conservatively based on the numerical limits
#' of \code{qchisq}, which works reliably up to approximately -log10(p) = 323.
#'
#' When only \code{p} is provided, the function protects against underflow by
#' ensuring p >= .Machine$double.xmin before calling qchisq.
#'
#' @return Numeric vector of chi-squared statistics (df = 1) corresponding to
#'   the input p-values
#'
#' @section Mathematical Details:
#' For large -log10(p), the chi-squared quantile can be approximated using the
#' asymptotic behavior of the chi-squared tail probability. The approximation:
#'
#' \deqn{\chi^2 \approx 2 \times \ln(10) \times p_{mlog10} \approx 4.60517 \times p_{mlog10}}
#'
#' becomes highly accurate when -log10(p) > 300 (corresponding to p < 10^-300).
#'
#' \strong{Derivation via Mills' Ratio:}
#'
#' For a chi-squared statistic with 1 degree of freedom, we have \eqn{\chi^2 = Z^2}
#' where Z is a standard normal variable. The two-sided p-value is
#' \eqn{p = 2\Phi(-|Z|)} where \eqn{\Phi} is the standard normal CDF.
#'
#' For large Z (extreme p-values), Mills' ratio provides the asymptotic expansion
#' for the tail probability:
#' \deqn{1 - \Phi(z) \sim \frac{\phi(z)}{z} = \frac{1}{z\sqrt{2\pi}} e^{-z^2/2}}
#'
#' Taking the logarithm:
#' \deqn{\ln(1 - \Phi(z)) \sim -\frac{z^2}{2} - \ln(z) - \frac{1}{2}\ln(2\pi)}
#'
#' For very large z, the dominant term is \eqn{-z^2/2}, while the \eqn{\ln(z)} and
#' constant terms become negligible. Therefore:
#' \deqn{\ln(p/2) \approx \ln(1 - \Phi(z)) \approx -\frac{z^2}{2}}
#'
#' This gives us:
#' \deqn{\chi^2 = z^2 \approx -2\ln(p/2) = -2\ln(p) + 2\ln(2)}
#'
#' For extremely small p, the \eqn{2\ln(2) \approx 1.39} term is negligible
#' compared to \eqn{-2\ln(p)}, yielding:
#' \deqn{\chi^2 \approx -2\ln(p)}
#'
#' Converting from natural log to base-10 log:
#' \deqn{-2\ln(p) = -2\ln(10^{-p_{mlog10}}) = 2 \times \ln(10) \times p_{mlog10}}
#'
#' where \eqn{\ln(10) \approx 2.302585}, giving the coefficient 4.60517.
#'
#' This derivation shows that \eqn{4.60517 \times (-\log_{10}(p))} is indeed the
#' dominant term in calculating the quantile for extreme p-values, with all
#' higher-order terms becoming negligible as p approaches zero.
#'
#' @references
#' Mills' ratio and asymptotic expansions: Abramowitz & Stegun (1972) Handbook of
#' Mathematical Functions, Section 26.2.12
#'
#' Chi-square tail probability: Johnson, Kotz & Balakrishnan (1994) Continuous
#' Univariate Distributions Vol 1, Sec 18.3
#'
#' @examples
#' \dontrun{
#' # Regular p-values
#' p <- c(0.05, 0.01, 1e-10)
#' chi2 <- .compute_chi2_from_pvalue(p = p)
#'
#' # Extremely small p-values using -log10 transformation
#' p_mlog10 <- c(10, 50, 100, 350)  # p-values from 1e-10 to 1e-350
#' chi2_log <- .compute_chi2_from_pvalue(p_mlog10 = p_mlog10)
#' }
#'
#' @keywords internal
#' @noRd
.compute_chi2_from_pvalue <- function(p = NULL, p_mlog10 = NULL) {
  # Validate input: need exactly one of p or p_mlog10
  if (is.null(p) && is.null(p_mlog10)) {
    stop("Either p or p_mlog10 must be provided")
  }
  if (!is.null(p) && !is.null(p_mlog10)) {
    stop("Only one of p or p_mlog10 should be provided, not both")
  }

  # If p_mlog10 is provided, use it preferentially for better numerical stability
  if (!is.null(p_mlog10)) {
    # Validate p_mlog10 values (should be positive)
    if (any(p_mlog10 < 0)) {
      stop("p_mlog10 values must be non-negative (representing -log10 of p-values)")
    }

    # Initialize result vector
    chi2 <- numeric(length(p_mlog10))

    # For moderate p-values (p_mlog10 <= 324), use exact qchisq
    moderate_idx <- p_mlog10 <= 324
    if (any(moderate_idx)) {
      p_values <- 10^(-p_mlog10[moderate_idx])
      chi2[moderate_idx] <- qchisq(p_values, df = 1, lower.tail = FALSE)
    }

    # For extreme p-values (p_mlog10 > 324), use asymptotic approximation
    # chi2 ≈ 2 * ln(10) * p_mlog10
    extreme_idx <- p_mlog10 > 324
    if (any(extreme_idx)) {
      # ln(10) ≈ 2.302585
      chi2[extreme_idx] <- 2 * log(10) * p_mlog10[extreme_idx]
    }

    return(chi2)

  } else {
    # Using p directly
    # Validate p values
    if (any(p < 0 | p > 1)) {
      stop("p values must be in the interval [0, 1]")
    }

    # Protect against underflow: ensure p >= machine epsilon
    p_protected <- pmax(p, .Machine$double.xmin)

    # Use qchisq with protected p-values
    chi2 <- qchisq(p_protected, df = 1, lower.tail = FALSE)

    return(chi2)
  }
}


#' Prepare Training Data for Model Fitting (Internal)
#'
#' @description
#' Internal helper function that prepares training data for B estimation model
#' fitting. Handles both direct beta method and p-value/Z-score method with
#' consistent data transformations and filtering.
#'
#' This function centralizes the data preparation logic to ensure consistency
#' across select_best_model() and compare_B_models().
#'
#' @param training_MAF Numeric vector of minor allele frequencies from training data
#' @param training_BETA Numeric vector of effect sizes (for beta method, optional)
#' @param training_P Numeric vector of p-values (for pvalue method, optional)
#' @param training_P_mlog10 Numeric vector of -log10 transformed p-values (optional)
#' @param training_N Numeric vector of sample sizes (for pvalue method, optional)
#' @param method Character string, either "beta" or "pvalue"
#' @param verbose Integer verbosity level (0=silent, 1=warnings, 2=info, 3=debug)
#'
#' @return List with elements:
#'   \itemize{
#'     \item X: Processed MAF values (numeric vector)
#'     \item Y: Processed response variable (numeric vector)
#'     \item n_removed: Number of observations removed due to invalid values
#'     \item n_total: Total number of observations before filtering
#'   }
#'
#' @details
#' For beta method:
#' \itemize{
#'   \item Y = BETA^2
#'   \item Removes observations with non-finite Y or invalid MAF
#' }
#'
#' For pvalue method:
#' \itemize{
#'   \item Converts P-values to chi-squared statistics using .compute_chi2_from_pvalue()
#'   \item Handles extremely small p-values (up to 10^-350) with proper numerical precision
#'   \item Uses asymptotic approximation for p-values < 10^-324 to avoid underflow
#'   \item Calculates trait-independent h^2 = Z^2 / (2 * N * MAF * (1-MAF))
#'   \item Removes observations with non-finite values
#' }
#'
#' @section Computational Complexity:
#' O(n) where n = length(training_MAF)
#'
#' @keywords internal
#' @noRd
.prepare_training_data_for_models <- function(
    training_MAF,
    training_BETA = NULL,
    training_P = NULL,
    training_P_mlog10 = NULL,
    training_N = NULL,
    method = "beta",
    verbose = 0) {

  # Validate method parameter
  if (!method %in% c("beta", "pvalue")) {
    stop("method must be either 'beta' or 'pvalue'")
  }

  if (method == "beta") {
    # Beta method: Y = BETA^2
    if (is.null(training_BETA)) {
      stop("training_BETA is required for beta method")
    }
    if (length(training_BETA) != length(training_MAF)) {
      stop("training_BETA and training_MAF must have the same length")
    }

    Y <- training_BETA^2

  } else {
    # P-value method: Y = h^2
    # Validate inputs
    if (is.null(training_P) && is.null(training_P_mlog10)) {
      stop("Either training_P or training_P_mlog10 is required for pvalue method")
    }
    if (!is.null(training_P) && length(training_P) != length(training_MAF)) {
      stop("training_P and training_MAF must have the same length")
    }
    if (!is.null(training_P_mlog10) && length(training_P_mlog10) != length(training_MAF)) {
      stop("training_P_mlog10 and training_MAF must have the same length")
    }
    if (is.null(training_N)) {
      stop("training_N is required for pvalue method")
    }
    if (length(training_N) != length(training_MAF)) {
      stop("training_N and training_MAF must have the same length")
    }

    # Convert p-values to chi-squared statistics using robust helper
    # This handles extreme p-values (up to 10^-350) with proper numerical precision
    Z2 <- .compute_chi2_from_pvalue(p = training_P, p_mlog10 = training_P_mlog10)

    # Calculate trait-independent h^2
    Y <- Z2 / (2 * training_N * training_MAF * (1 - training_MAF))
  }

  # Remove observations with non-finite values
  # This handles NA, NaN, Inf, -Inf in both Y and MAF
  valid_idx <- is.finite(Y) & is.finite(training_MAF) &
               training_MAF > 0 & training_MAF < 1

  n_removed <- sum(!valid_idx)
  n_total <- length(training_MAF)

  if (verbose >= 1 && n_removed > 0) {
    message(sprintf("Removed %d observations with invalid values (%.1f%% of data)",
                   n_removed, 100 * n_removed / n_total))
  }

  # Check if we have enough data left
  if (sum(valid_idx) < 3) {
    warning("Only ", sum(valid_idx), " valid observations remaining after filtering. ",
            "Results may be unreliable.")
  }

  return(list(
    X = training_MAF[valid_idx],
    Y = Y[valid_idx],
    n_removed = n_removed,
    n_total = n_total
  ))
}


#' Fit All Candidate Models (Internal)
#'
#' @description
#' Internal helper function that fits all candidate regression models for
#' B estimation. This is the ONLY place where models are fitted in the entire
#' B estimation framework, ensuring consistency across model selection and
#' model comparison functions.
#'
#' The function fits 8 standard models plus any custom models provided.
#' Standard models represent all combinations of transformations:
#' \itemize{
#'   \item Y transformations: Y, log(Y)
#'   \item X transformations: X, f(X)=X(1-X), log(X), log(f(X))
#' }
#'
#' @param X Numeric vector of MAF values (already filtered for valid range)
#' @param Y Numeric vector of response values (BETA^2 or h^2, already filtered)
#' @param custom_models Optional list of custom formula objects to fit in addition
#'   to the 8 standard models. Custom models can use X, Y, and their transformations
#'   (fX, logX, logfX, logY) which are available in the fitting environment
#' @param compute_cv Logical. If TRUE, computes cross-validated R^2 for each model.
#'   Default is FALSE for performance (CV is expensive). Only set to TRUE when
#'   CV_R2 criterion is being used
#' @param verbose Integer verbosity level (0=silent, 1=warnings, 2=info, 3=debug)
#'
#' @return Named list with one element per successfully fitted model. Each element
#'   is a list containing:
#'   \itemize{
#'     \item model: The fitted lm object
#'     \item formula: The formula used (as formula object)
#'     \item formula_text: Character representation of formula (for printing)
#'     \item R2: R-squared value from summary
#'     \item adj_R2: Adjusted R-squared value from summary
#'     \item CV_R2: Cross-validated R^2 (NA if compute_cv=FALSE or if CV failed)
#'     \item sigma: Residual standard error
#'     \item coefficients: Named vector of model coefficients
#'   }
#'
#'   Models that fail to fit are omitted from the returned list. An empty list
#'   indicates all models failed.
#'
#' @details
#' The 8 candidate models are:
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
#' Each model is fitted using lm() with error handling. Failed models (e.g.,
#' due to numerical issues with log transformations) are skipped with optional
#' warning messages.
#'
#' The transformations are computed directly in this function's environment
#' to ensure they are available for both model fitting and cross-validation.
#'
#' @section Computational Complexity:
#' \itemize{
#'   \item Without CV: O(k * n) where k=8 models, n=length(X)
#'   \item With CV: O(k * n^2) due to leave-one-out cross-validation
#' }
#'
#' For typical use: n=50, k=8 -> 400 model fits without CV, 20,000 with CV
#'
#' @section Performance Note:
#' Cross-validation is computationally expensive (n-fold increase). Only compute
#' CV_R2 when it will actually be used for model selection (criterion="CV_R2").
#'
#' @keywords internal
#' @noRd
.fit_all_candidate_models <- function(X, Y, custom_models = NULL, compute_cv = FALSE, verbose = 0) {
  # Validate inputs
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length")
  }
  if (length(X) < 3) {
    stop("Need at least 3 observations to fit models")
  }

  # Validate custom models if provided
  if (!is.null(custom_models)) {
    if (!is.list(custom_models)) {
      stop("custom_models must be a list of formula objects")
    }
    if (!all(sapply(custom_models, inherits, "formula"))) {
      stop("All elements of custom_models must be formula objects")
    }
  }

  # Create transformations in this environment
  # These will be available for both lm() and cv_r_squared()
  fX <- X * (1 - X)
  logX <- log(X)
  logfX <- log(fX)
  logY <- log(Y)

  # Define all 8 model formulas
  # Note: The formula names are used for identification
  model_formulas <- list(
    "Y ~ X" = Y ~ X,
    "Y ~ f(X)" = Y ~ fX,
    "Y ~ log(X)" = Y ~ logX,
    "Y ~ log(f(X))" = Y ~ logfX,
    "log(Y) ~ X" = logY ~ X,
    "log(Y) ~ f(X)" = logY ~ fX,
    "log(Y) ~ log(X)" = logY ~ logX,
    "log(Y) ~ log(f(X))" = logY ~ logfX
  )

  # Fit all models
  results <- list()

  for (i in seq_along(model_formulas)) {
    model_name <- names(model_formulas)[i]
    fmla <- model_formulas[[i]]

    tryCatch({
      # Fit model using data from this environment
      fit <- lm(fmla, data = environment())
      summ <- summary(fit)

      # Calculate CV R^2 if requested
      cv_r2 <- if (compute_cv) {
        tryCatch({
          # Pass environment() so cv_r_squared can access transformed variables
          cv_r_squared(X, Y, fmla, model_env = environment(), verbose = 0)
        }, error = function(e) {
          if (verbose >= 1) {
            warning(sprintf("CV calculation failed for model '%s': %s",
                           model_name, e$message))
          }
          NA_real_
        })
      } else {
        NA_real_
      }

      # Attach model_id attribute to the lm object (1-8 for standard models)
      attr(fit, "model_id") <- i

      # Store comprehensive model information
      results[[model_name]] <- list(
        model = fit,
        formula = fmla,
        formula_text = paste(deparse(fmla), collapse = " "),
        R2 = summ$r.squared,
        adj_R2 = summ$adj.r.squared,
        CV_R2 = cv_r2,
        sigma = summ$sigma,
        coefficients = coef(fit)
      )

    }, error = function(e) {
      if (verbose >= 1) {
        warning(sprintf("Model '%s' failed to fit: %s", model_name, e$message))
      }
      # Model will be omitted from results
    })
  }

  # Fit custom models if provided
  if (!is.null(custom_models)) {
    if (verbose >= 2) {
      message(sprintf("Fitting %d custom model(s)...", length(custom_models)))
    }

    for (j in seq_along(custom_models)) {
      # Generate unique name for custom model
      custom_name <- if (!is.null(names(custom_models)[j]) && names(custom_models)[j] != "") {
        names(custom_models)[j]
      } else {
        sprintf("Custom_%d", j)
      }

      fmla <- custom_models[[j]]

      tryCatch({
        # Fit custom model using data from this environment
        fit <- lm(fmla, data = environment())
        summ <- summary(fit)

        # Compute CV R^2 if requested
        cv_r2 <- NA_real_
        if (compute_cv) {
          tryCatch({
            # Pass this environment so cv_r_squared can access variables
            cv_r2 <- cv_r_squared(X, Y, fmla, model_env = environment(), verbose = 0)
          }, error = function(e) {
            if (verbose >= 1) {
              warning(sprintf("CV calculation failed for custom model '%s': %s",
                            custom_name, e$message))
            }
          })
        }

        # Attach model_id attribute (8 + j for custom models)
        attr(fit, "model_id") <- 8L + j

        # Store custom model information
        results[[custom_name]] <- list(
          model = fit,
          formula = fmla,
          formula_text = paste(deparse(fmla), collapse = " "),
          R2 = summ$r.squared,
          adj_R2 = summ$adj.r.squared,
          CV_R2 = cv_r2,
          sigma = summ$sigma,
          coefficients = coef(fit)
        )

      }, error = function(e) {
        if (verbose >= 1) {
          warning(sprintf("Custom model '%s' (%s) failed to fit: %s",
                        custom_name, deparse(fmla), e$message))
        }
        # Model will be omitted from results
      })
    }
  }

  # Check if any models succeeded
  if (length(results) == 0) {
    warning("All candidate models (including custom) failed to fit. Check input data validity.")
  }

  return(results)
}


#' Validate Trait Compatibility Between Training and Target Data
#'
#' @description
#' Internal helper function that validates trait compatibility between training
#' and target datasets. Returns validation result with appropriate severity level
#' (error, warning, or ok) and informative message.
#'
#' @param training_trait Character string or NULL: "binary", "continuous", "mixed", or NULL
#' @param target_trait Character string or NULL: "binary", "continuous", or NULL
#' @param method Character string: "beta" or "pvalue"
#'
#' @return List with elements:
#'   \itemize{
#'     \item valid: TRUE if validation passes (warning or ok), FALSE if error
#'     \item message: Character string with warning or error message (empty if ok)
#'     \item severity: "error", "warning", or "ok"
#'   }
#'
#' @details
#' Validation logic:
#' \itemize{
#'   \item P-value method: target_trait is REQUIRED (error if NULL)
#'   \item Beta method with both traits specified and different: ERROR
#'   \item Beta method with incomplete trait info: WARNING
#'   \item Beta method with matching traits: OK
#' }
#'
#' @keywords internal
#' @noRd
.validate_trait_compatibility <- function(training_trait, target_trait, method) {
  # For p-value method, target_trait is required
  if (method == "pvalue") {
    if (is.null(target_trait)) {
      return(list(
        valid = FALSE,
        message = "P-value/Z-score method requires target_trait to be specified ('binary' or 'continuous')",
        severity = "error"
      ))
    }
    # If target_trait is specified for p-value method, it's ok
    return(list(
      valid = TRUE,
      message = "",
      severity = "ok"
    ))
  }

  # For beta method, validate trait compatibility
  # Normalize NULL and "mixed" training traits for comparison
  training_is_unspecified <- is.null(training_trait) || training_trait == "mixed"
  target_is_unspecified <- is.null(target_trait)

  # Case 1: Both traits are NULL or mixed
  if (training_is_unspecified && target_is_unspecified) {
    return(list(
      valid = TRUE,
      message = paste0(
        "Beta method assumes training and target data have the same trait type, ",
        "but trait types are unspecified. The beta method is only valid when trait types match."
      ),
      severity = "warning"
    ))
  }

  # Case 2: Only training is specified (target is NULL)
  if (!training_is_unspecified && target_is_unspecified) {
    return(list(
      valid = TRUE,
      message = paste0(
        "Beta method assumes target data has trait type '", training_trait,
        "' (same as training). If target differs, results will be invalid."
      ),
      severity = "warning"
    ))
  }

  # Case 3: Only target is specified (training is NULL/mixed)
  if (training_is_unspecified && !target_is_unspecified) {
    training_desc <- if (is.null(training_trait)) "NULL" else training_trait
    return(list(
      valid = TRUE,
      message = paste0(
        "Cannot verify trait matching: training_trait is ", training_desc,
        " while target_trait is '", target_trait,
        "'. Beta method requires matching trait types."
      ),
      severity = "warning"
    ))
  }

  # Case 4: Both are specified
  if (!training_is_unspecified && !target_is_unspecified) {
    if (training_trait != target_trait) {
      # Different traits specified - this is an ERROR
      return(list(
        valid = FALSE,
        message = paste0(
          "Beta method requires matching trait types. Training is '", training_trait,
          "' but target is '", target_trait, "'. Use p-value method instead."
        ),
        severity = "error"
      ))
    } else {
      # Same traits - OK
      return(list(
        valid = TRUE,
        message = "",
        severity = "ok"
      ))
    }
  }

  # Should not reach here, but return error if we do
  return(list(
    valid = FALSE,
    message = "Unexpected trait configuration",
    severity = "error"
  ))
}


#' Extract Model Quality Metric for Comparison
#'
#' @description
#' Internal helper function that extracts the model quality metric used for
#' selection (R2, adj_R2, or CV_R2) from a fitted lm object. This is used to
#' compare models when both beta and pvalue methods are run.
#'
#' @param model A fitted lm object returned by select_best_model()
#' @param X Numeric vector of MAF values used to train the model (for CV calculation)
#' @param Y Numeric vector of response values used to train the model (for CV calculation)
#' @param criterion Character string: "R2", "adj_R2", or "CV_R2"
#' @param verbose Integer verbosity level
#'
#' @return Numeric value of the specified criterion for the model
#'
#' @keywords internal
#' @noRd
.extract_model_criterion <- function(model, X, Y, criterion, verbose = 1) {
  # Extract the requested criterion from the model
  value <- switch(criterion,
    R2 = {
      # Higher R2 is better
      summary(model)$r.squared
    },
    adj_R2 = {
      # Higher adjusted R2 is better
      summary(model)$adj.r.squared
    },
    CV_R2 = {
      # Higher CV R2 is better
      # Need to recalculate CV for this specific model
      fmla <- formula(model)
      tryCatch({
        cv_r_squared(X, Y, fmla, verbose = 0)
      }, error = function(e) {
        if (verbose >= 1) {
          warning("CV calculation failed for model comparison, using R2 instead")
        }
        summary(model)$r.squared  # Fallback to R2
      })
    }
  )

  return(value)
}


#' Train B Model Using Direct Beta Method
#'
#' @description
#' Internal helper function that trains a regression model for B estimation
#' using the direct beta method. This method fits BETA^2 ~ f(MAF) models
#' when training and target traits are the same type.
#'
#' This function performs model training only -- it does NOT produce B
#' predictions. Use \code{predict_B()} to generate predictions from the
#' returned model.
#'
#' @param training_MAF Numeric vector of MAF from training data
#' @param training_BETA Numeric vector of effect sizes from training data
#' @param selection_criterion Character string for model selection criterion
#' @param custom_models Optional list of custom model formulas
#' @param show_model_selection Logical. If TRUE, prints model selection table
#' @param verbose Integer verbosity level
#'
#' @return List with elements:
#'   \itemize{
#'     \item model: The best fitted lm object (with model_id attribute)
#'     \item all_models_info: Complete model info from select_best_model
#'       (list with models, summary_table, criterion_used, training_data)
#'   }
#'
#' @keywords internal
#' @noRd
.train_B_beta_method <- function(training_MAF, training_BETA,
                                 selection_criterion = "R2",
                                 custom_models = NULL,
                                 show_model_selection = TRUE,
                                 verbose = 1) {
  # Prepare training data using centralized utility
  # This handles validation and data cleaning consistently
  data_prepared <- .prepare_training_data_for_models(
    training_MAF = training_MAF,
    training_BETA = training_BETA,
    method = "beta",
    verbose = verbose
  )

  X <- data_prepared$X
  Y <- data_prepared$Y

  # Select best model with user-specified criterion
  # Always request full model info for the glow_B_model object
  model_result <- select_best_model(X = X, Y = Y,
                                     criterion = selection_criterion,
                                     custom_models = custom_models,
                                     return_all = TRUE,
                                     verbose = show_model_selection)

  # Extract best model and all_models_info from the return_all result
  best_model <- model_result$best_model
  all_models_info <- model_result$all_models_info

  # Return model info only -- no predictions
  return(list(
    model = best_model,
    all_models_info = all_models_info
  ))
}


#' Train B Model Using P-value/Z-score Method
#'
#' @description
#' Internal helper function that trains a regression model for B estimation
#' using the p-value/Z-score method. This method converts p-values to
#' chi-squared statistics, then to trait-independent h^2, and fits
#' h^2 ~ f(MAF) models.
#'
#' This function performs model training only -- it does NOT produce B
#' predictions or perform the h-squared to B conversion. That conversion
#' (which requires target_trait, target_SE/target_case_prop) is handled
#' by \code{predict_B()} in Phase 2.
#'
#' @param training_MAF Numeric vector of MAF from training data
#' @param training_P Numeric vector of p-values from training data (either this or training_P_mlog10)
#' @param training_P_mlog10 Numeric vector of -log10 transformed p-values from training data
#' @param training_N Numeric vector of sample sizes from training data
#' @param selection_criterion Character string for model selection criterion
#' @param custom_models Optional list of custom model formulas
#' @param show_model_selection Logical. If TRUE, prints model selection table
#' @param verbose Integer verbosity level
#'
#' @return List with elements:
#'   \itemize{
#'     \item model: The best fitted lm object (with model_id attribute)
#'     \item all_models_info: Complete model info from select_best_model
#'       (list with models, summary_table, criterion_used, training_data)
#'   }
#'
#' @keywords internal
#' @noRd
.train_B_pvalue_method <- function(training_MAF, training_P = NULL,
                                    training_P_mlog10 = NULL,
                                    training_N,
                                    selection_criterion = "R2",
                                    custom_models = NULL,
                                    show_model_selection = TRUE,
                                    verbose = 1) {
  # Prepare training data using centralized utility
  # This handles validation, p-value conversion, and data cleaning consistently
  data_prepared <- .prepare_training_data_for_models(
    training_MAF = training_MAF,
    training_P = training_P,
    training_P_mlog10 = training_P_mlog10,
    training_N = training_N,
    method = "pvalue",
    verbose = verbose
  )

  X <- data_prepared$X
  Y <- data_prepared$Y  # This is h^2

  # Select best model for h^2 ~ MAF with user-specified criterion
  # Always request full model info for the glow_B_model object
  model_result <- select_best_model(X = X, Y = Y,
                                     criterion = selection_criterion,
                                     custom_models = custom_models,
                                     return_all = TRUE,
                                     verbose = show_model_selection)

  # Extract best model and all_models_info from the return_all result
  best_model <- model_result$best_model
  all_models_info <- model_result$all_models_info

  # Return model info only -- no predictions or h^2-to-B conversion
  return(list(
    model = best_model,
    all_models_info = all_models_info
  ))
}


#' Compare Prediction Results from Both B Estimation Methods
#'
#' @description
#' Internal helper function that calculates comparison metrics between
#' beta method and p-value method B prediction results.
#'
#' Note: This function compares B *predictions* (not model fits). It requires
#' actual B vectors from both methods, which are only available after
#' prediction. Used by \code{get_B()} after calling \code{predict_B()}.
#'
#' @param B_beta_method Numeric vector of B estimates from the beta method
#' @param B_pvalue_method Numeric vector of B estimates from the p-value method
#' @param target_MAF Numeric vector of target MAF values
#'
#' @return List with comparison metrics:
#'   \itemize{
#'     \item correlation: Pearson correlation between methods
#'     \item rmse: Root mean squared error between methods
#'     \item mean_percent_diff: Mean percentage difference
#'     \item max_percent_diff: Maximum percentage difference
#'     \item summary_stats: Data frame with summary statistics
#'   }
#'
#' @keywords internal
#' @noRd
.compare_B_methods <- function(B_beta_method, B_pvalue_method, target_MAF) {
  # Calculate comparison metrics
  correlation <- cor(B_beta_method, B_pvalue_method, use = "complete.obs")
  rmse <- sqrt(mean((B_beta_method - B_pvalue_method)^2, na.rm = TRUE))

  # Calculate percentage difference
  # Use mean of two methods as denominator to avoid division by very small numbers
  percent_diff <- 100 * abs(B_beta_method - B_pvalue_method) /
    ((B_beta_method + B_pvalue_method) / 2)

  # Summary statistics
  summary_stats <- data.frame(
    Method = c("Beta Method", "P-value Method"),
    Mean_B = c(mean(B_beta_method, na.rm = TRUE),
               mean(B_pvalue_method, na.rm = TRUE)),
    SD_B = c(sd(B_beta_method, na.rm = TRUE),
             sd(B_pvalue_method, na.rm = TRUE)),
    Min_B = c(min(B_beta_method, na.rm = TRUE),
              min(B_pvalue_method, na.rm = TRUE)),
    Max_B = c(max(B_beta_method, na.rm = TRUE),
              max(B_pvalue_method, na.rm = TRUE))
  )

  return(list(
    correlation = correlation,
    rmse = rmse,
    mean_percent_diff = mean(percent_diff, na.rm = TRUE),
    max_percent_diff = max(percent_diff, na.rm = TRUE),
    summary_stats = summary_stats
  ))
}
