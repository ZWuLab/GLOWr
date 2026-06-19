########## Helper Functions for B Estimation ##########
#
# This file contains helper functions used by get_B and related functions

#' Calculate Cross-Validation R-squared
#'
#' @description
#' Internal function to calculate leave-one-out cross-validation R-squared
#' for model selection. Uses the analytically efficient PRESS (Predicted
#' Residual Sum of Squares) statistic to compute LOOCV predictions without
#' refitting the model n times.
#'
#' The PRESS statistic computes LOOCV predictions using the formula:
#' \deqn{\hat{y}_{(-i)} = y_i - \frac{r_i}{1 - h_{ii}}}
#' where \eqn{r_i} is the residual from the full model and \eqn{h_{ii}} is
#' the i-th diagonal element of the hat matrix.
#'
#' This approach is mathematically equivalent to explicitly refitting the model
#' n times, but is O(n) instead of O(n^2) in complexity, providing a 50-100x
#' speedup for typical training sample sizes (n=50-100).
#'
#' @param X Predictor variable(s) - can be a vector, matrix, or data frame
#' @param Y Response variable (numeric vector of length n)
#' @param model_formula Formula for the model (e.g., Y ~ X or Y ~ X + I(X^2))
#' @param model_env Environment containing variables referenced in model_formula.
#'   Default is parent.frame(), which allows access to variables in calling environment
#' @param verbose Verbosity level (0 = silent, 1+ = messages)
#'
#' @return Cross-validated R-squared value (numeric).
#'   Can be negative if predictions are worse than using the mean.
#'   Range: (-Inf, 1], where:
#'   \itemize{
#'     \item 1.0 = perfect predictions
#'     \item 0.0 = predictions as good as using the mean
#'     \item < 0 = predictions worse than using the mean (poor model fit)
#'   }
#'
#' @details
#' The function handles edge cases gracefully:
#' - Returns 0 if the model fails to fit
#' - Returns 0 if R^2 is NA or infinite
#' - Returns negative values if the model performs poorly (this is valid and expected)
#'
#' For high-leverage points (h_ii close to 1), the denominator is protected
#' against division by zero using a minimum threshold.
#'
#' IMPORTANT: When model_formula references transformed variables (e.g., logY, logX),
#' those variables must either be:
#' 1. Included in the data frame X (if X is a data frame), OR
#' 2. Available in model_env (the calling environment)
#'
#' @references
#' Allen, D. M. (1974). "The Relationship Between Variable Selection and Data
#' Augmentation and a Method for Prediction." Technometrics, 16(1), 125-127.
#'
#' @keywords internal
#' @noRd
cv_r_squared <- function(X, Y, model_formula, model_env = parent.frame(), verbose = 0) {

  # Get sample size
  n <- length(Y)

  # Create data frame for modeling
  # Handle different input types: vector, matrix, or data frame
  if (is.matrix(X) || is.data.frame(X)) {
    model_data <- cbind(Y = Y, X)
  } else {
    model_data <- data.frame(Y = Y, X = X)
  }

  # Fit full model once (not n times!)
  # Ensure the formula can find variables in model_env
  # by setting the formula's environment attribute
  environment(model_formula) <- model_env

  fit <- tryCatch({
    lm(model_formula, data = model_data)
  }, error = function(e) {
    # If model fails to fit, return NULL
    # This can happen with rank-deficient designs or other issues
    if (verbose >= 1) {
      message("  Warning: Model fitting failed in cv_r_squared: ", e$message)
    }
    return(NULL)
  })

  # If model fitting failed, return 0
  if (is.null(fit)) {
    return(0)
  }

  # Get the response variable from the fitted model
  # This ensures we use the same transformed Y that was used in fitting
  # (e.g., if formula is logY ~ X, this gets logY values)
  Y_response <- model.response(model.frame(fit))

  # Calculate LOOCV predictions analytically using PRESS statistic
  # Formula: y_hat_(-i) = y_i - r_i / (1 - h_ii)
  # where r_i = residual, h_ii = i-th diagonal of hat matrix

  # Get residuals from full model
  residuals <- residuals(fit)

  # Get diagonal of hat matrix H = X(X'X)^(-1)X'
  # hatvalues() is an efficient built-in R function for this
  hat_values <- hatvalues(fit)

  # Protect against division by zero for high-leverage points
  # (points with h_ii very close to 1)
  safe_hat <- pmax(hat_values, 1e-10)

  # Compute LOOCV predictions
  # This is mathematically equivalent to refitting n models,
  # but only requires fitting one model
  predictions <- Y_response - residuals / (1 - safe_hat)

  # Calculate R-squared from LOOCV predictions
  ss_res <- sum((Y_response - predictions)^2)   # Residual sum of squares
  ss_tot <- sum((Y_response - mean(Y_response))^2)       # Total sum of squares

  r_squared <- 1 - (ss_res / ss_tot) # Negative r_squared is possible by the formula - it indicates predictions are worse than using the mean. 

  # Verbose output for debugging
  if (verbose >= 1) {
    message("  CV R^2 calculation: ss_res = ", round(ss_res, 4),
            ", ss_tot = ", round(ss_tot, 4),
            ", raw CV R^2 = ", round(r_squared, 4))
  }

  # Handle edge cases
  # This can occur with perfect fits, rank-deficient models, etc.
  if (is.na(r_squared) || is.infinite(r_squared)) {
    if (verbose >= 1) {
      message("  Warning: R^2 calculation resulted in ", r_squared, ", returning 0")
    }
    return(0)
  }

  return(r_squared)
}
