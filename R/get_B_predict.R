########## Predict B Values from a Trained Model ##########
#
# This file contains predict_B() for applying a trained B estimation model
# to new MAF values without re-training. Handles both beta method models
# (direct B prediction) and p-value method models (h-squared to B conversion).
#
# EXPORTED FUNCTIONS:
#   predict_B() - Predict B values from a trained model for new MAF values
#
# INTERNAL HELPERS:
#   .predict_B_from_lm() - Core prediction + back-transformation from lm
#   .convert_h2_to_B()   - Convert h (sqrt of h-squared) predictions to B values

#################### EXPORTED MAIN FUNCTIONS ####################

#' Predict B Values from a Trained Model
#'
#' @description
#' Applies a trained B estimation model to new MAF values without re-training.
#' This enables the "train once, predict many times" workflow needed by the
#' GLOW testing integration pipeline.
#'
#' For models trained via the beta method, the function predicts BETA^2 from
#' MAF and back-transforms to B = sqrt(BETA^2).
#'
#' For models trained via the p-value method, the function predicts h-squared
#' from MAF and converts to B using trait-specific formulas:
#' \itemize{
#'   \item Continuous traits: \eqn{B = \sqrt{h^2} \times SE}
#'   \item Binary traits: \eqn{B = \sqrt{h^2} / \sqrt{p_0 (1 - p_0)}}
#'     where \eqn{p_0} is the case proportion
#' }
#'
#' @param model A \code{glow_B_model} object (from \code{\link{train_B_model}}),
#'   or a raw \code{lm} object with a \code{model_id} attribute indicating
#'   which of the 8 candidate models it represents. For \code{glow_B_model}
#'   objects, the function extracts the appropriate lm model based on the
#'   \code{method} parameter.
#' @param target_MAF Numeric vector of MAF values to predict B for. Values
#'   should be in \code{[0, 0.5]}. Values > 0.5 are automatically folded to
#'   \code{1 - MAF} with a warning. MAF = 0 returns B = 0.
#' @param target_trait Character string or NULL (default NULL). Required when
#'   using a p-value method model: \code{"binary"} or \code{"continuous"}.
#'   Ignored for beta method models.
#' @param target_SE Numeric scalar or vector of standard errors for the target
#'   data. Required for p-value method with continuous traits. If scalar,
#'   recycled to match length of \code{target_MAF}. If vector, must have
#'   same length as \code{target_MAF}.
#' @param target_case_prop Numeric scalar or vector of case proportions for
#'   binary traits. Required for p-value method with binary traits. Values
#'   must be in (0, 1). If scalar, recycled to match length of
#'   \code{target_MAF}. If vector, must have same length as \code{target_MAF}.
#' @param method Character (default \code{"auto"}). For \code{glow_B_model}
#'   objects: \code{"auto"} uses the primary method from the model,
#'   \code{"beta_method"} forces the beta method model,
#'   \code{"pvalue_method"} forces the p-value method model. Ignored for
#'   raw \code{lm} objects.
#'
#' @return Numeric vector of predicted B values (same length as
#'   \code{target_MAF}). Non-negative (floored at 0).
#'
#' @details
#' \strong{Beta Method Models:}
#'
#' The function handles all 8 candidate model transformations:
#' \itemize{
#'   \item Models 1-4: Y response (B = sqrt(abs(predicted)))
#'   \item Models 5-8: log(Y) response (B = sqrt(exp(predicted)))
#' }
#'
#' Each model uses one of 4 predictor transformations:
#' \itemize{
#'   \item X (MAF directly)
#'   \item fX = MAF * (1 - MAF)
#'   \item logX = log(MAF)
#'   \item logfX = log(MAF * (1 - MAF))
#' }
#'
#' \strong{P-value Method Models:}
#'
#' These models predict h-squared (heritability) rather than BETA^2.
#' The h-squared predictions must be converted to B values using trait-specific
#' formulas. This conversion was previously performed inside
#' \code{.estimate_B_pvalue()} but has been moved here to support the
#' train/predict separation.
#'
#' If \code{target_trait} is not provided for a p-value method model, the
#' function returns \code{sqrt(h^2)} with a warning (an approximation
#' that ignores trait-specific scaling).
#'
#' \strong{Method Selection for glow_B_model:}
#' \itemize{
#'   \item \code{method = "auto"}: Uses the primary method. If the model was
#'     trained with a single method, that method is used. If trained with
#'     "both", uses the method selected at training time by criterion
#'     comparison (\code{model$comparison$method_selected}).
#'   \item \code{method = "beta_method"}: Forces the beta method model.
#'   \item \code{method = "pvalue_method"}: Forces the p-value method model.
#' }
#'
#' @examples
#' \dontrun{
#' # Train a B model
#' b_model <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = known_snps$MAF,
#'   training_BETA = known_snps$BETA,
#'   method = "beta"
#' )
#'
#' # Predict B for new MAF values (beta method -- no conversion needed)
#' new_maf <- c(0.001, 0.01, 0.05, 0.1, 0.2)
#' B_predicted <- predict_B(b_model, new_maf)
#'
#' # P-value method with binary trait conversion
#' b_model_pval <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = known_snps$MAF,
#'   training_P = known_snps$P,
#'   training_N = known_snps$N,
#'   method = "pvalue"
#' )
#' B_binary <- predict_B(b_model_pval, new_maf,
#'                        target_trait = "binary",
#'                        target_case_prop = 0.3)
#'
#' # P-value method with continuous trait conversion
#' B_continuous <- predict_B(b_model_pval, new_maf,
#'                            target_trait = "continuous",
#'                            target_SE = 0.5)
#' }
#'
#' @seealso
#' \code{\link{train_B_model}} for training B estimation models
#' \code{\link{get_B}} for combined train + predict in one call
#' \code{\link{save_B_model}} for saving trained models
#' \code{\link{load_B_model}} for loading saved models
#'
#' @export
predict_B <- function(model, target_MAF,
                      target_trait = NULL,
                      target_SE = NULL,
                      target_case_prop = NULL,
                      method = "auto") {

  # ========== Input Validation ==========

  # Validate target_MAF
  if (!is.numeric(target_MAF)) {
    stop("target_MAF must be a numeric vector")
  }
  if (any(is.na(target_MAF))) {
    stop("target_MAF must not contain NA values")
  }

  # Fold MAF values > 0.5 to minor allele frequency convention
  # After folding, validates that all values are in [0, 0.5]
  # (allow_zero = TRUE because MAF = 0 is valid for prediction -> returns B = 0)
  target_MAF <- .fold_maf(target_MAF, allow_zero = TRUE,
                          context = "target_MAF")

  # Handle empty input
  if (length(target_MAF) == 0) {
    return(numeric(0))
  }

  # Validate method parameter
  valid_methods <- c("auto", "beta_method", "pvalue_method")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # ========== Extract lm Object and Determine Source Method ==========

  # source_method tracks whether the lm came from beta or pvalue training.

  # This determines whether h2-to-B conversion is needed after prediction.
  source_method <- NULL

  if (inherits(model, "glow_B_model")) {
    # --- glow_B_model input ---

    # Determine which method to use
    if (method == "auto") {
      resolved_method <- model$method_used
      # If method_used is "both", use the selected method from comparison
      if (resolved_method == "both") {
        resolved_method <- model$comparison$method_selected
      }
    } else {
      resolved_method <- method
    }

    lm_obj <- model$models[[resolved_method]]
    if (is.null(lm_obj)) {
      available <- names(Filter(Negate(is.null), model$models))
      stop("Method '", resolved_method, "' not available in this model. ",
           "Available methods: ", paste(available, collapse = ", "))
    }

    source_method <- resolved_method

  } else if (inherits(model, "lm")) {
    # --- Raw lm input ---
    lm_obj <- model
    if (is.null(attr(lm_obj, "model_id"))) {
      stop("Raw lm object must have a 'model_id' attribute (1-8). ",
           "Use a glow_B_model object from train_B_model() instead.")
    }

    # For raw lm objects, assume beta method (direct B prediction)
    # unless the caller explicitly specifies otherwise via method
    source_method <- "beta_method"

  } else {
    stop("model must be a 'glow_B_model' object or an 'lm' object, ",
         "got class: ", paste(class(model), collapse = ", "))
  }

  # ========== Identify Zero-MAF Entries ==========

  # MAF = 0 is a special case: no allele means no effect
  zero_maf <- target_MAF == 0
  nonzero_maf <- target_MAF[!zero_maf]

  # If all MAFs are zero, return all zeros
  if (length(nonzero_maf) == 0) {
    return(rep(0, length(target_MAF)))
  }

  # ========== Predict from lm Object ==========

  # Core prediction and back-transformation (shared for both methods)
  predicted_raw <- .predict_B_from_lm(lm_obj, nonzero_maf)

  # ========== Post-Processing Based on Source Method ==========

  if (source_method == "pvalue_method") {
    # P-value method: predicted_raw contains h = sqrt(h-squared) (already
    # back-transformed by .predict_B_from_lm). We need to convert h -> B.
    B_nonzero <- .convert_h2_to_B(
      h = predicted_raw,
      target_trait = target_trait,
      target_SE = target_SE,
      target_case_prop = target_case_prop,
      n_target = length(nonzero_maf)
    )
  } else {
    # Beta method: predicted_raw is already B values
    B_nonzero <- predicted_raw
  }

  # ========== Enforce Non-Negative ==========

  B_nonzero[B_nonzero < 0] <- 0

  # ========== Reassemble with Zero-MAF Entries ==========

  B <- numeric(length(target_MAF))
  B[!zero_maf] <- B_nonzero
  B[zero_maf] <- 0

  # Return unnamed numeric vector
  as.numeric(B)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Core Prediction and Back-Transformation from lm Object
#'
#' @description
#' Internal helper that performs the prediction from an lm object and
#' applies appropriate back-transformation based on whether the model
#' used a log(Y) response. For beta method models, returns B values
#' (sqrt of predicted Y). For p-value method models, returns h-squared
#' values (the raw predicted quantity before trait-specific conversion).
#'
#' @param lm_obj A fitted lm object with model_id attribute
#' @param nonzero_maf Numeric vector of non-zero MAF values
#'
#' @return Numeric vector of back-transformed predictions. For beta
#'   method: B = sqrt(|Y|) or sqrt(exp(logY)). For p-value method:
#'   h2 = |Y| or exp(logY) (before sqrt, since caller handles h2->B).
#'
#' @details
#' The back-transformation logic depends on the response variable in
#' the model formula:
#' \itemize{
#'   \item If "logY" is in the model terms: predicted value is log(Y),
#'     so result = sqrt(|exp(predicted)|).
#'   \item Otherwise: predicted value is Y directly,
#'     so result = sqrt(|predicted|).
#' }
#'
#' This function applies the same back-transformation regardless of
#' whether the model came from beta or pvalue training. For beta method
#' models, Y = BETA^2, so the result is B = sqrt(BETA^2). For p-value
#' method models, Y = h^2, so the result is h = sqrt(h^2). The caller
#' (predict_B) is responsible for the additional h-to-B conversion when
#' the source method is p-value.
#'
#' @keywords internal
#' @noRd
.predict_B_from_lm <- function(lm_obj, nonzero_maf) {
  # Construct transformed predictors
  # Provide all 4 predictor columns; predict() picks the right one
  newdata <- data.frame(
    X    = nonzero_maf,
    fX   = nonzero_maf * (1 - nonzero_maf),
    logX = log(nonzero_maf),
    logfX = log(nonzero_maf * (1 - nonzero_maf))
  )

  # Predict
  predicted <- predict(lm_obj, newdata = newdata)

  # Back-transform based on response variable
  model_terms <- all.vars(formula(lm_obj))

  if ("logY" %in% model_terms) {
    # Model predicted log(Y); Y = B^2 (or h^2), so B = sqrt(exp(predicted))
    B <- sqrt(abs(exp(predicted)))
  } else {
    # Model predicted Y directly; B = sqrt(abs(Y))
    B <- sqrt(abs(predicted))
  }

  as.numeric(B)
}


#' Convert h-squared Predictions to B Values
#'
#' @description
#' Internal helper that converts h-squared (heritability) predictions from
#' the p-value method to trait-specific B (effect size) values. This logic
#' was previously in \code{.estimate_B_pvalue()} and has been moved here
#' as part of the train/predict separation refactoring.
#'
#' @param h Numeric vector of sqrt(h-squared) values from .predict_B_from_lm().
#'   The function receives h = sqrt(h2), not h2 itself, because
#'   .predict_B_from_lm() applies sqrt() during back-transformation.
#' @param target_trait Character string or NULL: "binary" or "continuous"
#' @param target_SE Numeric scalar or vector. Required for continuous traits
#' @param target_case_prop Numeric scalar or vector. Required for binary traits
#' @param n_target Integer, the number of target MAF values (for recycling)
#'
#' @return Numeric vector of B values
#'
#' @details
#' The input \code{h2} is actually h = sqrt(h-squared), because
#' .predict_B_from_lm() already applies sqrt() to the predicted Y values.
#' For the p-value method, Y = h-squared, so sqrt(Y) = h.
#'
#' Conversion formulas (matching the old .estimate_B_pvalue() logic):
#' \itemize{
#'   \item Continuous trait: B = h * SE (where h = sqrt(h2))
#'   \item Binary trait: B = h / sqrt(p0 * (1 - p0))
#'     where p0 is the case proportion
#'   \item NULL trait: B = h (returns sqrt(h2) with a warning)
#' }
#'
#' @keywords internal
#' @noRd
.convert_h2_to_B <- function(h, target_trait, target_SE, target_case_prop,
                              n_target) {
  # h is sqrt(h-squared), produced by .predict_B_from_lm() which applies
  # sqrt(abs(predicted)) to the model's predicted h^2 values

  if (is.null(target_trait)) {
    # No trait specified: return h = sqrt(h2) as approximation
    warning("P-value method model requires target_trait for proper h2-to-B ",
            "conversion. Returning sqrt(h2) as approximation (no trait-specific ",
            "scaling applied).")
    return(h)
  }

  if (target_trait == "continuous") {
    # Continuous trait: B = h * SE
    if (is.null(target_SE)) {
      stop("target_SE is required for continuous target trait with p-value method")
    }
    # Allow target_SE to be scalar or vector
    if (length(target_SE) == 1) {
      target_SE <- rep(target_SE, n_target)
    } else if (length(target_SE) != n_target) {
      stop("target_SE must be either a single value or have the same length ",
           "as target_MAF")
    }
    B <- h * target_SE

  } else if (target_trait == "binary") {
    # Binary trait: B = h / sqrt(p0 * (1 - p0))
    # where p0 is case proportion in case-control study
    if (is.null(target_case_prop)) {
      stop("target_case_prop is required for binary target trait with ",
           "p-value method")
    }
    # Allow target_case_prop to be scalar or vector
    if (length(target_case_prop) == 1) {
      target_case_prop <- rep(target_case_prop, n_target)
    } else if (length(target_case_prop) != n_target) {
      stop("target_case_prop must be either a single value or have the same ",
           "length as target_MAF")
    }
    if (any(target_case_prop <= 0 | target_case_prop >= 1)) {
      stop("target_case_prop values must be in the interval (0, 1)")
    }
    B <- h / sqrt(target_case_prop * (1 - target_case_prop))

  } else {
    stop("target_trait must be 'binary', 'continuous', or NULL, got: '",
         target_trait, "'")
  }

  as.numeric(B)
}
