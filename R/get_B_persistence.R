########## Save and Load Trained B Estimation Models ##########
#
# This file contains save_B_model() and load_B_model() for persisting
# trained B estimation models with metadata validation.
#
# EXPORTED FUNCTIONS:
#   save_B_model()  - Save a glow_B_model (or glow_B_estimate) to RDS
#   load_B_model()  - Load and validate a saved glow_B_model from RDS

#################### EXPORTED MAIN FUNCTIONS ####################

#' Save a Trained B Estimation Model
#'
#' @description
#' Saves a trained \code{glow_B_model} object to an RDS file with a metadata
#' envelope for validation on load. The metadata includes package version,
#' save date, training size, and model selection information.
#'
#' Also accepts a \code{glow_B_estimate} object (from \code{\link{get_B}} with
#' \code{return_full = TRUE}), from which the nested \code{$model} is extracted
#' and saved.
#'
#' @param model A \code{glow_B_model} object (from \code{\link{train_B_model}})
#'   or a \code{glow_B_estimate} object (from \code{\link{get_B}} with
#'   \code{return_full = TRUE}). If a \code{glow_B_estimate} is provided,
#'   the nested \code{$model} (a \code{glow_B_model}) is extracted and saved.
#' @param path Character. File path for the saved model (.rds extension).
#' @param compress Logical (default \code{TRUE}). Use compression when saving.
#'
#' @return Invisible character: the path where the model was saved.
#'
#' @examples
#' \dontrun{
#' # Save a trained model directly
#' b_model <- train_B_model(
#'   training_trait = "binary",
#'   training_MAF = known_snps$MAF,
#'   training_BETA = known_snps$BETA,
#'   method = "beta"
#' )
#' save_B_model(b_model, "B_model_als.rds")
#'
#' # Or save from a full estimation result
#' B_result <- get_B(
#'   training_trait = "binary",
#'   training_MAF = known_snps$MAF,
#'   training_BETA = known_snps$BETA,
#'   target_MAF = seq(0.001, 0.3, length.out = 100),
#'   return_full = TRUE
#' )
#' save_B_model(B_result, "B_model_als.rds")
#' }
#'
#' @seealso
#' \code{\link{load_B_model}} for loading saved models
#' \code{\link{predict_B}} for applying loaded models to new data
#' \code{\link{train_B_model}} for training B estimation models
#' \code{\link{get_B}} for combined train + predict
#'
#' @export
save_B_model <- function(model, path, compress = TRUE) {

  # Accept glow_B_model directly or glow_B_estimate (extract $model)
  if (inherits(model, "glow_B_estimate")) {
    b_model <- model$model
  } else if (inherits(model, "glow_B_model")) {
    b_model <- model
  } else {
    stop("model must be a 'glow_B_model' or 'glow_B_estimate' object")
  }

  # Determine best method for metadata
  best_method <- b_model$method_used
  if (best_method == "both" && !is.null(b_model$comparison$method_selected)) {
    best_method <- b_model$comparison$method_selected
  }

  # Extract model_id from the best lm object
  best_lm <- b_model$models[[best_method]]
  best_model_id <- if (!is.null(best_lm)) attr(best_lm, "model_id") else NA_integer_

  # Extract R2 from the best model
  best_R2 <- if (!is.null(best_lm)) summary(best_lm)$r.squared else NA_real_

  # Build metadata
  metadata <- list(
    glowr_version = as.character(utils::packageVersion("GLOWr")),
    save_date     = Sys.Date(),
    training_n    = b_model$training_summary$n_used,
    best_method   = best_method,
    best_model_id = best_model_id,
    R2            = best_R2
  )

  # Wrap the glow_B_model (not the full estimate) in envelope
  envelope <- list(
    model          = b_model,
    metadata       = metadata,
    format_version = 1L
  )

  # Save
  saveRDS(envelope, file = path, compress = compress)
  message("Saved B model to ", path)

  invisible(path)
}


#' Load a Trained B Estimation Model
#'
#' @description
#' Loads a saved \code{glow_B_model} object and optionally validates its
#' structure, including checking the metadata envelope and testing prediction
#' capability.
#'
#' @param path Character. Path to a saved B model (.rds file).
#' @param validate Logical (default \code{TRUE}). Check object structure on load.
#'
#' @return A \code{glow_B_model} object. Use \code{\link{predict_B}} to
#'   apply it to new MAF values.
#'
#' @examples
#' \dontrun{
#' B_model <- load_B_model("B_model_als.rds")
#' B_values <- predict_B(B_model, target_MAF = c(0.01, 0.05, 0.1))
#' }
#'
#' @seealso
#' \code{\link{save_B_model}} for saving models
#' \code{\link{predict_B}} for applying loaded models to new data
#' \code{\link{train_B_model}} for training B estimation models
#'
#' @export
load_B_model <- function(path, validate = TRUE) {

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  envelope <- readRDS(path)

  if (validate) {
    # Check envelope structure
    required_fields <- c("model", "metadata", "format_version")
    missing <- setdiff(required_fields, names(envelope))
    if (length(missing) > 0) {
      stop("Invalid B model file: missing fields: ",
           paste(missing, collapse = ", "),
           ". Was this file created by save_B_model()?")
    }

    # Check model class: accept glow_B_model (new) or glow_B_estimate (legacy)
    loaded_model <- envelope$model
    if (inherits(loaded_model, "glow_B_estimate")) {
      # Legacy format: extract the nested glow_B_model
      loaded_model <- loaded_model$model
    }

    if (!inherits(loaded_model, "glow_B_model")) {
      stop("Loaded object is not a glow_B_model (class: ",
           paste(class(loaded_model), collapse = ", "), ")")
    }

    # Check that at least one lm model is present
    has_model <- !is.null(loaded_model$models$beta_method) ||
                 !is.null(loaded_model$models$pvalue_method)
    if (!has_model) {
      stop("Loaded model contains no fitted lm objects")
    }

    # Test prediction capability
    lm_obj <- if (!is.null(loaded_model$models$beta_method)) {
      loaded_model$models$beta_method
    } else {
      loaded_model$models$pvalue_method
    }
    tryCatch({
      test_newdata <- data.frame(
        X = 0.01, fX = 0.01 * 0.99,
        logX = log(0.01), logfX = log(0.01 * 0.99)
      )
      predict(lm_obj, newdata = test_newdata)
    }, error = function(e) {
      stop("Loaded model failed prediction test: ", e$message)
    })

    # Version check
    current <- as.character(utils::packageVersion("GLOWr"))
    if (!is.null(envelope$metadata$glowr_version) &&
        envelope$metadata$glowr_version != current) {
      warning("Model was saved with GLOWr ", envelope$metadata$glowr_version,
              ", current version is ", current)
    }

    # Return the extracted/validated glow_B_model
    return(loaded_model)
  }

  # No validation: extract appropriately
  loaded_model <- envelope$model
  if (inherits(loaded_model, "glow_B_estimate")) {
    loaded_model <- loaded_model$model
  }
  loaded_model
}
