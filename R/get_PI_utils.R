#' PI Model Utility Functions
#'
#' Functions for summarizing, evaluating, and visualizing trained PI
#' (variant-importance score) ensemble models. Supports both LASSO (glmnet)
#' and GLM models.
#'
#' File Log (reverse chronological order):
#' - 2026-06-13: Fixed by Claude Code (Opus 4.8, r-developer agent), prompted by
#'   ZWu for the GLOWr release-gate Rd hygiene. Rewrote the unit-interval @return
#'   range of the internal .compute_auc helper using \eqn{...} so the markup is
#'   text, not a markdown link (the helper is @noRd, so this is preventive).
#'   Docs-only; no logic change.
#' - 2026-02-04: Modified by Claude Code - plot_PI_coefficient_summary: removed
#'   "x." prefix from feature labels, replaced `las` param with `label_angle`
#'   (default 45) for flexible rotation angle control (useful for publications).
#' - 2026-02-04: Modified by Claude Code - GLM models: added p-value extraction
#'   and significance metrics (pct_significant, mean_neg_log10_pval). Removed
#'   coefficient boxplot (coefficients not comparable across features due to
#'   different scales). Simplified plot function to bar plot only. Removed
#'   duplicate .default_PI_features() (use version from get_PI_train.R).
#' - 2026-02-03: Created by Claude Code - Initial implementation with
#'   load_PI_models, summarize_PI_coefficients, predict_PI_ensemble,
#'   evaluate_PI_models, and plotting functions (base R, no dependencies)
#'
#' @note Uses \code{impute_na_median()} and \code{.default_PI_features()} from
#'   get_PI_train.R (accessible within package namespace).
#' @note Uses \code{glmnet} package for LASSO model predictions.
#' @note Uses \code{data.table::fread()} for efficient file loading.
#'
#' @name get_PI_utils
#' @docType package
NULL


#################### EXPORTED MAIN FUNCTIONS ####################


#' Load PI Models from Directory
#'
#' @description
#' Loads all PI model RDS files from a directory into a list. Automatically
#' detects whether models are LASSO (glmnet) or GLM type.
#'
#' @param model_dir Character. Path to directory containing model_*.rds files.
#' @param pattern Character. Regex pattern for model file names.
#'   Default: "^model_.*\\.rds$"
#'
#' @return List with components:
#'   \item{models}{List of fitted model objects}
#'   \item{model_type}{"LASSO" or "GLM" (auto-detected from first model)}
#'   \item{n_models}{Number of models loaded}
#'   \item{model_dir}{Source directory path}
#'
#' @examples
#' \dontrun{
#' # Load models from directory
#' result <- load_PI_models("path/to/piModels/")
#' cat("Loaded", result$n_models, result$model_type, "models\n")
#' }
#'
#' @export
load_PI_models <- function(model_dir, pattern = "^model_.*\\.rds$") {

  # Validate directory exists
  if (!dir.exists(model_dir)) {
    stop("Model directory does not exist: ", model_dir)
  }

  # List matching files
  model_files <- list.files(model_dir, pattern = pattern, full.names = TRUE)

  if (length(model_files) == 0) {
    stop("No model files matching pattern '", pattern, "' found in: ", model_dir)
  }

  # Sort by model number (extract numeric from filename)
  # e.g., model_1.rds, model_2.rds, ..., model_100.rds
  model_nums <- as.numeric(gsub(".*model_([0-9]+)\\.rds$", "\\1",
                                 basename(model_files)))
  model_files <- model_files[order(model_nums)]

  # Load all models
  models <- lapply(model_files, readRDS)

  # Detect model type from first model
  model_type <- .detect_PI_model_type(models[[1]])

  return(list(
    models = models,
    model_type = model_type,
    n_models = length(models),
    model_dir = normalizePath(model_dir)
  ))
}


#' Summarize PI Model Coefficients
#'
#' @description
#' Extracts and summarizes coefficients across all models in an ensemble.
#' For LASSO models, reports selection frequency. For GLM models, reports
#' significance frequency (proportion of models where p < 0.05) and mean
#' -log10(p-value). Both model types include coefficient statistics.
#'
#' @param models List of model objects, or output from \code{load_PI_models()}.
#' @param model_type Character. "auto" (detect), "LASSO", or "GLM".
#'
#' @return List with components:
#'   \item{summary}{Data frame with per-feature statistics. Columns differ by
#'     model type: LASSO has n_selected, pct_selected, mean_coef_nonzero; GLM
#'     has n_significant, pct_significant, mean_neg_log10_pval. Both have
#'     mean_coef, sd_coef, median_coef.}
#'   \item{coef_matrix}{Matrix (n_models x n_features) of coefficients}
#'   \item{pval_matrix}{Matrix (n_models x n_features) of p-values (GLM only,
#'     NULL for LASSO)}
#'   \item{model_type}{Detected or specified model type}
#'   \item{n_models}{Number of models summarized}
#'
#' @details
#' Raw coefficient magnitudes are not directly comparable across features with
#' different scales (e.g., CADD ~0-50 vs LINSIGHT 0-1). For cross-feature
#' importance comparisons, use the scale-invariant metrics:
#' \itemize{
#'   \item LASSO: \code{pct_selected} (selection frequency)
#'   \item GLM: \code{pct_significant} (significance frequency, p < 0.05)
#' }
#' Coefficient statistics (mean_coef, sd_coef, median_coef) are included for
#' understanding effect direction and within-feature variability.
#'
#' @examples
#' \dontrun{
#' models <- load_PI_models("path/to/piModels/")
#' coef_summary <- summarize_PI_coefficients(models)
#' print(coef_summary$summary)
#' }
#'
#' @export
summarize_PI_coefficients <- function(models, model_type = "auto") {

  # Handle input from load_PI_models() or raw list
  if (is.list(models) && "models" %in% names(models)) {
    if (model_type == "auto" && "model_type" %in% names(models)) {
      model_type <- models$model_type
    }
    model_list <- models$models
  } else {
    model_list <- models
  }

  n_models <- length(model_list)

  # Detect model type if auto
  if (model_type == "auto") {
    model_type <- .detect_PI_model_type(model_list[[1]])
  }

  # Extract coefficients (and p-values for GLM) from all models
  pval_matrix <- NULL

  if (model_type == "LASSO") {
    # LASSO: extract coefficients only
    coef_list <- lapply(model_list, .extract_coef_lasso)
    coef_matrix <- do.call(rbind, coef_list)
  } else {
    # GLM: extract both coefficients and p-values
    extracted <- lapply(model_list, .extract_coef_pval_glm)
    coef_matrix <- do.call(rbind, lapply(extracted, `[[`, "coef"))
    pval_matrix <- do.call(rbind, lapply(extracted, `[[`, "pval"))
  }

  # Compute summary statistics
  n_features <- ncol(coef_matrix)
  feature_names <- colnames(coef_matrix)

  # Mean and SD across all models
  mean_coef <- colMeans(coef_matrix)
  sd_coef <- apply(coef_matrix, 2, sd)
  median_coef <- apply(coef_matrix, 2, median)

  # Build summary data frame based on model type
  if (model_type == "LASSO") {
    # LASSO: selection frequency (nonzero coefficients)
    n_selected <- colSums(coef_matrix != 0)
    pct_selected <- round(100 * n_selected / n_models, 1)

    # Mean coefficient when nonzero
    mean_coef_nonzero <- sapply(seq_len(n_features), function(j) {
      vals <- coef_matrix[coef_matrix[, j] != 0, j]
      if (length(vals) == 0) NA_real_ else mean(vals)
    })

    summary_df <- data.frame(
      feature = feature_names,
      n_selected = n_selected,
      pct_selected = pct_selected,
      mean_coef = round(mean_coef, 6),
      sd_coef = round(sd_coef, 6),
      median_coef = round(median_coef, 6),
      mean_coef_nonzero = round(mean_coef_nonzero, 6),
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    # Sort by selection frequency (descending)
    summary_df <- summary_df[order(-summary_df$pct_selected), ]

  } else {
    # GLM: significance frequency (p < 0.05)
    n_significant <- colSums(pval_matrix < 0.05, na.rm = TRUE)
    pct_significant <- round(100 * n_significant / n_models, 1)

    # Mean -log10(p-value) for importance ranking
    # Use pmax to avoid -Inf from p=0 (set floor at 1e-300)
    neg_log10_pval <- -log10(pmax(pval_matrix, 1e-300))
    mean_neg_log10_pval <- colMeans(neg_log10_pval, na.rm = TRUE)

    summary_df <- data.frame(
      feature = feature_names,
      n_significant = n_significant,
      pct_significant = pct_significant,
      mean_neg_log10_pval = round(mean_neg_log10_pval, 3),
      mean_coef = round(mean_coef, 6),
      sd_coef = round(sd_coef, 6),
      median_coef = round(median_coef, 6),
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    # Sort by significance frequency (descending)
    summary_df <- summary_df[order(-summary_df$pct_significant), ]
  }

  rownames(summary_df) <- NULL

  return(list(
    summary = summary_df,
    coef_matrix = coef_matrix,
    pval_matrix = pval_matrix,
    model_type = model_type,
    n_models = n_models
  ))
}


#' Plot PI Coefficient Summary
#'
#' @description
#' Visualizes feature importance across the ensemble using scale-invariant
#' metrics. For LASSO models, shows selection frequency (percent of models
#' where coefficient is nonzero). For GLM models, shows significance frequency
#' (percent of models where p < 0.05).
#'
#' @param coef_summary Output from \code{summarize_PI_coefficients()}.
#' @param main Character. Plot title. NULL for auto-generated title based on
#'   model type.
#' @param label_angle Numeric. Rotation angle (degrees) for x-axis labels.
#'   Default: 45. Use 0 for horizontal, 90 for vertical. Useful for adjusting
#'   label readability in publication figures.
#' @param ... Additional parameters passed to \code{barplot()}.
#'
#' @return Invisible NULL. Plot is drawn to current graphics device.
#'
#' @details
#' Coefficient boxplots were intentionally excluded because raw coefficient
#' magnitudes are not comparable across features with different scales (e.g.,
#' CADD ~0-50 vs LINSIGHT 0-1). The bar plot uses scale-invariant importance
#' metrics that allow meaningful cross-feature comparison.
#'
#' @examples
#' \dontrun{
#' models <- load_PI_models("path/to/piModels/")
#' coef_summary <- summarize_PI_coefficients(models)
#' plot_PI_coefficient_summary(coef_summary)
#' }
#'
#' @export
plot_PI_coefficient_summary <- function(coef_summary, main = NULL, label_angle = 45, ...) {

  summary_df <- coef_summary$summary
  n_models <- coef_summary$n_models
  model_type <- coef_summary$model_type

  # Set up plot margins for rotated labels
  old_par <- par(mar = c(8, 4, 4, 2))
  on.exit(par(old_par))

  # Determine bar values and labels based on model type
  if (model_type == "LASSO") {
    bar_values <- summary_df$pct_selected
    bar_title <- if (is.null(main)) {
      sprintf("Selection Frequency (%d LASSO models)", n_models)
    } else main
    ylab_text <- "Selection Frequency (%)"
  } else {
    # GLM: use significance frequency (p < 0.05)
    bar_values <- summary_df$pct_significant
    bar_title <- if (is.null(main)) {
      sprintf("Significance Frequency (%d GLM models)", n_models)
    } else main
    ylab_text <- "% Models with p < 0.05"
  }

  # Clean feature names (remove "x." prefix added by R matrix conversion)
  feature_labels <- gsub("^x\\.", "", summary_df$feature)

  # Create bar plot without x-axis labels (we'll add them manually with angle)
  bp <- barplot(bar_values,
                names.arg = "",  # Suppress default labels
                col = "steelblue",
                main = bar_title,
                ylab = ylab_text,
                ...)

  # Add rotated x-axis labels
  text(x = bp, y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
       labels = feature_labels,
       srt = label_angle,
       adj = 1,   # Right-align at rotation point
       xpd = TRUE,
       cex = 0.8)

  # Add value labels on bars
  text(bp, bar_values + max(bar_values) * 0.02,
       labels = paste0(round(bar_values, 1), "%"),
       cex = 0.7, pos = 3)

  invisible(NULL)
}


#' Predict PI Using Ensemble of Models
#'
#' @description
#' Applies all models in the ensemble to new annotation data and computes
#' both per-model predictions and ensemble-averaged PI values.
#'
#' @param models List of model objects, or output from \code{load_PI_models()}.
#' @param newdata Matrix or data.frame of annotations (N variants x M features).
#'   Column names should match the features used in model training.
#' @param model_type Character. "auto" (detect), "LASSO", or "GLM".
#' @param na_impute Logical. If TRUE, impute NA values before prediction.
#'   Default: TRUE.
#' @param reference_medians Named numeric vector of reference medians for NA
#'   imputation (optional). Typically from \code{compute_annotation_medians()}.
#'   When provided, NAs are imputed with these chromosome-wide medians instead
#'   of per-gene column medians. Names must match column names in \code{newdata}.
#'
#' @return List with components:
#'   \item{predictions}{Matrix (N x n_models) of per-model predictions}
#'   \item{ensemble_pi}{Numeric vector (N) of averaged predictions}
#'   \item{n_models}{Number of models used}
#'   \item{na_imputed}{Logical, whether imputation was applied}
#'
#' @examples
#' \dontrun{
#' models <- load_PI_models("path/to/piModels/")
#' new_annot <- read.csv("new_variants_annotated.csv")
#' pred <- predict_PI_ensemble(models, new_annot)
#' hist(pred$ensemble_pi)
#' }
#'
#' @export
predict_PI_ensemble <- function(models, newdata, model_type = "auto",
                                 na_impute = TRUE,
                                 reference_medians = NULL) {

  # Handle input from load_PI_models() or raw list
  if (is.list(models) && "models" %in% names(models)) {
    if (model_type == "auto" && "model_type" %in% names(models)) {
      model_type <- models$model_type
    }
    model_list <- models$models
  } else {
    model_list <- models
  }

  n_models <- length(model_list)

  # Detect model type if auto
  if (model_type == "auto") {
    model_type <- .detect_PI_model_type(model_list[[1]])
  }

  # Convert to matrix
  if (is.data.frame(newdata)) {
    newdata <- as.matrix(newdata)
  }

  # Impute NA values if requested
  na_imputed <- FALSE
  if (na_impute && any(is.na(newdata))) {
    if (!is.null(reference_medians)) {
      # Use chromosome-wide reference medians for imputation
      for (j in seq_len(ncol(newdata))) {
        na_idx <- is.na(newdata[, j])
        if (any(na_idx)) {
          col_name <- colnames(newdata)[j]
          if (!is.null(col_name) && col_name %in% names(reference_medians) &&
              !is.na(reference_medians[[col_name]])) {
            newdata[na_idx, j] <- reference_medians[[col_name]]
          } else {
            # Fall back to per-gene median if reference not available
            newdata[na_idx, j] <- median(newdata[!na_idx, j], na.rm = TRUE)
          }
        }
      }
    } else {
      # Original behavior: per-gene column medians
      newdata <- impute_na_median(newdata)
    }
    na_imputed <- TRUE
  }

  # Add "x." prefix to column names if needed (glmnet expects this format)
  # Check if the first model's features have x. prefix
  if (!is.null(colnames(newdata))) {
    model_features <- rownames(model_list[[1]]$beta)
    if (!is.null(model_features) && any(grepl("^x\\.", model_features))) {
      # Model has x. prefix, check if data needs it
      if (!any(grepl("^x\\.", colnames(newdata)))) {
        colnames(newdata) <- paste0("x.", colnames(newdata))
      }
    }
  }

  n_variants <- nrow(newdata)

  # Preallocate prediction matrix
  predictions <- matrix(NA_real_, nrow = n_variants, ncol = n_models)

  # Generate predictions from each model
  for (i in seq_len(n_models)) {
    model <- model_list[[i]]

    if (model_type == "LASSO") {
      # For LASSO: use predict.glmnet
      predictions[, i] <- as.numeric(predict(
        model,
        newx = newdata,
        type = "response",
        s = model$lambda
      ))
    } else {
      # For GLM: use predict.glm
      predictions[, i] <- as.numeric(predict(
        model,
        newdata = list(x = newdata),
        type = "response"
      ))
    }
  }

  # Compute ensemble prediction (row means)
  ensemble_pi <- rowMeans(predictions, na.rm = TRUE)

  return(list(
    predictions = predictions,
    ensemble_pi = ensemble_pi,
    n_models = n_models,
    na_imputed = na_imputed
  ))
}


#' Evaluate PI Model Performance
#'
#' @description
#' Evaluates the discriminative performance of PI models using AUC (Area Under
#' the ROC Curve). Computes per-model AUC and ensemble AUC on case/control data.
#'
#' @param models List of model objects, or output from \code{load_PI_models()}.
#' @param cases Matrix, data.frame, or file path to case (positive) annotations.
#' @param controls Matrix, data.frame, file path, or directory for control
#'   (negative) annotations. If a directory, loads all chr*.csv files.
#' @param model_type Character. "auto" (detect), "LASSO", or "GLM".
#' @param features Character vector. Feature columns to use. NULL uses default
#'   11 PI features.
#' @param max_controls Integer. Maximum controls to load from directory (memory
#'   limit). Uses proportional sampling per chromosome. NULL loads all.
#' @param n_eval_controls Integer. Number of controls to use in evaluation after
#'   loading. NULL uses all loaded controls.
#' @param random_seed Integer. For reproducible control sampling.
#'
#' @return List with components:
#'   \item{per_model}{Data frame with model_id and auc columns}
#'   \item{ensemble}{List with auc and predictions}
#'   \item{summary}{List with mean_auc, sd_auc, min_auc, max_auc}
#'   \item{roc_data}{List with labels, per-model predictions, ensemble predictions}
#'   \item{metadata}{Evaluation settings and data info}
#'
#' @examples
#' \dontrun{
#' models <- load_PI_models("path/to/piModels/")
#' eval_result <- evaluate_PI_models(
#'   models,
#'   cases = "cases_annotated.csv",
#'   controls = "controls_dir/",
#'   max_controls = 10000
#' )
#' cat("Ensemble AUC:", eval_result$ensemble$auc, "\n")
#' }
#'
#' @export
evaluate_PI_models <- function(models, cases, controls,
                                model_type = "auto",
                                features = NULL,
                                max_controls = NULL,
                                n_eval_controls = NULL,
                                random_seed = NULL) {

  # Set seed if provided
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  # Handle models input
  if (is.list(models) && "models" %in% names(models)) {
    if (model_type == "auto" && "model_type" %in% names(models)) {
      model_type <- models$model_type
    }
    model_list <- models$models
  } else {
    model_list <- models
  }

  n_models <- length(model_list)

  # Detect model type if auto
  if (model_type == "auto") {
    model_type <- .detect_PI_model_type(model_list[[1]])
  }

  # Default features
  if (is.null(features)) {
    features <- .default_PI_features()
  }

  # Load case annotations
  case_mat <- .load_evaluation_annotations(cases, features, max_n = NULL)
  n_cases <- nrow(case_mat)

  # Load control annotations (with memory limit)
  ctrl_mat <- .load_evaluation_annotations(controls, features, max_n = max_controls)

  # Subsample controls if requested
  if (!is.null(n_eval_controls) && nrow(ctrl_mat) > n_eval_controls) {
    ctrl_idx <- sample(nrow(ctrl_mat), n_eval_controls, replace = FALSE)
    ctrl_mat <- ctrl_mat[ctrl_idx, , drop = FALSE]
  }
  n_controls <- nrow(ctrl_mat)

  # Combine case and control data
  eval_data <- rbind(case_mat, ctrl_mat)
  labels <- c(rep(1, n_cases), rep(0, n_controls))

  # Impute NA values
  eval_data <- impute_na_median(eval_data)

  # Get predictions using ensemble
  pred_result <- predict_PI_ensemble(model_list, eval_data,
                                      model_type = model_type,
                                      na_impute = FALSE)  # Already imputed

  # Compute per-model AUC
  per_model_auc <- numeric(n_models)
  for (i in seq_len(n_models)) {
    per_model_auc[i] <- .compute_auc(labels, pred_result$predictions[, i])
  }

  # Compute ensemble AUC
  ensemble_auc <- .compute_auc(labels, pred_result$ensemble_pi)

  # Build results
  per_model_df <- data.frame(
    model_id = seq_len(n_models),
    auc = round(per_model_auc, 4)
  )

  summary_stats <- list(
    mean_auc = mean(per_model_auc, na.rm = TRUE),
    sd_auc = sd(per_model_auc, na.rm = TRUE),
    min_auc = min(per_model_auc, na.rm = TRUE),
    max_auc = max(per_model_auc, na.rm = TRUE)
  )

  roc_data <- list(
    labels = labels,
    predictions = pred_result$predictions,
    ensemble_pi = pred_result$ensemble_pi
  )

  metadata <- list(
    n_cases = n_cases,
    n_controls = n_controls,
    n_models = n_models,
    model_type = model_type,
    features = features,
    max_controls = max_controls,
    n_eval_controls = n_eval_controls,
    random_seed = random_seed
  )

  return(list(
    per_model = per_model_df,
    ensemble = list(
      auc = ensemble_auc,
      predictions = pred_result$ensemble_pi
    ),
    summary = summary_stats,
    roc_data = roc_data,
    metadata = metadata
  ))
}


#' Plot ROC Curves for PI Models
#'
#' @description
#' Plots ROC (Receiver Operating Characteristic) curves for individual models
#' and the ensemble using base R graphics.
#'
#' @param evaluation_result Output from \code{evaluate_PI_models()}.
#' @param show_individual Logical. Show individual model curves. Default: TRUE.
#' @param highlight_ensemble Logical. Highlight ensemble curve. Default: TRUE.
#' @param col_individual Color for individual model curves. Default: "gray".
#' @param col_ensemble Color for ensemble curve. Default: "red".
#' @param alpha_individual Numeric (0-1). Transparency for individual curves.
#'   Default: 0.3.
#' @param main Character. Plot title. NULL for default.
#' @param ... Additional parameters passed to plot().
#'
#' @return Invisible list with ROC coordinates for further customization:
#'   \item{individual}{List of per-model ROC coordinates (fpr, tpr)}
#'   \item{ensemble}{Ensemble ROC coordinates (fpr, tpr)}
#'
#' @examples
#' \dontrun{
#' eval_result <- evaluate_PI_models(models, cases, controls)
#' plot_PI_roc(eval_result)
#' }
#'
#' @export
plot_PI_roc <- function(evaluation_result, show_individual = TRUE,
                         highlight_ensemble = TRUE,
                         col_individual = "gray",
                         col_ensemble = "red",
                         alpha_individual = 0.3,
                         main = NULL, ...) {

  labels <- evaluation_result$roc_data$labels
  predictions <- evaluation_result$roc_data$predictions
  ensemble_pi <- evaluation_result$roc_data$ensemble_pi
  n_models <- ncol(predictions)
  ensemble_auc <- evaluation_result$ensemble$auc
  mean_auc <- evaluation_result$summary$mean_auc

  # Default title
  if (is.null(main)) {
    main <- sprintf("ROC Curves (%d models)", n_models)
  }

  # Compute ROC coordinates for all models
  roc_individual <- lapply(seq_len(n_models), function(i) {
    .compute_roc_coords(labels, predictions[, i])
  })

  # Compute ensemble ROC
  roc_ensemble <- .compute_roc_coords(labels, ensemble_pi)

  # Create semi-transparent color for individual curves
  col_ind_alpha <- adjustcolor(col_individual, alpha.f = alpha_individual)

  # Set up plot
  plot(c(0, 1), c(0, 1), type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate",
       ylab = "True Positive Rate",
       main = main, ...)

  # Diagonal reference line
  abline(0, 1, lty = 2, col = "darkgray")

  # Plot individual model curves
  if (show_individual) {
    for (i in seq_len(n_models)) {
      lines(roc_individual[[i]]$fpr, roc_individual[[i]]$tpr,
            col = col_ind_alpha, lwd = 0.5)
    }
  }

  # Plot ensemble curve (bold)
  if (highlight_ensemble) {
    lines(roc_ensemble$fpr, roc_ensemble$tpr,
          col = col_ensemble, lwd = 2)
  }

  # Add legend
  legend_text <- c()
  legend_col <- c()
  legend_lwd <- c()

  if (highlight_ensemble) {
    legend_text <- c(legend_text, sprintf("Ensemble (AUC=%.3f)", ensemble_auc))
    legend_col <- c(legend_col, col_ensemble)
    legend_lwd <- c(legend_lwd, 2)
  }

  if (show_individual) {
    legend_text <- c(legend_text, sprintf("Individual (mean AUC=%.3f)", mean_auc))
    legend_col <- c(legend_col, col_individual)
    legend_lwd <- c(legend_lwd, 1)
  }

  if (length(legend_text) > 0) {
    legend("bottomright", legend = legend_text,
           col = legend_col, lwd = legend_lwd, bg = "white")
  }

  invisible(list(
    individual = roc_individual,
    ensemble = roc_ensemble
  ))
}


#################### INTERNAL HELPER FUNCTIONS ####################


#' Detect PI Model Type
#'
#' @description
#' Determines whether a model is LASSO (glmnet) or GLM based on its class.
#'
#' @param model A fitted model object.
#'
#' @return Character: "LASSO" or "GLM".
#'
#' @keywords internal
#' @noRd
.detect_PI_model_type <- function(model) {
  if (inherits(model, "elnet") || inherits(model, "glmnet")) {
    return("LASSO")
  } else if (inherits(model, "glm")) {
    return("GLM")
  } else {
    stop("Unknown model type: ", paste(class(model), collapse = ", "),
         ". Expected glmnet/elnet (LASSO) or glm (GLM).")
  }
}


#' Extract Coefficients from LASSO Model
#'
#' @description
#' Extracts coefficient vector from a glmnet LASSO model object.
#'
#' @param model A fitted glmnet model object.
#'
#' @return Named numeric vector of coefficients (excluding intercept).
#'
#' @keywords internal
#' @noRd
.extract_coef_lasso <- function(model) {
  # model$beta is a sparse matrix (features x 1)
  beta <- as.numeric(model$beta)
  names(beta) <- rownames(model$beta)
  return(beta)
}


#' Extract Coefficients from GLM Model
#'
#' @description
#' Extracts coefficient vector from a glm model object.
#'
#' @param model A fitted glm model object.
#'
#' @return Named numeric vector of coefficients (excluding intercept).
#'
#' @keywords internal
#' @noRd
.extract_coef_glm <- function(model) {
  coef_all <- coef(model)
  # Remove intercept
  coef_all <- coef_all[!names(coef_all) %in% c("(Intercept)")]
  # Remove "x." prefix if present (from glm fitting with matrix input)
  names(coef_all) <- gsub("^x\\.", "", names(coef_all))
  return(coef_all)
}


#' Extract Coefficients and P-values from GLM Model
#'
#' @description
#' Extracts both coefficient values and p-values from a glm model object.
#' Used for computing significance-based importance metrics for GLM models.
#'
#' @param model A fitted glm model object.
#'
#' @return List with components:
#'   \item{coef}{Named numeric vector of coefficients (excluding intercept)}
#'   \item{pval}{Named numeric vector of p-values (excluding intercept)}
#'
#' @keywords internal
#' @noRd
.extract_coef_pval_glm <- function(model) {
  # Get coefficient summary table (includes Estimate, Std. Error, z value, Pr(>|z|))
  coef_table <- summary(model)$coefficients

  # Identify non-intercept rows
  coef_rows <- !rownames(coef_table) %in% c("(Intercept)")

  # Extract coefficients and p-values
  coef_vals <- coef_table[coef_rows, "Estimate"]
  pval_vals <- coef_table[coef_rows, "Pr(>|z|)"]

  # Remove "x." prefix from names (from glm fitting with matrix input)
  clean_names <- gsub("^x\\.", "", names(coef_vals))
  names(coef_vals) <- clean_names
  names(pval_vals) <- clean_names

  return(list(coef = coef_vals, pval = pval_vals))
}


#' Compute AUC (Area Under ROC Curve)
#'
#' @description
#' Computes the Area Under the ROC Curve using the trapezoidal rule.
#' No external dependencies required.
#'
#' @param labels Numeric vector of true binary labels (0 or 1).
#' @param predictions Numeric vector of predicted probabilities.
#'
#' @return Numeric AUC value in \eqn{[0, 1]}.
#'
#' @keywords internal
#' @noRd
.compute_auc <- function(labels, predictions) {

  # Handle NA values
  valid_idx <- !is.na(labels) & !is.na(predictions)
  labels <- labels[valid_idx]
  predictions <- predictions[valid_idx]

  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)

  if (n_pos == 0 || n_neg == 0) {
    warning("AUC undefined: need both positive and negative cases")
    return(NA_real_)
  }

  # Sort by predictions descending
  ord <- order(predictions, decreasing = TRUE)
  labels <- labels[ord]

  # Compute cumulative TPR and FPR
  tpr <- cumsum(labels == 1) / n_pos
  fpr <- cumsum(labels == 0) / n_neg

  # Prepend (0, 0) for complete ROC curve
  fpr <- c(0, fpr)
  tpr <- c(0, tpr)

  # Trapezoidal integration
  # AUC = sum of trapezoids: 0.5 * (fpr[i+1] - fpr[i]) * (tpr[i+1] + tpr[i])
  n <- length(fpr)
  auc <- sum(diff(fpr) * (tpr[-1] + tpr[-n]) / 2)

  return(auc)
}


#' Compute ROC Curve Coordinates
#'
#' @description
#' Computes the coordinates (FPR, TPR) for plotting an ROC curve.
#'
#' @param labels Numeric vector of true binary labels (0 or 1).
#' @param predictions Numeric vector of predicted probabilities.
#'
#' @return List with components fpr (False Positive Rate) and tpr (True Positive
#'   Rate), both numeric vectors.
#'
#' @keywords internal
#' @noRd
.compute_roc_coords <- function(labels, predictions) {

  # Handle NA values
  valid_idx <- !is.na(labels) & !is.na(predictions)
  labels <- labels[valid_idx]
  predictions <- predictions[valid_idx]

  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)

  if (n_pos == 0 || n_neg == 0) {
    warning("ROC undefined: need both positive and negative cases")
    return(list(fpr = c(0, 1), tpr = c(0, 1)))
  }

  # Sort by predictions descending
  ord <- order(predictions, decreasing = TRUE)
  labels <- labels[ord]

  # Compute cumulative rates
  tpr <- c(0, cumsum(labels == 1) / n_pos)
  fpr <- c(0, cumsum(labels == 0) / n_neg)

  return(list(fpr = fpr, tpr = tpr))
}


#' Load Evaluation Annotations
#'
#' @description
#' Loads annotation data for evaluation from various sources: matrix, data.frame,
#' file path (CSV), or directory (multiple CSVs).
#'
#' @param source Matrix, data.frame, file path, or directory.
#' @param features Character vector of feature columns to extract.
#' @param max_n Integer. Maximum rows to load (for memory limit). NULL for all.
#'
#' @return Matrix of annotations (N x M).
#'
#' @keywords internal
#' @noRd
.load_evaluation_annotations <- function(source, features, max_n = NULL) {

  # If already a matrix
  if (is.matrix(source)) {
    mat <- source
  }
  # If data.frame
  else if (is.data.frame(source)) {
    mat <- as.matrix(source)
  }
  # If directory path
  else if (is.character(source) && dir.exists(source)) {
    # Use existing load_control_annotations if available
    if (exists("load_control_annotations", mode = "function")) {
      mat <- load_control_annotations(
        source = source,
        features = features,
        max_controls = max_n,
        format = "csv"
      )
    } else {
      # Fallback: load CSVs manually
      csv_files <- list.files(source, pattern = "^chr.*\\.csv$", full.names = TRUE)
      if (length(csv_files) == 0) {
        stop("No chr*.csv files found in directory: ", source)
      }

      # Load and combine
      data_list <- lapply(csv_files, function(f) {
        df <- data.table::fread(f, data.table = FALSE)
        if (!all(features %in% names(df))) {
          missing <- setdiff(features, names(df))
          stop("Missing features in ", f, ": ", paste(missing, collapse = ", "))
        }
        as.matrix(df[, features, drop = FALSE])
      })
      mat <- do.call(rbind, data_list)

      # Sample if max_n specified
      if (!is.null(max_n) && nrow(mat) > max_n) {
        mat <- mat[sample(nrow(mat), max_n, replace = FALSE), , drop = FALSE]
      }
    }
    return(mat)
  }
  # If file path
  else if (is.character(source) && file.exists(source)) {
    df <- data.table::fread(source, data.table = FALSE)
    if (!all(features %in% names(df))) {
      missing <- setdiff(features, names(df))
      stop("Missing features in ", source, ": ", paste(missing, collapse = ", "))
    }
    mat <- as.matrix(df[, features, drop = FALSE])
  }
  else {
    stop("Invalid source: must be matrix, data.frame, file path, or directory")
  }

  # Extract features if needed
  if (!is.null(features) && ncol(mat) != length(features)) {
    if (all(features %in% colnames(mat))) {
      mat <- mat[, features, drop = FALSE]
    }
  }

  # Sample if max_n specified
  if (!is.null(max_n) && nrow(mat) > max_n) {
    mat <- mat[sample(nrow(mat), max_n, replace = FALSE), , drop = FALSE]
  }

  return(mat)
}
