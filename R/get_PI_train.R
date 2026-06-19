#' PI Model Training Functions
#'
#' Functions for training PI (variant-importance score) ensemble models
#' from pre-annotated case and control variant data.
#'
#' File Log (reverse chronological order):
#' - 2026-06-11: Modified by Claude Code (Opus 4.8) via the r-developer agent,
#'   prompted by ZWu -- annotate_favor aGDS-format fix: .load_control_gds() now
#'   reads EITHER the STAARpipeline sub-node folder format (per-feature sub-nodes,
#'   the format annotate_favor now writes) OR the legacy single-matrix format
#'   (read.gdsn + feature_names attr). Branches on length(ls.gdsn(node)) > 0 before
#'   any read.gdsn (which errors on a folder node). Legacy matrix fixtures still load.
#' - 2026-02-03: Renamed by Claude Code - train_PI_models.R -> get_PI_train.R
#'   for consistent naming with other get_PI_*.R files
#' - 2026-02-02: Modified by Claude Code - Changed controls_per_model to default
#'   NULL (auto-calculate from n_cases); added controls_multiplier parameter
#' - 2026-02-02: Modified by Claude Code - Fixed .load_control_gds() to use
#'   gdsfmt::read.gdsn for matrix-based FunctionalAnnotation; added GDS
#'   proportional sampling (.load_gds_with_proportional_sampling)
#' - 2026-02-02: Modified by Claude Code - Implemented proportional sampling
#'   during file loading to limit memory; default max_controls=50000
#' - 2026-02-02: Modified by Claude Code - Reorganized per project rules:
#'   main functions first, helpers last; added max_controls parameter
#' - 2026-02-02: Created by Claude Code - Initial implementation
#'
#' @name get_PI_train
#' @docType package
NULL


#################### EXPORTED MAIN FUNCTIONS ####################


#' Train PI Ensemble Models
#'
#' @description
#' High-level function for training PI (variant-importance score) ensemble models.
#' Loads pre-annotated case and control data, imputes missing values, trains
#' models, and saves them to disk.
#'
#' @param case_csv Path to the annotated case variants CSV file.
#' @param control_source Path to control data. Either:
#'   \itemize{
#'     \item A directory containing \code{chr*.csv} or \code{chr*.gds} files
#'     \item A single CSV or GDS file
#'   }
#' @param output_dir Directory where model RDS files will be saved.
#' @param n_models Number of ensemble models to train. Default: 100.
#' @param controls_per_model Number of controls sampled per model. Default: NULL
#'   (automatically set to match the number of cases times \code{controls_multiplier}).
#'   If specified explicitly, overrides the automatic calculation.
#' @param controls_multiplier Multiplier for automatic \code{controls_per_model}
#'   calculation. Only used when \code{controls_per_model = NULL}. Default: 1
#'   (equal number of controls and cases). Set to 2 for twice as many controls
#'   as cases, etc.
#' @param max_controls Maximum number of control variants to load. If the source
#'   contains more variants, proportional sampling is performed during loading
#'   to limit memory usage. Default: 50000. Set to NULL to load all variants
#'   (not recommended for large datasets).
#' @param model_type Model type: "GLM" (default) or "LASSO".
#' @param features Character vector of annotation features to use.
#'   Default: NULL (uses 11 default FAVOR features).
#' @param chromosomes Integer vector of chromosomes to include for controls.
#'   Default: 1:22.
#' @param random_seed Random seed for reproducibility. Applied BEFORE any
#'   sampling operations (control loading and model training). Default: NULL.
#' @param verbose Verbosity level: 0 (silent), 1 (progress), 2 (detailed).
#'
#' @return A list with components:
#'   \describe{
#'     \item{models}{List of fitted model objects (length = n_models)}
#'     \item{metadata}{List with training metadata: n_cases, n_controls,
#'       n_models, controls_per_model, controls_multiplier, model_type,
#'       features, output_dir, random_seed}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Set random seed (if provided) for reproducibility
#'   \item Load case annotations from CSV
#'   \item Load control annotations from CSV/GDS (samples if max_controls set)
#'   \item Calculate controls_per_model if not specified (n_cases * multiplier)
#'   \item Impute NA values with column median
#'   \item Train ensemble using \code{model_PI()}
#'   \item Save models to output directory
#' }
#'
#' @section Controls Per Model:
#' By default, the number of controls sampled per model matches the number of
#' cases (\code{controls_multiplier = 1}). This creates balanced training sets.
#' Use \code{controls_multiplier = 2} or higher to oversample controls, which
#' may improve model stability when controls are abundant.
#'
#' @section Reproducibility:
#' When \code{random_seed} is set, all random operations are reproducible:
#' \itemize{
#'   \item Control variant sampling (if max_controls < total controls)
#'   \item Control selection within each model
#' }
#' Running with the same seed produces identical models.
#'
#' @section Memory Usage:
#' For large control datasets (millions of variants), use \code{max_controls}
#' to limit memory usage. A value of 100,000-500,000 is typically sufficient
#' since each model only samples \code{controls_per_model} variants.
#'
#' @examples
#' \dontrun{
#' # Default: controls_per_model = number of cases
#' result <- train_PI_models(
#'   case_csv = "als_cases_annotated.csv",
#'   control_source = "/path/to/favor_csv/",
#'   output_dir = "piModels/als/",
#'   n_models = 100,
#'   random_seed = 42
#' )
#'
#' # Use 2x controls per case
#' result <- train_PI_models(
#'   case_csv = "als_cases_annotated.csv",
#'   control_source = "/path/to/favor_csv/",
#'   output_dir = "piModels/als_2x/",
#'   controls_multiplier = 2,
#'   random_seed = 42
#' )
#'
#' # Explicit controls_per_model (overrides multiplier)
#' result <- train_PI_models(
#'   case_csv = "als_cases_annotated.csv",
#'   control_source = "/path/to/favor_csv/",
#'   output_dir = "piModels/als_500ctrls/",
#'   controls_per_model = 500,
#'   random_seed = 42
#' )
#' }
#'
#' @export
train_PI_models <- function(case_csv,
                            control_source,
                            output_dir,
                            n_models = 100,
                            controls_per_model = NULL,
                            controls_multiplier = 1,
                            max_controls = 50000,
                            model_type = "GLM",
                            features = NULL,
                            chromosomes = 1:22,
                            random_seed = NULL,
                            verbose = 1) {
  # Set random seed FIRST for reproducibility of all subsequent operations
  if (!is.null(random_seed)) {
    set.seed(random_seed)
    if (verbose >= 1) message("Random seed set to: ", random_seed)
  }

  # Use default features if not specified
  if (is.null(features)) {
    features <- .default_PI_features()
  }

  # Validate model_type
  if (!model_type %in% c("GLM", "LASSO")) {
    stop("model_type must be 'GLM' or 'LASSO'")
  }

  # Validate controls_multiplier
 if (!is.null(controls_multiplier)) {
    if (controls_multiplier <= 0) {
      stop("controls_multiplier must be a positive number")
    }
  }

  # Step 1: Load case annotations
  if (verbose >= 1) message("Loading case annotations...")
  case_mat <- load_case_annotations(case_csv, features)
  n_cases <- nrow(case_mat)
  if (verbose >= 1) message("  Loaded ", n_cases, " cases with ",
                            ncol(case_mat), " features.")

  # Step 2: Load control annotations (with optional sampling)
  if (verbose >= 1) message("Loading control annotations...")
  ctrl_mat <- load_control_annotations(
    source = control_source,
    features = features,
    chromosomes = chromosomes,
    max_controls = max_controls
  )
  if (verbose >= 1) message("  Loaded ", nrow(ctrl_mat), " controls.")

  # Calculate controls_per_model if not specified
  if (is.null(controls_per_model)) {
    controls_per_model <- as.integer(n_cases * controls_multiplier)
    if (verbose >= 1) {
      message("  Controls per model: ", controls_per_model,
              " (", n_cases, " cases x ", controls_multiplier, " multiplier)")
    }
  } else {
    if (verbose >= 1) {
      message("  Controls per model: ", controls_per_model, " (explicitly set)")
    }
  }

  # Validate controls_per_model
  if (controls_per_model > nrow(ctrl_mat)) {
    stop("controls_per_model (", controls_per_model, ") exceeds available controls (",
         nrow(ctrl_mat), ").")
  }

  # Step 3: Impute NA values
  if (verbose >= 1) message("Imputing missing values with column median...")
  na_before_case <- sum(is.na(case_mat))
  na_before_ctrl <- sum(is.na(ctrl_mat))

  case_mat <- impute_na_median(case_mat)
  ctrl_mat <- impute_na_median(ctrl_mat)

  if (verbose >= 2 && (na_before_case > 0 || na_before_ctrl > 0)) {
    message("  Imputed ", na_before_case, " case NAs and ",
            na_before_ctrl, " control NAs.")
  }

  # Check for remaining NAs (shouldn't happen unless entire column was NA)
  if (any(is.na(case_mat)) || any(is.na(ctrl_mat))) {
    warning("Some columns had all NA values and could not be imputed.")
  }

  # Step 4: Train ensemble models
  if (verbose >= 1) message("Training ", n_models, " ", model_type, " models...")

  models <- model_PI(
    caseAnnotation = case_mat,
    controlAnnotation = ctrl_mat,
    modelType = model_type,
    control_need_N = controls_per_model,
    model_need_N = n_models
  )

  if (verbose >= 1) message("  Training complete.")

  # Step 5: Save models
  if (verbose >= 1) message("Saving models to: ", output_dir)
  file_paths <- save_PI_models(models, output_dir)
  if (verbose >= 1) message("  Saved ", length(file_paths), " model files.")

  # Prepare metadata
  metadata <- list(
    n_cases = n_cases,
    n_controls = nrow(ctrl_mat),
    n_models = n_models,
    controls_per_model = controls_per_model,
    controls_multiplier = controls_multiplier,
    max_controls = max_controls,
    model_type = model_type,
    features = features,
    output_dir = normalizePath(output_dir),
    random_seed = random_seed,
    case_csv = normalizePath(case_csv),
    control_source = ifelse(dir.exists(control_source),
                            normalizePath(control_source),
                            normalizePath(control_source)),
    timestamp = Sys.time()
  )

  if (verbose >= 1) {
    message("\nSummary:")
    message("  Cases: ", metadata$n_cases)
    message("  Controls loaded: ", metadata$n_controls)
    message("  Controls per model: ", metadata$controls_per_model)
    message("  Models: ", metadata$n_models)
    message("  Type: ", metadata$model_type)
    message("  Output: ", metadata$output_dir)
  }

  list(
    models = models,
    metadata = metadata
  )
}


#' Load Case Annotations from CSV
#'
#' @description
#' Loads annotation features from a pre-annotated case variants CSV file.
#'
#' @param csv_path Path to the annotated CSV file. Must contain a VarInfo column
#'   and the requested feature columns.
#' @param features Character vector of annotation column names to extract.
#'   Default: 11 FAVOR features from \code{.default_PI_features()}.
#'
#' @return Numeric matrix with dimensions (N_variants x M_features).
#'   Row names are set to VarInfo values. May contain NA values.
#'
#' @details
#' The CSV file should have been generated by \code{annotate_favor()} or
#' similar FAVOR annotation tools. Expected columns include VarInfo (variant
#' identifier in CHR-POS-REF-ALT format) and the annotation features.
#'
#' @section NA Handling:
#' This function preserves NA values in the output. Use \code{impute_na_median()}
#' before passing to \code{model_PI()} if needed.
#'
#' @examples
#' \dontrun{
#' case_annot <- load_case_annotations("cases_annotated.csv")
#' dim(case_annot)  # N_cases x 11
#' }
#'
#' @export
load_case_annotations <- function(csv_path,
                                  features = .default_PI_features()) {
  # Validate input file exists
  if (!file.exists(csv_path)) {
    stop("File not found: ", csv_path)
  }

  # Read CSV using data.table for efficiency
  data <- data.table::fread(csv_path, data.table = FALSE)

  # Validate VarInfo column exists
  if (!"VarInfo" %in% names(data)) {
    stop("CSV must contain a 'VarInfo' column for variant identification.")
  }

  # Validate feature columns exist
  missing_features <- setdiff(features, names(data))
  if (length(missing_features) > 0) {
    stop("Missing annotation columns: ", paste(missing_features, collapse = ", "))
  }

  # Extract feature columns as matrix
  annot_mat <- as.matrix(data[, features, drop = FALSE])

  # Set row names to VarInfo for variant identification
  rownames(annot_mat) <- data$VarInfo

  # Report NA summary for user awareness
  na_count <- sum(is.na(annot_mat))
  if (na_count > 0) {
    na_pct <- round(100 * na_count / length(annot_mat), 1)
    message("Note: ", na_count, " NA values (", na_pct, "%) in case annotations.")
  }

  annot_mat
}


#' Load Control Annotations from CSV or GDS
#'
#' @description
#' Loads annotation features from pre-annotated control variants.
#' Supports single file or directory with per-chromosome files.
#'
#' @param source Path to either:
#'   \itemize{
#'     \item A single CSV or GDS file
#'     \item A directory containing \code{chr*.csv} or \code{chr*.gds} files
#'   }
#' @param features Character vector of annotation column names to extract.
#'   Default: 11 FAVOR features.
#' @param chromosomes Integer vector of chromosomes to include when loading
#'   from a directory. Default: 1:22 (all autosomes).
#' @param format File format: "auto" (detect from extension), "csv", or "gds".
#' @param max_controls Maximum number of control variants to return. If source
#'   contains more variants, proportional sampling is performed during file
#'   loading to limit memory usage. Default: NULL (all).
#'   Note: Random seed should be set before calling this function for
#'   reproducible sampling.
#'
#' @return Numeric matrix with dimensions (N_variants x M_features).
#'   Row names are set to VarInfo values. May contain NA values.
#'
#' @details
#' When \code{source} is a directory, the function scans for files matching
#' \code{chr{N}_*.csv} or \code{chr{N}_*.gds} patterns and loads only those
#' corresponding to the requested chromosomes.
#'
#' @section Memory Management:
#' When \code{max_controls} is set and loading from a directory, proportional
#' sampling is performed during loading (not after). Each chromosome file
#' contributes variants proportional to its size, so memory usage is bounded
#' by max_controls regardless of total dataset size.
#'
#' @examples
#' \dontrun{
#' # Load from directory with chr1-22 CSVs (sample 100K)
#' ctrl_annot <- load_control_annotations(
#'   "/path/to/favor_csv/",
#'   max_controls = 100000
#' )
#'
#' # Load specific chromosomes
#' ctrl_annot <- load_control_annotations(
#'   "/path/to/favor_csv/",
#'   chromosomes = c(21, 22)
#' )
#' }
#'
#' @export
load_control_annotations <- function(source,
                                     features = .default_PI_features(),
                                     chromosomes = 1:22,
                                     format = "auto",
                                     max_controls = NULL) {
  # Determine if source is directory or file and load accordingly
  if (dir.exists(source)) {
    # Pass max_controls to directory loader for proportional sampling during load
    annot_mat <- .load_control_from_directory(
      source, features, chromosomes, format, max_controls
    )
  } else if (file.exists(source)) {
    # For single file, load then sample if needed
    annot_mat <- .load_control_from_file(source, features, format)
    if (!is.null(max_controls) && nrow(annot_mat) > max_controls) {
      message("Sampling ", max_controls, " from ", nrow(annot_mat), " controls...")
      sample_idx <- sample(nrow(annot_mat), max_controls, replace = FALSE)
      annot_mat <- annot_mat[sample_idx, , drop = FALSE]
    }
  } else {
    stop("Source not found: ", source)
  }

  # Report NA summary for user awareness
  na_count <- sum(is.na(annot_mat))
  if (na_count > 0) {
    na_pct <- round(100 * na_count / length(annot_mat), 1)
    message("Note: ", na_count, " NA values (", na_pct, "%) in control annotations.")
  }

  annot_mat
}


#' Impute NA Values with Column Median
#'
#' @description
#' Replaces NA values in a numeric matrix with the column median.
#'
#' @param mat Numeric matrix that may contain NA values.
#'
#' @return Numeric matrix with NA values replaced by column medians.
#'
#' @details
#' For each column, NA values are replaced with the median of non-NA values.
#' If a column is entirely NA, it remains NA (a warning should be issued by
#' the caller).
#'
#' @examples
#' mat <- matrix(c(1, 2, NA, 4, NA, 6), ncol = 2)
#' impute_na_median(mat)
#'
#' @export
impute_na_median <- function(mat) {
  # Ensure input is matrix
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }

  # Impute each column with its median

for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) {
      col_median <- median(mat[!na_idx, j], na.rm = TRUE)
      mat[na_idx, j] <- col_median
    }
  }

  mat
}


#' Save PI Models to Directory
#'
#' @description
#' Saves a list of fitted PI models as individual RDS files.
#'
#' @param models List of fitted model objects.
#' @param output_dir Directory path for saving models.
#' @param prefix File name prefix. Default: "model_".
#'
#' @return Character vector of saved file paths (invisibly).
#'
#' @details
#' Creates the output directory if it doesn't exist. Each model is saved
#' as \code{{prefix}{i}.rds} (e.g., model_1.rds, model_2.rds, ...).
#'
#' @examples
#' \dontrun{
#' models <- list(model1, model2, model3)
#' save_PI_models(models, "output/piModels/")
#' }
#'
#' @export
save_PI_models <- function(models, output_dir, prefix = "model_") {
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save each model as separate RDS file
  file_paths <- character(length(models))
  for (i in seq_along(models)) {
    file_path <- file.path(output_dir, paste0(prefix, i, ".rds"))
    saveRDS(models[[i]], file_path)
    file_paths[i] <- file_path
  }

  invisible(file_paths)
}


#################### INTERNAL HELPER FUNCTIONS ####################


#' Default PI Training Features
#'
#' @description
#' Returns the default 16 FAVOR annotation features used for PI model training:
#' 13 Annotation Principal Components (latest version per category, epigenetics
#' expanded) + 3 integrative scores. Uses original FAVOR database column names.
#'
#' This is a curated subset of \code{.default_favor_features()} designed for
#' PI model training. Uses latest version for versioned aPCs and includes
#' epigenetics sub-features to let LASSO perform feature selection.
#'
#' @return Character vector of 16 feature names.
#'
#' @seealso \code{.default_favor_features()} for the full annotation set.
#'
#' @keywords internal
#' @noRd
.default_PI_features <- function() {
  c(
    # 13 Annotation Principal Components (latest version per category)
    "apc_conservation_v2",
    "apc_epigenetics",
    "apc_epigenetics_active",
    "apc_epigenetics_repressed",
    "apc_epigenetics_transcription",
    "apc_protein_function_v3",
    "apc_local_nucleotide_diversity_v3",
    "apc_mutation_density",
    "apc_transcription_factor",
    "apc_mappability",
    "apc_proximity_to_tsstes",
    "apc_proximity_to_coding_v2",
    "apc_micro_rna",
    # 3 Integrative Scores
    "cadd_phred",
    "linsight",
    "fathmm_xf"
  )
}


#' Load Control Annotations from Directory
#'
#' @description
#' Internal function to load and concatenate chromosome files from a directory.
#' Supports proportional sampling during loading to limit memory usage.
#'
#' @param dir_path Directory path containing chr*.csv or chr*.gds files.
#' @param features Feature columns to extract.
#' @param chromosomes Chromosomes to include.
#' @param format File format.
#' @param max_controls Maximum controls to load (NULL = all). When set,
#'   proportional sampling is performed during loading.
#'
#' @return Numeric matrix.
#'
#' @keywords internal
#' @noRd
.load_control_from_directory <- function(dir_path, features, chromosomes, format,
                                         max_controls = NULL) {
  # Auto-detect format from files in directory
  if (format == "auto") {
    csv_files <- list.files(dir_path, pattern = "^chr[0-9]+.*\\.csv$", full.names = TRUE)
    gds_files <- list.files(dir_path, pattern = "^chr[0-9]+.*\\.gds$", full.names = TRUE)

    if (length(csv_files) > 0 && length(gds_files) == 0) {
      format <- "csv"
    } else if (length(gds_files) > 0 && length(csv_files) == 0) {
      format <- "gds"
    } else if (length(csv_files) > 0 && length(gds_files) > 0) {
      message("Both CSV and GDS files found. Using CSV format.")
      format <- "csv"
    } else {
      stop("No chr*.csv or chr*.gds files found in directory: ", dir_path)
    }
  }

  # Get all files of the detected format
  if (format == "csv") {
    all_files <- list.files(dir_path, pattern = "^chr[0-9]+.*\\.csv$", full.names = TRUE)
  } else if (format == "gds") {
    all_files <- list.files(dir_path, pattern = "^chr[0-9]+.*\\.gds$", full.names = TRUE)
  }

  # Filter to requested chromosomes using regex
  chr_pattern <- paste0("chr(", paste(chromosomes, collapse = "|"), ")[^0-9]")
  selected_files <- all_files[grepl(chr_pattern, basename(all_files))]

  if (length(selected_files) == 0) {
    stop("No files found for chromosomes: ", paste(chromosomes, collapse = ", "))
  }

  message("Loading ", length(selected_files), " chromosome file(s)...")

  # Use proportional sampling during loading when max_controls is specified
  if (!is.null(max_controls)) {
    if (format == "csv") {
      combined <- .load_csv_with_proportional_sampling(
        selected_files, features, max_controls
      )
    } else if (format == "gds") {
      combined <- .load_gds_with_proportional_sampling(
        selected_files, features, max_controls
      )
    }
  } else {
    # Load all data from each file (no max_controls limit)
    data_list <- lapply(selected_files, function(f) {
      .load_control_from_file(f, features, format)
    })
    combined <- do.call(rbind, data_list)
  }

  message("Loaded ", nrow(combined), " control variants from ",
          length(selected_files), " files.")

  combined
}


#' Load CSV Files with Proportional Sampling
#'
#' @description
#' Internal function that samples proportionally from each CSV file during
#' loading to limit memory usage. Only reads sampled rows into memory.
#'
#' @param files Character vector of CSV file paths.
#' @param features Feature columns to extract.
#' @param max_controls Total number of controls to sample across all files.
#'
#' @return Numeric matrix with sampled rows.
#'
#' @keywords internal
#' @noRd
.load_csv_with_proportional_sampling <- function(files, features, max_controls) {
  # Step 1: Count rows in each file (fast - just count lines)
  row_counts <- vapply(files, function(f) {
    # Use wc -l for fast line counting (subtract 1 for header)
    result <- tryCatch({
      as.integer(system(paste("wc -l <", shQuote(f)), intern = TRUE)) - 1L
    }, error = function(e) {
      # Fallback: read with nrows to count
      nrow(data.table::fread(f, select = 1L, data.table = FALSE))
    })
    result
  }, integer(1))

  total_rows <- sum(row_counts)

  # Step 2: If total is less than max_controls, load all
  if (total_rows <= max_controls) {
    message("Total controls (", total_rows, ") <= max_controls. Loading all.")
    data_list <- lapply(files, function(f) {
      .load_control_csv(f, features)
    })
    return(do.call(rbind, data_list))
  }

  # Step 3: Calculate proportional allocation for each file
  message("Proportionally sampling ", max_controls, " from ", total_rows, " total controls...")
  proportions <- row_counts / total_rows
  allocations <- round(proportions * max_controls)

  # Adjust to ensure exact total (rounding may cause slight differences)
  diff <- max_controls - sum(allocations)
  if (diff != 0) {
    # Add/subtract from largest file
    largest_idx <- which.max(allocations)
    allocations[largest_idx] <- allocations[largest_idx] + diff
  }

  # Step 4: Load sampled rows from each file
  data_list <- mapply(function(f, n_rows, n_sample) {
    if (n_sample <= 0) {
      return(NULL)
    }
    if (n_sample >= n_rows) {
      # Load all rows from this file
      return(.load_control_csv(f, features))
    }

    # Sample row indices (1-based, excluding header)
    sample_rows <- sort(sample(n_rows, n_sample, replace = FALSE))

    # Read only sampled rows using data.table's skip/nrows is inefficient
    # Instead, read all and subset (still better than concatenating then sampling)
    # For very large files, could use fread with select rows, but complex
    all_data <- data.table::fread(f, data.table = FALSE)

    # Validate columns
    if (!"VarInfo" %in% names(all_data)) {
      stop("CSV must contain 'VarInfo' column: ", f)
    }
    missing_features <- setdiff(features, names(all_data))
    if (length(missing_features) > 0) {
      stop("Missing columns in ", basename(f), ": ",
           paste(missing_features, collapse = ", "))
    }

    # Extract sampled rows
    annot_mat <- as.matrix(all_data[sample_rows, features, drop = FALSE])
    rownames(annot_mat) <- all_data$VarInfo[sample_rows]

    annot_mat
  }, files, row_counts, allocations, SIMPLIFY = FALSE)

  # Remove NULLs and combine
  data_list <- data_list[!vapply(data_list, is.null, logical(1))]
  do.call(rbind, data_list)
}


#' Load Control Annotations from Single File
#'
#' @description
#' Internal function to load annotations from a single CSV or GDS file.
#'
#' @param file_path File path.
#' @param features Feature columns to extract.
#' @param format File format ("auto", "csv", or "gds").
#'
#' @return Numeric matrix.
#'
#' @keywords internal
#' @noRd
.load_control_from_file <- function(file_path, features, format) {
  # Auto-detect format from file extension
  if (format == "auto") {
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      format <- "csv"
    } else if (grepl("\\.gds$", file_path, ignore.case = TRUE)) {
      format <- "gds"
    } else {
      stop("Cannot determine file format. Please specify 'csv' or 'gds'.")
    }
  }

  # Dispatch to appropriate loader
  if (format == "csv") {
    .load_control_csv(file_path, features)
  } else if (format == "gds") {
    .load_control_gds(file_path, features)
  } else {
    stop("Unsupported format: ", format)
  }
}


#' Load Control Annotations from CSV File
#'
#' @description
#' Internal function to load annotations from a CSV file.
#'
#' @param file_path CSV file path.
#' @param features Feature columns to extract.
#'
#' @return Numeric matrix.
#'
#' @keywords internal
#' @noRd
.load_control_csv <- function(file_path, features) {
  # Read CSV efficiently
  data <- data.table::fread(file_path, data.table = FALSE)

  # Validate VarInfo column
  if (!"VarInfo" %in% names(data)) {
    stop("CSV must contain 'VarInfo' column: ", file_path)
  }

  # Validate feature columns
  missing_features <- setdiff(features, names(data))
  if (length(missing_features) > 0) {
    stop("Missing columns in ", basename(file_path), ": ",
         paste(missing_features, collapse = ", "))
  }

  # Extract matrix with VarInfo as row names
  annot_mat <- as.matrix(data[, features, drop = FALSE])
  rownames(annot_mat) <- data$VarInfo

  annot_mat
}


#' Load Control Annotations from GDS File
#'
#' @description
#' Internal function to load annotations from an aGDS file. Supports both the
#' STAARpipeline sub-node \emph{folder} format (one typed sub-node per feature
#' under \code{annotation/info/FunctionalAnnotation/}, as written by
#' \code{annotate_favor()}) and the legacy single-\emph{matrix} format (with a
#' \code{feature_names} attribute).
#'
#' @param file_path GDS file path.
#' @param features Feature columns to extract.
#' @param variant_indices Optional integer vector of variant indices to load
#'   (1-based). If NULL, loads all variants.
#'
#' @return Numeric matrix.
#'
#' @details
#' The node format is detected from the FunctionalAnnotation node:
#' \itemize{
#'   \item \strong{Folder} (\code{length(ls.gdsn(node)) > 0}): each requested
#'     feature's sub-node is read individually and assembled into a data.frame.
#'     \code{read.gdsn()} on a folder node would error ("no data field"), so the
#'     per-sub-node read is required. This is the format the per-feature
#'     variant-set scan also reads.
#'   \item \strong{Matrix} (no sub-nodes): the whole node is read via
#'     \code{read.gdsn()} and columns are named from the \code{feature_names}
#'     attribute (legacy back-compat).
#' }
#'
#' @keywords internal
#' @noRd
.load_control_gds <- function(file_path, features, variant_indices = NULL) {
  # Check if SeqArray and gdsfmt are available
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("Package 'SeqArray' required for GDS support. ",
         "Install with: BiocManager::install('SeqArray')")
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("Package 'gdsfmt' required for GDS support. ",
         "Install with: BiocManager::install('gdsfmt')")
  }

  # Open GDS file read-only
  gds <- SeqArray::seqOpen(file_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Build VarInfo from chromosome, position, ref, alt (all variants first)
  chr <- SeqArray::seqGetData(gds, "chromosome")
  pos <- SeqArray::seqGetData(gds, "position")
  ref <- SeqArray::seqGetData(gds, "$ref")
  alt <- SeqArray::seqGetData(gds, "$alt")
  varinfo <- paste(chr, pos, ref, alt, sep = "-")

  # Try to read FAVOR annotations using gdsfmt (more flexible than seqGetData)
  annot_path <- "annotation/info/FunctionalAnnotation"

  # Check if annotation node exists using gdsfmt
  annot_node <- tryCatch(
    gdsfmt::index.gdsn(gds, annot_path),
    error = function(e) NULL
  )

  if (is.null(annot_node)) {
    stop("FAVOR annotations not found in GDS file at path: ", annot_path)
  }

  # Detect the FunctionalAnnotation node format. A folder node (STAARpipeline
  # sub-node layout, as written by annotate_favor()) exposes one sub-node per
  # feature; a legacy matrix node has no sub-nodes. read.gdsn() on a folder node
  # errors ("no data field"), so the format must be branched on ls.gdsn() BEFORE
  # any read.
  sub_nodes <- gdsfmt::ls.gdsn(annot_node)

  if (length(sub_nodes) > 0) {
    # ---- Folder format: read each requested feature's sub-node ----
    missing_features <- setdiff(features, sub_nodes)
    if (length(missing_features) > 0) {
      stop("Missing annotation columns in GDS: ",
           paste(missing_features, collapse = ", "))
    }
    annot_df <- as.data.frame(
      lapply(features, function(feat) {
        gdsfmt::read.gdsn(gdsfmt::index.gdsn(annot_node, feat))
      }),
      stringsAsFactors = FALSE
    )
    names(annot_df) <- features
  } else {
    # ---- Legacy matrix format: read the whole node, name from attribute ----
    annot_data <- gdsfmt::read.gdsn(annot_node)
    if (is.matrix(annot_data)) {
      attr_list <- gdsfmt::get.attr.gdsn(annot_node)
      if ("feature_names" %in% names(attr_list)) {
        colnames(annot_data) <- attr_list$feature_names
      }
      annot_df <- as.data.frame(annot_data)
    } else if (is.list(annot_data)) {
      annot_df <- as.data.frame(annot_data)
    } else {
      stop("Unexpected annotation data structure in GDS file.")
    }

    # Validate requested features exist
    missing_features <- setdiff(features, names(annot_df))
    if (length(missing_features) > 0) {
      stop("Missing annotation columns in GDS: ",
           paste(missing_features, collapse = ", "))
    }
  }

  # Extract requested features as matrix
  annot_mat <- as.matrix(annot_df[, features, drop = FALSE])
  rownames(annot_mat) <- varinfo

  # Apply variant filter if indices specified (for proportional sampling)
  if (!is.null(variant_indices)) {
    annot_mat <- annot_mat[variant_indices, , drop = FALSE]
  }

  annot_mat
}


#' Load GDS Files with Proportional Sampling
#'
#' @description
#' Internal function that samples proportionally from each GDS file during
#' loading to limit memory usage. Uses SeqArray's variant filtering to read
#' only sampled variants.
#'
#' @param files Character vector of GDS file paths.
#' @param features Feature columns to extract.
#' @param max_controls Total number of controls to sample across all files.
#'
#' @return Numeric matrix with sampled rows.
#'
#' @keywords internal
#' @noRd
.load_gds_with_proportional_sampling <- function(files, features, max_controls) {
  # Check if SeqArray is available
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("Package 'SeqArray' required for GDS support. ",
         "Install with: BiocManager::install('SeqArray')")
  }

  # Step 1: Count variants in each GDS file (fast - just read variant.id length)
  variant_counts <- vapply(files, function(f) {
    gds <- SeqArray::seqOpen(f, readonly = TRUE)
    on.exit(SeqArray::seqClose(gds))
    length(SeqArray::seqGetData(gds, "variant.id"))
  }, integer(1))

  total_variants <- sum(variant_counts)

  # Step 2: If total is less than max_controls, load all
  if (total_variants <= max_controls) {
    message("Total controls (", total_variants, ") <= max_controls. Loading all.")
    data_list <- lapply(files, function(f) {
      .load_control_gds(f, features)
    })
    return(do.call(rbind, data_list))
  }

  # Step 3: Calculate proportional allocation for each file
  message("Proportionally sampling ", max_controls, " from ", total_variants,
          " total controls...")
  proportions <- variant_counts / total_variants
  allocations <- round(proportions * max_controls)

  # Adjust to ensure exact total (rounding may cause slight differences)
  diff <- max_controls - sum(allocations)
  if (diff != 0) {
    largest_idx <- which.max(allocations)
    allocations[largest_idx] <- allocations[largest_idx] + diff
  }

  # Step 4: Load sampled variants from each file using SeqArray filtering
  data_list <- mapply(function(f, n_variants, n_sample) {
    if (n_sample <= 0) {
      return(NULL)
    }
    if (n_sample >= n_variants) {
      # Load all variants from this file
      return(.load_control_gds(f, features))
    }

    # Sample variant indices (1-based)
    sample_indices <- sort(sample(n_variants, n_sample, replace = FALSE))

    # Load only sampled variants using SeqArray filter
    .load_control_gds(f, features, variant_indices = sample_indices)
  }, files, variant_counts, allocations, SIMPLIFY = FALSE)

  # Remove NULLs and combine
  data_list <- data_list[!vapply(data_list, is.null, logical(1))]
  do.call(rbind, data_list)
}
