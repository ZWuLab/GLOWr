########## Data Loader for PI Control Variant Preparation ##########
#
# This file contains functions to prepare control (neutral) variant lists
# for PI estimation training in GLOWr. These functions handle data loading
# from multiple formats (PLINK, GDS, CSV, data.frame) with GDS as the
# unified intermediate format for genotype data.

#################### EXPORTED MAIN FUNCTIONS ####################

#' Prepare Control Variant List for PI Estimation Training
#'
#' @description
#' Prepares control (neutral/non-causal) variant lists for PI estimation training.
#' Supports multiple input formats with GDS as the unified intermediate format
#' for genotype data. Control variants are used to estimate the null distribution
#' of variant-importance scores in the GLOW framework.
#'
#' @param source Character or data.frame. Input source:
#'   \itemize{
#'     \item PLINK prefix (without .bed/.bim/.fam extension)
#'     \item Path to GDS file (.gds)
#'     \item Path to CSV file (.csv)
#'     \item data.frame with variant information
#'   }
#' @param format Character. Input format: "auto" (default), "plink", "gds",
#'   "csv", "dataframe". When "auto", format is detected from file extension or
#'   object type.
#' @param strategy Character. Control selection strategy (Phase 1: only "provided"):
#'   \itemize{
#'     \item "provided" (default): Use all variants from source as controls
#'   }
#'   Future strategies: "synonymous", "intergenic", "common" (requires annotation)
#' @param n_controls Integer or NULL. If specified, randomly sample to this number
#'   of variants. If NULL (default), use all variants from source.
#' @param gds_output Character or NULL. When format="plink", path to save the
#'   intermediate GDS file. If NULL (default), creates a temporary file that is
#'   deleted after variant extraction. Set this to save GDS for reuse in pipeline
#'   workflows.
#' @param random_seed Integer or NULL. Seed for reproducible random sampling when
#'   \code{n_controls} is specified. Use this for reproducible training sets.
#' @param verbose Integer. Verbosity level: 0=silent, 1=messages (default), 2=detailed.
#'
#' @return An S3 object of class \code{glow_pi_control_data} containing:
#'   \itemize{
#'     \item \strong{data}: data.frame with columns:
#'       \itemize{
#'         \item VarInfo: FAVOR format "CHR-POS-REF-ALT"
#'         \item rsID: Variant identifier (from GDS annotation/id, or NA)
#'         \item CHR: Chromosome (character or numeric)
#'         \item POS: Position (integer)
#'         \item REF: Reference allele
#'         \item ALT: Alternate allele
#'       }
#'     \item \strong{metadata}: List with:
#'       \itemize{
#'         \item n_input: Number of variants in source
#'         \item n_output: Number of variants after sampling
#'         \item format: Detected or specified format
#'         \item gds_path: Path to GDS file (if created/used)
#'         \item gds_persistent: Whether GDS was saved permanently
#'         \item strategy: Control selection strategy used
#'         \item random_seed: Seed used (if sampling performed)
#'         \item date_prepared: Date of preparation
#'       }
#'   }
#'
#' @details
#' \strong{Processing Pipeline (format="plink"):}
#'
#' \enumerate{
#'   \item Convert PLINK to GDS using \code{\link{plink_to_gds}}
#'   \item Extract variant information from GDS (CHR, POS, REF, ALT, rsID)
#'   \item Generate VarInfo column in FAVOR format (CHR-POS-REF-ALT)
#'   \item Random sample if \code{n_controls} specified
#'   \item Delete temp GDS if \code{gds_output} is NULL
#' }
#'
#' \strong{Processing Pipeline (format="gds"):}
#'
#' \enumerate{
#'   \item Open GDS file (read-only)
#'   \item Extract variant information from GDS
#'   \item Generate VarInfo column
#'   \item Random sample if \code{n_controls} specified
#' }
#'
#' \strong{Processing Pipeline (format="csv" or "dataframe"):}
#'
#' \enumerate{
#'   \item Load data from CSV or accept data.frame
#'   \item Validate required columns (CHR, POS, REF, ALT)
#'   \item Generate VarInfo column
#'   \item Random sample if \code{n_controls} specified
#' }
#'
#' \strong{GDS as Canonical Format:}
#'
#' GDS (Genomic Data Structure) serves as the unified intermediate format in
#' GLOWpipeline. This design ensures:
#' \itemize{
#'   \item Consistent variant representation across input formats
#'   \item Efficient extraction of variant information
#'   \item Integration with downstream pipeline steps (aGDS, FAVOR)
#'   \item Optional persistence for pipeline reuse
#' }
#'
#' \strong{Control Variant Selection:}
#'
#' Phase 1 implementation uses "provided" strategy: all variants from the source
#' are treated as controls. Users are responsible for providing appropriate
#' control variants (e.g., common synonymous variants, intergenic variants).
#'
#' Future strategies (requires FAVOR annotation):
#' \itemize{
#'   \item "synonymous": Filter to synonymous variants
#'   \item "intergenic": Filter to intergenic variants
#'   \item "common": Filter to common variants (MAF > threshold)
#' }
#'
#' \strong{Random Sampling:}
#'
#' When \code{n_controls} is specified:
#' \itemize{
#'   \item Random sampling without replacement
#'   \item Set \code{random_seed} for reproducibility
#'   \item Sampling applied after format conversion
#'   \item If \code{n_controls} > available variants, uses all variants with warning
#' }
#'
#' \strong{Computational Complexity:}
#'
#' \itemize{
#'   \item PLINK format: O(n x m) for conversion + O(n) for extraction
#'   \item GDS format: O(n) for extraction
#'   \item CSV format: O(n) for reading
#'   \item Random sampling: O(n)
#' }
#' where n = number of variants, m = number of samples.
#'
#' @examples
#' \dontrun{
#' # Example 1: From PLINK files (GDS saved permanently for reuse)
#' controls <- prepare_PI_control_data(
#'   source = "data/reference_panel",
#'   format = "plink",
#'   gds_output = "data/reference_panel.gds",
#'   verbose = 1
#' )
#'
#' # Example 2: From PLINK files (temp GDS, deleted after use)
#' controls <- prepare_PI_control_data(
#'   source = "data/reference_panel",
#'   format = "plink",
#'   gds_output = NULL,  # Default: temp file
#'   verbose = 1
#' )
#'
#' # Example 3: From existing GDS
#' controls <- prepare_PI_control_data(
#'   source = "data/reference_panel.gds",
#'   format = "gds"
#' )
#'
#' # Example 4: With random sampling (reproducible)
#' controls <- prepare_PI_control_data(
#'   source = "data/reference_panel.gds",
#'   n_controls = 1000,
#'   random_seed = 42,
#'   verbose = 1
#' )
#'
#' # Example 5: From CSV file with variant list
#' controls <- prepare_PI_control_data(
#'   source = "data/control_variants.csv",
#'   format = "csv"
#' )
#'
#' # Example 6: From data.frame
#' control_df <- data.frame(
#'   CHR = c(1, 1, 2),
#'   POS = c(12345, 67890, 11111),
#'   REF = c("A", "G", "C"),
#'   ALT = c("G", "T", "T")
#' )
#' controls <- prepare_PI_control_data(
#'   source = control_df,
#'   format = "dataframe"
#' )
#'
#' # Access results
#' print(controls)
#' summary(controls)
#' head(controls$data)
#' str(controls$metadata)
#' }
#'
#' @seealso
#' \code{\link{plink_to_gds}} for PLINK to GDS conversion
#' \code{\link{get_PI}} for PI estimation using prepared control data
#' \code{\link{prepare_PI_case_data}} for case variant preparation
#'
#' @export
prepare_PI_control_data <- function(
  source,
  format = "auto",
  strategy = "provided",
  n_controls = NULL,
  gds_output = NULL,
  random_seed = NULL,
  verbose = 1
) {

  # ========== Step 1: Validate Strategy ==========

  # Phase 1: Only "provided" strategy supported
  if (!strategy %in% c("provided")) {
    stop("Currently only strategy='provided' is supported. Future: 'synonymous', 'intergenic', 'common'")
  }

  # ========== Step 2: Detect Format ==========

  if (verbose >= 1) {
    message("=== Preparing PI Control Data ===")
  }

  if (format == "auto") {
    detected_format <- .detect_control_format(source)
    if (verbose >= 1) {
      message(sprintf("Auto-detected format: %s", detected_format))
    }
  } else {
    detected_format <- format
  }

  # Validate format
  if (!detected_format %in% c("plink", "gds", "csv", "dataframe")) {
    stop("Unsupported format: ", detected_format)
  }

  # ========== Step 3: Load Variants Based on Format ==========

  gds_path <- NULL
  gds_persistent <- FALSE
  gds_was_temp <- FALSE

  if (detected_format == "plink") {
    # ========== PLINK Format: Convert to GDS, then extract ==========

    if (verbose >= 1) {
      message("Loading variants from PLINK files...")
    }

    # Determine GDS output path (temp or persistent)
    if (is.null(gds_output)) {
      # Create temporary GDS file
      gds_path <- tempfile(fileext = ".gds")
      gds_persistent <- FALSE
      gds_was_temp <- TRUE
      if (verbose >= 2) {
        message(sprintf("  Using temporary GDS: %s", gds_path))
      }
    } else {
      # Use user-specified path (persistent)
      gds_path <- gds_output
      gds_persistent <- TRUE
      if (verbose >= 2) {
        message(sprintf("  Saving GDS to: %s", gds_path))
      }
    }

    # Convert PLINK to GDS
    plink_to_gds(
      plink_prefix = source,
      output_gds = gds_path,
      verbose = verbose
    )

    # Extract variants from GDS
    control_data <- .extract_variants_from_gds(gds_path, verbose = verbose)

    # Delete temp GDS if not persistent
    if (gds_was_temp) {
      if (file.exists(gds_path)) {
        unlink(gds_path)
        if (verbose >= 2) {
          message("  Deleted temporary GDS file")
        }
      }
    }

  } else if (detected_format == "gds") {
    # ========== GDS Format: Extract directly ==========

    if (verbose >= 1) {
      message("Loading variants from GDS file...")
    }

    if (!file.exists(source)) {
      stop("GDS file not found: ", source)
    }

    gds_path <- source
    gds_persistent <- TRUE  # User-provided GDS is persistent

    # Extract variants from GDS
    control_data <- .extract_variants_from_gds(gds_path, verbose = verbose)

  } else if (detected_format == "csv") {
    # ========== CSV Format: Read directly ==========

    if (verbose >= 1) {
      message("Loading variants from CSV file...")
    }

    if (!file.exists(source)) {
      stop("CSV file not found: ", source)
    }

    # Read CSV file
    control_data <- utils::read.csv(source, stringsAsFactors = FALSE)

    # Validate required columns
    required_cols <- c("CHR", "POS", "REF", "ALT")
    missing_cols <- setdiff(required_cols, names(control_data))
    if (length(missing_cols) > 0) {
      stop(sprintf("CSV file missing required columns: %s",
                   paste(missing_cols, collapse = ", ")))
    }

    # Add rsID if not present
    if (!"rsID" %in% names(control_data)) {
      control_data$rsID <- NA_character_
    }

    if (verbose >= 2) {
      message(sprintf("  Loaded %d variants from CSV", nrow(control_data)))
    }

  } else if (detected_format == "dataframe") {
    # ========== Data.frame Format: Use directly ==========

    if (verbose >= 1) {
      message("Loading variants from data.frame...")
    }

    control_data <- source

    # Validate required columns
    required_cols <- c("CHR", "POS", "REF", "ALT")
    missing_cols <- setdiff(required_cols, names(control_data))
    if (length(missing_cols) > 0) {
      stop(sprintf("data.frame missing required columns: %s",
                   paste(missing_cols, collapse = ", ")))
    }

    # Add rsID if not present
    if (!"rsID" %in% names(control_data)) {
      control_data$rsID <- NA_character_
    }

    if (verbose >= 2) {
      message(sprintf("  Loaded %d variants from data.frame", nrow(control_data)))
    }
  }

  # Store original count
  n_input <- nrow(control_data)

  # ========== Step 4: Generate VarInfo Column ==========

  if (verbose >= 1) {
    message("Generating VarInfo column (CHR-POS-REF-ALT format)...")
  }

  # VarInfo already generated by .extract_variants_from_gds() for GDS/PLINK
  # For CSV/dataframe, generate now
  if (!("VarInfo" %in% names(control_data))) {
    control_data$VarInfo <- paste(
      control_data$CHR,
      control_data$POS,
      control_data$REF,
      control_data$ALT,
      sep = "-"
    )
  }

  if (verbose >= 2) {
    message(sprintf("  Generated VarInfo for %d variants", nrow(control_data)))
  }

  # ========== Step 5: Random Sampling (if requested) ==========

  if (!is.null(n_controls)) {
    if (verbose >= 1) {
      message(sprintf("Random sampling to %d variants...", n_controls))
    }

    # Set seed if provided
    if (!is.null(random_seed)) {
      set.seed(random_seed)
      if (verbose >= 2) {
        message(sprintf("  Using random seed: %d", random_seed))
      }
    }

    # Check if requested sample size is larger than available
    if (n_controls > n_input) {
      warning(sprintf(
        "Requested n_controls=%d but only %d variants available. Using all variants.",
        n_controls, n_input
      ), call. = FALSE)
      n_controls <- n_input
    }

    # Random sample
    sample_idx <- sample(seq_len(n_input), size = n_controls, replace = FALSE)
    control_data <- control_data[sample_idx, , drop = FALSE]

    if (verbose >= 1) {
      message(sprintf("  Sampled %d variants from %d", n_controls, n_input))
    }
  }

  n_output <- nrow(control_data)

  # ========== Step 6: Reorder Columns (VarInfo First) ==========

  if (verbose >= 1) {
    message("Reordering columns (VarInfo first)...")
  }

  # Define priority column order
  priority_cols <- c("VarInfo", "rsID", "CHR", "POS", "REF", "ALT")
  present_priority <- priority_cols[priority_cols %in% names(control_data)]
  remaining_cols <- setdiff(names(control_data), priority_cols)
  remaining_cols <- sort(remaining_cols)  # Alphabetically

  # Reorder columns
  final_order <- c(present_priority, remaining_cols)
  control_data <- control_data[, final_order, drop = FALSE]

  # Reset row names
  rownames(control_data) <- NULL

  # ========== Step 7: Build Metadata ==========

  metadata <- list(
    n_input = n_input,
    n_output = n_output,
    format = detected_format,
    gds_path = gds_path,
    gds_persistent = gds_persistent,
    strategy = strategy,
    random_seed = random_seed,
    date_prepared = Sys.Date()
  )

  # ========== Step 8: Create S3 Object ==========

  result <- list(
    data = control_data,
    metadata = metadata
  )

  class(result) <- c("glow_pi_control_data", "list")

  # ========== Summary ==========

  if (verbose >= 1) {
    message(sprintf("\n=== PI Control Data Preparation Complete ==="))
    message(sprintf("  Input:    %d variants", n_input))
    message(sprintf("  Output:   %d variants", n_output))
    message(sprintf("  Strategy: %s", strategy))
    if (!is.null(gds_path)) {
      message(sprintf("  GDS:      %s (%s)",
                      gds_path,
                      ifelse(gds_persistent, "persistent", "deleted")))
    }
    message("  Ready for PI estimation training")
  }

  return(result)
}


#################### S3 METHODS ####################

#' Print Method for glow_pi_control_data
#'
#' @description
#' Prints a summary of prepared PI control variant data.
#'
#' @param x Object of class \code{glow_pi_control_data}
#' @param ... Additional arguments (unused)
#'
#' @export
print.glow_pi_control_data <- function(x, ...) {
  cat("GLOW PI Control Data for Variant-Importance Score Estimation\n")
  cat("========================================================\n\n")

  cat("Metadata:\n")
  cat("  Data format:      ", x$metadata$format, "\n")
  cat("  Variants:         ", x$metadata$n_output, "\n")
  cat("  Strategy:         ", x$metadata$strategy, "\n")
  if (!is.null(x$metadata$gds_path)) {
    cat("  GDS file:         ", x$metadata$gds_path, "\n")
    cat("  GDS persistent:   ", x$metadata$gds_persistent, "\n")
  }
  if (!is.null(x$metadata$random_seed)) {
    cat("  Random seed:      ", x$metadata$random_seed, "\n")
  }
  cat("  Date prepared:    ", as.character(x$metadata$date_prepared), "\n")

  cat("\nFirst 5 variants:\n")
  print(head(x$data, 5))

  cat("\nReady for PI estimation training\n")

  invisible(x)
}


#' Summary Method for glow_pi_control_data
#'
#' @description
#' Provides a detailed summary of prepared PI control variant data.
#'
#' @param object Object of class \code{glow_pi_control_data}
#' @param ... Additional arguments (unused)
#'
#' @export
summary.glow_pi_control_data <- function(object, ...) {
  cat("GLOW PI Control Data Summary\n")
  cat("============================\n\n")

  cat("Source format:      ", object$metadata$format, "\n")
  cat("Input variants:     ", object$metadata$n_input, "\n")
  cat("Output variants:    ", object$metadata$n_output, "\n")
  cat("Strategy:           ", object$metadata$strategy, "\n")
  cat("Prepared:           ", as.character(object$metadata$date_prepared), "\n")

  if (!is.null(object$metadata$gds_path)) {
    cat("\nGDS file:           ", object$metadata$gds_path, "\n")
    cat("GDS persistent:     ", object$metadata$gds_persistent, "\n")
  }

  if (!is.null(object$metadata$random_seed)) {
    cat("\nRandom sampling:\n")
    cat("  Seed:             ", object$metadata$random_seed, "\n")
    cat("  Sample size:      ", object$metadata$n_output, "\n")
    cat("  From total:       ", object$metadata$n_input, "\n")
  }

  cat("\nChromosome distribution:\n")
  print(table(object$data$CHR))

  cat("\nData structure:\n")
  cat("  Columns:          ", paste(names(object$data), collapse = ", "), "\n")
  cat("  Rows:             ", nrow(object$data), "\n")

  invisible(object)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Extract Variant Information from GDS File
#'
#' @description
#' Opens a GDS file and extracts variant information (CHR, POS, REF, ALT, rsID)
#' for control variant preparation. Generates VarInfo column in FAVOR format.
#'
#' @param gds_path Character. Path to GDS file
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with columns: VarInfo, rsID, CHR, POS, REF, ALT
#'
#' @details
#' \strong{GDS Data Extraction:}
#'
#' The function extracts the following data nodes from the GDS file:
#' \itemize{
#'   \item chromosome: Chromosome values (numeric or character)
#'   \item position: Position values (integer)
#'   \item $ref: Reference allele
#'   \item $alt: Alternate allele (takes first if multi-allelic)
#'   \item annotation/id: rsID (or variant.id as fallback)
#' }
#'
#' \strong{rsID Fallback Strategy:}
#'
#' Different GDS files store rsID in different locations:
#' \enumerate{
#'   \item Try "annotation/id" (standard SeqArray location)
#'   \item Try "variant.id" (alternative location)
#'   \item Use NA if neither exists
#' }
#'
#' \strong{Multi-Allelic Variants:}
#'
#' For variants with multiple alternate alleles, only the first allele is used.
#' This matches FAVOR annotation format requirements.
#'
#' @keywords internal
#' @noRd
.extract_variants_from_gds <- function(gds_path, verbose = 1) {

  if (verbose >= 2) {
    message("  Opening GDS file for variant extraction...")
  }

  # Open GDS file (read-only)
  gds <- SeqArray::seqOpen(gds_path, readonly = TRUE)

  # Use tryCatch to ensure GDS is always closed
  variant_data <- tryCatch({

    # Extract chromosome
    chr <- SeqArray::seqGetData(gds, "chromosome")

    # Extract position
    pos <- SeqArray::seqGetData(gds, "position")

    # Extract reference allele
    ref <- SeqArray::seqGetData(gds, "$ref")

    # Extract alternate allele (may be multi-allelic)
    alt_raw <- SeqArray::seqGetData(gds, "$alt")

    # Handle multi-allelic variants (take first allele)
    # alt_raw may be a list or character vector
    if (is.list(alt_raw)) {
      # Extract first allele from each list element
      alt <- sapply(alt_raw, function(x) {
        if (length(x) > 0) x[1] else NA_character_
      })
    } else {
      # Character vector: split by comma and take first
      alt <- sapply(strsplit(as.character(alt_raw), ",", fixed = TRUE),
                    function(x) x[1])
    }

    # Try to get rsID (may be in different locations)
    rsid <- tryCatch(
      SeqArray::seqGetData(gds, "annotation/id"),
      error = function(e) {
        tryCatch(
          as.character(SeqArray::seqGetData(gds, "variant.id")),
          error = function(e2) {
            rep(NA_character_, length(chr))
          }
        )
      }
    )

    # Build data.frame
    df <- data.frame(
      CHR = chr,
      POS = pos,
      REF = ref,
      ALT = alt,
      rsID = rsid,
      stringsAsFactors = FALSE
    )

    # Generate VarInfo
    df$VarInfo <- paste(df$CHR, df$POS, df$REF, df$ALT, sep = "-")

    df

  }, finally = {
    SeqArray::seqClose(gds)
  })

  if (verbose >= 2) {
    message(sprintf("  Extracted %d variants from GDS", nrow(variant_data)))
  }

  return(variant_data)
}


#' Detect Control Data Format from Source
#'
#' @description
#' Auto-detects the format of control data source (PLINK, GDS, CSV, or data.frame).
#'
#' @param source Character or data.frame. Input source
#'
#' @return Character. Detected format: "plink", "gds", "csv", "dataframe"
#'
#' @details
#' \strong{Detection Logic:}
#'
#' \enumerate{
#'   \item If data.frame: return "dataframe"
#'   \item If file ends with .gds: return "gds"
#'   \item If file ends with .csv: return "csv"
#'   \item If PLINK files exist (.bed/.bim/.fam): return "plink"
#'   \item Otherwise: error
#' }
#'
#' \strong{PLINK File Detection:}
#'
#' Checks if all three PLINK files exist:
#' \itemize{
#'   \item {prefix}.bed
#'   \item {prefix}.bim
#'   \item {prefix}.fam
#' }
#'
#' Also handles case where user provides path with extension (e.g., "data.bed").
#'
#' @keywords internal
#' @noRd
.detect_control_format <- function(source) {

  # Check if data.frame
  if (is.data.frame(source)) {
    return("dataframe")
  }

  # Check if character string
  if (!is.character(source) || length(source) != 1) {
    stop("source must be a file path (character) or data.frame")
  }

  # Check file extension
  if (grepl("\\.gds$", source, ignore.case = TRUE)) {
    if (!file.exists(source)) {
      stop("GDS file not found: ", source)
    }
    return("gds")
  }

  if (grepl("\\.csv$", source, ignore.case = TRUE)) {
    if (!file.exists(source)) {
      stop("CSV file not found: ", source)
    }
    return("csv")
  }

  # Check if PLINK files exist
  plink_paths <- .resolve_plink_paths(source)
  plink_files_exist <- all(file.exists(c(
    plink_paths$bed,
    plink_paths$bim,
    plink_paths$fam
  )))

  if (plink_files_exist) {
    return("plink")
  }

  # Cannot detect format
  stop(sprintf(
    "Cannot auto-detect format for: %s\nPlease specify format explicitly ('plink', 'gds', 'csv', 'dataframe')",
    source
  ))
}


#' Resolve PLINK File Paths
#'
#' @description
#' Resolves PLINK file paths from a prefix, handling cases where user provides
#' path with extension.
#'
#' @param plink_prefix Character. Path to PLINK files (with or without extension)
#'
#' @return List with elements: prefix, bed, bim, fam
#'
#' @details
#' \strong{Path Resolution:}
#'
#' If user provides:
#' \itemize{
#'   \item "data/file" -> prefix="data/file", bed="data/file.bed", ...
#'   \item "data/file.bed" -> prefix="data/file", bed="data/file.bed", ...
#'   \item "data/file.bim" -> prefix="data/file", bed="data/file.bed", ...
#' }
#'
#' Strips .bed, .bim, or .fam extension if present.
#'
#' @keywords internal
#' @noRd
.resolve_plink_paths <- function(plink_prefix) {

  # Remove extension if user provided one
  prefix <- sub("\\.(bed|bim|fam)$", "", plink_prefix, ignore.case = TRUE)

  list(
    prefix = prefix,
    bed = paste0(prefix, ".bed"),
    bim = paste0(prefix, ".bim"),
    fam = paste0(prefix, ".fam")
  )
}
