########## Data Loader for PI Case Variant Preparation ##########
#
# This file contains functions to prepare trait-associated variant lists
# (case variants) for PI estimation training in GLOWr. These functions
# handle data loading, standardization, quality control, and VarInfo
# generation for FAVOR annotation.

#################### EXPORTED MAIN FUNCTIONS ####################

#' Prepare Case Variant List for PI Estimation Training
#'
#' @description
#' Standardizes and validates trait-associated variant lists for use in PI
#' (variant-importance score) estimation training. This function loads curated
#' lists of known disease/trait-associated variants from literature or GWAS
#' Catalog, applies quality control filters, and generates FAVOR-compatible
#' VarInfo identifiers for annotation retrieval.
#'
#' @param data Either a file path (CSV/TSV/Excel) or a data.frame containing
#'   trait-associated variant information. Supported formats: .csv, .tsv, .txt,
#'   .xlsx, .xls
#' @param format Character string specifying input format: "auto" (default),
#'   "gwas_catalog", "csv", "excel", or "dataframe". When "auto", auto-detects
#'   from file extension
#' @param column_mapping Named list for custom column name mapping. Specify
#'   mappings like \code{list(N = "Sample Size", BETA = "Beta_numeric")}
#' @param filter_include Named list for inclusion filters. Only rows matching
#'   specified values are kept. Example: \code{list(TRAIT = "ALS", CHR = 1:5)}.
#'   Takes precedence over filter_exclude
#' @param filter_exclude Named list for exclusion filters. Rows matching
#'   specified values are removed. Example: \code{list(TRAIT = c("trait1", "trait2"))}.
#'   Ignored for columns in filter_include
#' @param qc_filters List of quality control parameters:
#'   \itemize{
#'     \item remove_duplicates: Remove duplicate rsIDs (default: TRUE)
#'     \item duplicate_priority: Column for duplicate resolution (default: "N")
#'     \item filter_autosomes: Keep only chromosomes 1-22 (default: TRUE)
#'     \item filter_indels: Remove insertion/deletion variants (default: FALSE)
#'     \item remove_na: Remove rows with NA in required columns (default: TRUE)
#'     \item remove_na_cols: Columns requiring non-NA values (default: c("rsID", "CHR", "POS"))
#'   }
#' @param exclude_authors Character vector of author names to exclude (case-insensitive
#'   partial matching). Example: \code{c("Nicolas A", "Smith J")}. Useful for
#'   removing overlapping studies or specific publications
#' @param exclude_chromosomes Character or numeric vector of chromosomes to exclude.
#'   Example: \code{c("X", "Y", "MT")} or \code{c(23, 24, 25)}. Applied before
#'   autosome filtering
#' @param verbose Integer controlling verbosity: 0=silent, 1=messages (default),
#'   2=detailed messages
#'
#' @return An S3 object of class \code{glow_pi_case_data} containing:
#'   \itemize{
#'     \item data: data.frame with standardized columns including:
#'       \itemize{
#'         \item rsID: SNP identifier
#'         \item CHR: Chromosome (numeric 1-22 if filter_autosomes=TRUE)
#'         \item POS: Position
#'         \item REF: Reference allele (may be NA if not available)
#'         \item ALT: Alternate allele (may be NA if not available)
#'         \item VarInfo: FAVOR format "CHR-POS-REF-ALT" or "CHR-POS-NA-NA"
#'         \item MAF, P, N, BETA: If available
#'         \item TRAIT, STUDY, PMID, FIRST_AUTHOR: Metadata columns if present
#'       }
#'     \item metadata: List with processing information:
#'       \itemize{
#'         \item n_input: Number of variants in input data
#'         \item n_output: Number of variants after all filtering
#'         \item n_removed: Number of variants removed
#'         \item format: Detected or specified data format
#'         \item filters_applied: Record of all filters used
#'         \item column_mappings: Which columns were renamed
#'         \item date_prepared: Date of preparation
#'       }
#'   }
#'
#' @details
#' \strong{Processing Pipeline:}
#'
#' The function performs the following operations in sequence:
#' \itemize{
#'   \item \strong{Load Data}: Read from file (Excel/CSV/TSV) or accept data.frame
#'   \item \strong{Standardize Columns}: Map to standard names (rsID, CHR, POS, etc.)
#'   \item \strong{Fill Missing Positions}: Search alternative columns if POS has NAs
#'   \item \strong{Filter by Chromosome}: Remove specified chromosomes (if exclude_chromosomes provided)
#'   \item \strong{Convert CHR & Filter Autosomes}: Numeric CHR, keep 1-22 (if filter_autosomes=TRUE)
#'   \item \strong{Apply General Filters}: Include/exclude filters on any column
#'   \item \strong{Remove Duplicates}: By rsID, keep highest priority (if remove_duplicates=TRUE)
#'   \item \strong{Filter by Author}: Remove specified authors (if exclude_authors provided)
#'   \item \strong{Remove NA}: In required columns (if remove_na=TRUE)
#'   \item \strong{Generate VarInfo}: Create CHR-POS-REF-ALT format for FAVOR
#'   \item \strong{Sort by Position}: Order variants by CHR then POS
#'   \item \strong{Reorder Columns}: VarInfo first, then standard order
#' }
#'
#' \strong{VarInfo Generation:}
#'
#' The VarInfo column is generated in FAVOR-compatible format for annotation:
#' \itemize{
#'   \item Standard: "1-12345-A-G" (CHR-POS-REF-ALT)
#'   \item Missing alleles: "1-12345-NA-NA" (with warning)
#'   \item Multi-allelic: Takes first ALT allele
#' }
#'
#' Many literature-curated lists lack REF/ALT information. The function
#' accommodates this by using "NA" placeholders, which allows FAVOR to
#' attempt annotation by position.
#'
#' \strong{Duplicate Handling:}
#'
#' When \code{remove_duplicates = TRUE}:
#' \itemize{
#'   \item Duplicates identified by rsID
#'   \item Variant with highest value in \code{duplicate_priority} column kept
#'   \item Default: Keep variant with largest sample size (N)
#'   \item Alternative: Set \code{duplicate_priority = "P"} to keep lowest p-value
#' }
#'
#' \strong{Author Filtering:}
#'
#' The \code{exclude_authors} parameter provides flexible author-based filtering:
#' \itemize{
#'   \item Case-insensitive partial matching
#'   \item Handles name variations: "Nicolas A", "Nicolas A.", "Nicolas"
#'   \item Common use: Remove overlapping studies or low-quality publications
#'   \item Applied before other filters for efficiency
#' }
#'
#' \strong{Chromosome Filtering:}
#'
#' Two-stage chromosome filtering:
#' \enumerate{
#'   \item \code{exclude_chromosomes}: Remove specific chromosomes (e.g., "X")
#'   \item \code{filter_autosomes}: If TRUE, keep only CHR 1-22 (default)
#' }
#'
#' This allows flexible control: exclude X but keep Y, or exclude all non-autosomes.
#'
#' \strong{Comparison with prepare_B_training_data():}
#'
#' Key differences:
#' \itemize{
#'   \item PI case data: Curated trait-associated variants (literature, GWAS Catalog)
#'   \item B training data: GWAS summary statistics (all variants, genome-wide)
#'   \item PI requires: VarInfo generation for FAVOR annotation
#'   \item PI allows: Missing REF/ALT (common in literature)
#'   \item PI adds: Author filtering, more flexible column handling
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n x f) where n = number of variants, f = number of filters.
#' Dominant operations: duplicate detection O(n log n), filtering O(n x f).
#'
#' @examples
#' \dontrun{
#' # Example 1: ALS known variants (replicate legacy workflow)
#' als_data <- prepare_PI_case_data(
#'   data = "ALS-known-SNPs-raw.xlsx",
#'   format = "excel",
#'   column_mapping = list(
#'     N = "Sample Size",
#'     BETA = "Beta_numeric"
#'   ),
#'   exclude_authors = "Nicolas A",
#'   exclude_chromosomes = "X",
#'   qc_filters = list(
#'     remove_duplicates = TRUE,
#'     duplicate_priority = "N",
#'     filter_autosomes = FALSE,  # Already filtered by exclude_chromosomes
#'     remove_na = FALSE           # Keep variants with missing info
#'   ),
#'   verbose = 2
#' )
#' # Expected: 297 variants (352 input -> 297 after filtering)
#'
#' # Example 2: GWAS Catalog download
#' catalog_data <- prepare_PI_case_data(
#'   data = "gwas_catalog_v1.0-associations.tsv",
#'   format = "gwas_catalog",
#'   filter_include = list(
#'     `DISEASE/TRAIT` = "Type 2 diabetes"
#'   ),
#'   qc_filters = list(
#'     remove_duplicates = TRUE,
#'     duplicate_priority = "N",
#'     remove_na = TRUE,
#'     remove_na_cols = c("rsID", "CHR", "POS")
#'   ),
#'   verbose = 1
#' )
#'
#' # Example 3: Custom filtering (multiple criteria)
#' custom_data <- prepare_PI_case_data(
#'   data = my_variants_df,
#'   filter_include = list(
#'     TRAIT = c("ALS", "FTD"),  # Include only these traits
#'     CHR = 1:5                  # Only first 5 chromosomes
#'   ),
#'   filter_exclude = list(
#'     PMID = c("12345678", "87654321")  # Exclude specific studies
#'   ),
#'   exclude_authors = c("Smith J", "Jones K"),
#'   verbose = 2
#' )
#'
#' # Access results
#' print(custom_data)
#' summary(custom_data)
#' head(custom_data$data)
#' str(custom_data$metadata)
#' }
#'
#' @seealso
#' \code{\link{get_PI}} for PI estimation using prepared case/control data
#' \code{\link{prepare_B_training_data}} for B estimation data preparation
#'
#' @export
prepare_PI_case_data <- function(
  data,
  format = "auto",
  column_mapping = NULL,
  filter_include = NULL,
  filter_exclude = NULL,
  qc_filters = list(
    remove_duplicates = TRUE,
    duplicate_priority = "N",
    filter_autosomes = TRUE,
    filter_indels = FALSE,
    remove_na = TRUE,
    remove_na_cols = c("rsID", "CHR", "POS")
  ),
  exclude_authors = NULL,
  exclude_chromosomes = NULL,
  verbose = 1
) {

  # ========== Step 1: Load Data ==========

  if (verbose >= 1) {
    message("=== Preparing PI Case Data ===")
    message("Loading data...")
  }

  # Determine if input is file path or data.frame
  if (is.character(data)) {
    # File path provided
    if (!file.exists(data)) {
      stop("File not found: ", data)
    }

    # Detect format from extension if format="auto"
    if (format == "auto") {
      ext <- tolower(tools::file_ext(data))
      if (ext %in% c("csv", "tsv", "txt")) {
        format <- "csv"
      } else if (ext %in% c("xlsx", "xls")) {
        format <- "excel"
      } else {
        format <- "csv"  # Default
        if (verbose >= 1) {
          message("  Unknown file extension, assuming CSV format")
        }
      }
    }

    # Read file based on format (reuse helpers from B loader)
    if (format == "gwas_catalog") {
      case_data <- .read_gwas_catalog(data, verbose = verbose)
    } else if (format == "csv") {
      case_data <- .read_csv_file(data, verbose = verbose)
    } else if (format == "excel") {
      case_data <- .read_excel_file(data, verbose = verbose)
    } else {
      stop("Unsupported format: ", format)
    }

    detected_format <- format

  } else if (is.data.frame(data)) {
    # data.frame provided
    case_data <- data
    detected_format <- "dataframe"
    if (verbose >= 1) {
      message(sprintf("  Loaded data.frame with %d variants and %d columns",
                     nrow(case_data), ncol(case_data)))
    }

  } else {
    stop("data must be either a file path (character) or a data.frame")
  }

  # Store original count
  n_input <- nrow(case_data)
  if (verbose >= 1) {
    message(sprintf("  Input: %d variants", n_input))
  }

  # ========== Step 2: Standardize Column Names ==========

  if (verbose >= 1) {
    message("Standardizing column names...")
  }

  # Apply column standardization with ONLY user-specified mappings
  # PI-specific patterns will be handled as a second pass below
  case_data <- standardize_column_names(
    data = case_data,
    column_mapping = column_mapping,  # Only user mappings (can be NULL)
    format_hint = detected_format,
    verbose = verbose
  )

  # Extract mapping metadata from first pass
  column_mappings <- attr(case_data, "column_mappings")
  mapping_sources <- attr(case_data, "mapping_sources")
  if (is.null(column_mappings)) column_mappings <- list()
  if (is.null(mapping_sources)) mapping_sources <- list()

  # Second pass: Apply PI-specific patterns for columns not yet mapped
  # These should be marked as "auto" not "user"
  pi_mappings <- .get_pi_column_variations()

  for (std_name in names(pi_mappings)) {
    # Skip if already mapped
    if (std_name %in% names(case_data)) next

    # Try to find a matching column from PI patterns
    variations <- pi_mappings[[std_name]]
    for (var in variations) {
      if (var %in% names(case_data) && !(var %in% names(column_mappings))) {
        # Rename column
        names(case_data)[names(case_data) == var] <- std_name
        column_mappings[[std_name]] <- var
        mapping_sources[[std_name]] <- "auto-detected"

        if (verbose >= 2) {
          message(sprintf("    %s <- %s [auto]", std_name, var))
        }
        break
      }
    }
  }

  # Update attributes
  attr(case_data, "column_mappings") <- column_mappings
  attr(case_data, "mapping_sources") <- mapping_sources

  # ========== Step 2b: Fill Missing Positions from Alternative Columns ==========
  # If POS has missing values, search for alternatives (pos_GRCh37, pos_GRCh38, etc.)
  case_data <- .fill_missing_positions(case_data, verbose = verbose)

  # ========== Step 3: Filter by Chromosome ==========
  # NOTE: Apply chromosome filtering BEFORE deduplication to match legacy workflow

  if (verbose >= 1) {
    if (!is.null(exclude_chromosomes) && length(exclude_chromosomes) > 0) {
      message(sprintf("Filtering by chromosome (excluding: %s)...",
                      paste(exclude_chromosomes, collapse = ", ")))
    } else {
      message("No chromosome filtering requested, skipping...")
    }
  }

  if (!is.null(exclude_chromosomes) && length(exclude_chromosomes) > 0) {
    n_before_chr <- nrow(case_data)
    case_data <- .filter_by_chromosome(
      data = case_data,
      exclude_chromosomes = exclude_chromosomes,
      chr_col = "CHR",
      verbose = verbose
    )
    n_removed_chr <- n_before_chr - nrow(case_data)
  } else {
    n_removed_chr <- 0
  }

  # ========== Step 4: Convert CHR to Numeric & Filter to Autosomes ==========

  if (verbose >= 1) {
    if ("CHR" %in% names(case_data) && isTRUE(qc_filters$filter_autosomes)) {
      message("Converting CHR to numeric and filtering to autosomes (1-22)...")
    } else {
      message("Skipping autosome filtering (filter_autosomes = FALSE)...")
    }
  }

  if ("CHR" %in% names(case_data) && isTRUE(qc_filters$filter_autosomes)) {
    n_before_autosome <- nrow(case_data)
    case_data <- .ensure_chr_numeric(case_data, verbose = verbose)
    n_removed_autosome <- n_before_autosome - nrow(case_data)
  } else {
    n_removed_autosome <- 0
  }

  # ========== Step 5: Apply General Filters ==========

  if (verbose >= 1) {
    if (!is.null(filter_include) || !is.null(filter_exclude)) {
      message("Applying custom data filters...")
    } else {
      message("No custom filtering requested, skipping...")
    }
  }

  if (!is.null(filter_include) || !is.null(filter_exclude)) {
    n_before_general <- nrow(case_data)
    case_data <- .apply_general_filters(
      data = case_data,
      include = filter_include,
      exclude = filter_exclude,
      pattern_matching = TRUE,
      verbose = verbose
    )
    n_removed_general <- n_before_general - nrow(case_data)
  } else {
    n_removed_general <- 0
  }

  # ========== Step 6: Remove Duplicates ==========
  # NOTE: Apply deduplication BEFORE author filtering to match legacy workflow

  if (verbose >= 1) {
    if (isTRUE(qc_filters$remove_duplicates)) {
      message("Removing duplicate variants...")
    } else {
      message("Duplicate removal disabled, skipping...")
    }
  }

  if (isTRUE(qc_filters$remove_duplicates) && "rsID" %in% names(case_data)) {
    n_before_dup <- nrow(case_data)

    # Sort by rsID and priority column (descending)
    priority_col <- qc_filters$duplicate_priority
    if (!is.null(priority_col) && priority_col %in% names(case_data)) {
      case_data <- case_data[order(case_data$rsID, -case_data[[priority_col]], na.last = TRUE), ]
      if (verbose >= 2) {
        message(sprintf("  Sorting by duplicate_priority column '%s' (descending)", priority_col))
      }
    } else {
      case_data <- case_data[order(case_data$rsID), ]
    }

    # Remove duplicates (keep first = highest priority)
    case_data <- case_data[!duplicated(case_data$rsID), ]

    n_removed_dup <- n_before_dup - nrow(case_data)
    if (verbose >= 1 && n_removed_dup > 0) {
      message(sprintf("  Removed %d duplicate rsIDs (kept highest %s)",
                     n_removed_dup, priority_col))
    }
  } else {
    n_removed_dup <- 0
  }

  # ========== Step 7: Filter by Author ==========
  # NOTE: Apply author filtering AFTER deduplication to match legacy workflow

  if (verbose >= 1) {
    if (!is.null(exclude_authors) && length(exclude_authors) > 0) {
      message(sprintf("Filtering by author (excluding: %s)...",
                      paste(exclude_authors, collapse = ", ")))
    } else {
      message("No author filtering requested, skipping...")
    }
  }

  if (!is.null(exclude_authors) && length(exclude_authors) > 0) {
    n_before_author <- nrow(case_data)
    case_data <- .filter_by_author(
      data = case_data,
      exclude_authors = exclude_authors,
      author_col = "FIRST_AUTHOR",
      verbose = verbose
    )
    n_removed_author <- n_before_author - nrow(case_data)
  } else {
    n_removed_author <- 0
  }

  # ========== Step 8: Remove NA in Required Columns ==========

  if (verbose >= 1) {
    if (isTRUE(qc_filters$remove_na)) {
      message("Removing variants with NA in required columns...")
    } else {
      message("NA removal disabled, skipping...")
    }
  }

  if (isTRUE(qc_filters$remove_na)) {
    n_before_na <- nrow(case_data)

    # Get required columns from qc_filters
    remove_na_cols <- qc_filters$remove_na_cols
    if (is.null(remove_na_cols)) {
      remove_na_cols <- c("rsID", "CHR", "POS")  # Default
    }

    # Only check columns that exist
    cols_to_check <- intersect(remove_na_cols, names(case_data))

    if (length(cols_to_check) > 0) {
      # Remove rows with NA in required columns
      complete_idx <- complete.cases(case_data[, cols_to_check, drop = FALSE])
      case_data <- case_data[complete_idx, ]

      n_removed_na <- n_before_na - nrow(case_data)
      if (verbose >= 1 && n_removed_na > 0) {
        message(sprintf("  Removed %d variants with NA in required columns (%s)",
                       n_removed_na, paste(cols_to_check, collapse = ", ")))
      }
    }
  } else {
    n_removed_na <- 0
  }

  # ========== Step 9: Generate VarInfo Column ==========

  if (verbose >= 1) {
    message("Generating VarInfo column (CHR-POS-REF-ALT format)...")
  }

  # Extract columns for VarInfo generation
  chr_vec <- if ("CHR" %in% names(case_data)) case_data$CHR else NULL
  pos_vec <- if ("POS" %in% names(case_data)) case_data$POS else NULL
  ref_vec <- if ("REF" %in% names(case_data)) case_data$REF else NULL
  alt_vec <- if ("ALT" %in% names(case_data)) case_data$ALT else NULL

  # Check required columns
  if (is.null(chr_vec) || is.null(pos_vec)) {
    stop("CHR and POS columns are required for VarInfo generation")
  }

  # Generate VarInfo
  case_data$VarInfo <- .generate_varinfo_column(
    chr = chr_vec,
    pos = pos_vec,
    ref = ref_vec,
    alt = alt_vec,
    verbose = verbose
  )

  # Count variants with incomplete VarInfo (potential FAVOR annotation issues)
  n_na_pos <- sum(grepl("-NA-", case_data$VarInfo) & !grepl("-NA-NA$", case_data$VarInfo))
  n_na_alleles <- sum(grepl("-NA-NA$", case_data$VarInfo))

  if (verbose >= 1) {
    if (n_na_alleles > 0 || n_na_pos > 0) {
      message(sprintf("  Generated VarInfo for %d variants", nrow(case_data)))
      if (n_na_alleles > 0) {
        message(sprintf("    - %d variants with NA alleles (missing REF/ALT)", n_na_alleles))
      }
      if (n_na_pos > 0) {
        message(sprintf("    - %d variants with NA positions", n_na_pos))
      }
    } else {
      message(sprintf("  Generated VarInfo for %d variants (all complete)", nrow(case_data)))
    }
  }

  # Warn specifically about FAVOR annotation implications
  if (n_na_alleles > 0 || n_na_pos > 0) {
    warning(sprintf(
      paste0("VarInfo column contains %d entries with missing data (NA values).\n",
             "FAVOR annotation may fail or return incomplete results for these variants.\n",
             "Consider:\n",
             "  1. Providing REF/ALT columns in input data, OR\n",
             "  2. Using rsID-to-position lookup to fill missing alleles, OR\n",
             "  3. Filtering out variants with incomplete VarInfo before FAVOR annotation."),
      n_na_alleles + n_na_pos
    ), call. = FALSE)
  }

  # ========== Step 10: Sort by CHR and POS ==========

  if (verbose >= 1) {
    message("Sorting variants by chromosome and position...")
  }

  # Sort by CHR (numeric) then POS (numeric)
  if ("CHR" %in% names(case_data) && "POS" %in% names(case_data)) {
    # Convert CHR to numeric for proper sorting (handles "1", "2", ..., "22")
    chr_numeric <- suppressWarnings(as.numeric(case_data$CHR))
    pos_numeric <- suppressWarnings(as.numeric(case_data$POS))

    sort_order <- order(chr_numeric, pos_numeric, na.last = TRUE)
    case_data <- case_data[sort_order, , drop = FALSE]
    rownames(case_data) <- NULL  # Reset row names
  }

  # ========== Step 11: Reorder Columns (VarInfo First) ==========

  if (verbose >= 1) {
    message("Reordering columns (VarInfo first)...")
  }

  case_data <- .reorder_columns_pi(case_data)

  if (verbose >= 2) {
    message(sprintf("  Column order: %s",
                   paste(names(case_data)[1:min(8, ncol(case_data))], collapse = ", ")))
    if (ncol(case_data) > 8) {
      message(sprintf("    ... and %d more columns", ncol(case_data) - 8))
    }
  }

  # ========== Prepare Output Object ==========

  n_output <- nrow(case_data)
  n_removed <- n_input - n_output

  # Create metadata
  metadata <- list(
    n_input = n_input,
    n_output = n_output,
    n_removed = n_removed,
    format = detected_format,
    filters_applied = list(
      exclude_authors = exclude_authors,
      n_removed_author = n_removed_author,
      exclude_chromosomes = exclude_chromosomes,
      n_removed_chr = n_removed_chr,
      filter_autosomes = qc_filters$filter_autosomes,
      n_removed_autosome = n_removed_autosome,
      filter_include = filter_include,
      filter_exclude = filter_exclude,
      n_removed_general = n_removed_general,
      remove_duplicates = qc_filters$remove_duplicates,
      duplicate_priority = qc_filters$duplicate_priority,
      n_removed_dup = n_removed_dup,
      remove_na = qc_filters$remove_na,
      remove_na_cols = qc_filters$remove_na_cols,
      n_removed_na = n_removed_na
    ),
    column_mappings = column_mappings,
    mapping_sources = mapping_sources,
    date_prepared = Sys.Date(),
    columns_present = names(case_data)
  )

  result <- list(
    data = case_data,
    metadata = metadata
  )

  class(result) <- c("glow_pi_case_data", "list")

  # ========== Summary ==========

  if (verbose >= 1) {
    message(sprintf("\n=== PI Case Data Preparation Complete ==="))
    message(sprintf("  Input:  %d variants", n_input))
    message(sprintf("  Output: %d variants", n_output))
    if (n_removed > 0) {
      message(sprintf("  Removed: %d variants (%.1f%%)",
                     n_removed, 100 * n_removed / n_input))
      if (verbose >= 2) {
        message("  Breakdown:")
        if (n_removed_author > 0) {
          message(sprintf("    - Author filtering: %d", n_removed_author))
        }
        if (n_removed_chr > 0) {
          message(sprintf("    - Chromosome filtering: %d", n_removed_chr))
        }
        if (n_removed_autosome > 0) {
          message(sprintf("    - Autosome filtering: %d", n_removed_autosome))
        }
        if (n_removed_general > 0) {
          message(sprintf("    - Custom filtering: %d", n_removed_general))
        }
        if (n_removed_dup > 0) {
          message(sprintf("    - Duplicate removal: %d", n_removed_dup))
        }
        if (n_removed_na > 0) {
          message(sprintf("    - NA removal: %d", n_removed_na))
        }
      }
    }
    message(sprintf("  Ready for FAVOR annotation"))
  }

  return(result)
}


#################### S3 METHODS ####################

#' Print Method for glow_pi_case_data
#'
#' @description
#' Prints a summary of prepared PI case variant data.
#'
#' @param x Object of class \code{glow_pi_case_data}
#' @param ... Additional arguments (unused)
#'
#' @export
print.glow_pi_case_data <- function(x, ...) {
  cat("GLOW PI Case Data for Variant-Importance Score Estimation\n")
  cat("====================================================\n\n")

  cat("Metadata:\n")
  cat("  Data format:      ", x$metadata$format, "\n")
  cat("  Variants:         ", x$metadata$n_output,
      " (", x$metadata$n_removed, " removed during processing)\n", sep = "")
  cat("  Date prepared:    ", as.character(x$metadata$date_prepared), "\n")

  cat("\nAvailable columns: ", paste(x$metadata$columns_present, collapse = ", "), "\n")

  cat("\nData summary:\n")
  if ("CHR" %in% names(x$data)) {
    # Convert CHR to numeric for range calculation
    chr_numeric <- suppressWarnings(as.numeric(x$data$CHR))
    if (any(!is.na(chr_numeric))) {
      chr_range <- range(chr_numeric, na.rm = TRUE)
      cat("  CHR range:        ", sprintf("%.0f - %.0f", chr_range[1], chr_range[2]), "\n")
    } else {
      cat("  CHR values:       ", "Non-numeric\n")
    }
  }
  if ("VarInfo" %in% names(x$data)) {
    n_na_alleles <- sum(grepl("-NA-NA$", x$data$VarInfo))
    if (n_na_alleles > 0) {
      cat("  VarInfo:          ", "Generated (",
          n_na_alleles, " with NA alleles)\n", sep = "")
    } else {
      cat("  VarInfo:          ", "Generated (all with REF/ALT)\n")
    }
  }

  cat("\nFilters applied:\n")
  filters <- x$metadata$filters_applied
  if (!is.null(filters$exclude_authors) && length(filters$exclude_authors) > 0) {
    cat("  Excluded authors: ", paste(filters$exclude_authors, collapse = ", "),
        " (", filters$n_removed_author, " variants)\n", sep = "")
  }
  if (!is.null(filters$exclude_chromosomes) && length(filters$exclude_chromosomes) > 0) {
    cat("  Excluded chr:     ", paste(filters$exclude_chromosomes, collapse = ", "),
        " (", filters$n_removed_chr, " variants)\n", sep = "")
  }
  if (isTRUE(filters$filter_autosomes)) {
    cat("  Autosomes only:   ", "Yes (CHR 1-22)\n")
  }
  if (isTRUE(filters$remove_duplicates)) {
    cat("  Duplicates:       ", "Removed (", filters$n_removed_dup,
        " by ", filters$duplicate_priority, ")\n", sep = "")
  }

  cat("\nReady for FAVOR annotation and PI estimation\n")

  invisible(x)
}


#' Summary Method for glow_pi_case_data
#'
#' @description
#' Provides a detailed summary of prepared PI case variant data.
#'
#' @param object Object of class \code{glow_pi_case_data}
#' @param ... Additional arguments (unused)
#'
#' @export
summary.glow_pi_case_data <- function(object, ...) {
  cat("GLOW PI Case Data Summary\n")
  cat("=========================\n\n")

  # Print basic info
  print(object)

  # Add detailed filter breakdown if available
  if (!is.null(object$metadata$filters_applied)) {
    cat("\n\nDetailed Filter Breakdown:\n")
    cat("--------------------------\n")
    filters <- object$metadata$filters_applied
    cat(sprintf("  Input variants:           %d\n", object$metadata$n_input))
    if (filters$n_removed_author > 0) {
      cat(sprintf("  After author filter:      %d (-%d)\n",
                 object$metadata$n_input - filters$n_removed_author,
                 filters$n_removed_author))
    }
    if (filters$n_removed_chr > 0) {
      cat(sprintf("  After chr filter:         %d (-%d)\n",
                 object$metadata$n_input - filters$n_removed_author - filters$n_removed_chr,
                 filters$n_removed_chr))
    }
    if (filters$n_removed_autosome > 0) {
      cat(sprintf("  After autosome filter:    %d (-%d)\n",
                 object$metadata$n_input - filters$n_removed_author -
                   filters$n_removed_chr - filters$n_removed_autosome,
                 filters$n_removed_autosome))
    }
    if (filters$n_removed_general > 0) {
      cat(sprintf("  After custom filters:     %d (-%d)\n",
                 object$metadata$n_input - filters$n_removed_author -
                   filters$n_removed_chr - filters$n_removed_autosome -
                   filters$n_removed_general,
                 filters$n_removed_general))
    }
    if (filters$n_removed_dup > 0) {
      cat(sprintf("  After duplicate removal:  %d (-%d)\n",
                 object$metadata$n_output + filters$n_removed_na,
                 filters$n_removed_dup))
    }
    if (filters$n_removed_na > 0) {
      cat(sprintf("  After NA removal:         %d (-%d)\n",
                 object$metadata$n_output,
                 filters$n_removed_na))
    }
    cat(sprintf("  Final output:             %d\n", object$metadata$n_output))
  }

  # Add data structure info
  cat("\n\nData Structure:\n")
  cat("---------------\n")
  str(object$data, max.level = 1, vec.len = 3)

  invisible(object)
}
