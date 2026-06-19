########## Data Loader for B Estimation ##########
#
# This file contains functions to help users prepare their own GWAS summary
# statistics for B estimation. These functions handle data loading,
# standardization, quality control, and validation.

#################### EXPORTED MAIN FUNCTIONS ####################

#' Prepare User GWAS Summary Statistics for B Estimation (ENHANCED)
#'
#' @description
#' Standardizes and validates user-provided GWAS summary statistics for use
#' in B estimation via \code{\link{get_B}}. This enhanced version includes
#' smart data conversions, study/PMID filtering, position-based duplicate
#' removal, and comprehensive processing to produce publication-ready datasets.
#'
#' @param data Either a file path (CSV/TSV/Excel) or a data.frame containing GWAS
#'   summary statistics. Supported file formats: .csv, .tsv, .txt, .xlsx, .xls
#' @param format Character string specifying the input format: "auto" (default),
#'   "gwas_catalog", "csv", "excel", or "dataframe". When "auto", the function
#'   auto-detects the format from file extension
#' @param trait_type Character string specifying trait type: "binary",
#'   "continuous", or NULL (default). When NULL, the function performs per-row
#'   auto-detection based on OR/BETA columns. When specified, serves as a default
#'   that can be overridden by row-level indicators. A TRAIT_TYPE column will be
#'   added to the output data with per-variant trait type assignments.
#' @param column_mapping Named list for custom column name mapping
#' @param filter_include Named list for inclusion filters. Each element should be
#'   a column name with a vector of values to include. Only rows matching at least
#'   one value per specified column will be kept. Takes precedence over filter_exclude.
#'   Example: \code{list(STUDY = c("Study1", "Study2"), CHR = c(1, 2, 3))}
#' @param filter_exclude Named list for exclusion filters. Each element should be
#'   a column name with a vector of values to exclude. Rows matching any value in
#'   any specified column will be removed. Ignored for columns specified in filter_include.
#'   Example: \code{list(`Reported trait` = c("trait1", "trait2"), CHR = "X")}
#' @param filter_pattern Logical. If TRUE (default), uses pattern matching (grepl)
#'   for character columns. If FALSE, uses exact matching.
#' @param qc_filters List of quality control parameters:
#'   \itemize{
#'     \item min_sample_size: Minimum N (default: 500)
#'     \item max_pvalue: Maximum P-value (default: 1)
#'     \item min_maf: Minimum MAF (default: 0)
#'     \item remove_duplicates: Remove duplicates (default: TRUE)
#'     \item duplicate_priority: Column for duplicate resolution (default: "N")
#'     \item remove_na: Remove rows with NA in required columns (default: TRUE)
#'     \item remove_na_cols: Columns that must have non-NA values when remove_na=TRUE
#'       (default: c("MAF", "P", "N")). Users can add additional required columns
#'       as needed, e.g., c("MAF", "P", "N", "BETA")
#'   }
#' @param convert_or_to_beta Logical. If TRUE (default), converts OR to log(OR)
#'   for binary traits when BETA is not available
#' @param extract_text_values Logical. If TRUE (default), attempts to extract
#'   numeric values from text fields (P-values, sample sizes). \strong{WARNING:}
#'   Text extraction for sample sizes uses heuristics to detect when multiple
#'   groups should be summed (e.g., "men and women", "cases and controls"). While
#'   robust, this is NOT 100\% accurate. Always verify extracted sample sizes match
#'   the original study descriptions, especially for publication-quality analyses.
#' @param filter_autosomes Logical. If TRUE (default), converts CHR to numeric
#'   and filters to keep only autosomal chromosomes (CHR 1-22), removing sex
#'   chromosomes (X, Y), mitochondrial chromosome (MT), and variants with
#'   invalid CHR values. If FALSE, leaves CHR column completely unchanged (no
#'   conversion, no filtering), retaining all variants regardless of chromosome.
#'   Set to FALSE to preserve original CHR values unchanged.
#' @param verbose Integer controlling verbosity: 0=silent, 1=messages (default),
#'   2=detailed messages
#'
#' @return An S3 object of class \code{glow_training_data} containing:
#'   \itemize{
#'     \item data: data.frame with standardized columns in standard order:
#'       rsID, CHR, POS, MAF, P, N, BETA, TRAIT_TYPE, TRAIT, STUDY, PMID, FIRST_AUTHOR, ...
#'     \item metadata: List with comprehensive processing information including
#'       conversions applied, filters used, QC summary, column mappings, and
#'       mapping sources (user-specified vs auto-detected)
#'   }
#'
#' @details
#' \strong{Enhanced Processing Pipeline:}
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Load data from file or data.frame
#'   \item Standardize column names
#'   \item Convert columns to numeric types
#'   \item \strong{Smart Conversions}:
#'     \itemize{
#'       \item Extract sample size (N) from text descriptions
#'       \item Extract P-values from text (handles "5 x 10^-8", etc.)
#'       \item Preserve PVALUE_MLOG as P_mlog10 (not converted to P)
#'       \item Convert OR to log(OR) for binary traits
#'     }
#'   \item MAF harmonization (MAF > 0.5 becomes 1 - MAF)
#'   \item Add TRAIT_TYPE column
#'   \item Filter to autosomes (CHR 1-22) when filter_autosomes=TRUE; leave CHR unchanged when FALSE
#'   \item Filter by studies and/or PMIDs (optional)
#'   \item Apply QC filters (enhanced with position-based duplicate removal)
#'   \item Reorder columns to standard format
#'   \item Validate data
#' }
#'
#' \strong{CRITICAL WARNINGS - Smart Data Processing:}
#'
#' This function includes several "smart" automatic data processing features for
#' convenience. While robust, these are heuristic approaches with limitations.
#' \strong{Users are responsible for verifying that automatic operations are
#' appropriate for their specific data.}
#'
#' \strong{1. Sample Size Extraction (from text):}
#' \itemize{
#'   \item \strong{What it does:} Extracts numeric sample sizes from text descriptions
#'         like "76,067 women, 66,420 men" and intelligently sums or takes max
#'   \item \strong{Limitations:} Pattern-based heuristic; may misinterpret complex text
#'   \item \strong{Failure scenarios:} Non-standard formats, ambiguous descriptions,
#'         multiple cohorts with unclear relationships
#'   \item \strong{User action:} ALWAYS verify extracted N values against source data;
#'         for critical analyses, manually prepare clean N column
#' }
#'
#' \strong{2. P-value Text Extraction:}
#' \itemize{
#'   \item \strong{What it does:} Parses p-values from text like "5 x 10^-8", "< 0.001"
#'   \item \strong{Limitations:} Pattern-based; handles common formats but not all variations
#'   \item \strong{Failure scenarios:} Non-standard notation, unusual characters,
#'         ambiguous inequality symbols
#'   \item \strong{User action:} Check for NA values after extraction; verify converted
#'         p-values are reasonable; for non-standard formats, prepare numeric P column
#' }
#'
#' \strong{3. PVALUE_MLOG Preservation:}
#' \itemize{
#'   \item \strong{What it does:} Preserves PVALUE_MLOG column as P_mlog10 (does NOT convert to P)
#'   \item \strong{Purpose:} Retains -log10(p) values for better numerical stability with extreme p-values
#'   \item \strong{User action:} If you need regular p-values, use get_B() with training_P_mlog10 parameter
#'         which handles conversion internally with proper numerical precision
#' }
#'
#' \strong{4. OR to BETA Conversion:}
#' \itemize{
#'   \item \strong{What it does:} Converts odds ratios to log odds ratios: BETA = log(OR)
#'   \item \strong{Assumptions:} Column truly contains OR (not RR, HR, etc.); correct
#'         allele direction; consistent units
#'   \item \strong{Failure scenarios:} Column contains beta values mislabeled as OR;
#'         percentage increases stored incorrectly; mixed effect measures
#'   \item \strong{User action:} VERIFY OR column contains true odds ratios; inspect
#'         values (typically 0.5-2.0); check effect directions are consistent
#' }
#'
#' \strong{5. CHR Handling (when filter_autosomes = TRUE):}
#' \itemize{
#'   \item \strong{What it does:} Converts chromosome names to numeric and filters to autosomes (1-22)
#'   \item \strong{When filter_autosomes = FALSE:} Leaves CHR column completely unchanged (no conversion, no filtering)
#'   \item \strong{Limitations:} Non-standard naming may not be recognized; assumes
#'         human chromosomes; sex/MT chromosomes removed when filter_autosomes=TRUE
#'   \item \strong{User action:} Set filter_autosomes=FALSE to retain all chromosomes unchanged;
#'         for non-human data, manually prepare CHR column before processing
#' }
#'
#' \strong{6. MAF Harmonization:}
#' \itemize{
#'   \item \strong{What it does:} Automatically converts MAF > 0.5 to 1 - MAF
#'   \item \strong{Assumptions:} Column represents minor allele frequency (not major
#'         allele frequency); MAF definition is standard
#'   \item \strong{Failure scenarios:} Column is actually effect allele frequency (not
#'         necessarily minor); allele frequency for specific allele regardless of
#'         whether it's minor
#'   \item \strong{User action:} VERIFY that automatic conversion is appropriate;
#'         check if your data uses MAF or EAF convention; ensure consistency
#' }
#'
#' \strong{General Recommendations for Critical Analyses:}
#' \enumerate{
#'   \item \strong{Exploratory work:} Use automatic features with verbose=2 to see warnings
#'   \item \strong{Publication-quality:} Manually prepare clean columns (N, P, BETA, CHR, MAF)
#'   \item \strong{Non-standard data:} Pre-process data manually before calling this function
#'   \item \strong{Always:} Inspect result$metadata to see what conversions were applied
#'   \item \strong{Always:} Check result$data to verify conversions are sensible
#'   \item \strong{Document:} Record which automatic features were used and verified
#' }
#'
#' \strong{CRITICAL WARNING - Sample Size Extraction:}
#'
#' When \code{extract_text_values = TRUE} (default), the function attempts to
#' extract sample sizes from text descriptions. This uses intelligent heuristics
#' to detect when multiple groups should be summed:
#'
#' \itemize{
#'   \item \strong{Summed:} "76,067 women, 66,420 men" -> 142,487
#'   \item \strong{Summed:} "5,000 cases + 5,000 controls" -> 10,000
#'   \item \strong{Single:} "394,929 individuals" -> 394,929
#' }
#'
#' \strong{USER RESPONSIBILITY:}
#' \itemize{
#'   \item This extraction is heuristic and may not be 100\% accurate
#'   \item ALWAYS verify extracted sample sizes against original study descriptions
#'   \item Check warnings produced during extraction
#'   \item For critical analyses, manually specify sample sizes in a clean N column
#'   \item The function will warn you when it makes decisions about summing
#' }
#'
#' To inspect extracted values, check \code{result$data$N} after processing and
#' compare with original text in your input data.
#'
#' \strong{TRAIT_TYPE Column (Sophisticated Per-Row Detection):}
#'
#' The output always includes a TRAIT_TYPE column with per-variant trait type
#' assignment using multiple evidence sources. The function uses a conservative
#' approach, defaulting to "unknown" when confidence is low.
#'
#' \strong{Detection Signals (in decreasing confidence order):}
#' \enumerate{
#'   \item \strong{OR column with value} (excluding ambiguous names like "OR.or.BETA") -> binary
#'   \item \strong{Case Size + Control Size both present} -> binary
#'   \item \strong{Disease name in trait column} (e.g., "sclerosis", "diabetes") -> binary
#'   \item \strong{Quantitative trait name} (e.g., "BMI", "bone mineral density") -> continuous
#'   \item \strong{"Susceptibility" keyword in annotations} -> binary
#'   \item \strong{BETA present with no binary indicators} -> continuous (low confidence)
#' }
#'
#' \strong{Important Notes:}
#' \itemize{
#'   \item \strong{This detection is robust but NOT error-proof.} Manual verification
#'         is recommended for critical analyses.
#'   \item The function prioritizes avoiding false classifications - when uncertain,
#'         it assigns "unknown" rather than guessing.
#'   \item Trait type is less critical when using P-value-based (vs. BETA-based)
#'         methods for B estimation.
#'   \item Users can override auto-detection by explicitly specifying \code{trait_type}.
#' }
#'
#' This enables correct handling of mixed-trait datasets (e.g., GWAS Catalog)
#' where different variants come from different trait types.
#'
#' \strong{Duplicate Removal:}
#'
#' Duplicates are removed by both rsID and genomic position (CHR:POS). When
#' duplicates exist, the variant with the largest value in the priority column
#' (default: N = sample size) is kept.
#'
#' \strong{Flexible Filtering System:}
#'
#' The filtering system allows you to filter on ANY column in your data:
#' \itemize{
#'   \item \code{filter_include}: Keep only rows matching specified values
#'   \item \code{filter_exclude}: Remove rows matching specified values
#'   \item \code{filter_pattern}: Use pattern matching (TRUE) or exact matching (FALSE)
#'   \item \code{qc_filters$remove_na_cols}: Specify columns that must have non-NA values
#' }
#'
#' NA filtering is controlled by \code{qc_filters$remove_na} and \code{qc_filters$remove_na_cols}.
#' By default, only MAF, P, and N are required to be non-NA. You can add additional
#' required columns (e.g., BETA) by specifying them in \code{qc_filters$remove_na_cols}.
#'
#' This eliminates the need for pre-filtering and works with any data structure.
#' See Examples for common use cases (ALS, BMD, excluding overlapping studies).
#'
#' @examples
#' \dontrun{
#' # Example 1: Standard analysis (default - autosomes only)
#' result <- prepare_B_training_data(
#'   data = "gwas.csv",
#'   filter_autosomes = TRUE  # Default
#' )
#'
#' # Example 2: Keep all chromosomes (legacy behavior)
#' result <- prepare_B_training_data(
#'   data = "gwas.csv",
#'   filter_autosomes = FALSE  # Retain X, Y, MT, invalid CHR
#' )
#'
#' # Check which columns were renamed
#' print(result$metadata$column_mappings)
#' print(result$metadata$mapping_sources)
#'
#' # Example 3: ALS filtering (exclude specific traits and annotations)
#' als_data <- prepare_B_training_data(
#'   data = "ALS-known-SNPs.xlsx",
#'   trait_type = "binary",
#'   column_mapping = list(BETA = "Beta_numeric"),
#'   filter_exclude = list(
#'     `Reported trait` = c(
#'       "Amyotrophic lateral sclerosis (C9orf72 mutation interaction)",
#'       "Amyotrophic lateral sclerosis (age of onset)",
#'       "Survival in sporadic amyotrophic lateral sclerosis"
#'     ),
#'     `P-value annotation` = c("(age of onset)", "(survival)", "(site of onset)"),
#'     chr = "X"
#'   ),
#'   qc_filters = list(
#'     remove_na = TRUE,
#'     remove_na_cols = c("MAF", "P", "N", "BETA")  # BETA also required for ALS
#'   ),
#'   verbose = 2
#' )
#'
#' # Example 4: BMD filtering (include specific studies only)
#' bmd_data <- prepare_B_training_data(
#'   data = "bmd_gwas.xlsx",
#'   trait_type = "continuous",
#'   filter_include = list(
#'     STUDY = c("Atlas", "153 new loci", "613 new loci")
#'   ),
#'   filter_autosomes = FALSE,  # Keep all chromosomes (legacy behavior)
#'   qc_filters = list(
#'     remove_na = TRUE,
#'     remove_na_cols = c("MAF", "P", "N", "BETA")  # BETA also required for BMD
#'   ),
#'   verbose = 2
#' )
#'
#' # Example 5: Flexible filtering on any column
#' gwas_data <- prepare_B_training_data(
#'   data = "my_gwas.csv",
#'   trait_type = "binary",
#'   filter_exclude = list(
#'     PMID = c("12345678", "87654321"),  # Exclude overlapping studies
#'     Population = "European"             # Exclude specific population
#'   ),
#'   filter_include = list(
#'     CHR = 1:22  # Autosomes only
#'   ),
#'   filter_pattern = TRUE,  # Use pattern matching
#'   verbose = 2
#' )
#' }
#'
#' @seealso
#' \code{\link{get_B}} for B estimation using prepared data
#'
#' @export
prepare_B_training_data <- function(
  data,
  format = "auto",
  trait_type = NULL,
  column_mapping = NULL,
  filter_include = NULL,
  filter_exclude = NULL,
  filter_pattern = TRUE,
  qc_filters = list(
    min_sample_size = 500,
    max_pvalue = 1,
    min_maf = 0,
    remove_duplicates = TRUE,
    duplicate_priority = "N",
    remove_na = TRUE,
    remove_na_cols = c("MAF", "P", "N")
  ),
  convert_or_to_beta = TRUE,
  extract_text_values = TRUE,
  filter_autosomes = TRUE,
  verbose = 1
) {

  # Track conversions for metadata
  conversions_log <- list(
    OR_to_BETA = FALSE,
    text_to_N = FALSE,
    text_to_P = FALSE,
    CHR_to_numeric = FALSE
  )

  # ========== Step 1: Load Data ==========

  if (verbose >= 1) {
    message("Step 1: Loading data...")
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
        format <- "csv"  # Default to CSV
        if (verbose >= 1) {
          message("  Unknown file extension, assuming CSV format")
        }
      }
    }

    # Read file based on format
    if (format == "gwas_catalog") {
      gwas_data <- .read_gwas_catalog(data, verbose = verbose)
    } else if (format == "csv") {
      gwas_data <- .read_csv_file(data, verbose = verbose)
    } else if (format == "excel") {
      gwas_data <- .read_excel_file(data, verbose = verbose)
    } else {
      stop("Unsupported format: ", format, ". Use 'csv', 'gwas_catalog', 'excel', or 'auto'")
    }

    detected_format <- format

  } else if (is.data.frame(data)) {
    # data.frame provided
    gwas_data <- data
    detected_format <- "dataframe"
    if (verbose >= 1) {
      message(sprintf("  Loaded data.frame with %d variants and %d columns",
                     nrow(gwas_data), ncol(gwas_data)))
    }

  } else {
    stop("data must be either a file path (character) or a data.frame")
  }

  # Store original count
  n_original <- nrow(gwas_data)

  # ========== Step 2: Standardize Column Names ==========

  if (verbose >= 1) {
    message("Step 2: Standardizing column names...")
  }

  filters_log <- list(
    filter_include = filter_include,
    filter_exclude = filter_exclude,
    n_removed = 0
  )

  gwas_data <- standardize_column_names(
    data = gwas_data,
    column_mapping = column_mapping,
    format_hint = detected_format,
    verbose = verbose
  )

  # NEW: Extract mapping metadata from attributes
  column_mappings <- attr(gwas_data, "column_mappings")
  mapping_sources <- attr(gwas_data, "mapping_sources")

  # ========== Step 3: Smart Text Extraction (BEFORE numeric conversion) ==========

  if (verbose >= 1) {
    message("Step 3: Extracting values from text columns...")
  }

  if (extract_text_values) {

    # 3a. Extract sample size from N column if it contains text
    # After standardize_column_names(), the column is already named "N"
    if ("N" %in% names(gwas_data)) {
      # Check if N contains text (not purely numeric)
      if (is.character(gwas_data$N) && any(grepl("[a-zA-Z]", gwas_data$N, perl = TRUE))) {
        if (verbose >= 1) {
          message("  Extracting sample sizes from text in N column...")
        }

        gwas_data$N <- extract_sample_size_from_text(gwas_data$N)
        conversions_log$text_to_N <- TRUE

        if (verbose >= 2) {
          n_extracted <- sum(!is.na(gwas_data$N))
          message(sprintf("    Successfully extracted N for %d/%d variants",
                         n_extracted, nrow(gwas_data)))
        }
      }
    }

    # 3b. Extract P from text if needed (BEFORE numeric conversion)
    p_cols <- c("P", "P.VALUE", "P-VALUE", "PVALUE", "P_VALUE")
    p_col <- intersect(p_cols, names(gwas_data))

    if (length(p_col) > 0) {
      for (col in p_col) {
        # Check if P has special text format
        if (is.character(gwas_data[[col]]) && any(grepl("[xxE^]", gwas_data[[col]], ignore.case = TRUE, perl = TRUE))) {
          if (verbose >= 2) {
            message(sprintf("  Extracting P-values from text format in column: %s", col))
          }

          gwas_data[[col]] <- .extract_pvalue_from_text(gwas_data[[col]], verbose = (verbose >= 1))
          conversions_log$text_to_P <- TRUE
          break
        }
      }
    }
  }

  # ========== Step 4: Convert Remaining Columns to Numeric ==========

  if (verbose >= 1) {
    message("Step 4: Converting columns to numeric types...")
  }

  # Now convert columns to numeric (after text extraction)
  # NOTE: CHR is NOT included here - it's handled separately in Step 9 by .ensure_chr_numeric()
  # which does BOTH conversion AND autosomal filtering (removes X, Y, MT chromosomes)
  numeric_cols <- c("MAF", "P", "BETA", "N", "POS", "OR", "P_mlog10")

  # Track conversion status for reporting
  already_numeric <- character()
  converted_to_numeric <- character()
  converted_with_na <- character()
  not_present <- character()

  for (col in numeric_cols) {
    if (col %in% names(gwas_data)) {
      if (!is.numeric(gwas_data[[col]])) {
        # Try to convert to numeric
        original_values <- gwas_data[[col]]
        gwas_data[[col]] <- suppressWarnings(as.numeric(gwas_data[[col]]))

        # Count how many became NA due to conversion
        n_na_after <- sum(is.na(gwas_data[[col]]))
        n_na_before <- sum(is.na(original_values))
        n_converted_to_na <- n_na_after - n_na_before

        # Track conversion status
        if (n_converted_to_na > 0) {
          if (verbose >= 2) {
            converted_with_na <- c(converted_with_na,
                                  sprintf("%s (%d non-numeric values became NA)", col, n_converted_to_na))
          } else {
            converted_with_na <- c(converted_with_na, col)
          }
        } else {
          converted_to_numeric <- c(converted_to_numeric, col)
        }
      } else {
        # Already numeric
        already_numeric <- c(already_numeric, col)
      }
    } else {
      # Column not present
      not_present <- c(not_present, col)
    }
  }

  # Report conversion summary
  if (verbose >= 1) {
    if (length(already_numeric) > 0) {
      message(sprintf("  Already numeric: %s", paste(already_numeric, collapse = ", ")))
    }
    if (length(converted_to_numeric) > 0) {
      message(sprintf("  Converted to numeric: %s", paste(converted_to_numeric, collapse = ", ")))
    }
    if (length(converted_with_na) > 0) {
      if (verbose >= 2) {
        # Show detailed info with NA counts
        message("  Converted to numeric (with NA conversion):")
        for (item in converted_with_na) {
          message(sprintf("    - %s", item))
        }
      } else {
        # Just list columns
        message(sprintf("  Converted to numeric (with NA conversion): %s",
                       paste(converted_with_na, collapse = ", ")))
      }
    }
    if (verbose >= 2 && length(not_present) > 0) {
      message(sprintf("  Not present: %s", paste(not_present, collapse = ", ")))
    }
  }

  # ========== Step 5: Additional Smart Conversions ==========

  if (verbose >= 1) {
    message("Step 5: Checking for OR to BETA conversion...")
  }

  if (convert_or_to_beta && !is.null(trait_type) && trait_type == "binary") {
    if ("OR" %in% names(gwas_data) && (!"BETA" %in% names(gwas_data) || all(is.na(gwas_data$BETA)))) {
      if (verbose >= 1) {
        message("  Converting OR to BETA (log odds ratio) for binary traits")
      }

      gwas_data$BETA <- .convert_or_to_beta(gwas_data$OR, warn = (verbose >= 1))
      conversions_log$OR_to_BETA <- TRUE

      if (verbose >= 2) {
        n_converted <- sum(!is.na(gwas_data$BETA))
        message(sprintf("    Converted %d OR values to BETA", n_converted))
      }
    } else {
      if (verbose >= 1) {
        message("  No OR to BETA conversion needed")
      }
    }
  } else {
    if (verbose >= 1) {
      message("  No OR to BETA conversion needed")
    }
  }

  # ========== Step 6: MAF Harmonization ==========

  if (verbose >= 1) {
    message("Step 6: Harmonizing MAF values...")
  }

  if ("MAF" %in% names(gwas_data)) {
    idx_swap <- which(gwas_data$MAF > 0.5)
    n_swapped <- length(idx_swap)
    if (n_swapped > 0) {
      gwas_data$MAF[idx_swap] <- 1 - gwas_data$MAF[idx_swap]

      if (verbose >= 1) {
        message(sprintf("  Converted %d variants with MAF > 0.5 to 1 - MAF (%.1f%%)",
                       n_swapped, 100 * n_swapped / nrow(gwas_data)))
      }
    }
  }

  # ========== Step 7: Add TRAIT_TYPE Column ==========

  if (verbose >= 1) {
    message("Step 7: Adding TRAIT_TYPE column...")
  }

  gwas_data <- .add_trait_type_column(gwas_data, trait_type = trait_type, verbose = 0)

  # ========== Step 8: Apply General Filters (BEFORE CHR conversion) ==========

  if (verbose >= 1 && (!is.null(filter_include) || !is.null(filter_exclude))) {
    message("Step 8: Applying custom data filters...")
  } else if (verbose >= 1) {
    message("Step 8: No custom filtering requested, skipping...")
  }

  if (!is.null(filter_include) || !is.null(filter_exclude)) {
    n_before_filters <- nrow(gwas_data)
    gwas_data <- .apply_general_filters(
      data = gwas_data,
      include = filter_include,
      exclude = filter_exclude,
      pattern_matching = filter_pattern,
      verbose = verbose
    )
    filters_log$n_removed <- filters_log$n_removed + (n_before_filters - nrow(gwas_data))
  }

  # ========== Step 9: Ensure CHR is Numeric & Optionally Filter to Autosomes ==========

  if (verbose >= 1) {
    if (filter_autosomes) {
      message("Step 9: Ensuring CHR column is numeric and filtering to autosomes...")
    } else {
      message("Step 9: Skipping CHR processing (filter_autosomes = FALSE)...")
    }
  }

  if ("CHR" %in% names(gwas_data) && filter_autosomes) {
    # Filter to autosomes (default behavior - removes X, Y, MT, invalid CHR)
    # This function does TWO things:
    # 1. Convert character CHR ("X", "1", etc.) to numeric
    # 2. Filter to keep ONLY autosomes (CHR 1-22), removing X, Y, MT
    gwas_data <- .ensure_chr_numeric(gwas_data, verbose = verbose)
    conversions_log$CHR_to_numeric <- TRUE
  }
  # When filter_autosomes = FALSE, leave CHR column completely unchanged

  # ========== Step 10: Apply QC Filters ==========

  if (verbose >= 1) {
    message("Step 10: Applying quality control filters...")
  }

  # Merge user filters with defaults
  default_filters <- list(
    min_sample_size = 500,
    max_pvalue = 1,
    min_maf = 0,
    remove_duplicates = TRUE,
    duplicate_priority = "N",
    remove_na = TRUE,
    remove_na_cols = c("MAF", "P", "N")
  )

  for (name in names(default_filters)) {
    if (is.null(qc_filters[[name]])) {
      qc_filters[[name]] <- default_filters[[name]]
    }
  }

  # Apply QC filters using enhanced function (includes P-value cleaning)
  qc_result <- apply_qc_filters(
    data = gwas_data,
    filters = qc_filters,
    verbose = verbose
  )

  gwas_data <- qc_result$data
  qc_summary <- qc_result$qc_summary

  # Enhanced duplicate removal (by position)
  if (qc_filters$remove_duplicates) {
    n_before_dup <- nrow(gwas_data)
    gwas_data <- .remove_duplicates_by_position(
      data = gwas_data,
      priority = qc_filters$duplicate_priority,
      verbose = verbose
    )

    # Update QC summary if duplicates were removed
    if (nrow(gwas_data) < n_before_dup) {
      # The QC summary was already updated by apply_qc_filters for rsID duplicates
      # This captures additional position-based duplicate removal
    }
  }

  # ========== Step 11: Reorder Columns ==========

  if (verbose >= 1) {
    message("Step 11: Reordering columns to standard format...")
  }

  gwas_data <- .reorder_columns_standard(gwas_data)

  if (verbose >= 2) {
    message(sprintf("  Column order: %s",
                   paste(names(gwas_data)[1:min(10, ncol(gwas_data))], collapse = ", ")))
    if (ncol(gwas_data) > 10) {
      message(sprintf("    ... and %d more columns", ncol(gwas_data) - 10))
    }
  }

  # Report TRAIT_TYPE distribution after all filtering
  if (verbose >= 1 && "TRAIT_TYPE" %in% names(gwas_data)) {
    trait_counts <- table(gwas_data$TRAIT_TYPE)
    n_types <- length(trait_counts)

    if (n_types == 1) {
      # Single trait type
      message(sprintf("  TRAIT_TYPE: %s for all %d variants",
                     names(trait_counts)[1], nrow(gwas_data)))
    } else {
      # Mixed trait types
      message(sprintf("  TRAIT_TYPE distribution (mixed dataset):"))
      for (tt in names(trait_counts)) {
        pct <- 100 * trait_counts[tt] / nrow(gwas_data)
        message(sprintf("    - %s: %d variants (%.1f%%)",
                       tt, trait_counts[tt], pct))
      }
    }

    if ("unknown" %in% names(trait_counts) && verbose >= 2) {
      message(sprintf("  Note: %d variants could not be confidently classified (conservative approach)",
                     trait_counts["unknown"]))
    }
  }

  # ========== Final Validation ==========

  if (verbose >= 1) {
    message("\nValidating prepared data...")
  }

  validation_result <- validate_B_training_data(
    data = gwas_data,
    method = "auto",
    verbose = verbose
  )

  if (!validation_result$valid) {
    stop("Data validation failed: ", paste(validation_result$errors, collapse = "; "))
  }

  # ========== Prepare Output Object ==========

  n_final <- nrow(gwas_data)

  # Determine actual trait type from TRAIT_TYPE column if present
  if ("TRAIT_TYPE" %in% names(gwas_data)) {
    trait_types_present <- unique(gwas_data$TRAIT_TYPE)
    if (length(trait_types_present) == 1) {
      final_trait_type <- trait_types_present[1]
    } else {
      final_trait_type <- "mixed"  # Both binary and continuous
    }
  } else {
    final_trait_type <- trait_type
  }

  metadata <- list(
    trait_type = final_trait_type,
    n_variants_original = n_original,
    n_variants_final = n_final,
    n_variants_removed = n_original - n_final,
    format = detected_format,
    conversions_applied = conversions_log,
    filters_applied = filters_log,
    qc_applied = qc_filters,
    qc_summary = qc_summary,
    column_mappings = column_mappings,        # NEW: Store column mapping info
    mapping_sources = mapping_sources,        # NEW: Store mapping source info
    filter_autosomes = filter_autosomes,      # NEW: Record this setting
    date_prepared = Sys.Date(),
    columns_present = names(gwas_data)
  )

  result <- list(
    data = gwas_data,
    metadata = metadata
  )

  class(result) <- c("glow_training_data", "list")

  if (verbose >= 1) {
    message(sprintf("\n=== Data preparation complete ==="))
    message(sprintf("  Final: %d variants ready for B estimation", n_final))
    if (n_original > n_final) {
      message(sprintf("  Removed: %d variants (%.1f%%) during processing",
                     n_original - n_final,
                     100 * (n_original - n_final) / n_original))
    }

    # Report conversions applied
    conversions_applied <- names(conversions_log)[unlist(conversions_log)]
    if (length(conversions_applied) > 0) {
      message(sprintf("  Conversions: %s", paste(conversions_applied, collapse = ", ")))
    }
  }

  return(result)
}



#' Standardize Column Names in GWAS Data
#'
#' @description
#' Maps various column name formats to standard names used by GLOWr functions.
#' This function recognizes common variations in GWAS summary statistics files
#' and converts them to a consistent naming scheme. Returns data with attributes
#' containing mapping metadata to track which columns were renamed and whether
#' mappings were user-specified or auto-detected.
#'
#' @param data data.frame containing GWAS summary statistics
#' @param column_mapping Named list providing custom column mappings. Each
#'   element should be named with a standard column name (e.g., "MAF", "P")
#'   and contain the corresponding column name in your data
#' @param format_hint Character string suggesting the data format ("gwas_catalog",
#'   "csv", etc.). Used to help with ambiguous cases
#' @param verbose Integer verbosity level (0=silent, 1=messages)
#'
#' @return data.frame with standardized column names and two attributes:
#'   \itemize{
#'     \item column_mappings: Named list showing original -> standard mappings applied
#'     \item mapping_sources: Named list indicating "user-specified" or "auto-detected"
#'       for each mapping
#'   }
#'
#' @details
#' Standard column names:
#' \itemize{
#'   \item MAF: Minor allele frequency
#'   \item P: P-value
#'   \item BETA: Effect size
#'   \item N: Sample size
#'   \item rsID: SNP identifier
#'   \item CHR: Chromosome
#'   \item POS: Position
#' }
#'
#' Note: This function only renames columns. It does NOT harmonize MAF values.
#' MAF harmonization (converting MAF > 0.5 to 1 - MAF) is handled separately
#' by \code{\link{prepare_B_training_data}}.
#'
#' \strong{Mapping Source Tracking:}
#'
#' The function tracks whether each column mapping was user-specified (via
#' \code{column_mapping} parameter) or auto-detected. This metadata is stored
#' as attributes on the returned data.frame and can be accessed via:
#' \code{attr(data, "column_mappings")} and \code{attr(data, "mapping_sources")}.
#'
#' \strong{Computational Complexity:}
#'
#' O(n x c x v) where n = rows, c = columns, v = number of variation patterns.
#' In practice, very fast since c and v are small (typically c < 50, v < 10).
#'
#' @examples
#' \dontrun{
#' # Example with custom column names
#' my_data <- data.frame(
#'   snp = paste0("rs", 1:10),
#'   freq = runif(10, 0.1, 0.9),
#'   pval = runif(10, 1e-5, 0.1)
#' )
#'
#' standardized <- standardize_column_names(
#'   data = my_data,
#'   column_mapping = list(
#'     rsID = "snp",
#'     MAF = "freq",
#'     P = "pval"
#'   )
#' )
#'
#' # Check mapping sources
#' print(attr(standardized, "column_mappings"))
#' print(attr(standardized, "mapping_sources"))
#' }
#'
#' @export
standardize_column_names <- function(
  data,
  column_mapping = NULL,
  format_hint = "auto",
  verbose = 1
) {

  # Define standard column names and their common variations
  standard_mappings <- list(
    MAF = c("MAF", "RAF", "FREQ", "EAF", "AF",
            "RISK.ALLELE.FREQUENCY", "RISK_ALLELE_FREQUENCY",
            "EFFECT_ALLELE_FREQUENCY", "ALT_FREQS"),
    P = c("P", "PVALUE", "P_VALUE", "P.VALUE", "P-VALUE", "PVAL"),
    BETA = c("BETA", "EFFECT", "EFFECT_SIZE", "B", "BETA_NUMERIC", "Beta_numeric"),
    OR = c("OR", "ODDS_RATIO", "ODDSRATIO"),  # OR as separate column (not BETA)
    N = c("N", "SAMPLE_SIZE", "SAMPLE.SIZE", "SAMPLE SIZE",
          "INITIAL.SAMPLE.SIZE", "INITIAL SAMPLE SIZE",
          "N_TOTAL", "NEFF", "N_EFF"),
    rsID = c("SNP", "SNPS", "RSID", "RS_ID", "VARIANT_ID", "ID",
             "MARKERNAME", "MARKER"),
    CHR = c("CHR", "CHROM", "CHROMOSOME", "CHR_ID", "#CHR"),
    POS = c("POS", "POSITION", "BP", "CHR_POS", "BASE_PAIR_LOCATION"),
    P_mlog10 = c("P_MLOG10", "PVALUE_MLOG", "PVALUE_MLOG10", "MLOG10P", "MLOG_P", "NEG_LOG10_P")
  )

  # Get current column names (case-insensitive matching)
  current_names <- names(data)
  new_names <- current_names
  mapping_applied <- list()
  mapping_sources <- list()  # NEW: Track source of each mapping

  # Apply user-provided mappings first
  if (!is.null(column_mapping)) {
    for (std_name in names(column_mapping)) {
      user_col <- column_mapping[[std_name]]

      if (user_col %in% current_names) {
        # Find the index of the user's column
        idx <- which(current_names == user_col)
        new_names[idx] <- std_name
        mapping_applied[[std_name]] <- user_col
        mapping_sources[[std_name]] <- "user-specified"  # NEW: Mark as user-specified
      } else {
        warning("Column '", user_col, "' specified in column_mapping not found in data")
      }
    }
  }

  # Apply automatic mappings for remaining columns
  for (std_name in names(standard_mappings)) {
    # Skip if already mapped by user
    if (std_name %in% names(mapping_applied)) {
      next
    }

    variations <- standard_mappings[[std_name]]

    # Try case-insensitive matching
    for (var in variations) {
      idx <- which(toupper(current_names) == toupper(var))

      if (length(idx) > 0) {
        # Found a match
        new_names[idx[1]] <- std_name  # Use first match
        mapping_applied[[std_name]] <- current_names[idx[1]]
        mapping_sources[[std_name]] <- "auto-detected"  # NEW: Mark as auto-detected

        if (length(idx) > 1 && verbose >= 1) {
          message(sprintf("  Warning: Multiple matches for '%s', using '%s'",
                         std_name, current_names[idx[1]]))
        }

        break  # Stop after first variation match
      }
    }
  }

  # Apply new names
  names(data) <- new_names

  # Report mappings with source indicators
  if (verbose >= 1 && length(mapping_applied) > 0) {
    message("  Standardized column names:")
    for (std_name in names(mapping_applied)) {
      # Add source tag to message
      source_tag <- if (mapping_sources[[std_name]] == "user-specified") "[user]" else "[auto]"
      message(sprintf("    %s <- %s %s", std_name, mapping_applied[[std_name]], source_tag))
    }
  }

  # NEW: Attach mapping metadata as attributes
  attr(data, "column_mappings") <- mapping_applied
  attr(data, "mapping_sources") <- mapping_sources

  return(data)
}


#' Apply Quality Control Filters to GWAS Data
#'
#' @description
#' Applies a series of quality control filters to GWAS summary statistics,
#' removing variants that don't meet specified criteria.
#'
#' @param data data.frame with GWAS summary statistics
#' @param filters List of filter parameters with elements:
#'   \itemize{
#'     \item min_sample_size: Minimum N (default: 500)
#'     \item max_pvalue: Maximum p-value (default: 1)
#'     \item min_maf: Minimum MAF (default: 0)
#'     \item remove_duplicates: Remove duplicate rsIDs (default: TRUE)
#'     \item remove_na: Remove rows with NA in required columns (default: TRUE)
#'     \item remove_na_cols: Columns that must have non-NA values when remove_na=TRUE
#'       (default: c("MAF", "P", "N")). Can be customized by users.
#'   }
#'   \strong{Note:} MAF > 0.5 values should be converted to 1 - MAF before calling
#'   this function (done automatically in \code{\link{prepare_B_training_data}}).
#' @param verbose Integer verbosity level (0=silent, 1=messages)
#'
#' @return List with elements:
#'   \itemize{
#'     \item data: Filtered data.frame
#'     \item qc_summary: data.frame summarizing what was filtered
#'   }
#'
#' @details
#' Filters are applied in this order:
#' \enumerate{
#'   \item Clean invalid P-values (P<0 or P>1)
#'   \item Remove rows with NA in required columns
#'   \item Filter by sample size
#'   \item Filter by p-value threshold
#'   \item Filter by MAF range
#'   \item Remove duplicate rsIDs
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n x f) where n = number of variants, f = number of filters (typically f = 6).
#' Current implementation creates multiple data frame copies during filtering, which
#' could be optimized for very large datasets by using a single logical index vector.
#'
#' @examples
#' \dontrun{
#' # Apply strict QC filters with custom required columns
#' qc_result <- apply_qc_filters(
#'   data = my_gwas_data,
#'   filters = list(
#'     min_sample_size = 10000,
#'     max_pvalue = 0.001,
#'     min_maf = 0.01,
#'     remove_na_cols = c("MAF", "P", "N", "BETA")  # Also require BETA
#'   )
#' )
#'
#' filtered_data <- qc_result$data
#' print(qc_result$qc_summary)
#' }
#'
#' @export
apply_qc_filters <- function(
  data,
  filters = list(),
  verbose = 1
) {

  # Set default filter values if not provided
  default_filters <- list(
    min_sample_size = 500,
    max_pvalue = 1,
    min_maf = 0,
    remove_duplicates = TRUE,
    remove_na = TRUE,
    remove_na_cols = c("MAF", "P", "N")
  )

  # Merge user filters with defaults
  for (name in names(default_filters)) {
    if (is.null(filters[[name]])) {
      filters[[name]] <- default_filters[[name]]
    }
  }

  # Initialize QC summary
  qc_summary <- data.frame(
    filter = character(),
    n_before = integer(),
    n_removed = integer(),
    n_after = integer(),
    stringsAsFactors = FALSE
  )

  n_start <- nrow(data)

  # Helper function to record filter step
  .record_filter <- function(filter_name, n_before, n_after) {
    qc_summary <<- rbind(qc_summary, data.frame(
      filter = filter_name,
      n_before = n_before,
      n_removed = n_before - n_after,
      n_after = n_after,
      stringsAsFactors = FALSE
    ))
  }

  # ========== Filter 0: Clean invalid P-values ==========

  if ("P" %in% names(data)) {
    n_before <- nrow(data)

    # Remove invalid P-values: P<0, P>1, or Inf
    invalid_p <- data$P < 0 | data$P > 1 | is.infinite(data$P)
    n_invalid <- sum(invalid_p, na.rm = TRUE)

    if (n_invalid > 0) {
      # Keep only valid P-values (or NA)
      data <- data[!invalid_p | is.na(data$P), ]

      n_after <- nrow(data)
      .record_filter("Clean invalid P-values (P<0 or P>1)", n_before, n_after)

      if (verbose >= 1 && n_before > n_after) {
        message(sprintf("  Removed %d variants with invalid P-values (P<0 or P>1)",
                       n_before - n_after))
      }
    }
  }

  # ========== Filter 1: Remove NA in required columns ==========

  if (filters$remove_na) {
    n_before <- nrow(data)

    # Get required columns from filters (default: MAF, P, N)
    remove_na_cols <- filters$remove_na_cols

    # Only filter on columns that exist
    cols_to_check <- intersect(remove_na_cols, names(data))

    if (length(cols_to_check) > 0) {
      # Create logical vector for complete cases
      complete_idx <- complete.cases(data[, cols_to_check, drop = FALSE])
      data <- data[complete_idx, ]

      n_after <- nrow(data)

      # Create informative filter name showing which columns were checked
      filter_name <- sprintf("Remove NA in required columns (%s)",
                            paste(cols_to_check, collapse = ", "))
      .record_filter(filter_name, n_before, n_after)

      if (verbose >= 1 && n_before > n_after) {
        message(sprintf("  Removed %d variants with NA in required columns (%s)",
                       n_before - n_after, paste(cols_to_check, collapse = ", ")))
      }
    }
  }

  # ========== Filter 2: Sample size filter ==========

  if ("N" %in% names(data) && !is.null(filters$min_sample_size)) {
    n_before <- nrow(data)
    data <- data[!is.na(data$N) & data$N >= filters$min_sample_size, ]
    n_after <- nrow(data)

    .record_filter(sprintf("N >= %d", filters$min_sample_size), n_before, n_after)

    if (verbose >= 1 && n_before > n_after) {
      message(sprintf("  Removed %d variants with N < %d",
                     n_before - n_after, filters$min_sample_size))
    }
  }

  # ========== Filter 3: P-value filter ==========

  if ("P" %in% names(data) && !is.null(filters$max_pvalue) && filters$max_pvalue < 1) {
    n_before <- nrow(data)
    data <- data[!is.na(data$P) & data$P <= filters$max_pvalue, ]
    n_after <- nrow(data)

    .record_filter(sprintf("P <= %g", filters$max_pvalue), n_before, n_after)

    if (verbose >= 1 && n_before > n_after) {
      message(sprintf("  Removed %d variants with P > %g",
                     n_before - n_after, filters$max_pvalue))
    }
  }

  # ========== Filter 4: MAF range filter ==========

  if ("MAF" %in% names(data)) {
    n_before <- nrow(data)

    # Apply min_maf filter
    if (!is.null(filters$min_maf) && filters$min_maf > 0) {
      data <- data[!is.na(data$MAF) & data$MAF >= filters$min_maf, ]
    }

    # Ensure MAF is in valid range (0, 0.5]
    # Note: MAF > 0.5 should already be converted to 1 - MAF in prepare_B_training_data
    # This just removes any remaining invalid values (MAF = 0 or MAF > 0.5 due to missing conversion)
    data <- data[!is.na(data$MAF) & data$MAF > 0 & data$MAF <= 0.5, ]

    n_after <- nrow(data)

    if (filters$min_maf > 0) {
      .record_filter(
        sprintf("MAF >= %g and MAF in (0, 0.5]", filters$min_maf),
        n_before, n_after
      )
    } else {
      .record_filter("MAF in (0, 0.5]", n_before, n_after)
    }

    if (verbose >= 1 && n_before > n_after) {
      message(sprintf("  Removed %d variants outside MAF range",
                     n_before - n_after))
    }
  }

  # ========== Filter 5: Remove duplicates ==========

  if (filters$remove_duplicates && "rsID" %in% names(data)) {
    n_before <- nrow(data)

    # Sort by duplicate_priority column (descending) before removing duplicates
    # This ensures we keep the variant with the highest priority value (e.g., largest N)
    if (!is.null(filters$duplicate_priority) &&
        filters$duplicate_priority %in% names(data)) {
      # Sort by rsID first (for stability), then by priority column (descending)
      data <- data[order(data$rsID, -data[[filters$duplicate_priority]]), ]

      if (verbose >= 2) {
        message(sprintf("  Sorting by duplicate_priority column '%s' (descending)",
                       filters$duplicate_priority))
      }
    }

    # Remove duplicates, keeping first occurrence (which is now the highest priority)
    data <- data[!duplicated(data$rsID), ]

    n_after <- nrow(data)
    .record_filter("Remove duplicate rsIDs", n_before, n_after)

    if (verbose >= 1 && n_before > n_after) {
      message(sprintf("  Removed %d duplicate rsIDs",
                     n_before - n_after))
    }
  }

  # ========== Summary ==========

  n_final <- nrow(data)

  if (verbose >= 1) {
    message(sprintf("  QC complete: %d -> %d variants (%.1f%% retained)",
                   n_start, n_final, 100 * n_final / n_start))
  }

  return(list(
    data = data,
    qc_summary = qc_summary
  ))
}


#' Validate GWAS Data for B Estimation
#'
#' @description
#' Validates that GWAS summary statistics meet the requirements for B
#' estimation. Checks for required columns, valid value ranges, and
#' data consistency.
#'
#' @param data data.frame with GWAS summary statistics
#' @param method Character string specifying B estimation method: "auto",
#'   "beta", or "pvalue". Determines which columns are required
#' @param verbose Integer verbosity level (0=silent, 1=messages)
#'
#' @return List with elements:
#'   \itemize{
#'     \item valid: Logical indicating if data is valid
#'     \item errors: Character vector of error messages (empty if valid)
#'     \item warnings: Character vector of warning messages
#'     \item summary: Summary statistics about the data
#'   }
#'
#' @details
#' \strong{Validation Checks:}
#'
#' Required for all methods:
#' \itemize{
#'   \item MAF column present and values in (0, 1)
#'   \item P column present and values in (0, 1)
#'   \item N column present and values > 0
#' }
#'
#' Additional for beta method:
#' \itemize{
#'   \item BETA column present and numeric
#' }
#'
#' Data consistency checks:
#' \itemize{
#'   \item No duplicate rsIDs (warning if present)
#'   \item Reasonable sample sizes (warning if very small or very large)
#'   \item MAF harmonized to <= 0.5
#' }
#'
#' Computational complexity: O(n x c) where n = number of variants, c = number of columns to validate.
#' Current implementation makes multiple passes through the data.
#'
#' @examples
#' \dontrun{
#' # Validate data
#' validation <- validate_B_training_data(my_data, method = "pvalue")
#'
#' if (!validation$valid) {
#'   stop("Validation failed: ", paste(validation$errors, collapse = "; "))
#' }
#' }
#'
#' @export
validate_B_training_data <- function(
  data,
  method = "auto",
  verbose = 1
) {

  errors <- character()
  warnings <- character()

  # ========== Check Required Columns ==========

  # Columns required for all methods
  remove_na_cols <- c("MAF", "P", "N")

  # Additional columns for specific methods
  if (method == "beta") {
    remove_na_cols <- c(remove_na_cols, "BETA")
  }

  # Check if required columns exist
  missing_cols <- setdiff(remove_na_cols, names(data))

  if (length(missing_cols) > 0) {
    errors <- c(errors, paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    ))
    # Cannot proceed with other checks if columns missing
    return(list(
      valid = FALSE,
      errors = errors,
      warnings = warnings,
      summary = NULL
    ))
  }

  # ========== Validate MAF ==========

  if ("MAF" %in% names(data)) {
    # Check range
    if (any(data$MAF <= 0 | data$MAF > 0.5, na.rm = TRUE)) {
      errors <- c(errors, "MAF values must be in the interval (0, 0.5]")
    }

    # Check for NAs
    if (any(is.na(data$MAF))) {
      n_na <- sum(is.na(data$MAF))
      errors <- c(errors, sprintf("%d variants have NA for MAF", n_na))
    }

    # Check if all MAFs are very small or very large
    if (all(data$MAF < 0.001, na.rm = TRUE)) {
      warnings <- c(warnings, "All MAF values are < 0.001 (very rare variants)")
    }
    if (all(data$MAF > 0.4, na.rm = TRUE)) {
      warnings <- c(warnings, "All MAF values are > 0.4 (very common variants)")
    }
  }

  # ========== Validate P-values ==========

  if ("P" %in% names(data)) {
    # Check range
    if (any(data$P < 0 | data$P > 1, na.rm = TRUE)) {
      errors <- c(errors, "P-values must be in the interval [0, 1]")
    }

    # Check for NAs
    if (any(is.na(data$P))) {
      n_na <- sum(is.na(data$P))
      errors <- c(errors, sprintf("%d variants have NA for P", n_na))
    }

    # Check for p-values of exactly 0 or 1
    if (any(data$P == 0, na.rm = TRUE)) {
      warnings <- c(warnings, "Some P-values are exactly 0 (may cause numerical issues)")
    }
    if (any(data$P == 1, na.rm = TRUE)) {
      warnings <- c(warnings, "Some P-values are exactly 1")
    }
  }

  # ========== Validate Sample Size ==========

  if ("N" %in% names(data)) {
    # Check range
    if (any(data$N <= 0, na.rm = TRUE)) {
      errors <- c(errors, "Sample sizes (N) must be positive")
    }

    # Check for NAs
    if (any(is.na(data$N))) {
      n_na <- sum(is.na(data$N))
      errors <- c(errors, sprintf("%d variants have NA for N", n_na))
    }

    # Check for very small sample sizes
    if (any(data$N < 100, na.rm = TRUE)) {
      n_small <- sum(data$N < 100, na.rm = TRUE)
      warnings <- c(warnings, sprintf(
        "%d variants have N < 100 (very small sample size)",
        n_small
      ))
    }

    # Check for very large sample sizes (potential data entry errors)
    if (any(data$N > 10000000, na.rm = TRUE)) {
      n_large <- sum(data$N > 10000000, na.rm = TRUE)
      warnings <- c(warnings, sprintf(
        "%d variants have N > 10 million (check for data entry errors)",
        n_large
      ))
    }
  }

  # ========== Validate BETA (if present and required) ==========

  if ("BETA" %in% names(data)) {
    # Check for NAs
    if (any(is.na(data$BETA))) {
      n_na <- sum(is.na(data$BETA))
      if (method == "beta") {
        errors <- c(errors, sprintf("%d variants have NA for BETA", n_na))
      } else {
        warnings <- c(warnings, sprintf("%d variants have NA for BETA", n_na))
      }
    }

    # Check if all BETAs are 0
    if (all(data$BETA == 0, na.rm = TRUE)) {
      warnings <- c(warnings, "All BETA values are 0")
    }

    # Check for extreme BETAs (potential data entry errors)
    if (any(abs(data$BETA) > 10, na.rm = TRUE)) {
      n_extreme <- sum(abs(data$BETA) > 10, na.rm = TRUE)
      warnings <- c(warnings, sprintf(
        "%d variants have |BETA| > 10 (check for data entry errors)",
        n_extreme
      ))
    }
  }

  # ========== Check for Duplicates ==========

  if ("rsID" %in% names(data)) {
    duplicated_rsIDs <- data$rsID[duplicated(data$rsID)]

    if (length(duplicated_rsIDs) > 0) {
      warnings <- c(warnings, sprintf(
        "%d duplicate rsIDs found (consider using remove_duplicates=TRUE)",
        length(duplicated_rsIDs)
      ))
    }
  }

  # ========== Create Summary Statistics ==========

  summary_stats <- list(
    n_variants = nrow(data),
    n_complete_cases = sum(complete.cases(data)),
    maf_range = if ("MAF" %in% names(data)) range(data$MAF, na.rm = TRUE) else NULL,
    p_range = if ("P" %in% names(data)) range(data$P, na.rm = TRUE) else NULL,
    n_range = if ("N" %in% names(data)) range(data$N, na.rm = TRUE) else NULL,
    beta_range = if ("BETA" %in% names(data)) range(data$BETA, na.rm = TRUE) else NULL
  )

  # ========== Report Results ==========

  valid <- length(errors) == 0

  if (verbose >= 1) {
    if (valid) {
      message("  Validation passed")
      if (length(warnings) > 0) {
        message("  Warnings:")
        for (w in warnings) {
          message("    - ", w)
        }
      }
    } else {
      message("  Validation failed:")
      for (e in errors) {
        message("    - ", e)
      }
    }
  }

  return(list(
    valid = valid,
    errors = errors,
    warnings = warnings,
    summary = summary_stats
  ))
}


#################### INTERNAL HELPER FUNCTIONS ####################


#' Read Excel File (XLSX/XLS)
#'
#' @description
#' Internal helper to read Excel files using the readxl package.
#'
#' @param file_path Path to Excel file (.xlsx or .xls)
#' @param verbose Verbosity level
#'
#' @return data.frame with file contents
#'
#' @keywords internal
#' @noRd
.read_excel_file <- function(file_path, verbose = 1) {

  # Check if readxl is available
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read Excel files. ",
         "Install with: install.packages('readxl')")
  }

  if (verbose >= 2) {
    message("  Using readxl::read_excel() to read Excel file")
  }

  # Read Excel file
  # readxl automatically handles .xlsx and .xls formats
  data <- readxl::read_excel(
    path = file_path,
    na = c("NA", "na", "N/A", "n/a", "", " ", ".")
  )

  # Convert to data.frame (readxl returns tibble)
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (verbose >= 1) {
    message(sprintf("  Loaded Excel file with %d variants and %d columns",
                   nrow(data), ncol(data)))
  }

  return(data)
}


#' Read GWAS Catalog Format File
#'
#' @description
#' Internal helper to read files in GWAS Catalog format (tab-delimited).
#'
#' @param file_path Path to file
#' @param verbose Verbosity level
#'
#' @return data.frame with GWAS data
#'
#' @keywords internal
#' @noRd
.read_gwas_catalog <- function(file_path, verbose = 1) {
  # GWAS Catalog files are typically tab-delimited
  data <- .read_csv_file(file_path, delimiter = "\t", verbose = verbose)

  return(data)
}


#' Read CSV/TSV File with Flexible Delimiter Detection
#'
#' @description
#' Internal helper to read CSV or TSV files with automatic delimiter detection.
#'
#' @param file_path Path to file
#' @param delimiter Character delimiter (default: auto-detect)
#' @param verbose Verbosity level
#'
#' @return data.frame with file contents
#'
#' @keywords internal
#' @noRd
.read_csv_file <- function(file_path, delimiter = NULL, verbose = 1) {

  # Auto-detect delimiter if not specified
  if (is.null(delimiter)) {
    # Read first line
    first_line <- readLines(file_path, n = 1)

    # Count occurrences of common delimiters using lengths() for efficiency
    n_comma <- lengths(gregexpr(",", first_line, fixed = TRUE))
    n_tab <- lengths(gregexpr("\t", first_line, fixed = TRUE))
    n_space <- lengths(gregexpr(" ", first_line, fixed = TRUE))

    # Choose delimiter with most occurrences
    if (n_tab > n_comma && n_tab > n_space) {
      delimiter <- "\t"
      delim_name <- "tab"
    } else if (n_comma > n_space) {
      delimiter <- ","
      delim_name <- "comma"
    } else {
      delimiter <- " "
      delim_name <- "space"
    }

    if (verbose >= 2) {
      message(sprintf("  Auto-detected delimiter: %s", delim_name))
    }
  }

  # Try to use data.table::fread for 10-100x faster file reading
  # Falls back to read.table if data.table is not available
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (verbose >= 2) {
      message("  Using data.table::fread() for fast file reading")
    }

    data <- data.table::fread(
      file_path,
      sep = delimiter,
      header = TRUE,
      data.table = FALSE,  # Return data.frame not data.table
      na.strings = c("NA", "na", "N/A", "n/a", "", " "),
      showProgress = (verbose >= 2)
    )
  } else {
    # Fallback to read.table if data.table not available
    if (verbose >= 1) {
      message("  Note: Install 'data.table' package for 10-100x faster file reading")
    }

    data <- read.table(
      file_path,
      header = TRUE,
      sep = delimiter,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "\"",
      na.strings = c("NA", "na", "N/A", "n/a", "", " ")
    )
  }

  if (verbose >= 1) {
    message(sprintf("  Loaded file with %d variants and %d columns",
                   nrow(data), ncol(data)))
  }

  return(data)
}



#################### S3 METHODS ####################


#' Print Method for glow_training_data
#'
#' @description
#' Prints a summary of prepared GWAS training data.
#'
#' @param x Object of class \code{glow_training_data}
#' @param ... Additional arguments (unused)
#'
#' @export
print.glow_training_data <- function(x, ...) {
  cat("GLOW Training Data for B Estimation\n")
  cat("====================================\n\n")

  cat("Metadata:\n")
  cat("  Trait type:       ", x$metadata$trait_type, "\n")
  cat("  Data format:      ", x$metadata$format, "\n")
  cat("  Variants:         ", x$metadata$n_variants_final,
      " (", x$metadata$n_variants_removed, " removed during QC)\n", sep = "")
  cat("  Date prepared:    ", as.character(x$metadata$date_prepared), "\n")

  cat("\nAvailable columns: ", paste(x$metadata$columns_present, collapse = ", "), "\n")

  cat("\nData summary:\n")
  if ("MAF" %in% names(x$data)) {
    cat("  MAF range:        ", sprintf("%.4f - %.4f",
                                       min(x$data$MAF, na.rm = TRUE),
                                       max(x$data$MAF, na.rm = TRUE)), "\n")
  }
  if ("P" %in% names(x$data)) {
    cat("  P-value range:    ", sprintf("%.2e - %.2e",
                                       min(x$data$P, na.rm = TRUE),
                                       max(x$data$P, na.rm = TRUE)), "\n")
  }
  if ("N" %in% names(x$data)) {
    cat("  Sample size (N):  ", sprintf("%.0f - %.0f (median: %.0f)",
                                       min(x$data$N, na.rm = TRUE),
                                       max(x$data$N, na.rm = TRUE),
                                       median(x$data$N, na.rm = TRUE)), "\n")
  }
  if ("BETA" %in% names(x$data)) {
    cat("  BETA range:       ", sprintf("%.4f - %.4f",
                                       min(x$data$BETA, na.rm = TRUE),
                                       max(x$data$BETA, na.rm = TRUE)), "\n")
  }

  cat("\nQC filters applied:\n")
  if (!is.null(x$metadata$qc_applied$min_sample_size)) {
    cat("  Min sample size:  ", x$metadata$qc_applied$min_sample_size, "\n")
  }
  if (!is.null(x$metadata$qc_applied$min_maf) && x$metadata$qc_applied$min_maf > 0) {
    cat("  Min MAF:          ", sprintf("%.3f", x$metadata$qc_applied$min_maf), "\n")
  }
  if (!is.null(x$metadata$qc_applied$remove_duplicates) && x$metadata$qc_applied$remove_duplicates) {
    cat("  Duplicates:       ", "Removed\n")
  }

  cat("\nReady for use with get_B()\n")

  invisible(x)
}


#' Summary Method for glow_training_data
#'
#' @description
#' Provides a detailed summary of prepared GWAS training data.
#'
#' @param object Object of class \code{glow_training_data}
#' @param ... Additional arguments (unused)
#'
#' @export
summary.glow_training_data <- function(object, ...) {
  cat("GLOW Training Data Summary\n")
  cat("==========================\n\n")

  # Print basic info
  print(object)

  # Add QC summary details if available
  if (!is.null(object$metadata$qc_summary) && nrow(object$metadata$qc_summary) > 0) {
    cat("\n\nDetailed QC Summary:\n")
    cat("--------------------\n")
    print(object$metadata$qc_summary, row.names = FALSE)
  }

  # Add data structure info
  cat("\n\nData Structure:\n")
  cat("---------------\n")
  str(object$data, max.level = 1, vec.len = 3)

  invisible(object)
}
