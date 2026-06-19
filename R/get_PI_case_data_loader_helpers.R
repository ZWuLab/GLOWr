########## Helper Functions for prepare_PI_case_data() ##########
#
# This file contains internal helper functions specific to PI case data
# preparation, including column mappings, VarInfo generation, chromosome
# standardization, and author/chromosome filtering.

#################### EXPORTED MAIN FUNCTIONS ####################

# None - this file contains internal helper functions only

#################### INTERNAL HELPER FUNCTIONS ####################

# ---------- Column Mapping Dictionary ----------

#' Get PI-Specific Column Name Variations
#'
#' @description
#' Returns a named list of column name variations specific to PI case data
#' preparation. These mappings extend the standard mappings in
#' \code{standardize_column_names} with additional patterns commonly found
#' in trait-associated variant lists.
#'
#' @return Named list where each element is a standard column name containing
#'   a character vector of variations that should map to that standard name
#'
#' @details
#' \strong{PI-Specific Column Mappings:}
#'
#' This function provides column name variations commonly found in:
#' \itemize{
#'   \item GWAS Catalog downloads
#'   \item Literature-curated variant lists (e.g., ALS-known-SNPs.xlsx)
#'   \item Published supplementary tables
#'   \item Meta-analysis summary statistics
#' }
#'
#' \strong{Standard Column Names:}
#' \itemize{
#'   \item rsID: SNP identifier
#'   \item CHR: Chromosome
#'   \item POS: Position
#'   \item REF: Reference allele
#'   \item ALT: Alternate allele
#'   \item MAF: Minor allele frequency
#'   \item P: P-value
#'   \item N: Sample size
#'   \item BETA: Effect size
#'   \item TRAIT: Trait name
#'   \item STUDY: Study name
#'   \item PMID: PubMed ID
#'   \item FIRST_AUTHOR: First author name
#' }
#'
#' @keywords internal
#' @noRd
.get_pi_column_variations <- function() {

  list(
    # SNP identifier
    rsID = c("SNP", "SNPS", "RSID", "RS_ID", "VARIANT_ID", "ID",
             "MARKERNAME", "MARKER", "SNP_ID_CURRENT", "SNPID"),

    # Chromosome
    CHR = c("CHR", "CHROM", "CHROMOSOME", "CHR_ID", "#CHR", "chr"),

    # Position
    POS = c("POS", "POSITION", "BP", "CHR_POS", "BASE_PAIR_LOCATION",
            "CHR_POSITION", "POSITION_HG38", "POSITION_HG19"),

    # Reference allele
    REF = c("REF", "REFERENCE", "REF_ALLELE", "REFERENCE_ALLELE",
            "ALLELE1", "A1", "EFFECT_ALLELE", "EA"),

    # Alternate allele
    ALT = c("ALT", "ALTERNATE", "ALT_ALLELE", "ALTERNATE_ALLELE",
            "ALLELE2", "A2", "OTHER_ALLELE", "OA", "NEA"),

    # Minor allele frequency
    MAF = c("MAF", "RAF", "FREQ", "EAF", "AF",
            "RISK.ALLELE.FREQUENCY", "RISK_ALLELE_FREQUENCY",
            "EFFECT_ALLELE_FREQUENCY", "ALT_FREQS",
            "EFFECT ALLELE FREQUENCY"),

    # P-value
    P = c("P", "PVALUE", "P_VALUE", "P.VALUE", "P-VALUE", "PVAL",
          "P-VALUE (TEXT)", "P.VALUE.TEXT", "P_VALUE_TEXT"),

    # Sample size (many text-based variations)
    N = c("N", "SAMPLE_SIZE", "SAMPLE.SIZE", "SAMPLE SIZE", "Sample Size",
          "INITIAL.SAMPLE.SIZE", "INITIAL SAMPLE SIZE",
          "N_TOTAL", "NEFF", "N_EFF", "TOTAL_N",
          "DISCOVERY SAMPLE SIZE", "REPLICATION SAMPLE SIZE"),

    # Effect size
    BETA = c("BETA", "EFFECT", "EFFECT_SIZE", "B", "BETA_NUMERIC",
             "Beta_numeric", "OR (TEXT)", "OR.95..CI.TEXT."),

    # Odds ratio (for binary traits)
    OR = c("OR", "ODDS_RATIO", "ODDSRATIO", "OR or 95% CI"),

    # Trait name
    TRAIT = c("TRAIT", "REPORTED TRAIT", "Reported trait", "REPORTED.TRAIT",
              "DISEASE/TRAIT", "DISEASE_TRAIT", "MAPPED_TRAIT",
              "MAPPED TRAIT", "PHENOTYPE", "TRAIT(S)"),

    # Study name
    STUDY = c("STUDY", "STUDY_NAME", "STUDY.NAME", "STUDY ACCESSION",
              "STUDY_ACCESSION", "PUBMEDID", "JOURNAL"),

    # PubMed ID
    PMID = c("PMID", "PUBMEDID", "PUBMED_ID", "PUBMED ID",
             "STUDY ACCESSION"),

    # First author (PI-specific - common in curated lists)
    FIRST_AUTHOR = c("FIRST AUTHOR", "First Author", "FIRST.AUTHOR",
                     "FIRST_AUTHOR", "AUTHOR", "AUTHORS",
                     "LEAD_AUTHOR", "LEAD AUTHOR")
  )
}

# ---------- VarInfo Generation ----------

#' Generate VarInfo Column in CHR-POS-REF-ALT Format
#'
#' @description
#' Creates a VarInfo column in FAVOR-compatible format (CHR-POS-REF-ALT).
#' When REF or ALT alleles are missing, uses "NA" as placeholder with warning.
#'
#' @param chr Numeric or character vector of chromosome values
#' @param pos Numeric vector of positions
#' @param ref Character vector of reference alleles (optional)
#' @param alt Character vector of alternate alleles (optional)
#' @param verbose Integer verbosity level (0=silent, 1=messages)
#'
#' @return Character vector in format "CHR-POS-REF-ALT" or "CHR-POS-NA-NA"
#'
#' @details
#' \strong{FAVOR Format Requirements:}
#'
#' The VarInfo format is used by FAVOR annotation database for variant lookup:
#' \itemize{
#'   \item Standard format: "1-12345-A-G" (CHR-POS-REF-ALT)
#'   \item Missing alleles: "1-12345-NA-NA" (with warning)
#'   \item Multi-allelic: Takes first ALT allele only
#' }
#'
#' \strong{Handling Missing Alleles:}
#'
#' Many literature-curated variant lists lack REF/ALT information (only rsID
#' and position available). This function allows "NA" placeholders to enable
#' FAVOR annotation by position:
#' \itemize{
#'   \item If REF is NULL or all NA: Uses "NA" for all variants (with warning)
#'   \item If ALT is NULL or all NA: Uses "NA" for all variants (with warning)
#'   \item Individual NA values: Replaced with "NA" string
#' }
#'
#' \strong{Multi-Allelic Variants:}
#'
#' For ALT alleles containing multiple alleles (comma-separated):
#' \itemize{
#'   \item Takes only the first allele
#'   \item Example: "A,G" -> "A" (first allele used)
#'   \item No warning issued (common in GWAS data)
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n) where n = number of variants. Single-pass vectorized string operations.
#'
#' @examples
#' \dontrun{
#' # Standard VarInfo generation
#' .generate_varinfo_column(c(1, 2), c(12345, 67890), c("A", "G"), c("G", "T"))
#' # Returns: c("1-12345-A-G", "2-67890-G-T")
#'
#' # Missing alleles (with warning)
#' .generate_varinfo_column(c(1, 2), c(12345, 67890), NULL, NULL, verbose = 1)
#' # Returns: c("1-12345-NA-NA", "2-67890-NA-NA")
#'
#' # Multi-allelic (takes first)
#' .generate_varinfo_column(1, 12345, "A", "G,T", verbose = 0)
#' # Returns: "1-12345-A-G"
#' }
#'
#' @keywords internal
#' @noRd
.generate_varinfo_column <- function(chr, pos, ref = NULL, alt = NULL, verbose = 1) {

  # Validate required inputs
  if (is.null(chr) || is.null(pos)) {
    stop("CHR and POS are required for VarInfo generation")
  }

  if (length(chr) != length(pos)) {
    stop("CHR and POS must have the same length")
  }

  n <- length(chr)

  # Handle missing REF/ALT
  # If entire column is NULL or all NA, use "NA" for all
  if (is.null(ref) || all(is.na(ref))) {
    ref <- rep("NA", n)
    if (verbose >= 1) {
      warning("REF allele not available. Using 'NA' placeholder in VarInfo format.",
              call. = FALSE)
    }
  }

  if (is.null(alt) || all(is.na(alt))) {
    alt <- rep("NA", n)
    if (verbose >= 1) {
      warning("ALT allele not available. Using 'NA' placeholder in VarInfo format.",
              call. = FALSE)
    }
  }

  # Ensure vectors are same length
  if (length(ref) != n) {
    stop("REF must have same length as CHR/POS")
  }
  if (length(alt) != n) {
    stop("ALT must have same length as CHR/POS")
  }

  # Convert to character
  chr <- as.character(chr)
  pos <- as.character(pos)
  ref <- as.character(ref)
  alt <- as.character(alt)

  # Replace individual NA values with "NA" string
  ref[is.na(ref)] <- "NA"
  alt[is.na(alt)] <- "NA"

  # Handle multi-allelic variants (take first allele)
  # Split by comma and take first element
  alt <- sapply(strsplit(alt, ",", fixed = TRUE), function(x) x[1])

  # Generate VarInfo in format: CHR-POS-REF-ALT
  varinfo <- paste(chr, pos, ref, alt, sep = "-")

  return(varinfo)
}


#' Fill Missing Position Values from Alternative Columns
#'
#' @description
#' When the POS column has missing values, searches for alternative position
#' columns (e.g., pos_GRCh37, pos_GRCh38, position) and fills in missing values.
#'
#' @param data data.frame containing position columns
#' @param verbose Integer verbosity level
#'
#' @return data.frame with POS column filled from alternative sources
#'
#' @details
#' \strong{Search Order for Position Columns:}
#' \enumerate{
#'   \item POS (primary - already standardized)
#'   \item pos_GRCh38 (preferred build)
#'   \item position_GRCh38, POSITION_HG38
#'   \item pos_GRCh37, position_GRCh37, POSITION_HG19
#'   \item Any column containing "pos" (case-insensitive)
#' }
#'
#' \strong{IMPORTANT - GRCh38 Requirement:}
#' FAVOR annotation database requires GRCh38 coordinates. Users are responsible
#' for ensuring all position values are based on GRCh38. If positions from
#' older builds (e.g., GRCh37/hg19) are used, FAVOR annotation will fail or
#' return incorrect results.
#'
#' @keywords internal
#' @noRd
.fill_missing_positions <- function(data, verbose = 1) {

  # Check if POS column exists

  if (!"POS" %in% names(data)) {
    return(data)
  }

  # Count missing positions
  n_missing <- sum(is.na(data$POS))
  if (n_missing == 0) {
    return(data)
  }

  if (verbose >= 1) {
    message(sprintf("  POS column has %d missing values (%.1f%%). Searching for alternatives...",
                   n_missing, 100 * n_missing / nrow(data)))
  }

  # Define alternative position column patterns (in priority order)
  alt_patterns <- c(
    "^pos_GRCh38$", "^position_GRCh38$", "^POSITION_HG38$",  # GRCh38 preferred
    "^pos_GRCh37$", "^position_GRCh37$", "^POSITION_HG19$",  # GRCh37 fallback
    "^pos$", "^position$", "^bp$"  # Generic fallback
  )

  col_names <- names(data)
  filled_from <- NULL

  # Search for alternative columns
  for (pattern in alt_patterns) {
    matches <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
    for (alt_col in matches) {
      if (alt_col == "POS") next  # Skip the main column

      alt_values <- data[[alt_col]]
      if (is.null(alt_values) || all(is.na(alt_values))) next

      # Fill missing POS values from this column
      missing_idx <- is.na(data$POS)
      has_alt <- !is.na(alt_values)
      fill_idx <- missing_idx & has_alt

      if (any(fill_idx)) {
        data$POS[fill_idx] <- alt_values[fill_idx]
        n_filled <- sum(fill_idx)
        filled_from <- alt_col

        if (verbose >= 1) {
          message(sprintf("  Filled %d missing POS values from '%s'", n_filled, alt_col))

          # Warn if column name suggests non-GRCh38 build
          if (grepl("GRCh37|hg19|HG19", alt_col, ignore.case = TRUE)) {
            warning(sprintf(
              paste0("Position values filled from '%s' which appears to be GRCh37/hg19.\n",
                     "FAVOR annotation requires GRCh38 coordinates.\n",
                     "Please ensure all positions are lifted over to GRCh38 before FAVOR annotation."),
              alt_col
            ), call. = FALSE)
          }
        }
        break  # Stop after finding first useful alternative
      }
    }
    if (!is.null(filled_from)) break
  }

  # Report remaining missing positions
  n_still_missing <- sum(is.na(data$POS))
  if (n_still_missing > 0 && verbose >= 1) {
    message(sprintf("  %d positions still missing after fallback search", n_still_missing))
  }

  return(data)
}


#' Standardize Chromosome Values to Numeric 1-22
#'
#' @description
#' Converts chromosome values to numeric format (1-22 for autosomes).
#' Recognizes common chromosome naming conventions and removes non-autosomal
#' chromosomes (X, Y, MT).
#'
#' @param chr_values Character or numeric vector of chromosome values
#' @param verbose Integer verbosity level
#'
#' @return Numeric vector with autosomes (1-22) and NA for non-autosomes
#'
#' @details
#' \strong{Recognized Patterns:}
#' \itemize{
#'   \item Autosomes: "1", "chr1", "CHR1" -> 1-22
#'   \item X chromosome: "X", "chrX", "23" -> NA
#'   \item Y chromosome: "Y", "chrY", "24" -> NA
#'   \item Mitochondrial: "MT", "M", "chrM", "25" -> NA
#'   \item Invalid: Non-standard names -> NA
#' }
#'
#' \strong{Processing Steps:}
#' \enumerate{
#'   \item Convert to uppercase
#'   \item Remove "CHR" prefix if present
#'   \item Map X/Y/MT to NA
#'   \item Convert remaining to numeric
#'   \item Set values outside 1-22 to NA
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n) where n = number of chromosome values. Single-pass conversion.
#'
#' @examples
#' \dontrun{
#' # Standard conversion
#' .standardize_chromosome(c("1", "2", "X", "chr5", "Y"))
#' # Returns: c(1, 2, NA, 5, NA)
#'
#' # Already numeric
#' .standardize_chromosome(c(1, 2, 3))
#' # Returns: c(1, 2, 3)
#' }
#'
#' @keywords internal
#' @noRd
.standardize_chromosome <- function(chr_values, verbose = 1) {

  # Convert to character
  chr_char <- as.character(chr_values)

  # Convert to uppercase for matching
  chr_upper <- toupper(chr_char)

  # Remove "CHR" prefix if present
  chr_clean <- sub("^CHR", "", chr_upper)

  # Initialize output
  chr_numeric <- rep(NA_real_, length(chr_values))

  # Identify non-autosomal chromosomes
  is_x <- chr_clean %in% c("X", "23")
  is_y <- chr_clean %in% c("Y", "24")
  is_mt <- chr_clean %in% c("MT", "M", "25")

  # Set non-autosomes to NA
  chr_numeric[is_x | is_y | is_mt] <- NA

  # Convert remaining to numeric
  is_autosome <- !(is_x | is_y | is_mt)
  chr_numeric[is_autosome] <- suppressWarnings(as.numeric(chr_clean[is_autosome]))

  # Validate range (1-22 for autosomes)
  out_of_range <- !is.na(chr_numeric) & (chr_numeric < 1 | chr_numeric > 22)
  chr_numeric[out_of_range] <- NA

  # Report statistics if verbose
  if (verbose >= 2) {
    n_x <- sum(is_x)
    n_y <- sum(is_y)
    n_mt <- sum(is_mt)
    n_invalid <- sum(is.na(chr_numeric) & !is_x & !is_y & !is_mt)
    n_autosomes <- sum(!is.na(chr_numeric))

    if (n_x > 0) message(sprintf("  Identified %d X chromosome variants", n_x))
    if (n_y > 0) message(sprintf("  Identified %d Y chromosome variants", n_y))
    if (n_mt > 0) message(sprintf("  Identified %d mitochondrial variants", n_mt))
    if (n_invalid > 0) message(sprintf("  %d variants with invalid CHR", n_invalid))
    message(sprintf("  %d autosomal variants (CHR 1-22)", n_autosomes))
  }

  return(chr_numeric)
}


#' Filter Data by Author Name
#'
#' @description
#' Removes rows matching specified author names. Uses case-insensitive
#' partial matching to handle variations in author name formatting.
#'
#' @param data data.frame to filter
#' @param exclude_authors Character vector of author names to exclude
#' @param author_col Character string specifying author column name
#'   (default: "FIRST_AUTHOR")
#' @param verbose Integer verbosity level
#'
#' @return data.frame with specified authors removed
#'
#' @details
#' \strong{Matching Strategy:}
#'
#' Uses case-insensitive partial matching (grepl) to handle:
#' \itemize{
#'   \item Name variations: "Nicolas A", "Nicolas A.", "Nicolas"
#'   \item Formatting differences: "Smith J", "Smith, J.", "Smith"
#'   \item Multiple authors in single field: "Smith J, Jones K"
#' }
#'
#' \strong{Common Use Cases:}
#' \itemize{
#'   \item Removing overlapping studies from same research group
#'   \item Excluding variants from specific publications
#'   \item Filtering out low-quality studies by author
#' }
#'
#' \strong{Important Notes:}
#' \itemize{
#'   \item Partial matching may remove more variants than intended
#'   \item Verify filtered results to ensure correct removal
#'   \item Use specific patterns to avoid over-filtering
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n * m) where n = number of variants, m = number of authors to exclude.
#' Uses grepl for each author pattern.
#'
#' @examples
#' \dontrun{
#' # Exclude specific author
#' filtered <- .filter_by_author(data, exclude_authors = "Nicolas A")
#'
#' # Exclude multiple authors
#' filtered <- .filter_by_author(data, exclude_authors = c("Nicolas A", "Smith J"))
#' }
#'
#' @keywords internal
#' @noRd
.filter_by_author <- function(data, exclude_authors = NULL,
                              author_col = "FIRST_AUTHOR", verbose = 1) {

  # Return unchanged if no filtering requested
  if (is.null(exclude_authors) || length(exclude_authors) == 0) {
    return(data)
  }

  # Check if author column exists
  if (!author_col %in% names(data)) {
    if (verbose >= 1) {
      warning(sprintf("Author column '%s' not found in data. Skipping author filtering.",
                     author_col))
    }
    return(data)
  }

  n_before <- nrow(data)

  # Create logical vector for rows to remove
  remove_idx <- rep(FALSE, nrow(data))

  # Get author column
  author_values <- as.character(data[[author_col]])

  # Check each excluded author pattern (case-insensitive partial matching)
  for (author in exclude_authors) {
    # Use grepl with ignore.case for flexible matching
    matches <- grepl(author, author_values, ignore.case = TRUE, fixed = FALSE)
    remove_idx <- remove_idx | matches

    if (verbose >= 2) {
      n_matched <- sum(matches, na.rm = TRUE)
      if (n_matched > 0) {
        message(sprintf("  Excluding %d variants with author pattern '%s'",
                       n_matched, author))
      }
    }
  }

  # Remove matched rows
  data <- data[!remove_idx, , drop = FALSE]

  n_after <- nrow(data)
  n_removed <- n_before - n_after

  if (verbose >= 1 && n_removed > 0) {
    message(sprintf("  Removed %d variants from %d author(s)",
                   n_removed, length(exclude_authors)))
  }

  return(data)
}


#' Filter Data by Chromosome
#'
#' @description
#' Removes rows with specified chromosome values. Handles both numeric
#' and character chromosome representations.
#'
#' @param data data.frame to filter
#' @param exclude_chromosomes Character or numeric vector of chromosomes to exclude
#' @param chr_col Character string specifying chromosome column name
#'   (default: "CHR")
#' @param verbose Integer verbosity level
#'
#' @return data.frame with specified chromosomes removed
#'
#' @details
#' \strong{Matching Strategy:}
#'
#' Uses case-insensitive matching to handle:
#' \itemize{
#'   \item Numeric: 23, 24, 25 (for X, Y, MT)
#'   \item Character: "X", "Y", "MT", "M"
#'   \item With prefix: "chr23", "chrX"
#' }
#'
#' \strong{Common Use Cases:}
#' \itemize{
#'   \item Excluding X chromosome: exclude_chromosomes = "X"
#'   \item Excluding sex chromosomes: exclude_chromosomes = c("X", "Y")
#'   \item Excluding MT: exclude_chromosomes = c("MT", "M")
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n * m) where n = number of variants, m = number of chromosomes to exclude.
#'
#' @examples
#' \dontrun{
#' # Exclude X chromosome
#' filtered <- .filter_by_chromosome(data, exclude_chromosomes = "X")
#'
#' # Exclude sex chromosomes
#' filtered <- .filter_by_chromosome(data, exclude_chromosomes = c("X", "Y"))
#' }
#'
#' @keywords internal
#' @noRd
.filter_by_chromosome <- function(data, exclude_chromosomes = NULL,
                                  chr_col = "CHR", verbose = 1) {

  # Return unchanged if no filtering requested
  if (is.null(exclude_chromosomes) || length(exclude_chromosomes) == 0) {
    return(data)
  }

  # Check if chromosome column exists
  if (!chr_col %in% names(data)) {
    if (verbose >= 1) {
      warning(sprintf("Chromosome column '%s' not found in data. Skipping chromosome filtering.",
                     chr_col))
    }
    return(data)
  }

  n_before <- nrow(data)

  # Get chromosome column
  chr_values <- data[[chr_col]]

  # Convert to character for matching
  chr_char <- toupper(as.character(chr_values))

  # Standardize excluded chromosomes for matching
  exclude_std <- toupper(as.character(exclude_chromosomes))

  # Also handle "CHR" prefix removal
  chr_char_clean <- sub("^CHR", "", chr_char)

  # Create logical vector for rows to remove
  remove_idx <- rep(FALSE, nrow(data))

  # Check each excluded chromosome
  for (excl in exclude_std) {
    # Remove "CHR" prefix if present
    excl_clean <- sub("^CHR", "", excl)

    # Match both original and cleaned versions
    matches <- (chr_char == excl) | (chr_char_clean == excl_clean)
    remove_idx <- remove_idx | matches

    if (verbose >= 2) {
      n_matched <- sum(matches, na.rm = TRUE)
      if (n_matched > 0) {
        message(sprintf("  Excluding %d variants on chromosome '%s'",
                       n_matched, excl))
      }
    }
  }

  # Remove matched rows
  data <- data[!remove_idx, , drop = FALSE]

  n_after <- nrow(data)
  n_removed <- n_before - n_after

  if (verbose >= 1 && n_removed > 0) {
    message(sprintf("  Removed %d variants from %d chromosome(s)",
                   n_removed, length(exclude_chromosomes)))
  }

  return(data)
}


#' Reorder Columns for PI Case Data (VarInfo First)
#'
#' @description
#' Reorders columns to a standard format for PI case data with VarInfo
#' as the first column for easy FAVOR annotation.
#'
#' @param data data.frame to reorder
#'
#' @return data.frame with reordered columns
#'
#' @details
#' \strong{Column Order:}
#' \enumerate{
#'   \item VarInfo (first - for FAVOR lookup)
#'   \item rsID, CHR, POS (identifiers)
#'   \item REF, ALT, MAF (allele info)
#'   \item P, N, BETA (statistics)
#'   \item TRAIT, STUDY, PMID, FIRST_AUTHOR (metadata)
#'   \item Everything else (alphabetically)
#' }
#'
#' \strong{Rationale:}
#' VarInfo is placed first because it is the primary key for FAVOR annotation
#' lookup. This makes it easy to extract the column for annotation workflows.
#'
#' @keywords internal
#' @noRd
.reorder_columns_pi <- function(data) {

  # Define column order (VarInfo first, then standard columns)
  priority_cols <- c("VarInfo", "rsID", "CHR", "POS", "REF", "ALT", "MAF",
                     "P", "N", "BETA", "TRAIT_TYPE")
  metadata_cols <- c("TRAIT", "STUDY", "PMID", "FIRST_AUTHOR")

  # Get current column names
  current_cols <- names(data)

  # Identify which priority columns are present (in order)
  present_priority <- priority_cols[priority_cols %in% current_cols]
  present_metadata <- metadata_cols[metadata_cols %in% current_cols]

  # Identify remaining columns (not in priority or metadata lists)
  remaining_cols <- setdiff(current_cols, c(priority_cols, metadata_cols))
  remaining_cols <- sort(remaining_cols)  # Alphabetically

  # Create final column order
  final_order <- c(present_priority, present_metadata, remaining_cols)

  # Reorder
  data <- data[, final_order, drop = FALSE]

  return(data)
}
