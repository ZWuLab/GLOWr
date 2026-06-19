########## FAVOR Annotation Module ##########
#
# This file provides functions to annotate variants with functional annotation
# scores from the FAVOR database. Used in PI estimation to assign variant-importance
# scores based on functional features.
#
# EXPORTED FUNCTIONS:
#   - annotate_favor()   Annotate variants with FAVOR scores (CSV / aGDS output)
#
# INTERNAL HELPERS (selected):
#   - .create_agds_from_gds()                  aGDS output for GDS input
#   - .update_gds_with_annotations()           in-place GDS -> aGDS conversion
#   - .write_agds()                            annotation-only aGDS (non-GDS input)
#   - .write_functional_annotation_subnodes()  shared writer: FunctionalAnnotation
#                                              as a STAARpipeline sub-node folder

#################### EXPORTED MAIN FUNCTIONS ####################

#' Annotate Variants with FAVOR Database Scores
#'
#' @description
#' Annotates variants with functional annotation scores from the FAVOR
#' (Functional Annotation of Variants Online Resource) database. This function
#' is used in PI estimation to retrieve pathogenicity-relevant features for
#' trait-associated and control variants.
#'
#' @param variants One of the following:
#'   \itemize{
#'     \item A data.frame with VarInfo column (CHR-POS-REF-ALT format)
#'     \item An S3 object of class \code{glow_pi_case_data} or
#'       \code{glow_pi_control_data} from \code{\link{prepare_PI_case_data}} or
#'       \code{\link{prepare_PI_control_data}}
#'     \item A path (character) to a GDS file containing variant data (requires
#'       SeqArray package). Variants are extracted automatically.
#'   }
#' @param favor_db_path Character. Path to directory containing FAVOR CSV chunk
#'   files (e.g., chr1_1.csv, chr1_2.csv, etc.)
#' @param favor_split_file Character or NULL. Path to `FAVORdatabase_chrsplit.csv`
#'   file that maps positions to chunk files. If NULL (default), searches in
#'   this order: (1) \code{favor_db_path}, (2) package bundled file. The bundled
#'   file ensures the function works even if the split file is not present in
#'   the FAVOR database directory
#' @param features Character vector. FAVOR feature column names to extract.
#'   Default is 11 standard features: apc_conservation, apc_protein_function_v3,
#'   apc_epigenetics_active, apc_epigenetics_repressed,
#'   apc_epigenetics_transcription, apc_local_nucleotide_diversity,
#'   apc_mappability, apc_transcription_factor, cadd_phred, linsight, fathmm_xf
#' @param output_csv Character or NULL. If provided, saves annotated results
#'   to CSV file at this path
#' @param output_agds Character or NULL. If provided, saves annotated results
#'   to aGDS file at this path (requires SeqArray package)
#' @param update_gds Logical. If TRUE and input is a GDS file, write annotations
#'   back to the input GDS file, converting it to an aGDS format. Default is
#'   FALSE. Ignored if variants is not a GDS file path.
#' @param variant_filter List or NULL. Optional filter to select a subset of
#'   variants. Works with all input types (GDS, data.frame, S3 objects).
#'   Supported filters:
#'   \itemize{
#'     \item \code{list(chr = "21")}: Filter by chromosome (works with all inputs)
#'     \item \code{list(chr = "21", start = 1e6, end = 5e6)}: Filter by position
#'       range (works with all inputs)
#'     \item \code{list(variant_ids = c(1, 5, 10))}: Filter by SeqArray variant
#'       IDs (GDS input only)
#'   }
#'   Default is NULL (use all variants). For non-GDS input, chromosome and
#'   position are parsed from the VarInfo column.
#' @param match_method Character. Matching strategy for variants:
#'   \itemize{
#'     \item "exact" (default): Exact VarInfo match only (fast, STAAR-compatible)
#'     \item "flexible": Hierarchical matching for strand flips and multiallelic
#'       variants. Priority order: (1) exact match, (2) same REF/different ALT,
#'       (3) swapped alleles, (4) swapped REF match, (5) position average
#'   }
#' @param na_allele_method Character. How to handle variants with NA alleles
#'   (CHR-POS-NA-NA format) when multiple FAVOR entries match by position:
#'   \itemize{
#'     \item "average" (default): Average all numeric annotation columns
#'     \item "first": Keep first matching entry
#'   }
#' @param na_handling Character. How to handle missing annotation values:
#'   \itemize{
#'     \item "keep" (default): Keep NA values as-is
#'     \item "zero": Replace NA with 0
#'     \item "drop": Remove variants with any NA annotations
#'   }
#' @param use_xsv Logical. If TRUE (default) and xsv CLI tool is available,
#'   use xsv for faster CSV joining (auto-falls back to R if unavailable)
#' @param verbose Integer. Verbosity level: 0=silent, 1=messages (default),
#'   2=detailed messages
#'
#' @return data.frame with VarInfo column and annotation feature columns.
#'   Order of variants is preserved from input.
#'
#' @details
#' \strong{FAVOR Database Structure:}
#'
#' The FAVOR database stores functional annotations for all possible variants
#' in the human genome, split into chromosome-specific chunks for efficient access:
#' \itemize{
#'   \item FAVORdatabase_chrsplit.csv: Maps genomic positions to chunk file numbers.
#'     A copy is bundled with GLOWr for convenience (used if not found in
#'     \code{favor_db_path})
#'   \item chr{N}_{K}.csv: Chunk files with annotations (variant_vcf key column)
#'   \item variant_vcf format: "CHR-POS-REF-ALT" (e.g., "1-12345-A-G")
#' }
#'
#' \strong{Annotation Process:}
#'
#' \enumerate{
#'   \item Parse VarInfo to extract CHR and POS
#'   \item Use split file to identify relevant FAVOR chunks
#'   \item Load chunks and join by VarInfo (exact match)
#'   \item Handle NA alleles: match by CHR-POS if VarInfo has NA-NA suffix
#'   \item Apply NA handling strategy (keep/zero/drop)
#'   \item Write to CSV and/or aGDS if requested
#' }
#'
#' \strong{Position-Only Matching (NA Alleles):}
#'
#' Literature-curated variant lists often lack REF/ALT information, resulting
#' in VarInfo like "1-12345-NA-NA". For these variants:
#' \itemize{
#'   \item Match by CHR and POS only (ignoring alleles)
#'   \item If multiple FAVOR entries at same position:
#'     \itemize{
#'       \item "average": Average numeric columns (more robust for multi-allelic)
#'       \item "first": Keep first match (faster, less comprehensive)
#'     }
#'   \item Warning emitted showing count of affected variants
#' }
#'
#' \strong{Default 11 FAVOR Features:}
#'
#' The default feature set balances comprehensiveness with computational
#' efficiency, covering conservation, protein function, epigenetics,
#' mappability, and variant effect prediction:
#' \enumerate{
#'   \item apc_conservation: Conservation scores
#'   \item apc_protein_function_v3: Protein functional impact (version 3)
#'   \item apc_epigenetics_active: Active chromatin marks
#'   \item apc_epigenetics_repressed: Repressed chromatin marks
#'   \item apc_epigenetics_transcription: Transcription-related marks
#'   \item apc_local_nucleotide_diversity: Local sequence diversity
#'   \item apc_mappability: Read mappability
#'   \item apc_transcription_factor: TF binding sites
#'   \item cadd_phred: CADD deleteriousness score
#'   \item linsight: Conservation fitness score
#'   \item fathmm_xf: Pathogenicity prediction
#' }
#'
#' \strong{Performance Optimization:}
#'
#' For large variant sets, the function uses several optimization strategies:
#' \itemize{
#'   \item Only loads FAVOR chunks containing relevant positions
#'   \item Uses data.table for efficient CSV reading and joining
#'   \item Optional xsv CLI tool for faster joins 
#'   \item Only extracts requested feature columns
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n log m) where n = number of variants, m = FAVOR entries per chunk.
#' Dominant operations: chunk loading O(k x m), joining O(n log m) per chunk,
#' where k = number of chunks accessed.
#'
#' @examples
#' \dontrun{
#' # Example 1: Annotate control variants (complete VarInfo)
#' controls <- prepare_PI_control_data(
#'   source = "data/reference_panel.gds",
#'   n_controls = 500
#' )
#' annotated_controls <- annotate_favor(
#'   variants = controls,
#'   favor_db_path = "data/FAVOR"
#' )
#'
#' # Example 2: Annotate case variants (may have NA alleles)
#' cases <- prepare_PI_case_data(
#'   data = "ALS_known_variants.xlsx",
#'   exclude_authors = "Nicolas A"
#' )
#' annotated_cases <- annotate_favor(
#'   variants = cases,
#'   favor_db_path = "data/FAVOR",
#'   na_allele_method = "average",  # Average multi-allelic matches
#'   verbose = 2
#' )
#'
#' # Example 3: Custom feature subset
#' annotated <- annotate_favor(
#'   variants = my_variants_df,
#'   favor_db_path = "data/FAVOR",
#'   features = c("apc_conservation", "cadd_phred", "linsight"),
#'   verbose = 1
#' )
#'
#' # Example 4: Save to CSV and aGDS
#' annotated <- annotate_favor(
#'   variants = cases,
#'   favor_db_path = "data/FAVOR",
#'   output_csv = "results/annotated_cases.csv",
#'   output_agds = "results/annotated_cases.agds"
#' )
#'
#' # Example 5: Handle missing annotations by dropping
#' annotated <- annotate_favor(
#'   variants = variants,
#'   favor_db_path = "data/FAVOR",
#'   na_handling = "drop",  # Remove variants with any NA
#'   verbose = 1
#' )
#'
#' # Example 6: Annotate variants from a GDS file
#' annotated <- annotate_favor(
#'   variants = "data/genotypes.gds",
#'   favor_db_path = "data/FAVOR",
#'   verbose = 1
#' )
#'
#' # Example 7: GDS input with variant filter (chromosome 21 only)
#' annotated <- annotate_favor(
#'   variants = "data/genotypes.gds",
#'   favor_db_path = "data/FAVOR",
#'   variant_filter = list(chr = "21"),
#'   verbose = 1
#' )
#'
#' # Example 8: GDS input with position range filter
#' annotated <- annotate_favor(
#'   variants = "data/genotypes.gds",
#'   favor_db_path = "data/FAVOR",
#'   variant_filter = list(chr = "21", start = 1e6, end = 5e6),
#'   verbose = 1
#' )
#'
#' # Example 9: Update GDS file in-place with annotations (convert to aGDS)
#' annotated <- annotate_favor(
#'   variants = "data/genotypes.gds",
#'   favor_db_path = "data/FAVOR",
#'   update_gds = TRUE,  # Write annotations back to input GDS
#'   verbose = 1
#' )
#'
#' # Example 10: Filter data.frame input by chromosome
#' annotated <- annotate_favor(
#'   variants = my_variants_df,  # data.frame with VarInfo column
#'   favor_db_path = "data/FAVOR",
#'   variant_filter = list(chr = "21"),  # Only annotate chr21 variants
#'   verbose = 1
#' )
#' }
#'
#' @references
#' Zhou, H., Arapoglou, T., Li, X., et al. (2023). FAVOR: functional annotation of
#' variants online resource and annotator for variation across the human genome.
#' Nucleic Acids Research, 51(D1), D1300-D1311. doi:10.1093/nar/gkac966
#'
#' @seealso
#' \code{\link{prepare_PI_case_data}} for case variant preparation
#' \code{\link{prepare_PI_control_data}} for control variant preparation
#' \code{\link{get_PI}} for PI estimation using annotated variants
#'
#' @export
annotate_favor <- function(
  variants,
  favor_db_path,
  favor_split_file = NULL,
  features = .default_favor_features(),
  output_csv = NULL,
  output_agds = NULL,
  update_gds = FALSE,
  variant_filter = NULL,
  match_method = "exact",
  na_allele_method = "average",
  na_handling = "keep",
  use_xsv = TRUE,
  verbose = 1
) {

  # ========== Step 1: Validate Inputs ==========

  if (verbose >= 1) {
    message("=== FAVOR Annotation ===")
  }

  # Validate match_method
  if (!match_method %in% c("exact", "flexible")) {
    stop("match_method must be 'exact' or 'flexible'")
  }

  # Validate na_allele_method
  if (!na_allele_method %in% c("average", "first")) {
    stop("na_allele_method must be 'average' or 'first'")
  }

  # Validate na_handling
  if (!na_handling %in% c("keep", "zero", "drop")) {
    stop("na_handling must be 'keep', 'zero', or 'drop'")
  }

  # Validate favor_db_path
  if (!dir.exists(favor_db_path)) {
    stop("FAVOR database directory not found: ", favor_db_path)
  }

  # ========== Step 2: Extract Data from GDS, S3 Object, or data.frame ==========

  # Track whether input is GDS (for update_gds option)
  gds_input_path <- NULL

  if (is.character(variants) && length(variants) == 1) {
    # Check if it's a file path
    if (file.exists(variants)) {
      # Check if GDS file
      if (grepl("\\.gds$", variants, ignore.case = TRUE)) {
        if (verbose >= 1) {
          message(sprintf("Extracting variants from GDS file: %s", basename(variants)))
        }
        gds_input_path <- variants
        variant_data <- .extract_varinfo_from_gds(
          gds_path = variants,
          variant_filter = variant_filter,
          verbose = verbose
        )
      } else {
        stop("File must be a .gds file or provide a data.frame")
      }
    } else {
      stop("GDS file not found: ", variants)
    }
  } else if (inherits(variants, "glow_pi_case_data") || inherits(variants, "glow_pi_control_data")) {
    if (verbose >= 2) {
      message(sprintf("Extracting variants from %s object...",
                      class(variants)[1]))
    }
    variant_data <- variants$data
  } else if (is.data.frame(variants)) {
    variant_data <- variants
  } else {
    stop("variants must be a data.frame, glow_pi_case_data/glow_pi_control_data object, or GDS file path")
  }

  # Apply variant_filter for non-GDS input (GDS filtering handled in .extract_varinfo_from_gds)
  if (!is.null(variant_filter) && is.null(gds_input_path)) {
    # Parse CHR from VarInfo (first component: CHR-POS-REF-ALT)
    if (!is.null(variant_filter$chr)) {
      parsed_chr <- sub("^([^-]+)-.*", "\\1", variant_data$VarInfo)
      chr_filter <- as.character(variant_filter$chr)
      variant_data <- variant_data[parsed_chr %in% chr_filter, , drop = FALSE]

      if (verbose >= 1) {
        message(sprintf("  Filtered to chromosome %s: %d variants",
                        paste(chr_filter, collapse = ","), nrow(variant_data)))
      }
    }

    # Position range filtering
    if (!is.null(variant_filter$start) && !is.null(variant_filter$end)) {
      parsed_pos <- as.integer(sub("^[^-]+-([0-9]+)-.*", "\\1", variant_data$VarInfo))
      in_range <- parsed_pos >= variant_filter$start & parsed_pos <= variant_filter$end
      variant_data <- variant_data[in_range, , drop = FALSE]

      if (verbose >= 1) {
        message(sprintf("  Filtered to position range %d-%d: %d variants",
                        variant_filter$start, variant_filter$end, nrow(variant_data)))
      }
    }

    # Check if any variants remain after filtering
    if (nrow(variant_data) == 0) {
      warning("No variants remain after applying variant_filter", call. = FALSE)
    }
  }

  # Validate VarInfo column exists
  if (!"VarInfo" %in% names(variant_data)) {
    stop("variants must have a 'VarInfo' column (CHR-POS-REF-ALT format)")
  }

  n_input <- nrow(variant_data)
  if (verbose >= 1) {
    message(sprintf("Input: %d variants", n_input))
  }

  # ========== Step 3: Load FAVOR Split File ==========
  # Priority order:
  #   1. User-provided favor_split_file (explicit)
  #   2. FAVORdatabase_chrsplit.csv in favor_db_path (if exists)
  #   3. Package bundled file in inst/extdata (fallback)

  if (is.null(favor_split_file)) {
    # Priority 1: Look in favor_db_path
    favor_split_file <- file.path(favor_db_path, "FAVORdatabase_chrsplit.csv")

    if (!file.exists(favor_split_file)) {
      # Priority 2: Use package bundled file
      favor_split_file <- system.file("extdata", "FAVORdatabase_chrsplit.csv",
                                       package = "GLOWr")
      if (!nzchar(favor_split_file) || !file.exists(favor_split_file)) {
        stop("Could not find FAVORdatabase_chrsplit.csv in ", favor_db_path,
             " or in package extdata")
      }
      if (verbose >= 1) {
        message("Using package bundled FAVORdatabase_chrsplit.csv")
      }
    } else {
      if (verbose >= 2) {
        message(sprintf("Using split file from FAVOR directory: %s", favor_split_file))
      }
    }
  } else {
    if (!file.exists(favor_split_file)) {
      stop("Split file not found: ", favor_split_file)
    }
    if (verbose >= 2) {
      message(sprintf("Using user-provided split file: %s", favor_split_file))
    }
  }

  if (verbose >= 1) {
    message("Loading FAVOR split file...")
  }

  split_data <- data.table::fread(favor_split_file, data.table = FALSE)

  # Validate split file structure
  required_cols <- c("Chr", "File_No", "Start_Pos", "End_Pos")
  if (!all(required_cols %in% names(split_data))) {
    stop(sprintf("Split file missing required columns: %s",
                 paste(setdiff(required_cols, names(split_data)), collapse = ", ")))
  }

  # ========== Step 4: Identify Relevant FAVOR Chunks ==========

  if (verbose >= 1) {
    message("Identifying relevant FAVOR chunks...")
  }

  chunks_needed <- .identify_favor_chunks(variant_data, split_data, verbose = verbose)

  if (length(chunks_needed) == 0) {
    warning("No FAVOR chunks matched variant positions. Returning input with NA annotations.",
            call. = FALSE)
    # Add NA columns for requested features
    for (feat in features) {
      variant_data[[feat]] <- NA_real_
    }
    return(variant_data)
  }

  if (verbose >= 1) {
    message(sprintf("  Will load %d FAVOR chunk file(s)", length(chunks_needed)))
  }

  # ========== Step 5: Annotate Variants ==========

  if (verbose >= 1) {
    message("Annotating variants from FAVOR database...")
  }

  # Determine whether to use xsv for joining
  # xsv-based join is only used for exact matching of complete VarInfo (no NA alleles, no flexible matching)
  # Flexible matching and NA allele handling require R-based processing
  has_na_alleles <- any(grepl("-NA-NA$", variant_data$VarInfo))
  xsv_available <- .check_xsv_available()
  can_use_xsv_join <- use_xsv &&
                      match_method == "exact" &&
                      !has_na_alleles &&
                      xsv_available

  # Log processing method
  if (verbose >= 1) {
    if (can_use_xsv_join) {
      message("  Method: xsv (exact matching)")
    } else if (match_method == "flexible") {
      message("  Method: R (flexible matching requires allele comparison)")
    } else if (has_na_alleles) {
      message("  Method: R (NA allele handling)")
    } else if (!xsv_available) {
      message("  Method: R (xsv not available)")
    } else {
      message("  Method: R")
    }
  }

  if (can_use_xsv_join) {
    # xsv path: direct join for exact matching
    annotated <- .join_favor_xsv(
      variant_data = variant_data,
      favor_db_path = favor_db_path,
      chunks_needed = chunks_needed,
      features = features,
      verbose = verbose
    )

  } else {
    # R path: data.table join with position filtering
    # Required for: flexible matching, NA allele handling, or when xsv unavailable
    annotated <- .join_favor_r(
      variant_data = variant_data,
      favor_db_path = favor_db_path,
      chunks_needed = chunks_needed,
      features = features,
      match_method = match_method,
      na_allele_method = na_allele_method,
      verbose = verbose
    )
  }

  # ========== Step 6: Handle Missing Annotations ==========

  if (verbose >= 1) {
    message("Handling missing annotations...")
  }

  annotated <- .handle_na_annotations(
    data = annotated,
    features = features,
    na_handling = na_handling,
    verbose = verbose
  )

  n_output <- nrow(annotated)
  if (verbose >= 1 && n_output < n_input) {
    message(sprintf("  %d variants removed due to NA handling", n_input - n_output))
  }

  # ========== Step 7: Write Outputs ==========

  if (!is.null(output_csv)) {
    if (verbose >= 1) {
      message(sprintf("Writing CSV output to: %s", output_csv))
    }
    data.table::fwrite(annotated, output_csv)
  }

  if (!is.null(output_agds)) {
    if (!is.null(gds_input_path)) {
      # GDS input -> create proper SeqArray aGDS by copying + annotating
      if (verbose >= 1) {
        message(sprintf("Creating aGDS output: %s", output_agds))
      }
      .create_agds_from_gds(
        input_gds_path = gds_input_path,
        output_path = output_agds,
        annotations = annotated,
        features = features,
        verbose = verbose
      )
    } else {
      # Non-GDS input -> annotation-only GDS (existing behavior)
      if (verbose >= 1) {
        message(sprintf("Writing annotation-only GDS to: %s", output_agds))
      }
      .write_agds(annotated, output_agds, features = features, verbose = verbose)
    }
  }

  # ========== Step 8: Update Input GDS (if requested) ==========

  if (update_gds && !is.null(gds_input_path)) {
    if (verbose >= 1) {
      message(sprintf("Updating GDS file with annotations: %s", basename(gds_input_path)))
    }
    .update_gds_with_annotations(
      gds_path = gds_input_path,
      annotations = annotated,
      features = features,
      verbose = verbose
    )
  } else if (update_gds && is.null(gds_input_path)) {
    warning("update_gds=TRUE ignored: input was not a GDS file", call. = FALSE)
  }

  # ========== Summary ==========

  if (verbose >= 1) {
    message(sprintf("\n=== Annotation Complete ==="))
    message(sprintf("  Input:  %d variants", n_input))
    message(sprintf("  Output: %d variants", n_output))
    message(sprintf("  Features: %d (%s)",
                    length(features),
                    paste(features[1:min(3, length(features))], collapse = ", ")))
    if (length(features) > 3) {
      message(sprintf("            ... and %d more", length(features) - 3))
    }
  }

  return(annotated)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Extract VarInfo from GDS File
#'
#' @description
#' Extracts variant information from a SeqArray GDS file and creates a
#' data.frame with VarInfo column in CHR-POS-REF-ALT format.
#'
#' @param gds_path Character. Path to GDS file
#' @param variant_filter List or NULL. Optional filter:
#'   \itemize{
#'     \item \code{list(chr = "21")}: Filter by chromosome
#'     \item \code{list(chr = "21", start = 1e6, end = 5e6)}: Filter by position
#'     \item \code{list(variant_ids = c(1, 5, 10))}: Filter by variant IDs
#'   }
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with VarInfo column and variant.id for tracking
#'
#' @details
#' Requires SeqArray package. Opens GDS file, applies optional filter,
#' extracts chromosome, position, ref, alt, and creates VarInfo string.
#'
#' @keywords internal
#' @noRd
.extract_varinfo_from_gds <- function(gds_path, variant_filter = NULL, verbose = 1) {

  # Check if SeqArray is available

  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required for GDS input. Install with: BiocManager::install('SeqArray')")
  }

  # Open GDS file
  gds <- SeqArray::seqOpen(gds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Apply variant filter if provided
  if (!is.null(variant_filter)) {
    .apply_variant_filter(gds, variant_filter, verbose = verbose)
  }

  # Extract variant information
  chr <- SeqArray::seqGetData(gds, "chromosome")
  pos <- SeqArray::seqGetData(gds, "position")

  # Handle ref and alt - these can be complex in multi-allelic cases
  ref <- SeqArray::seqGetData(gds, "$ref")
  alt <- SeqArray::seqGetData(gds, "$alt")

  # Get variant IDs for tracking
  variant_id <- SeqArray::seqGetData(gds, "variant.id")

  # Create VarInfo: CHR-POS-REF-ALT format
  # For multi-allelic sites, alt may be comma-separated; take first ALT
  alt_first <- sapply(strsplit(alt, ",", fixed = TRUE), function(x) x[1])

  # Create VarInfo string
  VarInfo <- paste(chr, pos, ref, alt_first, sep = "-")

  # Create result data frame
  result <- data.frame(
    VarInfo = VarInfo,
    # chr = chr, #ZWu: chr, pos, ref, alt_first could be added to data to avoid parsing them later and speed up the process?  To be tested.
    # pos = pos,
    # ref = ref,
    # alt = alt_first,
    variant_id = variant_id,
    stringsAsFactors = FALSE
  )

  if (verbose >= 2) {
    message(sprintf("  Extracted %d variants from GDS", nrow(result)))
  }

  return(result)
}


#' Apply Variant Filter to GDS File
#'
#' @description
#' Applies a filter to select a subset of variants from a GDS file.
#'
#' @param gds SeqArray GDS object (opened)
#' @param variant_filter List with filter specifications
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (modifies GDS filter in place)
#'
#' @keywords internal
#' @noRd
.apply_variant_filter <- function(gds, variant_filter, verbose = 1) {

  if (is.null(variant_filter)) return(invisible(NULL))

  # Filter by variant IDs
  if (!is.null(variant_filter$variant_ids)) {
    SeqArray::seqSetFilter(gds, variant.id = variant_filter$variant_ids)
    if (verbose >= 2) {
      message(sprintf("  Applied variant ID filter: %d variants selected",
                      length(variant_filter$variant_ids)))
    }
    return(invisible(NULL))
  }

  # Filter by chromosome
  if (!is.null(variant_filter$chr)) {
    # Get chromosome filter
    chr_filter <- variant_filter$chr

    # Apply chromosome filter
    SeqArray::seqSetFilterChrom(gds, chr_filter)

    if (verbose >= 2) {
      n_filtered <- sum(SeqArray::seqGetFilter(gds)$variant.sel)
      message(sprintf("  Applied chromosome filter (chr=%s): %d variants selected",
                      chr_filter, n_filtered))
    }

    # Apply position range filter if provided
    if (!is.null(variant_filter$start) && !is.null(variant_filter$end)) {
      pos <- SeqArray::seqGetData(gds, "position")
      in_range <- pos >= variant_filter$start & pos <= variant_filter$end
      SeqArray::seqSetFilter(gds, variant.sel = in_range, action = "intersect")

      if (verbose >= 2) {
        n_filtered <- sum(SeqArray::seqGetFilter(gds)$variant.sel)
        message(sprintf("  Applied position filter (%d-%d): %d variants selected",
                        variant_filter$start, variant_filter$end, n_filtered))
      }
    }
  }

  return(invisible(NULL))
}


#' Update GDS File with FAVOR Annotations
#'
#' @description
#' Writes annotation data back to an existing GDS file, converting it to
#' aGDS format by adding a FunctionalAnnotation node.
#'
#' @param gds_path Character. Path to GDS file (will be modified)
#' @param annotations data.frame with VarInfo and annotation feature columns
#' @param features Character vector of feature column names
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (side effect: modifies GDS file)
#'
#' @details
#' Opens the GDS file in read-write mode and adds
#' \code{/annotation/info/FunctionalAnnotation} as a \emph{folder}
#' (\code{addfolder.gdsn}) holding one native-typed sub-node per feature
#' (numeric features stay numeric, string features stay character). This is the
#' STAARpipeline sub-node layout read by the per-feature variant-set scan via
#' \code{seqGetData(gds, ".../FunctionalAnnotation/<feature>")}; it is NOT a
#' single numeric matrix. A \code{feature_names} attribute is also placed on the
#' folder node for discoverability, but the sub-node names are the source of
#' truth.
#'
#' If annotation nodes already exist, they are overwritten with a warning.
#'
#' @keywords internal
#' @noRd
.update_gds_with_annotations <- function(gds_path, annotations, features, verbose = 1) {

  # Check if gdsfmt is available
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required to update GDS files. Install with: BiocManager::install('gdsfmt')")
  }

  # Open GDS file for read-write
  gds <- gdsfmt::openfn.gds(gds_path, readonly = FALSE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)

  # Check if annotation folder exists, create if not
  annotation_exists <- "annotation" %in% gdsfmt::ls.gdsn(gds)
  if (!annotation_exists) {
    annot_folder <- gdsfmt::addfolder.gdsn(gds, "annotation")
  } else {
    annot_folder <- gdsfmt::index.gdsn(gds, "annotation")
  }

  # Check if info folder exists, create if not
  info_exists <- "info" %in% gdsfmt::ls.gdsn(annot_folder)
  if (!info_exists) {
    info_folder <- gdsfmt::addfolder.gdsn(annot_folder, "info")
  } else {
    info_folder <- gdsfmt::index.gdsn(annot_folder, "info")
  }

  # Write FunctionalAnnotation as a STAARpipeline-style sub-node folder (one
  # typed sub-node per feature). The shared helper handles overwrite-with-warning
  # of any pre-existing node and preserves each feature's native type.
  .write_functional_annotation_subnodes(
    parent_node = info_folder,
    annotations = annotations,
    features = features,
    verbose = verbose
  )

  if (verbose >= 1) {
    message(sprintf("  Added FunctionalAnnotation node: %d variants x %d features",
                    nrow(annotations), length(features)))
  }

  return(invisible(NULL))
}


#' Create aGDS by Copying Input GDS and Adding Annotations
#'
#' @description
#' Creates a STAARpipeline-compatible aGDS file by copying the input SeqArray
#' GDS file and adding annotation data. The result is a valid SeqArray file
#' that can be opened with \code{seqOpen()}.
#'
#' @param input_gds_path Character. Path to input SeqArray GDS file
#' @param output_path Character. Path for output aGDS file
#' @param annotations data.frame with VarInfo and annotation feature columns
#' @param features Character vector of feature column names
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (side effect: creates aGDS file)
#'
#' @details
#' Follows STAARpipeline's gds2agds.R pattern:
#' \enumerate{
#'   \item Copy input GDS file to output path
#'   \item Open copy with seqOpen(readonly = FALSE)
#'   \item Add /annotation/info/FunctionalAnnotation node
#'   \item Close file
#' }
#'
#' \strong{Annotation node format:} \code{/annotation/info/FunctionalAnnotation}
#' is created as a \emph{folder} (\code{addfolder.gdsn}) holding one sub-node per
#' feature, each carrying that feature's native-typed vector in GDS variant order
#' (numeric features stay numeric, string features such as
#' \code{genecode_comprehensive_category} stay character). This is the
#' STAARpipeline sub-node layout that the per-feature variant-set scan reads via
#' \code{seqGetData(gds, ".../FunctionalAnnotation/<feature>")}. It is NOT a
#' single numeric matrix, so mixed numeric/string feature sets are carried
#' without coercion.
#'
#' This approach preserves the original input GDS file while creating a new
#' aGDS file with both genotype data and functional annotations.
#'
#' @keywords internal
#' @noRd
.create_agds_from_gds <- function(input_gds_path, output_path, annotations,
                                   features, verbose = 1) {

  # Check required packages
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required for aGDS output. Install with: BiocManager::install('SeqArray')")
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required for aGDS output. Install with: BiocManager::install('gdsfmt')")
  }

  # Step 1: Copy input GDS to output location
  if (verbose >= 1) {
    message(sprintf("  Copying GDS file to: %s", basename(output_path)))
  }

  # Ensure output directory exists
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Copy the file
  copy_success <- file.copy(input_gds_path, output_path, overwrite = TRUE)
  if (!copy_success) {
    stop("Failed to copy GDS file to: ", output_path)
  }

  # Step 2: Open the copy for writing
  gds <- SeqArray::seqOpen(output_path, readonly = FALSE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Step 3: Get VarInfo from the copied GDS to match annotation order
  # Extract directly without closing/reopening
  chr <- SeqArray::seqGetData(gds, "chromosome")
  pos <- SeqArray::seqGetData(gds, "position")
  ref <- SeqArray::seqGetData(gds, "$ref")
  alt <- SeqArray::seqGetData(gds, "$alt")

  # For multi-allelic sites, take first ALT
  alt_first <- sapply(strsplit(alt, ",", fixed = TRUE), function(x) x[1])

  # Create VarInfo string
  gds_varinfo <- paste(chr, pos, ref, alt_first, sep = "-")

  # Match annotations to GDS variant order
  annot_idx <- match(gds_varinfo, annotations$VarInfo)

  # Reorder the annotation table into GDS variant order. Reordering the
  # data.frame (not coercing to a matrix) preserves each feature's native type,
  # so string coding nodes (e.g. genecode_comprehensive_category) stay character.
  annot_ordered <- annotations[annot_idx, , drop = FALSE]

  # Step 4: Add FunctionalAnnotation as a STAARpipeline-style sub-node folder
  # (one typed sub-node per feature), so the per-feature scan reads it via
  # seqGetData(gds, ".../FunctionalAnnotation/<feature>").
  # Navigate to or create annotation/info folder
  anno_folder <- tryCatch(
    gdsfmt::index.gdsn(gds, "annotation/info"),
    error = function(e) NULL
  )

  if (is.null(anno_folder)) {
    # Create annotation/info folder structure
    anno_root <- tryCatch(
      gdsfmt::index.gdsn(gds, "annotation"),
      error = function(e) NULL
    )
    if (is.null(anno_root)) {
      anno_root <- gdsfmt::addfolder.gdsn(gds, "annotation")
    }
    anno_folder <- gdsfmt::addfolder.gdsn(anno_root, "info")
  }

  fa_node <- .write_functional_annotation_subnodes(
    parent_node = anno_folder,
    annotations = annot_ordered,
    features = features,
    verbose = verbose
  )

  if (verbose >= 1) {
    n_annotated <- sum(!is.na(annot_idx))
    message(sprintf("  Created aGDS: %d variants, %d annotated, %d features",
                    length(gds_varinfo), n_annotated, length(features)))
  }

  return(invisible(NULL))
}


#' Default FAVOR Annotation Features
#'
#' @description
#' Returns all 20 numerical annotation features from the FAVOR Essential DB:
#' 17 Annotation Principal Components (all versions and sub-features) +
#' 3 integrative scores. Used as the default for \code{annotate_favor()} to
#' annotate data broadly for downstream flexibility.
#'
#' For PI model training, use \code{.default_PI_features()} from
#' \code{get_PI_train.R}, which is a curated 16-feature subset.
#'
#' @return Character vector of 20 feature names.
#'
#' @seealso \code{.default_PI_features()} for the curated PI training set.
#'
#' @keywords internal
#' @noRd
.default_favor_features <- function() {
  c(
    # 17 Annotation Principal Components (all versions and sub-features)
    "apc_conservation",
    "apc_conservation_v2",
    "apc_epigenetics",
    "apc_epigenetics_active",
    "apc_epigenetics_repressed",
    "apc_epigenetics_transcription",
    "apc_local_nucleotide_diversity",
    "apc_local_nucleotide_diversity_v2",
    "apc_local_nucleotide_diversity_v3",
    "apc_mappability",
    "apc_micro_rna",
    "apc_mutation_density",
    "apc_protein_function_v3",
    "apc_proximity_to_coding",
    "apc_proximity_to_coding_v2",
    "apc_proximity_to_tsstes",
    "apc_transcription_factor",
    # 3 Integrative Scores
    "cadd_phred",
    "linsight",
    "fathmm_xf"
  )
}


#' Identify FAVOR Chunk Files Needed for Variants
#'
#' @description
#' Maps variants to FAVOR chunk files using the split file, based on
#' chromosome and position ranges.
#'
#' @param variant_data data.frame with VarInfo column
#' @param split_data data.frame with Chr, File_No, Start_Pos, End_Pos columns
#' @param verbose Integer. Verbosity level
#'
#' @return Character vector of chunk file names (e.g., "chr1_1.csv")
#'
#' @details
#' Parses VarInfo to extract CHR and POS, then matches against split file
#' position ranges to determine which chunk files contain annotations for
#' the input variants.
#'
#' @keywords internal
#' @noRd
.identify_favor_chunks <- function(variant_data, split_data, verbose = 1) {

  # Parse VarInfo to extract CHR and POS
  varinfo_split <- strsplit(variant_data$VarInfo, "-", fixed = TRUE)

  # Extract CHR (first element)
  chr_vec <- sapply(varinfo_split, function(x) {
    if (length(x) >= 1) x[1] else NA_character_
  })

  # Extract POS (second element)
  pos_vec <- as.numeric(sapply(varinfo_split, function(x) {
    if (length(x) >= 2) x[2] else NA_character_
  }))

  # Remove variants with NA CHR or POS
  valid_idx <- !is.na(chr_vec) & !is.na(pos_vec)
  if (sum(valid_idx) == 0) {
    warning("No variants with valid CHR-POS in VarInfo", call. = FALSE)
    return(character(0))
  }

  chr_vec <- chr_vec[valid_idx]
  pos_vec <- pos_vec[valid_idx]

  # Get unique CHR-POS ranges
  unique_chr <- unique(chr_vec)

  # Find matching chunks
  chunks_needed <- character(0)

  for (chr in unique_chr) {
    # Get positions for this chromosome
    chr_idx <- chr_vec == chr
    chr_positions <- pos_vec[chr_idx]

    # Get range
    min_pos <- min(chr_positions, na.rm = TRUE)
    max_pos <- max(chr_positions, na.rm = TRUE)

    # Find chunks that overlap this range
    chr_chunks <- split_data[split_data$Chr == chr, ]

    # A chunk overlaps if:
    # - chunk_start <= max_pos AND chunk_end >= min_pos
    overlap_idx <- (chr_chunks$Start_Pos <= max_pos) & (chr_chunks$End_Pos >= min_pos)
    overlapping_chunks <- chr_chunks[overlap_idx, ]

    if (nrow(overlapping_chunks) > 0) {
      # Build chunk file names
      chunk_names <- paste0("chr", overlapping_chunks$Chr, "_",
                            overlapping_chunks$File_No, ".csv")
      chunks_needed <- c(chunks_needed, chunk_names)
    }
  }

  chunks_needed <- unique(chunks_needed)

  if (verbose >= 2) {
    message(sprintf("    Chromosomes: %s", paste(unique_chr, collapse = ", ")))
    message(sprintf("    Chunks: %s", paste(chunks_needed, collapse = ", ")))
  }

  return(chunks_needed)
}


#' Filter FAVOR Chunk by Position
#'
#' @description
#' Filters a FAVOR chunk file to only rows matching specified positions.
#' Uses R's data.table with hash-based \code{\%in\%} for efficient O(n+m)
#' filtering.
#'
#' @param chunk_path Path to FAVOR chunk CSV file
#' @param positions Integer vector of positions to filter for
#' @param cols Character vector of columns to load
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with filtered rows and requested columns
#'
#' @details
#' For flexible matching, we only need FAVOR rows at positions present in the
#' query variants. This function filters during loading rather than after,
#' reducing memory footprint for subsequent R operations.
#'
#' Note: xsv-based filtering was removed in 2026-01 due to O(N) regex
#' complexity causing catastrophic performance for large position sets.
#' R's hash-based \code{\%in\%} provides O(1) lookup per row regardless
#' of position count.
#'
#' @keywords internal
#' @noRd
.filter_favor_by_position <- function(chunk_path, positions, cols, verbose = 1) {

  # Load chunk and filter by position using hash-based %in%
  if (verbose >= 2) {
    message("    Loading and filtering by position...")
  }

  chunk <- data.table::fread(
    chunk_path,
    select = cols,
    data.table = FALSE,
    showProgress = FALSE
  )

  # Filter by position (hash-based, O(1) per lookup)
  result <- chunk[chunk$position %in% positions, , drop = FALSE]

  if (verbose >= 2) {
    message(sprintf("    Position filter: %d -> %d rows", nrow(chunk), nrow(result)))
  }

  return(result)
}


#' Join Variants with FAVOR Annotations using R/data.table
#'
#' @description
#' Loads FAVOR chunk files and joins with variant data to retrieve annotations.
#' Primary annotation method using data.table for efficiency.
#'
#' @param variant_data data.frame with VarInfo column
#' @param favor_db_path Character. Path to FAVOR directory
#' @param chunks_needed Character vector of chunk file names
#' @param features Character vector of feature column names to extract
#' @param na_allele_method Character. "average" or "first" for position matching
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with VarInfo and annotation feature columns
#'
#' @details
#' \strong{Join Strategy:}
#'
#' \enumerate{
#'   \item Load each FAVOR chunk (select only needed columns)
#'   \item Join by VarInfo exact match (CHR-POS-REF-ALT)
#'   \item For variants with NA alleles (CHR-POS-NA-NA), match by position only
#'   \item Combine results from all chunks
#' }
#'
#' Uses data.table::fread() for fast CSV loading and merge for joining.
#'
#' @keywords internal
#' @noRd
.join_favor_r <- function(
  variant_data,
  favor_db_path,
  chunks_needed,
  features,
  match_method = "exact",
  na_allele_method = "average",
  verbose = 1
) {

  # Initialize result with input variants
  # Convert to plain data.frame to ensure reliable element-wise assignment.
  # data.table's reference semantics can cause result[[col]][idx] <- value
  # to silently fail (modified copy not written back to the table).
  result <- as.data.frame(variant_data)

  # Add annotation columns (initially NA)
  for (feat in features) {
    result[[feat]] <- NA_real_
  }

  # Separate variants with complete VarInfo vs NA alleles
  has_na_alleles <- grepl("-NA-NA$", result$VarInfo)
  complete_varinfo <- result[!has_na_alleles, , drop = FALSE]
  na_allele_varinfo <- result[has_na_alleles, , drop = FALSE]

  if (verbose >= 2) {
    message(sprintf("    %d variants with complete VarInfo", nrow(complete_varinfo)))
    message(sprintf("    %d variants with NA alleles", nrow(na_allele_varinfo)))
  }

  # Track which variants were annotated
  annotated_idx <- rep(FALSE, nrow(result))

  # Process each chunk
  for (chunk_file in chunks_needed) {
    chunk_path <- file.path(favor_db_path, chunk_file)

    if (!file.exists(chunk_path)) {
      if (verbose >= 1) {
        warning(sprintf("Chunk file not found, skipping: %s", chunk_file),
                call. = FALSE)
      }
      next
    }

    if (verbose >= 2) {
      message(sprintf("    Loading chunk: %s", chunk_file))
    }

    # Define columns to load
    # Include chromosome, position, ref_vcf, alt_vcf for flexible matching optimization
    # (these columns already exist in FAVOR - avoids parsing variant_vcf)
    cols_to_load <- c("variant_vcf", "chromosome", "position", "ref_vcf", "alt_vcf", features)

    # Load chunk with position filtering
    # Reduces memory footprint before R-based join operations
    query_positions <- unique(as.integer(sub("^[^-]+-([0-9]+)-.*", "\\1", variant_data$VarInfo)))
    chunk_data <- .filter_favor_by_position(
      chunk_path = chunk_path,
      positions = query_positions,
      cols = cols_to_load,
      verbose = verbose
    )

    # Join with complete VarInfo variants
    if (nrow(complete_varinfo) > 0) {
      # Step 1: Try exact match first
      matched <- merge(
        complete_varinfo[, "VarInfo", drop = FALSE],
        chunk_data,
        by.x = "VarInfo",
        by.y = "variant_vcf",
        all.x = FALSE,  # Only keep matches
        sort = FALSE
      )

      if (nrow(matched) > 0) {
        # Update result with annotations
        for (i in seq_len(nrow(matched))) {
          varinfo <- matched$VarInfo[i]
          result_idx <- which(result$VarInfo == varinfo)[1]
          if (!is.na(result_idx) && !annotated_idx[result_idx]) {
            has_any_annotation <- FALSE
            for (feat in features) {
              if (feat %in% names(matched) && !is.na(matched[[feat]][i])) {
                result[[feat]][result_idx] <- matched[[feat]][i]
                has_any_annotation <- TRUE
              }
            }
            # Only mark as annotated if at least one feature was non-NA
            if (has_any_annotation) {
              annotated_idx[result_idx] <- TRUE
            }
          }
        }
      }

      # Step 2: If flexible matching, try hierarchical matching for unmatched variants
      if (match_method == "flexible") {
        # Find variants not yet annotated
        unannotated_mask <- !has_na_alleles & !annotated_idx
        unannotated <- result[unannotated_mask, , drop = FALSE]

        if (nrow(unannotated) > 0) {
          # chunk_data is already position-filtered during loading (via .filter_favor_by_position)
          # so no additional filtering needed here

          matched_flex <- .match_flexible(
            variants = unannotated,
            favor_data = chunk_data,
            features = features,
            verbose = verbose
          )

          if (nrow(matched_flex) > 0) {
            # Update result with flexible match annotations
            for (i in seq_len(nrow(matched_flex))) {
              varinfo <- matched_flex$VarInfo[i]
              result_idx <- which(result$VarInfo == varinfo)[1]
              if (!is.na(result_idx) && !annotated_idx[result_idx]) {
                has_any_annotation <- FALSE
                for (feat in features) {
                  if (feat %in% names(matched_flex) && !is.na(matched_flex[[feat]][i])) {
                    result[[feat]][result_idx] <- matched_flex[[feat]][i]
                    has_any_annotation <- TRUE
                  }
                }
                # Only mark as annotated if at least one feature was non-NA
                if (has_any_annotation) {
                  annotated_idx[result_idx] <- TRUE
                }
              }
            }
          }
        }
      }
    }

    # Join with NA allele variants (position matching)
    if (nrow(na_allele_varinfo) > 0) {
      matched_na <- .match_by_position(
        na_allele_varinfo,
        chunk_data,
        features = features,
        method = na_allele_method,
        verbose = verbose
      )

      if (nrow(matched_na) > 0) {
        # Update result with annotations
        for (i in seq_len(nrow(matched_na))) {
          varinfo <- matched_na$VarInfo[i]
          result_idx <- which(result$VarInfo == varinfo)[1]
          if (!is.na(result_idx) && !annotated_idx[result_idx]) {
            has_any_annotation <- FALSE
            for (feat in features) {
              if (feat %in% names(matched_na) && !is.na(matched_na[[feat]][i])) {
                result[[feat]][result_idx] <- matched_na[[feat]][i]
                has_any_annotation <- TRUE
              }
            }
            # Only mark as annotated if at least one feature was non-NA
            if (has_any_annotation) {
              annotated_idx[result_idx] <- TRUE
            }
          }
        }
      }
    }
  }

  if (verbose >= 1) {
    n_annotated <- sum(annotated_idx)
    message(sprintf("  Annotated %d/%d variants (%.1f%%)",
                    n_annotated, nrow(result), 100 * n_annotated / nrow(result)))
  }

  return(result)
}


#' Match Variants by Position Only (for NA Alleles)
#'
#' @description
#' Handles variants with CHR-POS-NA-NA format by matching on position only.
#' If multiple FAVOR entries exist at the same position, either averages
#' numeric columns or takes the first match.
#'
#' @param variants_na data.frame with VarInfo containing NA-NA alleles
#' @param favor_data data.frame from FAVOR chunk (variant_vcf column)
#' @param features Character vector of feature columns
#' @param method Character. "average" or "first"
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with VarInfo and annotation columns
#'
#' @details
#' Creates position key (CHR-POS) from VarInfo and variant_vcf, then joins
#' by position. For multi-allelic sites with multiple matches:
#' \itemize{
#'   \item "average": Average all numeric feature columns across matches
#'   \item "first": Keep only the first match
#' }
#'
#' @keywords internal
#' @noRd
.match_by_position <- function(variants_na, favor_data, features, method, verbose) {

  if (nrow(variants_na) == 0) {
    return(variants_na)
  }

  # Create position key for variants (CHR-POS)
  variant_pos_keys <- sub("-NA-NA$", "", variants_na$VarInfo)

  # Parse FAVOR variant_vcf to extract CHR and POS more efficiently
  # Convert favor_data to data.table if not already
  if (!data.table::is.data.table(favor_data)) {
    favor_dt <- data.table::as.data.table(favor_data)
  } else {
    favor_dt <- data.table::copy(favor_data)
  }

  # Extract CHR and POS from variant_vcf using regex (faster than split+paste)
  # Format: CHR-POS-REF-ALT -> CHR-POS
  # Use set() instead of := to avoid namespace issues
  data.table::set(favor_dt, j = "pos_key",
                  value = sub("^([^-]+-[^-]+)-.*$", "\\1", favor_dt$variant_vcf))

  # Set key for fast lookup
  data.table::setkey(favor_dt, pos_key)

  # Initialize result list
  result_list <- vector("list", nrow(variants_na))
  n_multi <- 0

  # Process each variant
  for (i in seq_len(nrow(variants_na))) {
    varinfo <- variants_na$VarInfo[i]
    pos_key <- variant_pos_keys[i]

    # Fast lookup using data.table key
    matched_rows <- favor_dt[pos_key, nomatch = NULL]

    if (nrow(matched_rows) == 0) {
      # No match - create row with NAs
      result_row <- data.frame(VarInfo = varinfo, stringsAsFactors = FALSE)
      for (feat in features) {
        result_row[[feat]] <- NA_real_
      }
      result_list[[i]] <- result_row
      next
    }

    # Create result row
    result_row <- data.frame(VarInfo = varinfo, stringsAsFactors = FALSE)

    if (nrow(matched_rows) == 1) {
      # Single match
      for (feat in features) {
        if (feat %in% names(matched_rows)) {
          result_row[[feat]] <- matched_rows[[feat]][1]
        } else {
          result_row[[feat]] <- NA_real_
        }
      }
    } else {
      # Multiple matches
      n_multi <- n_multi + 1

      if (method == "average") {
        for (feat in features) {
          if (feat %in% names(matched_rows)) {
            values <- matched_rows[[feat]]
            if (is.numeric(values)) {
              result_row[[feat]] <- mean(values, na.rm = TRUE)
            } else {
              result_row[[feat]] <- values[1]
            }
          } else {
            result_row[[feat]] <- NA_real_
          }
        }
      } else {
        # Keep first
        for (feat in features) {
          if (feat %in% names(matched_rows)) {
            result_row[[feat]] <- matched_rows[[feat]][1]
          } else {
            result_row[[feat]] <- NA_real_
          }
        }
      }
    }

    result_list[[i]] <- result_row
  }

  # Combine all results
  result <- data.table::rbindlist(result_list, fill = TRUE)
  result <- as.data.frame(result)

  if (n_multi > 0 && verbose >= 1) {
    action <- ifelse(method == "average", "averaged", "kept first")
    message(sprintf("    Position matching: %d variants matched multiple FAVOR entries (%s)",
                    n_multi, action))
  }

  # Remove row names
  rownames(result) <- NULL

  return(result)
}


#' Flexible Hierarchical Matching for Variants
#'
#' @description
#' Implements hierarchical matching for variants with complete alleles that
#' didn't get exact matches. Handles strand flips and multiallelic sites.
#'
#' @param variants data.frame with VarInfo column (CHR-POS-REF-ALT format)
#' @param favor_data data.frame from FAVOR chunk (variant_vcf column)
#' @param features Character vector of feature columns
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with VarInfo and annotation columns for matched variants
#'
#' @details
#' Matching priority hierarchy:
#' \enumerate{
#'   \item Same REF (multiallelic): Same CHR-POS-REF, different ALT
#'   \item Swapped exact: CHR-POS with ALT-REF swapped
#'   \item Swapped REF (multiallelic): FAVOR's ALT matches input's REF
#'   \item Position average: Average all entries at CHR-POS
#' }
#'
#' Within each priority level:
#' \itemize{
#'   \item Single match: return it
#'   \item Multiple matches: average numeric columns
#' }
#'
#' @keywords internal
#' @noRd
.match_flexible <- function(variants, favor_data, features, verbose = 1) {

  if (nrow(variants) == 0 || nrow(favor_data) == 0) {
    # Return empty result with expected columns
    result <- variants[0, "VarInfo", drop = FALSE]
    for (feat in features) result[[feat]] <- numeric(0)
    return(result)
  }

  # Get chr/pos/ref/alt - use existing columns if available, else parse VarInfo
  # (Parsing 500 strings is fast <0.01s; this is just for cleaner code)
  if (all(c("CHR", "POS", "REF", "ALT") %in% names(variants))) {
    variants$chr <- as.character(variants$CHR)
    variants$pos <- as.character(variants$POS)
    variants$ref <- variants$REF
    variants$alt <- variants$ALT
  } else if (all(c("chr", "pos", "ref", "alt") %in% names(variants))) {
    variants$chr <- as.character(variants$chr)
    variants$pos <- as.character(variants$pos)
    # ref/alt already exist
  } else {
    # Fallback: parse VarInfo
    parsed <- do.call(rbind, strsplit(variants$VarInfo, "-", fixed = TRUE))
    variants$chr <- parsed[, 1]
    variants$pos <- parsed[, 2]
    variants$ref <- parsed[, 3]
    variants$alt <- parsed[, 4]
  }

  # Use existing FAVOR columns (NO parsing of variant_vcf)
  favor_data$chr <- as.character(favor_data$chromosome)
  favor_data$pos <- as.character(favor_data$position)
  favor_data$ref <- favor_data$ref_vcf
  favor_data$alt <- favor_data$alt_vcf

  # Initialize result
  result <- data.frame(VarInfo = variants$VarInfo, stringsAsFactors = FALSE)
  for (feat in features) result[[feat]] <- NA_real_

  matched_mask <- rep(FALSE, nrow(variants))
  n_same_ref <- 0
  n_swapped <- 0
  n_pos_avg <- 0

  # === Priority 2: Same REF (multiallelic) - vectorized ===
  unmatched <- variants[!matched_mask, , drop = FALSE]
  if (nrow(unmatched) > 0) {
    # Merge on chr-pos-ref
    merged <- merge(
      unmatched[, c("VarInfo", "chr", "pos", "ref")],
      favor_data[, c("chr", "pos", "ref", features), drop = FALSE],
      by = c("chr", "pos", "ref"),
      all.x = FALSE
    )
    if (nrow(merged) > 0) {
      # Aggregate multiple matches by averaging
      agg <- aggregate(merged[, features, drop = FALSE],
                       by = list(VarInfo = merged$VarInfo),
                       FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
      # Update result
      for (i in seq_len(nrow(agg))) {
        idx <- which(result$VarInfo == agg$VarInfo[i])
        if (length(idx) == 1 && !matched_mask[idx]) {
          for (feat in features) result[[feat]][idx] <- agg[[feat]][i]
          matched_mask[idx] <- TRUE
          n_same_ref <- n_same_ref + 1
        }
      }
    }
  }

  # === Priority 3: Swapped alleles (REF<->ALT) - vectorized ===
  unmatched <- variants[!matched_mask, , drop = FALSE]
  if (nrow(unmatched) > 0) {
    # Look for FAVOR entries where FAVOR.ref = variant.alt and FAVOR.alt = variant.ref
    merged <- merge(
      unmatched[, c("VarInfo", "chr", "pos", "ref", "alt")],
      favor_data[, c("chr", "pos", "ref", "alt", features), drop = FALSE],
      by.x = c("chr", "pos", "alt", "ref"),  # variant alt matches favor ref, etc.
      by.y = c("chr", "pos", "ref", "alt"),
      all.x = FALSE
    )
    if (nrow(merged) > 0) {
      agg <- aggregate(merged[, features, drop = FALSE],
                       by = list(VarInfo = merged$VarInfo),
                       FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
      for (i in seq_len(nrow(agg))) {
        idx <- which(result$VarInfo == agg$VarInfo[i])
        if (length(idx) == 1 && !matched_mask[idx]) {
          for (feat in features) result[[feat]][idx] <- agg[[feat]][i]
          matched_mask[idx] <- TRUE
          n_swapped <- n_swapped + 1
        }
      }
    }
  }

  # === Priority 5: Position average - vectorized ===
  unmatched <- variants[!matched_mask, , drop = FALSE]
  if (nrow(unmatched) > 0) {
    merged <- merge(
      unmatched[, c("VarInfo", "chr", "pos")],
      favor_data[, c("chr", "pos", features), drop = FALSE],
      by = c("chr", "pos"),
      all.x = FALSE
    )
    if (nrow(merged) > 0) {
      agg <- aggregate(merged[, features, drop = FALSE],
                       by = list(VarInfo = merged$VarInfo),
                       FUN = function(x) mean(as.numeric(x), na.rm = TRUE))
      for (i in seq_len(nrow(agg))) {
        idx <- which(result$VarInfo == agg$VarInfo[i])
        if (length(idx) == 1 && !matched_mask[idx]) {
          for (feat in features) result[[feat]][idx] <- agg[[feat]][i]
          matched_mask[idx] <- TRUE
          n_pos_avg <- n_pos_avg + 1
        }
      }
    }
  }

  if (verbose >= 1) {
    total <- n_same_ref + n_swapped + n_pos_avg
    if (total > 0) {
      message(sprintf("    Flexible matching: %d variants (same_ref=%d, swapped=%d, pos_avg=%d)",
                      total, n_same_ref, n_swapped, n_pos_avg))
    }
  }

  return(result[, c("VarInfo", features), drop = FALSE])
}


#' Handle Missing Annotation Values
#'
#' @description
#' Applies NA handling strategy to annotation feature columns: keep as-is,
#' replace with zero, or drop variants with any NA.
#'
#' @param data data.frame with annotation feature columns
#' @param features Character vector of feature column names
#' @param na_handling Character. "keep", "zero", or "drop"
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with NA handling applied
#'
#' @details
#' \strong{NA Handling Strategies:}
#'
#' \itemize{
#'   \item "keep": No changes, NA values remain
#'   \item "zero": Replace all NA with 0 (assumes neutral/no effect)
#'   \item "drop": Remove rows with any NA in feature columns
#' }
#'
#' @keywords internal
#' @noRd
.handle_na_annotations <- function(data, features, na_handling, verbose) {

  n_before <- nrow(data)

  if (na_handling == "keep") {
    # No changes
    if (verbose >= 2) {
      # Count NAs for reporting
      n_na <- sum(is.na(data[, features, drop = FALSE]))
      if (n_na > 0) {
        message(sprintf("    Keeping %d NA values in annotations", n_na))
      }
    }

  } else if (na_handling == "zero") {
    # Replace NA with 0
    for (feat in features) {
      if (feat %in% names(data)) {
        na_idx <- is.na(data[[feat]])
        n_na <- sum(na_idx)
        if (n_na > 0) {
          data[[feat]][na_idx] <- 0
          if (verbose >= 2) {
            message(sprintf("    Replaced %d NA values with 0 in %s", n_na, feat))
          }
        }
      }
    }

  } else if (na_handling == "drop") {
    # Remove rows with any NA in feature columns
    feature_cols <- features[features %in% names(data)]
    if (length(feature_cols) > 0) {
      complete_idx <- complete.cases(data[, feature_cols, drop = FALSE])
      data <- data[complete_idx, , drop = FALSE]

      n_removed <- n_before - nrow(data)
      if (verbose >= 1 && n_removed > 0) {
        message(sprintf("    Dropped %d variants with NA annotations", n_removed))
      }
    }
  }

  return(data)
}


#' Check if xsv CLI Tool is Available
#'
#' @description
#' Checks if the xsv command-line tool is installed and accessible in PATH.
#' xsv provides streaming CSV operations that may modestly speed up FAVOR
#' annotation joins for small datasets.
#'
#' @return Logical. TRUE if xsv is available, FALSE otherwise.
#'
#' @details
#' Uses \code{Sys.which()} to check if xsv is in the system PATH.
#' If xsv is not found, the function returns FALSE silently (no warning).
#'
#' Install xsv via: \code{cargo install xsv} (Rust) or Homebrew on macOS.
#'
#' @keywords internal
#' @noRd
.check_xsv_available <- function() {
  xsv_path <- Sys.which("xsv")
  return(nzchar(xsv_path))
}


#' Join Variants with FAVOR Annotations using xsv CLI
#'
#' @description
#' Annotation join using the xsv command-line tool. Performs direct hash-based
#' join between input VarInfo and FAVOR variant_vcf columns.
#'
#' @param variant_data data.frame with VarInfo column
#' @param favor_db_path Character. Path to FAVOR directory
#' @param chunks_needed Character vector of chunk file names
#' @param features Character vector of feature column names to extract
#' @param verbose Integer. Verbosity level
#'
#' @return data.frame with VarInfo and annotation feature columns
#'
#' @details
#' \strong{Algorithm:}
#'
#' \enumerate{
#'   \item Write input variants to temporary CSV (VarInfo column only)
#'   \item For each FAVOR chunk: \code{xsv join --left VarInfo input.csv variant_vcf chunk.csv}
#'   \item Concatenate chunk results: \code{xsv cat rows}
#'   \item Select needed columns: \code{xsv select}
#'   \item Read result back into R
#' }
#'
#' \strong{Performance Notes:}
#'
#' xsv uses memory-mapped I/O and streaming CSV parsing. May provide modest
#' speedup for small datasets due to avoiding full chunk loading into R memory.
#' For large datasets, performance is similar to R-based join as both are
#' I/O bound.
#'
#' Note: Position pre-filtering was removed in 2026-01 due to regex-based
#' implementation having O(N) complexity per row, causing catastrophic
#' performance for large position sets.
#'
#' @keywords internal
#' @noRd
.join_favor_xsv <- function(
  variant_data,
  favor_db_path,
  chunks_needed,
  features,
  verbose = 1
) {

  # Create temp directory for intermediate files
  temp_dir <- tempfile(pattern = "favor_xsv_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE) #this schedules deletion of that directory (and all files inside it) when .join_favor_xsv() finishes, whether it returns normally or errors.

  # Step 1: Write input VarInfo to temp CSV
  input_csv <- file.path(temp_dir, "input_varinfo.csv")
  input_df <- data.frame(VarInfo = variant_data$VarInfo, stringsAsFactors = FALSE)
  data.table::fwrite(input_df, input_csv)

  if (verbose >= 2) {
    message(sprintf("    xsv: Wrote %d variants to temp file", nrow(input_df)))
  }

  # Step 2: Join with each FAVOR chunk using xsv (direct join, no pre-filtering)
  chunk_results <- character(0)

  for (i in seq_along(chunks_needed)) {
    chunk_file <- chunks_needed[i]
    chunk_path <- file.path(favor_db_path, chunk_file)

    if (!file.exists(chunk_path)) {
      if (verbose >= 1) {
        warning(sprintf("Chunk file not found, skipping: %s", chunk_file), call. = FALSE)
      }
      next
    }

    if (verbose >= 2) {
      message(sprintf("    xsv: Joining with chunk %d/%d: %s",
                      i, length(chunks_needed), chunk_file))
    }

    # Output file for this chunk
    chunk_output <- file.path(temp_dir, sprintf("joined_%d.csv", i))

    # Direct join: xsv join --left VarInfo input.csv variant_vcf chunk.csv
    result <- tryCatch({
      system2(
        "xsv",
        args = c("join", "--left", "VarInfo", input_csv, "variant_vcf", chunk_path),
        stdout = chunk_output,
        stderr = FALSE
      )
      0L  # Success
    }, error = function(e) {
      if (verbose >= 2) {
        message(sprintf("      Error: %s", e$message))
      }
      1L  # Failure
    })

    if (result == 0 && file.exists(chunk_output) && file.size(chunk_output) > 0) {
      chunk_results <- c(chunk_results, chunk_output)
    }
  }

  if (length(chunk_results) == 0) {
    if (verbose >= 1) {
      warning("xsv: No chunks produced results", call. = FALSE)
    }
    # Return input with NA annotations
    for (feat in features) {
      variant_data[[feat]] <- NA_real_
    }
    return(variant_data)
  }

  # Step 3: Concatenate chunk results (if multiple)
  if (length(chunk_results) == 1) {
    combined_csv <- chunk_results[1]
  } else {
    combined_csv <- file.path(temp_dir, "combined.csv")

    # xsv cat rows file1.csv file2.csv ... > combined.csv
    # Use system2() for better cross-platform compatibility
    tryCatch({
      system2(
        "xsv",
        args = c("cat", "rows", chunk_results),
        stdout = combined_csv,
        stderr = FALSE
      )
    }, error = function(e) {
      warning("xsv cat rows failed: ", e$message, call. = FALSE)
    })
  }

  # Step 4: Select only needed columns
  # Columns: VarInfo + requested features
  cols_to_select <- c("VarInfo", features)

  # First, check which columns actually exist in the combined file
  # Read header only using R (more portable than head command)
  header_line <- readLines(combined_csv, n = 1)
  available_cols <- strsplit(header_line, ",")[[1]]

  # Filter to columns that exist
  cols_present <- cols_to_select[cols_to_select %in% available_cols]

  if (length(cols_present) < 2) {
    # Only VarInfo, no features found
    if (verbose >= 1) {
      warning("xsv: No requested feature columns found in FAVOR data", call. = FALSE)
    }
    for (feat in features) {
      variant_data[[feat]] <- NA_real_
    }
    return(variant_data)
  }

  # Select columns
  selected_csv <- file.path(temp_dir, "selected.csv")

  # Use system2() for better cross-platform compatibility
  tryCatch({
    system2(
      "xsv",
      args = c("select", paste(cols_present, collapse = ","), combined_csv),
      stdout = selected_csv,
      stderr = FALSE
    )
  }, error = function(e) {
    warning("xsv select failed: ", e$message, call. = FALSE)
  })

  # Step 5: Read result back into R
  if (verbose >= 2) {
    message("    xsv: Reading results back into R")
  }

  annotated <- data.table::fread(selected_csv, data.table = FALSE)

  # Handle duplicates from multiple chunks (keep first match per VarInfo)
  if (any(duplicated(annotated$VarInfo))) {
    annotated <- annotated[!duplicated(annotated$VarInfo), , drop = FALSE]
  }

  # Merge back with original variant_data to preserve order and handle unmatched
  result <- merge(
    variant_data,
    annotated,
    by = "VarInfo",
    all.x = TRUE,
    sort = FALSE
  )

  # Add any missing feature columns as NA
  for (feat in features) {
    if (!feat %in% names(result)) {
      result[[feat]] <- NA_real_
    }
  }

  if (verbose >= 1) {
    n_annotated <- sum(!is.na(result[[features[1]]]))
    message(sprintf("  xsv: Annotated %d/%d variants (%.1f%%)",
                    n_annotated, nrow(result), 100 * n_annotated / nrow(result)))
  }

  return(result)
}


#' Write Annotated Variants to aGDS Format
#'
#' @description
#' Writes annotated variants to aGDS (annotated GDS) format using gdsfmt
#' package. The aGDS format integrates genotype and annotation data for
#' efficient downstream analysis.
#'
#' @param data data.frame with VarInfo and annotation feature columns
#' @param output_path Character. Path for output aGDS file
#' @param features Character vector of feature column names
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (side effect: creates aGDS file)
#'
#' @details
#' \strong{aGDS Structure:}
#'
#' The function creates a GDS file with the following structure:
#' \itemize{
#'   \item /chromosome: Chromosome values
#'   \item /position: Position values
#'   \item /ref: Reference alleles
#'   \item /alt: Alternate alleles
#'   \item /VarInfo: FAVOR format identifiers
#'   \item /annotation/info/FunctionalAnnotation: a \emph{folder} holding one
#'     native-typed sub-node per feature (numeric features stay numeric, string
#'     features stay character), in input row order
#' }
#'
#' The FunctionalAnnotation node uses the STAARpipeline sub-node layout (NOT a
#' single numeric matrix), so it is read by the per-feature variant-set scan via
#' \code{seqGetData(gds, ".../FunctionalAnnotation/<feature>")} and can carry
#' mixed numeric/string feature sets without coercion. This format is compatible
#' with SeqArray and STAAR pipeline tools.
#'
#' @keywords internal
#' @noRd
.write_agds <- function(data, output_path, features, verbose) {

  # Check if gdsfmt package is available
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required for aGDS output. Install with: BiocManager::install('gdsfmt')")
  }

  # Parse VarInfo to extract components
  varinfo_split <- strsplit(data$VarInfo, "-", fixed = TRUE)

  chr_vec <- sapply(varinfo_split, function(x) if (length(x) >= 1) x[1] else NA_character_)
  pos_vec <- as.integer(sapply(varinfo_split, function(x) if (length(x) >= 2) x[2] else NA_character_))
  ref_vec <- sapply(varinfo_split, function(x) if (length(x) >= 3) x[3] else NA_character_)
  alt_vec <- sapply(varinfo_split, function(x) if (length(x) >= 4) x[4] else NA_character_)

  # Create GDS file
  gds_file <- gdsfmt::createfn.gds(output_path)

  tryCatch({
    # Add chromosome node
    gdsfmt::add.gdsn(gds_file, "chromosome", chr_vec, compress = "LZMA_RA", closezip = TRUE)

    # Add position node
    gdsfmt::add.gdsn(gds_file, "position", pos_vec, compress = "LZMA_RA", closezip = TRUE)

    # Add ref allele node
    gdsfmt::add.gdsn(gds_file, "ref", ref_vec, compress = "LZMA_RA", closezip = TRUE)

    # Add alt allele node
    gdsfmt::add.gdsn(gds_file, "alt", alt_vec, compress = "LZMA_RA", closezip = TRUE)

    # Add VarInfo node
    gdsfmt::add.gdsn(gds_file, "VarInfo", data$VarInfo, compress = "LZMA_RA", closezip = TRUE)

    # Create annotation group
    annot_group <- gdsfmt::addfolder.gdsn(gds_file, "annotation")
    info_group <- gdsfmt::addfolder.gdsn(annot_group, "info")

    # Write FunctionalAnnotation as a STAARpipeline-style sub-node folder (one
    # native-typed sub-node per feature). `data` is already in the intended row
    # order here, so no reordering is needed.
    .write_functional_annotation_subnodes(
      parent_node = info_group,
      annotations = data,
      features = features,
      verbose = verbose
    )

  }, finally = {
    gdsfmt::closefn.gds(gds_file)
  })

  if (verbose >= 2) {
    message(sprintf("    Created aGDS with %d variants and %d features",
                    nrow(data), length(features)))
  }
}


#' Write FunctionalAnnotation as a Sub-Node Folder (STAARpipeline Format)
#'
#' @description
#' Creates (or overwrites) a \code{FunctionalAnnotation} folder node under
#' \code{parent_node} and adds one sub-node per feature, each carrying that
#' feature's native-typed vector. This is the STAARpipeline aGDS layout that the
#' per-feature variant-set scan reads via
#' \code{seqGetData(gds, ".../FunctionalAnnotation/<feature>")}, and it mirrors
#' the \code{.add_anno} pattern used by the test helper
#' \code{create_test_agds()}.
#'
#' @param parent_node A \code{gdsn.class} folder node (typically
#'   \code{annotation/info}) under which the FunctionalAnnotation folder is
#'   created.
#' @param annotations data.frame with at least the \code{features} columns, in
#'   the desired variant (row) order. Callers are responsible for ordering the
#'   rows to match the GDS variant order before calling.
#' @param features Character vector of feature column names to write as
#'   sub-nodes.
#' @param verbose Integer. Verbosity level.
#'
#' @return Invisibly, the created \code{FunctionalAnnotation} folder node
#'   (\code{gdsn.class}).
#'
#' @details
#' \strong{Why a folder of sub-nodes (not a matrix):} a single numeric matrix
#' node exposes no per-feature sub-nodes, so the scan's
#' \code{seqGetData(gds, ".../FunctionalAnnotation/<feature>")} read fails; a
#' matrix also coerces mixed numeric/string features to one type, dropping string
#' coding nodes (e.g. \code{genecode_comprehensive_category}). Writing one typed
#' sub-node per feature avoids both problems.
#'
#' \strong{Native type preservation:} each sub-node is written from
#' \code{annotations[[feature]]} directly (no \code{as.matrix}), so numeric
#' features stay numeric and character features stay character.
#'
#' \strong{Overwrite:} if a \code{FunctionalAnnotation} node already exists under
#' \code{parent_node}, it is deleted (with a warning when \code{verbose >= 1})
#' before the new folder is created.
#'
#' A \code{feature_names} attribute is also placed on the folder node for
#' discoverability; the sub-node names remain the source of truth.
#'
#' @keywords internal
#' @noRd
.write_functional_annotation_subnodes <- function(parent_node, annotations,
                                                   features, verbose = 1) {

  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required for aGDS output. Install with: BiocManager::install('gdsfmt')")
  }

  # Overwrite any pre-existing FunctionalAnnotation node (folder or matrix).
  fa_exists <- "FunctionalAnnotation" %in% gdsfmt::ls.gdsn(parent_node)
  if (fa_exists) {
    if (verbose >= 1) {
      warning("Overwriting existing FunctionalAnnotation node", call. = FALSE)
    }
    gdsfmt::delete.gdsn(
      gdsfmt::index.gdsn(parent_node, "FunctionalAnnotation"),
      force = TRUE
    )
  }

  # Create the FunctionalAnnotation folder, then one typed sub-node per feature.
  fa_folder <- gdsfmt::addfolder.gdsn(parent_node, "FunctionalAnnotation")

  for (feat in features) {
    if (!feat %in% names(annotations)) {
      # Skip silently-absent features under verbose 0; warn otherwise. This keeps
      # the writer robust if a requested feature was not produced by the join.
      if (verbose >= 1) {
        warning(sprintf("Feature '%s' not present in annotations; skipping sub-node",
                        feat), call. = FALSE)
      }
      next
    }
    vec <- annotations[[feat]]
    # Preserve native type: numeric stays numeric, character stays character.
    # Factors (rare here) are flattened to character so the scan reads strings.
    if (is.factor(vec)) vec <- as.character(vec)
    gdsfmt::add.gdsn(fa_folder, feat, val = vec,
                     compress = "LZMA_RA", closezip = TRUE)
  }

  # Record feature_names on the folder for discoverability (sub-node names are
  # authoritative).
  gdsfmt::put.attr.gdsn(fa_folder, "feature_names", features)

  invisible(fa_folder)
}
