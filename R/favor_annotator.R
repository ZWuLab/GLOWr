########## FAVOR Annotation Module ##########
#
# This file provides functions to annotate variants with functional annotation
# scores from the FAVOR database. Used in PI estimation to assign variant-importance
# scores based on functional features, and to build the cohort aGDS that GLOW and
# STAAR scan.
#
# EXPORTED FUNCTIONS:
#   - annotate_favor()   Annotate variants with FAVOR scores (CSV / aGDS output)
#
# INTERNAL HELPERS (selected):
#   - .favor_annotate_rows()                   match every input row and extract features
#   - .favor_parse_keys()                      CHR-POS-REF-ALT -> lookup fields + input class
#   - .favor_classify()                        matching outcome per allele-bearing row
#   - .favor_rsid_check()                      input rsID against the matched row's (D9)
#   - .favor_position_values()                 position-only aggregation (SNV rows, D10)
#   - .favor_csv_rows() / .favor_parquet_rows() database rows at the cohort's positions
#   - .favor2_leaf_paths()                     Parquet leaf columns by dotted path
#   - .favor2_field_map() / .favor2_serialize() FAVOR 2.0 sources and v1-layout strings
#   - .favor_provenance()                      provenance attributes for the aGDS
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
#' (Functional Annotation of Variants Online Resource) database, read either
#' from the FAVOR v1 CSV chunk files or from the FAVOR 2.0 per-chromosome
#' Parquet files. Every input row receives a matching outcome
#' (\code{match_tier}), when matched the key of the database row that
#' supplied its annotation (\code{favor_key}), and the rsID evidence for that
#' row (\code{rsid_check}). Used in PI estimation to retrieve
#' pathogenicity-relevant features for trait-associated and control variants,
#' and to build the cohort aGDS that GLOW and STAAR scan.
#'
#' @param variants One of the following:
#'   \itemize{
#'     \item A data.frame with VarInfo column (CHR-POS-REF-ALT format), and
#'       optionally an \code{rsID} column (the input's rsID per row), which
#'       feeds the rsID check described in Details
#'     \item An S3 object of class \code{glow_pi_case_data} or
#'       \code{glow_pi_control_data} from \code{\link{prepare_PI_case_data}} or
#'       \code{\link{prepare_PI_control_data}}
#'     \item A path (character) to a GDS file containing variant data (requires
#'       SeqArray package). Variants are extracted automatically, and the
#'       GDS \code{annotation/id} (the rsID on a chip cohort) is read as the
#'       \code{rsID} column.
#'   }
#' @param favor_db_path Character. The database directory: either the FAVOR v1
#'   CSV chunk files (\code{chr1_1.csv}, \code{chr1_2.csv}, ...) or the FAVOR 2.0
#'   Parquet files (\code{chromosome_1.parquet}, ...).
#' @param favor_split_file Character or NULL. CSV backend only. Path to
#'   \code{FAVORdatabase_chrsplit.csv}, which maps positions to chunk files. If
#'   NULL (default), searches in this order: (1) \code{favor_db_path}, (2) the
#'   package's bundled file.
#' @param features Character vector. Annotation fields to extract, named as the
#'   FAVOR v1 columns. The default is the complete annotation content of the
#'   FAVOR Essential Database (30 fields): 20 numeric scores (17 annotation
#'   principal components, CADD PHRED, LINSIGHT, FATHMM-XF) plus 10 categorical
#'   fields (GENCODE category/info and exonic category/info, MetaSVM prediction,
#'   GeneHancer, CAGE, rDHS, rsID). The categorical fields are what the built-in
#'   coding masks read, so an aGDS built with the default is directly scannable
#'   by every built-in variant category. Fields the database does not carry are
#'   skipped and reported (FAVOR 2.0 lacks four older APC versions, so the
#'   default yields 26 fields there).
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
#' @param match_method Character. How allele-bearing variants are matched:
#'   \itemize{
#'     \item "exact" (default): the key as given only (STAAR-compatible)
#'     \item "flexible": for SNVs, the lookup key is also normalized against
#'       FAVOR's reference base by swapping REF and ALT or complementing both
#'       alleles, for cohorts whose alleles may differ from the reference in
#'       order or strand (genotyping-chip data). Indels always match by the
#'       exact key. See Details.
#'   }
#' @param na_allele_method Character. How position-only keys
#'   (CHR-POS-NA-NA) are annotated from the SNV rows at the position (indel
#'   rows are never included; a position with no SNV row is "uncovered"):
#'   \itemize{
#'     \item "average" (default): numeric features averaged over the SNV rows'
#'       non-missing values, string features take the first non-empty value
#'     \item "first": every feature from the first SNV row at the position
#'   }
#' @param na_handling Character. How to handle missing annotation values
#'   after matching (never changes \code{match_tier}):
#'   \itemize{
#'     \item "keep" (default): Keep NA values as-is
#'     \item "zero": Replace NA with 0 in numeric features
#'     \item "drop": Remove variants with any NA annotations from the
#'       returned table and the CSV. The aGDS writers keep every variant of
#'       the GDS (the file stays a copy of the cohort) and leave the dropped
#'       variants' annotation and tier entries empty.
#'   }
#' @param use_xsv Logical. CSV backend only. If TRUE (default) and the xsv CLI
#'   tool is available, xsv fetches the database rows at the input positions
#'   (a streaming join on \code{position}); otherwise each chunk is read with
#'   \code{data.table::fread()} and filtered. Both give the same rows.
#' @param verbose Integer. Verbosity level: 0=silent, 1=messages (default),
#'   2=detailed messages
#' @param favor_db_format Character. "auto" (default), "csv" or "parquet".
#'   Under "auto" the format is detected from the files in
#'   \code{favor_db_path}; a directory holding both kinds, or neither, is an
#'   error. (Placed after \code{verbose} so that positional calls written for
#'   GLOWr 0.1.1 keep their meaning.)
#' @param favor_release Character or NULL. The FAVOR release the database
#'   files come from, recorded in the aGDS provenance. NULL records "unknown".
#' @param rsid_policy Character. What the rsID check does with a match:
#'   \itemize{
#'     \item "require" (default): a \code{swapped}, \code{flipped} or
#'       \code{flipped_swapped} match whose check is \code{differs} or
#'       \code{favor_none} is withheld as the non-match \code{rsid_conflict}
#'       (its annotation stays missing; \code{favor_key} names the withheld
#'       row). Exact matches and matches checked \code{chip_none} are kept,
#'       so exact matching is unaffected and an input without rsIDs loses
#'       nothing.
#'     \item "record": every match is kept and the check is only recorded
#'       in \code{rsid_check}, for diagnostics.
#'   }
#'
#' @return data.frame with the input columns (for GDS input also
#'   \code{variant_id} and \code{rsID}), one column per served feature,
#'   \code{match_tier} (the matching outcome of each row), \code{favor_key}
#'   (the \code{variant_vcf} of the database row that supplied the
#'   annotation, or of the row withheld under \code{rsid_conflict}; NA
#'   otherwise) and \code{rsid_check} (the rsID evidence of every matched
#'   row; NA for non-matches and position-only keys). Order of variants is
#'   preserved from input.
#'
#' @details
#' \strong{Database formats:}
#' \itemize{
#'   \item FAVOR v1 CSV: \code{chr{N}_{K}.csv} chunk files keyed by
#'     \code{variant_vcf} ("CHR-POS-REF-ALT"), with
#'     \code{FAVORdatabase_chrsplit.csv} mapping positions to chunks (a copy is
#'     bundled with GLOWr).
#'   \item FAVOR 2.0 Parquet: \code{chromosome_{N}.parquet}, one file per
#'     chromosome (1-22, X, Y), sorted by position, with the annotation in
#'     nested struct columns, read leaf by leaf one row group at a time.
#'     Requires the \pkg{arrow} package. The map from our field names to the Parquet sources is
#'     internal (\code{.favor2_field_map()}); GeneHancer, GENCODE exonic
#'     information and GENCODE information are serialized into the v1 string
#'     layouts. FAVOR 2.0 stores scores as float32, which R reads as the exact
#'     double value of the float32.
#' }
#' The users obtain the database from the FAVOR team; GLOWr does not
#' redistribute it.
#'
#' \strong{Matching (the lookup key is normalized, never the data):}
#' Each input row is classified as an SNV (REF and ALT each one of A, C, G, T),
#' an indel (A/C/G/T strings, not both of length one), a position-only key
#' (both alleles "NA"), or unsupported (anything else, including symbolic
#' alleles and multiallelic records). For an SNV at position p with FAVOR
#' reference base r:
#' \enumerate{
#'   \item "exact": the key CHR-POS-REF-ALT exists.
#'   \item "swapped" (flexible): r equals ALT and CHR-POS-ALT-REF exists.
#'   \item "flipped" / "flipped_swapped" (flexible): r equals the complement
#'     of REF (or of ALT) and the complemented key exists.
#' }
#' Indels match by the exact key only, because a VCF indel's REF is the genome
#' sequence and its reciprocal key names a different variant. Non-matches carry
#' a reason: "unmatched_alt" (the position is covered and the reference base
#' is found, either by REF as given or, under flexible matching, after one of
#' the enabled transformations, but no row carries the resulting REF-ALT pair),
#' "unmatched_ref" (the position is covered, REF as given does not agree with
#' r, and no enabled transformation finds r), "uncovered" (no database row at the position),
#' "no_database" (the database has no file or chunk for the chromosome),
#' "unsupported_allele" and "rsid_conflict" (below). Position-only keys get
#' "position_only", with the features aggregated over the SNV rows at the
#' position only (\code{na_allele_method}); a position that carries no SNV row
#' is "uncovered". A tier records
#' the lookup transformation, not a verified strand: a non-palindromic SNV pair
#' at a covered position always resolves through the first four tiers,
#' including a pair of two non-reference alleles. The cohort's own alleles,
#' genotypes and keys are never modified.
#'
#' \strong{The rsID as evidence, never a key.} FAVOR carries dbSNP's rsID on
#' the row of the substitution dbSNP lists, and the other rows at the position
#' carry none, so agreement between the input's rsID and the matched row's
#' confirms the site and the allele pair. For every
#' matched row \code{rsid_check} records "same", "differs", "chip_none" (the
#' input carries no rsID of the form \code{rs<digits>}) or "favor_none" (the
#' database row carries none). The input's rsID is the GDS \code{annotation/id}
#' or the data.frame's \code{rsID} column; without either, every check is
#' "chip_none". Under \code{rsid_policy = "require"}, the default, a
#' transformed match checked "differs" or "favor_none" becomes
#' "rsid_conflict"; \code{"record"} keeps it. The check reads
#' the database's rsID column whether or not \code{rsid} is a requested
#' feature.
#'
#' \strong{Outputs:} the aGDS writers align rows to the GDS through
#' \code{variant.id}, write one native-typed sub-node per feature under
#' \code{annotation/info/FunctionalAnnotation}, the matching outcome as
#' \code{annotation/info/favor_match_tier}, the rsID check as
#' \code{annotation/info/favor_rsid_check}, and provenance (database format,
#' files, sizes and SHA-256 when a \code{SHA256SUMS} file is present, release,
#' GLOWr version, matching settings, the rsID policy and source, the rows the
#' position-only aggregation uses, outcome and check counts, date) as
#' attributes of the \code{FunctionalAnnotation} folder.
#'
#' \strong{Computational Complexity:} the matching is vectorized and linear in
#' the number of input rows plus the number of database rows at their
#' positions. Reading dominates: the CSV backend streams each chunk that covers
#' the input positions; the Parquet backend reads the position column of every
#' row group and the requested leaf columns of the groups that hold input
#' positions, one group at a time.
#'
#' @examples
#' \dontrun{
#' # Annotate a chip cohort's GDS from FAVOR 2.0, normalizing the lookup key for
#' # allele order and strand, and withholding transformed matches that the
#' # chip's rsID does not confirm
#' annotate_favor(
#'   variants = "data/chr22.gds",
#'   favor_db_path = "data/favor2-db/FAVOR2.0",
#'   match_method = "flexible",
#'   rsid_policy = "require",
#'   output_agds = "results/chr22_favor2.gds"
#' )
#'
#' # Annotate case variants given by position only, from FAVOR v1 CSV
#' annotated_cases <- annotate_favor(
#'   variants = cases,
#'   favor_db_path = "data/FAVOR",
#'   na_allele_method = "average"
#' )
#' table(annotated_cases$match_tier)
#'
#' # GDS input with a position range filter
#' annotated <- annotate_favor(
#'   variants = "data/genotypes.gds",
#'   favor_db_path = "data/FAVOR",
#'   variant_filter = list(chr = "21", start = 1e6, end = 5e6)
#' )
#' }
#'
#' @references
#' Zhou, H., Arapoglou, T., Li, X., et al. (2023). FAVOR: functional annotation of
#' variants online resource and annotator for variation across the human genome.
#' Nucleic Acids Research, 51(D1), D1300-D1311. doi:10.1093/nar/gkac966
#'
#' Zhou, H., Verma, V., Li, X., et al. (2026). FAVOR 2.0: A reengineered
#' functional annotation of variants online resource for interpreting genomic
#' variation. Nucleic Acids Research, 54(D1), D1405-D1414. doi:10.1093/nar/gkaf1217
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
  verbose = 1,
  favor_db_format = c("auto", "csv", "parquet"),
  favor_release = NULL,
  rsid_policy = c("require", "record")
) {

  # ========== Step 1: Validate Inputs ==========

  if (verbose >= 1) {
    message("=== FAVOR Annotation ===")
  }
  if (!match_method %in% c("exact", "flexible")) {
    stop("match_method must be 'exact' or 'flexible'")
  }
  if (!na_allele_method %in% c("average", "first")) {
    stop("na_allele_method must be 'average' or 'first'")
  }
  if (!na_handling %in% c("keep", "zero", "drop")) {
    stop("na_handling must be 'keep', 'zero', or 'drop'")
  }
  favor_db_format <- match.arg(favor_db_format)
  rsid_policy <- match.arg(rsid_policy)
  if (!dir.exists(favor_db_path)) {
    stop("FAVOR database directory not found: ", favor_db_path)
  }
  if (!is.null(favor_release) &&
      !(is.character(favor_release) && length(favor_release) == 1L)) {
    stop("favor_release must be NULL or a single character string")
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

  # Validate VarInfo column exists
  if (!"VarInfo" %in% names(variant_data)) {
    stop("variants must have a 'VarInfo' column (CHR-POS-REF-ALT format)")
  }
  # A plain data.frame: element-wise assignment into a data.table can fail
  # silently (see GLOWr-package.R), and every later step indexes by row.
  variant_data <- as.data.frame(variant_data, stringsAsFactors = FALSE)
  variant_data$VarInfo <- as.character(variant_data$VarInfo)

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
      parsed_pos <- suppressWarnings(
        as.integer(sub("^[^-]+-([0-9]+)-.*", "\\1", variant_data$VarInfo)))
      in_range <- !is.na(parsed_pos) & parsed_pos >= variant_filter$start &
        parsed_pos <= variant_filter$end
      variant_data <- variant_data[in_range, , drop = FALSE]
      if (verbose >= 1) {
        message(sprintf("  Filtered to position range %d-%d: %d variants",
                        variant_filter$start, variant_filter$end, nrow(variant_data)))
      }
    }
    if (nrow(variant_data) == 0) {
      warning("No variants remain after applying variant_filter", call. = FALSE)
    }
  }

  n_input <- nrow(variant_data)
  if (verbose >= 1) {
    message(sprintf("Input: %d variants", n_input))
  }
  # The input's rsIDs, the evidence of the rsID check (decision D9): the GDS
  # annotation/id (read into rsID by .extract_varinfo_from_gds) or a data.frame's
  # rsID column. Without either, every check is "chip_none".
  rsid_source <- if (!"rsID" %in% names(variant_data)) "none" else
    if (!is.null(gds_input_path)) "annotation/id" else "rsID column"
  rsid_in <- if (rsid_source == "none") NULL else variant_data$rsID

  # ========== Step 3: Choose the Database Backend ==========

  backend <- .favor_resolve_backend(favor_db_path, favor_db_format)
  split_data <- NULL
  if (backend == "csv") {
    split_data <- .favor_load_split_file(favor_db_path, favor_split_file, verbose)
  } else if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required for the FAVOR 2.0 Parquet backend.")
  }
  if (verbose >= 1) {
    message(sprintf("Database: %s (%s), match_method = %s", favor_db_path, backend, match_method))
  }

  # ========== Step 4: Match Every Row and Extract Features ==========

  res <- .favor_annotate_rows(
    variant_data     = variant_data,
    backend          = backend,
    favor_db_path    = favor_db_path,
    split_data       = split_data,
    features         = features,
    match_method     = match_method,
    na_allele_method = na_allele_method,
    use_xsv          = use_xsv,
    verbose          = verbose,
    rsid             = rsid_in,
    rsid_policy      = rsid_policy
  )
  annotated <- res$data
  features_served <- res$features

  # ========== Step 5: Handle Missing Annotations ==========

  annotated <- .handle_na_annotations(
    data = annotated,
    features = features_served,
    na_handling = na_handling,
    verbose = verbose
  )
  n_output <- nrow(annotated)
  if (verbose >= 1 && n_output < n_input) {
    message(sprintf("  %d variants removed due to NA handling", n_input - n_output))
  }

  provenance <- .favor_provenance(
    backend = backend, favor_db_path = favor_db_path, files = res$files,
    favor_release = favor_release, match_method = match_method,
    na_allele_method = na_allele_method, na_handling = na_handling,
    tiers = res$data$match_tier, features_missing = res$features_missing,
    rsid_policy = rsid_policy, rsid_source = rsid_source,
    rsid_checks = res$data$rsid_check
  )

  # ========== Step 6: Write Outputs ==========

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
        features = features_served,
        provenance = provenance,
        verbose = verbose
      )
    } else {
      # Non-GDS input -> annotation-only GDS
      if (verbose >= 1) {
        message(sprintf("Writing annotation-only GDS to: %s", output_agds))
      }
      .write_agds(annotated, output_agds, features = features_served,
                  provenance = provenance, verbose = verbose)
    }
  }

  if (update_gds && !is.null(gds_input_path)) {
    if (verbose >= 1) {
      message(sprintf("Updating GDS file with annotations: %s", basename(gds_input_path)))
    }
    .update_gds_with_annotations(
      gds_path = gds_input_path,
      annotations = annotated,
      features = features_served,
      provenance = provenance,
      verbose = verbose
    )
  } else if (update_gds && is.null(gds_input_path)) {
    warning("update_gds=TRUE ignored: input was not a GDS file", call. = FALSE)
  }

  # ========== Summary ==========

  if (verbose >= 1) {
    message("\n=== Annotation Complete ===")
    message(sprintf("  Input:  %d variants", n_input))
    message(sprintf("  Output: %d variants", n_output))
    message(sprintf("  Features: %d served (%s)%s", length(features_served),
                    paste(utils::head(features_served, 3), collapse = ", "),
                    if (length(features_served) > 3) ", ..." else ""))
    message("  Matching outcomes: ", .favor_tier_summary(res$data$match_tier))
    message(sprintf("  rsID check (policy %s, source %s): %s", rsid_policy, rsid_source,
                    .favor_rsid_summary(res$data$rsid_check)))
  }

  return(annotated)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#################### Matching core (both backends) ####################

#' Match Every Input Row Against the Database and Extract Its Features
#'
#' @description
#' Per chromosome: fetch every database row at the input positions (from the
#' CSV or Parquet backend), classify each allele-bearing row with
#' \code{.favor_classify()}, aggregate position-only keys, and write the
#' feature values, the outcome and the matched database key into the result.
#' Resolution is per input row, so duplicate keys each get the annotation and
#' the input order is kept.
#'
#' @param variant_data data.frame with a character VarInfo column.
#' @param backend "csv" or "parquet".
#' @param favor_db_path Database directory.
#' @param split_data The chunk table (CSV backend) or NULL.
#' @param features Requested feature names.
#' @param match_method "exact" or "flexible".
#' @param na_allele_method "average" or "first".
#' @param use_xsv Logical (CSV backend).
#' @param verbose Integer.
#' @param rsid The input's rsIDs (character, one per row) or NULL.
#' @param rsid_policy "record" or "require" (decision D9).
#'
#' @return list(data = data.frame with features, match_tier, favor_key,
#'   rsid_check; features = the features the database served;
#'   features_missing; files = database files read).
#' @keywords internal
#' @noRd
.favor_annotate_rows <- function(variant_data, backend, favor_db_path, split_data,
                                 features, match_method, na_allele_method,
                                 use_xsv, verbose, rsid = NULL, rsid_policy = "record") {
  n <- nrow(variant_data)
  keys <- .favor_parse_keys(variant_data$VarInfo)
  tier <- rep(NA_character_, n)
  fkey <- rep(NA_character_, n)
  rcheck <- rep(NA_character_, n)
  # The input's rsIDs as evidence (never a lookup key): NA where none is given.
  chip_rs <- .favor_norm_rsid(if (is.null(rsid)) rep(NA_character_, n) else rsid)
  str_feats <- .favor_string_features()
  # Typed NA per feature: an untyped (logical) all-NA column would be written to
  # the GDS as an integer node.
  vals <- lapply(features, function(f)
    if (f %in% str_feats) rep(NA_character_, n) else rep(NA_real_, n))
  names(vals) <- features
  served <- character(0)
  files_read <- character(0)
  missing_feats <- character(0)

  tier[keys$class == "unsupported"] <- "unsupported_allele"
  todo <- which(keys$class != "unsupported")

  for (chr in unique(keys$chr_lookup[todo])) {
    idx <- todo[keys$chr_lookup[todo] == chr]
    positions <- sort(unique(keys$pos[idx]))
    got <- if (backend == "csv") {
      .favor_csv_rows(chr, positions, favor_db_path, split_data, features,
                      use_xsv = use_xsv, verbose = verbose)
    } else {
      .favor_parquet_rows(chr, positions, favor_db_path, features, verbose = verbose)
    }
    if (is.null(got)) {
      # The database has no file or chunk for this chromosome at all.
      tier[idx] <- "no_database"
      warning(sprintf("FAVOR database at %s has no data for chromosome %s; %d variant(s) recorded as 'no_database'.",
                      favor_db_path, chr, length(idx)), call. = FALSE)
      next
    }
    fr <- got$rows
    files_read <- union(files_read, got$files)
    served <- union(served, got$served)
    missing_feats <- union(missing_feats, got$missing)

    # Allele-bearing rows: one database row per input row, or a reason.
    al <- idx[keys$class[idx] %in% c("snv", "indel")]
    if (length(al)) {
      cl <- .favor_classify(keys$pos[al], keys$ref[al], keys$alt[al],
                            keys$class[al] == "snv",
                            fr$position, fr$ref_vcf, fr$alt_vcf, match_method)
      tier[al] <- cl$tier
      hit <- !is.na(cl$row)
      if (any(hit)) {
        hi <- al[hit]; rows <- cl$row[hit]
        fkey[hi] <- paste(chr, fr$position[rows], fr$ref_vcf[rows], fr$alt_vcf[rows],
                          sep = "-")
        # The rsID evidence (decision D9): the input's rsID against the matched
        # row's. FAVOR writes dbSNP's rsID on the specific REF-ALT row dbSNP knows,
        # so "same" confirms the site and the allele pair.
        rcheck[hi] <- .favor_rsid_check(chip_rs[hi], fr$.rsid[rows])
        # Under "require", a transformed match that dbSNP does not confirm is
        # withheld: the outcome becomes rsid_conflict, favor_key keeps the withheld
        # row, and the annotation stays missing. Exact matches are never withheld,
        # nor is a match whose input carries no rsID (chip_none).
        keep <- rep(TRUE, length(hi))
        if (rsid_policy == "require") {
          keep <- !(tier[hi] %in% c("swapped", "flipped", "flipped_swapped") &
                    rcheck[hi] %in% c("differs", "favor_none"))
          tier[hi[!keep]] <- "rsid_conflict"
        }
        if (any(keep)) {
          fv <- .favor_feature_values(fr, rows[keep], got$served, backend)
          for (f in names(fv)) vals[[f]][hi[keep]] <- fv[[f]]
        }
      }
    }

    # Position-only keys: aggregate over the SNV rows at the position (decision
    # D10). Indel rows are never included, so a position covered only by indel
    # rows is "uncovered" for a position-only key.
    po <- idx[keys$class[idx] == "position_only"]
    if (length(po)) {
      snv_row <- grepl("^[ACGT]$", fr$ref_vcf) & grepl("^[ACGT]$", fr$alt_vcf)
      covered <- keys$pos[po] %in% fr$position[snv_row]
      tier[po] <- ifelse(covered, "position_only", "uncovered")
      if (any(covered)) {
        prow <- which(snv_row & fr$position %in% keys$pos[po][covered])
        ft <- .favor_feature_values(fr, prow, got$served, backend)
        pv <- .favor_position_values(keys$pos[po][covered], fr$position[prow], ft,
                                     na_allele_method)
        for (f in names(pv)) vals[[f]][po[covered]] <- pv[[f]]
      }
    }
  }

  # Features the database does not carry are dropped from the output and
  # reported, rather than written as all-missing nodes.
  missing_feats <- setdiff(missing_feats, served)
  if (length(missing_feats) && verbose >= 1) {
    message(sprintf("  Requested features not in this database, skipped: %s",
                    paste(missing_feats, collapse = ", ")))
  }
  # When no database file was read at all, keep the request as typed NA columns.
  keep <- if (length(files_read)) features[features %in% served] else features
  result <- variant_data
  for (f in keep) result[[f]] <- vals[[f]]
  result$match_tier <- tier
  result$favor_key <- fkey
  result$rsid_check <- rcheck
  rownames(result) <- NULL

  if (verbose >= 1) {
    present <- keep[keep %in% names(result)]
    any_annot <- if (length(present))
      Reduce(`|`, lapply(present, function(f) !is.na(result[[f]]))) else rep(FALSE, n)
    message(sprintf("  Annotated %d/%d variants (%.1f%%) with at least one non-missing feature",
                    sum(any_annot), n, if (n) 100 * sum(any_annot) / n else 0))
  }
  list(data = result, features = keep, features_missing = missing_feats,
       files = files_read)
}


#' Parse CHR-POS-REF-ALT Keys into Lookup Fields and an Input Class
#'
#' @description
#' Splits each key into chromosome, position, REF and ALT, and assigns one
#' input class: "snv" (REF and ALT each one of A, C, G, T, and different),
#' "indel" (strings over A, C, G, T, not both of length one, and different),
#' "position_only" (both alleles "NA") or "unsupported" (anything else:
#' symbolic or non-ACGT alleles, a single missing allele, identical alleles, a
#' multiallelic ALT with a comma, a malformed key). The lookup alleles are
#' upper-cased and the lookup chromosome drops a leading "chr"; the input's
#' own key is not modified.
#'
#' @param varinfo Character vector of keys.
#' @return data.frame(chr, chr_lookup, pos, ref, alt, class).
#' @keywords internal
#' @noRd
.favor_parse_keys <- function(varinfo) {
  n <- length(varinfo)
  sp <- data.table::tstrsplit(varinfo, "-", fixed = TRUE, fill = NA_character_)
  part <- function(k) if (length(sp) >= k) sp[[k]] else rep(NA_character_, n)
  n_parts <- nchar(varinfo) - nchar(gsub("-", "", varinfo, fixed = TRUE)) + 1L
  chr <- part(1)
  # A position is a positive integer written in decimal digits. as.integer()
  # alone would accept "101.9" (truncated to 101), "1e3" and "0", and the key
  # would then receive another coordinate's annotation, so such keys are refused.
  pos_txt <- part(2)
  pos_digits <- !is.na(pos_txt) & grepl("^[0-9]+$", pos_txt)
  pos <- suppressWarnings(as.integer(ifelse(pos_digits, pos_txt, NA_character_)))
  pos_bad <- n_parts == 4L & !is.na(pos_txt) & (!pos_digits | is.na(pos) | pos < 1L)
  if (any(pos_bad, na.rm = TRUE)) {
    bad <- varinfo[which(pos_bad)]
    stop(sprintf("%d VarInfo key(s) have a malformed position (a positive integer in decimal digits is required, at most %d): %s%s",
                 length(bad), .Machine$integer.max, paste(utils::head(bad, 5), collapse = ", "),
                 if (length(bad) > 5) ", ..." else ""), call. = FALSE)
  }
  ref <- toupper(part(3))
  alt <- toupper(part(4))
  acgt <- function(x) !is.na(x) & grepl("^[ACGT]+$", x)
  well_formed <- !is.na(varinfo) & n_parts == 4L & !is.na(chr) & nzchar(chr) & !is.na(pos)
  pos_only <- well_formed & ref %in% "NA" & alt %in% "NA"
  both <- well_formed & acgt(ref) & acgt(alt) & ref != alt
  snv <- both & nchar(ref) == 1L & nchar(alt) == 1L
  class <- ifelse(pos_only, "position_only",
           ifelse(snv, "snv", ifelse(both, "indel", "unsupported")))
  data.frame(chr = chr, chr_lookup = sub("^chr", "", chr, ignore.case = TRUE),
             pos = pos, ref = ref, alt = alt, class = class,
             stringsAsFactors = FALSE)
}


#' Matching Outcomes, in the Order the Design Lists Them
#' @keywords internal
#' @noRd
.favor_tier_levels <- function() {
  c("exact", "swapped", "flipped", "flipped_swapped", "position_only",
    "rsid_conflict", "unmatched_alt", "unmatched_ref", "uncovered",
    "no_database", "unsupported_allele")
}


#' Values of the rsID Check (Decision D9)
#' @keywords internal
#' @noRd
.favor_rsid_check_levels <- function() {
  c("same", "differs", "chip_none", "favor_none")
}


#' Normalize Input rsIDs: rs<digits> (any case) Counts, Anything Else Is None
#' @keywords internal
#' @noRd
.favor_norm_rsid <- function(x) {
  x <- trimws(as.character(x))
  ok <- !is.na(x) & grepl("^rs[0-9]+$", x, ignore.case = TRUE)
  ifelse(ok, tolower(x), NA_character_)
}


#' The rsID Check of Matched Rows
#'
#' @description "chip_none" when the input carries no rsID, else "favor_none"
#'   when the database row carries none, else "same" or "differs" (compared
#'   case-insensitively). The rsID is evidence about the matched row, never a
#'   lookup key.
#' @param chip Normalized input rsIDs (NA for none), one per matched row.
#' @param db The database rows' rsid values.
#' @keywords internal
#' @noRd
.favor_rsid_check <- function(chip, db) {
  db <- .favor_norm_rsid(db)
  ifelse(is.na(chip), "chip_none",
         ifelse(is.na(db), "favor_none",
                ifelse(chip == db, "same", "differs")))
}


#' Count per rsID Check Value, as One Line
#' @keywords internal
#' @noRd
.favor_rsid_summary <- function(check) {
  tab <- table(factor(check, levels = .favor_rsid_check_levels()))
  tab <- tab[tab > 0]
  if (!length(tab)) return("no matched row")
  paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
}


#' Count per Matching Outcome, as One Line
#' @keywords internal
#' @noRd
.favor_tier_summary <- function(tier) {
  tab <- table(factor(tier, levels = .favor_tier_levels()))
  tab <- tab[tab > 0]
  if (!length(tab)) return("none")
  paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", ")
}


#' Watson-Crick Complement of Single-Base Alleles (NA for Anything Else)
#' @keywords internal
#' @noRd
.favor_complement <- function(a) {
  unname(c(A = "T", C = "G", G = "C", T = "A")[a])
}


#' Classify Allele-Bearing Variants of One Chromosome
#'
#' @description
#' Implements the matching rule of the design
#' corrected design"). The lookup key is normalized, never the data.
#'
#' Let r be FAVOR's reference base at the position: the first base of the REF
#' of any database row there (SNV rows carry the reference base, indel rows are
#' anchored on it). For an SNV (a1, a2):
#' \enumerate{
#'   \item exact: key pos-a1-a2 exists;
#'   \item swapped (flexible): r == a2 and pos-a2-a1 exists;
#'   \item flipped (flexible): r == comp(a1) and pos-comp(a1)-comp(a2) exists;
#'   \item flipped_swapped (flexible): r == comp(a2) and pos-comp(a2)-comp(a1)
#'     exists.
#' }
#' An indel matches by the exact key only. Non-matches: uncovered (no row at
#' the position); unmatched_alt (REF agrees with r on its first base, or under
#' flexible matching an SNV whose ALT or complement is r, but no row carries
#' the pair); unmatched_ref (otherwise).
#'
#' For a non-palindromic SNV pair the pair and its complement together cover
#' A, C, G and T, so at a position carrying all three SNV rows the pair always
#' resolves through tiers 1 to 4; a tier is the lookup transformation, not a
#' verified strand.
#'
#' @param pos,ref,alt Input positions and upper-cased alleles (never modified).
#' @param is_snv Logical, TRUE for the SNV class.
#' @param fr_pos,fr_ref,fr_alt Every database row at the input positions.
#' @param match_method "exact" or "flexible".
#' @return list(tier = character, row = integer index into the database rows).
#' @keywords internal
#' @noRd
.favor_classify <- function(pos, ref, alt, is_snv, fr_pos, fr_ref, fr_alt, match_method) {
  n <- length(pos)
  tier <- rep(NA_character_, n)
  row  <- rep(NA_integer_, n)
  eq <- function(x, y) !is.na(x) & !is.na(y) & x == y   # element-wise, NA-safe

  fkey  <- paste(fr_pos, fr_ref, fr_alt, sep = "-")
  first <- !duplicated(fr_pos)
  r     <- substr(fr_ref[first], 1L, 1L)[match(pos, fr_pos[first])]  # NA: position not in FAVOR
  covered <- !is.na(r)
  look <- function(a1, a2) match(paste(pos, a1, a2, sep = "-"), fkey)

  # Tier 1: the key as given.
  hit <- look(ref, alt)
  ok  <- !is.na(hit)
  tier[ok] <- "exact"; row[ok] <- hit[ok]

  if (match_method == "flexible") {
    # Tier 2: FAVOR's REF is the input's ALT (SNVs only; an indel's reciprocal key
    # is a different variant).
    todo <- is.na(tier) & covered & is_snv & eq(alt, r)
    hit <- look(alt, ref)
    sel <- todo & !is.na(hit)
    tier[sel] <- "swapped"; row[sel] <- hit[sel]

    # Tiers 3 and 4: neither allele is FAVOR's REF, so complement both.
    cref <- .favor_complement(ref); calt <- .favor_complement(alt)
    todo <- is.na(tier) & covered & is_snv & !eq(ref, r) & !eq(alt, r)
    hit <- look(cref, calt)
    sel <- todo & eq(cref, r) & !is.na(hit)
    tier[sel] <- "flipped"; row[sel] <- hit[sel]
    hit <- look(calt, cref)
    sel <- todo & is.na(tier) & eq(calt, r) & !is.na(hit)
    tier[sel] <- "flipped_swapped"; row[sel] <- hit[sel]
  }

  # Non-matches carry a reason.
  rest <- is.na(tier)
  ref_ok <- eq(substr(ref, 1L, 1L), r)
  if (match_method == "flexible") {
    cref <- .favor_complement(ref); calt <- .favor_complement(alt)
    ref_ok <- ref_ok | (is_snv & (eq(alt, r) | eq(cref, r) | eq(calt, r)))
  }
  tier[rest & !covered] <- "uncovered"
  tier[rest & covered & ref_ok] <- "unmatched_alt"
  tier[rest & covered & !ref_ok] <- "unmatched_ref"
  list(tier = tier, row = row)
}


#' Aggregate the Database Rows at Each Position (Position-Only Keys)
#'
#' @description
#' "average": numeric features are the mean of their non-missing values at
#' the position (NA when none), string features the first non-empty value in
#' the database's row order. "first": every feature from the first row. The
#' caller passes the SNV rows at the position only (decision D10).
#'
#' @param pos Positions of the position-only inputs (all covered).
#' @param row_pos Position of each candidate database row.
#' @param ft Named list of feature vectors aligned with row_pos.
#' @param method "average" or "first".
#' @return Named list of feature vectors aligned with pos.
#' @keywords internal
#' @noRd
.favor_position_values <- function(pos, row_pos, ft, method) {
  upos <- unique(row_pos)
  grp  <- match(row_pos, upos)
  idx  <- match(pos, upos)
  first_row <- match(seq_along(upos), grp)
  out <- lapply(ft, function(v) {
    if (method == "average" && (is.numeric(v) || is.logical(v))) {
      v <- as.numeric(v)
      s <- rowsum(ifelse(is.na(v), 0, v), grp, reorder = TRUE)[, 1]
      k <- rowsum(as.numeric(!is.na(v)), grp, reorder = TRUE)[, 1]
      unname(ifelse(k > 0, s / k, NA_real_)[idx])
    } else if (method == "average") {
      ok <- !is.na(v) & nzchar(v)
      fne <- match(seq_along(upos), ifelse(ok, grp, NA_integer_))
      v[fne][idx]
    } else {
      v[first_row][idx]
    }
  })
  names(out) <- names(ft)
  out
}


#' Feature Values of Selected Database Rows
#'
#' @description CSV: the columns as read. Parquet: scalar columns as read, and
#'   the nested fields serialized into the v1 string layout. An empty string is
#'   returned as NA on both backends, because FAVOR v1 writes a missing string
#'   as an empty CSV field and FAVOR 2.0 as a null (the aGDS stores both as "").
#' @keywords internal
#' @noRd
.favor_feature_values <- function(fr, rows, features, backend) {
  out <- list()
  if (backend == "parquet") {
    nested <- names(.favor2_field_map()$nested)
    for (f in setdiff(features, nested)) out[[f]] <- fr[[f]][rows]
    out <- c(out, .favor2_serialize(fr, rows, intersect(features, nested)))
  } else {
    for (f in features) out[[f]] <- fr[[f]][rows]
  }
  lapply(out[features], function(v) {
    if (is.character(v)) v[!is.na(v) & !nzchar(v)] <- NA_character_
    v
  })
}


#' The String-Valued Default Features
#' @keywords internal
#' @noRd
.favor_string_features <- function() {
  c("genecode_comprehensive_category", "genecode_comprehensive_exonic_category",
    "genecode_comprehensive_info", "genecode_comprehensive_exonic_info",
    "metasvm_pred", "genehancer", "cage_tc", "cage_promoter", "rdhs", "rsid")
}


#' Choose the Database Backend
#'
#' @description "auto": Parquet when the directory holds chromosome_*.parquet
#'   and no chr*_*.csv chunks, CSV in the reverse case; an error when it holds
#'   both kinds or neither.
#' @keywords internal
#' @noRd
.favor_resolve_backend <- function(favor_db_path, favor_db_format) {
  if (favor_db_format != "auto") return(favor_db_format)
  has_pq  <- length(list.files(favor_db_path, pattern = "^chromosome_.*\\.parquet$")) > 0
  has_csv <- length(list.files(favor_db_path, pattern = "^chr[^_]+_[0-9]+\\.csv$")) > 0
  if (has_pq && has_csv) {
    stop("FAVOR database directory holds both Parquet and CSV chunk files: ",
         favor_db_path, ". Set favor_db_format explicitly.")
  }
  if (!has_pq && !has_csv) {
    stop("No FAVOR database files (chromosome_*.parquet or chr*_*.csv) in: ",
         favor_db_path)
  }
  if (has_pq) "parquet" else "csv"
}


#' Load the CSV Chunk Table (FAVORdatabase_chrsplit.csv)
#'
#' @description Priority: the user's file, then one in favor_db_path, then the
#'   package's bundled copy.
#' @keywords internal
#' @noRd
.favor_load_split_file <- function(favor_db_path, favor_split_file, verbose) {
  if (is.null(favor_split_file)) {
    favor_split_file <- file.path(favor_db_path, "FAVORdatabase_chrsplit.csv")
    if (!file.exists(favor_split_file)) {
      favor_split_file <- system.file("extdata", "FAVORdatabase_chrsplit.csv",
                                      package = "GLOWr")
      if (!nzchar(favor_split_file) || !file.exists(favor_split_file)) {
        stop("Could not find FAVORdatabase_chrsplit.csv in ", favor_db_path,
             " or in package extdata")
      }
      if (verbose >= 1) message("Using package bundled FAVORdatabase_chrsplit.csv")
    }
  } else if (!file.exists(favor_split_file)) {
    stop("Split file not found: ", favor_split_file)
  }
  split_data <- data.table::fread(favor_split_file, data.table = FALSE)
  required_cols <- c("Chr", "File_No", "Start_Pos", "End_Pos")
  if (!all(required_cols %in% names(split_data))) {
    stop(sprintf("Split file missing required columns: %s",
                 paste(setdiff(required_cols, names(split_data)), collapse = ", ")))
  }
  split_data
}


#################### CSV backend (FAVOR v1) ####################

#' Database Rows at the Given Positions of One Chromosome (CSV Backend)
#'
#' @description
#' Finds the chunks that cover the positions in the split table, and returns
#' every row of those chunks at the positions, with the key columns and the
#' requested features that the chunks carry. With xsv, the rows are fetched
#' by a streaming join on \code{position} (the positions file is the small,
#' hashed side); otherwise each chunk is read with fread() and filtered. A
#' chunk that the table lists but that is absent from disk is fatal (decision
#' D5); a chromosome absent from the table returns NULL.
#'
#' @return NULL, or list(rows = data.frame(position, ref_vcf, alt_vcf,
#'   features..., .rsid), files, served, missing). The column .rsid is the
#'   chunk's rsid (NA where the chunk has none), read for the rsID check
#'   whether or not rsid is a requested feature.
#' @keywords internal
#' @noRd
.favor_csv_rows <- function(chr, positions, favor_db_path, split_data, features,
                            use_xsv = TRUE, verbose = 1) {
  sd <- split_data[as.character(split_data$Chr) == chr, , drop = FALSE]
  if (!nrow(sd)) return(NULL)
  sd <- sd[sd$Start_Pos <= max(positions) & sd$End_Pos >= min(positions), , drop = FALSE]
  str_feats <- .favor_string_features()
  empty <- function(cols) {
    df <- data.frame(position = integer(0), ref_vcf = character(0), alt_vcf = character(0))
    for (f in cols) df[[f]] <- if (f %in% str_feats) character(0) else numeric(0)
    df$.rsid <- character(0)
    df
  }
  if (!nrow(sd)) {
    return(list(rows = empty(character(0)), files = character(0),
                served = character(0), missing = character(0)))
  }
  files <- file.path(favor_db_path, paste0("chr", sd$Chr, "_", sd$File_No, ".csv"))
  absent <- files[!file.exists(files)]
  if (length(absent)) {
    # Fatal: proceeding would silently leave every variant in the chunk's range
    # unannotated (the data-loss class of the 2026-09-01 multi-chunk bug).
    stop(sprintf(paste0("FAVOR chunk file not found: %s. All chunks covering the ",
                        "requested variants must be present."),
                 paste(absent, collapse = ", ")), call. = FALSE)
  }
  use_xsv <- isTRUE(use_xsv) && .check_xsv_available()
  pos_file <- NULL
  if (use_xsv) {
    pos_file <- tempfile("favor_positions_", fileext = ".csv")
    on.exit(unlink(pos_file), add = TRUE)
    # A distinct header, so the join output carries no duplicate column name.
    data.table::fwrite(data.frame(glowr_query_position = as.integer(positions)), pos_file)
  }
  out <- vector("list", length(files))
  served <- NULL
  for (i in seq_along(files)) {
    header <- names(data.table::fread(files[i], nrows = 0L, showProgress = FALSE))
    if (!all(c("position", "ref_vcf", "alt_vcf") %in% header)) {
      stop("FAVOR chunk lacks position/ref_vcf/alt_vcf columns: ", files[i], call. = FALSE)
    }
    have <- intersect(features, header)
    served <- if (is.null(served)) have else intersect(served, have)
    has_rs <- "rsid" %in% header   # the database rsID, for the rsID check
    cols <- unique(c("position", "ref_vcf", "alt_vcf", have, if (has_rs) "rsid"))
    cls <- list(integer = "position",
                character = unique(c("ref_vcf", "alt_vcf", intersect(have, str_feats),
                                     if (has_rs) "rsid")),
                numeric = intersect(setdiff(have, str_feats), .default_favor_features()))
    cls <- cls[lengths(cls) > 0]
    if (use_xsv) {
      joined <- tempfile("favor_join_", fileext = ".csv")
      stderr_f <- tempfile("favor_join_", fileext = ".stderr")
      status <- tryCatch(
        system2("xsv", args = c("join", "position", shQuote(files[i]),
                                "glowr_query_position", shQuote(pos_file)),
                stdout = joined, stderr = stderr_f),
        error = function(e) e$message)
      if (!identical(status, 0L)) {
        err_txt <- if (file.exists(stderr_f)) paste(readLines(stderr_f, warn = FALSE), collapse = " ") else ""
        unlink(c(joined, stderr_f))
        stop(sprintf("xsv join failed on FAVOR chunk %s (status %s)%s. Rerun with use_xsv = FALSE.",
                     basename(files[i]), paste(status, collapse = ","),
                     if (nzchar(err_txt)) paste0(": ", err_txt) else ""), call. = FALSE)
      }
      rows <- data.table::fread(joined, select = cols, colClasses = cls,
                                data.table = FALSE, showProgress = FALSE)
      unlink(c(joined, stderr_f))
    } else {
      rows <- data.table::fread(files[i], select = cols, colClasses = cls,
                                data.table = FALSE, showProgress = FALSE)
      rows <- rows[rows$position %in% positions, , drop = FALSE]
    }
    if (verbose >= 2) {
      message(sprintf("    %s: %d rows at the query positions", basename(files[i]), nrow(rows)))
    }
    rows$.rsid <- if (has_rs) as.character(rows$rsid) else rep(NA_character_, nrow(rows))
    out[[i]] <- rows
  }
  served <- as.character(served)
  rows <- data.table::rbindlist(lapply(out, function(d) d[, c("position", "ref_vcf", "alt_vcf", served, ".rsid"), drop = FALSE]),
                                use.names = TRUE)
  rows <- as.data.frame(rows)
  if (!nrow(rows)) rows <- empty(served)
  list(rows = rows, files = files, served = served, missing = setdiff(features, served))
}


#################### Parquet backend (FAVOR 2.0) ####################

#' FAVOR 2.0 Parquet Sources of Our Feature Names
#'
#' @description Scalar fields map to one Parquet leaf column, named by its
#'   dotted path. The three "nested" fields are read as their raw list or
#'   struct columns (all the leaves under `path`) and rebuilt into the v1 string
#'   layout by \code{.favor2_serialize()}. "absent" lists the FAVOR v1 fields
#'   that FAVOR 2.0 does not carry (its apc group holds only the thirteen
#'   versions mapped here). Plan:
#' @keywords internal
#' @noRd
.favor2_field_map <- function() {
  apc13 <- c("conservation_v2", "epigenetics", "epigenetics_active",
             "epigenetics_repressed", "epigenetics_transcription",
             "local_nucleotide_diversity_v3", "mappability", "micro_rna",
             "mutation_density", "protein_function_v3", "proximity_to_coding_v2",
             "proximity_to_tsstes", "transcription_factor")
  scalar <- c(
    cadd_phred = "main.cadd.phred",
    linsight   = "linsight",
    fathmm_xf  = "fathmm_xf",
    stats::setNames(paste0("apc.", apc13), paste0("apc_", apc13)),
    genecode_comprehensive_category        = "gencode.region_type",
    genecode_comprehensive_exonic_category = "gencode.consequence",
    metasvm_pred  = "dbnsfp.metasvm_pred",
    rsid          = "dbsnp.rsid",
    cage_tc       = "cage.cage_tc",
    cage_promoter = "cage.cage_promoter",
    rdhs          = "ccre.ids")
  # Each nested field: the raw columns it needs, as (internal column name ->
  # dotted path), and the leaves under each path.
  tx <- c("gene", "transcript_id", "location", "hgvsc", "hgvsp")
  nested <- list(
    genecode_comprehensive_info = list(
      .gi_genes = list(path = "gencode.genes", leaves = "gencode.genes")),
    genecode_comprehensive_exonic_info = list(
      .gei_tx = list(path = "gencode.transcripts", leaves = paste0("gencode.transcripts.", tx))),
    genehancer = list(
      .gh_id = list(path = "genehancer.id", leaves = "genehancer.id"),
      .gh_score = list(path = "genehancer.feature_score", leaves = "genehancer.feature_score"),
      .gh_targets = list(path = "genehancer.targets",
                         leaves = c("genehancer.targets.gene", "genehancer.targets.score"))))
  absent <- c("apc_conservation", "apc_local_nucleotide_diversity",
              "apc_local_nucleotide_diversity_v2", "apc_proximity_to_coding")
  list(scalar = scalar, nested = nested, absent = absent)
}


#' Database Rows at the Given Positions of One Chromosome (Parquet Backend)
#'
#' @description
#' Opens chromosome_<chr>.parquet alone (the 24 files do not share one schema:
#' `chromosome` is int64 in 1-22 and a string in X and Y) and walks its row
#' groups. For each group it reads the \code{position} leaf, keeps the rows at
#' the requested positions, and only then reads the key columns and the leaf
#' columns of the requested features for those rows. Reading leaves, not
#' whole struct columns, matters: a feature such as \code{main.cadd.phred}
#' sits inside a struct of about a hundred leaves. Memory is therefore bounded
#' by one row group, whatever the chromosome's size or the cohort's. Every row
#' at a position is returned, including rows of a position that spans two row
#' groups. A missing file for chromosome 1-22, X or Y is fatal (decision D5);
#' for any other chromosome it returns NULL.
#'
#' @return NULL, or list(rows, files, served, missing). The rows carry .rsid,
#'   the file's dbsnp.rsid (NA where the file has no such leaf), read for the
#'   rsID check whether or not rsid is a requested feature.
#' @keywords internal
#' @noRd
.favor_parquet_rows <- function(chr, positions, favor_db_path, features, verbose = 1) {
  f <- file.path(favor_db_path, paste0("chromosome_", chr, ".parquet"))
  if (!file.exists(f)) {
    if (chr %in% c(as.character(1:22), "X", "Y")) {
      stop("FAVOR 2.0 file not found: ", f, call. = FALSE)
    }
    return(NULL)
  }
  fm <- .favor2_field_map()
  scalar <- fm$scalar[intersect(names(fm$scalar), features)]
  nested_fields <- fm$nested[intersect(names(fm$nested), features)]

  pr <- arrow::ParquetFileReader$create(f)
  leaves <- .favor2_leaf_paths(pr$GetSchema())
  # The key columns are mandatory. A requested feature whose leaf columns the
  # file lacks is skipped and reported, like a feature absent from a CSV chunk.
  key_leaves <- c("position", "ref_vcf", "alt_vcf")
  if (!all(key_leaves %in% leaves)) {
    stop("FAVOR 2.0 file ", basename(f), " lacks the key columns: ",
         paste(setdiff(key_leaves, leaves), collapse = ", "), call. = FALSE)
  }
  lack_scalar <- names(scalar)[!unname(scalar) %in% leaves]
  lack_nested <- names(nested_fields)[!vapply(nested_fields, function(parts)
    all(unlist(lapply(parts, `[[`, "leaves"), use.names = FALSE) %in% leaves), TRUE)]
  if (length(c(lack_scalar, lack_nested)) && verbose >= 1) {
    message(sprintf("  FAVOR 2.0 file %s lacks the columns of %s; skipped",
                    basename(f), paste(c(lack_scalar, lack_nested), collapse = ", ")))
  }
  scalar <- scalar[setdiff(names(scalar), lack_scalar)]
  nested_fields <- nested_fields[setdiff(names(nested_fields), lack_nested)]
  nested <- unlist(unname(nested_fields), recursive = FALSE)
  served <- intersect(features, c(names(scalar), names(nested_fields)))
  rs_leaf <- if ("dbsnp.rsid" %in% leaves) "dbsnp.rsid" else NULL   # for the rsID check
  want <- unique(c(key_leaves, unname(scalar),
                   unlist(lapply(nested, `[[`, "leaves"), use.names = FALSE), rs_leaf))
  idx <- match(want, leaves)
  pos_leaf <- match("position", leaves) - 1L
  positions <- unique(as.integer(positions))

  parts <- list()
  n_groups <- pr$num_row_groups
  for (g in seq_len(n_groups) - 1L) {
    # A file reader accumulates memory over the row groups it has read, so it is
    # reopened every few groups (opening costs one footer read).
    if (g > 0L && g %% 8L == 0L) {
      rm(pr); invisible(gc(verbose = FALSE))
      pr <- arrow::ParquetFileReader$create(f)
    }
    pg <- as.vector(pr$ReadRowGroup(g, pos_leaf)$position)
    keep <- which(pg %in% positions)
    if (!length(keep)) next
    d <- as.data.frame(pr$ReadRowGroup(g, idx - 1L)$Take(keep - 1L))
    part <- data.frame(position = as.integer(d$position), ref_vcf = d$ref_vcf,
                       alt_vcf = d$alt_vcf, stringsAsFactors = FALSE)
    for (nm in names(scalar)) part[[nm]] <- .favor2_get_path(d, scalar[[nm]])
    for (nm in names(nested)) part[[nm]] <- .favor2_get_path(d, nested[[nm]]$path)
    part$.rsid <- if (is.null(rs_leaf)) rep(NA_character_, length(keep)) else
      as.character(.favor2_get_path(d, rs_leaf))
    parts[[length(parts) + 1L]] <- part
    if (verbose >= 2) {
      message(sprintf("    %s row group %d: %d rows at the query positions", basename(f), g, length(keep)))
    }
    rm(d)
  }
  rows <- .favor2_bind_parts(parts, c("position", "ref_vcf", "alt_vcf", names(scalar), names(nested), ".rsid"))
  list(rows = rows, files = f, served = served, missing = setdiff(features, served))
}


#' Dotted Paths of the Parquet Leaf Columns, in Leaf-Index Order
#'
#' @description Walks an arrow schema depth first. A struct contributes its
#'   children, a list of structs the children of its element, and anything else
#'   one leaf. The position of a path in the result is its Parquet leaf column
#'   index (0-based after subtracting one), which ReadRowGroup() accepts.
#' @keywords internal
#' @noRd
.favor2_leaf_paths <- function(schema) {
  walk <- function(field, prefix) {
    t <- field$type
    nm <- paste0(prefix, field$name)
    if (inherits(t, "StructType")) {
      unlist(lapply(t$fields(), walk, prefix = paste0(nm, ".")))
    } else if (inherits(t, c("ListType", "LargeListType")) &&
               inherits(t$value_type, "StructType")) {
      unlist(lapply(t$value_type$fields(), walk, prefix = paste0(nm, ".")))
    } else {
      nm
    }
  }
  unlist(lapply(schema$fields, walk, prefix = ""))
}


#' Extract a Column by Dotted Path from a Converted (Nested) data.frame
#' @keywords internal
#' @noRd
.favor2_get_path <- function(d, path) {
  x <- d
  for (p in strsplit(path, ".", fixed = TRUE)[[1]]) x <- x[[p]]
  if (is.list(x) && !is.data.frame(x)) x <- unclass(as.list(x))  # arrow_list -> plain list
  if (is.factor(x)) x <- as.character(x)
  x
}


#' Bind the Per-Row-Group Parts Column by Column (List Columns Included)
#' @keywords internal
#' @noRd
.favor2_bind_parts <- function(parts, cols) {
  if (!length(parts)) {
    out <- data.frame(position = integer(0), ref_vcf = character(0), alt_vcf = character(0))
    for (nm in setdiff(cols, names(out))) out[[nm]] <-
      if (nm == ".rsid") character(0) else if (startsWith(nm, ".")) list() else NA[0]
    return(out)
  }
  out <- data.frame(row.names = seq_len(sum(vapply(parts, nrow, 0L))))
  for (nm in cols) {
    v <- lapply(parts, `[[`, nm)
    out[[nm]] <- if (is.list(v[[1]])) do.call(c, v) else unlist(v, use.names = FALSE)
  }
  out
}


#' Rebuild the v1 String Layout of the Nested FAVOR 2.0 Fields
#'
#' @description Serialization contracts of plan 11 section 4.2:
#' \itemize{
#'   \item genehancer: "Name=<score>;genehancer_id=<id>" then
#'     ";connected_gene=<gene>;score=<score>" per target in FAVOR 2.0's order,
#'     numbers to two decimals (STAARpipeline reads the first connected gene).
#'   \item genecode_comprehensive_exonic_info: "gene:transcript:exon:c.:p." per
#'     transcript, comma-joined with v1's trailing comma; a list whose
#'     transcripts carry nothing but a gene name (v1's "UNKNOWN" sentinel) is
#'     written as the gene name alone.
#'   \item genecode_comprehensive_info: gene names comma-joined in FAVOR 2.0's
#'     order, dropping list elements that contain ":" (fragments of ANNOVAR's
#'     UTR transcript strings).
#' }
#' @param fr Rows from .favor_parquet_rows().
#' @param rows Integer indices of the rows to serialize.
#' @param features Which of the three nested fields to build.
#' @return Named list of character vectors aligned with rows.
#' @keywords internal
#' @noRd
.favor2_serialize <- function(fr, rows, features) {
  res <- list()
  blank <- function(x) is.na(x) | !nzchar(x)
  fmt2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
  if ("genehancer" %in% features) {
    id <- fr$.gh_id[rows]; sc <- fr$.gh_score[rows]; tg <- fr$.gh_targets[rows]
    res$genehancer <- vapply(seq_along(rows), function(k) {
      if (blank(id[k])) return(NA_character_)
      t <- tg[[k]]
      targets <- if (is.null(t) || !NROW(t)) "" else
        paste0(";connected_gene=", t$gene, ";score=", fmt2(t$score), collapse = "")
      paste0("Name=", fmt2(sc[k]), ";genehancer_id=", id[k], targets)
    }, character(1))
  }
  if ("genecode_comprehensive_exonic_info" %in% features) {
    tx <- fr$.gei_tx[rows]
    res$genecode_comprehensive_exonic_info <- vapply(tx, function(t) {
      if (is.null(t) || !NROW(t)) return(NA_character_)
      detail <- !(blank(t$transcript_id) & blank(t$location) & blank(t$hgvsc) & blank(t$hgvsp))
      if (!any(detail)) {
        g <- unique(t$gene[!blank(t$gene)])
        return(if (length(g)) paste(g, collapse = ",") else NA_character_)
      }
      paste0(paste(t$gene, t$transcript_id, t$location, t$hgvsc, t$hgvsp,
                   sep = ":", collapse = ","), ",")
    }, character(1))
  }
  if ("genecode_comprehensive_info" %in% features) {
    g <- fr$.gi_genes[rows]
    res$genecode_comprehensive_info <- vapply(g, function(x) {
      if (is.null(x) || !length(x)) return(NA_character_)
      x <- unique(x[!blank(x) & !grepl(":", x, fixed = TRUE)])
      if (!length(x)) NA_character_ else paste(x, collapse = ",")
    }, character(1))
  }
  lapply(res, unname)
}


#################### Provenance ####################

#' Provenance Record for the aGDS (Plan 11 Section 4.6)
#'
#' @description Database format and directory, the files read with byte
#'   counts and (when a SHA256SUMS file sits in the database directory or its
#'   parent) their SHA-256, the FAVOR release or "unknown", the GLOWr version,
#'   the matching settings, the rsID policy and source (decision D9), the rows
#'   the position-only aggregation uses (decision D10), the outcome and check
#'   counts and the date.
#' @return Named list of character or numeric vectors (GDS attributes).
#' @keywords internal
#' @noRd
.favor_provenance <- function(backend, favor_db_path, files, favor_release,
                              match_method, na_allele_method, na_handling,
                              tiers, features_missing, rsid_policy = "record",
                              rsid_source = "none", rsid_checks = NULL) {
  files <- as.character(files)
  sha <- rep("not computed", length(files))
  for (d in unique(c(favor_db_path, dirname(normalizePath(favor_db_path, mustWork = FALSE))))) {
    sf <- file.path(d, "SHA256SUMS")
    if (length(files) && file.exists(sf)) {
      lines <- readLines(sf, warn = FALSE)
      hash <- sub("^([0-9a-f]{64}).*$", "\\1", lines)
      name <- basename(sub("^[0-9a-f]{64}[ *]+", "", lines))
      m <- match(basename(files), name)
      sha[!is.na(m)] <- hash[m[!is.na(m)]]
    }
  }
  tab <- table(factor(tiers, levels = .favor_tier_levels()))
  rtab <- table(factor(rsid_checks, levels = .favor_rsid_check_levels()))
  list(
    favor_db_format        = backend,
    favor_db_path          = normalizePath(favor_db_path, mustWork = FALSE),
    favor_db_files         = if (length(files)) basename(files) else "none",
    favor_db_file_bytes    = if (length(files)) as.numeric(file.size(files)) else 0,
    favor_db_file_sha256   = if (length(files)) sha else "none",
    favor_release          = if (is.null(favor_release)) "unknown" else favor_release,
    favor_features_missing = if (length(features_missing)) features_missing else "none",
    glowr_version          = as.character(utils::packageVersion("GLOWr")),
    favor_match_method     = match_method,
    favor_position_only_method = na_allele_method,
    favor_position_only_rows = "snv",
    favor_na_handling      = na_handling,
    favor_rsid_policy      = rsid_policy,
    favor_rsid_source      = rsid_source,
    favor_tier_names       = names(tab),
    favor_tier_counts      = as.integer(tab),
    favor_rsid_check_names = names(rtab),
    favor_rsid_check_counts = as.integer(rtab),
    annotation_date        = format(Sys.Date(), "%Y-%m-%d")
  )
}

#################### aGDS writers ####################

#' Update GDS File with FAVOR Annotations
#'
#' @description
#' Writes annotation data back to an existing GDS file, converting it to
#' aGDS format by adding a FunctionalAnnotation node.
#'
#' @param gds_path Character. Path to GDS file (will be modified)
#' @param annotations data.frame with VarInfo (and, for GDS input, variant_id),
#'   the feature columns and match_tier
#' @param features Character vector of feature column names
#' @param provenance Named list of provenance attributes, or NULL
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (side effect: modifies GDS file)
#'
#' @details
#' Opens the GDS file in read-write mode and adds
#' \code{/annotation/info/FunctionalAnnotation} as a \emph{folder}
#' (\code{addfolder.gdsn}) holding one native-typed sub-node per feature
#' (numeric features stay numeric, string features stay character), plus
#' \code{/annotation/info/favor_match_tier}. Rows are aligned to the GDS
#' through \code{variant.id}, so duplicate keys and reordered or shortened
#' annotation tables still land on their own variants. If annotation nodes
#' already exist, they are overwritten with a warning.
#'
#' @keywords internal
#' @noRd
.update_gds_with_annotations <- function(gds_path, annotations, features,
                                         provenance = NULL, verbose = 1) {

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

  # Align the annotation rows to the GDS variant order before writing. The
  # in-place update writes columns verbatim, so a row-order mismatch (or a
  # shortened table from na_handling = "drop") would attach each variant its
  # neighbour's annotations. variant.id is the immutable identity; the key
  # string is only a fallback for tables that carry no variant_id.
  vid <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "variant.id"))
  if ("variant_id" %in% names(annotations)) {
    annot_idx <- match(vid, annotations$variant_id)
  } else {
    chrom <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "chromosome"))
    pos   <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "position"))
    alle  <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "allele"))
    ref   <- sub(",.*$", "", alle)
    alt   <- ifelse(grepl(",", alle, fixed = TRUE), sub("^[^,]*,", "", alle), "")
    annot_idx <- match(paste(chrom, pos, ref, alt, sep = "-"), annotations$VarInfo)
  }
  annotations <- annotations[annot_idx, , drop = FALSE]

  # Write FunctionalAnnotation as a STAARpipeline-style sub-node folder (one
  # typed sub-node per feature), then the matching outcome beside it.
  .write_functional_annotation_subnodes(
    parent_node = info_folder,
    annotations = annotations,
    features = features,
    provenance = provenance,
    verbose = verbose
  )
  .write_match_tier_node(info_folder, annotations$match_tier, verbose = verbose)
  .write_rsid_check_node(info_folder, annotations$rsid_check, verbose = verbose)

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
#' @param annotations data.frame with VarInfo, variant_id, the feature columns
#'   and match_tier
#' @param features Character vector of feature column names
#' @param provenance Named list of provenance attributes, or NULL
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (side effect: creates aGDS file)
#'
#' @details
#' Follows STAARpipeline's gds2agds.R pattern:
#' \enumerate{
#'   \item Copy input GDS file to output path
#'   \item Open copy with seqOpen(readonly = FALSE)
#'   \item Add /annotation/info/FunctionalAnnotation and
#'     /annotation/info/favor_match_tier
#'   \item Close file
#' }
#'
#' \strong{Annotation node format:} \code{/annotation/info/FunctionalAnnotation}
#' is created as a \emph{folder} (\code{addfolder.gdsn}) holding one sub-node per
#' feature, each carrying that feature's native-typed vector in GDS variant order
#' (numeric features stay numeric, string features such as
#' \code{genecode_comprehensive_category} stay character). This is the
#' STAARpipeline sub-node layout that the per-feature variant-set scan reads via
#' \code{seqGetData(gds, ".../FunctionalAnnotation/<feature>")}. Rows are aligned
#' through \code{variant.id}. The genotypes, sample IDs and variant fields of
#' the copy are the input's, untouched.
#'
#' @keywords internal
#' @noRd
.create_agds_from_gds <- function(input_gds_path, output_path, annotations,
                                   features, provenance = NULL, verbose = 1) {

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
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  copy_success <- file.copy(input_gds_path, output_path, overwrite = TRUE)
  if (!copy_success) {
    stop("Failed to copy GDS file to: ", output_path)
  }

  # Step 2: Open the copy for writing
  gds <- SeqArray::seqOpen(output_path, readonly = FALSE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Step 3: Align the annotation rows to the GDS variant order through
  # variant.id (the key string is only a fallback for tables without it).
  vid <- SeqArray::seqGetData(gds, "variant.id")
  if ("variant_id" %in% names(annotations)) {
    annot_idx <- match(vid, annotations$variant_id)
  } else {
    gds_key <- paste(SeqArray::seqGetData(gds, "chromosome"),
                     SeqArray::seqGetData(gds, "position"),
                     SeqArray::seqGetData(gds, "$ref"),
                     SeqArray::seqGetData(gds, "$alt"), sep = "-")
    annot_idx <- match(gds_key, annotations$VarInfo)
  }
  # Reordering the data.frame (not coercing to a matrix) preserves each
  # feature's native type.
  annot_ordered <- annotations[annot_idx, , drop = FALSE]

  # Step 4: Navigate to or create annotation/info, then write the
  # FunctionalAnnotation sub-node folder and the matching outcome.
  anno_folder <- tryCatch(
    gdsfmt::index.gdsn(gds, "annotation/info"),
    error = function(e) NULL
  )
  if (is.null(anno_folder)) {
    anno_root <- tryCatch(
      gdsfmt::index.gdsn(gds, "annotation"),
      error = function(e) NULL
    )
    if (is.null(anno_root)) {
      anno_root <- gdsfmt::addfolder.gdsn(gds, "annotation")
    }
    anno_folder <- gdsfmt::addfolder.gdsn(anno_root, "info")
  }

  .write_functional_annotation_subnodes(
    parent_node = anno_folder,
    annotations = annot_ordered,
    features = features,
    provenance = provenance,
    verbose = verbose
  )
  .write_match_tier_node(anno_folder, annot_ordered$match_tier, verbose = verbose)
  .write_rsid_check_node(anno_folder, annot_ordered$rsid_check, verbose = verbose)

  if (verbose >= 1) {
    n_annotated <- sum(!is.na(annot_idx))
    message(sprintf("  Created aGDS: %d variants, %d annotated, %d features",
                    length(vid), n_annotated, length(features)))
  }

  return(invisible(NULL))
}


#' Write Annotated Variants to aGDS Format
#'
#' @description
#' Writes annotated variants to an annotation-only GDS file (non-GDS input).
#'
#' @param data data.frame with VarInfo, the annotation feature columns and
#'   match_tier
#' @param output_path Character. Path for output aGDS file
#' @param features Character vector of feature column names
#' @param provenance Named list of provenance attributes, or NULL
#' @param verbose Integer. Verbosity level
#'
#' @return NULL (side effect: creates aGDS file)
#'
#' @details
#' \strong{aGDS Structure:}
#' \itemize{
#'   \item /chromosome, /position, /ref, /alt, /VarInfo: parsed from VarInfo
#'   \item /annotation/info/FunctionalAnnotation: a \emph{folder} holding one
#'     native-typed sub-node per feature, in input row order
#'   \item /annotation/info/favor_match_tier: the matching outcome per row
#' }
#'
#' @keywords internal
#' @noRd
.write_agds <- function(data, output_path, features, provenance = NULL, verbose = 1) {

  # Check if gdsfmt package is available
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required for aGDS output. Install with: BiocManager::install('gdsfmt')")
  }

  # Parse VarInfo to extract components
  varinfo_split <- strsplit(data$VarInfo, "-", fixed = TRUE)

  chr_vec <- sapply(varinfo_split, function(x) if (length(x) >= 1) x[1] else NA_character_)
  pos_vec <- suppressWarnings(as.integer(sapply(varinfo_split, function(x) if (length(x) >= 2) x[2] else NA_character_)))
  ref_vec <- sapply(varinfo_split, function(x) if (length(x) >= 3) x[3] else NA_character_)
  alt_vec <- sapply(varinfo_split, function(x) if (length(x) >= 4) x[4] else NA_character_)
  na_to_blank <- function(x) ifelse(is.na(x), "", x)

  # Create GDS file
  gds_file <- gdsfmt::createfn.gds(output_path)

  tryCatch({
    gdsfmt::add.gdsn(gds_file, "chromosome", na_to_blank(chr_vec), compress = "LZMA_RA", closezip = TRUE)
    gdsfmt::add.gdsn(gds_file, "position", pos_vec, compress = "LZMA_RA", closezip = TRUE)
    gdsfmt::add.gdsn(gds_file, "ref", na_to_blank(ref_vec), compress = "LZMA_RA", closezip = TRUE)
    gdsfmt::add.gdsn(gds_file, "alt", na_to_blank(alt_vec), compress = "LZMA_RA", closezip = TRUE)
    gdsfmt::add.gdsn(gds_file, "VarInfo", data$VarInfo, compress = "LZMA_RA", closezip = TRUE)

    # Create annotation group
    annot_group <- gdsfmt::addfolder.gdsn(gds_file, "annotation")
    info_group <- gdsfmt::addfolder.gdsn(annot_group, "info")

    # `data` is already in the intended row order here, so no reordering.
    .write_functional_annotation_subnodes(
      parent_node = info_group,
      annotations = data,
      features = features,
      provenance = provenance,
      verbose = verbose
    )
    .write_match_tier_node(info_group, data$match_tier, verbose = verbose)
    .write_rsid_check_node(info_group, data$rsid_check, verbose = verbose)

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
#' @param provenance Named list written as attributes of the folder, or NULL.
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
#' \strong{Native type preservation:} numeric features are written as
#' float64, character features as strings with a missing value stored as
#' \code{""} (GDS strings have no NA, and this is how the existing trees store
#' it).
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
                                                   features, provenance = NULL,
                                                   verbose = 1) {

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
    if (is.character(vec)) vec[is.na(vec)] <- ""
    if (is.logical(vec)) vec <- as.numeric(vec)   # an untyped all-NA column
    gdsfmt::add.gdsn(fa_folder, feat, val = vec,
                     compress = "LZMA_RA", closezip = TRUE)
  }

  # Record feature_names on the folder for discoverability (sub-node names are
  # authoritative), then the provenance of the annotation.
  gdsfmt::put.attr.gdsn(fa_folder, "feature_names", features)
  for (nm in names(provenance)) {
    gdsfmt::put.attr.gdsn(fa_folder, nm, provenance[[nm]])
  }

  invisible(fa_folder)
}


#' Write the Matching Outcome per Variant beside FunctionalAnnotation
#'
#' @description Adds (or overwrites) \code{favor_match_tier} under
#'   \code{parent_node} (typically \code{annotation/info}). A variant that was
#'   not annotated in this call (for example outside a variant filter) is
#'   stored as \code{""}.
#' @keywords internal
#' @noRd
.write_match_tier_node <- function(parent_node, tier, verbose = 1) {
  .write_outcome_node(parent_node, "favor_match_tier", tier,
                      "GLOWr::annotate_favor() matching outcome per variant")
}


#' Write the rsID Check per Variant beside FunctionalAnnotation (Decision D9)
#'
#' @description Adds (or overwrites) \code{favor_rsid_check} under
#'   \code{parent_node}: same, differs, chip_none or favor_none for every
#'   matched variant, \code{""} for the rest.
#' @keywords internal
#' @noRd
.write_rsid_check_node <- function(parent_node, check, verbose = 1) {
  .write_outcome_node(parent_node, "favor_rsid_check", check,
                      "GLOWr::annotate_favor() rsID check per matched variant (same, differs, chip_none, favor_none)")
}


#' Write One Per-Variant Character Outcome Node (Shared by the Two Above)
#' @keywords internal
#' @noRd
.write_outcome_node <- function(parent_node, name, values, description) {
  if (is.null(values)) return(invisible(NULL))
  if (name %in% gdsfmt::ls.gdsn(parent_node)) {
    gdsfmt::delete.gdsn(gdsfmt::index.gdsn(parent_node, name), force = TRUE)
  }
  values <- as.character(values)
  values[is.na(values)] <- ""
  node <- gdsfmt::add.gdsn(parent_node, name, val = values,
                           compress = "LZMA_RA", closezip = TRUE)
  gdsfmt::put.attr.gdsn(node, "Description", description)
  invisible(node)
}

#################### Input extraction and shared helpers ####################

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
#' @return data.frame with VarInfo column, variant_id for tracking and, when
#'   the GDS carries \code{annotation/id}, rsID (its value per variant; the
#'   rsID on a chip cohort), the input side of the rsID check.
#'
#' @details
#' Requires SeqArray package. Opens GDS file, applies optional filter,
#' extracts chromosome, position, ref, alt, and creates VarInfo string
#' (the whole ALT list for a multi-allelic site).
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

  # REF and ALT; a multi-allelic site's ALT is a comma-separated list
  ref <- SeqArray::seqGetData(gds, "$ref")
  alt <- SeqArray::seqGetData(gds, "$alt")

  # Get variant IDs for tracking (the writers align rows through them)
  variant_id <- SeqArray::seqGetData(gds, "variant.id")

  # Create VarInfo: CHR-POS-REF-ALT. A multi-allelic site keeps its whole ALT
  # list, which the annotator classifies as unsupported rather than annotating
  # the first ALT alone.
  VarInfo <- paste(chr, pos, ref, alt, sep = "-")

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
  # The GDS variant identifier (annotation/id): the rsID on a chip cohort, kept
  # as the evidence column of the rsID check. Absent node -> no column.
  ids <- tryCatch(SeqArray::seqGetData(gds, "annotation/id"), error = function(e) NULL)
  if (!is.null(ids) && length(ids) == nrow(result)) result$rsID <- as.character(ids)

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


#' Default FAVOR Annotation Features
#'
#' @description
#' Returns the complete annotation content of the FAVOR Essential DB, 30
#' fields: 17 Annotation Principal Components (all versions and sub-features),
#' 3 integrative scores, and 10 categorical fields. Used as the default for
#' \code{annotate_favor()} to annotate data broadly for downstream flexibility.
#'
#' For PI model training, use \code{.default_PI_features()} from
#' \code{get_PI_train.R}, which is a curated 16-feature subset.
#'
#' @return Character vector of 30 feature names.
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
    "fathmm_xf",
    # 10 categorical / string fields of the FAVOR Essential DB. Read by the
    # built-in coding masks (GENCODE category + exonic category, MetaSVM), by
    # the gene-centric noncoding masks (GENCODE info, CAGE, rDHS, GeneHancer)
    # and for variant identification (rsid). Written as native character
    # sub-nodes; never averaged.
    "genecode_comprehensive_category",
    "genecode_comprehensive_exonic_category",
    "genecode_comprehensive_info",
    "genecode_comprehensive_exonic_info",
    "metasvm_pred",
    "genehancer",
    "cage_tc",
    "cage_promoter",
    "rdhs",
    "rsid"
  )
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
    # Replace NA with 0 -- numeric features only. Character features (e.g. the
    # GENCODE category) keep NA: a "0" string would be a fabricated category.
    for (feat in features) {
      if (feat %in% names(data)) {
        if (!(is.numeric(data[[feat]]) || is.logical(data[[feat]]))) {
          if (verbose >= 2) {
            message(sprintf("    Keeping NA in non-numeric feature %s (na_handling='zero' applies to numeric features)", feat))
          }
          next
        }
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
#' The CSV backend uses it to fetch the database rows at the input positions
#' by a streaming join, without reading whole chunks into memory.
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
