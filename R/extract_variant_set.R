# ==============================================================================
# Variant Set Extraction from GDS/aGDS Files
# ==============================================================================
#
# Core function for extracting filtered variant sets from GDS files. Applies
# QC, variant type, positional, annotation-based, and MAF/MAC filters to
# produce ready-to-use genotype matrices and variant metadata. Also provides
# a chromosome-wide annotation median computation drawn from the same filter
# pool, used as a fallback reference for PI NA imputation.
#
# EXPORTED FUNCTIONS:
#   - extract_variant_set()           Extract filtered variant set from GDS
#   - compute_annotation_medians()    Compute chromosome-wide annotation medians
#
# INTERNAL HELPERS:
#   - .select_variants_by_spec()       Apply filter_spec to a candidate pool
#                                      (shared by the two exported functions)
#   - .normalize_region()              Normalize region input
#   - .evaluate_annotation_clauses()   Evaluate DNF annotation filter
#   - .resolve_annotation_path()       Map annotation name to GDS path
#   - .extract_annotation_scores()     Read numeric annotation columns
#   - .impute_geno_vset()              Impute missing genotypes + compute MAF/MAC


#################### EXPORTED MAIN FUNCTIONS ####################

#' Extract Filtered Variant Set from GDS
#'
#' Reads genotype data from a GDS/aGDS file, applies QC, variant type,
#' positional, and annotation-based filters, and returns a prepared
#' variant set for downstream processing and testing.
#'
#' @param gds An open \code{SeqVarGDSClass} object or character file path.
#'   If a path, the file is opened and closed internally.
#' @param region A list or single-row data.frame with elements/columns:
#'   chr, start, end. Optional: region_id, label.
#' @param filter_spec A \code{glow_variant_filter} from \code{variant_filter()}
#'   or \code{coding_filter()}.
#' @param sample_id Character vector of sample IDs for extraction and ordering.
#'   Must match a subset of GDS sample.id. If NULL, all samples are used
#'   in GDS order.
#' @param annotation_names Character vector of annotation names to extract
#'   for PI computation (e.g., c("CADD", "LINSIGHT")). Uses the annotation
#'   name catalog to resolve GDS paths. NULL means no annotations extracted.
#' @param Annotation_dir Character path prefix for annotations in GDS.
#'   Default "annotation/info/FunctionalAnnotation".
#' @param Annotation_name_catalog Data frame mapping annotation names to
#'   GDS sub-paths. NULL uses the bundled catalog.
#' @param impute_method Character: "mean" or "minor". Default "mean".
#'   \code{"mean"} replaces NA with the column mean (2 * AF), the expected
#'   dosage under Hardy-Weinberg. \code{"minor"} replaces NA with 0 when
#'   AF <= 0.5, or 2 when AF > 0.5; in both cases the imputed value assumes
#'   the individual carries zero copies of the minor allele. This matches
#'   the legacy STAARpipeline \code{matrix_flip_minor()} behavior.
#' @param verbose Integer verbosity level.
#'
#' @return A \code{glow_variant_set} S3 object, or NULL if the region
#'   has fewer than \code{filter_spec$min_variants} qualifying variants.
#'
#' @examples
#' \dontrun{
#' library(SeqArray)
#' gds <- seqOpen("chr22_essentialdb.gds")
#'
#' # Extract plof variants for a gene
#' region <- list(chr = "22", start = 17565844L, end = 17596578L,
#'                label = "CECR1")
#' spec <- coding_filter("plof")
#' vset <- extract_variant_set(gds, region, spec,
#'                              sample_id = sample_ids,
#'                              annotation_names = c("CADD", "LINSIGHT"))
#' dim(vset$G)          # n_samples x n_variants
#' vset$variant_info     # variant metadata
#' vset$annotations      # annotation score matrix
#'
#' seqClose(gds)
#' }
#'
#' @export
extract_variant_set <- function(gds,
                                 region,
                                 filter_spec,
                                 sample_id = NULL,
                                 annotation_names = NULL,
                                 Annotation_dir = "annotation/info/FunctionalAnnotation",
                                 Annotation_name_catalog = NULL,
                                 impute_method = "mean",
                                 verbose = 1) {

  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required. Install with: BiocManager::install('SeqArray')")
  }

  # ---- 1. Handle GDS connection ----
  opened_here <- FALSE
  if (is.character(gds)) {
    stopifnot(file.exists(gds))
    gds <- SeqArray::seqOpen(gds, readonly = TRUE)
    opened_here <- TRUE
    on.exit(SeqArray::seqClose(gds), add = TRUE)
  }

  # ---- 2. Normalize region ----
  region <- .normalize_region(region)

  # ---- 3. Resolve annotation name catalog ----
  if (is.null(Annotation_name_catalog)) {
    Annotation_name_catalog <- get("annotation_name_catalog",
                                    envir = asNamespace("GLOWr"))
  }

  # ---- 4. Read chromosome-level position and variant.id ----
  SeqArray::seqResetFilter(gds, verbose = FALSE)
  all_positions <- SeqArray::seqGetData(gds, "position")
  all_variant_ids <- SeqArray::seqGetData(gds, "variant.id")
  all_chr <- SeqArray::seqGetData(gds, "chromosome")

  # ---- 5. Positional filter ----
  chr_match <- as.character(all_chr) == as.character(region$chr)
  pos_match <- (all_positions >= region$start) & (all_positions <= region$end)
  in_region <- chr_match & pos_match
  n_total_in_region <- sum(in_region)

  if (n_total_in_region == 0) {
    if (verbose >= 1) message("No variants in region ", region$label)
    return(NULL)
  }

  region_variant_ids <- all_variant_ids[in_region]

  # ---- 6. Apply filter spec (QC + type + annotation clauses + MAF/MAC) ----
  # Delegates to .select_variants_by_spec() so per-gene (this function) and
  # chromosome-wide (compute_annotation_medians) filtering stay consistent.
  selected <- .select_variants_by_spec(
    gds = gds,
    candidate_variant_ids = region_variant_ids,
    filter_spec = filter_spec,
    sample_id = sample_id,
    Annotation_dir = Annotation_dir,
    Annotation_name_catalog = Annotation_name_catalog,
    impute_method = impute_method,
    return_genotypes = TRUE,
    reorder_samples = TRUE
  )

  n_passing <- length(selected$variant_ids)
  if (n_passing < filter_spec$min_variants) {
    if (verbose >= 2) {
      if (n_passing == 0) {
        message("Region ", region$label, ": 0 qualifying variants, skipping")
      } else {
        message("Region ", region$label, ": ", n_passing,
                " variants < min_variants (", filter_spec$min_variants,
                "), skipping")
      }
    }
    return(NULL)
  }

  final_variant_ids <- selected$variant_ids
  G <- selected$G
  MAF <- selected$MAF
  MAC <- selected$MAC

  # ---- 7. Build variant_info ----
  SeqArray::seqSetFilter(gds, variant.id = final_variant_ids, verbose = FALSE)
  positions <- SeqArray::seqGetData(gds, "position")
  alleles <- SeqArray::seqGetData(gds, "allele")
  allele_split <- strsplit(alleles, ",")
  ref <- vapply(allele_split, `[`, character(1), 1)
  alt <- vapply(allele_split, `[`, character(1), 2)

  # Pull rsid from annotation/id when present (typical dbSNP-annotated VCFs);
  # fall back to NA when the GDS has no such node.
  rsid <- tryCatch(
    SeqArray::seqGetData(gds, "annotation/id"),
    error = function(e) rep(NA_character_, length(final_variant_ids))
  )
  if (length(rsid) != length(final_variant_ids)) {
    rsid <- rep(NA_character_, length(final_variant_ids))
  }

  variant_info <- data.frame(
    variant_id = final_variant_ids,
    rsid = rsid,
    chr = rep(as.character(region$chr), length(final_variant_ids)),
    pos = positions,
    ref = ref,
    alt = alt,
    MAF = MAF,
    MAC = MAC,
    stringsAsFactors = FALSE
  )

  # ---- 8. Extract annotation scores (for PI) ----
  annotations <- NULL
  if (!is.null(annotation_names) && length(annotation_names) > 0) {
    annotations <- .extract_annotation_scores(
      gds, annotation_names, Annotation_dir, Annotation_name_catalog
    )
  }

  # ---- 9. Return glow_variant_set ----
  structure(
    list(
      G = G,
      variant_info = variant_info,
      annotations = annotations,
      region = region,
      filter_spec = filter_spec,
      n_samples = nrow(G),
      n_variants = ncol(G),
      n_total_in_region = n_total_in_region,
      n_after_annotation = selected$n_after_annotation
    ),
    class = "glow_variant_set"
  )
}


#' Compute Chromosome-Wide Annotation Medians
#'
#' @description
#' Computes the median of each annotation score across all variants in a GDS
#' file that pass the supplied \code{filter_spec}. These medians serve as
#' reference values for imputing NA annotations when computing PI scores
#' per gene, so the filter applied here must match the per-gene filter used
#' by \code{extract_variant_set()} — otherwise the global and per-gene
#' imputation pools are drawn from different populations.
#'
#' @param gds An open \code{SeqVarGDSClass} object or a file path to a GDS file.
#'   If a file path, the file is opened and closed automatically.
#' @param annotation_names Character vector of annotation column names to
#'   compute medians for (e.g., from \code{.default_PI_features()}).
#' @param filter_spec A \code{glow_variant_filter} object. Applied in full:
#'   QC status, variant type, annotation clauses (pLoF, missense, etc.),
#'   MAF cutoff, and MAC cutoff. The \code{min_variants} field is ignored —
#'   this is a reference-pool computation, not a per-gene test, so a small
#'   pool is acceptable.
#' @param sample_id Character vector of sample IDs (optional). Used for
#'   MAF/MAC computation. If NULL, all samples in the GDS are used.
#' @param Annotation_dir Character. GDS path prefix for annotations.
#'   Default: \code{"annotation/info/FunctionalAnnotation"}.
#' @param Annotation_name_catalog Data frame mapping annotation names to GDS
#'   paths. If NULL, uses the package-bundled catalog.
#' @param verbose Integer. Verbosity level (0 = silent, 1 = progress).
#'
#' @return Named numeric vector of per-annotation medians. Names correspond
#'   to \code{annotation_names}. An entry is \code{NA} if no non-NA values
#'   exist for that annotation across all matching variants.
#'
#' @details
#' Filter consistency with per-gene extraction is the point of this function.
#' When \code{prepare_glow_input()} falls back to per-gene median imputation
#' for an annotation that is non-NA in at least some of the gene's variants,
#' those variants were selected by \code{extract_variant_set()} under the
#' same \code{filter_spec}. For consistency, the global reference passed via
#' \code{reference_medians} must be drawn from the same filter pool —
#' chromosome-wide instead of gene-scoped. Both functions delegate filtering
#' to the internal \code{.select_variants_by_spec()} helper to guarantee
#' this pooling rule by construction.
#'
#' MAF and MAC are computed from \code{SeqArray::seqAlleleFreq()}, which
#' avoids materializing a whole-chromosome dosage matrix and lets this
#' function scale to real WGS cohorts.
#'
#' @examples
#' \dontrun{
#' spec <- coding_filter("plof")
#' medians <- compute_annotation_medians(
#'   "chr22.gds",
#'   annotation_names = c("cadd_phred", "linsight"),
#'   filter_spec = spec,
#'   sample_id = my_sample_ids
#' )
#' }
#'
#' @export
compute_annotation_medians <- function(gds,
                                        annotation_names,
                                        filter_spec,
                                        sample_id = NULL,
                                        Annotation_dir = "annotation/info/FunctionalAnnotation",
                                        Annotation_name_catalog = NULL,
                                        verbose = 1) {

  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required. Install with: BiocManager::install('SeqArray')")
  }

  if (missing(filter_spec) || is.null(filter_spec)) {
    stop("filter_spec is required so the global reference pool matches the ",
         "per-gene extraction filter")
  }

  # ---- 1. Handle GDS connection ----
  if (is.character(gds)) {
    stopifnot(file.exists(gds))
    gds <- SeqArray::seqOpen(gds, readonly = TRUE)
    on.exit(SeqArray::seqClose(gds), add = TRUE)
  }

  # ---- 2. Resolve annotation name catalog ----
  if (is.null(Annotation_name_catalog)) {
    Annotation_name_catalog <- get("annotation_name_catalog",
                                    envir = asNamespace("GLOWr"))
  }

  # ---- 3. Apply full filter spec via shared helper ----
  SeqArray::seqResetFilter(gds, verbose = FALSE)
  all_variant_ids <- SeqArray::seqGetData(gds, "variant.id")

  selected <- .select_variants_by_spec(
    gds = gds,
    candidate_variant_ids = all_variant_ids,
    filter_spec = filter_spec,
    sample_id = sample_id,
    Annotation_dir = Annotation_dir,
    Annotation_name_catalog = Annotation_name_catalog,
    return_genotypes = FALSE,
    reorder_samples = FALSE
  )

  if (length(selected$variant_ids) == 0) {
    if (verbose >= 1) message("No variants match filter criteria")
    result <- rep(NA_real_, length(annotation_names))
    names(result) <- annotation_names
    return(result)
  }

  if (verbose >= 1) {
    message("Computing annotation medians from ",
            length(selected$variant_ids), " variants")
  }

  # ---- 4. Compute per-annotation medians ----
  # Helper already set the filter to selected$variant_ids (with sample_id
  # if supplied). Read one annotation at a time to cap memory.
  medians <- numeric(length(annotation_names))
  names(medians) <- annotation_names

  for (i in seq_along(annotation_names)) {
    gds_path <- .resolve_annotation_path(
      annotation_names[i], Annotation_dir, Annotation_name_catalog
    )
    vals <- as.numeric(SeqArray::seqGetData(gds, gds_path))
    medians[i] <- median(vals, na.rm = TRUE)
  }

  medians
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Select Variants Matching a Filter Specification
#'
#' Applies a \code{glow_variant_filter} to a pre-specified pool of variant
#' IDs in an open GDS connection: QC status, variant type, DNF annotation
#' clauses, then MAF and MAC cutoffs. Positional scoping is the caller's
#' responsibility — this helper only filters an already-chosen candidate set.
#'
#' MAF and MAC are computed from \code{SeqArray::seqAlleleFreq()} rather than
#' by materializing a dosage matrix, so the helper is memory-efficient for
#' large candidate pools (whole-chromosome calls). The MAC formula
#' \code{round(2 * n_total * AF)} reproduces the integer MAC that
#' \code{.impute_geno_vset()} would compute from a mean-imputed dosage matrix.
#'
#' \code{min_variants} is intentionally NOT enforced here: that is a caller
#' concern (per-gene drop vs. chromosome-wide reference pool).
#'
#' @param gds An open \code{SeqVarGDSClass}.
#' @param candidate_variant_ids Integer or character vector of variant IDs
#'   to consider. Must already reflect any positional scoping.
#' @param filter_spec A \code{glow_variant_filter}.
#' @param sample_id Optional sample ID vector. If supplied, the GDS sample
#'   filter is set to this vector before AF/MAC computation.
#' @param Annotation_dir Character GDS path prefix for annotations.
#' @param Annotation_name_catalog Data frame mapping annotation names to
#'   GDS sub-paths.
#' @param impute_method Character: "mean" or "minor". Only used when
#'   \code{return_genotypes = TRUE}.
#' @param return_genotypes Logical. If TRUE, reads \code{$dosage} for the
#'   final variant set and imputes missing genotypes, returning \code{G}.
#'   If FALSE, \code{G} is \code{NULL} and no dosage matrix is materialized.
#' @param reorder_samples Logical. Only meaningful when
#'   \code{return_genotypes = TRUE} and \code{sample_id} is supplied.
#'   If TRUE, reorders \code{G} rows to match \code{sample_id}.
#'
#' @return A list with:
#'   \item{variant_ids}{Final variant IDs passing all filters.}
#'   \item{MAF}{Numeric vector of minor allele frequencies.}
#'   \item{MAC}{Integer vector of minor allele counts.}
#'   \item{G}{Imputed dosage matrix if \code{return_genotypes = TRUE},
#'     otherwise \code{NULL}.}
#'
#' @keywords internal
#' @noRd
.select_variants_by_spec <- function(gds,
                                      candidate_variant_ids,
                                      filter_spec,
                                      sample_id = NULL,
                                      Annotation_dir,
                                      Annotation_name_catalog,
                                      impute_method = "mean",
                                      return_genotypes = FALSE,
                                      reorder_samples = FALSE) {

  empty_result <- list(
    variant_ids      = candidate_variant_ids[0],
    MAF              = numeric(0),
    MAC              = integer(0),
    G                = NULL,
    n_after_annotation = 0L
  )

  if (length(candidate_variant_ids) == 0) return(empty_result)

  # ---- 1. Set filter to candidate variants + sample subset ----
  if (!is.null(sample_id)) {
    SeqArray::seqSetFilter(gds, variant.id = candidate_variant_ids,
                           sample.id = sample_id, verbose = FALSE)
  } else {
    SeqArray::seqSetFilter(gds, variant.id = candidate_variant_ids,
                           verbose = FALSE)
  }

  n_candidates <- length(candidate_variant_ids)

  # ---- 2. QC filter ----
  # Falls back to all-PASS if the QC node is absent (plain GDS without
  # VCF FILTER annotations).
  qc_values <- tryCatch(
    SeqArray::seqGetData(gds, filter_spec$qc_label),
    error = function(e) rep(filter_spec$qc_pass_value, n_candidates)
  )
  qc_pass <- (qc_values == filter_spec$qc_pass_value)

  # ---- 3. Variant type filter ----
  if (filter_spec$variant_type != "variant") {
    allele_strs <- SeqArray::seqGetData(gds, "allele")
    is_snv <- .is_snv_from_alleles(allele_strs)
    type_pass <- if (filter_spec$variant_type == "SNV") is_snv else !is_snv
  } else {
    type_pass <- rep(TRUE, n_candidates)
  }

  # ---- 4. Annotation clause filter ----
  if (!is.null(filter_spec$annotation_clauses)) {
    anno_pass <- .evaluate_annotation_clauses(
      gds, filter_spec$annotation_clauses,
      Annotation_dir, Annotation_name_catalog
    )
  } else {
    anno_pass <- rep(TRUE, n_candidates)
  }

  # ---- 5. Combine narrow filters ----
  narrow_pass <- qc_pass & type_pass & anno_pass
  if (sum(narrow_pass) == 0) return(empty_result)

  narrow_ids <- candidate_variant_ids[narrow_pass]

  # ---- 6. Set filter to narrowed variants (keep sample filter) ----
  if (!is.null(sample_id)) {
    SeqArray::seqSetFilter(gds, variant.id = narrow_ids,
                           sample.id = sample_id, verbose = FALSE)
  } else {
    SeqArray::seqSetFilter(gds, variant.id = narrow_ids, verbose = FALSE)
  }

  # ---- 7. Compute AF via seqAlleleFreq (memory-efficient, no dosage read) ----
  # ref.allele = 1L returns alt allele frequency, matching
  # colMeans($dosage, na.rm=TRUE) / 2.
  af_alt <- SeqArray::seqAlleleFreq(gds, ref.allele = 1L, verbose = FALSE)
  MAF <- pmin(af_alt, 1 - af_alt)

  # ---- 8. Compute MAC from AF and sample count ----
  # round(2 * n_total * af_alt) equals colSums(G_imputed) where
  # G_imputed has NAs replaced by 2*af_alt. See .impute_geno_vset().
  n_total <- length(SeqArray::seqGetData(gds, "sample.id"))
  MAC_raw <- as.integer(round(2 * n_total * af_alt))
  MAC <- pmin(MAC_raw, 2L * as.integer(n_total) - MAC_raw)

  # ---- 9. Apply MAF and MAC filter ----
  maf_pass <- MAF <= filter_spec$rare_maf_cutoff
  mac_pass <- MAC >= filter_spec$min_mac
  keep <- maf_pass & mac_pass

  if (sum(keep) == 0) {
    # narrow_ids passed QC+type+anno but all failed MAF/MAC.
    empty_result$n_after_annotation <- length(narrow_ids)
    return(empty_result)
  }

  final_ids <- narrow_ids[keep]
  MAF <- MAF[keep]
  MAC <- MAC[keep]

  # ---- 10. Set GDS filter to final_ids so the caller sees a consistent
  #         state. Callers that only want medians (return_genotypes = FALSE)
  #         rely on this filter being narrowed to final_ids, not narrow_ids.
  if (!is.null(sample_id)) {
    SeqArray::seqSetFilter(gds, variant.id = final_ids,
                           sample.id = sample_id, verbose = FALSE)
  } else {
    SeqArray::seqSetFilter(gds, variant.id = final_ids, verbose = FALSE)
  }

  # ---- 11. Read genotypes only when requested ----
  G <- NULL
  if (return_genotypes) {
    G <- SeqArray::seqGetData(gds, "$dosage")
    if (!is.matrix(G)) G <- matrix(G, ncol = 1)

    if (reorder_samples && !is.null(sample_id)) {
      gds_sample_ids <- SeqArray::seqGetData(gds, "sample.id")
      reorder_idx <- match(sample_id, gds_sample_ids)
      if (any(is.na(reorder_idx))) {
        n_missing <- sum(is.na(reorder_idx))
        stop(n_missing, " sample IDs not found in GDS file")
      }
      G <- G[reorder_idx, , drop = FALSE]
    }

    # Impute missing dosages. The AF values that drive imputation match
    # af_alt[keep] up to floating-point identity with colMeans(G, na.rm=TRUE)/2,
    # so the returned MAF/MAC (from steps 7-8) remain correct.
    geno_info <- .impute_geno_vset(G, method = impute_method)
    G <- geno_info$G
  }

  list(
    variant_ids        = final_ids,
    MAF                = MAF,
    MAC                = MAC,
    G                  = G,
    n_after_annotation = length(narrow_ids)
  )
}


#' Normalize region input to a standard list
#'
#' Accepts a data.frame row or a list and returns a standardized list
#' with chr, start, end, region_id, label.
#'
#' @param region A list or single-row data.frame.
#' @return Standardized list.
#' @keywords internal
#' @noRd
.normalize_region <- function(region) {
  if (is.data.frame(region)) {
    stopifnot(nrow(region) == 1)
    region <- as.list(region[1, ])
  }

  required <- c("chr", "start", "end")
  missing_fields <- setdiff(required, names(region))
  if (length(missing_fields) > 0) {
    stop("Region missing required fields: ", paste(missing_fields, collapse = ", "))
  }

  region$chr <- as.character(region$chr)
  region$start <- as.integer(region$start)
  region$end <- as.integer(region$end)

  if (is.null(region$region_id)) {
    region$region_id <- paste0("chr", region$chr, ":", region$start, "-", region$end)
  }
  if (is.null(region$label)) {
    region$label <- region$region_id
  }

  region
}


#' Evaluate DNF annotation clauses against current GDS filter
#'
#' Reads annotation values from the GDS file for the currently filtered
#' variants and evaluates the DNF (disjunctive normal form) clauses.
#'
#' @param gds An open SeqVarGDSClass with filter set.
#' @param clauses List of clauses (OR-ed). Each clause is a named list
#'   of field = accepted_values (AND-ed).
#' @param Annotation_dir Character GDS path prefix.
#' @param Annotation_name_catalog Data frame with name and dir columns.
#' @return Logical vector (one per filtered variant).
#' @keywords internal
#' @noRd
.evaluate_annotation_clauses <- function(gds, clauses, Annotation_dir,
                                          Annotation_name_catalog) {
  n <- length(SeqArray::seqGetData(gds, "variant.id"))

  # Collect all required fields and read annotation values
  all_fields <- unique(unlist(lapply(clauses, names)))
  anno_values <- list()
  for (field in all_fields) {
    gds_path <- .resolve_annotation_path(field, Annotation_dir,
                                          Annotation_name_catalog)
    anno_values[[field]] <- SeqArray::seqGetData(gds, gds_path)
  }

  # Evaluate each clause (OR of ANDs)
  result <- rep(FALSE, n)
  for (clause in clauses) {
    clause_result <- rep(TRUE, n)
    for (field_name in names(clause)) {
      accepted_values <- clause[[field_name]]
      field_values <- anno_values[[field_name]]
      clause_result <- clause_result & (field_values %in% accepted_values)
    }
    result <- result | clause_result
  }

  result
}


#' Resolve annotation field name to full GDS path
#'
#' @param field_name Character annotation name (e.g., "GENCODE.Category").
#' @param Annotation_dir Character GDS path prefix.
#' @param catalog Data frame with name and dir columns.
#' @return Character full GDS node path.
#' @keywords internal
#' @noRd
.resolve_annotation_path <- function(field_name, Annotation_dir, catalog) {
  # Try catalog name first (e.g., "CADD" -> "cadd_phred")
  idx <- which(catalog$name == field_name)
  if (length(idx) > 0) {
    return(paste0(Annotation_dir, "/", catalog$dir[idx[1]]))
  }
  # Try catalog dir column (e.g., "cadd_phred" maps to itself)
  idx_dir <- which(catalog$dir == field_name)
  if (length(idx_dir) > 0) {
    return(paste0(Annotation_dir, "/", catalog$dir[idx_dir[1]]))
  }
  # Fall back to using field_name directly as GDS sub-node name
  paste0(Annotation_dir, "/", field_name)
}


#' Extract numeric annotation scores for filtered variants
#'
#' @param gds An open SeqVarGDSClass with filter set.
#' @param annotation_names Character vector of annotation names.
#' @param Annotation_dir Character GDS path prefix.
#' @param catalog Data frame with name and dir columns.
#' @return Numeric matrix (n_variants x n_annotations).
#' @keywords internal
#' @noRd
.extract_annotation_scores <- function(gds, annotation_names, Annotation_dir,
                                        catalog) {
  n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
  scores <- matrix(NA_real_, nrow = n_variants, ncol = length(annotation_names))
  colnames(scores) <- annotation_names

  for (i in seq_along(annotation_names)) {
    gds_path <- .resolve_annotation_path(annotation_names[i], Annotation_dir,
                                          catalog)
    vals <- SeqArray::seqGetData(gds, gds_path)
    scores[, i] <- as.numeric(vals)
  }
  scores
}


#' Impute Missing Genotypes and Compute Allele Frequencies
#'
#' Same logic as .impute_geno() in marginal_scan.R, adapted for variant
#' set extraction. Replaces NA with column mean (2*AF) or minor allele,
#' computes MAF and MAC.
#'
#' @param G Genotype dosage matrix (n x m), may contain NAs.
#' @param method Character: "mean" or "minor".
#' @return List with G (imputed), MAF, MAC.
#' @keywords internal
#' @noRd
.impute_geno_vset <- function(G, method = "mean") {
  m <- ncol(G)
  n <- nrow(G)

  AF <- colMeans(G, na.rm = TRUE) / 2
  MAF <- pmin(AF, 1 - AF)

  if (method == "mean") {
    for (j in seq_len(m)) {
      na_idx <- is.na(G[, j])
      if (any(na_idx)) G[na_idx, j] <- 2 * AF[j]
    }
  } else if (method == "minor") {
    for (j in seq_len(m)) {
      na_idx <- is.na(G[, j])
      if (any(na_idx)) G[na_idx, j] <- ifelse(AF[j] <= 0.5, 0, 2)
    }
  }

  MAC <- as.integer(round(colSums(G)))
  MAC <- pmin(MAC, 2L * as.integer(n) - MAC)

  list(G = G, MAF = MAF, MAC = MAC)
}


#' Determine SNV Status from Allele Strings
#'
#' A variant is classified as SNV when both the reference and alternate
#' alleles are single nucleotides. Replaces SeqVarTools::isSNV() to
#' avoid an extra dependency.
#'
#' @param allele_strs Character vector of allele strings (format "REF,ALT").
#' @return Logical vector.
#' @keywords internal
#' @noRd
.is_snv_from_alleles <- function(allele_strs) {
  allele_split <- strsplit(allele_strs, ",")
  ref <- vapply(allele_split, `[`, character(1), 1)
  alt <- vapply(allele_split, `[`, character(1), 2)
  nchar(ref) == 1L & nchar(alt) == 1L
}
