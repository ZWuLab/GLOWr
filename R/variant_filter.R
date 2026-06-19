# ==============================================================================
# Variant Filter Specification Functions
# ==============================================================================
#
# Functions for creating declarative variant filter specifications for use
# with extract_variant_set(). Filters encode QC, variant type, annotation,
# and MAF/MAC criteria without accessing data.
#
# EXPORTED FUNCTIONS:
#   - variant_filter()              General-purpose filter constructor
#   - coding_filter()               Predefined coding category filter
#   - print.glow_variant_filter()   Print method for filter specs
#
# INTERNAL HELPERS:
#   - .coding_annotation_masks()    Predefined coding category masks

#################### EXPORTED MAIN FUNCTIONS ####################

#' Create a Variant Filter Specification
#'
#' Builds a declarative filter specification for use with
#' \code{\link{extract_variant_set}}. The spec encodes filtering rules
#' (QC, variant type, annotation masks, MAF) without accessing data.
#'
#' @param qc_label Character GDS node path for QC status.
#'   Default "annotation/filter".
#' @param qc_pass_value Character value that passes QC. Default "PASS".
#' @param variant_type Character: "SNV", "Indel", or "variant" (both).
#'   Default "SNV".
#' @param annotation_clauses List of clauses in DNF form. Each clause is a
#'   named list of field = accepted_values conditions (AND-ed within clause).
#'   Clauses are OR-ed together. NULL means no annotation filtering.
#' @param rare_maf_cutoff Numeric MAF ceiling. Variants with MAF > this
#'   are excluded. Default 0.01.
#' @param min_mac Integer minimum cohort minor allele count. Variants with
#'   cohort-MAC < \code{min_mac} are excluded. Default 1 — variants that
#'   are monomorphic in the analysis cohort are excluded because they
#'   carry no information for any association test. Set to 0 only if you
#'   have a specific reason to retain MAC=0 variants (rare).
#' @param min_variants Integer minimum variants after filtering to proceed.
#'   Default 2.
#'
#' @return A \code{glow_variant_filter} S3 object.
#'
#' @examples
#' # Default filter: rare SNVs passing QC
#' spec <- variant_filter()
#' spec
#'
#' # Custom annotation filter in DNF form (plof_ds equivalent):
#' # stopgain OR stoploss OR splicing OR (nonsynonymous AND MetaSVM=="D")
#' spec_custom <- variant_filter(
#'   annotation_clauses = list(
#'     list("GENCODE.EXONIC.Category" = "stopgain"),
#'     list("GENCODE.EXONIC.Category" = "stoploss"),
#'     list("GENCODE.Category" = "splicing"),
#'     list("GENCODE.EXONIC.Category" = "nonsynonymous SNV", "MetaSVM" = "D")
#'   ),
#'   rare_maf_cutoff = 0.01
#' )
#'
#' # All rare variants (SNV + Indel), no annotation filter
#' spec_all <- variant_filter(variant_type = "variant", rare_maf_cutoff = 0.05)
#'
#' @export
variant_filter <- function(qc_label = "annotation/filter",
                            qc_pass_value = "PASS",
                            variant_type = "SNV",
                            annotation_clauses = NULL,
                            rare_maf_cutoff = 0.01,
                            min_mac = 1L,
                            min_variants = 2L) {
  stopifnot(variant_type %in% c("SNV", "Indel", "variant"))
  stopifnot(rare_maf_cutoff > 0 && rare_maf_cutoff <= 0.5)

  # Validate annotation_clauses structure
  if (!is.null(annotation_clauses)) {
    stopifnot(is.list(annotation_clauses))
    for (i in seq_along(annotation_clauses)) {
      clause <- annotation_clauses[[i]]
      stopifnot(is.list(clause), length(clause) > 0)
      stopifnot(all(nchar(names(clause)) > 0))
    }
  }

  # Collect required annotation fields from clauses
  required_fields <- if (!is.null(annotation_clauses)) {
    unique(unlist(lapply(annotation_clauses, names)))
  } else {
    character(0)
  }

  structure(
    list(
      qc_label = qc_label,
      qc_pass_value = qc_pass_value,
      variant_type = variant_type,
      annotation_clauses = annotation_clauses,
      required_annotation_fields = required_fields,
      rare_maf_cutoff = rare_maf_cutoff,
      min_mac = as.integer(min_mac),
      min_variants = as.integer(min_variants)
    ),
    class = "glow_variant_filter"
  )
}


#' Create a Predefined Coding Category Filter
#'
#' Convenience constructor that returns a \code{glow_variant_filter} with
#' the standard STAARpipeline boolean mask for a coding category.
#'
#' @param category Character: one of "plof", "plof_ds", "missense",
#'   "disruptive_missense", "synonymous", "ptv", "ptv_ds".
#' @param rare_maf_cutoff Numeric MAF ceiling. Default 0.01.
#' @param min_mac Integer minimum cohort MAC. Default 1
#'   (exclude variants monomorphic in the analysis cohort).
#' @param min_variants Integer minimum variants. Default 2.
#'
#' @return A \code{glow_variant_filter} S3 object with a
#'   \code{category_name} attribute.
#'
#' @examples
#' # Predefined coding category filters
#' spec_plof <- coding_filter("plof")
#' spec_plof
#'
#' spec_missense <- coding_filter("missense", rare_maf_cutoff = 0.05)
#'
#' # Available categories:
#' # "plof", "plof_ds", "missense", "disruptive_missense",
#' # "synonymous", "ptv", "ptv_ds"
#'
#' @export
coding_filter <- function(category,
                           rare_maf_cutoff = 0.01,
                           min_mac = 1L,
                           min_variants = 2L) {
  masks <- .coding_annotation_masks()
  if (!category %in% names(masks)) {
    stop("Unknown coding category: '", category, "'. ",
         "Available: ", paste(names(masks), collapse = ", "))
  }

  spec <- variant_filter(
    variant_type = "SNV",
    annotation_clauses = masks[[category]],
    rare_maf_cutoff = rare_maf_cutoff,
    min_mac = min_mac,
    min_variants = min_variants
  )
  attr(spec, "category_name") <- category
  spec
}


#' Print Method for Variant Filter Specifications
#'
#' @param x A \code{glow_variant_filter} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.glow_variant_filter <- function(x, ...) {
  cat("GLOWr Variant Filter Specification\n")
  cat("  Variant type:", x$variant_type, "\n")
  cat("  MAF cutoff:", x$rare_maf_cutoff, "\n")
  cat("  Min MAC:", x$min_mac, "\n")
  cat("  Min variants:", x$min_variants, "\n")
  if (!is.null(attr(x, "category_name"))) {
    cat("  Category:", attr(x, "category_name"), "\n")
  }
  if (!is.null(x$annotation_clauses)) {
    cat("  Annotation clauses:", length(x$annotation_clauses), "clause(s)\n")
    cat("  Required fields:", paste(x$required_annotation_fields, collapse = ", "), "\n")
  } else {
    cat("  Annotation filter: none (all passing variants)\n")
  }
  invisible(x)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Predefined Coding Category Annotation Masks
#'
#' Returns the 7 standard STAARpipeline coding category masks in DNF form.
#' Each mask is a list of clauses (OR-ed); each clause is a named list of
#' field = accepted_value conditions (AND-ed within clause).
#'
#' @return Named list of 7 coding category masks.
#' @keywords internal
#' @noRd
.coding_annotation_masks <- function() {
  # pLoF: predicted loss-of-function (includes ncRNA splicing)
  plof_clauses <- list(
    list("GENCODE.EXONIC.Category" = "stopgain"),
    list("GENCODE.EXONIC.Category" = "stoploss"),
    list("GENCODE.Category" = "splicing"),
    list("GENCODE.Category" = "exonic;splicing"),
    list("GENCODE.Category" = "ncRNA_splicing"),
    list("GENCODE.Category" = "ncRNA_exonic;splicing")
  )

  # Disruptive missense: nonsynonymous + MetaSVM predicted damaging
  disruptive_missense_clauses <- list(
    list("GENCODE.EXONIC.Category" = "nonsynonymous SNV",
         "MetaSVM" = "D")
  )

  # Missense: all nonsynonymous SNVs
  missense_clauses <- list(
    list("GENCODE.EXONIC.Category" = "nonsynonymous SNV")
  )

  # Synonymous
  synonymous_clauses <- list(
    list("GENCODE.EXONIC.Category" = "synonymous SNV")
  )

  # PTV: protein-truncating variants (no ncRNA categories)
  ptv_clauses <- list(
    list("GENCODE.EXONIC.Category" = "stopgain"),
    list("GENCODE.EXONIC.Category" = "stoploss"),
    list("GENCODE.Category" = "splicing"),
    list("GENCODE.Category" = "exonic;splicing")
  )

  list(
    plof = plof_clauses,
    plof_ds = c(plof_clauses, disruptive_missense_clauses),
    missense = missense_clauses,
    disruptive_missense = disruptive_missense_clauses,
    synonymous = synonymous_clauses,
    ptv = ptv_clauses,
    ptv_ds = c(ptv_clauses, disruptive_missense_clauses)
  )
}
