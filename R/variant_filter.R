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
#'   named list of \code{field = condition} terms (AND-ed within a clause);
#'   clauses are OR-ed together. A condition is either an atomic vector of
#'   accepted values (set membership, e.g. \code{"GENCODE.Category" =
#'   c("UTR3", "UTR5")}) or a predicate built with
#'   \code{\link{annotation_predicate}} for tests membership cannot express:
#'   non-empty / empty fields (\code{annotation_predicate("nonempty")}) and
#'   numeric comparisons (\code{annotation_predicate("gt", 20)}). NULL means
#'   no annotation filtering.
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
#' # Predicates beyond set membership: nonsynonymous SNVs with CADD PHRED > 20,
#' # OR any variant inside a GeneHancer element (non-empty field)
#' spec_pred <- variant_filter(
#'   annotation_clauses = list(
#'     list("GENCODE.EXONIC.Category" = "nonsynonymous SNV",
#'          "CADD" = annotation_predicate("gt", 20)),
#'     list("GeneHancer" = annotation_predicate("nonempty"))
#'   ),
#'   rare_maf_cutoff = 0.5
#' )
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
      for (field in names(clause)) {
        .validate_clause_condition(clause[[field]], field, i)
      }
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


#' Build an Annotation Predicate for Variant Filtering
#'
#' Expresses a per-variant condition on an annotation field that set
#' membership cannot: non-empty / empty fields and numeric comparisons. Use
#' the result as a condition value inside the \code{annotation_clauses} of
#' \code{\link{variant_filter}}.
#'
#' @param op Character operator: \code{"nonempty"} (not NA and not the empty
#'   string), \code{"empty"} (NA or the empty string), \code{"in"} /
#'   \code{"not_in"} (set membership), or a numeric comparison \code{"gt"},
#'   \code{"ge"}, \code{"lt"}, \code{"le"} (greater / less than, or equal).
#' @param value Comparison value: a single non-missing number for the numeric
#'   operators; a non-empty atomic vector for \code{"in"} / \code{"not_in"};
#'   must be NULL for \code{"nonempty"} / \code{"empty"}.
#'
#' @details
#' Missing values never satisfy a positive test: \code{NA} is \code{FALSE}
#' under \code{"gt"}, \code{"ge"}, \code{"lt"}, \code{"le"} and \code{"in"},
#' and \code{TRUE} under \code{"empty"} and \code{"not_in"}. For the numeric
#' operators, character-stored numbers are coerced with \code{as.numeric()};
#' values that do not parse count as missing.
#'
#' @return A \code{glow_annotation_predicate} S3 object.
#'
#' @examples
#' annotation_predicate("nonempty")
#' annotation_predicate("gt", 20)
#' annotation_predicate("in", c("D", "T"))
#'
#' @seealso \code{\link{variant_filter}}
#' @export
annotation_predicate <- function(op, value = NULL) {
  ops <- c("in", "not_in", "nonempty", "empty", "gt", "ge", "lt", "le")
  if (!is.character(op) || length(op) != 1L || !op %in% ops) {
    stop("annotation_predicate(): op must be one of ", paste(ops, collapse = ", "))
  }
  if (op %in% c("nonempty", "empty")) {
    if (!is.null(value)) stop("annotation_predicate(): op '", op, "' takes no value")
  } else if (op %in% c("gt", "ge", "lt", "le")) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
      stop("annotation_predicate(): op '", op, "' needs a single non-missing numeric value")
    }
  } else {
    if (is.null(value) || !is.atomic(value) || length(value) == 0L) {
      stop("annotation_predicate(): op '", op, "' needs a non-empty atomic vector of values")
    }
  }
  structure(list(op = op, value = value), class = "glow_annotation_predicate")
}


#' Print Method for Annotation Predicates
#'
#' @param x A \code{glow_annotation_predicate} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.glow_annotation_predicate <- function(x, ...) {
  val <- if (is.null(x$value)) "" else paste(format(x$value), collapse = ", ")
  cat("GLOWr annotation predicate:", x$op, val, "\n")
  invisible(x)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Is this object an annotation predicate?
#' @keywords internal
#' @noRd
.is_annotation_predicate <- function(x) inherits(x, "glow_annotation_predicate")


#' Validate one clause condition (atomic accepted-values vector or predicate)
#' @keywords internal
#' @noRd
.validate_clause_condition <- function(cond, field, clause_index) {
  if (.is_annotation_predicate(cond)) return(invisible(TRUE))
  if (!is.null(cond) && is.atomic(cond) && length(cond) > 0L) return(invisible(TRUE))
  stop(sprintf(paste0("variant_filter(): clause %d, field '%s': a condition must be ",
                      "an atomic vector of accepted values or an annotation_predicate()"),
               clause_index, field))
}

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
