# ==============================================================================
# Prepare Data for GLOW Testing
# ==============================================================================
#
# Transforms a raw glow_variant_set (from extract_variant_set()) into processed
# data ready for score statistics and GLOW variant-set tests. Handles allele
# flipping, B/PI computation, LD filtering, rare variant collapsing, and
# B/PI aggregation.
#
# EXPORTED FUNCTIONS:
#   - prepare_glow_input()  Transform variant set into GLOW test-ready data
#   - print.glow_input()    S3 print method for glow_input objects
#
# INTERNAL HELPERS:
#   - .validate_prepare_inputs()  Input validation for prepare_glow_input
#
# Note: .default_PI_features() is defined in get_PI_train.R and accessible
#   within the package namespace. This file does NOT duplicate it.


#################### EXPORTED MAIN FUNCTIONS ####################

#' Prepare Data for GLOW Testing
#'
#' Transforms a raw \code{glow_variant_set} (from
#' \code{\link{extract_variant_set}}) into processed data ready for score
#' statistics and GLOW variant-set tests. Handles allele flipping, B/PI
#' computation, LD filtering, rare variant collapsing, and B/PI aggregation.
#'
#' @param variant_set A \code{glow_variant_set} object from
#'   \code{\link{extract_variant_set}}.
#' @param B Numeric vector (length = n_variants) of pre-computed B values, or
#'   NULL. Exactly one of \code{B}, \code{B_func}, \code{B_model} must be
#'   provided.
#' @param B_func Function taking a numeric MAF vector and returning B values.
#'   For example:
#'   \code{function(maf) sqrt(-0.14307 * log(maf * (1 - maf)) - 0.19522)}.
#' @param B_model A \code{glow_B_model} object. Calls
#'   \code{\link{predict_B}} internally.
#' @param PI Numeric vector (length = n_variants) of pre-computed PI values, or
#'   NULL. Exactly one of \code{PI} and \code{PI_models} must be provided.
#' @param PI_models Either a list from \code{\link{load_PI_models}} or a
#'   character path to a model directory.
#' @param PI_features Character vector of annotation column names to use for PI
#'   prediction. Required when using \code{PI_models}. If NULL, uses the
#'   default 16 curated PI training features from
#'   \code{.default_PI_features()}.
#' @param reference_medians Named numeric vector of chromosome-wide annotation
#'   medians for NA imputation (optional). From
#'   \code{\link{compute_annotation_medians}}. When provided, NA annotation
#'   values are imputed with these medians instead of per-gene column medians.
#' @param ld_threshold Numeric (default 0.9). Correlation threshold for LD
#'   pruning.
#' @param ld_prefer_keep Character (default "lower_maf"). Tie-breaking for LD
#'   pruning: "lower_maf", "higher_maf", or "first".
#' @param remove_lindep Logical (default TRUE). Remove linearly dependent
#'   columns via QR decomposition.
#' @param mac_threshold Integer (default 10). MAC threshold for rare variant
#'   collapsing. Set to 0 to disable collapsing.
#' @param collapse_method Character (default "mean"). B/PI aggregation method
#'   for collapsed groups: "mean", "max", or "sum".
#' @param max_group_size Numeric (default Inf). Maximum variants per collapsed
#'   group.
#' @param verbose Integer (default 1). Verbosity level.
#'
#' @return A \code{glow_input} S3 object (see Details), or NULL if too few
#'   variants remain after filtering.
#'
#' @details
#' The processing pipeline applies the following steps in order:
#' \enumerate{
#'   \item Validate inputs
#'   \item Flip alleles to minor allele coding via \code{\link{flip_alleles}}
#'   \item Compute B values (from vector, function, or model)
#'   \item Compute PI values (from vector or ensemble models)
#'   \item LD filtering via \code{\link{filter_variants_ld}}
#'   \item Rare variant collapsing via \code{\link{collapse_rare_variants}}
#'   \item Aggregate B/PI for collapsed groups via \code{\link{aggregate_B_PI}}
#' }
#'
#' The returned \code{glow_input} object contains:
#' \describe{
#'   \item{G}{Numeric matrix (n x p_final) of processed genotypes}
#'   \item{B}{Numeric vector (length p_final) of B values}
#'   \item{PI}{Numeric vector (length p_final) of PI values}
#'   \item{is_collapsed}{Logical vector indicating collapsed groups}
#'   \item{col_mapping}{List mapping final columns to original column indices}
#'   \item{ld_keep_idx}{Integer vector (length \code{n_after_ld}); element
#'     \code{j} is the post-filter variant index (i.e., 1-based row index into
#'     the input variant set's \code{variant_info}) occupying post-LD column
#'     \code{j}. Chained with \code{col_mapping}, this traces every post-collapse
#'     test unit back to its post-filter raw variants.}
#'   \item{region}{Region information from the input variant set}
#'   \item{n_original}{Raw variant count in the region before any filtering
#'     (from \code{variant_set$n_total_in_region})}
#'   \item{n_after_filter}{Variants remaining after Component 3's filter spec
#'     (QC, variant type, annotation mask, MAF/MAC), i.e., the variants
#'     entering \code{prepare_glow_input}}
#'   \item{n_after_ld}{Variants remaining after LD filtering}
#'   \item{n_after_collapse}{Variants remaining after rare-variant collapsing
#'     (= final number of columns in \code{G})}
#'   \item{cMAC}{Cumulative MAC across all final columns}
#'   \item{processing_log}{Character vector of processing step descriptions}
#' }
#'
#' @examples
#' \dontrun{
#' # Extract variant set from GDS
#' vset <- extract_variant_set(gds, region, filter_spec)
#'
#' # Method 1: Pre-computed B and PI vectors
#' result <- prepare_glow_input(vset, B = B_vec, PI = PI_vec)
#'
#' # Method 2: B from function, PI from models
#' B_formula <- function(maf) sqrt(-0.14307 * log(maf * (1 - maf)) - 0.19522)
#' result <- prepare_glow_input(vset, B_func = B_formula,
#'                               PI_models = "path/to/models/")
#'
#' # Method 3: B from trained model
#' result <- prepare_glow_input(vset, B_model = b_model, PI = PI_vec)
#' }
#'
#' @seealso \code{\link{extract_variant_set}} for creating the input variant set,
#'   \code{\link{flip_alleles}}, \code{\link{filter_variants_ld}},
#'   \code{\link{collapse_rare_variants}}, \code{\link{aggregate_B_PI}}
#'
#' @export
prepare_glow_input <- function(variant_set,
                               B = NULL, B_func = NULL, B_model = NULL,
                               PI = NULL, PI_models = NULL,
                               PI_features = NULL,
                               reference_medians = NULL,
                               ld_threshold = 0.9,
                               ld_prefer_keep = "lower_maf",
                               remove_lindep = TRUE,
                               mac_threshold = 10L,
                               collapse_method = "mean",
                               max_group_size = Inf,
                               verbose = 1) {

  # ---- 1. Validate inputs ----
  .validate_prepare_inputs(variant_set, B, B_func, B_model, PI, PI_models)

  # ---- 2. Extract data from variant set ----
  G <- variant_set$G
  annotations <- variant_set$annotations
  # Four-stage variant count tracking. n_original is the raw positional count
  # from the region (pre-Component-3-filter); n_after_filter is the count
  # entering this function (post-filter-spec from extract_variant_set).
  n_original <- if (is.null(variant_set$n_total_in_region)) ncol(G) else variant_set$n_total_in_region
  n_after_filter <- ncol(G)
  log_entries <- character()

  # ---- 3. Flip alleles to minor allele coding ----
  G <- flip_alleles(G)

  # Recompute MAF from flipped genotypes (used for B computation)
  maf <- colMeans(G, na.rm = TRUE) / 2

  log_entries <- c(log_entries,
                   paste0("Flipped alleles: ", n_after_filter, " variants"))

  # ---- 4. Compute B ----
  if (!is.null(B)) {
    # Pre-computed B vector
    if (length(B) != n_after_filter) {
      stop("B vector length (", length(B), ") must equal number of variants (",
           n_after_filter, ")")
    }
    B_vec <- B
  } else if (!is.null(B_func)) {
    # Compute B from user-supplied function
    B_vec <- B_func(maf)
  } else {
    # Compute B from trained model
    B_vec <- predict_B(B_model, target_MAF = maf)
  }

  # Validate B_vec output
  if (!is.numeric(B_vec) || length(B_vec) != n_after_filter) {
    stop("B computation must produce a numeric vector of length ", n_after_filter,
         ", got length ", length(B_vec))
  }

  # Handle edge cases: NAs and zero-MAF variants
  if (any(is.na(B_vec))) {
    warning("B contains ", sum(is.na(B_vec)), " NA values; replacing with 0")
    B_vec[is.na(B_vec)] <- 0
  }
  B_vec[maf == 0] <- 0

  log_entries <- c(log_entries,
                   paste0("Computed B: range [", round(min(B_vec), 4),
                          ", ", round(max(B_vec), 4), "]"))

  # ---- 5. Compute PI ----
  if (!is.null(PI)) {
    # Pre-computed PI vector
    if (length(PI) != n_after_filter) {
      stop("PI vector length (", length(PI), ") must equal number of variants (",
           n_after_filter, ")")
    }
    PI_vec <- PI
  } else {
    # Compute PI from ensemble models
    if (is.character(PI_models)) {
      PI_models <- load_PI_models(PI_models)
    }
    if (is.null(PI_features)) {
      PI_features <- .default_PI_features()
    }
    if (is.null(annotations)) {
      stop("PI_models requires annotations in variant_set. ",
           "Set annotation_names in extract_variant_set().")
    }
    # Check that requested PI features exist in annotations
    missing_feats <- setdiff(PI_features, colnames(annotations))
    if (length(missing_feats) > 0) {
      stop("PI features not found in annotations: ",
           paste(missing_feats, collapse = ", "))
    }
    annot_matrix <- annotations[, PI_features, drop = FALSE]
    pi_result <- predict_PI_ensemble(PI_models, annot_matrix,
                                      reference_medians = reference_medians)
    PI_vec <- pi_result$ensemble_pi
  }

  # Validate PI_vec output
  if (!is.numeric(PI_vec) || length(PI_vec) != n_after_filter) {
    stop("PI computation must produce a numeric vector of length ", n_after_filter,
         ", got length ", length(PI_vec))
  }

  log_entries <- c(log_entries,
                   paste0("Computed PI: mean=", round(mean(PI_vec), 4)))

  # ---- 6. LD filter ----
  ld_keep_idx <- filter_variants_ld(G, ld_threshold = ld_threshold,
                                     remove_lindep = remove_lindep,
                                     prefer_keep = ld_prefer_keep)

  # Subset all data to kept variants
  G <- G[, ld_keep_idx, drop = FALSE]
  B_vec <- B_vec[ld_keep_idx]
  PI_vec <- PI_vec[ld_keep_idx]
  if (!is.null(annotations)) {
    annotations <- annotations[ld_keep_idx, , drop = FALSE]
  }

  n_after_ld <- ncol(G)
  log_entries <- c(log_entries,
                   paste0("LD filter: ", n_after_filter, " -> ", n_after_ld,
                          " (removed ", n_after_filter - n_after_ld, ")"))

  # Check: enough variants remaining?
  if (n_after_ld < 1) {
    if (verbose >= 1) message("No variants after LD filtering")
    return(NULL)
  }

  # ---- 7. Collapse rare variants ----
  if (mac_threshold > 0) {
    collapse_result <- collapse_rare_variants(G,
                                              mac_threshold = mac_threshold,
                                              max_group_size = max_group_size,
                                              agg_method = collapse_method)
    G <- collapse_result$G_collapsed
    agg <- aggregate_B_PI(B_vec, PI_vec, collapse_result)
    B_vec <- agg$B_collapsed
    PI_vec <- agg$PI_collapsed
    is_collapsed <- collapse_result$is_collapsed
    col_mapping <- collapse_result$col_mapping
    n_collapsed <- sum(is_collapsed)
    log_entries <- c(log_entries,
                     paste0("Collapse: ", n_after_ld, " -> ", ncol(G),
                            " (", n_collapsed, " collapsed groups)"))
  } else {
    is_collapsed <- rep(FALSE, ncol(G))
    col_mapping <- as.list(seq_len(ncol(G)))
    log_entries <- c(log_entries, "Collapse: disabled")
  }

  # ---- 8. Final check ----
  n_after_collapse <- ncol(G)
  if (n_after_collapse < 1) {
    if (verbose >= 1) message("No variants after collapsing")
    return(NULL)
  }

  # ---- 9. Compute cumulative MAC ----
  cMAC <- sum(colSums(G))

  # ---- 10. Assemble and return glow_input object ----
  structure(
    list(
      G = G,
      B = B_vec,
      PI = PI_vec,
      is_collapsed = is_collapsed,
      col_mapping = col_mapping,
      ld_keep_idx = ld_keep_idx,
      region = variant_set$region,
      n_original = n_original,
      n_after_filter = n_after_filter,
      n_after_ld = n_after_ld,
      n_after_collapse = n_after_collapse,
      cMAC = cMAC,
      processing_log = log_entries
    ),
    class = "glow_input"
  )
}


#' Print Method for glow_input Objects
#'
#' @param x A \code{glow_input} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisible \code{x}.
#'
#' @export
print.glow_input <- function(x, ...) {
  cat("GLOW Input (processed variant set)\n")
  region_label <- if (!is.null(x$region$label)) x$region$label else "unnamed"
  cat("  Region: ", region_label, "\n", sep = "")
  cat("  Variants: ", x$n_original, " original -> ",
      x$n_after_filter, " after filter -> ",
      x$n_after_ld, " after LD -> ",
      x$n_after_collapse, " after collapse\n", sep = "")
  cat("  B range: [", round(min(x$B), 4), ", ",
      round(max(x$B), 4), "]\n", sep = "")
  cat("  PI range: [", round(min(x$PI), 4), ", ",
      round(max(x$PI), 4), "]\n", sep = "")
  cat("  cMAC: ", x$cMAC, "\n", sep = "")
  invisible(x)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Validate inputs for prepare_glow_input
#'
#' Checks class of variant_set and exactly-one-of constraints for B and PI
#' source parameters.
#'
#' @param variant_set The variant_set argument.
#' @param B The B argument.
#' @param B_func The B_func argument.
#' @param B_model The B_model argument.
#' @param PI The PI argument.
#' @param PI_models The PI_models argument.
#'
#' @return NULL (invisible). Stops with informative error on failure.
#' @keywords internal
#' @noRd
.validate_prepare_inputs <- function(variant_set, B, B_func, B_model,
                                      PI, PI_models) {
  # Validate variant_set class

  if (!inherits(variant_set, "glow_variant_set")) {
    stop("variant_set must be a 'glow_variant_set' object (from extract_variant_set())")
  }

  # Exactly one of B, B_func, B_model must be non-NULL
  n_B_sources <- sum(!is.null(B), !is.null(B_func), !is.null(B_model))
  if (n_B_sources != 1) {
    stop("Exactly one of B, B_func, B_model must be provided (got ",
         n_B_sources, ")")
  }

  # Exactly one of PI, PI_models must be non-NULL
  n_PI_sources <- sum(!is.null(PI), !is.null(PI_models))
  if (n_PI_sources != 1) {
    stop("Exactly one of PI, PI_models must be provided (got ",
         n_PI_sources, ")")
  }

  invisible(NULL)
}
