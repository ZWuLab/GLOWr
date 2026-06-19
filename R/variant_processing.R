# ==============================================================================
# Post-Extraction Variant Processing
# ==============================================================================
#
# Processing functions for variant sets after extraction from GDS. These
# functions operate on genotype matrices (G) or glow_variant_set objects
# to prepare data for downstream B/PI computation and GLOW testing.
#
# EXPORTED FUNCTIONS:
#   - flip_alleles()            Flip genotype dosage to minor allele coding
#   - filter_variants_ld()      LD pruning + linear dependence removal
#   - collapse_rare_variants()  Collapse ultra-rare variants by rowSums
#   - aggregate_B_PI()          Aggregate B/PI vectors after collapsing
#
# INTERNAL HELPERS:
#   - .greedy_ld_prune()        Greedy LD pruning with MAF-based tie-breaking
#   - .group_rare_variants()    Group rare indices into chunks
#   - .group_rare_spatial()     Group rare variants between adjacent common variants
#   - .build_collapsed_result() Build collapsed genotype matrix and mapping


#################### EXPORTED MAIN FUNCTIONS ####################

#' Flip Genotype Dosage to Minor Allele Coding
#'
#' For columns where the allele frequency exceeds 0.5 (i.e., the coded
#' allele is the major allele), flips dosage values: 0 becomes 2 and 2
#' becomes 0, with 1 unchanged. This ensures all variants are coded
#' relative to the minor allele.
#'
#' @param G Numeric matrix (n samples x p variants) of dosage values (0, 1, 2).
#'   Missing values (NA) should be imputed before calling this function.
#'
#' @return Numeric matrix of same dimensions with flipped columns where needed.
#'
#' @examples
#' # Variant 2 has AF > 0.5, so it gets flipped
#' G <- matrix(c(0, 0, 2, 2,   # variant 1: AF = 0.5 (no flip)
#'               2, 2, 0, 0),   # variant 2: AF = 0.5 (no flip)
#'             nrow = 4, ncol = 2)
#' flip_alleles(G)
#'
#' # Variant with AF = 0.75 gets flipped: 2→0, 0→2
#' G2 <- matrix(c(2, 2, 2, 0), ncol = 1)  # AF = 0.75
#' flip_alleles(G2)  # becomes c(0, 0, 0, 2)
#'
#' @export
flip_alleles <- function(G) {
  stopifnot(is.matrix(G), is.numeric(G))

  AF <- colMeans(G, na.rm = TRUE) / 2
  to_flip <- which(AF > 0.5)

  if (length(to_flip) > 0) {
    G[, to_flip] <- 2 - G[, to_flip]
  }

  G
}


#' Filter Variants by Linkage Disequilibrium and Linear Dependence
#'
#' Identifies pairs of variants with absolute correlation exceeding a
#' threshold and removes one from each pair using a greedy algorithm.
#' Optionally also removes linearly dependent columns via QR decomposition.
#'
#' @param G Numeric matrix (n x p) of genotype dosage.
#' @param ld_threshold Numeric correlation threshold. Pairs with
#'   |r| > ld_threshold are considered redundant. Default 0.9.
#' @param remove_lindep Logical. If TRUE (default), also remove columns
#'   that are exact linear combinations of other columns after LD pruning.
#'   Uses QR decomposition. Prevents singularity in downstream Z-score
#'   computation and correlation matrices.
#' @param prefer_keep Character controlling which variant to keep when
#'   breaking ties in LD pruning. One of:
#'   \itemize{
#'     \item \code{"lower_maf"} (default): keep the rarer variant
#'     \item \code{"higher_maf"}: keep the more common variant
#'     \item \code{"first"}: keep the earlier column (positional order)
#'   }
#'
#' @return Integer vector of column indices to KEEP.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' # Three variants: v1 and v2 highly correlated, v3 independent
#' v1 <- rbinom(n, 2, 0.1)
#' v2 <- v1  # perfect correlation
#' v2[1:3] <- 2 - v2[1:3]  # introduce small differences
#' v3 <- rbinom(n, 2, 0.05)
#' G <- cbind(v1, v2, v3)
#'
#' keep <- filter_variants_ld(G, ld_threshold = 0.9)
#' keep  # one of v1/v2 removed, v3 kept
#' G_pruned <- G[, keep, drop = FALSE]
#'
#' @export
filter_variants_ld <- function(G, ld_threshold = 0.9, remove_lindep = TRUE,
                                prefer_keep = "lower_maf") {
  stopifnot(is.matrix(G), is.numeric(G))
  stopifnot(ld_threshold > 0 && ld_threshold <= 1)
  prefer_keep <- match.arg(prefer_keep, c("lower_maf", "higher_maf", "first"))

  p <- ncol(G)
  if (p <= 1) return(seq_len(p))

  # Compute MAF for tie-breaking (needed for "lower_maf" and "higher_maf")
  MAF <- if (prefer_keep != "first") {
    af <- colMeans(G, na.rm = TRUE) / 2
    pmin(af, 1 - af)
  } else {
    NULL
  }

  # Handle zero-variance columns (monomorphic after imputation)
  col_vars <- apply(G, 2, var, na.rm = TRUE)
  zero_var <- col_vars < .Machine$double.eps

  if (all(zero_var)) return(integer(0))

  # Only compute correlation among variable columns
  var_idx <- which(!zero_var)
  if (length(var_idx) <= 1) {
    keep <- var_idx
  } else {
    R <- cor(G[, var_idx, drop = FALSE], use = "pairwise.complete.obs")

    # Greedy removal with MAF-based tie-breaking
    # (replaces caret::findCorrelation)
    var_MAF <- if (!is.null(MAF)) MAF[var_idx] else NULL
    to_remove <- .greedy_ld_prune(R, ld_threshold, var_MAF, prefer_keep)

    keep_among_var <- setdiff(seq_along(var_idx), to_remove)
    keep <- var_idx[keep_among_var]
  }

  # Linear dependence removal via QR decomposition
  # After LD pruning, some columns may still be exact linear combinations
  # of others (e.g., after ultra-rare collapsing). QR decomposition
  # identifies the maximal linearly independent subset.
  if (remove_lindep && length(keep) > 1) {
    G_sub <- G[, keep, drop = FALSE]
    qr_result <- qr(G_sub)
    lindep_keep <- qr_result$pivot[seq_len(qr_result$rank)]
    keep <- keep[sort(lindep_keep)]
  }

  keep
}


#' Collapse Ultra-Rare Variants
#'
#' Groups ultra-rare variants (MAC below threshold) and collapses each
#' group into a single "super-variant" by summing genotype dosages.
#' Common variants (MAC >= threshold) are left unchanged.
#'
#' @param G Numeric matrix (n x p) of genotype dosage.
#' @param mac_threshold Integer MAC threshold. Variants with MAC < this
#'   are candidates for collapsing. Default 10.
#' @param max_group_size Integer maximum variants per collapsed group.
#'   Inf means all rare variants collapse into one. Default Inf.
#' @param spatial_grouping Logical. If TRUE, groups rare variants only
#'   between pairs of adjacent common variants (respects genomic order).
#'   If FALSE, groups all rare variants sequentially. Default FALSE.
#' @param agg_method Character method for aggregating B/PI values of
#'   collapsed variants. One of "mean", "max", "sum". Default "mean".
#'
#' @return A list with components:
#'   \describe{
#'     \item{G_collapsed}{Numeric matrix (n x p') of collapsed genotypes}
#'     \item{col_mapping}{List of length p'. Each element is an integer
#'       vector of original column indices that were merged into this column.}
#'     \item{is_collapsed}{Logical vector of length p'. TRUE for output
#'       columns formed by merging 2 or more variants. A singleton rare
#'       variant that did not get merged with anything is treated as a
#'       passthrough column and reported as FALSE.}
#'     \item{agg_method}{Character aggregation method (for downstream B/PI computation).}
#'   }
#'
#' @examples
#' # 50 samples, 5 variants with varying MAC
#' set.seed(1)
#' G <- cbind(
#'   rbinom(50, 2, 0.3),   # common (MAC ~ 30)
#'   rbinom(50, 2, 0.02),  # rare   (MAC ~ 2)
#'   rbinom(50, 2, 0.01),  # rare   (MAC ~ 1)
#'   rbinom(50, 2, 0.25),  # common (MAC ~ 25)
#'   rbinom(50, 2, 0.03)   # rare   (MAC ~ 3)
#' )
#' result <- collapse_rare_variants(G, mac_threshold = 10)
#' ncol(result$G_collapsed)  # fewer columns: rare variants grouped
#' result$is_collapsed       # TRUE for collapsed groups
#'
#' @export
collapse_rare_variants <- function(G,
                                    mac_threshold = 10L,
                                    max_group_size = Inf,
                                    spatial_grouping = FALSE,
                                    agg_method = "mean") {
  stopifnot(is.matrix(G), is.numeric(G))
  stopifnot(agg_method %in% c("mean", "max", "sum"))

  p <- ncol(G)
  n <- nrow(G)
  mac_threshold <- as.integer(mac_threshold)

  # Compute MAC per column
  MAC <- pmin(colSums(G), 2L * n - colSums(G))
  is_rare <- MAC < mac_threshold
  is_common <- !is_rare

  # If no rare variants, return unchanged
  if (!any(is_rare)) {
    return(list(
      G_collapsed = G,
      col_mapping = as.list(seq_len(p)),
      is_collapsed = rep(FALSE, p),
      agg_method = agg_method
    ))
  }

  # If no common variants, all rare
  if (!any(is_common)) {
    groups <- .group_rare_variants(which(is_rare), max_group_size)
    return(.build_collapsed_result(G, groups, integer(0), agg_method))
  }

  # Group rare variants
  rare_idx <- which(is_rare)
  common_idx <- which(is_common)

  if (spatial_grouping) {
    groups <- .group_rare_spatial(rare_idx, common_idx, max_group_size)
  } else {
    groups <- .group_rare_variants(rare_idx, max_group_size)
  }

  .build_collapsed_result(G, groups, common_idx, agg_method)
}


#' Aggregate B and PI Vectors After Collapsing
#'
#' After collapsing ultra-rare variants, the per-variant B and PI values
#' need to be aggregated for collapsed groups. This function uses the
#' col_mapping from collapse_rare_variants() to aggregate B and PI.
#'
#' @param B Numeric vector of per-variant B values (length p, original).
#' @param PI Numeric vector of per-variant PI values (length p, original).
#' @param collapse_result Output from \code{collapse_rare_variants()}.
#'
#' @return List with B_collapsed and PI_collapsed (length p').
#'
#' @examples
#' # Suppose we have 5 original variants with B and PI values
#' B <- c(0.5, 0.8, 0.9, 0.3, 0.7)
#' PI <- c(0.1, 0.4, 0.6, 0.2, 0.5)
#'
#' # After collapsing, variants 2-3 and 5 were grouped (see collapse_rare_variants)
#' collapse_result <- list(
#'   col_mapping = list(1L, c(2L, 3L), 4L, 5L),
#'   is_collapsed = c(FALSE, TRUE, FALSE, FALSE),
#'   agg_method = "mean"
#' )
#' agg <- aggregate_B_PI(B, PI, collapse_result)
#' agg$B_collapsed   # c(0.5, mean(0.8,0.9), 0.3, 0.7)
#' agg$PI_collapsed  # c(0.1, mean(0.4,0.6), 0.2, 0.5)
#'
#' @export
aggregate_B_PI <- function(B, PI, collapse_result) {
  agg_fn <- switch(collapse_result$agg_method,
    "mean" = mean,
    "max" = max,
    "sum" = sum
  )

  B_new <- vapply(collapse_result$col_mapping, function(idx) agg_fn(B[idx]),
                  numeric(1))
  PI_new <- vapply(collapse_result$col_mapping, function(idx) agg_fn(PI[idx]),
                   numeric(1))

  list(B_collapsed = B_new, PI_collapsed = PI_new)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Greedy LD pruning with MAF-based tie-breaking
#'
#' Repeatedly removes the variant with the most high-LD partners until
#' no pair exceeds the threshold. When multiple variants tie for most
#' partners, prefer_keep determines which to remove:
#' - "lower_maf": remove the one with higher MAF (keep rarer)
#' - "higher_maf": remove the one with lower MAF (keep more common)
#' - "first": remove the one with higher index (keep earlier)
#'
#' @param R Correlation matrix (p x p).
#' @param threshold Numeric LD threshold.
#' @param MAF Numeric vector of MAFs (length p), or NULL.
#' @param prefer_keep Character tie-breaking strategy.
#' @return Integer vector of indices (within R) to REMOVE.
#' @keywords internal
#' @noRd
.greedy_ld_prune <- function(R, threshold, MAF = NULL, prefer_keep = "first") {
  p <- ncol(R)
  diag(R) <- 0  # ignore self-correlation
  removed <- logical(p)

  repeat {
    # Count high-LD pairs per variant (among non-removed)
    high_ld <- abs(R) > threshold
    high_ld[removed, ] <- FALSE
    high_ld[, removed] <- FALSE

    n_partners <- colSums(high_ld)
    if (max(n_partners) == 0) break

    # Find candidates with the most high-LD partners
    max_partners <- max(n_partners)
    candidates <- which(n_partners == max_partners)

    # Tie-breaking: decide which candidate to remove
    if (length(candidates) == 1) {
      worst <- candidates
    } else if (prefer_keep == "lower_maf" && !is.null(MAF)) {
      # Remove the candidate with the HIGHEST MAF (keep rarer)
      worst <- candidates[which.max(MAF[candidates])]
    } else if (prefer_keep == "higher_maf" && !is.null(MAF)) {
      # Remove the candidate with the LOWEST MAF (keep more common)
      worst <- candidates[which.min(MAF[candidates])]
    } else {
      # "first" or fallback: remove the last one (keep earlier)
      worst <- candidates[length(candidates)]
    }

    removed[worst] <- TRUE
  }

  which(removed)
}


#' Group rare variant indices into chunks of max_group_size
#'
#' @param rare_idx Integer vector of rare variant column indices.
#' @param max_group_size Integer or Inf.
#' @return List of integer vectors (groups of indices).
#' @keywords internal
#' @noRd
.group_rare_variants <- function(rare_idx, max_group_size) {
  if (length(rare_idx) == 0) return(list())
  if (is.infinite(max_group_size)) return(list(rare_idx))

  split(rare_idx, ceiling(seq_along(rare_idx) / max_group_size))
}


#' Group rare variants between adjacent common variants
#'
#' Rare variants are grouped into segments defined by the positions of
#' common (MAC >= threshold) variants. Each segment contains rare variants
#' between two adjacent common variants (or before the first / after the
#' last common variant). Segments are further split by max_group_size.
#'
#' @param rare_idx Integer vector of rare variant column indices.
#' @param common_idx Integer vector of common variant column indices.
#' @param max_group_size Integer or Inf.
#' @return List of integer vectors (groups of indices).
#' @keywords internal
#' @noRd
.group_rare_spatial <- function(rare_idx, common_idx, max_group_size) {
  # Sort both by position (column index = genomic order)
  rare_idx <- sort(rare_idx)
  common_idx <- sort(common_idx)

  groups <- list()

  # Define boundaries: before first common, between each pair, after last common
  boundaries <- c(0, common_idx, max(c(rare_idx, common_idx)) + 1)

  for (i in seq_len(length(boundaries) - 1)) {
    lo <- boundaries[i]
    hi <- boundaries[i + 1]
    in_segment <- rare_idx[rare_idx > lo & rare_idx < hi]
    if (length(in_segment) > 0) {
      # Further split by max_group_size
      segment_groups <- .group_rare_variants(in_segment, max_group_size)
      groups <- c(groups, segment_groups)
    }
  }

  groups
}


#' Build the collapsed genotype matrix and mapping
#'
#' Takes rare variant groups and common variant indices and builds the
#' final collapsed output in positional order.
#'
#' @param G Original genotype matrix.
#' @param rare_groups List of integer vectors (rare variant groups).
#' @param common_idx Integer vector of common variant column indices.
#' @param agg_method Character aggregation method.
#' @return List with G_collapsed, col_mapping, is_collapsed, agg_method.
#' @keywords internal
#' @noRd
.build_collapsed_result <- function(G, rare_groups, common_idx, agg_method) {
  # Build output columns in positional order
  # Each output column is either: a single common variant, or a collapsed group

  # Determine position key for each output column
  # Common variants: position = their column index
  # Collapsed groups: position = min column index in group
  items <- list()
  for (ci in common_idx) {
    items[[length(items) + 1]] <- list(type = "common", idx = ci, pos = ci)
  }
  for (grp in rare_groups) {
    items[[length(items) + 1]] <- list(type = "group", idx = grp, pos = min(grp))
  }

  # Sort by position
  positions <- vapply(items, function(x) x$pos, numeric(1))
  items <- items[order(positions)]

  # Build output
  G_cols <- vector("list", length(items))
  col_mapping <- vector("list", length(items))
  is_collapsed <- logical(length(items))

  for (i in seq_along(items)) {
    item <- items[[i]]
    if (item$type == "common") {
      G_cols[[i]] <- G[, item$idx]
      col_mapping[[i]] <- item$idx
      is_collapsed[i] <- FALSE
    } else {
      G_cols[[i]] <- rowSums(G[, item$idx, drop = FALSE])
      col_mapping[[i]] <- item$idx
      # is_collapsed flags TRUE mergers only — a singleton rare group is a
      # passthrough (no aggregation occurred). See File Log 2026-05-14.
      is_collapsed[i] <- length(item$idx) > 1L
    }
  }

  G_collapsed <- do.call(cbind, G_cols)

  list(
    G_collapsed = G_collapsed,
    col_mapping = col_mapping,
    is_collapsed = is_collapsed,
    agg_method = agg_method
  )
}
