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
#   - .independent_columns()   In-order linear-dependence check (replaces qr())
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
#' Optionally also removes linearly dependent columns, in column order.
#'
#' @param G Numeric matrix (n x p) of genotype dosage.
#' @param ld_threshold Numeric correlation threshold. Pairs with
#'   |r| > ld_threshold are considered redundant. Default 0.9.
#' @param remove_lindep Logical. If TRUE (default), also remove columns
#'   that are exact linear combinations of other columns after LD pruning.
#'   The columns are visited in order and a column is dropped when its
#'   residual after projection onto the kept columns before it has norm
#'   below \code{1e-07} times its own norm, which is the rule of \code{qr()}
#'   with its default tolerance. Prevents singularity in downstream Z-score
#'   computation and correlation matrices. Requires a matrix without
#'   missing values.
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
#' @details
#' Cost: for \eqn{N} samples and \eqn{m} variable columns, the correlation is
#' one BLAS cross-product, \eqn{O(N m^2)}, on the centered and scaled columns
#' (the two-pass column variance, equal to \code{var()}, decides the
#' zero-variance exclusion and scales the columns). The greedy prune thresholds
#' the correlation once and keeps each variant's count of above-threshold
#' partners up to date as variants are removed, \eqn{O(m^2 + k m)} for
#' \eqn{k} removals. Memory holds the \eqn{m \times m} correlation and its
#' logical threshold mask (about 130 MB and 65 MB at \eqn{m = 4000}). A matrix
#' with missing values takes the pairwise-complete \code{cor()} path instead,
#' which is far slower; \code{\link{extract_variant_set}} imputes before this
#' step, so the scan never does.
#'
#' The linear-dependence step is a blocked in-order Gram-Schmidt applied twice
#' (\code{.independent_columns()}), \eqn{O(N r m)} in BLAS calls for \eqn{r}
#' kept columns among the \eqn{m} visited, and it stops once \eqn{r} reaches
#' the number of rows with a nonzero entry. It keeps the columns that
#' LINPACK \code{qr()} keeps, which shifted every remaining column once per
#' dropped column and so dominated genes with more variants than samples.
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
  n <- nrow(G)
  if (p <= 1) return(seq_len(p))

  # Compute MAF for tie-breaking (needed for "lower_maf" and "higher_maf")
  MAF <- if (prefer_keep != "first") {
    af <- colMeans(G, na.rm = TRUE) / 2
    pmin(af, 1 - af)
  } else {
    NULL
  }

  has_na <- anyNA(G)

  # Handle zero-variance columns (monomorphic after imputation). The two-pass
  # column variance from the centered matrix equals var(); a constant column
  # gives exactly 0 either way. The centered matrix is reused below for the
  # correlation, so the NA-free path centers once.
  if (has_na) {
    col_vars <- apply(G, 2, var, na.rm = TRUE)
  } else {
    col_means <- colMeans(G)
    Gc <- G - rep(col_means, each = n)     # column-major: each mean repeated n times
    col_vars <- colSums(Gc * Gc) / (n - 1)
  }
  zero_var <- col_vars < .Machine$double.eps

  if (all(zero_var)) return(integer(0))

  # Only compute correlation among variable columns
  var_idx <- which(!zero_var)
  if (length(var_idx) <= 1) {
    keep <- var_idx
  } else {
    if (has_na) {
      # Missing values: the pairwise-complete path (slow, kept for direct callers
      # that did not impute). The scan never reaches it.
      R <- cor(G[, var_idx, drop = FALSE], use = "pairwise.complete.obs")
    } else {
      # Correlation as the cross-product of the standardized columns,
      # Z_ij = (G_ij - mean_j) / (s_j * sqrt(n - 1)), so R = t(Z) %*% Z. One BLAS
      # call; equals cor() to the rounding of the summation order.
      Z <- Gc[, var_idx, drop = FALSE]
      rm(Gc)
      Z <- Z / rep(sqrt(col_vars[var_idx] * (n - 1)), each = n)
      R <- crossprod(Z)
      rm(Z)
    }

    # Greedy removal with MAF-based tie-breaking
    # (replaces caret::findCorrelation)
    var_MAF <- if (!is.null(MAF)) MAF[var_idx] else NULL
    to_remove <- .greedy_ld_prune(R, ld_threshold, var_MAF, prefer_keep)
    rm(R)

    keep_among_var <- setdiff(seq_along(var_idx), to_remove)
    keep <- var_idx[keep_among_var]
  }

  # Linear dependence removal, in column order
  # After LD pruning, some columns may still be exact linear combinations
  # of others (e.g., after ultra-rare collapsing). .independent_columns()
  # keeps the first maximal linearly independent subset in column order,
  # the set LINPACK qr() returned here before (same rule, same tolerance).
  if (remove_lindep && length(keep) > 1) {
    G_sub <- G[, keep, drop = FALSE]
    if (anyNA(G_sub)) {
      stop("filter_variants_ld: remove_lindep = TRUE requires a matrix without missing values")
    }
    keep <- keep[.independent_columns(G_sub)]
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
#' The adjacency (|r| above the threshold, diagonal FALSE) is formed once. Each
#' variant's count of partners among the variants not yet removed is kept up to
#' date: removing variant w zeroes its own count and decrements the count of
#' every non-removed partner of w. That count equals the column sum of the
#' adjacency with the removed rows and columns zeroed, which is what a full
#' recomputation per removal gave, so the candidate set and the choice at every
#' step are the same, at O(m^2 + k m) instead of O(k m^2).
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
  high_ld <- abs(R) > threshold   # the adjacency, thresholded once
  diag(high_ld) <- FALSE          # ignore self-correlation
  n_partners <- colSums(high_ld)  # partners among the not-yet-removed variants
  removed <- logical(p)

  repeat {
    max_partners <- max(n_partners)
    if (max_partners == 0) break

    # Find candidates with the most high-LD partners (removed variants hold 0)
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
    # The removed variant leaves every partner's count; its own count goes to 0.
    partners <- which(high_ld[, worst] & !removed)
    n_partners[partners] <- n_partners[partners] - 1L
    n_partners[worst] <- 0L
  }

  which(removed)
}


#' Linearly Independent Columns, in Column Order
#'
#' Returns the indices of the columns of \code{G} that LINPACK \code{qr()} with
#' its limited column pivoting keeps: the columns are visited in order, and
#' column \eqn{j} is kept when its residual after projection onto the kept
#' columns before it has norm at least \code{tol} times the column's own norm.
#' The projection is a blocked classical Gram-Schmidt applied twice (CGS2):
#' a block of 256 columns is projected onto the kept basis in BLAS level-3
#' calls, then sub-blocks of 16 onto the block's own kept columns, then each
#' column onto its sub-block's, so the residual norms are accurate to rounding.
#' Rows with no nonzero entry are dropped first, which changes no residual, and
#' the visit stops once the kept count reaches the remaining row count (the
#' rank bound), after which every column is dependent.
#'
#' \code{qr()} applies the same rule, but cycles each dropped column to the
#' end by shifting every remaining column one slot, a cost of order
#' \eqn{N \times} dropped \eqn{\times} remaining that reached 367 s on a gene
#' with 14,038 columns and 3,195 kept at \eqn{N = 3{,}202}. Here the cost is
#' \eqn{O(N r m)} for \eqn{r} kept columns among the \eqn{m} visited.
#'
#' Both procedures estimate the same residual norm to rounding error, so they
#' differ only for a column whose residual sits within rounding of
#' \code{tol} times its norm. Integer dosage columns are either exactly
#' dependent (residual at machine precision) or far from it.
#'
#' @param G Numeric matrix without missing values whose columns all have a
#'   nonzero norm (\code{filter_variants_ld()} excludes zero-variance columns
#'   before this step). A zero column, should one arrive, is dropped.
#' @param tol Relative tolerance, \code{qr()}'s default \code{1e-07}.
#'
#' @return Integer vector of the kept column indices, increasing.
#'
#' @keywords internal
#' @noRd
.independent_columns <- function(G, tol = 1e-7) {
  block <- 256L   # columns projected onto the kept basis in one BLAS call
  sub <- 16L      # columns projected onto the block's kept columns in one call
  norm0 <- sqrt(colSums(G * G))          # the original column norms (qr()'s reference)
  nz <- rowSums(G != 0) > 0
  if (!all(nz)) G <- G[nz, , drop = FALSE]   # zero rows carry no dependence information
  n <- nrow(G); m <- ncol(G)
  keep <- logical(m)
  if (n == 0L || m == 0L) return(integer(0))

  # CGS2 against a list of orthonormal blocks: two passes keep the residual
  # orthogonal to the basis to rounding, whatever the conditioning of G.
  project <- function(B, blocks) {
    if (length(blocks)) for (pass in 1:2) for (Qk in blocks) B <- B - Qk %*% crossprod(Qk, B)
    B
  }

  basis <- list()   # orthonormal blocks of the kept columns, in visiting order
  r <- 0L           # number of kept columns so far
  start <- 1L
  while (start <= m && r < n) {         # r == n: rank bound, the rest is dependent
    idx <- start:min(start + block - 1L, m)
    B <- project(G[, idx, drop = FALSE], basis)
    inblock <- list(); nb <- 0L; s0 <- 1L
    while (s0 <= length(idx) && r + nb < n) {
      sidx <- s0:min(s0 + sub - 1L, length(idx))
      Bs <- project(B[, sidx, drop = FALSE], inblock)
      qs <- NULL                         # kept columns of this sub-block, orthonormal
      for (jj in seq_along(sidx)) {      # in order within the sub-block
        v <- Bs[, jj]
        if (!is.null(qs)) for (pass in 1:2) v <- v - qs %*% crossprod(qs, v)
        nv <- sqrt(sum(v * v))
        j <- idx[sidx[jj]]
        if (norm0[j] > 0 && nv >= tol * norm0[j]) {   # the test of LINPACK dqrdc2
          qs <- cbind(qs, v / nv)
          keep[j] <- TRUE
          nb <- nb + 1L
          if (r + nb >= n) break
        }
      }
      if (!is.null(qs)) inblock[[length(inblock) + 1L]] <- qs
      s0 <- sidx[length(sidx)] + 1L
    }
    if (nb > 0L) {
      basis[[length(basis) + 1L]] <- do.call(cbind, inblock)
      r <- r + nb
    }
    start <- idx[length(idx)] + 1L
  }
  which(keep)
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
