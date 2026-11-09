# ==============================================================================
# In-sample LD Scores (for LD Score Regression)
# ==============================================================================
#
# Per-variant LD scores ell_j = sum_k r2tilde_jk over biallelic SNVs k within
# +/- `window` bp of j on the same chromosome, using the analysis cohort's own
# genotypes. r2tilde = r^2 - (1 - r^2)/(n - 2) is the LDSC-unbiased r^2
# estimator. LD scores are genotype-only (GWAS-config-independent), so they are
# computed once and reused by ldsc_regression() across covariate configs.
#
# EXPORTED FUNCTIONS:
#   - compute_ld_scores()   Single-chromosome in-sample LD scores
#
# This is a packaged port of the validated research script
# ldsc-compute-ldscores.R. It uses SeqArray (a Suggested dependency, guarded by
# requireNamespace) and GLOWr's own SNV predicate .is_snv_from_alleles() rather
# than SeqVarTools; the LD score is a correlation, so it is invariant to whether
# the reference- or alternate-allele dosage is used.


#' In-sample LD Scores for One Chromosome
#'
#' Computes, for each biallelic SNV on a single chromosome, the LD score
#' \eqn{\ell_j = \sum_k \tilde r^2_{jk}} summed over biallelic SNVs \eqn{k}
#' within \eqn{\pm}\code{window} bp of \eqn{j}, from the cohort's own genotypes.
#' This is the LD-score input to \code{\link{ldsc_regression}}.
#'
#' @details
#' The unbiased \eqn{r^2} estimator (Bulik-Sullivan et al. 2015) is used,
#' \deqn{\tilde r^2_{jk} = r^2_{jk} - \frac{1 - r^2_{jk}}{n - 2},}
#' so the LD score is \eqn{\ell_j = \sum_{k:\,|\mathrm{pos}_k-\mathrm{pos}_j|\le w}
#' \tilde r^2_{jk}}. Variants are processed in contiguous core segments of
#' \code{segment} variants; for each segment the genotypes of
#' \eqn{[\mathrm{segment} \pm \mathrm{window}]} are read, mean-imputed,
#' column-standardised, and a single BLAS cross-product gives the
#' window-by-core correlation block, from which each core variant's LD score is
#' a \eqn{\pm}\code{window} prefix-sum over \eqn{\tilde r^2}, the window inclusive
#' at both ends. The file's storage order need not be position order: the genotype
#' columns are put in position order after they are read. Peak memory is bounded by one
#' segment plus its flanks, never a whole chromosome: per worker, the standardized
#' genotypes of the segment and its two flanks (8 bytes per sample per variant) plus
#' three matrices of the flanked window by the segment. With \code{n_cores > 1} the
#' segments are split into contiguous blocks and each block runs in a forked worker on
#' its own GDS handle, each holding its own matrices; a segment always reads its full
#' flank, so the result does not depend on \code{n_cores}.
#'
#' \strong{Single-chromosome contract.} Positions repeat across chromosomes, so
#' the \eqn{\pm}\code{window} mask is only meaningful within one chromosome;
#' the function therefore requires \code{gds_file} to hold a single chromosome
#' (and errors otherwise). For a genome-wide table, call it per chromosome and
#' \code{rbind} (mirrors how \code{\link{marginal_scan}} is driven per-chr).
#'
#' \strong{LD scores are a population property} estimated from genotypes only;
#' they do not depend on the phenotype or covariates, so one cached table serves
#' every covariate config. \code{sample_n > 0} uses a fixed (seeded) random
#' reference subsample of individuals for speed; the unbiased estimator absorbs
#' the resulting small-\eqn{n} bias.
#'
#' \strong{Ancestry caveat.} In an ancestry-mixed cohort, in-sample LD scores
#' carry some cross-ancestry LD. The LDSC \emph{intercept} (the confounding
#' read-out that \code{\link{ldsc_regression}} returns) is robust to LD-score
#' mis-scaling, which the slope absorbs.
#'
#' @param gds_file Path to a GDS/aGDS file for \strong{one} chromosome. Opened
#'   read-only and closed on exit.
#' @param variant_id Optional integer vector of \code{variant.id}s to restrict
#'   the variant set (both targets and LD-score neighbours). Default \code{NULL}
#'   uses all biallelic SNVs in the file.
#' @param window Integer LD-score window half-width in bp (default \code{1e6},
#'   i.e. \eqn{\pm}1 Mb).
#' @param segment Integer number of core variants processed per BLAS block
#'   (default \code{2000L}). A memory/speed knob only; results are invariant.
#' @param sample_n Integer. If \code{> 0}, a fixed random subsample of that many
#'   individuals is used as the LD reference (default \code{0L} = all samples).
#' @param seed Integer seed for the \code{sample_n} subsample (default
#'   \code{1L}); unused when \code{sample_n = 0L}.
#' @param verbose Integer verbosity (0 = silent, 1 = a one-line summary).
#' @param n_cores Integer (default \code{1L}). Forked workers over which the core segments
#'   are spread (\code{parallel::mclapply}). Each worker opens its own read-only handle on
#'   \code{gds_file} and computes whole segments with their flanking windows, so the table
#'   is the one the serial loop writes. Not available on Windows (no fork): a value above 1
#'   warns and runs with one worker.
#'
#' @return A \code{data.frame} with one row per biallelic SNV:
#'   \code{chr}, \code{pos}, \code{ref}, \code{alt}, \code{variant_id}, and
#'   \code{ld} (the LD score), ordered by \code{pos}.
#'
#' @seealso \code{\link{ldsc_regression}}, \code{\link{marginal_scan}}
#'
#' @references
#' Bulik-Sullivan B. et al. (2015). LD Score regression distinguishes
#' confounding from polygenicity in GWAS. \emph{Nat Genet} 47:291.
#'
#' @examples
#' \dontrun{
#' ld22 <- compute_ld_scores("chr22_essentialdb.gds")
#' # Genome-wide: loop per-chromosome GDS files and bind.
#' ld <- do.call(rbind, lapply(sprintf("chr%d.gds", 1:22), compute_ld_scores))
#' }
#'
#' @export
compute_ld_scores <- function(gds_file, variant_id = NULL, window = 1e6,
                              segment = 2000L, sample_n = 0L, seed = 1L,
                              verbose = 1, n_cores = 1L) {
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray is required for compute_ld_scores(). Install with: ",
         "BiocManager::install('SeqArray')", call. = FALSE)
  }
  stopifnot(file.exists(gds_file), window > 0, segment >= 1L)
  WIN <- as.numeric(window)
  SEG <- as.integer(segment)

  gds <- SeqArray::seqOpen(gds_file, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Optional fixed reference subsample of individuals (LD is a population
  # property; the unbiased r^2 estimator handles the small-n bias).
  if (sample_n > 0L) {
    sid <- SeqArray::seqGetData(gds, "sample.id")
    if (sample_n < length(sid)) {
      set.seed(seed)
      SeqArray::seqSetFilter(
        gds, sample.id = sid[sort(sample(length(sid), sample_n))],
        verbose = FALSE)
    }
  }

  # Biallelic SNVs only (2 alleles, both single-nucleotide), optionally
  # restricted to `variant_id`.
  na <- SeqArray::seqGetData(gds, "$num_allele")
  SeqArray::seqSetFilter(gds, variant.sel = (na == 2L), verbose = FALSE)
  allele_str <- SeqArray::seqGetData(gds, "allele")
  snv <- .is_snv_from_alleles(allele_str)
  SeqArray::seqSetFilter(gds, variant.sel = snv, action = "intersect",
                         verbose = FALSE)
  if (!is.null(variant_id)) {
    SeqArray::seqSetFilter(gds, variant.id = variant_id, action = "intersect",
                           verbose = FALSE)
  }

  chr_all <- SeqArray::seqGetData(gds, "chromosome")
  if (length(unique(chr_all)) != 1L) {
    stop("compute_ld_scores() requires a single-chromosome GDS; got ",
         length(unique(chr_all)), " chromosomes (",
         paste(utils::head(unique(chr_all), 5L), collapse = ", "),
         "). Call per chromosome and rbind.", call. = FALSE)
  }

  vid <- SeqArray::seqGetData(gds, "variant.id")
  pos <- SeqArray::seqGetData(gds, "position")
  alleles <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
  ref <- vapply(alleles, `[`, character(1L), 1L)
  alt <- vapply(alleles, `[`, character(1L), 2L)

  o <- order(pos)
  vid <- vid[o]; pos <- pos[o]; ref <- ref[o]; alt <- alt[o]
  m <- length(vid)
  ld <- numeric(m)
  n_samp <- length(SeqArray::seqGetData(gds, "sample.id"))

  if (m == 0L) {
    return(data.frame(chr = character(0), pos = integer(0), ref = character(0),
                      alt = character(0), variant_id = integer(0),
                      ld = numeric(0), stringsAsFactors = FALSE))
  }

  # The reference-sample ids after the optional subsample, so that a worker applies the
  # same filter on its own handle.
  sid_keep <- if (sample_n > 0L) SeqArray::seqGetData(gds, "sample.id") else NULL

  # One core segment [s, e] on the handle `g`: read [pos[s] - WIN, pos[e] + WIN], inclusive
  # at both ends (left.open = TRUE makes the lower bound inclusive for integer and
  # fractional windows alike; the per-target bound below uses the same rule, so a target's
  # score does not depend on where its segment starts), mean-impute, standardize, one BLAS
  # cross-product against the core columns, and sum the adjusted r^2 over each core
  # variant's own window by prefix sums.
  ld_segment <- function(g, s) {
    e   <- min(s + SEG - 1L, m)
    lo  <- findInterval(pos[s] - WIN, pos, left.open = TRUE) + 1L   # first pos >= pos[s]-WIN
    hi  <- findInterval(pos[e] + WIN, pos)                           # last  pos <= pos[e]+WIN
    idx <- lo:hi
    SeqArray::seqSetFilter(g, variant.id = vid[idx], verbose = FALSE)
    D  <- SeqArray::seqGetData(g, "$dosage")     # n_samp x |idx|, in the file's STORAGE
                                                 # order; ref/alt immaterial (r^2 is invariant)
    # The filter selects variants, it does not order them: the columns come in storage
    # order, which is position order only for a position-sorted file. Put them in the
    # order of vid[idx] (position order), which the core columns and the window bounds
    read_ids <- SeqArray::seqGetData(g, "variant.id")
    if (!identical(read_ids, vid[idx])) D <- D[, match(vid[idx], read_ids), drop = FALSE]
    storage.mode(D) <- "double"                  # once; the steps below then work in place
    cm <- colMeans(D, na.rm = TRUE)
    if (anyNA(D)) {                              # impute NA -> col mean
      nas <- which(is.na(D), arr.ind = TRUE)
      D[nas] <- cm[nas[, 2L]]
    }
    csd <- sqrt((colSums(D * D) - n_samp * cm * cm) / (n_samp - 1))
    csd[!is.finite(csd) | csd == 0] <- 1
    # Standardize column by column, in place. sweep() would hold two more n_samp x |idx|
    # copies, and a worker's memory is this matrix times the number of copies (8 bytes per
    # entry: 7,743 samples by 20,000 variants is 1.2 GB per copy).
    for (j in seq_along(idx)) D[, j] <- (D[, j] - cm[j]) / csd[j]
    core_cols <- (s - lo + 1L):(e - lo + 1L)     # cores within the read window
    R  <- crossprod(D, D[, core_cols, drop = FALSE]) / (n_samp - 1)   # |idx| x |core|
    R  <- R * R - (1 - R * R) / (n_samp - 2)     # the LDSC-adjusted r^2, in place
    cs <- apply(R, 2L, cumsum)                   # |idx| x |core| prefix sums
    if (!is.matrix(cs)) cs <- matrix(cs, nrow = nrow(R))   # apply() drops a 1-row result
    wpos  <- pos[idx]; cpos <- pos[s:e]
    a <- findInterval(cpos - WIN, wpos, left.open = TRUE) + 1L  # first window row >= cpos-WIN
    b <- findInterval(cpos + WIN, wpos)                          # last  window row <= cpos+WIN
    q <- seq_along(core_cols)
    # the sum over window rows a..b of column q is cs[b, q] - cs[a - 1, q], with cs[0, q] = 0
    cs[cbind(b, q)] - ifelse(a > 1L, cs[cbind(pmax(a - 1L, 1L), q)], 0)
  }

  seg_starts <- seq(1L, m, by = SEG)
  n_cores <- as.integer(n_cores)
  stopifnot(length(n_cores) == 1L, !is.na(n_cores), n_cores >= 1L)
  if (n_cores > 1L && .Platform$OS.type == "windows") {
    warning("n_cores > 1 needs fork(), which Windows lacks; running with one worker.")
    n_cores <- 1L
  }
  n_workers <- min(n_cores, length(seg_starts))
  if (n_workers == 1L) {
    for (s in seg_starts) {
      e <- min(s + SEG - 1L, m)
      ld[s:e] <- ld_segment(gds, s)
    }
  } else {
    # Contiguous blocks of whole segments, one forked worker each, every worker on its own
    # read-only handle (a forked copy of the parent's handle must not be used). A segment
    # always reads its full flank, so the values do not depend on the split.
    blocks <- split(seg_starts, cut(seq_along(seg_starts), n_workers, labels = FALSE))
    run_block <- function(starts) {
      g <- SeqArray::seqOpen(gds_file, readonly = TRUE, allow.duplicate = TRUE)
      on.exit(SeqArray::seqClose(g), add = TRUE)
      if (!is.null(sid_keep)) {
        SeqArray::seqSetFilter(g, sample.id = sid_keep, verbose = FALSE)
      }
      lapply(starts, function(s) ld_segment(g, s))
    }
    parts <- parallel::mclapply(blocks, run_block, mc.cores = n_workers)
    for (k in seq_along(blocks)) {
      starts <- blocks[[k]]
      if (inherits(parts[[k]], "try-error") || length(parts[[k]]) != length(starts)) {
        stop("compute_ld_scores(): a forked worker failed: ",
             paste(as.character(parts[[k]]), collapse = " "), call. = FALSE)
      }
      for (j in seq_along(starts)) {
        s <- starts[j]; e <- min(s + SEG - 1L, m)
        ld[s:e] <- parts[[k]][[j]]
      }
    }
  }

  if (verbose >= 1) {
    message(sprintf("  chr%s: %d biallelic SNVs, mean LD score %.2f (median %.2f)",
                    chr_all[1L], m, mean(ld), stats::median(ld)))
  }
  data.frame(chr = as.character(chr_all[1L]), pos = pos, ref = ref, alt = alt,
             variant_id = vid, ld = ld, stringsAsFactors = FALSE)
}
