########## PC Computation from GDS Genotype Data ##########
#
# This file provides functions to compute principal components from GDS files
# for population stratification adjustment. Uses SNPRelate for LD pruning and
# PCA computation.
#
# EXPORTED FUNCTIONS:
#   - compute_pcs_gds(): Compute PCs from per-chromosome GDS files
#   - print.glow_pcs(): S3 print method for glow_pcs objects
#   - plot.glow_pcs(): S3 plot method (dispatches to scree or pairs)
#
# INTERNAL FUNCTIONS:
#   - .plot_scree(): Scree plot helper (eigenvalue bars + cumulative variance)
#   - .plot_pc_pairs(): PC scatter matrix helper
#   - .write_pcs_csv(): CSV output helper with metadata header
#   - .check_pca_dependencies(): Dependency checker for SNPRelate and SeqArray

#################### EXPORTED MAIN FUNCTIONS ####################

#' Compute Principal Components from GDS Genotype Data
#'
#' @description
#' Computes principal components from per-chromosome SeqArray GDS files for
#' population stratification adjustment in association analyses. Variants are
#' filtered by MAF and missing rate, LD pruned, then PCA is performed via
#' SNPRelate.
#'
#' @param gds_files Character vector. Paths to one or more SeqArray GDS files,
#'   typically one per chromosome. All files must exist and be readable.
#' @param n_pcs Integer. Number of principal components to compute. Must be
#'   a positive integer (\eqn{\geq 1}). Default: 10.
#' @param maf_threshold Numeric. Minimum minor allele frequency for variant
#'   inclusion. Must be in \eqn{[0, 0.5]}. Default: 0.05.
#' @param missing_rate Numeric. Maximum per-variant missing rate. Variants
#'   with missing rate above this threshold are excluded. Must be in
#'   \eqn{[0, 1]}. Default: 0.05.
#' @param ld_threshold Numeric. LD pruning \eqn{r^2} threshold. Variants
#'   with pairwise \eqn{r^2} above this threshold are pruned. Must be in
#'   (0, 1]. Default: 0.2.
#' @param output_file Character or NULL. Optional path to save results.
#'   Extension determines format: \code{.rds} uses \code{saveRDS()} (preserves
#'   S3 class and all attributes), \code{.csv} writes a CSV with metadata
#'   comment header. Default: NULL (no file output).
#' @param num_thread Integer. Number of threads for SNPRelate operations.
#'   Default: 1L.
#' @param algorithm Character. PCA algorithm: \code{"randomized"} (default) uses
#'   randomized truncated SVD—fast and memory-efficient for large datasets
#'   (\eqn{O(npk)} time, \eqn{O(nk)} memory). \code{"exact"} computes the full
#'   \eqn{n \times n} GRM then extracts eigenvalues—requires \eqn{O(n^2)} memory
#'   (e.g., ~18 GB for 48K samples). Use \code{"exact"} only for small datasets
#'   or when exact eigendecomposition is needed. Default: \code{"randomized"}.
#' @param seed Integer or NULL. Random seed for reproducibility of LD pruning.
#'   Default: NULL (no seed set).
#' @param verbose Integer. Verbosity level: 0 = silent, 1 = progress messages,
#'   2 = detailed messages. Default: 1.
#'
#' @return An object of class \code{c("glow_pcs", "data.frame")} with columns:
#'   \describe{
#'     \item{sample.id}{Character. Sample identifiers from the GDS file(s).}
#'     \item{PC1, PC2, ..., PC\{n_pcs\}}{Numeric. Principal component scores.}
#'   }
#'   The following attributes are attached:
#'   \describe{
#'     \item{eigenvalues}{Numeric vector of length \code{n_pcs}.}
#'     \item{variance_proportion}{Numeric vector. Proportion of variance
#'       explained by each PC.}
#'     \item{total_variance_explained}{Numeric scalar. Sum of variance
#'       proportions for all returned PCs.}
#'     \item{n_variants_used}{Integer. Number of variants after filtering and
#'       LD pruning.}
#'     \item{n_variants_per_chr}{Named integer vector. Pruned variant count
#'       per chromosome.}
#'     \item{call_args}{List. Arguments used for this call (for
#'       reproducibility).}
#'   }
#'
#' @details
#' \strong{Processing Pipeline:}
#'
#' \enumerate{
#'   \item Validate inputs and check SNPRelate/SeqArray availability
#'   \item Merge chromosome GDS files via \code{SeqArray::seqMerge()} (if
#'     multiple files; skipped for single file)
#'   \item Convert SeqArray GDS to SNP GDS via \code{SeqArray::seqGDS2SNP()}
#'   \item LD prune via \code{SNPRelate::snpgdsLDpruning()} with MAF and
#'     missing rate filtering
#'   \item Compute PCA via \code{SNPRelate::snpgdsPCA()} on pruned variants
#'   \item Build \code{glow_pcs} S3 object with eigenvalue diagnostics
#'   \item Save to disk if \code{output_file} is specified
#'   \item Clean up temporary files (via \code{on.exit()})
#' }
#'
#' \strong{GDS Format Compatibility:}
#'
#' SNPRelate functions (\code{snpgdsPCA}, \code{snpgdsLDpruning}) require
#' SNP GDS format (\code{SNPGDSFileClass}), not SeqArray format
#' (\code{SeqVarGDSClass}). This function handles the conversion
#' automatically via \code{seqGDS2SNP()}.
#'
#' \strong{Computational Complexity:}
#'
#' For \eqn{n} samples and \eqn{p} variants after pruning, PCA computation
#' is \eqn{O(n^2 p)} for the covariance matrix. LD pruning is
#' \eqn{O(p \cdot w)} where \eqn{w} is the sliding window size.
#'
#' \strong{Dependencies:}
#'
#' SNPRelate and SeqArray are optional dependencies (in \code{Suggests}).
#' The function checks availability at runtime and provides installation
#' instructions if missing.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic PC computation from per-chromosome GDS files
#' gds_files <- sprintf("data/gds/chr%d.gds", 1:22)
#' pcs <- compute_pcs_gds(gds_files, n_pcs = 10)
#' print(pcs)
#'
#' # Example 2: Compute and save with seed for reproducibility
#' pcs <- compute_pcs_gds(
#'   gds_files = gds_files,
#'   n_pcs = 20,
#'   output_file = "results/pcs.rds",
#'   seed = 42L,
#'   num_thread = 4L,
#'   verbose = 1
#' )
#'
#' # Example 3: Scree plot to determine optimal number of PCs
#' plot(pcs, type = "scree")
#'
#' # Example 4: PC scatter matrix for population structure
#' plot(pcs, type = "pairs", n_pairs = 6)
#'
#' # Example 5: Single GDS file (skips merge step)
#' pcs <- compute_pcs_gds("data/merged.gds", n_pcs = 10)
#' }
#'
#' @references
#' Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012).
#' A high-performance computing toolset for relatedness and principal
#' component analysis of SNP data. Bioinformatics, 28(24), 3326-3328.
#' DOI: 10.1093/bioinformatics/bts606
#'
#' @seealso
#' \code{\link{plink_to_gds}} for converting PLINK files to GDS format
#'
#' @export
compute_pcs_gds <- function(gds_files,
                            n_pcs = 10L,
                            maf_threshold = 0.05,
                            missing_rate = 0.05,
                            ld_threshold = 0.2,
                            output_file = NULL,
                            num_thread = 1L,
                            algorithm = c("randomized", "exact"),
                            seed = NULL,
                            verbose = 1) {

  # ========== Step 1: Validate Inputs ==========
  # (Validate before dependency check so invalid args are caught even
  #  without SNPRelate installed)

  # gds_files: character vector, all files must exist

  if (!is.character(gds_files) || length(gds_files) < 1L) {
    stop("gds_files must be a character vector of one or more GDS file paths.",
         call. = FALSE)
  }
  missing_files <- gds_files[!file.exists(gds_files)]
  if (length(missing_files) > 0L) {
    stop("GDS file(s) not found:\n  ",
         paste(missing_files, collapse = "\n  "),
         call. = FALSE)
  }

  # n_pcs: positive integer
  n_pcs <- as.integer(n_pcs)
  if (is.na(n_pcs) || n_pcs < 1L) {
    stop("n_pcs must be a positive integer (>= 1), got: ", n_pcs,
         call. = FALSE)
  }

  # maf_threshold: numeric in [0, 0.5]
  if (!is.numeric(maf_threshold) || length(maf_threshold) != 1L ||
      is.na(maf_threshold) || maf_threshold < 0 || maf_threshold > 0.5) {
    stop("maf_threshold must be a numeric value in [0, 0.5], got: ",
         maf_threshold, call. = FALSE)
  }

  # missing_rate: numeric in [0, 1]
  if (!is.numeric(missing_rate) || length(missing_rate) != 1L ||
      is.na(missing_rate) || missing_rate < 0 || missing_rate > 1) {
    stop("missing_rate must be a numeric value in [0, 1], got: ",
         missing_rate, call. = FALSE)
  }

  # ld_threshold: numeric in (0, 1]
  if (!is.numeric(ld_threshold) || length(ld_threshold) != 1L ||
      is.na(ld_threshold) || ld_threshold <= 0 || ld_threshold > 1) {
    stop("ld_threshold must be a numeric value in (0, 1], got: ",
         ld_threshold, call. = FALSE)
  }

  # output_file: if provided, directory must exist and extension must be .rds or .csv
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1L) {
      stop("output_file must be a single character string.", call. = FALSE)
    }
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir)) {
      stop("Output directory does not exist: ", out_dir, call. = FALSE)
    }
    out_ext <- tolower(tools::file_ext(output_file))
    if (!out_ext %in% c("rds", "csv")) {
      stop("output_file must have extension .rds or .csv, got: .", out_ext,
           call. = FALSE)
    }
  }

  # num_thread: positive integer
  num_thread <- as.integer(num_thread)
  if (is.na(num_thread) || num_thread < 1L) {
    stop("num_thread must be a positive integer (>= 1).", call. = FALSE)
  }

  # algorithm: match.arg
  algorithm <- match.arg(algorithm)

  # seed: integer or NULL
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (is.na(seed)) {
      stop("seed must be an integer or NULL.", call. = FALSE)
    }
  }

  # verbose: integer 0, 1, or 2
  verbose <- as.integer(verbose)
  if (is.na(verbose) || !verbose %in% c(0L, 1L, 2L)) {
    stop("verbose must be 0, 1, or 2.", call. = FALSE)
  }

  # ========== Step 2: Check Dependencies ==========
  .check_pca_dependencies()

  # ========== Step 3: Initialize temp file tracking ==========

  # Track temp files for cleanup (these will be populated below)
  temp_merged <- NULL
  temp_snp <- NULL
  snp_gds <- NULL

  # Register cleanup on exit (runs even on error)
  on.exit({
    # Close SNP GDS connection if open
    if (!is.null(snp_gds)) {
      try(SNPRelate::snpgdsClose(snp_gds), silent = TRUE)
    }
    # Remove temp merged SeqArray file
    if (!is.null(temp_merged)) {
      unlink(temp_merged, force = TRUE)
    }
    # Remove temp SNP GDS file
    if (!is.null(temp_snp)) {
      unlink(temp_snp, force = TRUE)
    }
  }, add = TRUE)

  # ========== Step 4: Merge Chromosome Files (if multiple) ==========

  if (verbose >= 1L) {
    message("=== Computing Principal Components ===")
    message(sprintf("  Input: %d GDS file(s)", length(gds_files)))
  }

  if (length(gds_files) > 1L) {
    if (verbose >= 1L) {
      message("  Step 1/4: Merging chromosome files...")
    }
    temp_merged <- tempfile(pattern = "glow_pcs_merged_", fileext = ".gds")
    SeqArray::seqMerge(gds_files, temp_merged, verbose = (verbose >= 2L))
    seq_file <- temp_merged

    if (verbose >= 1L) {
      message("    Merge complete.")
    }
  } else {
    if (verbose >= 1L) {
      message("  Step 1/4: Single GDS file, skipping merge.")
    }
    seq_file <- gds_files[1L]
  }

  # ========== Step 5: Convert SeqArray to SNP GDS ==========

  if (verbose >= 1L) {
    message("  Step 2/4: Converting to SNP GDS format...")
  }

  temp_snp <- tempfile(pattern = "glow_pcs_snp_", fileext = ".gds")
  SeqArray::seqGDS2SNP(seq_file, temp_snp, verbose = (verbose >= 2L))
  snp_gds <- SNPRelate::snpgdsOpen(temp_snp)

  if (verbose >= 2L) {
    message("    SNP GDS opened successfully.")
  }

  # ========== Step 6: LD Pruning ==========

  if (verbose >= 1L) {
    message("  Step 3/4: LD pruning (MAF >= ", maf_threshold,
            ", missing <= ", missing_rate,
            ", LD r2 < ", ld_threshold, ")...")
  }

  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  pruned <- SNPRelate::snpgdsLDpruning(
    snp_gds,
    maf = maf_threshold,
    missing.rate = missing_rate,
    ld.threshold = ld_threshold,
    method = "corr",
    slide.max.bp = 500000L,
    num.thread = num_thread,
    verbose = (verbose >= 2L)
  )

  # Collect pruned variant IDs
  pruned_ids <- unlist(pruned, use.names = FALSE)
  n_per_chr <- vapply(pruned, length, integer(1L))

  if (verbose >= 1L) {
    message(sprintf("    %d variants retained after pruning (%d chromosomes)",
                    length(pruned_ids), length(n_per_chr)))
  }
  if (verbose >= 2L) {
    for (chr_name in names(n_per_chr)) {
      message(sprintf("      %s: %d variants", chr_name, n_per_chr[chr_name]))
    }
  }

  # Validate n_pcs against pruned variant count
  if (n_pcs >= length(pruned_ids)) {
    stop("n_pcs (", n_pcs, ") must be less than the number of pruned variants (",
         length(pruned_ids), "). Reduce n_pcs or relax filtering thresholds.",
         call. = FALSE)
  }

  # ========== Step 7: Compute PCA ==========

  if (verbose >= 1L) {
    message(sprintf("  Step 4/4: Computing %d PCs (algorithm: %s)...",
                    n_pcs, algorithm))
  }

  pca <- SNPRelate::snpgdsPCA(
    snp_gds,
    snp.id = pruned_ids,
    eigen.cnt = n_pcs,
    algorithm = algorithm,
    num.thread = num_thread,
    autosome.only = FALSE,
    verbose = (verbose >= 2L)
  )

  if (verbose >= 1L) {
    message("    PCA complete.")
  }

  # ========== Step 8: Build glow_pcs Object ==========

  # Build data.frame with sample IDs and PC scores
  pc_df <- data.frame(
    sample.id = pca$sample.id,
    pca$eigenvect[, seq_len(n_pcs), drop = FALSE],
    stringsAsFactors = FALSE
  )
  colnames(pc_df) <- c("sample.id", paste0("PC", seq_len(n_pcs)))

  # Attach S3 class
  class(pc_df) <- c("glow_pcs", "data.frame")

  # Attach attributes
  attr(pc_df, "eigenvalues") <- pca$eigenval[seq_len(n_pcs)]
  attr(pc_df, "variance_proportion") <- pca$varprop[seq_len(n_pcs)]
  attr(pc_df, "total_variance_explained") <- sum(pca$varprop[seq_len(n_pcs)])
  attr(pc_df, "n_variants_used") <- length(pruned_ids)
  attr(pc_df, "n_variants_per_chr") <- n_per_chr
  attr(pc_df, "call_args") <- list(
    gds_files = gds_files,
    n_pcs = n_pcs,
    maf_threshold = maf_threshold,
    missing_rate = missing_rate,
    ld_threshold = ld_threshold,
    algorithm = algorithm,
    seed = seed
  )

  # ========== Step 9: Save to Disk (if requested) ==========

  if (!is.null(output_file)) {
    ext <- tolower(tools::file_ext(output_file))
    if (ext == "rds") {
      saveRDS(pc_df, file = output_file)
      if (verbose >= 1L) {
        message(sprintf("  Saved to RDS: %s", output_file))
      }
    } else if (ext == "csv") {
      .write_pcs_csv(pc_df, output_file)
      if (verbose >= 1L) {
        message(sprintf("  Saved to CSV: %s", output_file))
      }
    }
  }

  # ========== Step 10: Summary ==========

  if (verbose >= 1L) {
    message(sprintf(
      "\n=== PC Computation Complete ===\n  Samples: %d\n  PCs: %d\n  Variants used: %d\n  Variance explained: %.1f%%",
      nrow(pc_df), n_pcs, length(pruned_ids),
      sum(pca$varprop[seq_len(n_pcs)]) * 100
    ))
  }

  return(pc_df)
}


#' Print Method for glow_pcs Objects
#'
#' @description
#' Displays a compact summary of a \code{glow_pcs} object including sample
#' count, number of PCs, variants used, and variance explained.
#'
#' @param x A \code{glow_pcs} object returned by \code{\link{compute_pcs_gds}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' pcs <- compute_pcs_gds(gds_files, n_pcs = 10)
#' print(pcs)
#' }
#'
#' @export
print.glow_pcs <- function(x, ...) {
  n_pcs <- ncol(x) - 1L
  vp <- attr(x, "variance_proportion")
  cat("GLOW Principal Components\n")
  cat("  Samples:   ", nrow(x), "\n")
  cat("  PCs:       ", n_pcs, "\n")
  cat("  Variants:  ", attr(x, "n_variants_used"),
      "(after MAF/missing/LD filtering)\n")
  cat("  Variance:  ",
      sprintf("%.1f%%", attr(x, "total_variance_explained") * 100),
      sprintf("(cumulative, %d PCs)\n", n_pcs))
  # Show top 3 PCs
  top_n <- min(3L, n_pcs)
  if (!is.null(vp) && length(vp) >= top_n) {
    top_str <- paste(
      sprintf("PC%d=%.1f%%", seq_len(top_n), vp[seq_len(top_n)] * 100),
      collapse = ", "
    )
    cat("  Top PCs:   ", top_str, "\n")
  }
  invisible(x)
}


#' Plot Method for glow_pcs Objects
#'
#' @description
#' Produces diagnostic plots for a \code{glow_pcs} object. Dispatches to
#' either a scree plot or a PC scatter matrix based on the \code{type}
#' argument.
#'
#' @param x A \code{glow_pcs} object returned by \code{\link{compute_pcs_gds}}.
#' @param type Character. Plot type: \code{"scree"} for a scree plot with
#'   individual eigenvalues (bars) and cumulative variance (line), or
#'   \code{"pairs"} for a scatter plot matrix of the first \code{n_pairs}
#'   PCs. Default: \code{"scree"}.
#' @param n_pairs Integer. Number of PCs to include in the pairs plot.
#'   Only used when \code{type = "pairs"}. Default: 4.
#' @param ... Additional graphical parameters passed to the underlying
#'   plot functions.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \dontrun{
#' pcs <- compute_pcs_gds(gds_files, n_pcs = 20)
#'
#' # Scree plot
#' plot(pcs, type = "scree")
#'
#' # PC scatter matrix
#' plot(pcs, type = "pairs", n_pairs = 6)
#' }
#'
#' @export
plot.glow_pcs <- function(x, type = c("scree", "pairs"), n_pairs = 4L, ...) {
  type <- match.arg(type)
  if (type == "scree") {
    .plot_scree(x, ...)
  } else {
    .plot_pc_pairs(x, n_pairs = n_pairs, ...)
  }
  invisible(x)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Scree Plot for glow_pcs Object
#'
#' Creates a scree plot with individual variance proportions (bars) and
#' cumulative variance (line overlay). Useful for determining the optimal
#' number of PCs to retain.
#'
#' @param x A \code{glow_pcs} object.
#' @param ... Additional graphical parameters passed to \code{barplot()}.
#'
#' @return Invisibly returns NULL.
#'
#' @keywords internal
#' @noRd
.plot_scree <- function(x, ...) {
  eigenvalues <- attr(x, "eigenvalues")
  vp <- attr(x, "variance_proportion")
  cum_var <- cumsum(vp)
  n <- length(eigenvalues)
  pc_idx <- seq_len(n)

  # Bar plot: individual variance proportion
  bp <- barplot(
    vp * 100,
    names.arg = paste0("PC", pc_idx),
    ylab = "Variance Explained (%)",
    xlab = "Principal Component",
    main = "PCA Scree Plot",
    col = "steelblue",
    border = NA,
    ylim = c(0, max(vp * 100) * 1.3),
    ...
  )

  # Overlay: cumulative variance line
  lines(bp, cum_var * 100, type = "b", pch = 19, col = "darkred", lwd = 2)

  # Legend
  legend("topright",
    legend = c("Individual", "Cumulative"),
    col = c("steelblue", "darkred"),
    pch = c(15, 19), lty = c(NA, 1), lwd = c(NA, 2),
    bty = "n"
  )

  # Annotate cumulative total at last PC
  text(bp[n], cum_var[n] * 100,
       labels = sprintf("%.1f%%", cum_var[n] * 100),
       pos = 3, col = "darkred", cex = 0.8)

  invisible(NULL)
}


#' PC Scatter Matrix Plot
#'
#' Creates a scatter plot matrix of the first \code{n_pairs} PCs, with
#' variance explained shown in axis labels.
#'
#' @param x A \code{glow_pcs} object.
#' @param n_pairs Integer. Number of PCs to include (capped at available PCs).
#' @param ... Additional graphical parameters passed to \code{pairs()}.
#'
#' @return Invisibly returns NULL.
#'
#' @keywords internal
#' @noRd
.plot_pc_pairs <- function(x, n_pairs = 4L, ...) {
  n_pcs <- ncol(x) - 1L
  n_pairs <- min(n_pairs, n_pcs)
  pc_cols <- paste0("PC", seq_len(n_pairs))
  vp <- attr(x, "variance_proportion")

  # Labels with variance explained
  labels <- sprintf("PC%d (%.1f%%)", seq_len(n_pairs),
                    vp[seq_len(n_pairs)] * 100)

  pairs(
    x[, pc_cols, drop = FALSE],
    labels = labels,
    pch = ".",
    col = grDevices::adjustcolor("steelblue", alpha.f = 0.3),
    main = "PC Scatter Matrix",
    ...
  )

  invisible(NULL)
}


#' Write PCs to CSV with Metadata Header
#'
#' Writes PC scores to a CSV file with metadata comment lines (starting
#' with \code{#}) above the data rows. Comment lines include date, sample
#' count, number of PCs, variants used, total variance explained, and
#' eigenvalues.
#'
#' @param pc_df A \code{glow_pcs} object (data.frame subclass).
#' @param path Character. Output file path.
#'
#' @return Invisibly returns NULL. Side effect: writes file to disk.
#'
#' @keywords internal
#' @noRd
.write_pcs_csv <- function(pc_df, path) {
  # Write a clean CSV (no comment headers) so that read.csv() works with
  # default settings. Metadata is saved to a companion _meta.txt file.
  write.csv(pc_df, file = path, row.names = FALSE)

  # Write metadata sidecar
  meta_path <- sub("\\.csv$", "_meta.txt", path)
  meta <- c(
    sprintf("GLOW Principal Components"),
    sprintf("Date: %s", Sys.Date()),
    sprintf("Samples: %d", nrow(pc_df)),
    sprintf("PCs: %d", ncol(pc_df) - 1L),
    sprintf("Variants used: %d", attr(pc_df, "n_variants_used")),
    sprintf("Total variance explained: %.4f",
            attr(pc_df, "total_variance_explained")),
    sprintf("Eigenvalues: %s",
            paste(round(attr(pc_df, "eigenvalues"), 6), collapse = ", "))
  )
  writeLines(meta, meta_path)

  invisible(NULL)
}


#' Check PCA Dependencies (SNPRelate and SeqArray)
#'
#' Verifies that SNPRelate and SeqArray packages are available. Provides
#' informative installation instructions if any are missing.
#'
#' @return Invisibly returns TRUE if all dependencies are present. Throws
#'   an error with installation instructions otherwise.
#'
#' @keywords internal
#' @noRd
.check_pca_dependencies <- function() {
  missing_pkgs <- character()
  if (!requireNamespace("SNPRelate", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "SNPRelate")
  }
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "SeqArray")
  }
  if (length(missing_pkgs) > 0L) {
    stop(
      "Required package(s) not available: ",
      paste(missing_pkgs, collapse = ", "), "\n",
      "Install with: BiocManager::install(c(",
      paste0("'", missing_pkgs, "'", collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
