########## Visualization Functions for Marginal Analysis Results ##########
#
# Manhattan plot and QQ plot for genome-wide single-variant association results.
# Includes S3 dispatch for glow_marginal_scan objects.
#
# EXPORTED FUNCTIONS:
#   - plot_manhattan()           Manhattan plot for -log10(p) by genomic position
#   - plot_qq()                  QQ plot with genomic inflation factor (lambda_GC)
#   - plot.glow_marginal_scan()  S3 dispatch: manhattan, qq, or both
#
# INTERNAL HELPERS:
#   - .cumulative_positions()    Compute x-axis positions with chromosome gaps
#   - .chr_to_numeric()          Standardize chromosome labels to integer


#################### EXPORTED MAIN FUNCTIONS ####################

#' Manhattan Plot for Marginal Analysis Results
#'
#' Creates a Manhattan plot showing \eqn{-\log_{10}(p)} for each variant,
#' arranged by genomic position across chromosomes.
#'
#' @param x A glow_marginal_scan data.frame, or any data.frame with
#'   columns for chromosome, position, and p-value.
#' @param pvalue_col Character name of the p-value column. Default: "pvalue".
#' @param chr_col Character name of the chromosome column. Default: "chr".
#' @param pos_col Character name of the position column. Default: "pos".
#' @param suggestive_line Numeric p-value threshold for suggestive significance
#'   line (blue dashed). Default: 1e-5. Set to NULL to suppress.
#' @param genome_wide_line Numeric p-value threshold for genome-wide
#'   significance line (red solid). Default: 5e-8. Set to NULL to suppress.
#' @param highlight Integer or character vector of variant_id values to
#'   highlight in a distinct color. Requires a "variant_id" column in x.
#'   Default: NULL.
#' @param chr_colors Character vector of length 2 for alternating chromosome
#'   colors. Default: c("grey30", "grey60").
#' @param highlight_color Color for highlighted variants. Default: "red".
#' @param title Character plot title. Default: "Manhattan Plot".
#' @param cex Numeric point size. Default: 0.5.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return \code{invisible(NULL)}.
#'
#' @details
#' Chromosome labels are standardized to numeric by stripping the "chr" prefix
#' (e.g., "chr22" becomes 22). Non-numeric chromosomes (X, Y, MT) are assigned
#' the next integer after the maximum numeric chromosome found. Within each
#' chromosome, variants are positioned by their base-pair coordinate. A 2%
#' inter-chromosome gap is added relative to the total genomic span.
#'
#' Variants with NA or zero p-values are silently excluded. The y-axis is
#' \eqn{-\log_{10}(p)}. The x-axis tick marks are placed at chromosome
#' midpoints.
#'
#' Complexity: O(m log m) where m is the number of variants (dominated by
#' the sort to compute cumulative positions).
#'
#' @examples
#' \dontrun{
#' # Simulate results
#' sim <- data.frame(
#'   chr    = rep(1:22, each = 100),
#'   pos    = rep(seq(1e6, 1e8, length.out = 100), 22),
#'   pvalue = runif(2200)
#' )
#' class(sim) <- c("glow_marginal_scan", "data.frame")
#' plot_manhattan(sim)
#'
#' # With highlights
#' sim$variant_id <- seq_len(nrow(sim))
#' plot_manhattan(sim, highlight = c(1, 500, 1000))
#' }
#'
#' @seealso \code{\link{plot_qq}}, \code{\link{plot.glow_marginal_scan}},
#'   \code{\link{marginal_scan}}
#'
#' @export
plot_manhattan <- function(x, pvalue_col = "pvalue", chr_col = "chr",
                            pos_col = "pos", suggestive_line = 1e-5,
                            genome_wide_line = 5e-8, highlight = NULL,
                            chr_colors = c("grey30", "grey60"),
                            highlight_color = "red", title = "Manhattan Plot",
                            cex = 0.5, ...) {
  # Validate required columns
  missing_cols <- setdiff(c(pvalue_col, chr_col, pos_col), names(x))
  if (length(missing_cols) > 0) {
    stop("Columns not found in x: ", paste(missing_cols, collapse = ", "))
  }

  pvals <- x[[pvalue_col]]
  chr   <- x[[chr_col]]
  pos   <- x[[pos_col]]

  # Remove NA or non-positive p-values
  valid <- !is.na(pvals) & pvals > 0
  if (sum(valid) == 0) {
    stop("No valid p-values found (all NA or <= 0)")
  }
  pvals <- pvals[valid]
  chr   <- chr[valid]
  pos   <- pos[valid]

  # Preserve variant_id for highlighting if present
  has_vid <- "variant_id" %in% names(x)
  if (has_vid) {
    vid <- x[["variant_id"]][valid]
  }

  # Standardize chromosome names to integer
  chr_numeric <- .chr_to_numeric(chr)

  # Compute cumulative x-positions with 2% gap between chromosomes
  cum_pos <- .cumulative_positions(chr_numeric, pos)

  # -log10(p) transform
  log10p <- -log10(pvals)

  # Assign alternating colors per chromosome (by sorted unique chromosome order)
  chr_levels <- sort(unique(chr_numeric))
  n_chr <- length(chr_levels)
  chr_rank <- match(chr_numeric, chr_levels)  # 1-based rank for each point
  col_vec <- chr_colors[(chr_rank %% 2) + 1]

  # Override color for highlighted variants
  if (!is.null(highlight) && has_vid) {
    highlight_idx <- vid %in% highlight
    col_vec[highlight_idx] <- highlight_color
  } else if (!is.null(highlight) && !has_vid) {
    warning("'highlight' ignored: 'variant_id' column not found in x")
  }

  # Compute chromosome midpoints for x-axis labels
  mid_x <- vapply(chr_levels, function(ch) {
    idx <- chr_numeric == ch
    mean(range(cum_pos[idx]))
  }, numeric(1))

  # Main plot
  plot(cum_pos, log10p,
       pch    = 20,
       col    = col_vec,
       cex    = cex,
       xaxt   = "n",
       xlab   = "Chromosome",
       ylab   = expression(-log[10](italic(p))),
       main   = title,
       ...)

  # Custom x-axis: chromosome labels at midpoints
  chr_labels <- as.character(chr_levels)
  axis(1, at = mid_x, labels = chr_labels, las = 1, tick = FALSE, cex.axis = 0.75)

  # Validate significance thresholds (must be valid p-values if non-NULL)
  if (!is.null(suggestive_line) && (suggestive_line <= 0 || suggestive_line > 1)) {
    stop("'suggestive_line' must be a p-value in (0, 1], got: ", suggestive_line)
  }
  if (!is.null(genome_wide_line) && (genome_wide_line <= 0 || genome_wide_line > 1)) {
    stop("'genome_wide_line' must be a p-value in (0, 1], got: ", genome_wide_line)
  }

  # Significance lines
  if (!is.null(suggestive_line)) {
    abline(h = -log10(suggestive_line), col = "blue", lty = 2, lwd = 1)
  }
  if (!is.null(genome_wide_line)) {
    abline(h = -log10(genome_wide_line), col = "red", lty = 1, lwd = 1.2)
  }

  invisible(NULL)
}


#' QQ Plot for Marginal Analysis P-values
#'
#' Creates a quantile-quantile plot comparing observed p-values against the
#' expected uniform distribution, with an optional 95% confidence band and
#' genomic inflation factor (\eqn{\lambda_{GC}}).
#'
#' @param x A glow_marginal_scan data.frame, or a numeric vector of p-values.
#' @param pvalue_col Character name of the p-value column (used when x is a
#'   data.frame). Default: "pvalue".
#' @param ci Logical: draw a 95% pointwise confidence band? Default: TRUE.
#' @param title Character plot title. Default: "QQ Plot".
#' @param cex Numeric point size. Default: 0.5.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @return The genomic inflation factor \eqn{\lambda_{GC}} (invisibly).
#'
#' @details
#' The expected p-values under the null are
#' \eqn{e_i = (i - 0.5) / n}, for \eqn{i = 1, \ldots, n}, where n is the
#' number of valid p-values. Observed p-values are sorted in ascending order
#' before plotting.
#'
#' The genomic inflation factor is computed as:
#' \deqn{\lambda_{GC} = \frac{\mathrm{median}(\chi^2_1)}{\chi^2_{1, 0.5}}}
#' where \eqn{\chi^2_1 = \mathrm{qchisq}(1-p, 1)} for all valid p-values
#' (before sorting).
#'
#' The 95% confidence band is derived from the Beta distribution:
#' lower bound is \eqn{\mathrm{qbeta}(0.025, i, n-i+1)} and upper bound is
#' \eqn{\mathrm{qbeta}(0.975, i, n-i+1)}.
#'
#' Complexity: O(n log n) due to sorting.
#'
#' @examples
#' \dontrun{
#' # Numeric vector input
#' set.seed(42)
#' pvals <- runif(10000)
#' lambda <- plot_qq(pvals)
#' message("lambda_GC = ", round(lambda, 3))
#'
#' # data.frame input
#' df <- data.frame(pvalue = runif(500))
#' class(df) <- c("glow_marginal_scan", "data.frame")
#' lambda <- plot_qq(df)
#' }
#'
#' @seealso \code{\link{plot_manhattan}}, \code{\link{plot.glow_marginal_scan}},
#'   \code{\link{marginal_scan}}
#'
#' @export
plot_qq <- function(x, pvalue_col = "pvalue", ci = TRUE,
                     title = "QQ Plot", cex = 0.5, ...) {
  # Extract p-values from numeric vector or data.frame
  if (is.numeric(x)) {
    pvals <- x
  } else if (is.data.frame(x)) {
    if (!pvalue_col %in% names(x)) {
      stop("Column '", pvalue_col, "' not found in x")
    }
    pvals <- x[[pvalue_col]]
  } else {
    stop("x must be a numeric vector or a data.frame")
  }

  # Filter to valid p-values: no NA, positive, at most 1
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
  n <- length(pvals)
  if (n == 0) {
    stop("No valid p-values after filtering (all NA, <= 0, or > 1)")
  }

  # Compute lambda_GC from ALL valid p-values (unsorted)
  chisq_obs <- stats::qchisq(1 - pvals, df = 1)
  lambda_gc <- median(chisq_obs) / stats::qchisq(0.5, df = 1)

  # Sort observed and compute expected p-values
  obs_sorted <- sort(pvals)
  expected   <- (seq_len(n) - 0.5) / n

  # Transform to -log10 scale
  obs_log    <- -log10(obs_sorted)
  exp_log    <- -log10(expected)

  # Axis limits: use the range of expected + a little headroom
  xlim_val <- range(exp_log)
  ylim_val <- c(0, max(obs_log, xlim_val[2]))

  # Initial plot (points may be redrawn on top of CI band below)
  plot(exp_log, obs_log,
       pch  = 20,
       cex  = cex,
       xlab = expression(Expected ~ -log[10](italic(p))),
       ylab = expression(Observed ~ -log[10](italic(p))),
       main = title,
       xlim = xlim_val,
       ylim = ylim_val,
       type = if (ci) "n" else "p",   # suppress points initially when drawing CI
       ...)

  # 95% confidence band using the Beta distribution
  if (ci) {
    i_seq <- seq_len(n)
    ci_lo <- -log10(stats::qbeta(0.975, i_seq, n - i_seq + 1))
    ci_hi <- -log10(stats::qbeta(0.025, i_seq, n - i_seq + 1))

    # Draw polygon (band in x = expected, y = ci_hi/ci_lo)
    graphics::polygon(
      x   = c(exp_log, rev(exp_log)),
      y   = c(ci_hi,   rev(ci_lo)),
      col = grDevices::rgb(0, 0, 1, 0.15),
      border = NA
    )

    # Re-draw points on top of band
    points(exp_log, obs_log, pch = 20, cex = cex)
  }

  # Identity line
  abline(0, 1, col = "red", lwd = 1.5)

  # Lambda annotation
  legend("topleft",
         legend = bquote(lambda[GC] == .(round(lambda_gc, 3))),
         bty = "n",
         cex = 0.9)

  invisible(lambda_gc)
}


#' Plot Method for Marginal Scan Results
#'
#' S3 \code{plot} method for objects of class \code{"glow_marginal_scan"}
#' (returned by \code{\link{marginal_scan}}). Dispatches to
#' \code{\link{plot_manhattan}}, \code{\link{plot_qq}}, or both.
#'
#' @param x A \code{glow_marginal_scan} object.
#' @param type Character: \code{"manhattan"}, \code{"qq"}, or \code{"both"}.
#'   Default: \code{"manhattan"}.
#' @param ... Additional arguments forwarded to the underlying plot function.
#'
#' @return For \code{type = "qq"} or \code{type = "both"}: the genomic
#'   inflation factor \eqn{\lambda_{GC}} (invisibly).
#'   For \code{type = "manhattan"}: \code{invisible(NULL)}.
#'
#' @details
#' When \code{type = "both"}, the graphics device is split into a 1x2 panel
#' using \code{par(mfrow = c(1, 2))}. The original \code{par} settings are
#' restored on exit via \code{on.exit()}.
#'
#' All further arguments in \code{...} are passed to \emph{both} underlying
#' plot functions when \code{type = "both"}.
#'
#' @examples
#' \dontrun{
#' results <- marginal_scan("chr22.gds", nm, use_SPA = FALSE)
#' plot(results)                        # Manhattan
#' plot(results, type = "qq")           # QQ
#' plot(results, type = "both")         # side by side
#' }
#'
#' @seealso \code{\link{plot_manhattan}}, \code{\link{plot_qq}},
#'   \code{\link{marginal_scan}}
#'
#' @export
plot.glow_marginal_scan <- function(x, type = "manhattan", ...) {
  type <- match.arg(type, c("manhattan", "qq", "both"))

  if (type == "manhattan") {
    plot_manhattan(x, ...)
    return(invisible(NULL))
  }

  if (type == "qq") {
    return(invisible(plot_qq(x, ...)))
  }

  # type == "both": split device into 1x2 panel
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2))

  plot_manhattan(x, ...)
  lambda <- plot_qq(x, ...)

  invisible(lambda)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Convert Chromosome Labels to Integer
#'
#' Strips "chr" prefix and converts to integer. Non-numeric chromosomes
#' (X, Y, MT) are assigned consecutive integers after the max numeric value.
#'
#' @param chr Character or integer vector of chromosome labels.
#' @return Integer vector of the same length.
#' @keywords internal
#' @noRd
.chr_to_numeric <- function(chr) {
  chr_str <- as.character(chr)
  # Strip leading "chr" (case-insensitive)
  chr_stripped <- gsub("^[Cc][Hh][Rr]", "", chr_str)

  chr_int <- suppressWarnings(as.integer(chr_stripped))
  na_idx  <- is.na(chr_int)

  if (any(na_idx)) {
    # Assign non-numeric chromosomes (X, Y, MT) integers beyond the max
    max_num <- if (any(!na_idx)) max(chr_int[!na_idx], na.rm = TRUE) else 0L
    non_num_labels <- unique(chr_stripped[na_idx])
    mapping <- setNames(
      seq_along(non_num_labels) + max_num,
      non_num_labels
    )
    chr_int[na_idx] <- mapping[chr_stripped[na_idx]]
  }

  chr_int
}


#' Compute Cumulative Genomic Positions with Inter-Chromosome Gaps
#'
#' Converts per-chromosome base-pair positions to a continuous x-axis by
#' appending chromosomes sequentially. A gap of 2% of the total raw span
#' (excluding gaps) is inserted between chromosomes.
#'
#' @param chr_numeric Integer vector of chromosome labels (from .chr_to_numeric).
#' @param pos Numeric vector of base-pair positions.
#' @return Numeric vector of cumulative x-positions, same length as \code{pos}.
#' @keywords internal
#' @noRd
.cumulative_positions <- function(chr_numeric, pos) {
  chr_levels <- sort(unique(chr_numeric))
  n_chr      <- length(chr_levels)
  n_pts      <- length(pos)

  # Per-chromosome max positions
  chr_max <- vapply(chr_levels, function(ch) {
    max(pos[chr_numeric == ch], na.rm = TRUE)
  }, numeric(1))

  # Gap = 2% of total raw span
  total_span <- sum(chr_max)
  gap        <- 0.02 * total_span

  # Compute offset for each chromosome (cumulative end of previous + gap)
  offsets <- numeric(n_chr)
  for (k in seq_len(n_chr)[-1]) {
    offsets[k] <- offsets[k - 1L] + chr_max[k - 1L] + gap
  }

  cum_pos <- numeric(n_pts)
  for (k in seq_len(n_chr)) {
    idx <- chr_numeric == chr_levels[k]
    cum_pos[idx] <- pos[idx] + offsets[k]
  }

  cum_pos
}
