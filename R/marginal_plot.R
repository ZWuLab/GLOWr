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
#   - .thin_cells()              Thinning grid resolution from the marker size
#   - .thin_points()             One drawn point per grid cell, the highest kept
#   - .qq_band_ranks()           Ranks at which the QQ confidence band is evaluated


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
#' @param thin Logical (default TRUE). Draw at most one point per cell of a grid at the
#'   marker's scale (see \code{thin_cells_per_inch}), which bounds the drawn points at a
#'   few tens of thousands however many variants there are, so a whole-genome PDF stays
#'   small. Every highlighted variant is drawn, and the point kept in a cell is its highest,
#'   so the extremes are always visible. The axes and colors come from the complete input,
#'   and a point outside a cropped axis range is kept and clipped as usual. At the marker's
#'   size the drawing matches the full set; a vector PDF enlarged well beyond that scale
#'   shows the thinning (\code{thin = FALSE} draws every point).
#' @param thin_cells_per_inch Numeric or NULL (default). The thinning grid's cells per inch
#'   of the plot region. NULL sets the cell to half the marker's drawn diameter
#'   (\code{pch = 20} at the given \code{cex}: about 0.035 inch at \code{cex = 0.5},
#'   measured at 1200 dpi), so two points in one cell overlap by most of their area.
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
                            cex = 0.5, thin = TRUE, thin_cells_per_inch = NULL, ...) {
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

  # Highlighted variants: their own color, drawn last and never thinned away
  is_high <- rep(FALSE, length(cum_pos))
  if (!is.null(highlight) && has_vid) {
    is_high <- vid %in% highlight
    col_vec[is_high] <- highlight_color
  } else if (!is.null(highlight) && !has_vid) {
    warning("'highlight' ignored: 'variant_id' column not found in x")
  }

  # Compute chromosome midpoints for x-axis labels
  mid_x <- vapply(chr_levels, function(ch) {
    idx <- chr_numeric == ch
    mean(range(cum_pos[idx]))
  }, numeric(1))

  # Main plot: the frame from the complete input, then the points, at most one per grid
  # cell when `thin` (the cell is sized from the marker once the frame fixes the scale)
  plot(range(cum_pos), range(log10p),
       type   = "n",
       xaxt   = "n",
       xlab   = "Chromosome",
       ylab   = expression(-log[10](italic(p))),
       main   = title,
       ...)
  keep <- if (isTRUE(thin)) {
    .thin_points(cum_pos, log10p, .thin_cells(thin_cells_per_inch, cex))
  } else {
    rep(TRUE, length(cum_pos))
  }
  draw <- keep & !is_high
  points(cum_pos[draw], log10p[draw], pch = 20, col = col_vec[draw], cex = cex)
  if (any(is_high)) {
    points(cum_pos[is_high], log10p[is_high], pch = 20, col = highlight_color, cex = cex)
  }

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
#' @param thin Logical (default TRUE). Draw at most one point per cell of a grid at the
#'   marker's scale (see \code{thin_cells_per_inch}); the lambda, the expected quantiles
#'   and the confidence band are computed from the complete input.
#' @param thin_cells_per_inch Numeric or NULL (default). As in \code{\link{plot_manhattan}}.
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
#' \eqn{\mathrm{qbeta}(0.975, i, n-i+1)}, with the original \eqn{n}. It is evaluated at
#' every rank among the 200 smallest p-values and at about 2,000 geometrically spaced
#' ranks beyond, which draws the same smooth band with a few thousand vertices instead
#' of \eqn{2n}.
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
                     title = "QQ Plot", cex = 0.5,
                     thin = TRUE, thin_cells_per_inch = NULL, ...) {
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

  # The frame from the complete input; the band and the points are drawn on it below
  plot(xlim_val, ylim_val,
       type = "n",
       xlab = expression(Expected ~ -log[10](italic(p))),
       ylab = expression(Observed ~ -log[10](italic(p))),
       main = title,
       xlim = xlim_val,
       ylim = ylim_val,
       ...)

  # 95% confidence band from the Beta distribution of the i-th order statistic, with the
  # original n, at a bounded selection of ranks: the same smooth band with a few thousand
  # vertices instead of 2n
  if (ci) {
    i_sel <- .qq_band_ranks(n)
    ci_lo <- -log10(stats::qbeta(0.975, i_sel, n - i_sel + 1))
    ci_hi <- -log10(stats::qbeta(0.025, i_sel, n - i_sel + 1))
    graphics::polygon(
      x   = c(exp_log[i_sel], rev(exp_log[i_sel])),
      y   = c(ci_hi,          rev(ci_lo)),
      col = grDevices::rgb(0, 0, 1, 0.15),
      border = NA
    )
  }

  # Points on top of the band, at most one per grid cell when `thin`
  keep <- if (isTRUE(thin)) {
    .thin_points(exp_log, obs_log, .thin_cells(thin_cells_per_inch, cex))
  } else {
    rep(TRUE, n)
  }
  points(exp_log[keep], obs_log[keep], pch = 20, cex = cex)

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

#' Thinning Grid Resolution From the Marker Size
#'
#' Cells per inch of the plot region: the caller's value when given, else half the drawn
#' diameter of \code{pch = 20} at \code{cex}. The diameter grows about linearly with
#' \code{cex}: 0.035 inch at 0.5 and 0.060 inch at 1 (measured on a 1200 dpi raster), so
#' it is taken as \eqn{0.010 + 0.050\,\mathrm{cex}} inch.
#'
#' @param cells_per_inch Numeric or NULL.
#' @param cex Numeric marker expansion.
#' @return A positive number, cells per inch.
#' @keywords internal
#' @noRd
.thin_cells <- function(cells_per_inch, cex) {
  if (!is.null(cells_per_inch)) {
    stopifnot(is.numeric(cells_per_inch), length(cells_per_inch) == 1L,
              is.finite(cells_per_inch), cells_per_inch > 0)
    return(cells_per_inch)
  }
  diam_in <- 0.010 + 0.050 * cex
  2 / diam_in
}


#' One Drawn Point per Grid Cell, the Highest Kept
#'
#' Maps the points to the cells of a grid over the plot region of the current frame
#' (\code{par("usr")} and \code{par("pin")}) and keeps, per cell, the point with the
#' largest y. Two points in one cell overlap at the marker's size, so at that scale the
#' picture is the full set's; a vector PDF enlarged well beyond it shows the thinning. A
#' logarithmic axis is mapped as \code{par("usr")} is, and a point outside the plot region
#' is kept without taking a cell, so a cropped axis hides nothing the full set would show.
#'
#' @param x,y Numeric user coordinates (untransformed; a log axis is handled here).
#' @param cells_per_inch Numeric grid resolution.
#' @return Logical vector, TRUE for the points to draw.
#' @keywords internal
#' @noRd
.thin_points <- function(x, y, cells_per_inch) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  # par("usr") is on the log10 scale for a logarithmic axis; map the coordinates the same way
  if (isTRUE(graphics::par("xlog"))) x <- log10(x)
  if (isTRUE(graphics::par("ylog"))) y <- log10(y)
  nx <- max(1, ceiling(pin[1L] * cells_per_inch))
  ny <- max(1, ceiling(pin[2L] * cells_per_inch))
  # A point outside the plot region (a cropped axis) takes no cell and is kept: the device
  # clips it as it would without thinning, and it cannot hide a visible point in its cell.
  inside <- is.finite(x) & is.finite(y) &
            x >= usr[1L] & x <= usr[2L] & y >= usr[3L] & y <= usr[4L]
  keep <- !inside
  ix <- floor((x[inside] - usr[1L]) / (usr[2L] - usr[1L]) * nx)
  iy <- floor((y[inside] - usr[3L]) / (usr[4L] - usr[3L]) * ny)
  ix <- pmin(ix, nx - 1)                   # the upper edge belongs to the last cell
  iy <- pmin(iy, ny - 1)
  cell <- ix * ny + iy                     # one number per (ix, iy); exact in double
  yi <- y[inside]
  o <- order(yi, decreasing = TRUE)        # the highest point of a cell is the one kept
  k <- logical(length(yi))
  k[o[!duplicated(cell[o])]] <- TRUE
  keep[inside] <- k
  keep
}


#' Ranks at Which the QQ Confidence Band Is Evaluated
#'
#' Every rank among the \code{dense} smallest p-values, then \code{geometric} ranks spaced
#' geometrically up to \code{n}, and \code{n} itself; all ranks when \code{n} is small.
#'
#' @param n Integer number of p-values.
#' @param dense,geometric Integers.
#' @return Sorted integer vector of distinct ranks in \code{1:n}.
#' @keywords internal
#' @noRd
.qq_band_ranks <- function(n, dense = 200L, geometric = 2000L) {
  if (n <= dense + 1L) return(seq_len(n))
  geo <- round(exp(seq(log(dense), log(n), length.out = geometric)))
  sort(unique(c(seq_len(dense), as.integer(geo), as.integer(n))))
}


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
