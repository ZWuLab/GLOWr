# ==============================================================================
# Highlight-Region Selection for Manhattan Plots
# ==============================================================================
#
# Region-agnostic helper that picks the regions (genes, windows, loci) to
# highlight on a Manhattan plot: every region at or below a p-value threshold,
# topped up to a minimum count so the plot is never empty. Used for both genes
# and sliding-window loci, hence the region-agnostic `label_col`.
#
# EXPORTED FUNCTIONS:
#   - select_highlight_regions()   Regions to highlight (threshold + top-up floor)


#################### EXPORTED MAIN FUNCTIONS ####################

#' Select Regions to Highlight on a Manhattan Plot
#'
#' Picks the regions (genes, windows, loci) to highlight on a Manhattan plot.
#' Returns every region whose p-value is at or below \code{threshold}; if fewer
#' than \code{min_n} regions pass, the selection is topped up with the
#' next-smallest p-values so the plot is never empty (useful when no region
#' reaches genome-wide significance).
#'
#' @details
#' The returned labels are taken in ascending p-value order. The selection size
#' is \code{min(max(n_pass, min_n), n_nonNA)}, where \code{n_pass} is the count
#' at or below \code{threshold} and \code{n_nonNA} is the number of non-\code{NA}
#' p-values; the result is empty only when there are no usable p-values.
#' Complexity is \eqn{O(n \log n)} from the p-value sort.
#'
#' @param df data.frame of region-level results.
#' @param pvalue_col Character. Name of the p-value column (e.g.,
#'   \code{"GLOW_Omni"}).
#' @param threshold Numeric. P-value cutoff; typically \code{0.05 / nrow(df)}.
#' @param label_col Character. Name of the region-label column. Default
#'   \code{"label"}.
#' @param min_n Integer. Minimum number of regions to return when fewer pass
#'   \code{threshold}. Default \code{10L}.
#'
#' @return Character vector of region labels sorted by ascending p-value, of
#'   length \code{max(n_pass, min_n)} (capped at the number of non-\code{NA}
#'   p-values).
#'
#' @examples
#' df <- data.frame(label = paste0("g", 1:5),
#'                  GLOW_Omni = c(1e-9, 0.2, 0.4, 0.6, 0.8))
#' select_highlight_regions(df, "GLOW_Omni", threshold = 1e-6, min_n = 3L)
#'
#' @export
select_highlight_regions <- function(df, pvalue_col, threshold,
                                     label_col = "label", min_n = 10L) {
  p <- df[[pvalue_col]]
  ok <- !is.na(p)
  n_pass <- sum(p[ok] <= threshold)
  n_pick <- min(max(n_pass, as.integer(min_n)), sum(ok))
  if (n_pick == 0L) return(character(0))
  ord <- order(p)
  df[[label_col]][ord][seq_len(n_pick)]
}
