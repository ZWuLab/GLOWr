# ==============================================================================
# File: glow_defaults.R
# Purpose: The one place that names GLOW's variant-preprocessing defaults (the LD-pruning
#   threshold, the MAC collapsing threshold and the B/PI aggregation method), so that GLOWr's
#   functions, GLOWpipeline's scan engine and the GLOWanalyses base configuration read the
#   same values and a Terra or cluster run cannot diverge from a direct call.
#
# EXPORTED FUNCTIONS:
#   - glow_preprocessing_defaults()   The three defaults as a named list

#################### EXPORTED MAIN FUNCTIONS ####################

#' GLOW's Variant-Preprocessing Defaults
#'
#' Returns the default LD-pruning threshold, rare-variant collapsing threshold and B/PI
#' aggregation method that \code{\link{prepare_glow_input}}, \code{\link{glow_region_test}},
#' \code{\link{filter_variants_ld}} and \code{\link{collapse_rare_variants}} use when the
#' caller does not set them. GLOWpipeline's scan engine and the GLOWanalyses base
#' configuration read the same values, so every route applies one convention unless a run
#' overrides it explicitly.
#'
#' @details
#' \describe{
#'   \item{\code{ld_threshold = 0.9}}{LD pruning removes one variant of every pair whose
#'     absolute genotype correlation \eqn{|r|} exceeds 0.9 (the threshold is on \eqn{|r|},
#'     not on \eqn{r^2}).}
#'   \item{\code{mac_threshold = 10L}}{Variants with cohort minor allele count strictly below
#'     10 are collapsed into one unit; 0 disables collapsing.}
#'   \item{\code{collapse_method = "mean"}}{The B and PI weights of a collapsed group are
#'     aggregated by their mean. The genotypes of the group are always summed.}
#' }
#' These are the settings of the methods paper (Zhang, Liu, Landers and Wu, 2026): its
#' simulations excluded variants with MAC below 10 and pairs with \eqn{|r|} above 0.9, and
#' its UK Biobank analysis collapsed variants with MAC below 10.
#'
#' @return A named list with elements \code{ld_threshold} (numeric), \code{mac_threshold}
#'   (integer) and \code{collapse_method} (character).
#'
#' @examples
#' glow_preprocessing_defaults()
#' glow_preprocessing_defaults()$ld_threshold
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. (2026). Integrated weighted
#' association test with application to genetic association studies. The Annals
#' of Applied Statistics, 20(3), 2074-2097. \doi{10.1214/26-AOAS2223}
#'
#' @seealso \code{\link{prepare_glow_input}}, \code{\link{filter_variants_ld}},
#'   \code{\link{collapse_rare_variants}}
#' @export
glow_preprocessing_defaults <- function() {
  list(ld_threshold = 0.9, mac_threshold = 10L, collapse_method = "mean")
}
