# ==============================================================================
# Genomic-Control p-value Calibration
# ==============================================================================
#
# Scalar genomic-control recalibration of p-values: rescale each 1-df
# chi-square by a single calibration factor, diminishing a potentially
# inflated association test at genome-wide scale. 
# The factor is GENERAL -- `method` only selects its source:
#   - "lambda_gc"      : the genomic-control lambda, estimated from `p` itself.
#   - "ldsc_intercept" : the LD Score regression intercept (supplied), which is designed to
#                        removes confounding while preserving polygenic signal
#                        (see ldsc_regression()).
# The current rescaling itself is identical across methods; only the source of the
# factor differs. This separation (a general `calibration_factor` value + a
# `method` that names where it came from) mirrors how a single scalar GC
# correction is applied throughout GLOW (e.g. the M2 z-scaling z_scale =
# 1/sqrt(factor) in glow_region_test()).
#
# EXPORTED FUNCTIONS:
#   - calibrate_pvalues()   Rescale p-values by a calibration factor
#   - compute_lambda_gc()   Genomic-control inflation factor lambda_GC


#' Calibrate p-values by a Scalar Genomic-Control Factor
#'
#' Recalibrates a vector of p-values by rescaling each on the 1-df chi-square
#' scale by a single \code{calibration_factor}, mapping a test whose inflation
#' factor is \eqn{f} toward \eqn{\lambda = 1}. The same rescaling underlies both
#' classical genomic control and the LD Score regression intercept correction;
#' \code{method} selects only where \eqn{f} comes from.
#'
#' @details
#' On the 1-df chi-square scale the calibration is
#' \deqn{\chi^2_{\mathrm{cal}} = \chi^2 / f, \qquad
#'       p_{\mathrm{cal}} = \Pr\!\left(\chi^2_1 > \chi^2_{\mathrm{obs}} / f\right),}
#' computed in the upper tail for accuracy at small \eqn{p}. The rescaling is
#' identical across methods; \code{method} only selects the source of the factor
#' \eqn{f}:
#' \describe{
#'   \item{\code{"lambda_gc"}}{When \code{calibration_factor} is \code{NULL} it
#'     defaults to \code{compute_lambda_gc(p)} (the genomic-control \eqn{\lambda}
#'     estimated from \code{p} itself). Here \code{calibration_factor} \emph{is}
#'     \eqn{\lambda_{\mathrm{GC}}}, which removes \emph{all} inflation
#'     (confounding and polygenic).}
#'   \item{\code{"ldsc_intercept"}}{\code{calibration_factor} is \strong{required}
#'     (the intercept from \code{\link{ldsc_regression}}); it is not estimable
#'     from the p-values being calibrated. The intercept removes confounding
#'     while \emph{preserving} polygenic signal.}
#' }
#' A non-finite or non-positive \code{calibration_factor} passes \code{p} through
#' unchanged (e.g. \code{"lambda_gc"} on fewer than 50 usable p-values, where
#' \code{compute_lambda_gc} returns \code{NA}). \code{NA} / out-of-\eqn{(0,1]}
#' entries flow through as \code{NA}.
#'
#' @param p Numeric vector of p-values.
#' @param method Calibration method -- the source of the factor:
#'   \code{"lambda_gc"} (default) or \code{"ldsc_intercept"}.
#' @param calibration_factor Optional numeric scalar, the chi-square divisor
#'   \eqn{f}. For \code{"lambda_gc"} it defaults to \code{compute_lambda_gc(p)};
#'   for \code{"ldsc_intercept"} it is required (the supplied LDSC intercept).
#'
#' @return A list with components \code{p} (the calibrated vector, same length
#'   and order as the input), \code{calibration_factor} (the factor used), and
#'   \code{method}.
#'
#' @seealso \code{\link{compute_lambda_gc}}, \code{\link{ldsc_regression}}
#'
#' @references
#' Bulik-Sullivan B. et al. (2015). LD Score regression distinguishes
#' confounding from polygenicity in GWAS. \emph{Nat Genet} 47:291.
#'
#' @examples
#' set.seed(1)
#' chi2 <- rchisq(2000, df = 1) * 1.3                  # an inflated null
#' p    <- pchisq(chi2, df = 1, lower.tail = FALSE)
#' compute_lambda_gc(p)                                 # ~1.3
#' cal <- calibrate_pvalues(p, method = "lambda_gc")    # factor estimated from p
#' compute_lambda_gc(cal$p)                             # ~1.0
#' # LDSC-intercept correction with an externally estimated factor:
#' calibrate_pvalues(p, method = "ldsc_intercept", calibration_factor = 1.08)$p[1:3]
#'
#' @export
calibrate_pvalues <- function(p,
                              method = c("lambda_gc", "ldsc_intercept"),
                              calibration_factor = NULL) {
  method <- match.arg(method)

  if (is.null(calibration_factor)) {
    if (method == "lambda_gc") {
      calibration_factor <- compute_lambda_gc(p)
    } else {
      stop("method = \"ldsc_intercept\" requires `calibration_factor` (e.g. ",
           "ldsc_regression()$intercept); it cannot be estimated from `p`.",
           call. = FALSE)
    }
  }

  pc <- rep(NA_real_, length(p))
  # Undefined / non-positive factor (e.g. lambda_gc on < 50 p-values, or a
  # supplied NA): pass `p` through unchanged.
  if (length(calibration_factor) != 1L || is.na(calibration_factor) ||
      calibration_factor <= 0) {
    return(list(p = p, calibration_factor = calibration_factor, method = method))
  }

  ok <- is.finite(p) & p >= 0 & p <= 1
  # Work in the upper tail throughout for accuracy at small p.
  x <- stats::qchisq(p[ok], df = 1L, lower.tail = FALSE)   # p -> chi-square_1 (upper)
  pc[ok] <- stats::pchisq(x / calibration_factor, df = 1L, lower.tail = FALSE)
  list(p = pc, calibration_factor = calibration_factor, method = method)
}


#' Genomic-Control Inflation Factor (lambda_GC)
#'
#' The genomic-control inflation factor
#' \eqn{\lambda_{\mathrm{GC}} = \mathrm{median}(\chi^2) / \mathrm{qchisq}(0.5, 1)},
#' with the chi-squares recovered from the p-values as \eqn{\mathrm{qchisq}(1-p, 1)}.
#'
#' @details
#' Only finite p-values in \eqn{(0, 1]} are used. Returns \code{NA} when fewer
#' than 50 such values are available (too few to interpret).
#'
#' @param pvalues Numeric vector of p-values.
#'
#' @return The inflation factor (numeric scalar), or \code{NA_real_} when there
#'   are fewer than 50 usable p-values.
#'
#' @seealso \code{\link{calibrate_pvalues}}
#'
#' @examples
#' set.seed(1)
#' compute_lambda_gc(pchisq(rchisq(2000, 1) * 1.2, 1, lower.tail = FALSE))  # ~1.2
#'
#' @export
compute_lambda_gc <- function(pvalues) {
  p <- pvalues[is.finite(pvalues) & pvalues > 0 & pvalues <= 1]
  if (length(p) < 50L) return(NA_real_)
  stats::median(stats::qchisq(1 - p, df = 1L)) / stats::qchisq(0.5, df = 1L)
}
