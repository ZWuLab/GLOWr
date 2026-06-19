# ==============================================================================
# LD Score Regression (Confounding vs Polygenicity)
# ==============================================================================
#
# Regresses single-variant chi-square statistics on their LD scores to split
# genomic-control inflation into a CONFOUNDING component (the intercept) and a
# POLYGENIC component (the slope). The INTERCEPT is the principled
# confounding-only inflation factor used to calibrate GLOW scans (see
# calibrate_pvalues(method = "ldsc_intercept") and the M2 z-scaling).
#
# EXPORTED FUNCTIONS:
#   - ldsc_regression()   Weighted two-pass LDSC regression + chr jackknife SEs
#
# INTERNAL HELPERS:
#   - .ldsc_fit()             One weighted two-pass fit -> c(intercept, slope)
#   - .ldsc_jackknife_se()    Chromosome block-jackknife SE of a statistic
#
# This function performs the regression itself; preparing and QC-filtering the
# per-variant chi-square statistics and LD scores it takes as input is the
# caller's responsibility.


#' LD Score Regression: Confounding vs Polygenicity
#'
#' Fits the LD Score regression of single-variant chi-square statistics on their
#' LD scores, returning the \strong{intercept} (the confounding-only inflation
#' factor), the slope (polygenicity), and the confounding ratio, each with a
#' chromosome block-jackknife standard error. The intercept is the principled
#' calibration factor for \code{\link{calibrate_pvalues}} (\code{method =
#' "ldsc_intercept"}) and the GLOW M2 Z-scaling.
#'
#' @details
#' The model (Bulik-Sullivan et al. 2015) is
#' \deqn{E[\chi^2_j] = 1 + N a + \frac{N h^2}{M}\,\ell_j,}
#' where \eqn{\ell_j} is the LD score. The fitted \strong{intercept} estimates
#' \eqn{1 + Na} (here \eqn{a} is the model's confounding parameter -- \emph{not}
#' the intercept itself); \eqn{(\mathrm{intercept} - 1)} is the confounding
#' contribution and the \strong{slope} \eqn{\propto h^2/M} is polygenicity. The
#' \strong{confounding ratio}
#' \deqn{\mathrm{ratio} = \frac{\mathrm{intercept} - 1}{\overline{\chi^2} - 1}}
#' is the fraction of the inflation due to confounding rather than polygenicity
#' -- a more robust read-out than \eqn{\lambda_{\mathrm{GC}}}, which itself grows
#' with \eqn{N} under pure polygenicity.
#'
#' The fit is weighted two-pass least squares: pass 1 weights by
#' \eqn{1/\max(\ell, 1)} (LD double-counting); pass 2 multiplies in a
#' \eqn{(1 + \max(\widehat{\mathrm{slope}}, 0)\,\ell)^{-2}} heteroskedasticity
#' weight. Standard errors (intercept, slope, ratio) come from a
#' \strong{chromosome block jackknife} over \code{unique(chr)}. \eqn{\chi^2} is
#' winsorised at \code{winsorize} inside the fit for robustness; the reported
#' \code{lambda_gc} uses the un-winsorised \eqn{\chi^2}.
#'
#' \strong{QC is the caller's responsibility.} Pass already-QC'd, matched
#' \code{(chi2, ld, chr)} (standard LDSC QC: common variants, e.g.
#' \eqn{\mathrm{MAF} \ge 0.01}; exclude the MHC; keys matched to the LD-score
#' table). Non-finite \code{(chi2, ld)} pairs are dropped.
#'
#' @param chi2 Numeric vector of single-variant chi-square statistics
#'   (\eqn{Z^2}).
#' @param ld Numeric vector of LD scores (from \code{\link{compute_ld_scores}}),
#'   aligned to \code{chi2}.
#' @param chr Vector (character/integer) of chromosomes, aligned to \code{chi2};
#'   defines the jackknife blocks.
#' @param winsorize Numeric cap applied to \eqn{\chi^2} inside the fit (default
#'   \code{80}).
#' @param two_pass Logical; if \code{TRUE} (default) use the heteroskedasticity-
#'   weighted second pass, else stop after pass 1.
#'
#' @return A one-row \code{data.frame}: \code{n_variants}, \code{lambda_gc},
#'   \code{mean_chi2}, \code{intercept}, \code{intercept_se}, \code{slope},
#'   \code{slope_se}, \code{confounding_ratio}, \code{ratio_se}. The
#'   \code{intercept} is the confounding-only inflation factor.
#'
#' @seealso \code{\link{compute_ld_scores}}, \code{\link{calibrate_pvalues}}
#'
#' @references
#' Bulik-Sullivan B. et al. (2015). LD Score regression distinguishes
#' confounding from polygenicity in GWAS. \emph{Nat Genet} 47:291.
#'
#' @examples
#' set.seed(1)
#' ld  <- stats::rgamma(5000, shape = 2, scale = 30)            # LD scores
#' chr <- rep(1:22, length.out = 5000)
#' chi2 <- stats::rchisq(5000, df = 1, ncp = 0) + 0.08 + 5e-4 * ld  # intercept ~1.08
#' fit <- ldsc_regression(chi2, ld, chr)
#' fit$intercept            # ~1.08 (the confounding-only factor)
#'
#' @export
ldsc_regression <- function(chi2, ld, chr, winsorize = 80, two_pass = TRUE) {
  stopifnot(length(chi2) == length(ld), length(ld) == length(chr))
  ok <- is.finite(chi2) & is.finite(ld)
  chi2 <- chi2[ok]; ld <- ld[ok]; chr <- as.character(chr)[ok]
  n <- length(chi2)
  if (n < 100L) {
    stop("ldsc_regression(): need at least 100 finite (chi2, ld) pairs after ",
         "QC; got ", n, ".", call. = FALSE)
  }

  chi2w   <- pmin(chi2, winsorize)
  fit     <- .ldsc_fit(chi2w, ld, two_pass)
  meanchi <- mean(chi2w)
  ratio   <- (fit[["intercept"]] - 1) / (meanchi - 1)

  # Chromosome block jackknife: one set of leave-one-block-out fits yields the
  # SEs of intercept, slope, and ratio together.
  blocks <- unique(chr)
  K <- length(blocks)
  if (K < 2L) {
    stop("ldsc_regression(): need >= 2 distinct chr blocks for the jackknife; ",
         "got ", K, ".", call. = FALSE)
  }
  jk <- vapply(blocks, function(b) {
    keep <- chr != b
    e <- .ldsc_fit(chi2w[keep], ld[keep], two_pass)
    c(intercept = e[["intercept"]], slope = e[["slope"]],
      ratio = (e[["intercept"]] - 1) / (mean(chi2w[keep]) - 1))
  }, c(intercept = 0, slope = 0, ratio = 0))
  jk_se <- function(x) sqrt((K - 1) / K * sum((x - mean(x))^2))

  data.frame(
    n_variants        = n,
    lambda_gc         = stats::median(chi2) / stats::qchisq(0.5, df = 1L),
    mean_chi2         = meanchi,
    intercept         = fit[["intercept"]],
    intercept_se      = jk_se(jk["intercept", ]),
    slope             = fit[["slope"]],
    slope_se          = jk_se(jk["slope", ]),
    confounding_ratio = ratio,
    ratio_se          = jk_se(jk["ratio", ]),
    row.names = NULL, stringsAsFactors = FALSE)
}


#################### INTERNAL HELPERS ####################

#' One weighted (two-pass) LDSC fit
#'
#' @param chi2w Winsorised chi-square vector.
#' @param ld LD-score vector (aligned to \code{chi2w}).
#' @param two_pass Logical; add the heteroskedasticity-weighted second pass.
#' @return Named numeric \code{c(intercept, slope)}.
#' @keywords internal
#' @noRd
.ldsc_fit <- function(chi2w, ld, two_pass = TRUE) {
  reg <- function(w) {
    cf <- as.numeric(stats::coef(stats::lm(chi2w ~ ld, weights = w)))
    c(intercept = cf[1L], slope = cf[2L])
  }
  e0 <- reg(1 / pmax(ld, 1))                        # pass 1: LD double-counting
  if (!isTRUE(two_pass)) return(e0)
  het <- (1 + pmax(e0[["slope"]], 0) * ld)^2        # heteroskedasticity weight
  reg(1 / (pmax(ld, 1) * het))                      # pass 2
}
