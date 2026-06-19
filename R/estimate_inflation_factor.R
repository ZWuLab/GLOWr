# ==============================================================================
# Estimate the Single-Variant Inflation Factor from a Scan
# ==============================================================================
#
# One reusable "scan -> calibration factor" tool: given a single-variant scan,
# estimate the genomic-control inflation factor by the chosen method, applying
# the standard LDSC QC (common variants; exclude the MHC). The factor is the
# input to both calibration methods:
#   - M1: calibrate_pvalues(method = "ldsc_intercept", calibration_factor = factor)
#   - M2: z_scale = 1 / sqrt(factor) into glow_region_test()/glow_scan_chunk()
# so the single-variant calibration step, the set-scan run config (Stage 1), and
# the diagnostic study all obtain the factor from one tested function.
#
# EXPORTED FUNCTIONS:
#   - estimate_inflation_factor()   Scan -> inflation factor (ldsc_intercept | lambda_gc)


#' Estimate the Single-Variant Inflation Factor from a Scan
#'
#' Given a single-variant association scan, estimate its genomic-control
#' inflation factor by one of two methods, after the standard LDSC QC
#' (common variants; the MHC excluded). The returned factor is the input to the
#' GLOW calibration methods: \code{\link{calibrate_pvalues}} (M1) and the
#' \code{z_scale = 1/sqrt(factor)} hook in \code{\link{glow_region_test}} (M2).
#'
#' @details
#' \describe{
#'   \item{\code{"ldsc_intercept"} (the principled, confounding-only factor)}{The
#'     scan is merged to \code{ld_scores} by \code{chr:pos:ref:alt}, QC'd
#'     (\eqn{\mathrm{MAF} \ge}\code{maf_min}; the MHC dropped; positive finite LD
#'     scores), and \code{\link{ldsc_regression}} is run on
#'     \eqn{\chi^2 = Z^2} vs the LD score; the \strong{intercept} is the factor.}
#'   \item{\code{"lambda_gc"} (removes all inflation)}{The genomic-control
#'     \eqn{\lambda_{\mathrm{GC}}} via \code{\link{compute_lambda_gc}} on the
#'     scan's p-values (after the same MAF QC). No LD scores needed.}
#' }
#' Applying the same MAF QC to both methods keeps them comparable as
#' common-variant inflation estimates; set \code{maf_min = 0} (or omit
#' \code{maf_col}) to use all variants.
#'
#' @param scan A single-variant scan \code{data.frame} (the
#'   \code{\link{marginal_scan}} shape). \code{"ldsc_intercept"} needs columns
#'   \code{chr}, \code{pos}, \code{ref}, \code{alt} and \code{z_col};
#'   \code{"lambda_gc"} needs \code{pvalue_col}. \code{maf_col} is used for QC
#'   when present.
#' @param method \code{"ldsc_intercept"} (default) or \code{"lambda_gc"}.
#' @param ld_scores LD-score table from \code{\link{compute_ld_scores}}
#'   (\code{chr}, \code{pos}, \code{ref}, \code{alt}, \code{ld}); required for
#'   \code{"ldsc_intercept"}, ignored for \code{"lambda_gc"}.
#' @param pvalue_col,z_col Column names for the p-value (\code{"lambda_gc"}) and
#'   the Z-score (\code{"ldsc_intercept"}, \eqn{\chi^2 = Z^2}). Defaults
#'   \code{"pvalue"}, \code{"Z"}.
#' @param maf_col,maf_min Minor-allele-frequency column and threshold for the
#'   common-variant QC (default \code{"MAF"}, \code{0.01}). QC is skipped if
#'   \code{maf_col} is absent or \code{maf_min <= 0}.
#' @param exclude_mhc,mhc_chr,mhc_start,mhc_end Exclude the MHC region (default
#'   TRUE; chr \code{"6"}, \code{25e6}-\code{35e6}, hg38).
#' @param winsorize \eqn{\chi^2} cap passed to \code{\link{ldsc_regression}}
#'   (default 80).
#'
#' @return A one-row \code{data.frame}: \code{method}, \code{factor},
#'   \code{n_variants}, and \code{intercept}, \code{intercept_se}, \code{slope},
#'   \code{slope_se}, \code{confounding_ratio}, \code{ratio_se}, \code{lambda_gc},
#'   \code{mean_chi2} (the LDSC fields are \code{NA} for \code{"lambda_gc"}). The
#'   \code{factor} column is the calibration factor; the row doubles as a
#'   provenance record to persist alongside a run.
#'
#' @seealso \code{\link{ldsc_regression}}, \code{\link{compute_ld_scores}},
#'   \code{\link{compute_lambda_gc}}, \code{\link{calibrate_pvalues}}
#'
#' @examples
#' \dontrun{
#' scan <- read.csv("marginal_all.csv")
#' ld   <- readRDS("ldscores.rds")
#' est  <- estimate_inflation_factor(scan, "ldsc_intercept", ld_scores = ld)
#' est$factor                       # the LDSC intercept 
#' estimate_inflation_factor(scan, "lambda_gc")$factor
#' }
#'
#' @export
estimate_inflation_factor <- function(scan,
                                      method = c("ldsc_intercept", "lambda_gc"),
                                      ld_scores = NULL,
                                      pvalue_col = "pvalue", z_col = "Z",
                                      maf_col = "MAF", maf_min = 0.01,
                                      exclude_mhc = TRUE, mhc_chr = "6",
                                      mhc_start = 25e6, mhc_end = 35e6,
                                      winsorize = 80) {
  method <- match.arg(method)
  if (!is.data.frame(scan)) stop("`scan` must be a data.frame.", call. = FALSE)

  # Common-variant QC (shared by both methods; skipped if no MAF column).
  if (maf_col %in% names(scan) && maf_min > 0) {
    scan <- scan[is.finite(scan[[maf_col]]) & scan[[maf_col]] >= maf_min, ,
                 drop = FALSE]
  }

  # NA-filled template so both methods return the same columns.
  out <- data.frame(method = method, factor = NA_real_, n_variants = 0L,
                    intercept = NA_real_, intercept_se = NA_real_,
                    slope = NA_real_, slope_se = NA_real_,
                    confounding_ratio = NA_real_, ratio_se = NA_real_,
                    lambda_gc = NA_real_, mean_chi2 = NA_real_,
                    stringsAsFactors = FALSE)

  if (method == "lambda_gc") {
    if (!pvalue_col %in% names(scan)) {
      stop("scan has no '", pvalue_col, "' column for method 'lambda_gc'.",
           call. = FALSE)
    }
    p <- scan[[pvalue_col]]
    lam <- compute_lambda_gc(p)
    out$factor     <- lam
    out$lambda_gc  <- lam
    out$n_variants <- sum(is.finite(p) & p > 0 & p <= 1)
    return(out)
  }

  # method == "ldsc_intercept"
  if (is.null(ld_scores)) {
    stop("method 'ldsc_intercept' requires `ld_scores` (from compute_ld_scores()).",
         call. = FALSE)
  }
  need <- c("chr", "pos", "ref", "alt", z_col)
  miss <- setdiff(need, names(scan))
  if (length(miss)) stop("scan is missing column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
  if (!all(c("chr", "pos", "ref", "alt", "ld") %in% names(ld_scores))) {
    stop("`ld_scores` must have chr, pos, ref, alt, ld (compute_ld_scores() shape).",
         call. = FALSE)
  }

  if (isTRUE(exclude_mhc)) {
    in_mhc <- as.character(scan$chr) == as.character(mhc_chr) &
              scan$pos >= mhc_start & scan$pos <= mhc_end
    scan <- scan[!in_mhc, , drop = FALSE]
  }

  key <- function(d) paste(as.character(d$chr), d$pos, d$ref, d$alt, sep = ":")
  scan$.key <- key(scan)
  ld <- ld_scores[is.finite(ld_scores$ld) & ld_scores$ld > 0, , drop = FALSE]
  ld$.key <- key(ld)
  m <- merge(data.frame(.key = scan$.key, chr = as.character(scan$chr),
                        chi2 = scan[[z_col]]^2, stringsAsFactors = FALSE),
             data.frame(.key = ld$.key, ld = ld$ld, stringsAsFactors = FALSE),
             by = ".key")

  fit <- ldsc_regression(m$chi2, m$ld, m$chr, winsorize = winsorize)
  out$factor            <- fit$intercept
  out$n_variants        <- fit$n_variants
  out$intercept         <- fit$intercept
  out$intercept_se      <- fit$intercept_se
  out$slope             <- fit$slope
  out$slope_se          <- fit$slope_se
  out$confounding_ratio <- fit$confounding_ratio
  out$ratio_se          <- fit$ratio_se
  out$lambda_gc         <- fit$lambda_gc
  out$mean_chi2         <- fit$mean_chi2
  out
}
