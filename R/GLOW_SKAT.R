# This file implements the GLOW SKAT test

#' GLOW SKAT Test with Optimal Weights
#'
#' @description
#' Performs the SKAT (Sequence Kernel Association Test) using optimal weights
#' derived from the GLOW framework. SKAT is most powerful when causal variants
#' have effects in mixed directions (some positive, some negative) or when
#' effect sizes are heterogeneous.
#'
#' The test statistic is:
#' \deqn{S_{SKAT} = \sum_{i=1}^{m} w_i Z_i^2}
#'
#' where \eqn{w_i} are optimal weights and \eqn{Z_i} are marginal score Z-scores.
#'
#' **CRITICAL**: The SKAT test statistic is a QUADRATIC form (sum of SQUARED
#' Z-scores). The p-value is calculated using p.GFisher from the GFisher package
#' with df=1. This is NOT the Liu method. This is NOT liu.mod. The test uses
#' the GFisher framework which properly accounts for correlation among Z-scores.
#'
#' @details
#' **Algorithm** (EXACT port from legacy lines 26-52):
#' 1. Compute scaled effect sizes: \eqn{B^* = \sqrt{diag(M_s)} \times B / s_0}
#' 2. Define transformation function: g(x) = x^2 (SQUARING function)
#' 3. Set degrees of freedom: df = 1 (SKAT uses df=1)
#' 4. Calculate optimal weights using Optimal_Weights_M with:
#'    - Transformation g(x) = x^2
#'    - Effect size estimates B*
#'    - Variant-importance scores PI
#'    - Correlation matrix M
#'    - is.posi.wts = TRUE (forces non-negative weights)
#' 5. Construct weight matrix with optimal weights and equal weights
#' 6. Call omni_SgZ_test which routes to skat_test() for df=1
#' 7. skat_test() calls p.GFisher with df=1 (NOT Liu method)
#'
#' **Optimal Weight Schemes**:
#' - **wts_BE_N**: Best Estimator using normal approximation
#' - **wts_APE_N**: Asymptotically Powerful Estimator using normal approximation
#' - **wts_BE_sparse**: Best Estimator using sparse mixture approximation
#' - **wts_APE_sparse**: APE using sparse mixture approximation
#' - **wts_equ**: Equal weights baseline (simple sum of Z^2)
#'
#' The final SKAT p-value is from the equal weights row (last row).
#'
#' **When SKAT Test is Optimal**:
#' - Causal variants have effects in mixed directions (some protective, some deleterious)
#' - Effect sizes are heterogeneous across variants
#' - No strong prior on direction of effects
#' - Variance component model is appropriate
#'
#' **P-value Calculation** (MANDATORY - DO NOT MODIFY):
#' The p-value is computed using p.GFisher from the GFisher package with df=1:
#' \deqn{p = p.GFisher(q = S, df = 1, w = weights, M = correlation, ...)}
#'
#' This is a quadratic form test. The statistic S is a weighted sum of SQUARED
#' Z-scores. The GFisher framework with df=1 provides the correct null distribution
#' accounting for the correlation structure.
#'
#' **What this is NOT**:
#' - NOT Liu's method for weighted chi-square mixtures
#' - NOT liu.mod method
#' - NOT Davies method applied directly
#' - The transformation is g(x) = x^2 (NOT g(x) = x as in Burden)
#' - The df = 1 (NOT df = Inf as in Burden, NOT df = 2 as in Fisher)
#'
#' **Computational Complexity**: O(m^2) for weight calculation, O(m^2) for p-value
#'
#' **Assumptions**:
#' - Z-scores are standard normal under null hypothesis
#' - M is the correlation matrix of Z-scores
#' - Effect sizes B are on the standardized scale
#' - Variant-importance scores PI are in \eqn{[0,1]}
#'
#' @param marg_score_stats List output from \code{\link{getZ_marg_score}} or
#'   \code{\link{getZ_marg_score_binary_SPA}} containing:
#'   \describe{
#'     \item{Zscores}{Numeric vector of marginal score Z-scores}
#'     \item{M_Z}{Correlation matrix of Z-scores}
#'     \item{M_s}{Covariance matrix of score statistics}
#'     \item{s0}{Scalar standard deviation under null}
#'   }
#' @param B Numeric vector of effect size estimates (length m). These can be
#'   from external studies, meta-analysis, or prior information. Should be on
#'   the standardized scale (e.g., log odds ratio, standardized beta).
#' @param PI Numeric vector of variant-importance scores (length m, values in
#'   \eqn{[0,1]}), each variant's annotation-derived relative importance. Can be
#'   computed using \code{\link{get_PI}} based on annotation data.
#' @param ... Additional arguments passed to the internal \code{omni_SgZ_test}
#'   and eventually to p.GFisher (e.g., method, nsim for GFisher package options)
#'
#' @return A list with two elements:
#' \describe{
#'   \item{STAT}{Matrix (k x 1) of test statistics, where k is the number of
#'     weight schemes (4 optimal + 1 equal). Last row is labeled "GLOW_SKAT"
#'     and corresponds to equal weights.}
#'   \item{PVAL}{Matrix (k x 1) of p-values corresponding to each test statistic.
#'     All p-values are computed using p.GFisher with df=1. Last row is the
#'     final GLOW SKAT test p-value.}
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., & Lin, X. (2011).
#' Rare-variant association testing for sequencing data with the sequence kernel
#' association test. The American Journal of Human Genetics, 89(1), 82-93.
#'
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and
#' accurate p-value calculation under dependence. Biometrics, 79(2), 1159-1172.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n <- 500  # sample size
#' m <- 20   # number of variants
#'
#' # Genotype matrix (n x m)
#' G <- matrix(rbinom(n*m, 2, 0.1), n, m)
#'
#' # Covariates
#' X <- matrix(rnorm(n*2), n, 2)
#'
#' # Binary trait
#' Y <- rbinom(n, 1, 0.3)
#'
#' # Compute marginal score statistics
#' marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
#'
#' # External effect sizes (e.g., from meta-analysis)
#' B <- rnorm(m, mean=0, sd=0.2)
#'
#' # Variant-importance scores (e.g., from annotation)
#' PI <- runif(m, 0.1, 0.9)
#'
#' # Run GLOW SKAT test
#' result <- GLOW_SKAT(marg_score_stats=marg_stats, B=B, PI=PI)
#'
#' # View results
#' print(result$STAT)
#' print(result$PVAL)
#' }
#'
#' @note For new code, consider using \code{\link{glow_test}} with appropriate
#'   \code{test_specs} instead. The \code{glow_test()} interface supports
#'   flexible test composition via \code{\link{default_test_specs}} and is the
#'   preferred entry point for the GLOW test system. This function remains
#'   available for direct use and backward compatibility.
#'
#' @seealso
#' \code{\link{glow_test}} for the preferred high-level interface
#' \code{\link{getZ_marg_score}} for computing marginal score statistics
#' \code{\link{getZ_marg_score_binary_SPA}} for binary traits with SPA
#' \code{\link{get_B}} for estimating effect sizes from external data
#' \code{\link{get_PI}} for computing variant-importance scores from annotations
#' \code{\link{Optimal_Weights_M}} for optimal weight calculation
#' \code{\link{GLOW_Burden}} for Burden test
#' \code{\link{GLOW_Fisher}} for Fisher test
#' \code{\link{GLOW_Omni}} for omnibus test combining all methods
#'
#' @export
GLOW_SKAT <- function(marg_score_stats, B, PI, ...) {
  # Input validation
  m <- nrow(marg_score_stats$M_Z)  # Number of variants

  if (length(B) != m) {
    stop("B must have length equal to the number of variants (", m, ")")
  }
  if (length(PI) != m) {
    stop("PI must have length equal to the number of variants (", m, ")")
  }

  # Extract correlation matrix of Z-scores
  M <- marg_score_stats$M_Z

  # Equal weights baseline
  wts_equ <- rep(1, nrow(M))

  # Compute scaled effect sizes: B* = sqrt(diag(M_s)) * B / s0
  # This standardizes B to be on the same scale as the score statistics
  Bstar <- sqrt(diag(marg_score_stats$M_s)) * B / marg_score_stats$s0

  # Extract Z-scores from marginal score statistics
  Zscores <- marg_score_stats$Zscores

  # ========== SKAT TEST CONFIGURATION ==========
  # CRITICAL: This is where we define SKAT
  # Transformation function: g(x) = x^2 (SQUARING, not identity)
  # This is the defining characteristic of the SKAT test
  g <- function(x) x^2

  # Degrees of freedom: 1 (NOT Inf, NOT 2)
  # This signals to omni_SgZ_test to use p.GFisher with df=1
  # df=1 is specifically for SKAT (quadratic form)
  stat_df <- 1

  # Calculate optimal weights for SKAT test
  # Because g(x) = x^2 (NOT identity), this will return 4 weight schemes:
  # wts_BE_N, wts_APE_N, wts_BE_sparse, wts_APE_sparse
  wts_opt_skat <- Optimal_Weights_M(
    g = g,
    Bstar = Bstar,
    PI = PI,
    M = M,
    is.posi.wts = TRUE
  )

  # Construct weight matrix: bind optimal weights and equal weights
  # Each row is a different weighting scheme
  WT_opt_skat <- rbind(t(do.call(cbind, wts_opt_skat)), wts_equ)
  rownames(WT_opt_skat) <- c(names(wts_opt_skat), "wts_equ")

  # Degrees of freedom matrix: all rows use df = 1 (SKAT test)
  DF_opt_skat <- matrix(rep(stat_df, nrow(WT_opt_skat)), ncol = 1)

  # Call omnibus test function
  # This will route to calcu_SgZ_p() with df=1
  # calcu_SgZ_p() will call p.GFisher with df=1 (NOT Liu method)
  omni_opt <- omni_SgZ_test(
    Zscores = Zscores,
    DF = DF_opt_skat,
    W = WT_opt_skat,
    M = M,
    ...
  )

  # Rename the last row (equal weights) to "GLOW_SKAT"
  # This is the final reported SKAT test result
  rownames(omni_opt$STAT)[nrow(omni_opt$STAT)] <- "GLOW_SKAT"
  rownames(omni_opt$PVAL)[nrow(omni_opt$PVAL)] <- "GLOW_SKAT"

  # Return test statistics and p-values
  return(list(STAT = omni_opt$STAT, PVAL = omni_opt$PVAL))
}
