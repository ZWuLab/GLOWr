# This file implements the GLOW Fisher test

#' GLOW Fisher Test with Optimal Weights
#'
#' @description
#' Performs the Fisher combination test using optimal weights derived from the
#' GLOW framework. The Fisher test is most powerful when causal variants have
#' mixed effect directions but moderate heterogeneity, providing a middle ground
#' between Burden and SKAT tests.
#'
#' The Fisher test statistic is:
#' \deqn{S_{Fisher} = \sum_{i=1}^{m} w_i g(Z_i)}
#'
#' where \eqn{w_i} are optimal weights, \eqn{Z_i} are marginal score Z-scores,
#' and \eqn{g(\cdot)} is the GFisher transformation function with df=2.
#'
#' **CRITICAL**: The Fisher test uses the g_GFisher transformation with df=2.
#' The p-value is calculated using p.GFisher from the GFisher package with df=2.
#' This is NOT the Liu method. This is NOT the Davies method applied directly.
#' The transformation converts Z-scores to chi-square statistics on the df=2 scale.
#'
#' @details
#' **Algorithm** (EXACT port from legacy lines 26-52):
#' 1. Compute scaled effect sizes: \eqn{B^* = \sqrt{diag(M_s)} \times B / s_0}
#' 2. Define transformation function: g(x) = g_GFisher(x, df=2, p.type)
#' 3. Set degrees of freedom: df = 2 (Fisher test uses df=2)
#' 4. Calculate optimal weights using Optimal_Weights_M with:
#'    - Transformation g(x) = g_GFisher(x, df=2)
#'    - Effect size estimates B*
#'    - Variant-importance scores PI
#'    - Correlation matrix M
#'    - is.posi.wts = TRUE (forces non-negative weights)
#' 5. Construct weight matrix with optimal weights and equal weights
#' 6. Call omni_SgZ_test which routes to calcu_SgZ_p() with df=2
#' 7. calcu_SgZ_p() calls p.GFisher with df=2
#'
#' **Optimal Weight Schemes**:
#' - **wts_BE_N**: Best Estimator using normal approximation
#' - **wts_APE_N**: Asymptotically Powerful Estimator using normal approximation
#' - **wts_BE_sparse**: Best Estimator using sparse mixture approximation
#' - **wts_APE_sparse**: APE using sparse mixture approximation
#' - **wts_equ**: Equal weights baseline (simple sum)
#'
#' The final Fisher p-value is from the equal weights row (last row).
#'
#' **When Fisher Test is Optimal**:
#' - Causal variants have mixed effect directions (some positive, some negative)
#' - Effect sizes have moderate heterogeneity (less extreme than SKAT scenario)
#' - Combines benefits of Burden and SKAT approaches
#' - More robust than Burden when effect directions are mixed
#' - More powerful than SKAT when effects are moderately homogeneous
#'
#' **P-value Calculation** (MANDATORY - DO NOT MODIFY):
#' The p-value is computed using p.GFisher from the GFisher package with df=2:
#' \deqn{p = p.GFisher(q = S, df = 2, w = weights, M = correlation, p.type, ...)}
#'
#' The transformation function is:
#' \deqn{g(Z_i) = F_{\chi^2_2}^{-1}(\log(2) + \log(\bar{\Phi}(|Z_i|)))}
#'
#' for two-sided tests, where \eqn{F_{\chi^2_2}^{-1}} is the inverse CDF of the
#' chi-square distribution with df=2, and \eqn{\bar{\Phi}} is the upper tail of
#' the standard normal.
#'
#' **What this is NOT**:
#' - NOT Liu's method for weighted chi-square mixtures
#' - NOT Davies method applied directly
#' - The transformation uses g_GFisher with df=2 (NOT df=1 as in SKAT, NOT df=Inf as in Burden)
#' - The df = 2 (NOT df=1, NOT df=Inf)
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
#'   and eventually to p.GFisher (e.g., p.type="two" or "one", method, nsim for
#'   GFisher package options)
#'
#' @return A list with two elements:
#' \describe{
#'   \item{STAT}{Matrix (k x 1) of test statistics, where k is the number of
#'     weight schemes (4 optimal + 1 equal = 5 rows). Last row is labeled
#'     "GLOW_Fisher" and corresponds to equal weights.}
#'   \item{PVAL}{Matrix (k x 1) of p-values corresponding to each test statistic.
#'     All p-values are computed using p.GFisher with df=2. Last row is the
#'     final GLOW Fisher test p-value.}
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' Fisher, R. A. (1932). Statistical Methods for Research Workers (4th ed.).
#' Oliver and Boyd.
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
#' # Run GLOW Fisher test (two-sided, default)
#' result <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI)
#'
#' # Run one-sided Fisher test
#' result_one <- GLOW_Fisher(marg_score_stats=marg_stats, B=B, PI=PI, p.type="one")
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
#' \code{\link{GLOW_SKAT}} for SKAT test
#' \code{\link{GLOW_Omni}} for omnibus test combining all methods
#'
#' @export
GLOW_Fisher <- function(marg_score_stats, B, PI, ...) {
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

  # ========== FISHER TEST CONFIGURATION ==========
  # CRITICAL: This is where we define Fisher test
  # Transformation function: g(x) = g_GFisher(x, df=2, p.type)
  # This is the defining characteristic of the Fisher test
  # NOTE: The legacy code uses g_GFisher_two directly, which corresponds to
  # g_GFisher with p.type="two" (two-sided test)
  g <- function(x, df = 2) g_GFisher_two(x, df)

  # Degrees of freedom: 2 (NOT 1, NOT Inf)
  # This signals to omni_SgZ_test to use p.GFisher with df=2
  # df=2 is specifically for Fisher test
  stat_df <- 2

  # Calculate optimal weights for Fisher test
  # Because g uses g_GFisher transformation with df=2,
  # this will return weight schemes: wts_BE, wts_APR
  wts_opt_fisher <- Optimal_Weights_M(
    g = g,
    Bstar = Bstar,
    PI = PI,
    M = M,
    is.posi.wts = TRUE
  )

  # Construct weight matrix: bind optimal weights and equal weights
  # Each row is a different weighting scheme
  WT_opt_fisher <- rbind(t(do.call(cbind, wts_opt_fisher)), wts_equ)
  rownames(WT_opt_fisher) <- c(names(wts_opt_fisher), "wts_equ")

  # Degrees of freedom matrix: all rows use df = 2 (Fisher test)
  DF_opt_fisher <- matrix(rep(stat_df, nrow(WT_opt_fisher)), ncol = 1)

  # Call omnibus test function
  # This will route to calcu_SgZ_p() with df=2
  # calcu_SgZ_p() will call p.GFisher with df=2
  omni_opt <- omni_SgZ_test(
    Zscores = Zscores,
    DF = DF_opt_fisher,
    W = WT_opt_fisher,
    M = M,
    ...
  )

  # Rename the last row (equal weights) to "GLOW_Fisher"
  # This is the final reported Fisher test result
  rownames(omni_opt$STAT)[nrow(omni_opt$STAT)] <- "GLOW_Fisher"
  rownames(omni_opt$PVAL)[nrow(omni_opt$PVAL)] <- "GLOW_Fisher"

  # Return test statistics and p-values
  return(list(STAT = omni_opt$STAT, PVAL = omni_opt$PVAL))
}
