# This file implements the GLOW Burden test

#' GLOW Burden Test with Optimal Weights
#'
#' @description
#' Performs the Burden test using optimal weights derived from the GLOW framework.
#' The Burden test combines Z-scores linearly and is most powerful when all causal
#' variants have effects in the same direction (all positive or all negative).
#'
#' The test statistic is:
#' \deqn{S_{Burden} = \sum_{i=1}^{m} w_i Z_i}
#'
#' where \eqn{w_i} are optimal weights and \eqn{Z_i} are marginal score Z-scores.
#'
#' **CRITICAL**: The Burden test statistic is a LINEAR combination of Z-scores.
#' Under the null hypothesis, S follows a normal distribution with variance
#' \eqn{\sigma^2 = w^T M w}, where M is the correlation matrix of Z-scores.
#' The p-value is calculated using the normal distribution (NOT Davies method,
#' NOT chi-square distribution).
#'
#' @details
#' **Algorithm**:
#' 1. Compute scaled effect sizes: \eqn{B^* = \sqrt{diag(M_s)} \times B / s_0}
#' 2. Calculate optimal weights using Optimal_Weights_M with:
#'    - Transformation g(x) = x (identity function)
#'    - Degrees of freedom df = Inf (indicates Burden test)
#'    - Effect size estimates B*
#'    - Variant-importance scores PI
#'    - Correlation matrix M
#' 3. Construct weight matrix with optimal weights (BE, APE) and equal weights
#' 4. Call omni_SgZ_test which routes to burden_test for each weight scheme
#' 5. burden_test uses normal distribution: \eqn{p = 2 \Phi(-|S|/\sigma)}
#'
#' **Optimal Weight Schemes**:
#' - **wts_BE** (Best Estimator): Optimal under normal approximation assuming
#'   random effect sizes \eqn{\theta_i \sim PI_i \delta_{B_i} + (1-PI_i) \delta_0}
#' - **wts_APE** (Asymptotically Powerful Estimator): Accounts for correlation
#'   and variance structure, optimal for sparse signals
#' - **wts_equ**: Equal weights baseline (simple sum of Z-scores)
#'
#' The final Burden p-value is the minimum of the equal weights p-value (last row).
#'
#' **When Burden Test is Optimal**:
#' - All causal variants have effect in the same direction
#' - Effect sizes are relatively homogeneous
#' - Strong directional hypothesis (e.g., all deleterious or all protective)
#'
#' **P-value Calculation** (MANDATORY - DO NOT MODIFY):
#' The p-value is computed using the normal distribution:
#' \deqn{p = 2 \times P(N(0, \sigma^2) > |S|) = 2 \times \Phi(-|S|/\sigma)}
#' where \eqn{\sigma^2 = w^T M w} is the variance of S under the null.
#'
#' This is NOT a quadratic form - it is a linear combination. The statistic S
#' is NOT squared. The Davies method is NOT used. The chi-square distribution
#' is NOT used.
#'
#' **Computational Complexity**: O(m^2) for weight calculation and p-value
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
#'   (e.g., method, nsim for GFisher package options)
#'
#' @return A list with two elements:
#' \describe{
#'   \item{STAT}{Matrix (k x 1) of test statistics, where k is the number of
#'     weight schemes (wts_BE, wts_APE, wts_equ). Last row is labeled
#'     "GLOW_Burden" and corresponds to equal weights.}
#'   \item{PVAL}{Matrix (k x 1) of p-values corresponding to each test statistic.
#'     All p-values are computed using the normal distribution. Last row is the
#'     final GLOW Burden test p-value.}
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
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
#' # Run GLOW Burden test
#' result <- GLOW_Burden(marg_score_stats=marg_stats, B=B, PI=PI)
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
#' \code{\link{GLOW_SKAT}} for SKAT test
#' \code{\link{GLOW_Fisher}} for Fisher test
#' \code{\link{GLOW_Omni}} for omnibus test combining all methods
#'
#' @export
GLOW_Burden <- function(marg_score_stats, B, PI, ...) {
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

  # Equal weights baseline (simple sum)
  wts_equ <- rep(1, nrow(M))

  # Compute scaled effect sizes: B* = sqrt(diag(M_s)) * B / s0
  # This standardizes B to be on the same scale as the score statistics
  Bstar <- sqrt(diag(marg_score_stats$M_s)) * B / marg_score_stats$s0

  # Extract Z-scores from marginal score statistics
  Zscores <- marg_score_stats$Zscores

  # ========== BURDEN TEST CONFIGURATION ==========
  # Transformation function: g(x) = x (identity)
  # This is the defining characteristic of the Burden test
  g <- function(x) x

  # Degrees of freedom: Inf
  # This signals to omni_SgZ_test to use normal distribution (NOT Davies)
  stat_df <- Inf

  # Calculate optimal weights for Burden test
  # This uses the special case in Optimal_Weights_M for g(x) = x
  # Returns: wts_BE and wts_APE
  wts_opt_burden <- Optimal_Weights_M(
    g = g,
    Bstar = Bstar,
    PI = PI,
    M = M,
    is.posi.wts = TRUE
  )

  # Construct weight matrix: bind optimal weights and equal weights
  # Each row is a different weighting scheme
  WT_opt_burden <- rbind(t(do.call(cbind, wts_opt_burden)), wts_equ)
  rownames(WT_opt_burden) <- c(names(wts_opt_burden), "wts_equ")

  # Degrees of freedom matrix: all rows use df = Inf (Burden test)
  DF_opt_burden <- matrix(rep(stat_df, nrow(WT_opt_burden)), ncol = 1)

  # Call omnibus test function
  # This will route to burden_test() for each weight scheme
  # burden_test() uses normal distribution for p-values (NOT Davies)
  omni_opt <- omni_SgZ_test(
    Zscores = Zscores,
    DF = DF_opt_burden,
    W = WT_opt_burden,
    M = M,
    ...
  )

  # Rename the last row (equal weights) to "GLOW_Burden"
  # This is the final reported Burden test result
  rownames(omni_opt$STAT)[nrow(omni_opt$STAT)] <- "GLOW_Burden"
  rownames(omni_opt$PVAL)[nrow(omni_opt$PVAL)] <- "GLOW_Burden"

  # Return test statistics and p-values
  return(list(STAT = omni_opt$STAT, PVAL = omni_opt$PVAL))
}
