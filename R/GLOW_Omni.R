# This file implements the GLOW Omni test

#' GLOW Omnibus Test Combining Burden, SKAT, and Fisher
#'
#' @description
#' Performs an omnibus test combining Burden, SKAT, and Fisher tests using the
#' Cauchy Combination Test (CCT). The omnibus test adapts to different genetic
#' architectures by combining multiple test statistics, providing robust power
#' across various scenarios.
#'
#' GLOW_Omni computes:
#' 1. Burden test (optimal + equal weights)
#' 2. SKAT test (optimal + equal weights)
#' 3. Fisher test (optimal + equal weights)
#' 4. CCT of SNV-level p-values (variant-level minimum p-value combination)
#' 5. CCT combining all the above tests (final omnibus p-value)
#'
#' The final omnibus p-value provides protection against model misspecification
#' and is powerful across diverse genetic architectures.
#'
#' @details
#' **Algorithm**:
#'
#' Internally delegates to \code{run_bsf_tests()}, which:
#' 1. For each test family in \code{test_specs}, computes optimal and equal
#'    weights, then evaluates the weighted test statistic.
#' 2. Combines all individual test p-values via CCT (BSF omnibus).
#' 3. Computes SNV-level CCT from marginal p-values (if
#'    \code{include_snv_cct = TRUE}).
#' 4. Combines individual + SNV CCT for final omnibus p-value.
#'
#' **Transformation Details**:
#' - Scaled effect sizes: \eqn{B^* = \sqrt{diag(M_s)} \times B / s_0}
#' - Burden uses g(x) = x, df = Inf
#' - SKAT uses g(x) = x^2, df = 1
#' - Fisher uses g(x) = g_GFisher_two(x, df=2), df = 2
#'
#' **Output Structure** (default 3-family spec, 16 rows):
#' - Rows 1-5: SKAT tests (BE_N, APR_N, BE_sparse, APR_sparse, equ)
#' - Rows 6-8: Burden tests (BE, APR, equ)
#' - Rows 9-13: Fisher tests (BE_N, APR_N, BE_sparse, APR_sparse, equ)
#' - Row 14: BSF Omni (CCT of rows 1-13)
#' - Row 15: SNV-level CCT
#' - Row 16: Final omnibus CCT
#'
#' With custom \code{test_specs}, the number of rows varies.
#'
#' **Computational Complexity**: O(k * m^2) where k is number of weight
#' schemes, m is number of variants
#'
#' @param marg_score_stats List output from \code{\link{getZ_marg_score}} or
#'   \code{\link{getZ_marg_score_binary_SPA}} containing:
#'   \describe{
#'     \item{Zscores}{Numeric vector of marginal score Z-scores}
#'     \item{M_Z}{Correlation matrix of Z-scores}
#'     \item{M_s}{Covariance matrix of score statistics}
#'     \item{s0}{Scalar standard deviation under null}
#'   }
#' @param B Numeric vector of effect size estimates (length m).
#' @param PI Numeric vector of variant-importance scores (length m, values in
#'   \eqn{[0,1]}).
#' @param test_specs List of test specification lists (default: NULL uses
#'   \code{\link{default_test_specs}}, the standard 3-family configuration).
#'   See \code{\link{default_test_specs}} for the spec format.
#' @param include_snv_cct Logical. If TRUE (default), append SNV-level CCT
#'   and final omnibus CCT rows to the output.
#' @param return_weights Logical (default FALSE). If TRUE, attach the
#'   internally-assembled weight matrix to the result (forwarded from
#'   \code{run_bsf_tests}). Rows: weight schemes; columns: variants.
#' @param ... Additional arguments passed to \code{omni_SgZ_test} and weight
#'   functions (e.g., method, nsim).
#'
#' @return A list with:
#' \describe{
#'   \item{STAT}{Matrix (n x 1) of test statistics. With default specs, n=16.}
#'   \item{PVAL}{Matrix (n x 1) of p-values. Last row is the final omnibus.}
#'   \item{test_names}{Character vector of row names.}
#'   \item{weights}{(Only when \code{return_weights = TRUE}.) Weight matrix
#'     (n_schemes x p).}
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' Liu Y, Xie J (2020). "Cauchy Combination Test: A Powerful Test With Analytic
#' p-Value Calculation Under Arbitrary Dependency Structures."
#' JASA, 115(529), 393-402.
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n <- 500; m <- 20
#' G <- matrix(rbinom(n*m, 2, 0.1), n, m)
#' X <- matrix(rnorm(n*2), n, 2)
#' Y <- rbinom(n, 1, 0.3)
#' marg_stats <- getZ_marg_score(G, X, Y, trait="binary")
#' B <- rnorm(m, mean=0, sd=0.2)
#' PI <- runif(m, 0.1, 0.9)
#'
#' # Standard omnibus (16 rows)
#' result <- GLOW_Omni(marg_stats, B, PI)
#' result$PVAL[nrow(result$PVAL), ]
#'
#' # Custom 2-family spec
#' my_specs <- list(
#'   list(family = "SKAT", g = function(x) x^2, df = 1),
#'   list(family = "Burden", g = function(x) x, df = Inf)
#' )
#' result2 <- GLOW_Omni(marg_stats, B, PI, test_specs = my_specs)
#' }
#'
#' @seealso
#' \code{\link{default_test_specs}} for the test specification format
#' \code{\link{glow_test}} for the high-level test runner
#' \code{\link{GLOW_Omni_byP}} for omnibus test from pre-computed p-values
#'
#' @export
GLOW_Omni <- function(marg_score_stats, B, PI,
                       test_specs = NULL,
                       include_snv_cct = TRUE,
                       return_weights = FALSE,
                       ...) {
  # Input validation
  m <- nrow(marg_score_stats$M_Z)

  if (length(B) != m) {
    stop("B must have length equal to the number of variants (", m, ")")
  }
  if (length(PI) != m) {
    stop("PI must have length equal to the number of variants (", m, ")")
  }

  # Resolve test specs
  specs <- if (is.null(test_specs)) default_test_specs() else test_specs

  # Compute scaled effect sizes: B* = sqrt(diag(M_s)) * B / s0
  Bstar <- sqrt(diag(marg_score_stats$M_s)) * B / marg_score_stats$s0

  # Delegate to data-driven composition layer
  run_bsf_tests(
    Zscores = marg_score_stats$Zscores,
    M = marg_score_stats$M_Z,
    Bstar = Bstar,
    PI = PI,
    test_specs = specs,
    include_snv_cct = include_snv_cct,
    return_weights = return_weights,
    ...
  )
}


#' GLOW Omnibus Test from Pre-Computed P-values
#'
#' @description
#' Performs the GLOW omnibus test when individual variant p-values are already
#' computed. Reconstructs Z-scores from p-values and their signs, then runs
#' the omnibus test pipeline via \code{run_bsf_tests()}.
#'
#' @details
#' **Algorithm**:
#'
#' 1. Reconstruct Z-scores: \eqn{Z = \Phi^{-1}(1 - p/2) \times sign}
#' 2. Call \code{run_bsf_tests()} with reconstructed Z-scores
#'
#' **CRITICAL**: The correlation matrix M must reflect the Z-score correlation
#' structure, not the adjusted p-value structure.
#'
#' @param Pvalues Numeric vector of variant-level p-values (length m, in
#'   \eqn{[0,1]}).
#' @param Zsigns Numeric vector of effect direction signs (length m, +1 or -1).
#' @param M Correlation matrix of Z-scores (m x m).
#' @param Bstar Numeric vector of scaled effect sizes (length m).
#' @param PI Numeric vector of variant-importance scores (length m, in \eqn{[0,1]}).
#' @param test_specs List of test spec lists (default: NULL uses
#'   \code{\link{default_test_specs}}).
#' @param include_snv_cct Logical (default TRUE). Append SNV-level and final
#'   omnibus CCT rows.
#' @param return_weights Logical (default FALSE). If TRUE, attach the
#'   internally-assembled weight matrix to the result (forwarded from
#'   \code{run_bsf_tests}).
#' @param ... Additional arguments passed to \code{run_bsf_tests}.
#'
#' @return Same structure as \code{\link{GLOW_Omni}}.
#'
#' @seealso
#' \code{\link{GLOW_Omni}} for omnibus test with Z-scores
#' \code{\link{getZ_marg_score_binary_SPA}} for SPA-adjusted p-values
#'
#' @export
GLOW_Omni_byP <- function(Pvalues, Zsigns, M, Bstar, PI,
                            test_specs = NULL,
                            include_snv_cct = TRUE,
                            return_weights = FALSE,
                            ...) {
  # Input validation
  m <- length(Pvalues)

  if (length(Zsigns) != m) {
    stop("Zsigns must have the same length as Pvalues (", m, ")")
  }
  if (!all(Zsigns %in% c(-1, 1))) {
    warning("Zsigns should be +1 or -1; other values may give unexpected results")
  }
  if (nrow(M) != m || ncol(M) != m) {
    stop("M must be a ", m, " x ", m, " correlation matrix")
  }
  if (length(Bstar) != m) {
    stop("Bstar must have length equal to the number of variants (", m, ")")
  }
  if (length(PI) != m) {
    stop("PI must have length equal to the number of variants (", m, ")")
  }

  # Resolve test specs
  specs <- if (is.null(test_specs)) default_test_specs() else test_specs

  # Reconstruct Z-scores from p-values
  Zscores <- qnorm(Pvalues / 2, lower.tail = FALSE) * Zsigns

  # Delegate to data-driven composition layer
  run_bsf_tests(
    Zscores = Zscores,
    M = M,
    Bstar = Bstar,
    PI = PI,
    test_specs = specs,
    include_snv_cct = include_snv_cct,
    return_weights = return_weights,
    ...
  )
}


