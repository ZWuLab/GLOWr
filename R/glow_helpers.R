# This file contains helper functions for GLOW variant-set tests.

#' Cauchy Combination Test
#'
#' Combines multiple p-values using the Cauchy Combination Test (CCT).
#' This is an internal helper function used in GLOW_Omni test to combine
#' p-values from different variant-set tests (Burden, SKAT, Fisher).
#'
#' @description
#' The Cauchy Combination Test transforms p-values using the Cauchy distribution
#' and combines them via averaging. The method is robust to correlation among
#' tests and maintains valid type I error control.
#'
#' The CCT statistic is computed as:
#' \deqn{T_{CCT} = \frac{1}{n} \sum_{i=1}^{n} \tan\left(\pi(0.5 - p_i)\right)}
#'
#' where \eqn{p_i} are the individual p-values. Under the null hypothesis,
#' \eqn{T_{CCT}} follows a standard Cauchy distribution. The combined p-value
#' is obtained as \eqn{P(C > T_{CCT})} where \eqn{C \sim Cauchy(0,1)}.
#'
#' For numerical stability:
#' - Very small p-values (< thr.smallp): Use special formula \eqn{1/(\pi \cdot p_i)}
#' - Large p-values (> thr.largp): Truncate to thr.largp
#'
#' @details
#' The Cauchy transformation handles extreme p-values gracefully:
#'
#' 1. **Regular p-values** (thr.smallp \eqn{\le} p \eqn{\le} thr.largp):
#'    Transform using \eqn{\tan(\pi(0.5 - p))}
#'
#' 2. **Very small p-values** (p < thr.smallp):
#'    Use approximation \eqn{1/(\pi \cdot p)} to avoid numerical issues
#'    with tan() near the asymptote
#'
#' 3. **Large p-values** (p > thr.largp):
#'    Truncate to thr.largp to prevent negative contributions that could
#'    mask significant signals
#'
#' The CCT is particularly useful for combining dependent tests because it
#' maintains valid type I error control without requiring knowledge of the
#' correlation structure among tests.
#'
#' **Computational Complexity**: O(n) where n is the number of p-values
#'
#' **Assumptions**:
#' - Input p-values are valid (between 0 and 1)
#' - P-values represent valid hypothesis tests
#'
#' @param PVAL Numeric vector of p-values to combine (must be in [0,1])
#' @param thr.largp Numeric threshold for large p-values. P-values exceeding
#'   this threshold are truncated to thr.largp. Default: 0.9
#' @param thr.smallp Numeric threshold for very small p-values. P-values below
#'   this threshold use a special stable formula. Default: 1e-15
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cct}{The CCT test statistic (numeric)}
#'   \item{pval_cct}{The combined p-value from the Cauchy distribution (numeric)}
#' }
#'
#' @references
#' Liu Y, Xie J (2020). "Cauchy Combination Test: A Powerful Test With Analytic
#' p-Value Calculation Under Arbitrary Dependency Structures."
#' Journal of the American Statistical Association, 115(529), 393-402.
#' \doi{10.1080/01621459.2018.1554485}
#'
#' @examples
#' # Combine three independent p-values
#' p_vals <- c(0.01, 0.05, 0.10)
#' result <- cct_test(p_vals)
#' result$pval_cct  # Combined p-value
#'
#' # Handle very small p-values
#' p_vals_small <- c(1e-20, 1e-18, 0.05)
#' result_small <- cct_test(p_vals_small)
#'
#' # Handle large p-values (will be truncated)
#' p_vals_large <- c(0.001, 0.95, 0.99)
#' result_large <- cct_test(p_vals_large)
#'
#' @noRd
cct_test <- function(PVAL, thr.largp = 0.9, thr.smallp = 1e-15) {
  # Truncate large p-values to threshold
  # This prevents negative contributions that could mask significant signals
  PVAL[PVAL > thr.largp] <- thr.largp

  # Identify very small p-values that need special handling
  is.small <- (PVAL < thr.smallp)

  # Initialize statistic vector
  CCTSTAT <- PVAL

  # Case 1: No very small p-values - use standard formula for all
  if (sum(is.small) == 0) {
    # Standard Cauchy transformation: tan(pi * (0.5 - p))
    # Mean gives equal weight (1/n) to each p-value
    cct <- mean(tan(pi * (0.5 - CCTSTAT)))
  }

  # Case 2: Some very small p-values - use mixed formula
  if (sum(is.small) > 0) {
    # Regular p-values: standard Cauchy transformation
    CCTSTAT[!is.small] <- tan((0.5 - CCTSTAT[!is.small]) * pi)

    # Very small p-values: use stable approximation 1/(pi * p)
    # This avoids numerical issues with tan() near the asymptote at pi/2
    CCTSTAT[is.small] <- 1 / CCTSTAT[is.small] / pi

    # Combine all transformed values
    cct <- mean(CCTSTAT)
  }

  # Compute combined p-value from Cauchy distribution
  # Upper tail: P(Cauchy > cct)
  pval_cct <- pcauchy(cct, lower.tail = FALSE)

  # Return both statistic and p-value
  return(list(cct = cct, pval_cct = pval_cct))
}


#' Generic S(g(Z)) statistic calculation
#'
#' Computes the weighted sum statistic S = sum(w_i * g(Z_i)) and optionally
#' calculates the p-value. This is the core computation underlying the Burden,
#' SKAT, and Fisher tests.
#'
#' @description
#' This function computes a general weighted combination statistic of transformed
#' Z-scores: \deqn{S = \sum_{i=1}^{n} w_i g(Z_i)}
#'
#' where:
#' - \eqn{Z_i} are standard normal Z-scores
#' - \eqn{w_i} are weights (optionally constrained to be non-negative)
#' - \eqn{g(\cdot)} is a transformation function
#'
#' The function routes to different p-value calculation methods based on the
#' transformation function and degrees of freedom:
#'
#' **Burden test (df=Inf or g(x)=x)**:
#' - Uses normal distribution
#' - Variance: \eqn{\sigma^2 = w^T M w}
#' - P-value: \eqn{P(|N(0,\sigma^2)| > |S|)}
#'
#' **GFisher tests (df=1,2,etc)**:
#' - Uses p.GFisher from GFisher package
#' - Accounts for correlation via M
#' - Can handle one-sided or two-sided p-values
#'
#' @details
#' **CRITICAL ROUTING LOGIC** (this determines p-value calculation):
#'
#' 1. **Burden test identification**:
#'    - If df = Inf OR
#'    - If g is identity function (g(x) = x) OR
#'    - If formals(g)$df = Inf
#'    - Then: Use normal distribution with variance = t(wts) %*% M %*% wts
#'
#' 2. **GFisher test identification**:
#'    - If df is finite (not Inf) OR
#'    - If formals(g)$df is finite
#'    - Then: Call p.GFisher with specified df
#'
#' **Weight handling**:
#' - If is.posi.wts=TRUE: negative weights are set to 0
#' - Weights are always scaled to sum to 1 (if non-zero)
#' - For GFisher tests, weights are re-forced to be non-negative
#'
#' **Computational Complexity**: O(n) for statistic, O(n^2) for p-value (correlation)
#'
#' **Assumptions**:
#' - Zscores follow standard normal under null
#' - M is the correlation matrix of Zscores
#' - For GFisher: weights must be non-negative
#'
#' @param g Function to transform Z-scores. Common choices:
#'   - function(x) x: Burden test
#'   - function(x) x^2: SKAT test
#'   - function(x, df=2) g_GFisher(x, df, p.type): Fisher test
#' @param Zscores Numeric vector of Z-scores (standard normal under null)
#' @param wts Numeric vector of weights (same length as Zscores)
#' @param calc_p Logical; if TRUE, calculate p-value. Default: FALSE
#' @param M Correlation matrix of the Z-scores. Required if calc_p=TRUE
#' @param df Degrees of freedom. Use Inf for burden, 1 for SKAT, 2 for Fisher
#' @param p.type Character string: "two" for two-sided (default), "one" for one-sided
#' @param is.posi.wts Logical; if TRUE, force weights to be non-negative. Default: TRUE
#' @param method Method for p.GFisher: "HYB" (default), "MR", or "GB"
#' @param nsim Number of simulations for method="MR". Default: NULL
#' @param seed Random seed for simulations. Default: NULL
#'
#' @return A list with elements:
#' \describe{
#'   \item{S}{The test statistic (numeric)}
#'   \item{p}{The p-value (numeric, only if calc_p=TRUE)}
#' }
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and
#' accurate p-value calculation under dependence. Biometrics, 79(2), 1159-1172.
#'
#' @examples
#' # Burden test (df=Inf, uses normal distribution)
#' Z <- rnorm(10)
#' M <- diag(10)
#' result <- calcu_SgZ_p(function(x) x, Z, rep(1,10), calc_p=TRUE, M=M, df=Inf)
#'
#' # SKAT test (df=1, uses p.GFisher)
#' result <- calcu_SgZ_p(function(x) x^2, Z, rep(1,10), calc_p=TRUE, M=M, df=1)
#'
#' @noRd
calcu_SgZ_p <- function(g, Zscores, wts, calc_p = FALSE, M = NULL, df = NULL,
                        p.type = "two", is.posi.wts = TRUE,
                        method = "HYB", nsim = NULL, seed = NULL) {

  # Force weights to be non-negative if requested
  if (is.posi.wts) {
    wts <- pmax(wts, 0)
  }

  # Always scale the weights to sum to 1 (if non-zero)
  if (mean(abs(wts)) > 0) {
    wts <- wts / sum(abs(wts))
  } else {
    # If all weights are zero, return S=0 and p=1
    if (calc_p) {
      return(list(S = 0, p = 1))
    } else {
      return(list(S = 0))
    }
  }

  # Calculate the statistic if g is provided
  if (!is.null(g)) {
    S <- sum(wts * g(Zscores))
  }

  # Calculate the p-value if requested
  if (calc_p) {
    # Separate the calculation for burden test (df=Inf) and GFisher test

    # Check for burden test: df=Inf OR g is identity OR formals(g)$df=Inf
    if (identical(df, Inf) ||
        identical(body(g), body(function(x) x)) ||
        identical(formals(g)$df, Inf)) {

      # BURDEN TEST: Use normal distribution (NOT Davies method)
      # S = sum(wts * Zscores) is a linear combination
      S <- sum(wts * Zscores)

      # Variance of S under null: Var(S) = t(wts) %*% M %*% wts
      S_sd <- sqrt(t(wts) %*% M %*% wts)

      # Two-sided p-value: P(|N(0, sigma^2)| > |S|)
      p <- pnorm(abs(S), mean = 0, sd = S_sd, lower.tail = FALSE) * 2

    } else if (!is.null(df) || !is.null(formals(g)$df)) {

      # GFISHER STATISTICS: Use p.GFisher from GFisher package

      # Get the degrees of freedom from the function g if not provided
      if (!is.null(formals(g)$df)) {
        df <- formals(g)$df
      }

      # Get the p-value type from the function g if available
      if (!is.null(formals(g)$p.type)) {
        p.type <- formals(g)$p.type
      }

      # Force the weights to be non-negative for GFisher
      wts <- pmax(wts, 0)

      # Re-scale the weights to sum to 1
      wts <- wts / sum(abs(wts))

      # Re-calculate the statistic to be consistent with GFisher package
      S <- sum(wts * g(Zscores))

      # Call p.GFisher from GFisher package
      p <- GFisher::p.GFisher(
        q = S,
        df = df,
        w = wts,
        M = M,
        p.type = p.type,
        method = method,
        nsim = nsim,
        seed = seed
      )

    } else {
      stop("The p-value is not available for the given g function.")
    }

    return(list(S = S, p = p))

  } else {
    return(list(S = S))
  }
}


#' Burden test
#'
#' Computes the Burden test statistic and p-value. The Burden test uses a
#' linear combination of Z-scores with specified weights.
#'
#' @description
#' The Burden test statistic is:
#' \deqn{S = \sum_{i=1}^{n} w_i Z_i}
#'
#' Under the null hypothesis, S follows a normal distribution with mean 0 and
#' variance \eqn{\sigma^2 = w^T M w}, where M is the correlation matrix.
#'
#' The p-value is calculated as:
#' \deqn{p = 2 \cdot P(N(0, \sigma^2) > |S|)}
#'
#' This is a two-sided test by definition in genetic association studies.
#'
#' @details
#' **CRITICAL**: This function uses the normal distribution for p-value calculation.
#' It does NOT use Davies method or chi-square approximation. The statistic is
#' a linear combination (NOT squared), and its null distribution is normal.
#'
#' **Algorithm**:
#' 1. Compute S = sum(wts * Zscores)
#' 2. Compute variance: sigma^2 = t(wts) %*% M %*% wts
#' 3. Compute p-value: 2 * P(|N(0, sigma^2)| > |S|)
#'
#' **Weight handling**:
#' - If is.posi.wts=TRUE: negative weights are set to 0
#' - Weights are scaled to sum to 1 (if non-zero)
#'
#' **Computational Complexity**: O(n^2) due to matrix multiplication
#'
#' **Assumptions**:
#' - Zscores follow standard normal under null
#' - M is the correlation matrix of Zscores (diagonal = 1)
#'
#' @param scores Numeric vector of Z-scores or score statistics
#' @param M Correlation matrix of the scores. Default: identity matrix
#' @param wts Numeric vector of weights. Default: equal weights
#' @param calc_p Logical; if TRUE, calculate p-value. Default: TRUE
#' @param is.posi.wts Logical; if TRUE, force weights to be non-negative. Default: TRUE
#'
#' @return A list with elements:
#' \describe{
#'   \item{S}{The Burden test statistic (numeric)}
#'   \item{p}{The two-sided p-value (numeric, only if calc_p=TRUE)}
#' }
#'
#' @examples
#' # Simple example with independent scores
#' Z <- rnorm(10)
#' M <- diag(10)
#' result <- burden_test(Z, M, wts=rep(1,10), calc_p=TRUE)
#'
#' # Example with correlation
#' library(MASS)
#' M <- matrix(0.3, 5, 5) + diag(0.7, 5)
#' Z <- mvrnorm(1, mu=rep(0,5), Sigma=M)
#' result <- burden_test(Z, M, wts=c(1,2,3,2,1), calc_p=TRUE)
#'
#' @noRd
burden_test <- function(scores, M = diag(length(scores)),
                        wts = rep(1, length(scores)),
                        calc_p = TRUE, is.posi.wts = TRUE) {

  # Use calcu_SgZ_p with identity transformation g(x) = x
  # This will route to normal distribution (NOT Davies method)
  calcu_SgZ_p(
    function(x) x,
    Zscores = scores,
    wts = wts,
    calc_p = calc_p,
    M = M,
    is.posi.wts = is.posi.wts
  )
}


#' SKAT test
#'
#' Computes the SKAT (Sequence Kernel Association Test) statistic and p-value.
#' SKAT uses a quadratic form of weighted squared Z-scores.
#'
#' @description
#' The SKAT test statistic is:
#' \deqn{S = \sum_{i=1}^{n} w_i Z_i^2}
#'
#' This is equivalent to a weighted chi-square mixture. The p-value is calculated
#' using the GFisher framework with df=1, which properly accounts for the correlation
#' structure among the Z-scores.
#'
#' @details
#' **CRITICAL**: This function uses p.GFisher with df=1 for p-value calculation.
#' It does NOT use the Liu method or direct Davies method. The GFisher framework
#' with df=1 provides accurate p-values under dependence.
#'
#' **Algorithm**:
#' 1. If scores are not Z-scores (diagonal of M_s != 1), standardize them:
#'    - Adjust weights: wts = wts * diag(M_s)
#'    - Standardize scores: scores = scores / sqrt(diag(M_s))
#'    - Convert M_s to correlation matrix
#' 2. Compute S = sum(wts * scores^2)
#' 3. Call p.GFisher with df=1 to get p-value
#'
#' **Weight handling**:
#' - Negative weights are automatically set to 0 in calcu_SgZ_p
#' - Weights are scaled to sum to 1
#'
#' **Computational Complexity**: O(n^2) for p-value calculation
#'
#' **Assumptions**:
#' - Scores follow standard normal under null (after standardization)
#' - M_s is the covariance/correlation matrix of scores
#'
#' @param scores Numeric vector of Z-scores or score statistics
#' @param M_s Covariance or correlation matrix of the scores. Default: identity matrix
#' @param wts Numeric vector of weights. Default: equal weights
#' @param calc_p Logical; if TRUE, calculate p-value. Default: TRUE
#' @param ... Additional arguments passed to calcu_SgZ_p (e.g., method, nsim)
#'
#' @return A list with elements:
#' \describe{
#'   \item{S}{The SKAT test statistic (numeric)}
#'   \item{p}{The one-sided p-value (numeric, only if calc_p=TRUE)}
#' }
#'
#' @references
#' Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., & Lin, X. (2011).
#' Rare-variant association testing for sequencing data with the sequence kernel
#' association test. The American Journal of Human Genetics, 89(1), 82-93.
#'
#' @examples
#' # Simple example with independent scores
#' Z <- rnorm(10)
#' M <- diag(10)
#' result <- skat_test(Z, M_s=M, wts=rep(1,10), calc_p=TRUE)
#'
#' # Example with correlation and non-standard scores
#' library(MASS)
#' M <- matrix(0.3, 5, 5) + diag(0.7, 5)
#' Z <- mvrnorm(1, mu=rep(0,5), Sigma=M)
#' result <- skat_test(Z, M_s=M, wts=c(1,2,3,2,1), calc_p=TRUE)
#'
#' @noRd
skat_test <- function(scores, M_s = diag(length(scores)),
                      wts = rep(1, length(scores)),
                      calc_p = TRUE, ...) {

  # If scores are not Z-scores (diagonal of M_s is not all 1's),
  # adjust weights and standardize scores
  if (!identical(diag(M_s), rep(1, length(scores)))) {
    # Get standard deviations from diagonal of M_s
    SDs <- sqrt(diag(M_s))

    # Adjust weights: multiply by variance (SD^2)
    wts <- wts * SDs^2

    # Standardize scores to Z-scores
    scores <- scores / SDs

    # Convert covariance matrix to correlation matrix
    M_s <- cov2cor(M_s)
  }

  # Use calcu_SgZ_p with quadratic transformation g(x) = x^2
  # This will route to p.GFisher with df=1 (NOT Liu method)
  calcu_SgZ_p(
    g = function(x) x^2,
    Zscores = scores,
    wts = wts,
    calc_p = calc_p,
    M = M_s,
    df = 1,
    ...
  )
}


#' Fisher test based on Z-scores
#'
#' Computes the Fisher combination test statistic and p-value using Z-scores.
#' This uses the GFisher transformation with df=2.
#'
#' @description
#' The Fisher test statistic is computed via the transformation:
#' \deqn{S = \sum_{i=1}^{n} w_i g(Z_i)}
#'
#' where \eqn{g(Z_i)} transforms Z-scores to chi-square statistics using the
#' GFisher framework. For two-sided tests:
#' \deqn{g(Z_i) = F_{\chi^2_2}^{-1}(\log(2) + \log(\bar{\Phi}(|Z_i|)))}
#'
#' For one-sided tests, the log(2) term is omitted.
#'
#' @details
#' **CRITICAL**: This function uses p.GFisher with df=2 for p-value calculation.
#' The transformation g_GFisher handles the conversion from Z-scores to the
#' appropriate chi-square scale.
#'
#' **Algorithm**:
#' 1. Define transformation function g based on p.type (one-sided or two-sided)
#' 2. Compute S = sum(wts * g(Zscores))
#' 3. Call p.GFisher with df=2 to get p-value
#'
#' **Weight handling**:
#' - Negative weights are automatically set to 0 in calcu_SgZ_p
#' - Weights are scaled to sum to 1
#'
#' **Computational Complexity**: O(n^2) for p-value calculation
#'
#' **Assumptions**:
#' - Zscores follow standard normal under null
#' - M is the correlation matrix of Zscores
#'
#' @param Zscores Numeric vector of Z-scores (standard normal under null)
#' @param M Correlation matrix of the Z-scores. Default: identity matrix
#' @param wts Numeric vector of weights. Default: equal weights
#' @param calc_p Logical; if TRUE, calculate p-value. Default: TRUE
#' @param p.type Character string: "two" for two-sided (default), "one" for one-sided
#' @param ... Additional arguments passed to calcu_SgZ_p (e.g., method, nsim)
#'
#' @return A list with elements:
#' \describe{
#'   \item{S}{The Fisher test statistic (numeric)}
#'   \item{p}{The p-value (numeric, only if calc_p=TRUE)}
#' }
#'
#' @references
#' Zhang, H., & Wu, Z. (2023). The generalized Fisher's combination and
#' accurate p-value calculation under dependence. Biometrics, 79(2), 1159-1172.
#'
#' @examples
#' # Two-sided Fisher test
#' Z <- rnorm(10)
#' M <- diag(10)
#' result <- fisher_test_Z(Z, M, wts=rep(1,10), calc_p=TRUE, p.type="two")
#'
#' # One-sided Fisher test
#' result <- fisher_test_Z(Z, M, wts=rep(1,10), calc_p=TRUE, p.type="one")
#'
#' @noRd
fisher_test_Z <- function(Zscores, M = diag(length(Zscores)),
                          wts = rep(1, length(Zscores)),
                          calc_p = TRUE, p.type = "two", ...) {

  # Define the transformation function based on p.type
  if (p.type == "two") {
    # Two-sided Fisher test: g_fisher(x) uses absolute value
    g_fisher <- function(x, df = 2, p.type = "two") g_GFisher(x, df, p.type)
  } else if (p.type == "one") {
    # One-sided Fisher test: g_fisher(x) preserves sign
    g_fisher <- function(x, df = 2, p.type = "one") g_GFisher(x, df, p.type)
  }

  # Use calcu_SgZ_p with GFisher transformation
  # This will route to p.GFisher with df=2
  calcu_SgZ_p(
    g = g_fisher,
    Zscores = Zscores,
    wts = wts,
    calc_p = calc_p,
    M = M,
    df = 2,
    p.type = p.type,
    ...
  )
}


#' Multiple S(g(Z)) tests with various weights and degrees of freedom
#'
#' Performs multiple variant-set tests with different weight schemes and
#' transformation functions. This is used to compute Burden, SKAT, and Fisher
#' tests with various optimal and equal weight configurations.
#'
#' @description
#' Computes multiple tests of the form \eqn{S = \sum w_i g(Z_i)} where different
#' rows of the weight matrix W and different degrees of freedom in DF define
#' different tests. Each test may use:
#' - Burden (df=Inf): Linear combination, normal distribution
#' - SKAT (df=1): Quadratic form, p.GFisher with df=1
#' - Fisher (df=2): GFisher transformation, p.GFisher with df=2
#'
#' @details
#' **Algorithm**:
#' For each test (row i in DF and W matrices):
#' 1. If DF[i] = Inf: Call burden_test
#' 2. Else: Create g function with g_GFisher and call calcu_SgZ_p
#' 3. Store statistic and p-value
#'
#' Row names are automatically generated as "df_[df]_[w_name]".
#'
#' **Computational Complexity**: O(k * n^2) where k is number of tests
#'
#' **Assumptions**:
#' - Zscores follow standard normal under null
#' - M is the correlation matrix of Zscores
#'
#' @param Zscores Numeric vector of Z-scores
#' @param DF Matrix or vector of degrees of freedom (one per test). Use Inf for Burden
#' @param W Matrix of weights (one row per test, columns correspond to variants)
#' @param M Correlation matrix of the Z-scores
#' @param p.type Character string: "two" for two-sided, "one" for one-sided
#' @param calcu_p Logical; if TRUE, calculate p-values. Default: TRUE
#' @param is.posi.wts Logical; if TRUE, force weights to be non-negative. Default: TRUE
#' @param w_names Optional character vector of weight scheme names. Default: rownames(W)
#' @param ... Additional arguments passed to burden_test or calcu_SgZ_p
#'
#' @return A list with elements:
#' \describe{
#'   \item{STAT}{Matrix of test statistics (one row per test)}
#'   \item{PVAL}{Matrix of p-values (one row per test, only if calcu_p=TRUE)}
#' }
#'
#' @examples
#' # Multiple tests with different df and weights
#' Z <- rnorm(10)
#' M <- diag(10)
#' DF <- matrix(c(Inf, 1, 2), ncol=1)  # Burden, SKAT, Fisher
#' W <- rbind(rep(1,10), rep(1,10), rep(1,10))  # Equal weights
#' result <- multi_SgZ_test(Z, DF, W, M, p.type="two", calcu_p=TRUE)
#'
#' @noRd
multi_SgZ_test <- function(Zscores, DF, W, M, p.type = "two", calcu_p = TRUE,
                           is.posi.wts = TRUE, w_names = rownames(W), ...) {

  # Number of tests
  testN <- nrow(DF)

  # Initialize result matrices
  STAT <- PVAL <- matrix(rep(NA, testN), ncol = 1)
  rownames(STAT) <- rownames(PVAL) <- paste("df", as.character(DF), w_names, sep = "_")

  # Loop through each test
  for (i in 1:testN) {

    if (DF[i, 1] == Inf) {
      # Burden test: df=Inf, use normal distribution
      result <- burden_test(
        scores = Zscores,
        M = M,
        wts = W[i, ],
        calc_p = calcu_p,
        is.posi.wts = is.posi.wts
      )

    } else {
      # GFisher-based tests (SKAT, Fisher, etc.)
      # Create transformation function with appropriate df and p.type
      g <- function(x) g_GFisher(x, df = DF[i, ], p.type = p.type)

      result <- calcu_SgZ_p(
        g = g,
        Zscores = Zscores,
        wts = W[i, ],
        calc_p = calcu_p,
        M = M,
        df = DF[i, ],
        p.type = p.type,
        is.posi.wts = is.posi.wts,
        ...
      )
    }

    # Store results
    if (calcu_p) {
      STAT[i, 1] <- result$S
      PVAL[i, 1] <- result$p
    } else {
      STAT[i, 1] <- result$S
    }
  }

  return(list(STAT = STAT, PVAL = PVAL))
}


#' Omnibus test combining multiple S(g(Z)) statistics via CCT
#'
#' Computes an omnibus test by combining multiple variant-set tests (Burden, SKAT,
#' Fisher) using the Cauchy Combination Test (CCT). This provides a robust test
#' that adapts to different genetic architectures.
#'
#' @description
#' The omnibus test performs the following:
#' 1. Compute multiple tests via multi_SgZ_test (Burden, SKAT, Fisher with various weights)
#' 2. Combine their p-values using CCT (Cauchy Combination Test)
#' 3. Return individual test results and the combined CCT statistic/p-value
#'
#' The CCT combination is robust to correlation among the individual tests and
#' provides a powerful omnibus test without requiring knowledge of the genetic
#' architecture.
#'
#' @details
#' **Algorithm**:
#' 1. Call multi_SgZ_test to get statistics and p-values for all tests
#' 2. If calcu_p=TRUE:
#'    - Call cct_test to combine p-values
#'    - Append CCT statistic and p-value to results
#' 3. Return combined results
#'
#' **Computational Complexity**: O(k * n^2) where k is number of tests
#'
#' **Assumptions**:
#' - Individual tests are valid
#' - CCT maintains type I error control under dependence
#'
#' @param Zscores Numeric vector of Z-scores
#' @param DF Matrix or vector of degrees of freedom (one per test). Use Inf for Burden
#' @param W Matrix of weights (one row per test, columns correspond to variants)
#' @param M Correlation matrix of the Z-scores
#' @param p.type Character string: "two" for two-sided, "one" for one-sided
#' @param calcu_p Logical; if TRUE, calculate p-values. Default: TRUE
#' @param is.posi.wts Logical; if TRUE, force weights to be non-negative. Default: TRUE
#' @param ... Additional arguments passed to multi_SgZ_test
#'
#' @return A list with elements:
#' \describe{
#'   \item{STAT}{Matrix of test statistics (individual tests + CCT)}
#'   \item{PVAL}{Matrix of p-values (individual tests + CCT p-value)}
#'   \item{cct}{The CCT statistic (numeric)}
#'   \item{pval_cct}{The CCT p-value (numeric)}
#' }
#'
#' @references
#' Liu Y, Xie J (2020). "Cauchy Combination Test: A Powerful Test With Analytic
#' p-Value Calculation Under Arbitrary Dependency Structures."
#' JASA, 115(529), 393-402.
#'
#' @examples
#' # Omnibus test combining Burden, SKAT, and Fisher
#' Z <- rnorm(10)
#' M <- diag(10)
#' DF <- matrix(c(Inf, 1, 2), ncol=1)  # Burden, SKAT, Fisher
#' W <- rbind(rep(1,10), rep(1,10), rep(1,10))  # Equal weights
#' result <- omni_SgZ_test(Z, DF, W, M, p.type="two", calcu_p=TRUE)
#' result$pval_cct  # Combined p-value
#'
#' @noRd
omni_SgZ_test <- function(Zscores, DF, W, M, p.type = "two", calcu_p = TRUE,
                          is.posi.wts = TRUE, ...) {

  # Compute multiple tests
  multi_tests <- multi_SgZ_test(
    Zscores = Zscores,
    DF = DF,
    W = W,
    M = M,
    p.type = p.type,
    calcu_p = calcu_p,
    is.posi.wts = is.posi.wts,
    ...
  )

  # Combine p-values using CCT if requested
  if (calcu_p) {
    cct_result <- cct_test(multi_tests$PVAL)
    cct <- cct_result$cct
    pval_cct <- cct_result$pval_cct
  } else {
    cct <- NA
    pval_cct <- NA
  }

  # Return individual test results plus CCT
  return(list(
    STAT = rbind(multi_tests$STAT, cct),
    PVAL = rbind(multi_tests$PVAL, pval_cct),
    cct = cct,
    pval_cct = pval_cct
  ))
}
