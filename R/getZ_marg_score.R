# ==============================================================================
# Score Calculation Functions for GLOW Tests
# ==============================================================================
#
# This file contains functions for computing standardized marginal score
# statistics (Z-scores) and their correlation structure under generalized
# linear models for both continuous and binary traits.
#
# EXPORTED FUNCTIONS:
#   - getZ_marg_score()          Main score calculation (continuous/binary traits)
#   - getZ_marg_score_binary_SPA()  SPA-adjusted scores for binary traits
#   - fit_null_model()           Pre-compute null model for repeated use (optimization)
#
# INTERNAL HELPERS:
#   - .validate_score_inputs()   Input validation
#   - .compute_projection_matrices()  H^(1/2) matrix computation
#   - .handle_zero_variance_snps()   Zero-variance SNP detection


#################### EXPORTED MAIN FUNCTIONS ####################

#' Get Standardized Marginal Score Statistics
#'
#' @useDynLib GLOWr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @description
#' Calculates the standardized marginal score statistics for each variant under
#' a generalized linear model framework. This is the foundation function for all
#' GLOW tests, providing the score statistics and their correlation structure.
#'
#' For continuous traits, there is an option to use t-statistics from the \code{lm()}
#' function of Y on each column of G and X. This t-statistic is close to the
#' null-model defined marginal score statistics but not exactly the same, partially
#' due to different ways of estimating the error variance.
#'
#' @details
#' The function implements the following algorithm:
#'
#' \strong{For Binary Traits (Logistic Regression):}
#' \enumerate{
#'   \item Fit null model: \eqn{\text{logit}(P(Y=1)) = X\beta}
#'   \item Compute fitted probabilities: \eqn{Y_0 = \hat{P}(Y=1|X)}
#'   \item Define weighted design matrix: \eqn{\tilde{X} = \sqrt{Y_0(1-Y_0)} \cdot X}
#'   \item Compute projection: \eqn{H^{1/2} = \tilde{X}({\tilde{X}}'\tilde{X})^{-1/2}}
#'   \item Define weighted genotype matrix: \eqn{\tilde{G} = \sqrt{Y_0(1-Y_0)} \cdot G}
#'   \item Compute score variance: \eqn{M_s = {\tilde{G}}'{\tilde{G}} - {\tilde{G}}'H^{1/2}{H^{1/2}}'{\tilde{G}}}
#'   \item Compute score statistics: \eqn{S = G'(Y - Y_0)}
#'   \item Compute Z-scores: \eqn{Z_j = S_j / \sqrt{M_{s,jj}}}
#'   \item Set dispersion: \eqn{s_0 = 1}
#' }
#'
#' \strong{For Continuous Traits (Linear Regression):}
#' \enumerate{
#'   \item Fit null model: \eqn{Y = X\beta + \epsilon}
#'   \item Compute projection: \eqn{H^{1/2} = X(X'X)^{-1/2}}
#'   \item Compute score variance: \eqn{M_s = G'G - G'H^{1/2}{H^{1/2}}'G}
#'   \item Obtain residuals: \eqn{r = Y - X\hat{\beta}}
#'   \item Estimate dispersion: \eqn{s_0 = \text{SD}(r)}
#'   \item Compute score statistics: \eqn{S = G'r / s_0}
#'   \item If \code{use_lm_t = TRUE}: Extract t-statistics from marginal regressions
#'   \item If \code{use_lm_t = FALSE}: Compute Z-scores as \eqn{Z_j = S_j / \sqrt{M_{s,jj}}}
#' }
#'
#' The correlation matrix \eqn{M_Z} is obtained by standardizing \eqn{M_s}:
#' \deqn{M_Z = \text{cor}(M_s)}
#'
#' \strong{Mathematical Properties:}
#' \itemize{
#'   \item The Z-scores are asymptotically standard normal under the null hypothesis
#'   \item \eqn{M_Z} represents the correlation of Z-scores when effects are fixed (including H0)
#'   \item \eqn{M_s} is the covariance matrix of score statistics under fixed effects
#'   \item The function properly accounts for correlation induced by covariates through projection
#' }
#'
#' \strong{Computational Complexity:}
#' \itemize{
#'   \item Time: \eqn{O(n \cdot p \cdot k + p \cdot k^2)} where n = sample size, p = number of SNPs, k = number of covariates
#'   \item Space: \eqn{O(p^2 + n \cdot p + n \cdot k)}
#' }
#'
#' @param G A numeric matrix of genotypes (n x p), where n is the number of individuals
#'   and p is the number of SNPs. Missing values (NA) are not allowed.
#' @param X A numeric matrix of covariates (n x k), where k is the number of covariates.
#'   Required if \code{null_model = NULL}. Ignored if \code{null_model} is provided.
#' @param Y A numeric vector of response variable (length n). For binary traits, must be
#'   coded as 0/1. For continuous traits, any numeric values are allowed.
#'   Required if \code{null_model = NULL}. Ignored if \code{null_model} is provided.
#' @param trait Character string indicating trait type. Must be either "binary" (logistic
#'   regression) or "continuous" (linear regression). Default is "binary".
#'   Ignored if \code{null_model} is provided (trait is extracted from null_model).
#' @param null_model Optional. A null_model object from \code{fit_null_model()}.
#'   If provided, the function reuses the pre-computed null model, which can improve
#'   performance for repeated calls with different SNP sets. When \code{null_model} is
#'   provided, \code{X}, \code{Y}, and \code{trait} are ignored. Default is \code{NULL}
#'   (standard mode - fit null model from X and Y).
#' @param use_lm_t Logical. For continuous traits, whether to use t-statistics from
#'   \code{lm()} function as Z-scores. If \code{FALSE} (default), uses null-model
#'   defined Z-scores. Ignored for binary traits.
#' @param use_cpp Logical. Whether to use C++ implementation for GHG matrix computation.
#'   Default is \code{FALSE} (use R). The C++ version may be faster on some platforms
#'   (Linux with standard BLAS), but is not necessarily so on Apple Silicon due to highly
#'   optimized ARM64 BLAS. Set to \code{TRUE} to force C++ usage for portability testing.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Zscores}{Numeric vector (length p) of standardized marginal score statistics
#'     (Z-scores) for each variant}
#'   \item{scores}{Numeric vector (length p) of marginal score statistics for each variant.
#'     For the standard method, \eqn{Z_j = S_j / (\sqrt{M_{s,jj}} \cdot s_0)} (exact
#'     relationship). Compare with the SPA method where Z-scores are derived from
#'     SPA-adjusted p-values.}
#'   \item{M_Z}{Numeric matrix (p x p) giving the correlation matrix of Z-scores when
#'     effects are fixed (including H0), given X and G fixed}
#'   \item{M_s}{Numeric matrix (p x p) giving the covariance matrix of score statistics
#'     when effects are fixed (including H0), given X and G fixed}
#'   \item{s0}{Numeric scalar giving the estimated dispersion parameter. For binary traits,
#'     s0 = 1. For continuous traits, s0 is the residual standard deviation under the
#'     null model}
#' }
#'
#' @examples
#' # Example 1: Continuous trait with null-model Z-scores (default)
#' set.seed(123)
#' n <- 500
#' p <- 10
#' # Generate genotypes with varying MAF
#' G <- matrix(0, n, p)
#' for (j in 1:p) {
#'   maf <- 0.05 * j / p  # MAF from 0.005 to 0.05
#'   G[, j] <- rbinom(n, 2, maf)
#' }
#' # Covariates: intercept will be added automatically
#' X <- matrix(rnorm(n * 2), n, 2)
#' Y <- rnorm(n)  # Continuous trait
#' result_cont <- getZ_marg_score(G, X, Y, trait = "continuous")
#' head(result_cont$Zscores)
#'
#' # Example 2: Continuous trait with lm() t-statistics
#' result_lm_t <- getZ_marg_score(G, X, Y, trait = "continuous", use_lm_t = TRUE)
#' head(result_lm_t$Zscores)
#'
#' # Example 3: Binary trait
#' Y_binary <- sample(c(0, 1), n, replace = TRUE)
#' result_binary <- getZ_marg_score(G, X, Y_binary, trait = "binary")
#' head(result_binary$Zscores)
#' result_binary$s0  # Should be 1 for binary traits
#'
#' # Example 4: Examine correlation structure
#' result <- getZ_marg_score(G, X, Y, trait = "continuous")
#' image(result$M_Z, main = "Z-score Correlation Matrix")
#'
#' # Example 5: Optimized workflow for multiple SNP sets (RECOMMENDED FOR GWAS)
#' # Create multiple gene sets
#' gene_sets <- lapply(1:100, function(i) {
#'   matrix(rbinom(n * 5, 2, runif(1, 0.01, 0.1)), n, 5)
#' })
#'
#' # Step 1: Fit null model once
#' null_obj <- fit_null_model(X, Y_binary, trait = "binary")
#'
#' # Step 2: Reuse for all gene sets (much faster!)
#' results_optimized <- lapply(gene_sets, function(G_gene) {
#'   getZ_marg_score(G_gene, null_model = null_obj)
#' })
#'
#' # This is ~10-100x faster than refitting null model for each gene set
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @export
getZ_marg_score <- function(G, X = NULL, Y = NULL, trait = "binary", null_model = NULL,
                            use_lm_t = FALSE, use_cpp = FALSE) {
  # ========== Input Mode Selection ==========
  # Two modes: (1) Standard - provide X and Y, (2) Optimized - provide null_model

  if (is.null(null_model)) {
    # ========== Mode 1: Standard (fit null model) ==========
    if (is.null(X) || is.null(Y)) {
      stop("Either provide (X, Y, trait) or null_model. ",
           "For standard mode: provide X and Y. ",
           "For optimized mode: provide null_model from fit_null_model().")
    }

    # Validate inputs using helper function
    validated <- .validate_score_inputs(G, X, Y, trait = trait)
    G <- validated$G
    X <- validated$X
    Y <- validated$Y

    # If the first column of X is not all 1's, we need to add it for intercept.
    if(ncol(X) == 0 || !all(X[,1] == 1)){
      X = cbind(1,X)
    }

  } else {
    # ========== Mode 2: Optimized (reuse null_model) ==========
    if (!inherits(null_model, "glow_null_model")) {
      stop("null_model must be created by fit_null_model()")
    }

    if (!is.null(X) || !is.null(Y)) {
      warning("null_model provided; ignoring X and Y parameters")
    }

    # Validate G dimensions only
    if (!is.matrix(G)) {
      stop("G must be a matrix")
    }
    if (any(is.na(G))) {
      stop("G contains missing values (NA). Please impute or remove before analysis")
    }
    if (nrow(G) != null_model$n) {
      stop("G has ", nrow(G), " rows but null_model was fit with ",
           null_model$n, " samples. Dimensions must match.")
    }

    # Extract components from null_model
    trait <- null_model$trait
    X <- null_model$X

    # Extract trait-specific components
    if (trait == "binary") {
      Y0 <- null_model$fitted_probs
      Y <- null_model$Y
    } else if (trait == "continuous") {
      res <- null_model$residuals
      s0 <- null_model$s0
      Y <- null_model$Y  # Needed for use_lm_t option
    }
  }

  # ========== Main Computation (same for both modes) ==========

  # Binary trait (logistic regression)
  if (trait == "binary") {
    # Get fitted probabilities (from null_model or fit new)
    if (is.null(null_model)) {
      # Mode 1: Fit null model
      mod0 <- glm(Y ~ X, family = "binomial")
      Y0 <- mod0$fitted.values
    }
    # Mode 2: Y0 already extracted from null_model above

    # Compute weighted matrices using helper function
    w <- sqrt(Y0 * (1 - Y0))
    proj_matrices <- .compute_projection_matrices(X, G, weights = w)
    Hhalf <- proj_matrices$Hhalf
    Gtilde <- proj_matrices$Gtilde

    # Compute score variance matrix: M_s = Gtilde'Gtilde - Gtilde'H^(1/2)(H^(1/2))'Gtilde
    # Option to use C++ or R implementation
    if (use_cpp) {
      # C++ implementation (may be faster on Linux/standard BLAS)
      ghg_result <- compute_GHG_cpp(Gtilde, Hhalf)
      GHG <- ghg_result$GHG
      GHhalf <- ghg_result$GHhalf
    } else {
      # R implementation (faster on Apple Silicon with optimized BLAS)
      GHhalf <- t(Gtilde) %*% Hhalf
      GHG <- t(Gtilde) %*% Gtilde - GHhalf %*% t(GHhalf)
    }

    # Score statistics: S = G' * (Y - Y0)
    score <- t(G) %*% (Y - Y0)

    # Handle zero-variance SNPs using helper function
    diag_GHG <- diag(GHG)
    zero_var <- .handle_zero_variance_snps(diag_GHG)
    if (zero_var$has_zero_var) {
      warning("Detected ", length(zero_var$zero_var_idx), " SNP(s) with zero variance ",
              "(constant genotype). Z-scores set to 0 for these variants.")
    }

    # Z-scores: Z_j = S_j / sqrt(M_s[j,j])
    Zscore <- score / sqrt(zero_var$diag_GHG_safe)
    # Set Z-scores to 0 for zero-variance SNPs
    if (zero_var$has_zero_var) {
      Zscore[zero_var$zero_var_idx] <- 0
    }

    # Dispersion parameter for binary trait is 1
    s0 <- 1

  } else if (trait == "continuous") {
    # Continuous trait (linear regression)

    # Compute projection matrices using helper function (no weights)
    proj_matrices <- .compute_projection_matrices(X, G, weights = NULL)
    Hhalf <- proj_matrices$Hhalf
    Gtilde <- proj_matrices$Gtilde  # Same as G for unweighted case

    # Compute score variance matrix: M_s = G'G - G'H^(1/2)(H^(1/2))'G
    # Option to use C++ or R implementation
    if (use_cpp) {
      # C++ implementation (may be faster on Linux/standard BLAS)
      ghg_result <- compute_GHG_cpp(G, Hhalf)
      GHG <- ghg_result$GHG
      GHhalf <- ghg_result$GHhalf
    } else {
      # R implementation (faster on Apple Silicon with optimized BLAS)
      GHhalf <- t(G) %*% Hhalf
      GHG <- t(G) %*% G - GHhalf %*% t(GHhalf)
    }

    # Get residuals and dispersion (from null_model or fit new)
    if (is.null(null_model)) {
      # Mode 1: Fit null model
      mod0 <- glm(Y ~ X, family = "gaussian")
      res <- residuals(mod0)
      s0 <- sd(res)
    }
    # Mode 2: res and s0 already extracted from null_model above

    # Score statistics: S = G' * res / s0
    score <- t(G) %*% res / s0

    # Handle zero-variance SNPs using helper function
    diag_GHG <- diag(GHG)
    zero_var <- .handle_zero_variance_snps(diag_GHG)
    if (zero_var$has_zero_var) {
      warning("Detected ", length(zero_var$zero_var_idx), " SNP(s) with zero variance ",
              "(constant genotype). Z-scores set to 0 for these variants.")
    }

    if (use_lm_t) {
      # Option: Use marginal t-statistics from lm(Y ~ g + X) for each SNP
      # This extracts the t-statistic for the coefficient of g (position 2 in coefficients)
      Zscore <- apply(G, 2, function(g) {
        summary(lm(Y ~ g + X))$coefficients[2, 3]
      })
    } else {
      # Default: Use null-model defined Z-scores
      # Z_j = S_j / sqrt(M_s[j,j])
      Zscore <- score / sqrt(zero_var$diag_GHG_safe)
      # Set Z-scores to 0 for zero-variance SNPs
      if (zero_var$has_zero_var) {
        Zscore[zero_var$zero_var_idx] <- 0
      }
    }
  }

  # Compute correlation matrix of Z-scores from covariance matrix of scores
  # M_Z = cor(M_s). A zero-variance (constant genotype) column has a zero
  # diagonal in GHG, which makes cov2cor() emit a "non-positive diagonal"
  # warning and return NaN for that row/column. Such columns are normally
  # filtered upstream (min_mac >= 1); guard here so a stray constant column
  # still yields a finite, valid correlation matrix (treated as uncorrelated,
  # diagonal 1) and only the intended "zero variance" warning is raised.
  if (zero_var$has_zero_var) {
    GHG_cor <- GHG
    diag(GHG_cor)[zero_var$zero_var_idx] <- 1  # avoid div-by-zero in cov2cor
    M <- cov2cor(GHG_cor)
    M[zero_var$zero_var_idx, ] <- 0
    M[, zero_var$zero_var_idx] <- 0
    diag(M)[zero_var$zero_var_idx] <- 1
  } else {
    M <- cov2cor(GHG)
  }

  # Return results
  return(list(
    Zscores = as.vector(Zscore),
    scores = score,
    M_Z = M,
    M_s = GHG,
    s0 = s0
  ))
}


#' Get SPA-Adjusted Standardized Marginal Score Statistics for Binary Traits
#'
#' @description
#' Calculates the saddlepoint approximation (SPA)-adjusted standardized marginal
#' score statistics for each variant under a logistic regression model for binary
#' traits. This function is particularly important for handling rare variants and
#' unbalanced case-control data where standard asymptotic p-values may be inflated.
#'
#' The SPA method provides more accurate p-values than the standard normal
#' approximation, especially when:
#' \itemize{
#'   \item Minor allele frequencies are very low (< 1%)
#'   \item Case-control ratios are highly unbalanced
#'   \item Sample sizes are moderate
#' }
#'
#' @details
#' \strong{Algorithm Overview:}
#'
#' The function implements the following steps:
#'
#' \enumerate{
#'   \item \strong{Null Model Fitting:}
#'   Fit logistic regression under the null hypothesis (no genetic effect):
#'   \deqn{\text{logit}(P(Y=1)) = X\beta}
#'   where \eqn{X} is the covariate matrix. Obtain fitted probabilities \eqn{Y_0 = P(Y=1|X)}.
#'
#'   \item \strong{Residual Genotype Calculation:}
#'   Compute genotypes adjusted for covariates:
#'   \deqn{G_{\text{scale}} = G - X(X'VX)^{-1}X'VG}
#'   where \eqn{V = \text{diag}(Y_0(1-Y_0))} is the variance matrix.
#'
#'   \item \strong{Score Statistics:}
#'   Calculate score statistics for each variant:
#'   \deqn{S_j = G_{\text{scale},j}'Y}
#'   where \eqn{G_{\text{scale},j}} is the j-th column of \eqn{G_{\text{scale}}}.
#'
#'   \item \strong{SPA-Adjusted P-values:}
#'   For each variant, apply the saddlepoint approximation to obtain p-values
#'   that account for the discrete nature of genotypes and the skewness of the
#'   score statistic distribution. Uses the "fastSPA" method from the SPAtest
#'   package, which employs a partially normal approximation for improved
#'   efficiency with rare variants.
#'
#'   \item \strong{SPA-Adjusted Z-scores:}
#'   Convert SPA p-values to Z-scores while preserving the direction of effect:
#'   \deqn{Z_{\text{SPA},j} = \Phi^{-1}(1 - p_{\text{SPA},j}/2) \cdot \text{sign}(S_j)}
#'   where \eqn{\Phi^{-1}} is the inverse standard normal CDF.
#'
#'   \item \strong{Correlation Structure:}
#'   Compute the correlation matrix \eqn{M_Z} and covariance matrix \eqn{M_s} using
#'   weighted genotypes:
#'   \deqn{\tilde{G} = \sqrt{Y_0(1-Y_0)} \cdot G}
#'   \deqn{M_s = \tilde{G}'(\tilde{G} - H^{1/2}{H^{1/2}}'\tilde{G})}
#'   where \eqn{H^{1/2} = \tilde{X}(\tilde{X}'\tilde{X})^{-1/2}} with
#'   \eqn{\tilde{X} = \sqrt{Y_0(1-Y_0)} \cdot X}.
#' }
#'
#' \strong{Saddlepoint Approximation (SPA):}
#'
#' The SPA method provides a more accurate approximation to the null distribution
#' of the score statistic by using the cumulant generating function (CGF). For a
#' score statistic \eqn{S}, the SPA approximates the tail probability as:
#'
#' \deqn{P(S \geq s) \approx \Phi(\hat{w} + \frac{1}{\hat{w}}\log(\frac{\hat{u}}{\hat{w}}))}
#'
#' where \eqn{\hat{w}} and \eqn{\hat{u}} are derived from the CGF and its derivatives
#' evaluated at the saddlepoint. The "fastSPA" method uses a hybrid approach:
#' \itemize{
#'   \item For common variants: uses standard normal approximation
#'   \item For rare variants: uses full SPA correction
#'   \item Cutoff determined by "BE" (Bonferroni-equivalent) threshold
#' }
#'
#' \strong{Mathematical Properties:}
#' \itemize{
#'   \item SPA-adjusted p-values are more accurate than asymptotic p-values for rare variants
#'   \item The correlation structure \eqn{M_Z} is identical to the non-SPA version
#'   \item Only the marginal Z-scores are adjusted; the correlation remains based on variance
#'   \item Dispersion parameter \eqn{s_0 = 1} for binary traits (canonical link)
#' }
#'
#' \strong{Computational Complexity:}
#' \itemize{
#'   \item Time: \eqn{O(n \cdot p \cdot k + p \cdot k^2 + p \cdot t_{\text{SPA}})}
#'     where \eqn{t_{\text{SPA}}} is the SPA computation time per variant
#'   \item Space: \eqn{O(p^2 + n \cdot p + n \cdot k)}
#'   \item The "fastSPA" method is optimized for large-scale GWAS with many rare variants
#' }
#'
#' \strong{When to Use SPA vs. Standard Method:}
#' \itemize{
#'   \item Use SPA when: MAF < 1%, case-control ratio < 0.1 or > 10, or moderate sample size (n < 5000)
#'   \item Standard method acceptable when: MAF > 5%, balanced design, large sample size (n > 10000)
#'   \item SPA adds computational cost but provides Type I error control for rare variants
#' }
#'
#' @param G A numeric matrix of genotypes (n x p), where n is the number of individuals
#'   and p is the number of SNPs. Each column represents a variant coded as 0, 1, or 2
#'   (number of minor alleles). Missing values (NA) are not allowed.
#' @param X A numeric matrix of covariates (n x k), where k is the number of covariates.
#'   Required if \code{null_model = NULL}. Ignored if \code{null_model} is provided.
#'   Should NOT include an intercept column; the function adds it automatically.
#' @param Y A numeric vector of binary response variable (length n). Must be coded as
#'   0 (controls) or 1 (cases). No missing values allowed.
#'   Required if \code{null_model = NULL}. Ignored if \code{null_model} is provided.
#' @param null_model Optional. A binary trait null_model object from \code{fit_null_model()}.
#'   If provided, the function reuses the pre-computed null model, dramatically improving
#'   performance for repeated calls with different SNP sets. When \code{null_model} is
#'   provided, \code{X} and \code{Y} are ignored. Default is \code{NULL}
#'   (standard mode - fit null model from X and Y).
#'   **Note**: Must be a binary trait null_model (from \code{fit_null_model(..., trait = "binary")}).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Zscores}{Numeric vector (length p) of SPA-adjusted standardized marginal
#'     score statistics (Z-scores) for each variant. These account for the discrete
#'     nature of genotypes and provide better Type I error control for rare variants.}
#'   \item{scores}{Numeric vector (length p) of marginal score statistics
#'     \eqn{S_j = G_j'(Y - Y_0)}. These are the raw (unstandardized) association
#'     signals. Note: unlike the standard method where
#'     \eqn{Z_j = S_j / \sqrt{M_{s,jj}}}, the SPA Z-scores are derived from
#'     SPA-adjusted p-values: \eqn{Z_{SPA,j} = \Phi^{-1}(1 - p_{SPA,j}/2) \cdot sign(S_j)}.
#'     The scores are returned for completeness and diagnostic purposes.}
#'   \item{M_Z}{Numeric matrix (p x p) giving the correlation matrix of Z-scores when
#'     effects are fixed (including H0), given X and G fixed. Identical to non-SPA version.}
#'   \item{M_s}{Numeric matrix (p x p) giving the covariance matrix of score statistics
#'     when effects are fixed (including H0), given X and G fixed.}
#'   \item{s0}{Numeric scalar equal to 1 (dispersion parameter for binary traits with
#'     canonical logistic link).}
#' }
#'
#' @examples
#' # Example 1: Balanced case-control with common and rare variants
#' set.seed(123)
#' n <- 500
#' p <- 10
#'
#' # Generate genotypes with varying MAF
#' G <- matrix(0, n, p)
#' for (j in 1:p) {
#'   maf <- 0.05 * j / p  # MAF from 0.005 to 0.05
#'   G[, j] <- rbinom(n, 2, maf)
#' }
#'
#' # Covariates: intercept will be added automatically
#' X <- matrix(rnorm(n * 2), n, 2)
#'
#' # Binary trait (balanced)
#' Y <- sample(c(0, 1), n, replace = TRUE)
#'
#' # SPA-adjusted analysis
#' result_spa <- getZ_marg_score_binary_SPA(G, X, Y)
#' head(result_spa$Zscores)
#'
#' # Compare with standard method (for common variants, should be similar)
#' result_std <- getZ_marg_score(G, cbind(1, X), Y, trait = "binary")
#' plot(result_std$Zscores, result_spa$Zscores,
#'      xlab = "Standard Z-scores", ylab = "SPA-adjusted Z-scores",
#'      main = "Comparison of Methods")
#' abline(0, 1, col = "red")
#'
#' # Example 2: Unbalanced case-control with very rare variants
#' set.seed(456)
#' n_case <- 100
#' n_ctrl <- 900
#' n <- n_case + n_ctrl
#' p <- 5
#'
#' # Very rare variants (MAF < 0.01)
#' G <- matrix(0, n, p)
#' for (j in 1:p) {
#'   maf <- 0.001 * j  # MAF from 0.001 to 0.005
#'   G[, j] <- rbinom(n, 2, maf)
#' }
#'
#' # Unbalanced binary trait
#' Y <- c(rep(1, n_case), rep(0, n_ctrl))
#'
#' # Covariates
#' X <- matrix(rnorm(n * 3), n, 3)
#'
#' # SPA-adjusted analysis (critical for this scenario)
#' result_rare <- getZ_marg_score_binary_SPA(G, X, Y)
#' print(result_rare$Zscores)
#'
#' # Examine correlation structure
#' image(result_rare$M_Z, main = "Z-score Correlation Matrix")
#'
#' # Example 3: Optimized workflow for multiple SNP sets (RECOMMENDED FOR GWAS)
#' # Create multiple gene sets
#' gene_sets_spa <- lapply(1:100, function(i) {
#'   matrix(rbinom(n * 5, 2, runif(1, 0.001, 0.01)), n, 5)  # Very rare variants
#' })
#'
#' # Step 1: Fit null model once
#' null_obj_spa <- fit_null_model(X, Y, trait = "binary")
#'
#' # Step 2: Reuse for all gene sets (much faster!)
#' results_spa_optimized <- lapply(gene_sets_spa, function(G_gene) {
#'   getZ_marg_score_binary_SPA(G_gene, null_model = null_obj_spa)
#' })
#'
#' # This is ~10-100x faster than refitting null model for each gene set
#'
#' @references
#' Dey, R., Schmidt, E. M., Abecasis, G. R., and Lee, S. (2017). A fast and accurate
#' algorithm to test for binary phenotypes and its application to PheWAS. American
#' Journal of Human Genetics, 101(1), 37-49. doi:10.1016/j.ajhg.2017.05.014
#'
#' Zhou, W., Nielsen, J. B., Fritsche, L. G., et al. (2018).
#' Efficiently controlling for case-control imbalance and sample relatedness in
#' large-scale genetic association studies. \emph{Nature Genetics}, 50(9), 1335-1341.
#' \doi{10.1038/s41588-018-0184-y}
#'
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @seealso
#' \code{\link{getZ_marg_score}} for the standard (non-SPA) version suitable for
#' common variants and balanced designs.
#'
#' \code{\link[SPAtest]{ScoreTest_SPA}} for details on the SPA implementation.
#'
#' @import SPAtest
#' @export
getZ_marg_score_binary_SPA <- function(G, X = NULL, Y = NULL, null_model = NULL) {
  # ========== Input Mode Selection ==========
  # Two modes: (1) Standard - provide X and Y, (2) Optimized - provide null_model

  if (is.null(null_model)) {
    # ========== Mode 1: Standard (fit null model) ==========
    if (is.null(X) || is.null(Y)) {
      stop("Either provide (X, Y) or null_model. ",
           "For standard mode: provide X and Y. ",
           "For optimized mode: provide null_model from fit_null_model().")
    }

    # Validate inputs using helper function (binary trait only)
    validated <- .validate_score_inputs(G, X, Y, trait = "binary")
    G <- validated$G
    X <- validated$X
    Y <- validated$Y

    # If the first column of X is not all 1's, we need to add it for intercept.
    if(ncol(X) == 0 || !all(X[,1] == 1)){
      X = cbind(1,X)
    }

  } else {
    # ========== Mode 2: Optimized (reuse null_model) ==========
    if (!inherits(null_model, "glow_null_model")) {
      stop("null_model must be created by fit_null_model()")
    }

    if (null_model$trait != "binary") {
      stop("getZ_marg_score_binary_SPA requires a binary trait null_model. ",
           "The provided null_model has trait = '", null_model$trait, "'")
    }

    if (!is.null(X) || !is.null(Y)) {
      warning("null_model provided; ignoring X and Y parameters")
    }

    # Validate G dimensions only
    if (!is.matrix(G)) {
      stop("G must be a matrix")
    }
    if (any(is.na(G))) {
      stop("G contains missing values (NA). Please impute or remove before analysis")
    }
    if (nrow(G) != null_model$n) {
      stop("G has ", nrow(G), " rows but null_model was fit with ",
           null_model$n, " samples. Dimensions must match.")
    }

    # Extract components from null_model
    X <- null_model$X
    Y0 <- null_model$fitted_probs
    Y <- null_model$Y
  }

  # ========== Main Computation (same for both modes) ==========

  # Get fitted probabilities (from null_model or fit new)
  if (is.null(null_model)) {
    # Mode 1: Fit null model
    mod0 <- glm(Y ~ X, family = "binomial")
    Y0 <- mod0$fitted.values
  }
  # Mode 2: Y0 already extracted from null_model above

  # Compute score statistics
  score = t(G)%*%(Y - Y0)
  n_var <- ncol(G)

  # Step 1: Variance structure (computed BEFORE the SPA call so degenerate
  # variants can be flagged first). Not SPA-adjusted; identical to the non-SPA
  # path. Uses the helper for the weighted projection matrices.
  w <- sqrt(Y0 * (1 - Y0))
  proj_matrices <- .compute_projection_matrices(X, G, weights = w)
  Hhalf <- proj_matrices$Hhalf
  Gtilde <- proj_matrices$Gtilde

  # Compute GHG (its diagonal is each variant's adjusted score variance).
  GHhalf <- t(Gtilde) %*% Hhalf
  GHG <- t(Gtilde) %*% Gtilde - GHhalf %*% t(GHhalf)

  # Flag degenerate variants: (near-)zero adjusted score variance, e.g. a variant
  # monomorphic in this sample or collinear with the covariates. SPAtest's
  # saddlepoint routine divides by sqrt(var1) and fails on var1 == 0 ("missing
  # value where TRUE/FALSE needed"); such a variant carries no single-variant
  # signal. diag(GHG) IS that adjusted variance ($var1$), so this is the exact
  # (not heuristic) failure condition. We exclude these from the SPA call and
  # report them as NA. The "BE" cutoff is computed per variant inside SPAtest, so
  # excluding degenerate variants does not change the p-values of the rest.
  var_adj <- diag(GHG)
  pos_var <- var_adj[is.finite(var_adj) & var_adj > 0]
  var_scale <- if (length(pos_var)) stats::median(pos_var) else 1
  ok_var <- is.finite(var_adj) & var_adj > 1e-8 * var_scale

  # Step 2: SPA p-values for the well-conditioned variants; NA for degenerate.
  pval_spa <- rep(NA_real_, n_var)
  if (any(ok_var)) {
    pval_spa[ok_var] <- SPAtest::ScoreTest_SPA(
      genos = t(G[, ok_var, drop = FALSE]),  # variants in rows, samples in cols
      pheno = Y,          # Binary phenotype (0/1)
      cov = X,  # Covariates (NULL if none, intercept added internally by SPAtest)
      method = "fastSPA", # Use fast SPA with partially normal approximation
      minmac = 0,         # Include all variants regardless of MAC
      Cutoff = "BE"       # Bonferroni-equivalent cutoff for SPA application
    )$p.value
  }
  if (!all(ok_var)) {
    warning(sprintf(
      paste0("getZ_marg_score_binary_SPA: %d of %d variant(s) had zero adjusted ",
             "variance (e.g. monomorphic in this sample) and were set to NA."),
      sum(!ok_var), n_var))
  }

  # Step 3: Convert SPA p-values to Z-scores, preserving direction of effect.
  # Z = Phi^(-1)(1 - p/2) * sign(score); degenerate variants stay NA.
  Zscores_spa <- qnorm(pval_spa / 2, lower.tail = FALSE) * sign(score)

  # Step 4: Correlation matrix. Degenerate variants have zero variance, so their
  # correlations are undefined: report them as NA rather than the NaN (with a
  # "non-positive diag" warning) that cov2cor() would produce. cov2cor() on the
  # well-conditioned submatrix is identical to the full result for those rows.
  M <- matrix(NA_real_, n_var, n_var)
  if (any(ok_var)) {
    M[ok_var, ok_var] <- cov2cor(GHG[ok_var, ok_var, drop = FALSE])
  }
  s0 <- 1  #dispersion parameter for binary trait

  # Return results
  return(list(
    Zscores = as.vector(Zscores_spa),
    scores = score,
    M_Z = M,
    M_s = GHG,
    s0 = s0
  ))
}


#' Fit Null Model for Efficient Repeated Score Calculations
#'
#' @description
#' Pre-computes the null model (without genetic effects) for efficient repeated
#' calls to \code{getZ_marg_score()} or \code{getZ_marg_score_binary_SPA()} with
#' different SNP sets but the same phenotype and covariates. This optimization is
#' critical for genome-wide analysis where the same null model is used across
#' thousands of gene/region tests.
#'
#' @details
#' In genome-wide association studies, researchers typically test many SNP sets
#' (genes, pathways, regions) against the same phenotype and covariates. The
#' standard workflow repeatedly fits the null model (expensive GLM fitting),
#' even though it's identical across all tests.
#'
#' This function enables a two-step optimized workflow:
#' \enumerate{
#'   \item Fit null model once: \code{null_obj <- fit_null_model(X, Y, trait)}
#'   \item Reuse for all SNP sets: \code{getZ_marg_score(G_i, null_model = null_obj)}
#' }
#'
#' \strong{For Binary Traits:}
#'
#' Fits logistic regression \eqn{\text{logit}(P(Y=1)) = X\beta} and stores:
#' \itemize{
#'   \item Fitted probabilities \eqn{Y_0 = P(Y=1|X)} (for weighting)
#'   \item Covariate matrix with intercept (for projection)
#'   \item Response vector (needed for SPA calculations)
#'   \item Sample size (for validation)
#' }
#'
#' \strong{For Continuous Traits:}
#'
#' Fits linear regression \eqn{Y = X\beta + \epsilon} and stores:
#' \itemize{
#'   \item Residuals \eqn{r = Y - X\hat{\beta}}
#'   \item Residual standard deviation \eqn{s_0 = \text{SD}(r)}
#'   \item Covariate matrix with intercept (for projection)
#'   \item Sample size (for validation)
#' }
#'
#' \strong{Performance Impact:}
#'
#' For testing 1000 genes with 500 samples:
#' \itemize{
#'   \item Standard mode: ~50 seconds (fit null model 1000 times)
#'   \item Optimized mode: ~5 seconds (fit null model once, reuse 1000 times)
#'   \item Speedup: 10x faster
#' }
#'
#' Speedup increases with more tests and larger samples.
#'
#' @param X A numeric matrix of covariates (n x k), where k is the number of covariates.
#'   If the first column is not all 1s (intercept), it will be added automatically.
#' @param Y A numeric vector of response variable (length n). For binary traits, must be
#'   coded as 0/1. For continuous traits, any numeric values are allowed.
#' @param trait Character string indicating trait type. Must be either "binary" (logistic
#'   regression) or "continuous" (linear regression). Default is "binary".
#' @param sample_id Optional character vector of sample IDs (length n). Stored in the
#'   null model for downstream sample alignment in \code{marginal_scan()}.
#'   Default is \code{NULL}.
#'
#' @return A null_model object (list) with class "glow_null_model" containing:
#' \describe{
#'   \item{trait}{Character, "binary" or "continuous"}
#'   \item{n}{Integer, sample size (number of observations)}
#'   \item{X}{Numeric matrix, covariate matrix with intercept included}
#'   \item{fitted_probs}{(Binary only) Numeric vector of fitted probabilities P(Y=1|X)}
#'   \item{Y}{Numeric vector, response (needed for SPA and use_lm_t)}
#'   \item{residuals}{(Continuous only) Numeric vector of residuals}
#'   \item{s0}{(Continuous only) Numeric scalar, residual standard deviation}
#'   \item{H_half}{Numeric matrix (n x k+1), precomputed projection matrix for
#'     efficient marginal scan. For binary: X_tilde (X_tilde' X_tilde)^{-1/2};
#'     for continuous: X (X' X)^{-1/2}.}
#'   \item{weights}{(Binary only) Numeric vector of length n, sqrt(mu_0 * (1 - mu_0))}
#'   \item{sample_id}{Character vector of sample IDs (if provided)}
#' }
#'
#' @examples
#' # Example 1: Binary trait - standard vs. optimized workflow
#' set.seed(123)
#' n <- 500
#' X <- matrix(rnorm(n * 2), n, 2)
#' Y <- sample(c(0, 1), n, replace = TRUE)
#'
#' # Create 10 different SNP sets
#' gene_sets <- lapply(1:10, function(i) {
#'   matrix(rbinom(n * 5, 2, 0.1), n, 5)
#' })
#'
#' # Standard workflow (slow - fits null model 10 times)
#' \dontrun{
#' results_slow <- lapply(gene_sets, function(G) {
#'   getZ_marg_score(G, X, Y, trait = "binary")
#' })
#' }
#'
#' # Optimized workflow (fast - fits null model once)
#' null_obj <- fit_null_model(X, Y, trait = "binary")
#' results_fast <- lapply(gene_sets, function(G) {
#'   getZ_marg_score(G, null_model = null_obj)
#' })
#'
#' # Example 2: Continuous trait
#' Y_cont <- rnorm(n)
#' null_obj_cont <- fit_null_model(X, Y_cont, trait = "continuous")
#' result <- getZ_marg_score(gene_sets[[1]], null_model = null_obj_cont)
#'
#' # Example 3: Using with SPA-adjusted method
#' null_obj_spa <- fit_null_model(X, Y, trait = "binary")
#' result_spa <- getZ_marg_score_binary_SPA(gene_sets[[1]], null_model = null_obj_spa)
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @seealso
#' \code{\link{getZ_marg_score}} for standard score calculation
#' \code{\link{getZ_marg_score_binary_SPA}} for SPA-adjusted score calculation
#'
#' @export
fit_null_model <- function(X, Y, trait = "binary", sample_id = NULL) {
  # Validate inputs using existing helper
  validated <- .validate_score_inputs(
    G = matrix(0, nrow = length(Y), ncol = 1),  # Dummy G for validation
    X = X,
    Y = Y,
    trait = trait
  )
  X <- validated$X
  Y <- validated$Y
  n <- validated$n

  # Add intercept if not present
  if (ncol(X) == 0 || !all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }

  # Fit null model based on trait type
  if (trait == "binary") {
    # Fit logistic regression null model
    mod0 <- glm(Y ~ X, family = "binomial")

    # Extract fitted probabilities
    fitted_probs <- mod0$fitted.values

    # Precompute projection matrix H_half (depends only on X, not G)
    w <- sqrt(fitted_probs * (1 - fitted_probs))
    X_tilde <- X * w
    XtX_chol <- chol(crossprod(X_tilde))
    H_half <- X_tilde %*% backsolve(XtX_chol, diag(ncol(XtX_chol)))

    # Create null_model object
    null_model <- list(
      trait = "binary",
      n = n,
      X = X,
      fitted_probs = fitted_probs,
      Y = Y,  # Needed for SPA calculations
      H_half = H_half,
      weights = w,  # sqrt(mu_0 * (1 - mu_0))
      sample_id = sample_id
    )

  } else if (trait == "continuous") {
    # Fit linear regression null model
    mod0 <- glm(Y ~ X, family = "gaussian")

    # Extract residuals and estimate dispersion
    residuals <- residuals(mod0)
    s0 <- sd(residuals)

    # Precompute projection matrix H_half (depends only on X, not G)
    XtX_chol <- chol(crossprod(X))
    H_half <- X %*% backsolve(XtX_chol, diag(ncol(XtX_chol)))

    # Create null_model object
    null_model <- list(
      trait = "continuous",
      n = n,
      X = X,
      residuals = residuals,
      s0 = s0,
      Y = Y,  # Needed for use_lm_t option
      H_half = H_half,
      sample_id = sample_id
    )

  } else {
    stop("trait must be either 'binary' or 'continuous'")
  }

  # Add class for method dispatch and validation
  class(null_model) <- c("glow_null_model", "list")

  return(null_model)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Validate Input Data for Score Calculation
#'
#' @description
#' Internal helper function to validate input dimensions, types, and check for
#' missing values. Used by both getZ_marg_score() and getZ_marg_score_binary_SPA().
#'
#' @param G Genotype matrix (n x p)
#' @param X Covariate matrix (n x k)
#' @param Y Response vector (length n)
#' @param trait Character string: "binary" or "continuous" (optional check)
#'
#' @return List with validated and potentially converted inputs:
#'   \item{G}{Validated genotype matrix}
#'   \item{X}{Validated covariate matrix}
#'   \item{Y}{Validated response vector}
#'   \item{n}{Number of observations}
#'   \item{p}{Number of SNPs}
#'   \item{k}{Number of covariates}
#'
#' @keywords internal
#' @noRd
.validate_score_inputs <- function(G, X, Y, trait = NULL) {
  # Check matrix types and convert if necessary
  if (!is.matrix(G)) {
    stop("G must be a matrix")
  }
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.vector(Y) && !is.matrix(Y)) {
    stop("Y must be a vector or single-column matrix")
  }
  if (is.matrix(Y)) {
    if (ncol(Y) != 1) {
      stop("Y must be a single column if provided as matrix")
    }
    Y <- as.vector(Y)
  }

  # Check dimensions
  n_G <- nrow(G)
  n_X <- nrow(X)
  n_Y <- length(Y)

  if (n_G != n_X || n_G != n_Y) {
    stop("G, X, and Y must have the same number of observations. ",
         "Currently: G has ", n_G, " rows, X has ", n_X, " rows, Y has ", n_Y, " elements")
  }

  # Check for missing values
  if (any(is.na(G))) {
    stop("G contains missing values (NA). Please impute or remove before analysis")
  }
  if (any(is.na(X))) {
    stop("X contains missing values (NA). Please remove before analysis")
  }
  if (any(is.na(Y))) {
    stop("Y contains missing values (NA). Please remove before analysis")
  }

  # Validate trait-specific requirements
  if (!is.null(trait)) {
    if (!trait %in% c("binary", "continuous")) {
      stop("trait must be either 'binary' or 'continuous'")
    }

    # For binary trait, check Y is 0/1
    if (trait == "binary") {
      if (!all(Y %in% c(0, 1))) {
        stop("For binary trait, Y must contain only 0 and 1")
      }
    }
  }

  # Return validated inputs with dimensions
  list(
    G = G,
    X = X,
    Y = Y,
    n = n_G,
    p = ncol(G),
    k = ncol(X)
  )
}


#' Compute Projection Matrices for Score Calculation
#'
#' @description
#' Internal helper function to compute weighted design matrices and projection
#' matrices used in score variance calculations. Handles both weighted (binary trait)
#' and unweighted (continuous trait) cases.
#'
#' @param X Covariate matrix (n x k)
#' @param G Genotype matrix (n x p)
#' @param weights Optional weight vector (length n). If NULL, unweighted computation.
#'
#' @return List with computed matrices:
#'   \item{Hhalf}{Projection matrix H^(1/2)}
#'   \item{Gtilde}{Weighted genotype matrix (or G if unweighted)}
#'   \item{Xtilde}{Weighted covariate matrix (or X if unweighted)}
#'
#' @details
#' Computes H^(1/2) = Xtilde * (Xtilde' * Xtilde)^(-1/2) where:
#' - For weighted case: Xtilde = sqrt(weights) * X
#' - For unweighted case: Xtilde = X
#'
#' Uses Cholesky decomposition for numerical stability:
#' (Xtilde' * Xtilde)^(-1/2) = t(chol(solve(Xtilde' * Xtilde)))
#'
#' @keywords internal
#' @noRd
.compute_projection_matrices <- function(X, G, weights = NULL) {
  # Apply weights if provided
  if (!is.null(weights)) {
    # Weighted case (binary trait)
    Xtilde <- X * weights
    Gtilde <- G * weights
  } else {
    # Unweighted case (continuous trait)
    Xtilde <- X
    Gtilde <- G
  }

  # Compute H^(1/2) = Xtilde * (Xtilde' * Xtilde)^(-1/2)
  # Check for numerical stability before Cholesky
  XtX <- t(Xtilde) %*% Xtilde
  XtX_inv <- tryCatch(
    solve(XtX),
    error = function(e) {
      if (!is.null(weights)) {
        stop("Covariate matrix X is singular or near-singular after weighting. ",
             "This may occur with separation in logistic regression or collinear covariates. ",
             "Check for redundant covariates or perfectly predicted outcomes.")
      } else {
        stop("Covariate matrix X is singular or near-singular. ",
             "Check for collinear covariates or constant columns.")
      }
    }
  )
  Hhalf <- Xtilde %*% t(chol(XtX_inv))

  list(
    Hhalf = Hhalf,
    Gtilde = Gtilde,
    Xtilde = Xtilde
  )
}


#' Handle Zero-Variance SNPs in Score Calculation
#'
#' @description
#' Internal helper function to identify zero-variance SNPs and prevent division
#' by zero in Z-score calculation.
#'
#' @param diag_GHG Diagonal elements of GHG matrix (variance terms)
#' @param Zscore Vector of Z-scores (will be modified in place)
#'
#' @return List with:
#'   \item{zero_var_idx}{Indices of zero-variance SNPs}
#'   \item{diag_GHG_safe}{Modified diagonal with 1 replacing zero-variance elements}
#'   \item{has_zero_var}{Logical indicating if any zero-variance SNPs found}
#'
#' @keywords internal
#' @noRd
.handle_zero_variance_snps <- function(diag_GHG) {
  # Identify zero-variance SNPs
  zero_var_idx <- which(diag_GHG < .Machine$double.eps)
  has_zero_var <- length(zero_var_idx) > 0

  # Create safe version for division
  diag_GHG_safe <- diag_GHG
  if (has_zero_var) {
    # Replace with 1 to prevent division by zero
    diag_GHG_safe[zero_var_idx] <- 1
  }

  list(
    zero_var_idx = zero_var_idx,
    diag_GHG_safe = diag_GHG_safe,
    has_zero_var = has_zero_var
  )
}
