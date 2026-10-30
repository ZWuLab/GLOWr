########## Optimal Weights Calculation Functions ##########
#
# This file contains the core optimal weighting methodology for the GLOW
# framework. The main user-facing function is Optimal_Weights_M(), which
# computes optimal weights for variant-set tests.
#
# MAIN FUNCTION: Optimal_Weights_M()
#   This is a core user-facing function for calculating optimal weights.
#   It supports different test types (Burden, SKAT, Fisher) and approximation
#   methods (normal, sparse mixture).



#################### Main Optimal Weights Function ####################

#' Calculate Optimal Weights for GLOW Tests
#'
#' @description
#' Main function for computing optimal weights for variant-set tests in the
#' GLOW framework. Computes weights under different hypotheses (H0 vs H1) and
#' approximation methods (normal vs sparse), tailored to the specific test
#' transformation function g.
#'
#' @param g Function, transformation applied to Z-scores (e.g., identity for
#'   Burden, GFisher transformation for Fisher-type tests)
#' @param Bstar Numeric vector, estimated effect sizes for each variant
#' @param PI Numeric vector, estimated variant-importance scores
#' @param M Numeric matrix, correlation matrix of marginal Z-scores
#' @param is.posi.wts Logical, if TRUE, negative weights are set to zero (default TRUE)
#'
#' @return A list with optimal weights. The structure depends on the transformation function:
#'
#' For Burden/Liptak tests (when g is identity):
#' \describe{
#'   \item{wts_BE}{Weights optimized for Best Effect under H0}
#'   \item{wts_APE}{Weights optimized for Average Power under Ensemble under H1}
#' }
#'
#' For other tests (Fisher, SKAT, etc.):
#' \describe{
#'   \item{wts_BE_N}{BE weights using normal approximation}
#'   \item{wts_APE_N}{APE weights using normal approximation}
#'   \item{wts_BE_sparse}{BE weights using sparse mixture approximation}
#'   \item{wts_APE_sparse}{APE weights using sparse mixture approximation}
#' }
#'
#' @details
#' This function orchestrates the optimal weight calculation for GLOW tests.
#' The procedure differs based on the transformation function g:
#'
#' \strong{For Burden/Liptak tests (g is identity):}
#'
#' The weights have closed-form solutions:
#' \deqn{wts_{BE} = MU}
#' \deqn{wts_{APE} = (I + diag(V)M)^{-1} MU}
#' where \eqn{MU = Bstar * PI} and \eqn{V = Bstar^2 * PI * (1-PI)}.
#'
#' \strong{For other tests (GFisher, SKAT, etc.):}
#'
#' Two approximation methods are used:
#'
#' 1. \emph{Normal Approximation}: Assumes Z-scores are approximately normal
#'    under the alternative, \eqn{\tilde{Z} \sim N(MMU, M + MVM)}, where:
#'    \itemize{
#'      \item MMU = M \%*\% MU (marginal means)
#'      \item MVM = t(M * V) \%*\% M (marginal variance contribution from random effects)
#'    }
#'    The moments of \eqn{g(\tilde{Z}_i)} are computed on the standardized scale. With
#'    \eqn{s_i = \sqrt{(M + MVM)_{ii}}} and the correlation matrix
#'    \eqn{R = \mathrm{cov2cor}(M + MVM)}, \eqn{\tilde{Z}_i = MMU_i + s_i U_i} with
#'    \eqn{U \sim N(0, R)}: \code{get_r_tilde} receives the standard deviations
#'    \eqn{s_i}, and \code{CovM_gXgY} receives \eqn{R} together with the scales
#'    \eqn{s_i}, because Mehler's expansion is an expansion in a correlation (bounded by
#'    one, which the truncation at order 64 relies on) with the Hermite coefficients of
#'    \eqn{g(MMU_i + s_i U)}. GLOWr up to 0.3.2, like the legacy code, passed the
#'    covariance \eqn{M + MVM} into the expansion and its diagonal, a variance, as the
#'    standard deviation; the \code{_N} weights of those versions differ from the current ones.
#'
#' 2. \emph{Sparse Mixture Approximation}: Explicitly models the sparse signal
#'    structure Z_i = Z_0i + Bstar_i * C_i where C_i ~ Bernoulli(PI_i).
#'
#' \strong{Numerical integration and the expansion order:}
#'
#' Every expectation of \eqn{g} is a one-dimensional integral against the standard normal
#' density. It is computed by a composite 10-point Gauss-Legendre rule on \eqn{[-12, 12]}
#' with panels of width 0.1, cut at the point where the argument of \eqn{g} is zero (the
#' kink of the two-sided transforms, which depend on \eqn{|x|}). The Hermite coefficients
#' \eqn{b_k = E[g(\mu + sU)\,he_k(U)]} come from one evaluation of \eqn{g} per distinct
#' \eqn{(\mu, s)} pair, and Mehler's expansion
#' \eqn{\mathrm{Cov}\{g(\mu_i + s_iU_i), g(\mu_j + s_jU_j)\} = \sum_k b_k(i)b_k(j)r_{ij}^k}
#' runs to order 64. The truncated sum is the covariance of the projected transform, and its
#' diagonal is the truncated variance; the sparse covariance carries the exact mixture
#' variance on its diagonal instead, because the mixture formula does not apply to a variant
#' with itself. Against a direct integration of every pairwise covariance on a 20-variant
#' window with a near-duplicate pair, the weights of the normal approximation agree to
#' \eqn{10^{-3}} of the largest weight and the sparse weights to \eqn{5 \times 10^{-3}}
#' (the package tests). GLOWr up to 0.3.3 used \code{integrate()} at its default tolerance
#' and order 8; on that window its \code{_N} weights were off by \eqn{9 \times 10^{-2}}.
#'
#' The rule is validated for the Fisher and polynomial transforms over the scales GLOW meets.
#' A custom \code{g} must be vectorized, square-integrable against the normal density and
#' smooth except possibly at zero; a transform with another kink, a cusp, a discontinuity or
#' rapid growth needs its own convergence check.
#'
#' For each approximation, weights are computed under two hypotheses:
#' \itemize{
#'   \item BE (Best under H0): Assumes no signal, optimizes for power under worst-case
#'   \item APE (Average Power under Ensemble): Accounts for signal, optimizes for average power
#' }
#'
#' \strong{Transformation Function Check:}
#'
#' The function checks if g is the identity by comparing function bodies:
#' \code{identical(body(g), body(function(x) x))}
#'
#' @section Computational Complexity:
#' \itemize{
#'   \item Burden/Liptak: O(n^3) for matrix inversion
#'   \item Other tests: O(n_pairs * n_nodes * K + 4 * n^2 * K + n^3) where:
#'     \itemize{
#'       \item n = number of variants
#'       \item n_pairs = distinct (mean, scale) pairs, at most 2n + 1
#'       \item n_nodes = quadrature nodes per pair (2,410)
#'       \item K = number of Hermite polynomial orders (64)
#'       \item Factor of 4 from the four Mehler sums (the null, the normal-approximation
#'         alternative, and two sparse components)
#'     }
#' }
#'
#' @section Assumptions:
#' \itemize{
#'   \item Bstar, PI, and M have compatible dimensions
#'   \item M is a valid correlation matrix (positive definite, diagonal = 1)
#'   \item 0 <= PI_i <= 1 for all i
#'   \item Bstar values are positive (standardized effect sizes)
#'   \item For non-Burden tests, g(mu + s U) is square-integrable under the normal density
#'     and smooth except possibly at zero (see the details)
#' }
#'
#' @section Weight Interpretation:
#' \itemize{
#'   \item Larger weights indicate variants that contribute more to the test statistic
#'   \item All weights are normalized by mean(abs(weights)) for comparability
#'   \item BE weights are conservative (good when signal is uncertain)
#'   \item APE weights are adaptive (good when signal estimates are reliable)
#'   \item Sparse weights are more accurate for rare causal variants
#'   \item Normal weights are computationally simpler approximations
#' }
#'
#' @section Performance Optimization:
#' The function uses several optimizations for improved performance:
#' \itemize{
#'   \item H0 covariance matrices are cached to speed up repeated calls
#'   \item Fast path computation for independent variants (identity correlation)
#'   \item Caching can be controlled: \code{options(GLOWr.use_cache = TRUE/FALSE)}
#'   \item Cache can be cleared: \code{clear_glow_cache()}
#'   \item Cache status: \code{glow_cache_info()}
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Burden test (identity transformation)
#' g_identity <- function(x) x
#' Bstar <- c(0.5, 1.0, 0.3)
#' PI <- c(0.01, 0.05, 0.02)
#' M <- diag(3)
#' weights_burden <- Optimal_Weights_M(g_identity, Bstar, PI, M)
#' names(weights_burden)  # wts_BE, wts_APE
#'
#' # Example 2: GFisher test
#' g_GFisher_two <- function(x, df = 2) {
#'   qchisq(log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
#'          df = df, lower.tail = FALSE, log.p = TRUE)
#' }
#' weights_fisher <- Optimal_Weights_M(
#'   function(x) g_GFisher_two(x, df = 2),
#'   Bstar, PI, M
#' )
#' names(weights_fisher)  # wts_BE_N, wts_APE_N, wts_BE_sparse, wts_APE_sparse
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. (2026). Integrated weighted
#' association test with application to genetic association studies. The Annals
#' of Applied Statistics, 20(3), 2074-2097. \doi{10.1214/26-AOAS2223}
#'
#' @export
Optimal_Weights_M <- function(g, Bstar, PI, M, is.posi.wts = TRUE) {
  # Input validation for NA/Inf
  # Note: Bstar can be negative (protective alleles), so we allow negative values
  validate_numeric_input(Bstar, "Bstar (effect sizes)", allow_negative = TRUE)
  validate_numeric_input(PI, "PI (variant-importance scores)", allow_negative = FALSE)
  validate_numeric_input(M, "M (correlation matrix)")

  # Additional validation
  if (length(Bstar) != length(PI)) {
    stop("Bstar and PI must have the same length")
  }
  if (length(Bstar) != nrow(M) || length(Bstar) != ncol(M)) {
    stop("Bstar, PI, and M must have compatible dimensions")
  }
  if (any(PI < 0 | PI > 1)) {
    stop("PI values must be in [0, 1]")
  }

  # Check that M is positive semi-definite
  check_positive_definite(M, "Correlation matrix M")

  # Calculate mean and variance of random effects (theta)
  MU <- as.matrix(Bstar * PI)  # Mean of theta: E[theta_i] = Bstar_i * PI_i
  V <- c(Bstar^2 * PI * (1 - PI))  # Variance of theta: Var[theta_i] = Bstar_i^2 * PI_i * (1 - PI_i)

  # Calculate marginal means and covariances when effects are random
  MMU <- M %*% MU  # Marginal means: E[Z_i] under random effects
  # Marginal variance: Var[Z] = M + M*diag(V)*M (fast computation for diagonal V)
  MVM <- t(M * V) %*% M  # M * diag(V) * M, exploiting diagonal structure of V

  # Check if this is a Burden/Liptak test (g is identity function)
  if (identical(body(g), body(function(x) x))) {
    # --- Special case: Burden/Liptak test ---
    # For Burden, we have closed-form solutions

    # BE weights: Simply the mean vector MU
    wts_BE <- MU

    # APE weights: Account for correlation and variance
    # wts_APE = (I + diag(V) * M)^{-1} * MU
    wts_APE <- safe_matrix_inverse(
      diag(length(MU)) + diag(V) %*% M,
      "APE weight calculation for Burden test"
    ) %*% MU

    # Apply positive weights constraint if requested
    if (is.posi.wts) {
      wts_BE[wts_BE < 0] <- 0
      wts_APE[wts_APE < 0] <- 0
    }

    # Normalize weights by mean absolute value with division by zero protection
    mean_abs_BE <- mean(abs(wts_BE))
    mean_abs_APE <- mean(abs(wts_APE))

    if (mean_abs_BE < .Machine$double.eps) {
      warning("All BE weights are effectively zero. Returning equal weights.",
              call. = FALSE)
      wts_BE <- rep(1/length(wts_BE), length(wts_BE))
    } else {
      wts_BE <- wts_BE / mean_abs_BE
    }

    if (mean_abs_APE < .Machine$double.eps) {
      warning("All APE weights are effectively zero. Returning equal weights.",
              call. = FALSE)
      wts_APE <- rep(1/length(wts_APE), length(wts_APE))
    } else {
      wts_APE <- wts_APE / mean_abs_APE
    }

    return(list(wts_BE = wts_BE, wts_APE = wts_APE))

  } else {
    # --- General case: Non-Burden tests (GFisher, SKAT, etc.) ---

    ## 1 and 2. The normal approximation and the sparse mixture, from one coefficient
    ## table over the null pair (0, 1), the normal-approximation pairs (MMU_i, s_i) and the
    ## sparse pairs (Bstar_i, 1): see .optimal_weight_moments()
    mom <- .optimal_weight_moments(g, Bstar, PI, M)

    # Weights under the normal approximation (BE: the null covariance; APE: the alternative's)
    wts_BE_N <- get_wts(mom$Sigma_BE_N, mom$r_N, is.posi.wts)
    wts_APE_N <- get_wts(mom$Sigma_APE_N, mom$r_N, is.posi.wts)

    # Weights under the sparse mixture approximation
    wts_BE_sparse <- get_wts(mom$Sigma_BE_sparse, mom$r_sparse, is.posi.wts)
    wts_APE_sparse <- get_wts(mom$Sigma_APE_sparse, mom$r_sparse, is.posi.wts)

    # Return all four sets of weights
    return(list(
      wts_BE_N = wts_BE_N$w_normalized,
      wts_APE_N = wts_APE_N$w_normalized,
      wts_BE_sparse = wts_BE_sparse$w_normalized,
      wts_APE_sparse = wts_APE_sparse$w_normalized
    ))
  }
}



#################### Weight Calculation Functions ####################

#' Calculate Optimal Weights from Covariance and Mean Vector
#'
#' @description
#' Computes optimal weights for combining test statistics based on their
#' covariance matrix and mean vector. Optionally constrains weights to be
#' non-negative.
#'
#' @param Sigma Numeric matrix, covariance matrix of test statistics
#' @param r Numeric vector, mean difference vector (E_1 - E_0)
#' @param is.posi.wts Logical, if TRUE, negative weights are set to zero
#'
#' @return A list with elements:
#' \describe{
#'   \item{w}{Raw optimal weights}
#'   \item{w_normalized}{Weights normalized by mean absolute value}
#' }
#'
#' @details
#' The optimal weights maximize power for detecting a signal with mean r
#' and covariance Sigma. The weights are computed as:
#' \deqn{w = \Sigma^{-1} r}
#'
#' This is the solution to maximizing the non-centrality parameter:
#' \deqn{\lambda = w^T r / \sqrt{w^T \Sigma w}}
#'
#' When is.posi.wts = TRUE, negative weights are truncated to zero:
#' \deqn{w_i = max(0, w_i)}
#'
#' Weights are then normalized by dividing by the mean absolute value:
#' \deqn{w_{normalized} = w / mean(|w|)}
#'
#' This normalization ensures that the weights are on a comparable scale
#' across different scenarios.
#'
#' @section Computational Complexity:
#' O(n^3) for matrix inversion where n = length(r)
#'
#' @section Assumptions:
#' \itemize{
#'   \item Sigma is positive definite (invertible)
#'   \item r and Sigma have compatible dimensions
#' }
#'
#' @examples
#' \dontrun{
#' # Simple example with independent statistics
#' Sigma <- diag(3)
#' r <- c(1, 2, 0.5)
#' wts <- get_wts(Sigma, r, is.posi.wts = FALSE)
#' wts$w_normalized
#'
#' # With positive weights constraint
#' r_mixed <- c(1, -0.5, 2)
#' wts_pos <- get_wts(Sigma, r_mixed, is.posi.wts = TRUE)
#' wts_pos$w_normalized  # Negative weight set to 0
#' }
#'
#' @keywords internal
#' @noRd
get_wts <- function(Sigma, r, is.posi.wts = TRUE) {
  # Compute optimal weights: w = Sigma^{-1} * r
  w <- safe_matrix_inverse(Sigma, "optimal weight calculation") %*% r

  # Apply positive weights constraint if requested
  if (is.posi.wts) {
    w <- w * (w > 0)  # Set negative weights to 0
  }

  # Normalize weights by mean absolute value
  # This makes weights comparable across different scenarios
  # Add protection against division by zero
  mean_abs_w <- mean(abs(w))
  if (mean_abs_w < .Machine$double.eps) {
    warning("All weights are effectively zero in optimal weight calculation. ",
            "Returning equal weights. This may indicate issues with parameter estimates.",
            call. = FALSE)
    w_normalized <- rep(1/length(w), length(w))
  } else {
    w_normalized <- w / mean_abs_w
  }

  return(list(w = w, w_normalized = w_normalized))
}


#' Calculate r-tilde for Normal Approximation
#'
#' @description
#' Computes the vector of mean differences for the normal approximation
#' approach to optimal weights. This is used when approximating the distribution
#' of transformed Z-scores under normal theory.
#'
#' @param g Function, transformation to apply to Z-scores
#' @param MU Numeric vector, means of the normal distributions
#' @param SD Numeric vector, standard deviations of the normal distributions
#'
#' @return Numeric vector of mean differences (same length as MU)
#'
#' @details
#' This function calculates:
#' \deqn{r_{tilde,i} = E[g(Z_i)] - E[g(Z_0)]}
#' where \eqn{Z_i \sim N(MU_i, SD_i^2)} and \eqn{Z_0 \sim N(0, 1)}.
#'
#' The normal approximation assumes that under the alternative hypothesis,
#' each Z-score follows a normal distribution with mean MU_i and standard
#' deviation SD_i, rather than a mixture distribution.
#'
#' This is computationally simpler than the sparse approximation but may
#' be less accurate when the true distribution is a sparse mixture.
#'
#' @section Computational Complexity:
#' O(n * cost_per_integration) where n = length(MU)
#'
#' @examples
#' \dontrun{
#' # Identity transformation
#' g_identity <- function(x) x
#' MU <- c(0.5, 1.0, 0.2)
#' SD <- c(1.1, 1.2, 1.05)
#' r_tilde <- get_r_tilde(g_identity, MU, SD)
#' }
#'
#' @keywords internal
#' @noRd
get_r_tilde <- function(g, MU, SD) {
  # E[g(Z_0)] once, then E[g(Z_i)] - E[g(Z_0)] with Z_i ~ N(MU[i], SD[i]^2)
  e0 <- E_gX_p(g, mu = 0, p = 1, sigma = 1)
  vapply(seq_along(MU), function(i) E_gX_p(g, mu = MU[i], p = 1, sigma = SD[i]),
         numeric(1)) - e0
}


#' Calculate r for Sparse Approximation
#'
#' @description
#' Computes the vector of mean differences for the sparse mixture approximation
#' approach to optimal weights. This accounts for the sparse signal structure
#' where only a fraction of variants have non-zero effects.
#'
#' @param g Function, transformation to apply to Z-scores
#' @param MU Numeric vector, effect sizes (B-star values)
#' @param PI Numeric vector, variant-importance scores
#'
#' @return Numeric vector of mean differences (same length as MU)
#'
#' @details
#' This function calculates:
#' \deqn{r_i = E_1[T_i] - E_0[T_i]}
#' where \eqn{T_i = g(Z_i)} and:
#' \itemize{
#'   \item Under H1 (alternative): \eqn{Z_i = Z_{0i} + MU_i * C_i}, \eqn{C_i \sim Bern(PI_i)}
#'   \item Under H0 (null): \eqn{Z_i = Z_{0i}}, \eqn{Z_{0i} \sim N(0, 1)}
#' }
#'
#' The sparse approximation explicitly models that:
#' \itemize{
#'   \item Only a fraction PI_i of variants are causal
#'   \item Causal variants have effect size MU_i
#'   \item Non-causal variants have effect size 0
#' }
#'
#' This is more accurate than the normal approximation when the true
#' genetic architecture is sparse.
#'
#' @section Computational Complexity:
#' O(n * cost_per_integration) where n = length(MU)
#'
#' @examples
#' \dontrun{
#' # GFisher transformation
#' g_GFisher <- function(x, df = 2) {
#'   qchisq(log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
#'          df = df, lower.tail = FALSE, log.p = TRUE)
#' }
#' MU <- c(0.5, 1.0, 0.2)
#' PI <- c(0.01, 0.05, 0.001)
#' r <- get_r(function(x) g_GFisher(x, df = 2), MU, PI)
#' }
#'
#' @keywords internal
#' @noRd
get_r <- function(g, MU, PI) {
  n <- length(MU)
  # For each i, compute E[g(Z_i)] under mixture model - E[g(Z_0)] under null
  r <- sapply(1:n, function(i) {
    E_T_mix(g, MU[i], PI[i]) - E_T_mix(g, 0, 0)
  })
  return(r)
}


#################### Covariance Matrix Functions ####################

#' Covariance Matrix via Hermite Expansion
#'
#' @description
#' Computes the covariance matrix of g(X) and g(Y), where X and Y are
#' correlated multivariate normal random vectors, using Hermite polynomial
#' expansion (Mehler's formula). This is the core function for computing
#' covariances under the normal approximation.
#'
#' @param g Function, transformation to apply to normal random variables
#' @param MU1 Numeric vector, means of X components
#' @param MU2 Numeric vector, means of Y components
#' @param M Numeric matrix, correlation matrix of the shared standard component Z
#'   (unit diagonal, must be positive definite). A covariance matrix is not accepted
#'   here: pass its correlation matrix and its standard deviations through SD1 and SD2.
#' @param ORD Integer vector, the orders of the Hermite expansion summed (default 1:64;
#'   \code{1:8} reproduces the order-8 expansion of GLOWr 0.3.3)
#' @param SD1 Numeric vector, standard deviations of X components (default all 1)
#' @param SD2 Numeric vector, standard deviations of Y components (default all 1)
#'
#' @return Numeric matrix, covariance matrix where element (i,j) is Cov[g(X_i), g(Y_j)]
#'
#' @details
#' This function computes:
#' \deqn{Cov[g(X_i), g(Y_j)]}
#' where \eqn{X_i = MU1_i + SD1_i Z_i} and \eqn{Y_j = MU2_j + SD2_j Z_j}, with
#' Z ~ MVN(0, M) the shared standard component, so that
#' \eqn{Cov(X_i, Y_j) = SD1_i SD2_j M_{ij}}. With the default scales, X ~ MVN(MU1, M)
#' and Y ~ MVN(MU2, M).
#'
#' The computation uses the Hermite expansion of the transformed standard variable:
#' \deqn{g(\mu + s Z) \approx \sum_{k=0}^{K} \frac{a_k(\mu, s)}{k!} He_k(Z), \quad
#'   a_k(\mu, s) = E[g(\mu + s Z) He_k(Z)], \quad Z \sim N(0, 1),}
#' where He_k are probabilist's Hermite polynomials.
#'
#' The key property used (Mehler's formula) is:
#' \deqn{Cov[He_k(Z_i), He_l(Z_j)] = k! \, M_{ij}^k \, 1(k = l),}
#' which needs the correlation \eqn{M_{ij}}, bounded by one in absolute value: the terms
#' \eqn{M_{ij}^k / k!} then decrease and the truncation at order 64 discards the smallest
#' ones. The expansion is not valid in a covariance with entries above one. The truncated
#' sum is the covariance of the projection of \eqn{g} onto the polynomials of degree at
#' most \eqn{K}, so its diagonal is the truncated variance \eqn{\sum_k a_k^2 / k!}.
#'
#' The algorithm:
#' \enumerate{
#'   \item Compute the normalized Hermite coefficients \eqn{b_k = a_k / \sqrt{k!}} of
#'     \eqn{g(\mu + sU)} for each distinct (mean, scale) pair among the X and Y
#'     components, by the composite Gauss-Legendre rule of
#'     \code{optimal_weights_quadrature.R} (one evaluation of g per pair)
#'   \item Sum the contributions of each order k in ORD: b1[k,i] * b2[k,j] * M_ij^k
#' }
#'
#' @section Computational Complexity:
#' O(n^2 * K * n_eval) for correlated variants, O(n * K * n_eval) for independent variants:
#' \itemize{
#'   \item n = length(MU1) = length(MU2)
#'   \item K = max(ORD) (64 by default)
#'   \item n_eval = cost per numerical integration (2,410 evaluations of g per distinct pair)
#'   \item Fast path optimization when M is identity matrix
#' }
#'
#' @section Numerical Considerations:
#' \itemize{
#'   \item The Hermite functions come from their three-term recurrence, stable to order 64
#'     and beyond on the quadrature domain
#'   \item Integration over [-12, 12] standard deviations, cut where the argument of g is
#'     zero (the kink of the two-sided transforms)
#'   \item Assumes g(mu + s U) square-integrable and smooth away from zero
#'   \item Fast path for independent variants (identity correlation matrix)
#' }
#'
#' @section Mathematical Background:
#' The Hermite expansion leverages the orthogonality of Hermite polynomials
#' with respect to the normal distribution. For correlated normals with
#' correlation rho, the covariance of k-th order Hermite polynomials is
#' exactly k! * rho^k. This allows semi-analytical computation of covariances
#' of arbitrary smooth transformations.
#'
#' @examples
#' \dontrun{
#' # Identity transformation
#' g_identity <- function(x) x
#' M <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' MU1 <- c(0, 0)
#' MU2 <- c(0, 0)
#' CovM_gXgY(g_identity, MU1, MU2, M)  # Should equal M
#' }
#'
#' @references
#' Lin, D. Y., and Tang, Z. Z. (2011). A general framework for detecting
#' disease associations with rare variants in sequencing studies.
#' The American Journal of Human Genetics, 89(3), 354-367.
#'
#' @keywords internal
#' @noRd
CovM_gXgY <- function(g, MU1, MU2, M, ORD = 1:64, SD1 = NULL, SD2 = NULL) {
  # Check that M is positive semi-definite
  check_positive_definite(M, "Correlation matrix M")

  n <- length(MU1)
  if (is.null(SD1)) SD1 <- rep(1, n)
  if (is.null(SD2)) SD2 <- rep(1, n)
  if (length(MU2) != n || length(SD1) != n || length(SD2) != n) {
    stop("MU1, MU2, SD1 and SD2 must have the same length")
  }
  if (any(SD1 <= 0) || any(SD2 <= 0)) {
    stop("SD1 and SD2 must be positive")
  }
  if (any(abs(diag(M) - 1) > 1e-8)) {
    stop("M must be a correlation matrix (unit diagonal); pass cov2cor(Sigma) with SD1 and SD2")
  }
  ORD <- .check_ord(ORD)

  # Normalized Hermite coefficients b_k = E[g(mu + s U) he_k(U)] of every distinct
  # (mean, scale) pair among the X and Y components, one evaluation of g per pair
  tab <- .hermite_table(g, c(MU1, MU2), c(SD1, SD2), max(ORD))
  B1 <- tab$b[, seq_len(n), drop = FALSE]          # row k + 1 holds order k, column i: X_i
  B2 <- tab$b[, n + seq_len(n), drop = FALSE]      # column j: Y_j

  # Fast path for independent variants (identity correlation matrix): only the
  # diagonal survives, Cov[g(X_i), g(Y_i)] = sum_k b_k(X_i) b_k(Y_i)
  is_identity <- (nrow(M) == ncol(M)) &&
                 all(abs(diag(M) - 1) < 1e-10) &&
                 all(abs(M[row(M) != col(M)]) < 1e-10)
  if (is_identity) {
    diag_vals <- colSums(B1[ORD + 1L, , drop = FALSE] * B2[ORD + 1L, , drop = FALSE])
    return(diag(diag_vals, nrow = n, ncol = n))
  }

  # Mehler's expansion: Cov[g(X_i), g(Y_j)] = sum_{k in ORD} b_k(X_i) b_k(Y_j) M_ij^k
  .mehler_sum(B1, B2, M, ORD)
}


#' Covariance Matrix for Mixture Distributions
#'
#' @description
#' Computes the covariance matrix of T_i = g(Z_i) and T_j = g(Z_j), where
#' Z follows a multivariate mixture distribution with sparse signals.
#' This is used for the sparse approximation in optimal weight calculation.
#'
#' @param g Function, transformation to apply to Z-scores
#' @param MU Numeric vector, effect sizes (B-star values)
#' @param PI Numeric vector, variant-importance scores
#' @param M Numeric matrix, correlation matrix of null Z-scores
#' @param ORD Integer vector, the orders of the Hermite expansion summed (default 1:64)
#'
#' @return Numeric matrix, covariance matrix where element (i,j) is Cov[g(Z_i), g(Z_j)]
#'
#' @details
#' This function computes the covariance matrix under the sparse mixture model:
#' \deqn{Z_i = Z_{0i} + MU_i * C_i}
#' where:
#' \itemize{
#'   \item Z_0 ~ MVN(0, M) is the null Z-score vector
#'   \item C_i ~ Bernoulli(PI_i) independently across i
#'   \item MU_i is the effect size when C_i = 1
#' }
#'
#' The covariance is:
#' \deqn{Cov[T_i, T_j] = E[T_i T_j] - E[T_i]E[T_j]}
#'
#' For the mixture model with independent causal indicators:
#' \deqn{E[T_i T_j] = \sum_{c_i,c_j \in \{0,1\}} P(C_i=c_i, C_j=c_j) E[g(Z_i) g(Z_j) | C_i=c_i, C_j=c_j]}
#'
#' This expands to four terms:
#' \enumerate{
#'   \item \eqn{P(C_i=1, C_j=1) = PI_i * PI_j}: Both causal
#'   \item \eqn{P(C_i=1, C_j=0) = PI_i * (1-PI_j)}: Only i causal
#'   \item \eqn{P(C_i=0, C_j=1) = (1-PI_i) * PI_j}: Only j causal
#'   \item \eqn{P(C_i=0, C_j=0) = (1-PI_i) * (1-PI_j)}: Both null
#' }
#'
#' Each conditional covariance is computed using CovM_gXgY with appropriate means.
#'
#' For diagonal elements, the variance is computed exactly using Var_T_mix to
#' avoid numerical issues.
#'
#' @section Computational Complexity:
#' O(4 * n^2 * K * n_eval + n * n_eval) where:
#' \itemize{
#'   \item Factor of 4 from the four mixture components
#'   \item n = length(MU)
#'   \item K = length(ORD)
#'   \item Additional O(n * n_eval) for diagonal variance calculations
#'   \item H0 covariance matrices (Cov_M3) are cached to improve performance
#' }
#'
#' @section Assumptions:
#' \itemize{
#'   \item Causal indicators C_i are independent across variants
#'   \item M is a valid correlation matrix (positive definite)
#'   \item 0 <= PI_i <= 1 for all i
#' }
#'
#' @examples
#' \dontrun{
#' # GFisher transformation
#' g_GFisher <- function(x, df = 2) {
#'   qchisq(log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
#'          df = df, lower.tail = FALSE, log.p = TRUE)
#' }
#' M <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' MU <- c(1.0, 0.5)
#' PI <- c(0.01, 0.05)
#' Sigma <- CovMT_mix(function(x) g_GFisher(x, 2), MU, PI, M)
#' }
#'
#' @keywords internal
#' @noRd
CovMT_mix <- function(g, MU, PI, M, ORD = 1:64) {
  # Check that M is positive semi-definite
  check_positive_definite(M, "Correlation matrix M")

  n <- length(MU)
  ORD <- .check_ord(ORD)
  K <- max(ORD)

  # One coefficient table: the causal pairs (MU_i, 1) and the null pair (0, 1)
  tab <- .hermite_table(g, c(MU, 0), rep(1, n + 1L), K)
  B1 <- tab$b[, seq_len(n), drop = FALSE]      # coefficients of g(MU_i + U)
  b0 <- tab$b[, n + 1L]                        # coefficients of g(U)
  B0 <- matrix(b0, K + 1L, n)

  # Covariances of the four causal combinations (i causal / j causal, i causal / j null,
  # and the null pair; the null / causal case is the transpose)
  Cov_M  <- .mehler_sum(B1, B1, M, ORD)
  Cov_M1 <- .mehler_sum(B1, B0, M, ORD)

  # The null covariance Cov_M3 is identical for every call with the same null
  # coefficients, M and orders, so it is cached under a key of those numerical inputs
  # (a key from the body of g would not see the values a closure captures, and a key
  # from a probe of g would not tell two transforms apart that agree at the probes)
  use_cache <- getOption("GLOWr.use_cache", default = TRUE)
  Cov_M3 <- NULL
  if (use_cache) {
    cache_key <- tryCatch({
      digest::digest(list(b0 = b0, M = M, ORD = ORD, type = "H0_covariance"),
                     algo = "xxhash64")
    }, error = function(e) {
      warning("Cache key generation failed: ", e$message,
              ". Falling back to no caching.", call. = FALSE)
      NULL
    })
    if (!is.null(cache_key)) {
      if (exists(cache_key, envir = .glow_cache, inherits = FALSE)) {
        Cov_M3 <- get(cache_key, envir = .glow_cache)
      } else {
        Cov_M3 <- .mehler_sum(B0, B0, M, ORD)
        assign(cache_key, Cov_M3, envir = .glow_cache)
      }
    }
  }
  if (is.null(Cov_M3)) Cov_M3 <- .mehler_sum(B0, B0, M, ORD)

  # Weighted sum over the causal combinations, P(C_i = a, C_j = b) = pi_i^a (1-pi_i)^(1-a) ...
  M_out <- outer(PI, PI) * Cov_M + outer(PI, 1 - PI) * Cov_M1 +
           outer(1 - PI, PI) * t(Cov_M1) + outer(1 - PI, 1 - PI) * Cov_M3

  # The diagonal is the mixture's own variance (the sum above would treat C_i as
  # independent of itself), by the law of total variance:
  # (1 - pi) Var_0 + pi Var_1 + pi (1 - pi) (E_1 - E_0)^2, with centered variances
  E1 <- B1[1L, ]
  E0 <- b0[1L]
  diag(M_out) <- (1 - PI) * tab$Var[n + 1L] + PI * tab$Var[seq_len(n)] +
                 PI * (1 - PI) * (E1 - E0)^2

  M_out
}


#' Get Covariance Matrix Under Specified Hypothesis
#'
#' @description
#' Wrapper function to compute the covariance matrix of transformed Z-scores
#' under either the null hypothesis (H0) or alternative hypothesis (H1) using
#' the sparse mixture approximation.
#'
#' @param g Function, transformation to apply to Z-scores
#' @param MU Numeric vector, effect sizes (B-star values)
#' @param PI Numeric vector, variant-importance scores
#' @param M Numeric matrix, correlation matrix of null Z-scores
#' @param hypo Character string, "H0" for null or "H1" for alternative
#'
#' @return Numeric matrix, covariance matrix under the specified hypothesis
#'
#' @details
#' Under H0 (null hypothesis):
#' \itemize{
#'   \item All variants are assumed to be non-causal
#'   \item Equivalent to setting MU = 0 and PI = 0 for all variants
#'   \item Covariance depends only on the correlation structure M
#' }
#'
#' Under H1 (alternative hypothesis):
#' \itemize{
#'   \item Uses the actual MU and PI values
#'   \item Accounts for sparse signals and their correlations
#'   \item Covariance depends on M, MU, and PI
#' }
#'
#' This function is used to compute the appropriate covariance matrix for
#' different weight calculation strategies:
#' \itemize{
#'   \item Sigma_BE (Best under H0): Uses H0 covariance
#'   \item Sigma_APE (Average Power under Ensemble): Uses H1 covariance
#' }
#'
#' @section Computational Complexity:
#' Same as CovMT_mix: O(4 * n^2 * K * n_eval) where n = length(MU)
#'
#' @examples
#' \dontrun{
#' # GFisher transformation
#' g_GFisher <- function(x, df = 2) {
#'   qchisq(log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
#'          df = df, lower.tail = FALSE, log.p = TRUE)
#' }
#' M <- diag(3)
#' MU <- c(1.0, 0.5, 0.8)
#' PI <- c(0.01, 0.05, 0.02)
#'
#' # Under null: all variants non-causal
#' Sigma_H0 <- get_Sigma(function(x) g_GFisher(x, 2), MU, PI, M, hypo = "H0")
#'
#' # Under alternative: use actual MU and PI
#' Sigma_H1 <- get_Sigma(function(x) g_GFisher(x, 2), MU, PI, M, hypo = "H1")
#' }
#'
#' @keywords internal
#' @noRd
get_Sigma <- function(g, MU, PI, M, hypo = "H1") {
  n <- length(MU)

  if (hypo == "H1") {
    # Under alternative: use actual MU and PI
    Sigma <- CovMT_mix(g, MU, PI, M)
  } else if (hypo == "H0") {
    # Under null: all variants are non-causal (MU = 0, PI = 0)
    Sigma <- CovMT_mix(g, rep(0, n), rep(0, n), M)
  } else {
    stop("hypo must be either 'H0' or 'H1'")
  }

  return(Sigma)
}


#################### Core Mathematical Functions ####################

#' Hermite Polynomials (Probabilist's Version)
#'
#' @description
#' Computes probabilist's Hermite polynomials of orders 1 through 8, shifted
#' by location parameter mu. These polynomials are orthogonal with respect to
#' the standard normal distribution and are used in Hermite series expansions
#' for approximating expectations involving transformed normal random variables.
#'
#' @param x Numeric vector of evaluation points
#' @param mu Numeric scalar, shift parameter (typically the mean)
#' @param ord Integer, order of the Hermite polynomial (a positive integer)
#'
#' @return Numeric vector of Hermite polynomial values evaluated at x
#'
#' @details
#' The probabilist's Hermite polynomials \eqn{He_n(x)} are defined such that:
#' \deqn{E[He_n(Z)] = 0} for \eqn{Z \sim N(0,1)}
#' \deqn{E[He_n(Z)He_m(Z)] = n! \delta_{nm}}
#'
#' This function implements the shifted versions:
#' \deqn{He_n(x - \mu)}
#'
#' The first 8 polynomials are:
#' \itemize{
#'   \item \eqn{He_1(x-\mu) = x - \mu}
#'   \item \eqn{He_2(x-\mu) = (x-\mu)^2 - 1}
#'   \item \eqn{He_3(x-\mu) = (x-\mu)^3 - 3(x-\mu)}
#'   \item \eqn{He_4(x-\mu) = (x-\mu)^4 - 6(x-\mu)^2 + 3}
#'   \item \eqn{He_5(x-\mu) = (x-\mu)^5 - 10(x-\mu)^3 + 15(x-\mu)}
#'   \item \eqn{He_6(x-\mu) = (x-\mu)^6 - 15(x-\mu)^4 + 45(x-\mu)^2 - 15}
#'   \item \eqn{He_7(x-\mu) = (x-\mu)^7 - 21(x-\mu)^5 + 105(x-\mu)^3 - 105(x-\mu)}
#'   \item \eqn{He_8(x-\mu) = (x-\mu)^8 - 28(x-\mu)^6 + 210(x-\mu)^4 - 420(x-\mu)^2 + 105}
#' }
#'
#' These polynomials are used in the Hermite series expansion to compute
#' covariances of transformed correlated normal random variables.
#'
#' @section Computational Complexity:
#' O(length(x)) - polynomial evaluation is vectorized
#'
#' @section Mathematical Background:
#' The Hermite expansion is used to approximate:
#' \deqn{Cov[g(X), g(Y)]}
#' where X and Y are correlated normal random variables. The expansion uses
#' the fact that if X and Y have correlation \eqn{\rho}, then:
#' \deqn{Cov[He_n(X), He_m(Y)] = n! \rho^n \delta_{nm}}
#'
#' @examples
#' \dontrun{
#' # Evaluate Hermite polynomials at standard normal quantiles
#' x <- seq(-3, 3, by = 0.5)
#' mu <- 0
#'
#' # First few orders
#' he1 <- hermite(x, mu, 1)
#' he2 <- hermite(x, mu, 2)
#' he3 <- hermite(x, mu, 3)
#'
#' # Verify orthogonality property (approximately)
#' z <- rnorm(10000)
#' mean(hermite(z, 0, 2) * hermite(z, 0, 3))  # Should be close to 0
#' }
#'
#' @references
#' Abramowitz, M. and Stegun, I. A. (1972). Handbook of Mathematical Functions.
#' Dover Publications.
#'
#' @keywords internal
#' @noRd
hermite <- function(x, mu, ord) {
  if (!is.numeric(ord) || length(ord) != 1L || !is.finite(ord) || ord < 1 ||
      ord != round(ord)) {
    stop("Hermite polynomial order must be a positive integer")
  }
  # Three-term recurrence He_{k+1}(y) = y He_k(y) - k He_{k-1}(y) from He_0 = 1, He_1 = y
  y <- x - mu
  h_prev <- rep(1, length(y))
  h <- y
  if (ord >= 2) {
    for (k in seq_len(ord - 1)) {
      h_next <- y * h - k * h_prev
      h_prev <- h
      h <- h_next
    }
  }
  h
}


#' Calculate E[g(X)^p] for Normal X
#'
#' @description
#' Computes the expectation of g(X) raised to power p, where X follows a
#' normal distribution with mean mu and standard deviation sigma. The
#' expectation is approximated using numerical integration.
#'
#' @param g Function to apply to X before taking expectation
#' @param mu Numeric scalar, mean of the normal distribution
#' @param p Numeric scalar, power to raise g(X) to (typically 1 or 2)
#' @param sigma Numeric scalar, standard deviation (default = 1)
#'
#' @return Numeric scalar, the approximate value of E[g(X)^p]
#'
#' @details
#' This function computes:
#' \deqn{E[g(X)^p] = \int g(x)^p f(x) dx}
#' where \eqn{f(x)} is the normal density with mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}.
#'
#' The integration uses the composite Gauss-Legendre rule of
#' \code{optimal_weights_quadrature.R}: the standardized range \eqn{[-12, 12]} in panels
#' of width 0.1, cut where the argument of g is zero (the kink of the two-sided
#' transforms).
#'
#' Common use cases:
#' \itemize{
#'   \item p = 1: Expected value of transformed variable
#'   \item p = 2: Expected squared value (used for variance calculations)
#' }
#'
#' @section Computational Complexity:
#' O(n_eval) with n_eval = 2,410 evaluations of g (one vectorized call)
#'
#' @section Numerical Considerations:
#' \itemize{
#'   \item A fixed composite Gauss-Legendre rule, deterministic and exact for polynomials
#'   \item Validated for the Fisher and polynomial transforms; a transform with a kink away
#'     from zero, a cusp, a discontinuity or rapid growth needs its own convergence check
#'   \item A non-finite value of g, or a non-finite moment, is an error
#' }
#'
#' @examples
#' \dontrun{
#' # For X ~ N(0, 1), E[X] should be 0
#' g_identity <- function(x) x
#' E_gX_p(g_identity, mu = 0, p = 1, sigma = 1)  # Should be close to 0
#'
#' # For X ~ N(0, 1), E[X^2] should be 1
#' E_gX_p(g_identity, mu = 0, p = 2, sigma = 1)  # Should be close to 1
#'
#' # For a transformation like squaring
#' g_square <- function(x) x^2
#' E_gX_p(g_square, mu = 0, p = 1, sigma = 1)  # E[X^2] = 1
#' }
#'
#' @keywords internal
#' @noRd
E_gX_p <- function(g, mu, p, sigma = 1) {
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be a single finite positive number")
  }
  # Nodes on the standardized scale, cut where mu + sigma * u = 0
  nd <- .normal_nodes(cut = -mu / sigma)
  x <- mu + sigma * nd$x
  gx <- g(x)
  if (length(gx) != length(x)) {
    stop("The transformation g must return one value per input (", length(x),
         " inputs, ", length(gx), " values).", call. = FALSE)
  }
  bad <- which(!is.finite(gx))
  if (length(bad)) {
    stop("The transformation g returned a non-finite value at x = ", format(x[bad[1L]]),
         " (mu = ", format(mu), ", sigma = ", format(sigma), ").", call. = FALSE)
  }
  value <- sum(gx^p * stats::dnorm(nd$x) * nd$w)
  if (!is.finite(value)) {
    stop("E[g(X)^", p, "] is not finite (mu = ", format(mu), ", sigma = ", format(sigma),
         "): the transformation grows too fast for the fixed quadrature domain.",
         call. = FALSE)
  }
  value
}


#' Calculate E[T_i] for Mixture Distribution
#'
#' @description
#' Computes the expectation of T_i = g(Z_i), where Z_i follows a mixture
#' distribution: Z_i = Z_0i + mu*C_i, with Z_0i ~ N(0,1) and C_i ~ Bernoulli(pi).
#'
#' @param g Function, the transformation to apply
#' @param mu Numeric scalar, the shift parameter (effect size)
#' @param pi Numeric scalar, the mixing probability P(C_i = 1)
#'
#' @return Numeric scalar, E[g(Z_i)]
#'
#' @details
#' This function calculates:
#' \deqn{E[T_i] = E[g(Z_i)] = \pi E[g(Z_0 + \mu)] + (1-\pi) E[g(Z_0)]}
#'
#' The mixture distribution arises in the GLOW framework where:
#' \itemize{
#'   \item Z_0i is the null Z-score (no effect)
#'   \item C_i is a binary indicator of whether variant i is causal
#'   \item mu is the effect size when variant is causal
#'   \item pi is the prior probability that variant is causal
#' }
#'
#' This represents a sparse signal model where only a fraction pi of variants
#' have non-zero effects.
#'
#' @section Computational Complexity:
#' O(n_eval) where n_eval is the cost of two calls to E_gX_p()
#'
#' @section Assumptions:
#' \itemize{
#'   \item 0 <= pi <= 1
#'   \item g is a well-behaved function
#' }
#'
#' @examples
#' \dontrun{
#' # Identity transformation
#' g_identity <- function(x) x
#' E_T_mix(g_identity, mu = 2, pi = 0.5)  # Should be close to 1
#'
#' # Square transformation
#' g_square <- function(x) x^2
#' E_T_mix(g_square, mu = 2, pi = 0.1)
#' }
#'
#' @keywords internal
#' @noRd
E_T_mix <- function(g, mu, pi) {
  # Expectation under mixture:
  # E[g(Z)] = pi * E[g(Z_0 + mu)] + (1 - pi) * E[g(Z_0)]
  # where Z_0 ~ N(0, 1)
  result <- pi * E_gX_p(g, mu, 1) + (1 - pi) * E_gX_p(g, 0, 1)
  return(result)
}


#' Calculate Var[T_i] for Mixture Distribution
#'
#' @description
#' Computes the variance of T_i = g(Z_i), where Z_i follows a mixture
#' distribution: Z_i = Z_0i + mu*C_i, with Z_0i ~ N(0,1) and C_i ~ Bernoulli(pi).
#'
#' @param g Function, the transformation to apply
#' @param mu Numeric scalar, the shift parameter (effect size)
#' @param pi Numeric scalar, the mixing probability P(C_i = 1)
#'
#' @return Numeric scalar, Var[g(Z_i)]
#'
#' @details
#' This function calculates:
#' \deqn{Var[T_i] = E[T_i^2] - (E[T_i])^2}
#' where:
#' \deqn{E[T_i^2] = \pi E[g(Z_0 + \mu)^2] + (1-\pi) E[g(Z_0)^2]}
#'
#' The variance accounts for both:
#' \itemize{
#'   \item Variability within each component (causal vs non-causal)
#'   \item Variability between components (mixture variability)
#' }
#'
#' @section Computational Complexity:
#' O(n_eval) where n_eval is the cost of three calls to E_gX_p() plus one
#' call to E_T_mix()
#'
#' @section Assumptions:
#' \itemize{
#'   \item 0 <= pi <= 1
#'   \item g is a well-behaved function
#' }
#'
#' @examples
#' \dontrun{
#' # Identity transformation
#' g_identity <- function(x) x
#' Var_T_mix(g_identity, mu = 2, pi = 0.5)
#'
#' # For pi = 0 or pi = 1, should reduce to variance of single component
#' Var_T_mix(g_identity, mu = 2, pi = 0)  # Var of N(0,1)
#' Var_T_mix(g_identity, mu = 2, pi = 1)  # Var of N(2,1)
#' }
#'
#' @keywords internal
#' @noRd
Var_T_mix <- function(g, mu, pi) {
  # Law of total variance over the causal indicator C ~ Bernoulli(pi):
  # Var[T] = pi Var_1 + (1 - pi) Var_0 + pi (1 - pi) (E_1 - E_0)^2, where E_a and Var_a
  # are the mean and the (centered) variance of g(a mu + Z_0), so that a large mean does
  # not cancel against a large second moment
  m1 <- .hermite_moments(g, mu, 1, 0L)
  m0 <- .hermite_moments(g, 0, 1, 0L)
  pi * m1$Var + (1 - pi) * m0$Var + pi * (1 - pi) * (m1$b[1L] - m0$b[1L])^2
}


#' Validate the orders of the Hermite expansion
#'
#' @param ORD Integer vector of positive integers (the orders summed).
#' @return The orders as a sorted integer vector without duplicates.
#' @keywords internal
#' @noRd
.check_ord <- function(ORD) {
  if (!is.numeric(ORD) || length(ORD) == 0L || any(!is.finite(ORD)) || any(ORD < 1) ||
      any(ORD != round(ORD))) {
    stop("ORD must be a vector of positive integers (the orders of the Hermite expansion)")
  }
  sort(unique(as.integer(ORD)))
}

#' The moments behind the four optimal-weight schemes, from one coefficient table
#'
#' @description
#' For a non-Burden transformation \code{g}, computes everything \code{Optimal_Weights_M()}
#' needs: the mean-difference vectors and the covariance matrices of the normal
#' approximation (\code{_N}) and of the sparse mixture (\code{_sparse}), under the null
#' (BE) and the alternative (APE). The Hermite coefficients of every distinct
#' (mean, scale) pair among the null pair (0, 1), the normal-approximation pairs
#' \eqn{(MMU_i, s_i)} and the sparse pairs \eqn{(Bstar_i, 1)} are computed once, and four
#' Mehler sums follow. The values equal those of \code{get_r_tilde()},
#' \code{CovM_gXgY()}, \code{get_r()} and \code{get_Sigma()} called one by one.
#'
#' @param g Function, the transformation.
#' @param Bstar,PI,M As in \code{Optimal_Weights_M()}.
#' @param ORD Integer vector, the orders of the Hermite expansion summed (default 1:64).
#' @return A list: \code{r_N}, \code{Sigma_BE_N}, \code{Sigma_APE_N}, \code{r_sparse},
#'   \code{Sigma_BE_sparse}, \code{Sigma_APE_sparse}.
#' @keywords internal
#' @noRd
.optimal_weight_moments <- function(g, Bstar, PI, M, ORD = 1:64) {
  n <- length(Bstar)
  ORD <- .check_ord(ORD)
  K <- max(ORD)

  # The normal approximation: Z~ ~ N(MMU, M + MVM), handled as MMU_i + s_i U_i with
  # U ~ N(0, R), R the correlation matrix of M + MVM (the expansion needs a correlation)
  MU <- Bstar * PI
  V <- Bstar^2 * PI * (1 - PI)
  MMU <- as.vector(M %*% MU)
  Sigma_tilde <- M + t(M * V) %*% M
  SD_tilde <- sqrt(diag(Sigma_tilde))
  R_tilde <- stats::cov2cor(Sigma_tilde)

  # One table: column 1 the null pair, then the normal-approximation pairs, then the
  # sparse pairs (each distinct pair integrated once)
  tab <- .hermite_table(g, c(0, MMU, Bstar), c(1, SD_tilde, rep(1, n)), K)
  b0 <- tab$b[, 1L]
  B_N <- tab$b[, 1L + seq_len(n), drop = FALSE]
  B_S <- tab$b[, 1L + n + seq_len(n), drop = FALSE]
  B0 <- matrix(b0, K + 1L, n)

  # Normal approximation: r_i = E[g(Z~_i)] - E[g(Z_0)]; the null covariance (truncated
  # series in M, the same matrix that the sparse scheme uses as its null component) and
  # the alternative covariance (the series in R with the stretched coefficients)
  r_N <- B_N[1L, ] - b0[1L]
  Sigma_null <- .mehler_sum(B0, B0, M, ORD)
  Sigma_APE_N <- .mehler_sum(B_N, B_N, R_tilde, ORD)

  # Sparse mixture: Z_i = Z_0i + Bstar_i C_i, C_i ~ Bernoulli(PI_i); the off-diagonals
  # from the four causal combinations, the diagonal the mixture's own variance (the law
  # of total variance with centered variances), as in CovMT_mix()
  Cov_11 <- .mehler_sum(B_S, B_S, M, ORD)
  Cov_10 <- .mehler_sum(B_S, B0, M, ORD)
  Sigma_APE_sparse <- outer(PI, PI) * Cov_11 + outer(PI, 1 - PI) * Cov_10 +
                      outer(1 - PI, PI) * t(Cov_10) + outer(1 - PI, 1 - PI) * Sigma_null
  E1 <- B_S[1L, ]
  E0 <- b0[1L]
  Var1 <- tab$Var[1L + n + seq_len(n)]
  Var0 <- tab$Var[1L]
  diag(Sigma_APE_sparse) <- (1 - PI) * Var0 + PI * Var1 + PI * (1 - PI) * (E1 - E0)^2
  Sigma_BE_sparse <- Sigma_null
  diag(Sigma_BE_sparse) <- Var0
  r_sparse <- PI * (E1 - E0)

  list(r_N = r_N, Sigma_BE_N = Sigma_null, Sigma_APE_N = Sigma_APE_N,
       r_sparse = r_sparse, Sigma_BE_sparse = Sigma_BE_sparse,
       Sigma_APE_sparse = Sigma_APE_sparse)
}


#################### Robustness Helper Functions ####################

#' Safely Invert Matrix with Condition Number Check
#'
#' @description
#' Performs safe matrix inversion with automatic detection and correction of
#' near-singular matrices using nearPD adjustment.
#'
#' @param Sigma Numeric matrix to invert
#' @param context Character string describing the context (for error messages)
#' @param threshold Numeric, condition number threshold above which nearPD is used (default 1e10)
#'
#' @return Inverted matrix
#'
#' @details
#' This function checks the condition number of Sigma before inversion. If the
#' condition number exceeds the threshold, the matrix is adjusted using
#' Matrix::nearPD() to ensure numerical stability.
#'
#' Condition number kappa(Sigma) = max(eigenvalues) / min(eigenvalues)
#' High condition numbers (> 1e10) indicate near-singularity.
#'
#' @keywords internal
#' @noRd
safe_matrix_inverse <- function(Sigma, context = "matrix inversion", threshold = 1e10) {
  # Check condition number
  kappa_val <- kappa(Sigma)

  if (kappa_val > threshold) {
    warning("Matrix is near-singular in ", context,
            " (condition number = ", sprintf("%.2e", kappa_val), "). ",
            "Using nearPD adjustment for numerical stability.",
            call. = FALSE)
    Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE)$mat)
    kappa_val <- kappa(Sigma)
  }

  # Attempt inversion with informative error
  tryCatch(
    solve(Sigma),
    error = function(e) {
      stop("Matrix inversion failed in ", context, ".\n",
           "  Condition number: ", sprintf("%.2e", kappa_val), "\n",
           "  Matrix dimensions: ", nrow(Sigma), " x ", ncol(Sigma), "\n",
           "  Minimum eigenvalue: ", sprintf("%.2e", min(eigen(Sigma, only.values = TRUE)$values)), "\n",
           "  Original error: ", e$message,
           call. = FALSE)
    }
  )
}


#' Validate Numeric Inputs for NA/Inf
#'
#' @description
#' Validates that numeric inputs do not contain NA or infinite values, with
#' optional checks for non-negativity and positivity.
#'
#' @param x Numeric vector or matrix to validate
#' @param name Character string, name of the variable (for error messages)
#' @param allow_negative Logical, whether negative values are allowed (default TRUE)
#' @param allow_zero Logical, whether zero values are allowed (default TRUE)
#'
#' @return Invisible TRUE if all checks pass
#'
#' @keywords internal
#' @noRd
validate_numeric_input <- function(x, name, allow_negative = TRUE, allow_zero = TRUE) {
  if (any(is.na(x))) {
    stop(name, " contains NA values. Please remove or impute missing data.",
         call. = FALSE)
  }
  if (any(is.infinite(x))) {
    stop(name, " contains infinite values.",
         call. = FALSE)
  }
  if (!allow_negative && any(x < 0)) {
    stop(name, " must be non-negative. Found negative values.",
         call. = FALSE)
  }
  if (!allow_zero && any(x == 0)) {
    stop(name, " must be positive (non-zero). Found zero values.",
         call. = FALSE)
  }
  invisible(TRUE)
}


#' Check if Matrix is Positive Semi-Definite
#'
#' @description
#' Validates that a matrix is positive semi-definite by checking eigenvalues.
#'
#' @param M Numeric matrix to check
#' @param name Character string, name of the matrix (for error/warning messages)
#' @param tol Numeric, tolerance for considering eigenvalues as zero (default 1e-10)
#'
#' @return Invisible TRUE if matrix is positive semi-definite
#'
#' @keywords internal
#' @noRd
check_positive_definite <- function(M, name = "matrix", tol = 1e-10) {
  eigenvalues <- eigen(M, only.values = TRUE)$values
  min_eigen <- min(eigenvalues)

  if (min_eigen < -tol) {
    stop(name, " must be positive semi-definite. ",
         "Minimum eigenvalue: ", sprintf("%.2e", min_eigen),
         call. = FALSE)
  }

  if (min_eigen < tol) {
    warning(name, " is nearly singular. ",
            "Minimum eigenvalue: ", sprintf("%.2e", min_eigen),
            call. = FALSE)
  }

  invisible(TRUE)
}




########## Cache Infrastructure ##########

# Package-level cache for H0 covariance matrices
# This is a private environment not exported to users
.glow_cache <- new.env(parent = emptyenv())

#' Clear GLOW Covariance Cache
#'
#' @description
#' Clears the internal cache used for storing H0 covariance matrices.
#' This can be useful to free memory or force recomputation.
#'
#' @return Invisible NULL
#'
#' @details
#' The GLOWr package caches H0 (null hypothesis) covariance matrices to improve
#' performance when \code{Optimal_Weights_M()} is called multiple times with the
#' same correlation structure. This function clears that cache.
#'
#' Caching can be controlled via \code{options(GLOWr.use_cache = TRUE/FALSE)}.
#' The default is TRUE (caching enabled).
#'
#' @examples
#' # Clear the cache
#' clear_glow_cache()
#'
#' # Disable caching entirely
#' options(GLOWr.use_cache = FALSE)
#'
#' @export
clear_glow_cache <- function() {
  rm(list = ls(envir = .glow_cache), envir = .glow_cache)
  invisible(NULL)
}

#' Get GLOW Cache Information
#'
#' @description
#' Returns information about the current state of the GLOW covariance cache.
#'
#' @return A list with elements:
#' \describe{
#'   \item{enabled}{Logical, whether caching is currently enabled}
#'   \item{cached_items}{Integer, number of items currently in cache}
#'   \item{cache_keys}{Character vector of cache keys (for debugging)}
#' }
#'
#' @examples
#' # Check cache status
#' glow_cache_info()
#'
#' @export
glow_cache_info <- function() {
  list(
    enabled = getOption("GLOWr.use_cache", default = TRUE),
    cached_items = length(ls(envir = .glow_cache)),
    cache_keys = ls(envir = .glow_cache)
  )
}
