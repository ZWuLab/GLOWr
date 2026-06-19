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
#'    under the alternative with means MMU and variances M + MVM, where:
#'    \itemize{
#'      \item MMU = M \%*\% MU (marginal means)
#'      \item MVM = t(M * V) \%*\% M (marginal variance contribution from random effects)
#'    }
#'
#' 2. \emph{Sparse Mixture Approximation}: Explicitly models the sparse signal
#'    structure Z_i = Z_0i + Bstar_i * C_i where C_i ~ Bernoulli(PI_i).
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
#'   \item Other tests: O(4 * n^2 * K * n_eval + n^3) where:
#'     \itemize{
#'       \item n = number of variants
#'       \item K = number of Hermite polynomial orders (8)
#'       \item n_eval = cost per numerical integration
#'       \item Factor of 4 from normal + sparse approximations, each under H0 and H1
#'     }
#' }
#'
#' @section Assumptions:
#' \itemize{
#'   \item Bstar, PI, and M have compatible dimensions
#'   \item M is a valid correlation matrix (positive definite, diagonal = 1)
#'   \item 0 <= PI_i <= 1 for all i
#'   \item Bstar values are positive (standardized effect sizes)
#'   \item For non-Burden tests, g is a smooth function suitable for Hermite expansion
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
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
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

    ## 1. Normal Approximation
    # Calculate r_tilde: mean differences under normal approximation
    r_tilde <- get_r_tilde(g, MU = MMU, SD = diag(M + MVM))

    # Covariance matrices under normal approximation
    # BE: Under H0 (no signal), means = 0
    Sigma_BE_tilde <- CovM_gXgY(g, MU1 = rep(0, length(MMU)), MU2 = rep(0, length(MMU)), M = M)
    # APE: Under H1 (with signal), means = MMU, covariance = M + MVM
    Sigma_APE_tilde <- CovM_gXgY(g, MU1 = MMU, MU2 = MMU, M = M + MVM)

    # Calculate weights using normal approximation
    wts_BE_N <- get_wts(Sigma_BE_tilde, r_tilde, is.posi.wts)
    wts_APE_N <- get_wts(Sigma_APE_tilde, r_tilde, is.posi.wts)

    ## 2. Sparse Mixture Approximation
    # Calculate r: mean differences under sparse mixture model
    r <- get_r(g, Bstar, PI)

    # Covariance matrices under sparse mixture approximation
    # BE: Under H0 (no signal)
    Sigma_BE <- get_Sigma(g, Bstar, PI, M, hypo = "H0")
    # APE: Under H1 (with signal)
    Sigma_APE <- get_Sigma(g, Bstar, PI, M, hypo = "H1")

    # Calculate weights using sparse approximation
    wts_BE_sparse <- get_wts(Sigma_BE, r, is.posi.wts)
    wts_APE_sparse <- get_wts(Sigma_APE, r, is.posi.wts)

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
  n <- length(MU)
  # For each i, compute E[g(Z_i)] - E[g(Z_0)]
  # where Z_i ~ N(MU[i], SD[i]^2) and Z_0 ~ N(0, 1)
  r_tilde <- sapply(1:n, function(i) {
    E_gX_p(g, mu = MU[i], p = 1, sigma = SD[i]) - E_gX_p(g, mu = 0, p = 1, sigma = 1)
  })
  return(r_tilde)
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
#' expansion. This is the core function for computing covariances under
#' the normal approximation.
#'
#' @param g Function, transformation to apply to normal random variables
#' @param MU1 Numeric vector, means of X components
#' @param MU2 Numeric vector, means of Y components
#' @param M Numeric matrix, correlation matrix between X and Y (must be positive definite)
#' @param ORD Integer vector, orders of Hermite polynomials to use (default 1:8)
#'
#' @return Numeric matrix, covariance matrix where element (i,j) is Cov[g(X_i), g(Y_j)]
#'
#' @details
#' This function computes:
#' \deqn{Cov[g(X_i), g(Y_j)]}
#' where X = Z + MU1 ~ MVN(MU1, M) and Y = Z + MU2 ~ MVN(MU2, M), with
#' Z ~ MVN(0, M) being the shared component.
#'
#' The computation uses the Hermite expansion:
#' \deqn{g(X) \approx \sum_{k=1}^{K} a_k He_k(X)}
#' where He_k are probabilist's Hermite polynomials.
#'
#' The key property used is:
#' \deqn{Cov[He_k(X_i), He_k(Y_j)] = k! M_{ij}^k}
#'
#' The algorithm:
#' \enumerate{
#'   \item Compute Hermite coefficients for g(X_i) around each mean MU1[i]
#'   \item Compute Hermite coefficients for g(Y_j) around each mean MU2[j]
#'   \item Sum contributions from each order k: coef1[k,i] * coef2[k,j] * M_ij^k / k!
#' }
#'
#' @section Computational Complexity:
#' O(n^2 * K * n_eval) for correlated variants, O(n * K * n_eval) for independent variants:
#' \itemize{
#'   \item n = length(MU1) = length(MU2)
#'   \item K = length(ORD) (typically 8)
#'   \item n_eval = cost per numerical integration
#'   \item Fast path optimization when M is identity matrix
#' }
#'
#' @section Numerical Considerations:
#' \itemize{
#'   \item Higher order terms (ORD > 8) may be unstable
#'   \item Integration performed over +/- 8 standard deviations
#'   \item Assumes g is sufficiently smooth for Hermite expansion
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
CovM_gXgY <- function(g, MU1, MU2, M, ORD = 1:8) {
  # Check that M is positive semi-definite
  check_positive_definite(M, "Correlation matrix M")

  n <- length(MU1)

  # Fast path for independent variants (identity correlation matrix)
  # When M is identity, off-diagonal covariances are zero
  # This avoids expensive Hermite expansion for off-diagonals
  is_identity <- (nrow(M) == ncol(M)) &&
                 all(abs(diag(M) - 1) < 1e-10) &&
                 all(abs(M[row(M) != col(M)]) < 1e-10)

  if (is_identity) {
    # Only compute diagonal variances
    diag_vals <- sapply(seq_len(n), function(i) {
      # Compute Hermite coefficients for g(X) around MU1[i]
      coef_i <- sapply(ORD, function(ord) {
        integrand <- function(x) {
          g(x) * dnorm(x - MU1[i]) * hermite(x, MU1[i], ord)
        }
        integrate(integrand, MU1[i] - 8, MU1[i] + 8)$value
      })

      # Variance = sum of squared coefficients / factorial(ord)
      # This uses the orthogonality property: E[He_k(X)^2] = k!
      sum(coef_i^2 / factorial(ORD))
    })

    # Handle single variant case (diag() on scalar creates 1x1 matrix)
    M_out <- diag(diag_vals, nrow = n, ncol = n)

    return(M_out)
  }

  # Continue with regular computation for correlated variants

  # Compute Hermite coefficients for g(X) around each mean MU1[i]
  # coef1[k, i] = integral of g(x) * phi(x - MU1[i]) * He_k(x - MU1[i]) dx
  coef1 <- sapply(MU1, function(mu) {
    sapply(ORD, function(ord) {
      # Integrand: g(x) * density * Hermite polynomial
      integrand <- function(x) {
        g(x) * dnorm(x - mu) * hermite(x, mu, ord)
      }
      # Integrate over +/- 8 standard deviations from mu
      integrate(integrand, mu - 8, mu + 8)$value
    })
  })

  # Compute Hermite coefficients for g(Y) around each mean MU2[j]
  coef2 <- sapply(MU2, function(mu) {
    sapply(ORD, function(ord) {
      integrand <- function(x) {
        g(x) * dnorm(x - mu) * hermite(x, mu, ord)
      }
      integrate(integrand, mu - 8, mu + 8)$value
    })
  })

  # Initialize output covariance matrix
  M_out <- matrix(0, ncol = n, nrow = n)

  # Sum over Hermite polynomial orders
  for (i in seq_along(ORD)) {
    ord <- ORD[i]
    # Contribution from order 'ord':
    # Cov[He_ord(X_i), He_ord(Y_j)] = ord! * M_ij^ord
    # So Cov[g(X_i), g(Y_j)] gets contribution: coef1[ord, i] * coef2[ord, j] * M_ij^ord / ord!
    M_out <- M_out + as.matrix(coef1[ord, ]) %*% t(as.matrix(coef2[ord, ])) * M^ord / factorial(ord)
  }

  return(M_out)
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
#' @param ORD Integer vector, orders of Hermite polynomials to use (default 1:8)
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
CovMT_mix <- function(g, MU, PI, M, ORD = 1:8) {
  # Check that M is positive semi-definite
  check_positive_definite(M, "Correlation matrix M")

  n <- length(MU)

  # Create probability matrices for the four combinations
  PI_M <- PI %*% t(PI)             # P(C_i=1, C_j=1) = pi_i * pi_j
  PI_M1 <- PI %*% t(1 - PI)        # P(C_i=1, C_j=0) = pi_i * (1-pi_j)
  PI_M2 <- t(PI_M1)                # P(C_i=0, C_j=1) = (1-pi_i) * pi_j
  PI_M3 <- (1 - PI) %*% t(1 - PI)  # P(C_i=0, C_j=0) = (1-pi_i) * (1-pi_j)

  # Compute covariance matrices for each combination of causal status
  # Cov_M: Both i and j are causal (means = MU for both)
  Cov_M <- CovM_gXgY(g, MU, MU, M, ORD = 1:8)

  # Cov_M1: i is causal (mean = MU), j is null (mean = 0)
  Cov_M1 <- CovM_gXgY(g, MU, rep(0, n), M, ORD = 1:8)

  # Cov_M2: i is null (mean = 0), j is causal (mean = MU)
  Cov_M2 <- t(Cov_M1)

  # Cov_M3: Both i and j are null (means = 0 for both)
  # This is expensive to compute but identical for all calls with same M and g
  # Use caching if enabled
  use_cache <- getOption("GLOWr.use_cache", default = TRUE)

  if (use_cache) {
    # Generate cache key based on g function body and M
    cache_key <- tryCatch({
      digest::digest(list(
        g_body = body(g),
        g_formals = formals(g),
        M = M,
        type = "H0_covariance"
      ), algo = "xxhash64")
    }, error = function(e) {
      warning("Cache key generation failed: ", e$message,
              ". Falling back to no caching.", call. = FALSE)
      NULL
    })

    # Try to use cache
    if (!is.null(cache_key)) {
      if (exists(cache_key, envir = .glow_cache)) {
        # Cache hit
        Cov_M3 <- get(cache_key, envir = .glow_cache)
      } else {
        # Cache miss - compute and store
        Cov_M3 <- tryCatch({
          result <- CovM_gXgY(g, rep(0, n), rep(0, n), M, ORD = 1:8)
          assign(cache_key, result, envir = .glow_cache)
          result
        }, error = function(e) {
          warning("Caching failed: ", e$message,
                  ". Computing without caching.", call. = FALSE)
          CovM_gXgY(g, rep(0, n), rep(0, n), M, ORD = 1:8)
        })
      }
    } else {
      # Cache key generation failed, compute without caching
      Cov_M3 <- CovM_gXgY(g, rep(0, n), rep(0, n), M, ORD = 1:8)
    }
  } else {
    # Caching disabled by user
    Cov_M3 <- CovM_gXgY(g, rep(0, n), rep(0, n), M, ORD = 1:8)
  }

  # Weighted sum of covariance matrices
  M_out <- PI_M * Cov_M + PI_M1 * Cov_M1 + PI_M2 * Cov_M2 + PI_M3 * Cov_M3

  # Replace diagonal elements with exact variance calculations
  # This avoids numerical issues from Hermite expansion on the diagonal
  diag(M_out) <- sapply(1:n, function(i) Var_T_mix(g, MU[i], PI[i]))

  return(M_out)
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
#' @param ord Integer, order of the Hermite polynomial (1 through 8)
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
  # Compute shifted variable
  y <- x - mu

  # Return appropriate Hermite polynomial based on order
  if (ord == 1) {
    return(y)
  } else if (ord == 2) {
    return(y^2 - 1)
  } else if (ord == 3) {
    return(y^3 - 3*y)
  } else if (ord == 4) {
    return(y^4 - 6*y^2 + 3)
  } else if (ord == 5) {
    return(y^5 - 10*y^3 + 15*y)
  } else if (ord == 6) {
    return(y^6 - 15*y^4 + 45*y^2 - 15)
  } else if (ord == 7) {
    return(y^7 - 21*y^5 + 105*y^3 - 105*y)
  } else if (ord == 8) {
    return(y^8 - 28*y^6 + 210*y^4 - 420*y^2 + 105)
  } else {
    stop("Hermite polynomial order must be between 1 and 8")
  }
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
#' The integration is performed over the range \eqn{[\mu - 8\sigma, \mu + 8\sigma]},
#' which captures virtually all probability mass (> 99.9999%) for a normal
#' distribution.
#'
#' Common use cases:
#' \itemize{
#'   \item p = 1: Expected value of transformed variable
#'   \item p = 2: Expected squared value (used for variance calculations)
#' }
#'
#' @section Computational Complexity:
#' O(n_eval) where n_eval is the number of function evaluations needed by
#' the adaptive quadrature algorithm (typically 100-500)
#'
#' @section Numerical Considerations:
#' \itemize{
#'   \item Uses adaptive quadrature via \code{integrate()}
#'   \item Integration bounds chosen to capture > 99.9999% of normal mass
#'   \item May have issues if g(x) has sharp discontinuities
#'   \item Default tolerances from \code{integrate()} are used
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
  # Define the integrand: g(x)^p * f(x) where f is normal density
  # We integrate over the standardized scale and adjust the density
  integrand <- function(x) {
    g(x)^p * dnorm((x - mu) / sigma)
  }

  # Integration bounds: mu +/- 8*sigma captures > 99.9999% of mass
  lower <- mu - 8 * sigma
  upper <- mu + 8 * sigma

  # Perform numerical integration
  # Divide by sigma to account for the scaling in the density
  result <- integrate(integrand, lower, upper)$value / sigma

  return(result)
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
  # E[T^2] for the mixture
  E_T2 <- pi * E_gX_p(g, mu, 2) + (1 - pi) * E_gX_p(g, 0, 2)

  # E[T] for the mixture
  E_T <- E_T_mix(g, mu, pi)

  # Var[T] = E[T^2] - (E[T])^2
  variance <- E_T2 - E_T^2

  return(variance)
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
