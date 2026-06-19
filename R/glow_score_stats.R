# ==============================================================================
# Compute Marginal Score Statistics (Dispatcher)
# ==============================================================================
#
# Thin dispatcher selecting between standard and SPA-corrected score statistics.
#
# EXPORTED FUNCTIONS:
#   - compute_score_stats()  Dispatch to appropriate score stat function


#################### EXPORTED MAIN FUNCTIONS ####################

#' Compute Marginal Score Statistics
#'
#' Dispatches to the appropriate score statistic function based on trait type
#' and SPA preference. For binary traits, defaults to SPA-corrected statistics.
#' For continuous traits, uses standard score statistics.
#'
#' @param G Numeric matrix (n_samples x p_variants). Processed genotype matrix.
#'   Must not contain NA values.
#' @param null_model A \code{glow_null_model} object from
#'   \code{\link{fit_null_model}}.
#' @param use_spa Logical or NULL (default NULL). NULL = auto-detect (SPA for
#'   binary traits, standard for continuous). TRUE = force SPA (binary only).
#'   FALSE = force standard (any trait).
#' @param verbose Integer. Verbosity level (default 1). Set to 0 to suppress
#'   messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{Zscores}{Numeric vector (length p) of standardized marginal score
#'     statistics.}
#'   \item{M_Z}{Numeric matrix (p x p) correlation matrix of Z-scores.}
#'   \item{M_s}{Numeric matrix (p x p) covariance matrix of score statistics.}
#'   \item{s0}{Numeric scalar. Dispersion parameter (1 for binary traits).}
#' }
#'
#' @details
#' This is a thin dispatcher that unifies the API for score statistic
#' computation. It selects between \code{\link{getZ_marg_score}} and
#' \code{\link{getZ_marg_score_binary_SPA}} based on trait type (extracted from
#' the null model) and the \code{use_spa} parameter.
#'
#' The default behavior (\code{use_spa = NULL}) auto-detects: SPA for binary
#' traits (recommended for rare variants and unbalanced case-control), standard
#' for continuous traits. Users can override with explicit \code{TRUE} or
#' \code{FALSE}.
#'
#' All computation is delegated to the underlying functions; this function
#' performs only input validation and dispatch logic.
#'
#' @examples
#' # Binary trait (auto-selects SPA)
#' set.seed(42)
#' n <- 200; p <- 5
#' G <- matrix(rbinom(n * p, 2, 0.1), n, p)
#' X <- cbind(1, rnorm(n))
#' Y <- rbinom(n, 1, 0.3)
#' nm <- fit_null_model(X, Y, trait = "binary")
#' res <- compute_score_stats(G, nm)
#'
#' # Continuous trait (auto-selects standard)
#' Y_cont <- rnorm(n)
#' nm_cont <- fit_null_model(X, Y_cont, trait = "continuous")
#' res_cont <- compute_score_stats(G, nm_cont)
#'
#' @seealso
#' \code{\link{getZ_marg_score}} for standard score statistics,
#' \code{\link{getZ_marg_score_binary_SPA}} for SPA-corrected statistics,
#' \code{\link{fit_null_model}} for null model fitting.
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' Dey, R., Schmidt, E. M., Abecasis, G. R., and Lee, S. (2017). A fast and accurate
#' algorithm to test for binary phenotypes and its application to PheWAS. American
#' Journal of Human Genetics, 101(1), 37-49. doi:10.1016/j.ajhg.2017.05.014
#'
#' @export
compute_score_stats <- function(G, null_model, use_spa = NULL, verbose = 1) {
  # ========== Input Validation ==========
  # Validate null_model class
  if (!inherits(null_model, "glow_null_model")) {
    stop("null_model must be a 'glow_null_model' object from fit_null_model()")
  }

  # Validate G is a numeric matrix without NAs
  if (!is.matrix(G) || !is.numeric(G)) {
    stop("G must be a numeric matrix")
  }
  if (any(is.na(G))) {
    stop("G contains NA values. Please impute or remove before analysis.")
  }

  # Validate dimension match (common user error)
  if (nrow(G) != null_model$n) {
    stop("Dimension mismatch: G has ", nrow(G), " rows but null_model has ",
         null_model$n, " samples. Row count of G must equal sample size in null_model.")
  }

  # ========== Determine SPA Usage ==========
  if (is.null(use_spa)) {
    # Auto-detect: SPA for binary, standard for continuous
    use_spa <- (null_model$trait == "binary")
  } else {
    # Validate explicit use_spa
    if (!is.logical(use_spa) || length(use_spa) != 1) {
      stop("use_spa must be NULL, TRUE, or FALSE")
    }
    if (use_spa && null_model$trait != "binary") {
      stop("SPA correction is only available for binary traits. ",
           "The null_model has trait = '", null_model$trait, "'.")
    }
  }

  # ========== Dispatch ==========
  if (use_spa) {
    result <- getZ_marg_score_binary_SPA(G, null_model = null_model)
  } else {
    result <- getZ_marg_score(G, null_model = null_model, trait = null_model$trait)
  }

  # ========== Log dispatch choice ==========
  if (verbose >= 1) {
    message("Score stats: trait=", null_model$trait,
            ", SPA=", use_spa, ", p=", ncol(G))
  }

  return(result)
}
