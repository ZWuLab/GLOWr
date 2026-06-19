
#' GLOWr: inteGrative anaLysis using Optimized Weights for genetic Association Analysis - R Package
#'
#' @description
#' The GLOWr package implements the GLOW (inteGrative anaLysis using Optimized Weights)
#' methodology for rare variant association testing in whole-genome sequencing studies.
#' The package provides optimal weighting schemes for multiple test statistics (Burden,
#' SKAT, Fisher, and Omnibus tests) that adaptively combine information from external
#' data sources and functional annotations to maximize statistical power.
#'
#' @section Main Functions:
#'
#' **Core GLOW Tests**:
#' \itemize{
#'   \item \code{\link{GLOW_Burden}}: Burden test with optimal weights
#'   \item \code{\link{GLOW_SKAT}}: SKAT test with optimal weights
#'   \item \code{\link{GLOW_Fisher}}: Fisher combination test with optimal weights
#'   \item \code{\link{GLOW_Omni}}: Omnibus test combining multiple test statistics
#'   \item \code{\link{GLOW_Omni_byP}}: Omnibus test from pre-computed p-values
#' }
#'
#' **Marginal Score Calculation**:
#' \itemize{
#'   \item \code{\link{getZ_marg_score}}: Calculate marginal Z-scores and correlation matrix
#'   \item \code{\link{getZ_marg_score_binary_SPA}}: Z-scores for binary traits with SPA correction
#' }
#'
#' **Parameter Estimation**:
#' \itemize{
#'   \item \code{\link{get_B}}: Estimate effect sizes from external data or pilot studies
#'   \item \code{\link{get_PI}}: Estimate variant-importance scores from functional annotations
#' }
#'
#' **Optimal Weights**:
#' \itemize{
#'   \item \code{\link{Optimal_Weights_M}}: Calculate optimal weights for any test statistic
#' }
#'
#' **Cache Management**:
#' \itemize{
#'   \item \code{\link{clear_glow_cache}}: Clear cached covariance matrices
#'   \item \code{\link{glow_cache_info}}: Check cache status
#' }
#'
#' @section Performance Features:
#'
#' GLOWr includes several performance optimizations:
#' \itemize{
#'   \item \strong{Robust Matrix Inversion}: Automatic correction for ill-conditioned
#'     correlation matrices using nearPD (handles condition numbers up to 10^13)
#'   \item \strong{Covariance Caching}: H0 covariance matrices are cached to speed up
#'     omnibus tests and permutation procedures
#'   \item \strong{Fast Paths}: Optimized computation for independent variants and
#'     simple transformations
#'   \item \strong{Efficient Scaling}: Typical genes (500 variants) analyze in < 1 second;
#'     genome-wide scans feasible in hours
#' }
#'
#' Performance can be controlled via:
#' \itemize{
#'   \item \code{options(GLOWr.use_cache = TRUE/FALSE)} - Enable/disable caching
#' }
#'
#' @section Typical Workflow:
#'
#' 1. **Calculate marginal scores** from genotype and phenotype data:
#'    \code{scores <- getZ_marg_score(G, X, Y)}
#'
#' 2. **Estimate parameters** for optimal weights:
#'    \code{B <- get_B(...)}
#'    \code{PI <- get_PI(...)}
#'
#' 3. **Run GLOW tests**:
#'    \code{burden_result <- GLOW_Burden(scores$Zscores, scores$M_Z, B, PI)}
#'    \code{skat_result <- GLOW_SKAT(scores$Zscores, scores$M_Z, B, PI)}
#'    \code{omni_result <- GLOW_Omni(scores$Zscores, scores$M_Z, B, PI)}
#'
#' @section References:
#'
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @docType package
#' @name GLOWr-package
#' @aliases GLOWr
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats cor pnorm qnorm integrate dnorm
#' @importFrom GFisher p.GFisher
#' @importFrom SPAtest ScoreTest_SPA
#' @importFrom glmnet cv.glmnet
#' @importFrom Matrix nearPD
#' @importFrom data.table data.table
#' @import graphics
#' @import grDevices
## usethis namespace: end
NULL

# Required for data.table's [.data.table dispatch within this package.
# Without this, data.table falls back to [.data.frame semantics when called
# from functions in this namespace, breaking key-based lookups in
# .match_by_position() and other internal functions.
.datatable.aware <- TRUE
