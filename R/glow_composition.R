# ==============================================================================
# GLOW Test Composition Layer
# ==============================================================================
#
# Data-driven composition functions that build and execute GLOW omnibus tests
# from test specifications. These functions replace hard-coded test assembly
# with a flexible, spec-driven approach while maintaining numerical equivalence
# with the GLOW_Omni pipeline.
#
# INTERNAL FUNCTIONS:
#   - call_weight_fn()          Dispatch a weight function with matched args
#   - run_bsf_tests()           Build and run the BSF omnibus from test specs
#   - generate_test_names()     Produce test names without running computations


#################### INTERNAL HELPER FUNCTIONS ####################

#' Dispatch a weight function with argument matching
#'
#' Calls a weight function with only the arguments it accepts. If the function
#' has a \code{...} formal, all arguments are passed through. Otherwise, only
#' arguments whose names match the function's formals are passed.
#'
#' @param fn A weight-generating function (e.g., \code{Optimal_Weights_M} or
#'   \code{.default_equal_weights}).
#' @param all_args A named list of all available arguments.
#'
#' @return The return value of \code{fn}, typically a named list of weight
#'   vectors.
#'
#' @keywords internal
#' @noRd
call_weight_fn <- function(fn, all_args) {
  fn_formals <- names(formals(fn))
  if ("..." %in% fn_formals) {
    # Function accepts ..., pass everything
    do.call(fn, all_args)
  } else {
    # Only pass arguments that match the function's formals
    do.call(fn, all_args[intersect(names(all_args), fn_formals)])
  }
}


#' Run BSF tests from test specifications
#'
#' Data-driven composition function that builds and executes the GLOW BSF
#' omnibus test from a list of test specifications. For each family in
#' \code{test_specs}, computes weights via the spec's weight functions, stacks
#' them into a combined weight matrix, and calls \code{omni_SgZ_test} to
#' produce statistics and p-values. Optionally appends SNV-level CCT and a
#' final omnibus CCT row.
#'
#' @details
#' The algorithm proceeds in 6 steps:
#' \enumerate{
#'   \item Validate specs and check family name uniqueness.
#'   \item Per-family weight computation: for each spec, iterate over
#'     \code{weight_fns} and collect weight vectors.
#'   \item Stack weights and degrees of freedom across families.
#'   \item Call \code{omni_SgZ_test} on the combined weight/df matrices.
#'   \item If \code{include_snv_cct = TRUE}, compute SNV-level CCT and
#'     final omnibus CCT.
#'   \item Return structured result.
#' }
#'
#' Row names follow the pattern \code{Family_scheme}, where \code{scheme}
#' is derived from the weight function output names with the \code{wts_}
#' prefix stripped and \code{APE} replaced by \code{APR}.
#'
#' **Computational Complexity**: Same as \code{omni_SgZ_test}: O(k * p^2)
#' where k is the total number of weight schemes and p is the number of
#' variants.
#'
#' @param Zscores Numeric vector of marginal Z-scores (length p).
#' @param M Numeric matrix (p x p), correlation matrix of Z-scores.
#' @param Bstar Numeric vector (length p), scaled effect sizes.
#' @param PI Numeric vector (length p), variant-importance scores in \eqn{[0, 1]}.
#' @param test_specs List of test spec lists (see \code{\link{default_test_specs}}).
#' @param include_snv_cct Logical. If TRUE (default), append SNV-level CCT
#'   and final omnibus CCT rows.
#' @param return_weights Logical (default FALSE). If TRUE, attach the
#'   internally-assembled weight matrix \code{WT_combined} (rows = weight
#'   schemes with rownames like \code{Burden_BE}, \code{SKAT_BE_N}, ...;
#'   columns = variants) to the result list as \code{weights}. The matrix
#'   is the exact set of weights consumed by the test, so callers can
#'   persist it for drilldown without risk of recomputation drift.
#' @param ... Additional arguments passed to weight functions (via
#'   \code{call_weight_fn}) and the test routine (\code{omni_SgZ_test}, which
#'   forwards them to \code{p.GFisher}). Weight functions only receive
#'   \code{...} arguments that match their formals (unless they define
#'   \code{...} themselves). To avoid conflicts, custom weight function
#'   parameters should not reuse names used by \code{p.GFisher}
#'   (e.g., \code{method}, \code{nsim}, \code{seed}).
#'
#' @note
#' **Relationship between spec \code{g} and spec \code{df}**:
#'
#' The spec's \code{g} function is used for **weight computation** via
#' \code{Optimal_Weights_M}, where it defines the transformation used in
#' the expected signal-to-noise calculation. The \code{df} value determines
#' the **test statistic transformation**: \code{multi_SgZ_test}
#' calls \code{g_GFisher(x, df=..., p.type=...)} to compute the statistic and
#' routes the p-value calculation accordingly.
#'
#' For the built-in families this is consistent:
#' \itemize{
#'   \item SKAT: \code{g(x) = x^2}, \code{df = 1} --
#'     \code{g_GFisher(x, 1, "two") = x^2}
#'   \item Burden: \code{g(x) = x}, \code{df = Inf} -- uses the normal
#'     distribution (no \code{g_GFisher} call)
#'   \item Fisher: \code{g(x) = g_GFisher_two(x, 2)}, \code{df = 2} --
#'     \code{g_GFisher(x, 2, "two") = g_GFisher_two(x, 2)}
#' }
#'
#' Custom families must ensure \code{g} and \code{df} are consistent:
#' \code{g(x)} should equal \code{g_GFisher(x, df, p.type)} for correct
#' results. Future work may modify the test computation to accept \code{g}
#' directly, removing this constraint.
#'
#' @return A list with:
#' \describe{
#'   \item{STAT}{Numeric matrix of test statistics (one column).}
#'   \item{PVAL}{Numeric matrix of p-values (one column).}
#'   \item{test_names}{Character vector of row names.}
#'   \item{weights}{(Only when \code{return_weights = TRUE}.) Numeric matrix
#'     of weight schemes (n_schemes x p), rows named by scheme (e.g.,
#'     \code{Burden_BE}, \code{SKAT_BE_N}, ...).}
#' }
#'
#' @keywords internal
#' @noRd
run_bsf_tests <- function(Zscores, M, Bstar, PI, test_specs,
                           include_snv_cct = TRUE,
                           return_weights = FALSE, ...) {

  # ========== STEP 1: Validate specs, check uniqueness ==========
  test_specs <- lapply(test_specs, validate_test_spec)
  families <- vapply(test_specs, `[[`, character(1), "family")
  if (anyDuplicated(families)) {
    stop("Duplicate family names in test_specs: ",
         paste(families[duplicated(families)], collapse = ", "))
  }

  # Validate all specs use the same p.type (omni_SgZ_test requires it)
  p_types <- vapply(test_specs, `[[`, character(1), "p.type")
  if (length(unique(p_types)) > 1) {
    stop("All test specs must use the same p.type. Found: ",
         paste(unique(p_types), collapse = ", "))
  }

  p <- length(Zscores)

  # ========== STEP 2: Per-family weight computation ==========
  all_W <- list()
  all_DF <- list()
  all_names <- list()

  for (spec in test_specs) {
    # Build argument pool for weight functions
    # Includes standard Optimal_Weights_M args plus p for .default_equal_weights
    all_args <- list(
      g = spec$g,
      Bstar = Bstar,
      PI = PI,
      M = M,
      p = p,
      is.posi.wts = TRUE
    )

    # Merge ... args (do not override standard args)
    dots <- list(...)
    all_args <- c(all_args, dots[!names(dots) %in% names(all_args)])

    W_rows <- NULL
    scheme_names <- character(0)

    for (wt_fn in spec$weight_fns) {
      wt_result <- call_weight_fn(wt_fn, all_args)
      if (!is.null(wt_result)) {
        for (wt_name in names(wt_result)) {
          # Optimal_Weights_M returns column matrices (p x 1); equal weights
          # returns plain vectors. Coerce to row vector (1 x p) for rbind.
          wt_vec <- as.vector(wt_result[[wt_name]])
          W_rows <- rbind(W_rows, wt_vec)
          # Clean name: strip wts_ prefix, map APE -> APR
          clean_name <- gsub("APE", "APR", sub("^wts_", "", wt_name))
          scheme_names <- c(scheme_names, clean_name)
        }
      }
    }

    # Build df column for this family
    DF_family <- matrix(rep(spec$df, nrow(W_rows)), ncol = 1)

    # Hierarchical row names: Family_scheme
    row_names <- paste0(spec$family, "_", scheme_names)

    all_W[[spec$family]] <- W_rows
    all_DF[[spec$family]] <- DF_family
    all_names[[spec$family]] <- row_names
  }

  # ========== STEP 3: Stack across families ==========
  WT_combined <- do.call(rbind, all_W)
  DF_combined <- do.call(rbind, all_DF)
  combined_names <- unlist(all_names)
  rownames(WT_combined) <- combined_names

  # ========== STEP 4: Call omni_SgZ_test ==========
  omni_result <- omni_SgZ_test(
    Zscores = Zscores,
    DF = DF_combined,
    W = WT_combined,
    M = M,
    p.type = test_specs[[1]]$p.type,
    ...
  )

  # Rename the BSF CCT row (last row from omni_SgZ_test is the CCT)
  n_individual <- length(combined_names)
  rownames(omni_result$STAT)[n_individual + 1] <- "BSF_Omni"
  rownames(omni_result$PVAL)[n_individual + 1] <- "BSF_Omni"

  # ========== STEP 5: SNV-level CCT and final omnibus ==========
  if (include_snv_cct) {
    snv_pvals <- 2 * pnorm(-abs(Zscores))
    snv_cct <- cct_test(snv_pvals)

    # Final omnibus: individual p-values (not the BSF CCT) + SNV CCT
    individual_pvals <- omni_result$PVAL[1:n_individual, 1]
    final_cct <- cct_test(c(individual_pvals, snv_cct$pval_cct))

    result_stat <- rbind(omni_result$STAT, snv_cct$cct, final_cct$cct)
    result_pval <- rbind(omni_result$PVAL, snv_cct$pval_cct, final_cct$pval_cct)
    rownames(result_stat)[(n_individual + 2):(n_individual + 3)] <-
      c("SNV_CCT", "Omni")
    rownames(result_pval)[(n_individual + 2):(n_individual + 3)] <-
      c("SNV_CCT", "Omni")
  } else {
    result_stat <- omni_result$STAT
    result_pval <- omni_result$PVAL
  }

  # ========== STEP 6: Return ==========
  result <- list(
    STAT = result_stat,
    PVAL = result_pval,
    test_names = rownames(result_pval)
  )

  if (isTRUE(return_weights)) {
    # WT_combined has scheme rownames set above; rows = schemes, cols = variants
    result$weights <- WT_combined
  }

  result
}


#' Generate test names from test specifications (without computation)
#'
#' Produces the full vector of test row names that \code{run_bsf_tests} would
#' generate, without actually computing weights or test statistics. This is
#' useful for the single-variant edge case in \code{glow_test} where we need
#' the name vector but skip the variant-set test machinery.
#'
#' @details
#' The function calls each weight function with minimal dummy arguments
#' (p=1, scalar Bstar/PI/M) to extract the output names. For
#' \code{Optimal_Weights_M}, the computation with p=1 is trivial and fast.
#'
#' @param test_specs List of test spec lists.
#' @param include_snv_cct Logical. If TRUE (default), append "BSF_Omni",
#'   "SNV_CCT", and "Omni" names.
#'
#' @return Character vector of test names.
#'
#' @keywords internal
#' @noRd
generate_test_names <- function(test_specs, include_snv_cct = TRUE) {
  test_specs <- lapply(test_specs, validate_test_spec)
  families <- vapply(test_specs, `[[`, character(1), "family")
  if (anyDuplicated(families)) {
    stop("Duplicate family names in test_specs: ",
         paste(families[duplicated(families)], collapse = ", "))
  }

  all_names <- character(0)

  for (spec in test_specs) {
    # Use p=2 dummy arguments for name extraction (p=1 causes issues with
    # diag() of a scalar in Optimal_Weights_M). With p=2, the computation
    # is trivial and fast.
    dummy_args <- list(
      g = spec$g,
      Bstar = c(0.5, 0.5),
      PI = c(0.5, 0.5),
      M = diag(2),
      p = 2,
      is.posi.wts = TRUE
    )

    scheme_names <- character(0)
    for (wt_fn in spec$weight_fns) {
      wt_result <- call_weight_fn(wt_fn, dummy_args)
      if (!is.null(wt_result)) {
        for (wt_name in names(wt_result)) {
          clean_name <- gsub("APE", "APR", sub("^wts_", "", wt_name))
          scheme_names <- c(scheme_names, clean_name)
        }
      }
    }

    # Build hierarchical names: Family_scheme
    family_names <- paste0(spec$family, "_", scheme_names)
    # multi_SgZ_test prepends "df_X_" to row names, so replicate that
    row_names <- paste0("df_", spec$df, "_", family_names)
    all_names <- c(all_names, row_names)
  }

  # Add CCT rows (these are appended by omni_SgZ_test and run_bsf_tests)
  all_names <- c(all_names, "BSF_Omni")
  if (include_snv_cct) {
    all_names <- c(all_names, "SNV_CCT", "Omni")
  }

  all_names
}
