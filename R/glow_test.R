# ==============================================================================
# GLOW Variant-Set Test Runner
# ==============================================================================
#
# Runs GLOW association tests and returns structured S3 result objects.
#
# EXPORTED FUNCTIONS:
#   - glow_test()                    Run GLOW variant-set tests
#   - print.glow_test_result()       Print method for results
#   - as.data.frame.glow_test_result() Flatten to data.frame
#
# INTERNAL HELPER FUNCTIONS:
#   - .validate_glow_test_inputs()   Input validation
#   - .extract_individual_pvalues()  Extract p-values from individual tests
#   - .strip_df_prefix()             Strip df_X_ prefix from row names


#################### EXPORTED MAIN FUNCTIONS ####################

#' Run GLOW Variant-Set Tests
#'
#' Executes GLOW association tests on a set of variants using pre-computed
#' marginal score statistics and per-variant B (effect size) and PI
#' (variant-importance score) parameters. Returns a structured S3 result object that can be
#' printed, converted to a data.frame, or further processed.
#'
#' @details
#' This function is the primary user-facing entry point for running GLOW tests.
#' It wraps the underlying \code{\link{GLOW_Omni}}, \code{\link{GLOW_Burden}},
#' \code{\link{GLOW_SKAT}}, and \code{\link{GLOW_Fisher}} functions, adding
#' input validation, structured output, and region metadata.
#'
#' **Test Modes**:
#' \describe{
#'   \item{\code{tests = "omni"} (default)}{Runs \code{GLOW_Omni}, which
#'     computes all test statistics. With default specs: 5 SKAT, 3 Burden,
#'     5 Fisher, BSF omnibus, SNV-level CCT, and final omnibus CCT (16 total).}
#'   \item{\code{tests = c("burden", "skat", "fisher")}}{Runs individual tests
#'     only. Each test returns its own set of p-values (varying number of weight
#'     schemes). Useful when only specific test types are needed.}
#' }
#'
#' **P-value Naming**:
#' In omnibus mode, p-value names are derived dynamically from the row names
#' produced by \code{run_bsf_tests()}, with a \code{GLOW_} prefix added.
#' With default test specs, this produces 16 names: \code{GLOW_SKAT_BE_N},
#' \code{GLOW_Burden_BE}, \code{GLOW_Fisher_BE_N}, \code{GLOW_BSF_Omni},
#' \code{GLOW_SNV_CCT}, \code{GLOW_Omni}, etc.
#'
#' **Custom Test Specs**:
#' Pass \code{test_specs} to configure which test families are run and how
#' they are weighted. See \code{\link{default_test_specs}} for the format.
#'
#' **Single-Variant Handling**:
#' When only one variant remains (\code{p = 1}), all variant-set tests reduce
#' to the marginal single-SNV test. The marginal p-value is assigned to all
#' output slots, with names derived from \code{generate_test_names()}.
#'
#' @param marg_score_stats List from \code{\link{compute_score_stats}} or
#'   \code{\link{getZ_marg_score}}, containing \code{Zscores}, \code{M_Z},
#'   \code{M_s}, and \code{s0}.
#' @param B Numeric vector (length p). Per-variant effect size parameter.
#' @param PI Numeric vector (length p). Per-variant variant-importance score in
#'   \eqn{[0, 1]}.
#' @param region_info List with optional elements \code{chr}, \code{start},
#'   \code{end}, \code{label}. Region metadata attached to output.
#' @param variant_summary List with variant count information (optional).
#'   Common fields: \code{n_original}, \code{n_after_filter},
#'   \code{n_after_ld}, \code{n_after_collapse}, \code{cMAC}.
#' @param tests Character. Which tests to run: \code{"omni"} (default) runs
#'   the full omnibus pipeline, or a character vector subset of
#'   \code{c("burden", "skat", "fisher")} for individual tests.
#' @param test_specs List of test specification lists (default: NULL uses
#'   \code{\link{default_test_specs}}). Only applies when
#'   \code{tests = "omni"}.
#' @param include_snv_cct Logical. If TRUE (default), include SNV-level CCT
#'   and final omnibus CCT in the output. Only applies when
#'   \code{tests = "omni"}.
#' @param return_weights Logical (default FALSE). If TRUE and
#'   \code{tests = "omni"}, attach the per-test-unit weight matrix used by
#'   the GLOW test to the result as \code{weights} (rows: weight schemes;
#'   columns: variants). In the single-variant edge case, attaches a
#'   \code{n_schemes x 1} matrix of 1s with the schema's scheme rownames so
#'   downstream consumers can rely on a uniform shape.
#' @param verbose Integer. Verbosity level (default 1). Set to 0 for silence.
#' @param ... Additional arguments passed to underlying test functions
#'   (e.g., \code{method}, \code{nsim} for GFisher).
#'
#' @return A \code{glow_test_result} S3 object (list) with components:
#' \describe{
#'   \item{pvalues}{Named numeric vector of p-values.}
#'   \item{statistics}{Named numeric vector of test statistics, same names
#'     as \code{pvalues}.}
#'   \item{region_info}{The \code{region_info} argument, passed through.}
#'   \item{variant_summary}{The \code{variant_summary} argument, passed through.}
#'   \item{settings}{List with \code{tests_run} recording which tests were run.}
#'   \item{raw}{List of raw output from the underlying test functions.}
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @examples
#' \dontrun{
#' # Simulate inputs
#' set.seed(42)
#' p <- 10
#' marg_stats <- list(
#'   Zscores = rnorm(p),
#'   M_Z = diag(p),
#'   M_s = diag(p),
#'   s0 = 1
#' )
#' B <- rep(0.5, p)
#' PI <- runif(p, 0.1, 0.9)
#'
#' # Run omnibus test (default: 16 p-values)
#' result <- glow_test(marg_stats, B, PI,
#'                     region_info = list(label = "GENE1", chr = "22"))
#' print(result)
#' as.data.frame(result)
#'
#' # Custom 2-family spec
#' my_specs <- list(
#'   list(family = "SKAT", g = function(x) x^2, df = 1),
#'   list(family = "Burden", g = function(x) x, df = Inf)
#' )
#' result2 <- glow_test(marg_stats, B, PI, test_specs = my_specs)
#' }
#'
#' @seealso
#' \code{\link{GLOW_Omni}} for the underlying omnibus test
#' \code{\link{default_test_specs}} for the test specification format
#' \code{\link{GLOW_Burden}}, \code{\link{GLOW_SKAT}}, \code{\link{GLOW_Fisher}}
#'   for individual tests
#'
#' @export
glow_test <- function(marg_score_stats, B, PI,
                      region_info = NULL,
                      variant_summary = NULL,
                      tests = "omni",
                      test_specs = NULL,
                      include_snv_cct = TRUE,
                      return_weights = FALSE,
                      verbose = 1,
                      ...) {

  # ------ Step 1: Validate inputs ------
  .validate_glow_test_inputs(marg_score_stats, B, PI, tests)

  p <- length(marg_score_stats$Zscores)

  # Resolve test specs (used in omni and single-variant paths)
  specs <- if (is.null(test_specs)) default_test_specs() else test_specs

  # ------ Handle single-variant edge case ------
  if (p == 1) {
    Z <- marg_score_stats$Zscores[1]
    marginal_pval <- 2 * pnorm(-abs(Z))

    # Derive names dynamically from test specs
    all_names <- generate_test_names(specs, include_snv_cct = include_snv_cct)
    glow_names <- paste0("GLOW_", .strip_df_prefix(all_names))

    # All p-values are the marginal p-value
    pvalues <- stats::setNames(rep(marginal_pval, length(glow_names)),
                               glow_names)

    # Statistics: use the spec's g(Z) for each family, CCT stat for
    # omnibus rows. Build a lookup from family name to g(Z) value.
    cct_stat <- tan((0.5 - marginal_pval) * pi)
    family_gZ <- vapply(specs, function(s) {
      s$g(Z)
    }, numeric(1))
    names(family_gZ) <- vapply(specs, `[[`, character(1), "family")

    statistics <- vapply(all_names, function(nm) {
      if (grepl("^df_", nm)) {
        # Extract family name: "df_X_Family_scheme" -> "Family"
        # Strip "df_<number>_" or "df_Inf_" prefix, then take first token
        stripped <- sub("^df_(\\d+|Inf)_", "", nm)
        fam <- sub("_.*", "", stripped)
        if (fam %in% names(family_gZ)) {
          family_gZ[[fam]]
        } else {
          Z^2  # fallback (should not happen with valid specs)
        }
      } else {
        cct_stat  # CCT rows (BSF_Omni, SNV_CCT, Omni)
      }
    }, numeric(1))
    statistics <- stats::setNames(as.vector(statistics), glow_names)

    if (verbose >= 1) {
      message("GLOW test: p=1 (single variant), marginal p-value=",
              format(marginal_pval, digits = 3))
    }

    sv_result <- list(
      pvalues = pvalues,
      statistics = statistics,
      region_info = region_info,
      variant_summary = variant_summary,
      settings = list(tests_run = tests, single_variant = TRUE),
      raw = list(single_variant = TRUE, Z = Z, marginal_pval = marginal_pval)
    )

    if (isTRUE(return_weights) && "omni" %in% tests) {
      # Build a uniform-shape weight matrix (n_schemes x 1) of 1s with rownames
      # matching the per-family schemes that run_bsf_tests would emit. We
      # derive the rownames from generate_test_names() by dropping df_X_
      # prefix and the trailing CCT rows.
      raw_names <- generate_test_names(specs, include_snv_cct = FALSE)
      # generate_test_names appends "BSF_Omni" at the end; drop it.
      raw_names <- raw_names[raw_names != "BSF_Omni"]
      scheme_names <- .strip_df_prefix(raw_names)
      w_mat <- matrix(1, nrow = length(scheme_names), ncol = 1L,
                      dimnames = list(scheme_names, NULL))
      sv_result$weights <- w_mat
    }

    return(structure(sv_result, class = "glow_test_result"))
  }

  # ------ Step 2: Run tests ------
  raw <- list()
  pvalues <- numeric(0)
  statistics <- numeric(0)

  if ("omni" %in% tests) {
    # Full omnibus test via GLOW_Omni (delegates to run_bsf_tests)
    raw$omni <- GLOW_Omni(marg_score_stats, B, PI,
                           test_specs = specs,
                           include_snv_cct = include_snv_cct,
                           return_weights = return_weights,
                           ...)

    # Name-based extraction: strip df_X_ prefix, add GLOW_ prefix
    raw_names <- rownames(raw$omni$PVAL)
    glow_names <- paste0("GLOW_", .strip_df_prefix(raw_names))

    pvalues <- stats::setNames(as.vector(raw$omni$PVAL[, 1]), glow_names)
    statistics <- stats::setNames(as.vector(raw$omni$STAT[, 1]), glow_names)

  } else {
    # Individual tests: run each requested test and concatenate results

    if ("burden" %in% tests) {
      raw$burden <- GLOW_Burden(marg_score_stats, B, PI, ...)
      burden_pv <- .extract_individual_pvalues(raw$burden)
      pvalues <- c(pvalues, burden_pv$pvalues)
      statistics <- c(statistics, burden_pv$statistics)
    }

    if ("skat" %in% tests) {
      raw$skat <- GLOW_SKAT(marg_score_stats, B, PI, ...)
      skat_pv <- .extract_individual_pvalues(raw$skat)
      pvalues <- c(pvalues, skat_pv$pvalues)
      statistics <- c(statistics, skat_pv$statistics)
    }

    if ("fisher" %in% tests) {
      raw$fisher <- GLOW_Fisher(marg_score_stats, B, PI, ...)
      fisher_pv <- .extract_individual_pvalues(raw$fisher)
      pvalues <- c(pvalues, fisher_pv$pvalues)
      statistics <- c(statistics, fisher_pv$statistics)
    }
  }

  # ------ Step 3: Log summary ------
  if (verbose >= 1) {
    message(
      "GLOW test: p=", p,
      ", tests=", paste(tests, collapse = ","),
      ", min_pval=", format(min(pvalues), digits = 3)
    )
  }

  # ------ Step 4: Construct S3 result object ------
  result <- list(
    pvalues = pvalues,
    statistics = statistics,
    region_info = region_info,
    variant_summary = variant_summary,
    settings = list(tests_run = tests),
    raw = raw
  )

  if (!is.null(raw$omni$weights)) {
    result$weights <- raw$omni$weights
  }

  structure(result, class = "glow_test_result")
}


#' Print Method for GLOW Test Results
#'
#' Prints a concise summary of a \code{glow_test_result} object, showing
#' region label, variant counts, tests run, and top p-values.
#'
#' @param x A \code{glow_test_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.glow_test_result <- function(x, ...) {
  cat("GLOW Test Result")
  if (!is.null(x$region_info$label)) {
    cat(": ", x$region_info$label, sep = "")
  }
  cat("\n")

  # Variant summary
  if (!is.null(x$variant_summary)) {
    if (!is.null(x$variant_summary$n_after_collapse)) {
      cat("  Variants: ", x$variant_summary$n_after_collapse,
          " (from ", x$variant_summary$n_original, " original)\n", sep = "")
    }
    if (!is.null(x$variant_summary$cMAC)) {
      cat("  Cumulative MAC: ", x$variant_summary$cMAC, "\n", sep = "")
    }
  }

  # Tests run
  cat("  Tests: ", paste(x$settings$tests_run, collapse = ", "), "\n", sep = "")

  # Top 5 p-values
  n_show <- min(5, length(x$pvalues))
  top <- sort(x$pvalues)[seq_len(n_show)]
  cat("  Top p-values:\n")
  for (nm in names(top)) {
    cat("    ", nm, ": ", format(top[nm], digits = 3), "\n", sep = "")
  }

  invisible(x)
}


#' Convert GLOW Test Results to Data Frame
#'
#' Flattens a \code{glow_test_result} object into a single-row data.frame
#' suitable for binding across regions (e.g., \code{rbind} or
#' \code{do.call(rbind, list_of_dfs)}).
#'
#' @param x A \code{glow_test_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A single-row \code{data.frame} with columns for region info,
#'   variant summary, and all p-values (dynamic number of columns).
#'
#' @export
as.data.frame.glow_test_result <- function(x, ...) {
  # NULL-safe accessor (do NOT use %||%)
  null_or <- function(a, b) if (is.null(a)) b else a

  # Region info columns
  ri <- x$region_info
  region_df <- data.frame(
    label = null_or(ri$label, NA_character_),
    chr   = null_or(ri$chr, NA_character_),
    start = null_or(ri$start, NA_integer_),
    end   = null_or(ri$end, NA_integer_),
    stringsAsFactors = FALSE
  )

  # Variant summary columns
  vs <- x$variant_summary
  var_df <- data.frame(
    n_variants         = null_or(vs$n_original, NA_integer_),
    n_after_annotation = null_or(vs$n_after_annotation, NA_integer_),
    n_after_filter     = null_or(vs$n_after_filter, NA_integer_),
    n_after_ld         = null_or(vs$n_after_ld, NA_integer_),
    n_after_collapse   = null_or(vs$n_after_collapse, NA_integer_),
    cMAC               = null_or(vs$cMAC, NA_real_),
    stringsAsFactors = FALSE
  )

  # P-value columns (each p-value becomes its own column, handles any count)
  pval_df <- as.data.frame(as.list(x$pvalues))

  cbind(region_df, var_df, pval_df)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Strip df_X_ prefix from test row names
#'
#' Removes the \code{df_<number>_} or \code{df_Inf_} prefix that
#' \code{multi_SgZ_test()} prepends to row names. Names without this
#' prefix (e.g., "BSF_Omni", "SNV_CCT", "Omni") are returned unchanged.
#'
#' @param names Character vector of row names.
#'
#' @return Character vector with prefixes stripped.
#'
#' @keywords internal
#' @noRd
.strip_df_prefix <- function(names) {
  sub("^df_(\\d+|Inf)_", "", names)
}


#' Validate inputs for glow_test
#'
#' Checks that marg_score_stats, B, PI, and tests arguments are valid.
#'
#' @param marg_score_stats List with Zscores, M_Z, M_s, s0.
#' @param B Numeric vector of effect sizes.
#' @param PI Numeric vector of variant-importance scores.
#' @param tests Character vector of test names.
#'
#' @return NULL (invisible). Throws informative errors on invalid input.
#'
#' @keywords internal
#' @noRd
.validate_glow_test_inputs <- function(marg_score_stats, B, PI, tests) {
  # Check marg_score_stats structure
  required_fields <- c("Zscores", "M_Z", "M_s", "s0")
  missing_fields <- setdiff(required_fields, names(marg_score_stats))
  if (length(missing_fields) > 0) {
    stop("marg_score_stats is missing required fields: ",
         paste(missing_fields, collapse = ", "))
  }

  p <- length(marg_score_stats$Zscores)

  # Check B length matches
  if (length(B) != p) {
    stop("B must have length ", p, " (matching Zscores), but has length ",
         length(B))
  }

  # Check PI length matches
  if (length(PI) != p) {
    stop("PI must have length ", p, " (matching Zscores), but has length ",
         length(PI))
  }

  # Check for NAs
  if (any(is.na(B))) {
    stop("B must not contain NA values")
  }
  if (any(is.na(PI))) {
    stop("PI must not contain NA values")
  }
  if (any(is.na(marg_score_stats$Zscores))) {
    stop("marg_score_stats$Zscores must not contain NA values")
  }

  # Check PI range
  if (any(PI < 0) || any(PI > 1)) {
    stop("PI values must be in [0, 1]")
  }

  # Check B and PI are numeric
  if (!is.numeric(B)) {
    stop("B must be a numeric vector")
  }
  if (!is.numeric(PI)) {
    stop("PI must be a numeric vector")
  }

  # Check M_Z is p x p
  if (!is.matrix(marg_score_stats$M_Z) ||
      nrow(marg_score_stats$M_Z) != p ||
      ncol(marg_score_stats$M_Z) != p) {
    stop("marg_score_stats$M_Z must be a ", p, " x ", p, " matrix")
  }

  # Validate tests argument
  valid_tests <- c("omni", "burden", "skat", "fisher")
  invalid_tests <- setdiff(tests, valid_tests)
  if (length(invalid_tests) > 0) {
    stop("Invalid test(s): ", paste(invalid_tests, collapse = ", "),
         ". Must be one or more of: ", paste(valid_tests, collapse = ", "))
  }

  invisible(NULL)
}


#' Extract p-values and statistics from individual test output
#'
#' Converts the matrix output from GLOW_Burden/SKAT/Fisher into named
#' vectors using the row names from the matrices.
#'
#' @param test_result List with STAT and PVAL matrices (from individual tests).
#'
#' @return List with named vectors \code{pvalues} and \code{statistics}.
#'
#' @keywords internal
#' @noRd
.extract_individual_pvalues <- function(test_result) {
  # Use row names from the matrices (these are unique within each test)
  rnames <- rownames(test_result$PVAL)

  pvalues <- stats::setNames(
    as.vector(test_result$PVAL[, 1]),
    rnames
  )
  statistics <- stats::setNames(
    as.vector(test_result$STAT[, 1]),
    rnames
  )

  list(pvalues = pvalues, statistics = statistics)
}
