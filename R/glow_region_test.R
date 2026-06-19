# ==============================================================================
# Single-Region GLOW Variant-Set Test
# ==============================================================================
#
# One region in, one structured result out. A gene-centric scan and a
# sliding-window scan share this single per-region code path, so both analyses
# test each region identically:
#   extract_variant_set -> prepare_glow_input (B-source dispatch)
#   -> compute_score_stats (+ standard-Z for binary) -> glow_test
#   -> optional STAAR-O on the post-filter genotype matrix
#   -> optional evidence packet.
#
# EXPORTED FUNCTIONS:
#   - glow_region_test()             Single-region GLOW (+ optional STAAR) test
#
# INTERNAL HELPER FUNCTIONS:
#   - .build_region_evidence_packet()  Per-region evidence packet
#   - .extract_staar_family_pvalues()  Extract the 8 STAAR family-omnibus p-values
#   - .run_region_staar()              Median-impute + STAAR call + skip/err logic
#   - .staar_family_cols()             The 8 STAAR family-omnibus column names
#   - .null_or()                       NULL-coalescing helper (local)


#################### EXPORTED MAIN FUNCTION ####################

#' Run the GLOW Test on a Single Genomic Region
#'
#' Tests one genomic region (a gene or a sliding window) with the GLOW
#' variant-set procedure, optionally adding a STAAR-O test and an evidence
#' packet. Both a gene-centric scan and a sliding-window scan call this function
#' once per region, so the two analyses run identical per-region code.
#'
#' @details
#' The processing pipeline, for one region, is:
#' \enumerate{
#'   \item \code{\link{extract_variant_set}} on the open \code{gds} (skip if the
#'     region has no qualifying variants).
#'   \item \code{\link{prepare_glow_input}} with B-source dispatch (exactly one
#'     of \code{b_func} / \code{b_model}); LD-prune + rare-variant collapse
#'     (skip if nothing survives).
#'   \item \code{\link{compute_score_stats}}, honoring \code{use_spa} (default
#'     NULL = auto: SPA for binary, standard for continuous). For a binary trait
#'     the complementary single-variant Z is also computed for the evidence
#'     packet (the SPA complement only when \code{return_evidence = TRUE}).
#'   \item \code{\link{glow_test}} (omnibus), with weights returned when
#'     \code{return_evidence = TRUE}.
#'   \item If \code{staar} is non-NULL, \code{STAAR::STAAR} on the
#'     post-filter / pre-LD-prune / pre-collapse genotype matrix
#'     (\code{vset$G}), with missing annotations median-imputed against
#'     \code{reference_medians} (STAAR does its own LD handling and collapsing).
#'   \item If \code{return_evidence = TRUE}, a per-region evidence packet.
#' }
#'
#' \strong{B-source contract.} Exactly one of \code{b_func}, \code{b_model} must
#' be non-NULL (mirrors \code{\link{prepare_glow_input}}). Supplying both or
#' neither is a caller error and throws.
#'
#' \strong{Status semantics.} The caller reads \code{$status} rather than
#' string-matching a message:
#' \describe{
#'   \item{\code{"ok"}}{GLOW produced a result row (\code{$result} non-NULL).
#'     STAAR, if requested, either ran (p-values in \code{$staar}) or failed for
#'     a non-skip reason (\code{$staar} all-NA, message in \code{$message});
#'     either way GLOW succeeded.}
#'   \item{\code{"skip_empty"}}{\code{extract_variant_set} or
#'     \code{prepare_glow_input} returned NULL (no qualifying variants).
#'     \code{$result} is NULL.}
#'   \item{\code{"staar_skip_single"}}{GLOW produced a row but STAAR could not run
#'     (fewer than 2 rare variants — STAAR's design-intended exclusion).
#'     \code{$result} is the GLOW row; \code{$staar} is all-NA.}
#'   \item{\code{"error"}}{An unexpected error in the GLOW path. \code{$result}
#'     is NULL and \code{$message} holds the error message.}
#' }
#'
#' @param gds An \strong{open} \code{SeqVarGDSClass} connection. The orchestration
#'   layer opens the GDS once per chunk and reuses it across windows, so this
#'   function does not open or close it. (A file path also works because
#'   \code{\link{extract_variant_set}} accepts one, but the intended contract is
#'   an open connection to avoid re-opening the GDS in the hot loop.)
#' @param region A list or single-row data.frame with \code{chr}, \code{start},
#'   \code{end}, and (optionally) \code{label} — the row shape produced by
#'   \code{\link{define_regions_window}} / \code{\link{define_regions_gene}}.
#' @param filter_spec A \code{glow_variant_filter} from
#'   \code{\link{variant_filter}} or \code{\link{coding_filter}}.
#' @param null_model A fitted \code{glow_null_model} from
#'   \code{\link{fit_null_model}}. Its \code{sample_id} drives sample
#'   extraction/ordering and its \code{trait} drives the score-stat dispatch.
#' @param pi_models A loaded PI ensemble (from \code{\link{load_PI_models}}) or a
#'   model-directory path, passed to \code{prepare_glow_input} as \code{PI_models}.
#' @param pi_features Character vector of annotation column names for PI
#'   prediction. Also the annotations extracted from the GDS, so any
#'   \code{staar$anno_cols} must be a subset of these.
#' @param reference_medians Named numeric vector of chromosome-wide annotation
#'   medians (from \code{\link{compute_annotation_medians}}) for NA imputation in
#'   both PI prediction and STAAR annotation imputation.
#' @param b_func Function mapping a numeric MAF vector to B values. Mutually
#'   exclusive with \code{b_model}.
#' @param b_model A \code{glow_B_model} (from \code{\link{train_B_model}} /
#'   \code{\link{load_B_model}}). Mutually exclusive with \code{b_func}.
#' @param ld_threshold Numeric (default 0.9). LD-pruning correlation threshold.
#' @param mac_threshold Integer (default 10). MAC threshold for rare-variant
#'   collapsing.
#' @param collapse_method Character (default "mean"). B/PI aggregation method
#'   for collapsed groups.
#' @param use_spa Logical or NULL (default NULL). Selects GLOW's single-variant
#'   score statistics: NULL auto-detects (SPA for binary, standard for
#'   continuous); TRUE forces SPA (binary only — a non-binary trait errors);
#'   FALSE forces standard for any trait. GLOW-only: every GLOW family test
#'   combines the chosen Z, while STAAR (if requested) runs its own, always
#'   non-SPA, path. SPA is recommended for rare variants / unbalanced
#'   case-control; see the guidance in \code{\link{getZ_marg_score}}.
#' @param z_scale Numeric scalar (default 1). Single-variant confounding
#'   correction (method M2) applied to every unit's Z \emph{before} the GLOW
#'   combination: \code{z_scale = 1/sqrt(inflation_factor)}, where
#'   \code{inflation_factor} is the single-variant genomic inflation (e.g. the
#'   LD Score regression intercept from \code{\link{ldsc_regression}}). It
#'   multiplies the score-statistic Z-scores that \code{\link{glow_test}}
#'   combines (and the evidence-packet \code{Z_standard}/\code{Z_SPA}); the
#'   correlation \code{M_Z} is invariant to a uniform Z rescaling and
#'   \code{M_s}/\code{s0} are left untouched. \code{z_scale = 1} is a no-op.
#' @param staar STAAR context, or NULL (default) for GLOW only. When non-NULL, a
#'   list with \code{null_model} (a \code{STAAR::fit_null_glm} object),
#'   \code{anno_cols} (character vector of annotation column names, a subset of
#'   \code{pi_features}), and \code{rare_maf_cutoff} (numeric).
#' @param return_evidence Logical (default FALSE). If TRUE, build and return the
#'   per-region evidence packet (post-filter variant table + per-unit
#'   B/PI/Z/weights), and request weights from \code{glow_test}.
#' @param region_summary Logical (default FALSE). If TRUE, splice 7 region-level
#'   diagnostic columns into \code{$result} (after \code{cMAC}, before the
#'   \code{GLOW_*} p-values): \code{n_collapsed_units} (number of collapsed test
#'   units), \code{min_MAF}/\code{max_MAF} and \code{min_MAC}/\code{max_MAC}
#'   (MAF/MAC range over the post-collapse test units), and
#'   \code{min_p_standard}/\code{min_p_SPA} (smallest per-unit single-variant
#'   p-value on the standard / SPA channel; \code{min_p_SPA} is \code{NA} unless
#'   the SPA Z exists, i.e. a binary trait). These are standard per-region
#'   diagnostic columns. They are cheap (computed from \code{input$G} and
#'   the score-stat Z this function already holds) and do not require
#'   \code{return_evidence}. Default FALSE leaves the canonical row unchanged.
#' @param verbose Integer (default 0). Verbosity threaded to the underlying calls.
#'
#' @return A list with components:
#' \describe{
#'   \item{result}{A single-row \code{data.frame} (the
#'     \code{as.data.frame.glow_test_result} shape: region columns, variant
#'     counts, and the GLOW p-value columns; plus the 7 region-summary columns
#'     spliced after \code{cMAC} when \code{region_summary = TRUE}), or NULL when
#'     \code{status \%in\% c("skip_empty", "error")}.}
#'   \item{evidence}{The evidence packet, or NULL (non-NULL only when
#'     \code{return_evidence = TRUE} and \code{status} kept a row).}
#'   \item{staar}{A named list of the 8 STAAR family-omnibus p-values when a
#'     \code{staar} context was supplied (all-NA on \code{staar_skip_single} or a
#'     STAAR error), or NULL when STAAR was off.}
#'   \item{staar_detail}{The full \code{STAAR::STAAR} output, or NULL.}
#'   \item{status}{One of \code{"ok"}, \code{"skip_empty"}, \code{"staar_skip_single"},
#'     \code{"error"} (see Details).}
#'   \item{message}{A diagnostic message: the error message when
#'     \code{status == "error"}, a STAAR error note when STAAR failed for a
#'     non-skip reason (GLOW still succeeded), else NULL.}
#' }
#'
#' @seealso \code{\link{extract_variant_set}}, \code{\link{prepare_glow_input}},
#'   \code{\link{compute_score_stats}}, \code{\link{glow_test}},
#'   \code{\link{define_regions_window}}
#'
#' @examples
#' \dontrun{
#' library(SeqArray)
#' gds <- seqOpen("chr22_essentialdb.gds")
#' on.exit(seqClose(gds))
#'
#' region <- list(chr = "22", start = 16570000L, end = 16670000L,
#'                label = "chr22:16570000-16670000")
#' spec <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV")
#' pi   <- load_PI_models("path/to/piModels/")
#' feats <- c("cadd_phred", "linsight", "fathmm_xf")
#' meds <- compute_annotation_medians(gds, feats, spec,
#'                                    sample_id = null_model$sample_id)
#'
#' res <- glow_region_test(
#'   gds, region, spec, null_model,
#'   pi_models = pi, pi_features = feats, reference_medians = meds,
#'   b_func = function(maf) sqrt(-0.14307 * log(maf * (1 - maf)) - 0.19522),
#'   ld_threshold = 0.95, mac_threshold = 11L)
#' res$status
#' res$result
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' Li, X., Li, Z., Zhou, H., et al. (2020). Dynamic incorporation of multiple in
#' silico functional annotations empowers rare variant association analysis of large
#' whole-genome sequencing studies at scale. Nature Genetics, 52, 969-983.
#' doi:10.1038/s41588-020-0676-4
#'
#' Li, Z., Li, X., Zhou, H., et al. (2022). A framework for detecting noncoding
#' rare-variant associations of large-scale whole-genome sequencing studies.
#' Nature Methods, 19, 1599-1611. doi:10.1038/s41592-022-01640-x
#'
#' @export
glow_region_test <- function(gds, region, filter_spec,
                             null_model,
                             pi_models, pi_features, reference_medians,
                             b_func = NULL, b_model = NULL,
                             ld_threshold = 0.9, mac_threshold = 10L,
                             collapse_method = "mean",
                             use_spa = NULL,
                             z_scale = 1,
                             staar = NULL,
                             return_evidence = FALSE,
                             region_summary = FALSE,
                             verbose = 0) {

  # ---- Caller-contract validation (fail fast on misuse) ----
  n_b_sources <- sum(!is.null(b_func), !is.null(b_model))
  if (n_b_sources != 1L) {
    stop("Exactly one of `b_func`, `b_model` must be provided (got ",
         n_b_sources, ").")
  }
  if (!is.numeric(z_scale) || length(z_scale) != 1L || !is.finite(z_scale) ||
      z_scale <= 0) {
    stop("`z_scale` must be a single finite positive number (got ",
         deparse(z_scale), ").", call. = FALSE)
  }
  if (!is.null(use_spa) &&
      (!is.logical(use_spa) || length(use_spa) != 1L || is.na(use_spa))) {
    stop("`use_spa` must be NULL, TRUE, or FALSE (got ", deparse(use_spa), ").",
         call. = FALSE)
  }
  if (isTRUE(use_spa) && null_model$trait != "binary") {
    stop("`use_spa = TRUE` requires a binary trait (null_model$trait = '",
         null_model$trait, "').", call. = FALSE)
  }
  if (!is.null(staar)) {
    if (!requireNamespace("STAAR", quietly = TRUE)) {
      stop("A `staar` context was supplied but the STAAR package is not ",
           "installed. Install STAAR (e.g. remotes::install_github(",
           "'xihaoli/STAAR')) or pass staar = NULL for GLOW only.")
    }
    .validate_staar_context(staar)
  }

  # NULL-/error-return template keeps the contract uniform.
  empty_result <- function(status, message = NULL) {
    list(result = NULL, evidence = NULL, staar = NULL, staar_detail = NULL,
         status = status, message = message)
  }

  # The GLOW path (steps 1-4, 6) is wrapped in tryCatch so an unexpected
  # failure becomes status = "error" rather than an uncaught throw (mirrors
  # 02-run-glow-chr.R:416-419). NULL returns from extract/prepare are normal
  # "skip_empty" outcomes, not errors, and return directly.
  tryCatch({

    # 1. Extract variants (02-run-glow-chr.R:300-304)
    vset <- extract_variant_set(gds, region, filter_spec,
                                sample_id = null_model$sample_id,
                                annotation_names = pi_features,
                                verbose = verbose)
    if (is.null(vset)) return(empty_result("skip_empty"))

    # 2. Prepare GLOW input — B-source dispatch (02-run-glow-chr.R:307-322)
    prep_args <- list(
      variant_set       = vset,
      PI_models         = pi_models,
      PI_features       = pi_features,
      reference_medians = reference_medians,
      ld_threshold      = ld_threshold,
      mac_threshold     = mac_threshold,
      collapse_method   = collapse_method,
      verbose           = verbose
    )
    if (!is.null(b_func)) prep_args$B_func <- b_func else prep_args$B_model <- b_model
    input <- do.call(prepare_glow_input, prep_args)
    if (is.null(input)) return(empty_result("skip_empty"))

    # 3. Score statistics. The dispatched stat (what glow_test() combines) honors
    #    `use_spa` (NULL = auto: SPA for binary, standard for continuous).
    #    `used_spa` is the resolved choice, reused for the evidence flag (step 6).
    #    (02-run-glow-chr.R:326-338)
    used_spa <- if (is.null(use_spa)) (null_model$trait == "binary") else use_spa
    stats_dispatch <- compute_score_stats(input$G, null_model,
                                          use_spa = use_spa, verbose = verbose)

    # The evidence packet wants both single-variant p-values for a binary trait;
    # one of the pair is the dispatch, so compute only the complement — and the
    # expensive SPA complement only when the evidence packet needs it, so a
    # non-SPA scan (return_evidence = FALSE) pays no SPA cost.
    if (null_model$trait == "binary") {
      if (used_spa) {                       # default binary path — unchanged
        Z_SPA      <- stats_dispatch$Zscores
        Z_standard <- compute_score_stats(input$G, null_model, use_spa = FALSE,
                                          verbose = verbose)$Zscores
      } else {                              # use_spa = FALSE override
        Z_standard <- stats_dispatch$Zscores
        Z_SPA      <- if (isTRUE(return_evidence))
          compute_score_stats(input$G, null_model, use_spa = TRUE,
                              verbose = verbose)$Zscores else NULL
      }
    } else {
      Z_standard <- stats_dispatch$Zscores
      Z_SPA <- NULL
    }

    # unit's Z by `z_scale` (= 1/sqrt(single-variant inflation factor)) BEFORE the
    # GLOW combination. glow_test() reads `stats_dispatch$Zscores`; the correlation
    # M_Z is invariant to a uniform Z rescaling and M_s/s0 are left untouched.
    # z_scale = 1 is a no-op; Z_standard/Z_SPA are scaled too so the evidence
    # packet reflects the corrected inputs.
    if (z_scale != 1) {
      stats_dispatch$Zscores <- stats_dispatch$Zscores * z_scale
      Z_standard <- Z_standard * z_scale
      if (!is.null(Z_SPA)) Z_SPA <- Z_SPA * z_scale
    }

    # 4. GLOW test (02-run-glow-chr.R:341-351)
    glow_result <- glow_test(
      stats_dispatch, input$B, input$PI,
      region_info = region,
      variant_summary = list(
        n_original         = input$n_original,
        n_after_annotation = vset$n_after_annotation,
        n_after_filter     = input$n_after_filter,
        n_after_ld         = input$n_after_ld,
        n_after_collapse   = input$n_after_collapse,
        cMAC               = input$cMAC
      ),
      return_weights = return_evidence,
      verbose = verbose
    )
    result_row <- as.data.frame(glow_result)

    # 4b. Optional region-summary columns. Splice the 7 diagnostic columns in
    #     after cMAC: the post-collapse MAF/MAC ranges (input$G), the
    #     collapsed-unit count, and the per-unit single-variant min p-values
    #     2*pnorm(-abs(Z)). The SPA column needs Z_SPA, which is non-NULL only on
    #     the binary path (default SPA dispatch); a continuous trait gets NA.
    #     Cheap and independent of return_evidence.
    if (isTRUE(region_summary)) {
      result_row <- .splice_region_summary(result_row, input,
                                           Z_standard = Z_standard,
                                           Z_SPA = Z_SPA)
    }

    # 5. STAAR-O on the post-filter vset$G (02-run-glow-chr.R:364-402)
    status        <- "ok"
    staar_pvalues <- NULL
    staar_detail  <- NULL
    message       <- NULL
    if (!is.null(staar)) {
      staar_pvalues <- .extract_staar_family_pvalues(NULL)   # all-NA default
      so <- .run_region_staar(vset, staar, reference_medians)
      if (isTRUE(so$skipped)) {
        # < 2 rare variants: STAAR's design-intended exclusion, not an error.
        status <- "staar_skip_single"
      } else if (!is.null(so$error)) {
        # A non-skip STAAR failure: GLOW already succeeded, so keep the row,
        # leave STAAR NA, and surface the message (does not fail the region).
        message <- so$error
      } else {
        staar_pvalues <- .extract_staar_family_pvalues(so$staar_out)
        staar_detail  <- so$staar_out
      }
    }

    # 6. Evidence packet (02-run-glow-chr.R:405-409)
    evidence <- NULL
    if (isTRUE(return_evidence)) {
      evidence <- .build_region_evidence_packet(
        vset = vset, input = input,
        Z_standard = Z_standard, Z_SPA = Z_SPA,
        result = glow_result, null_model = null_model, used_spa = used_spa
      )
    }

    # 7. Assemble return list
    list(result = result_row, evidence = evidence,
         staar = staar_pvalues, staar_detail = staar_detail,
         status = status, message = message)

  }, error = function(e) {
    empty_result("error", conditionMessage(e))
  })
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' NULL-coalescing helper
#'
#' Returns \code{a} unless it is NULL, in which case \code{b}. Defined here so
#' this file is self-contained.
#'
#' @keywords internal
#' @noRd
.null_or <- function(a, b) if (is.null(a)) b else a


#' Validate a STAAR context list
#'
#' @param staar A STAAR context (see \code{glow_region_test}'s \code{staar}).
#' @keywords internal
#' @noRd
.validate_staar_context <- function(staar) {
  if (!is.list(staar)) stop("`staar` must be a list (the STAAR context) or NULL.")
  missing_fields <- setdiff(c("null_model", "anno_cols", "rare_maf_cutoff"),
                            names(staar))
  if (length(missing_fields) > 0) {
    stop("STAAR context is missing required field(s): ",
         paste(missing_fields, collapse = ", "),
         ". Expected: null_model, anno_cols, rare_maf_cutoff.")
  }
  invisible(NULL)
}


#' The 8 STAAR family-omnibus column names
#'
#' The fixed names (and order) of the STAAR family-omnibus p-values collected at
#' the region level.
#'
#' @keywords internal
#' @noRd
.staar_family_cols <- function() {
  c("STAAR_O", "ACAT_O",
    "STAAR_S_1_25", "STAAR_S_1_1",
    "STAAR_B_1_25", "STAAR_B_1_1",
    "STAAR_A_1_25", "STAAR_A_1_1")
}


#' Project the 8 STAAR family-omnibus p-values from a STAAR() output
#'
#' Returns an all-NA named list when \code{staar_out} is NULL.
#'
#' @param staar_out The full output of \code{STAAR::STAAR}, or NULL.
#' @return A named list of 8 numeric p-values (the \code{.staar_family_cols()}).
#' @keywords internal
#' @noRd
.extract_staar_family_pvalues <- function(staar_out) {
  cols <- .staar_family_cols()
  na_row <- stats::setNames(as.list(rep(NA_real_, length(cols))), cols)
  if (is.null(staar_out)) return(na_row)
  get_omnibus <- function(df, col_name) {
    if (is.null(df) || !is.data.frame(df) || !(col_name %in% names(df))) {
      return(NA_real_)
    }
    as.numeric(df[[col_name]][1L])
  }
  list(
    STAAR_O      = .null_or(staar_out$results_STAAR_O, NA_real_),
    ACAT_O       = .null_or(staar_out$results_ACAT_O,  NA_real_),
    STAAR_S_1_25 = get_omnibus(staar_out$results_STAAR_S_1_25, "STAAR-S(1,25)"),
    STAAR_S_1_1  = get_omnibus(staar_out$results_STAAR_S_1_1,  "STAAR-S(1,1)"),
    STAAR_B_1_25 = get_omnibus(staar_out$results_STAAR_B_1_25, "STAAR-B(1,25)"),
    STAAR_B_1_1  = get_omnibus(staar_out$results_STAAR_B_1_1,  "STAAR-B(1,1)"),
    STAAR_A_1_25 = get_omnibus(staar_out$results_STAAR_A_1_25, "STAAR-A(1,25)"),
    STAAR_A_1_1  = get_omnibus(staar_out$results_STAAR_A_1_1,  "STAAR-A(1,1)")
  )
}


#' Run STAAR-O on a region's post-filter genotype matrix
#'
#' Median-imputes missing annotations (so STAAR's kernel eigendecomposition does
#' not fail on unannotated variants), then calls \code{STAAR::STAAR} inside a
#' tryCatch that distinguishes STAAR's design-intended "< 2 rare variants" skip
#' (and a NULL return for the same reason) from a genuine error. Mirrors
#' 02-run-glow-chr.R:366-401.
#'
#' @param vset A \code{glow_variant_set} (uses \code{$G} and \code{$annotations}).
#' @param staar The validated STAAR context.
#' @param reference_medians Named numeric vector of annotation medians.
#' @return A list with \code{staar_out} (or NULL), \code{skipped} (logical), and
#'   \code{error} (a message string for a non-skip failure, else NULL).
#' @keywords internal
#' @noRd
.run_region_staar <- function(vset, staar, reference_medians) {
  anno_cols <- staar$anno_cols
  staar_anno <- as.matrix(vset$annotations[, anno_cols, drop = FALSE])

  # Impute missing annotation values with the chromosome-wide medians (the same
  # values prepare_glow_input uses for PI prediction).
  for (col_name in anno_cols) {
    na_idx <- is.na(staar_anno[, col_name])
    if (any(na_idx)) {
      staar_anno[na_idx, col_name] <- reference_medians[[col_name]]
    }
  }

  tryCatch({
    staar_out <- STAAR::STAAR(
      genotype         = vset$G,
      obj_nullmodel    = staar$null_model,
      annotation_phred = staar_anno,
      rare_maf_cutoff  = staar$rare_maf_cutoff
    )
    # NULL when STAAR has < rv_num_cutoff (default 2) qualifying variants.
    if (is.null(staar_out)) {
      list(staar_out = NULL, skipped = TRUE, error = NULL)
    } else {
      list(staar_out = staar_out, skipped = FALSE, error = NULL)
    }
  }, error = function(e) {
    em <- conditionMessage(e)
    # "Number of rare variant in the set is less than 2!" is STAAR's
    # design-intended exclusion, not a code error.
    if (grepl("less than 2", em, fixed = TRUE)) {
      list(staar_out = NULL, skipped = TRUE, error = NULL)
    } else {
      list(staar_out = NULL, skipped = FALSE, error = em)
    }
  })
}


#' Splice the 7 region-summary columns into a result row
#'
#' The per-region diagnostic columns, derived from the
#' \code{glow_input} this function already holds:
#' \itemize{
#'   \item \code{n_collapsed_units = sum(input$is_collapsed)};
#'   \item \code{min_MAF}/\code{max_MAF} over post-collapse MAF
#'     (\code{colMeans(input$G)/2});
#'   \item \code{min_MAC}/\code{max_MAC} over post-collapse MAC
#'     (\code{colSums(input$G)});
#'   \item \code{min_p_standard = min(2*pnorm(-abs(Z_standard)))} (NA when no
#'     units);
#'   \item \code{min_p_SPA = min(2*pnorm(-abs(Z_SPA)))} when \code{Z_SPA} is
#'     non-NULL (binary trait), else NA.
#' }
#' Spliced after \code{cMAC} and before the GLOW p-value columns.
#'
#' @param result_row A 1-row \code{$result} data.frame.
#' @param input A \code{glow_input} (uses \code{$G}, \code{$is_collapsed}).
#' @param Z_standard Standard-Z scores (length \code{n_after_collapse}).
#' @param Z_SPA SPA Z-scores (binary trait) or NULL (continuous / no SPA).
#' @return The \code{result_row} with 7 columns spliced after \code{cMAC}.
#' @keywords internal
#' @noRd
.splice_region_summary <- function(result_row, input, Z_standard, Z_SPA) {
  MAF <- colMeans(input$G) / 2
  MAC <- colSums(input$G)
  pvalue_standard <- 2 * stats::pnorm(-abs(Z_standard))
  pvalue_SPA <- if (!is.null(Z_SPA)) 2 * stats::pnorm(-abs(Z_SPA)) else NULL

  summary_df <- data.frame(
    n_collapsed_units = sum(input$is_collapsed),
    min_MAF           = min(MAF),
    max_MAF           = max(MAF),
    min_MAC           = min(MAC),
    max_MAC           = max(MAC),
    min_p_standard    = if (length(pvalue_standard) > 0L)
      min(pvalue_standard) else NA_real_,
    min_p_SPA         = if (!is.null(pvalue_SPA))
      min(pvalue_SPA) else NA_real_,
    stringsAsFactors = FALSE
  )

  insert_after <- match("cMAC", names(result_row))
  cbind(result_row[, seq_len(insert_after), drop = FALSE],
        summary_df,
        result_row[, seq.int(insert_after + 1L, ncol(result_row)),
                   drop = FALSE])
}


#' Build the per-region evidence packet
#'
#' Reconstructs the filter -> LD -> collapse traceability (which post-filter
#' variants survived LD pruning and which collapse group they joined) and
#' records per-unit B/PI/Z/weights.
#'
#' @param vset A \code{glow_variant_set}.
#' @param input A \code{glow_input} (from \code{prepare_glow_input}).
#' @param Z_standard Standard-Z scores (length \code{n_after_collapse}).
#' @param Z_SPA SPA Z-scores (binary trait) or NULL (continuous).
#' @param result A \code{glow_test_result} carrying \code{$weights}.
#' @param null_model The fitted \code{glow_null_model}.
#' @param used_spa Logical. The resolved SPA choice GLOW actually used for this
#'   region (recorded as \code{glow_test_used_spa}).
#' @return A named list (the evidence packet).
#' @keywords internal
#' @noRd
.build_region_evidence_packet <- function(vset, input, Z_standard, Z_SPA,
                                          result, null_model, used_spa) {
  n_filter <- input$n_after_filter

  ld_kept <- seq_len(n_filter) %in% input$ld_keep_idx
  collapse_group <- rep(NA_integer_, n_filter)
  for (k in seq_along(input$col_mapping)) {
    post_ld_idx     <- input$col_mapping[[k]]
    post_filter_idx <- input$ld_keep_idx[post_ld_idx]
    collapse_group[post_filter_idx] <- k
  }

  post_filter <- cbind(
    vset$variant_info,
    data.frame(ld_kept = ld_kept,
               collapse_group = collapse_group,
               stringsAsFactors = FALSE)
  )

  component_post_filter_idx <- lapply(input$col_mapping, function(post_ld_idx) {
    input$ld_keep_idx[post_ld_idx]
  })

  test_units <- list(
    is_collapsed              = input$is_collapsed,
    component_post_filter_idx = component_post_filter_idx,
    MAF                       = colMeans(input$G) / 2,
    MAC                       = colSums(input$G),
    B                         = input$B,
    PI                        = input$PI,
    Z_standard                = Z_standard,
    pvalue_standard           = 2 * stats::pnorm(-abs(Z_standard)),
    Z_SPA                     = Z_SPA,
    pvalue_SPA                = if (!is.null(Z_SPA))
      2 * stats::pnorm(-abs(Z_SPA)) else NULL,
    weights                   = result$weights
  )

  list(
    region             = input$region,
    n_total_in_region  = input$n_original,
    n_after_filter     = n_filter,
    n_after_ld         = input$n_after_ld,
    n_after_collapse   = input$n_after_collapse,
    cMAC               = input$cMAC,
    trait              = null_model$trait,
    glow_test_used_spa = used_spa,
    post_filter        = post_filter,
    test_units         = test_units
  )
}
