# ==============================================================================
# Tests for glow_region_test() — the single-region GLOW primitive
# ==============================================================================
#
# Uses the synthetic test aGDS (helper-create-test-agds.R) plus a small,
# prediction-compatible GLM PI ensemble, a binary null model, and (where
# installed) a STAAR context. Regions are chosen deterministically from the
# fixture so statuses are reproducible:
#   - GENE_A [1000, 5000]  : multi-variant (n_after_collapse = 4)  -> "ok"
#   - [1000, 1000]         : exactly one post-filter variant       -> single-SNV
#   - [999990, 999999]     : no variants                           -> "skip_empty"


# ---- Shared fixtures (built once for the file) ----

test_agds_path <- tempfile(fileext = ".gds")
create_test_agds(test_agds_path)
withr::defer(unlink(test_agds_path), teardown_env())

# Open one connection reused across tests (the primitive's intended contract).
.gds <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
withr::defer(SeqArray::seqClose(.gds), teardown_env())

.sids   <- SeqArray::seqGetData(.gds, "sample.id")
.n      <- length(.sids)
.feats  <- c("cadd_phred", "linsight", "fathmm_xf")
.spec   <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV",
                          min_mac = 1L, min_variants = 1L)

# Binary null model on the GDS sample order.
set.seed(20260526)
.X  <- cbind(1, rnorm(.n))
.Y  <- rbinom(.n, 1, 0.3)
.nm <- fit_null_model(.X, .Y, trait = "binary", sample_id = .sids)

# Chromosome-wide annotation medians (for PI + STAAR NA imputation).
.meds <- compute_annotation_medians(.gds, .feats, .spec,
                                    sample_id = .sids, verbose = 0)

# Small prediction-compatible GLM PI ensemble (matches load_PI_models() shape).
.make_glm_pi <- function(seed) {
  set.seed(seed)
  xx <- matrix(runif(300 * length(.feats)), 300, length(.feats))
  yy <- rbinom(300, 1, plogis(xx %*% rep(0.4, length(.feats))))
  glm(y ~ x, data = list(x = xx, y = yy), family = "binomial")
}
.pi_models <- list(models = list(.make_glm_pi(11), .make_glm_pi(12)),
                   model_type = "GLM", n_models = 2L)

# B sources: a closed-form function and an equivalent trained-model stand-in.
.b_func <- function(maf) rep(0.3, length(maf))
.make_b_model <- function() {
  set.seed(7)
  maf  <- runif(60, 0.01, 0.3)
  logX <- log(maf)
  Y    <- -0.14 * logX + 0.5 + rnorm(60, 0, 0.02)
  lm_obj <- lm(Y ~ logX, data = data.frame(Y = Y, logX = logX))
  attr(lm_obj, "model_id") <- 3L  # Y ~ log(MAF); predict_B() reconstructs design
  lm_obj
}

# STAAR context (only used by skip_if_not_installed tests).
.make_staar_ctx <- function() {
  sdat  <- data.frame(phenotype = .Y, cov1 = .X[, 2])
  snull <- STAAR::fit_null_glm(phenotype ~ cov1, data = sdat,
                               family = stats::binomial(link = "logit"))
  list(null_model = snull, anno_cols = .feats, rare_maf_cutoff = 0.5)
}

# Convenience wrapper with the fixed fixtures bound.
.call_region <- function(region, staar = NULL, return_evidence = FALSE,
                         b_func = .b_func, b_model = NULL, null_model = .nm,
                         z_scale = 1, use_spa = NULL, region_summary = FALSE) {
  glow_region_test(.gds, region, .spec, null_model,
                   pi_models = .pi_models, pi_features = .feats,
                   reference_medians = .meds,
                   b_func = b_func, b_model = b_model,
                   ld_threshold = 0.95, mac_threshold = 11L,
                   collapse_method = "mean", use_spa = use_spa, z_scale = z_scale,
                   staar = staar, return_evidence = return_evidence,
                   region_summary = region_summary, verbose = 0)
}

.region_multi  <- list(chr = "22", start = 1000L,   end = 5000L,   label = "GENE_A")
.region_single <- list(chr = "22", start = 1000L,   end = 1000L,   label = "SINGLE")
.region_empty  <- list(chr = "22", start = 999990L, end = 999999L, label = "EMPTY")

.GLOW_COLS <- c("GLOW_SKAT_BE_N", "GLOW_SKAT_APR_N", "GLOW_SKAT_BE_sparse",
                "GLOW_SKAT_APR_sparse", "GLOW_SKAT_equ",
                "GLOW_Burden_BE", "GLOW_Burden_APR", "GLOW_Burden_equ",
                "GLOW_Fisher_BE_N", "GLOW_Fisher_APR_N", "GLOW_Fisher_BE_sparse",
                "GLOW_Fisher_APR_sparse", "GLOW_Fisher_equ",
                "GLOW_BSF_Omni", "GLOW_SNV_CCT", "GLOW_Omni")


# ---- Happy path, GLOW only ----

test_that("multi-variant region returns status 'ok' with the canonical row", {
  res <- .call_region(.region_multi)

  expect_equal(res$status, "ok")
  expect_null(res$staar)
  expect_null(res$staar_detail)
  expect_null(res$evidence)          # return_evidence = FALSE
  expect_null(res$message)

  expect_s3_class(res$result, "data.frame")
  expect_equal(nrow(res$result), 1L)
  id_count_cols <- c("label", "chr", "start", "end", "n_variants",
                     "n_after_annotation", "n_after_filter",
                     "n_after_ld", "n_after_collapse", "cMAC")
  expect_true(all(id_count_cols %in% names(res$result)))
  expect_true(all(.GLOW_COLS %in% names(res$result)))
  expect_equal(res$result$label, "GENE_A")
  expect_true(res$result$n_after_collapse >= 2L)   # genuinely multi-variant
  # Ordering invariant: n_variants >= n_after_annotation >= n_after_filter.
  expect_true(res$result$n_variants >= res$result$n_after_annotation)
  expect_true(res$result$n_after_annotation >= res$result$n_after_filter)
  # No spliced gene-driver diagnostic columns in the primitive's canonical row.
  expect_false("min_MAF" %in% names(res$result))
  expect_false("n_collapsed_units" %in% names(res$result))
  # GLOW p-values are valid probabilities.
  pv <- unlist(res$result[.GLOW_COLS])
  expect_true(all(pv >= 0 & pv <= 1))
})


# ---- return_evidence = TRUE ----

test_that("return_evidence = TRUE attaches a packet of the gene shape", {
  res <- .call_region(.region_multi, return_evidence = TRUE)

  expect_equal(res$status, "ok")
  expect_false(is.null(res$evidence))
  ev_fields <- c("region", "n_after_filter", "n_after_ld", "n_after_collapse",
                 "cMAC", "post_filter", "test_units")
  expect_true(all(ev_fields %in% names(res$evidence)))
  expect_s3_class(res$evidence$post_filter, "data.frame")
  expect_true(all(c("ld_kept", "collapse_group") %in%
                  names(res$evidence$post_filter)))
  # Per-unit vectors are length n_after_collapse.
  k <- res$result$n_after_collapse
  expect_length(res$evidence$test_units$B, k)
  expect_length(res$evidence$test_units$PI, k)
  # Binary trait -> SPA channel populated.
  expect_true(res$evidence$glow_test_used_spa)
  expect_false(is.null(res$evidence$test_units$Z_SPA))
})


# ---- Empty region ----

test_that("region with no qualifying variants is skip_empty", {
  res <- .call_region(.region_empty)
  expect_equal(res$status, "skip_empty")
  expect_null(res$result)
  expect_null(res$evidence)
  expect_null(res$staar)
})


# ---- Single-variant region: GLOW p==1 path, and staar_skip_single ----

test_that("single-variant region: all 16 GLOW p-values are equal (GLOW only)", {
  res <- .call_region(.region_single)
  expect_equal(res$status, "ok")
  expect_equal(res$result$n_after_collapse, 1L)
  pv <- unlist(res$result[.GLOW_COLS])
  expect_equal(length(unique(signif(pv, 12))), 1L)  # marginal p assigned to all
})

test_that("single-variant region with STAAR context yields staar_skip_single + NA STAAR", {
  skip_if_not_installed("STAAR")
  res <- .call_region(.region_single, staar = .make_staar_ctx())

  expect_equal(res$status, "staar_skip_single")
  expect_false(is.null(res$result))                 # GLOW row IS kept
  expect_equal(res$result$n_after_collapse, 1L)
  expect_false(is.null(res$staar))                  # STAAR enabled -> list of 8
  expect_length(res$staar, 8L)
  expect_true(all(is.na(unlist(res$staar))))        # all NA on skip
  expect_null(res$staar_detail)
})


# ---- B-source dispatch (sanity, not bit-equality across sources) ----

test_that("b_func and b_model paths both produce a valid result", {
  res_func  <- .call_region(.region_multi, b_func = .b_func, b_model = NULL)
  res_model <- .call_region(.region_multi, b_func = NULL,
                            b_model = .make_b_model())

  expect_equal(res_func$status,  "ok")
  expect_equal(res_model$status, "ok")
  expect_true(is.finite(res_func$result$GLOW_Omni))
  expect_true(is.finite(res_model$result$GLOW_Omni))
  # Same region geometry / variant counts regardless of B source.
  expect_equal(res_func$result$n_after_collapse,
               res_model$result$n_after_collapse)
})

test_that("supplying both or neither B source is a caller error", {
  expect_error(.call_region(.region_multi, b_func = .b_func,
                            b_model = .make_b_model()),
               "Exactly one")
  expect_error(.call_region(.region_multi, b_func = NULL, b_model = NULL),
               "Exactly one")
})


# ---- STAAR context happy path ----

test_that("multi-variant region with STAAR context returns non-NA STAAR_O", {
  skip_if_not_installed("STAAR")
  res <- .call_region(.region_multi, staar = .make_staar_ctx())

  expect_equal(res$status, "ok")
  expect_false(is.null(res$staar))
  expect_length(res$staar, 8L)
  expect_named(res$staar, c("STAAR_O", "ACAT_O",
                            "STAAR_S_1_25", "STAAR_S_1_1",
                            "STAAR_B_1_25", "STAAR_B_1_1",
                            "STAAR_A_1_25", "STAAR_A_1_1"))
  expect_false(is.na(res$staar$STAAR_O))
  expect_true(res$staar$STAAR_O >= 0 && res$staar$STAAR_O <= 1)
  expect_false(is.null(res$staar_detail))
})

test_that("STAAR context is rejected when malformed", {
  skip_if_not_installed("STAAR")
  bad_ctx <- list(null_model = NULL, anno_cols = .feats)  # missing rare_maf_cutoff
  expect_error(.call_region(.region_multi, staar = bad_ctx),
               "rare_maf_cutoff")
})


# ---- Error surface ----

test_that("a malformed null_model surfaces as status 'error', not a throw", {
  nm_bad <- .nm
  nm_bad$n <- 999L   # forces a dimension-mismatch stop() inside compute_score_stats

  res <- expect_no_error(.call_region(.region_multi, null_model = nm_bad))
  expect_equal(res$status, "error")
  expect_null(res$result)
  expect_true(is.character(res$message) && nzchar(res$message))
})


test_that("z_scale = 1 is a bit-for-bit no-op (GLOW p-values and evidence Z)", {
  base <- .call_region(.region_multi, return_evidence = TRUE)              # default z_scale = 1
  z1   <- .call_region(.region_multi, return_evidence = TRUE, z_scale = 1)
  expect_equal(z1$result[, .GLOW_COLS], base$result[, .GLOW_COLS])
  expect_equal(z1$evidence$test_units$Z_standard, base$evidence$test_units$Z_standard)
})

test_that("z_scale multiplies each unit's Z by exactly z_scale (evidence packet)", {
  cc   <- 0.6
  base <- .call_region(.region_multi, return_evidence = TRUE)
  zc   <- .call_region(.region_multi, return_evidence = TRUE, z_scale = cc)
  expect_equal(zc$evidence$test_units$Z_standard,
               cc * base$evidence$test_units$Z_standard)
  # Binary trait: the SPA Z is the test-driving Z and is scaled too.
  if (!is.null(base$evidence$test_units$Z_SPA)) {
    expect_equal(zc$evidence$test_units$Z_SPA, cc * base$evidence$test_units$Z_SPA)
  }
})

test_that("shrinking Z (z_scale < 1) makes GLOW tests no more significant", {
  base <- .call_region(.region_multi)
  shr  <- .call_region(.region_multi, z_scale = 0.5)   # halve every unit's Z
  expect_gte(shr$result$GLOW_Omni,      base$result$GLOW_Omni)
  expect_gte(shr$result$GLOW_Burden_equ, base$result$GLOW_Burden_equ)
})

test_that("invalid z_scale is rejected", {
  expect_error(.call_region(.region_multi, z_scale = 0),         "z_scale")
  expect_error(.call_region(.region_multi, z_scale = -1),        "z_scale")
  expect_error(.call_region(.region_multi, z_scale = c(1, 2)),   "z_scale")
  expect_error(.call_region(.region_multi, z_scale = NA_real_),  "z_scale")
})


test_that("use_spa = NULL (default) equals use_spa = TRUE on a binary trait", {
  # The binary auto-default IS SPA, so NULL and TRUE must be bit-identical.
  d  <- .call_region(.region_multi)                 # default NULL -> auto SPA
  tt <- .call_region(.region_multi, use_spa = TRUE)
  expect_equal(d$result[, .GLOW_COLS], tt$result[, .GLOW_COLS])
})

test_that("use_spa = FALSE switches the binary dispatch to standard score stats", {
  spa <- .call_region(.region_multi, use_spa = TRUE,  return_evidence = TRUE)
  std <- .call_region(.region_multi, use_spa = FALSE, return_evidence = TRUE)

  # Precondition: SPA and standard single-variant Z genuinely differ here.
  expect_false(isTRUE(all.equal(std$evidence$test_units$Z_standard,
                                std$evidence$test_units$Z_SPA)))
  # Consequence: the combined GLOW test now reads the non-SPA Z -> p differs.
  expect_false(isTRUE(all.equal(spa$result[, .GLOW_COLS],
                                std$result[, .GLOW_COLS])))

  # The evidence flag records the actual choice (not the trait).
  expect_true(spa$evidence$glow_test_used_spa)
  expect_false(std$evidence$glow_test_used_spa)

  # Both runs carry the SAME dual single-variant Z regardless of which fed the
  # test (the dispatch + its complement are the same two vectors).
  expect_equal(std$evidence$test_units$Z_standard,
               spa$evidence$test_units$Z_standard)
  expect_equal(std$evidence$test_units$Z_SPA,
               spa$evidence$test_units$Z_SPA)
})

test_that("evidence packet keeps both single-variant p-values when use_spa = FALSE", {
  std <- .call_region(.region_multi, use_spa = FALSE, return_evidence = TRUE)
  expect_false(is.null(std$evidence$test_units$pvalue_standard))
  expect_false(is.null(std$evidence$test_units$pvalue_SPA))   # SPA complement built
})

test_that("use_spa = FALSE composes with z_scale (Z_standard non-NULL invariant)", {
  std    <- .call_region(.region_multi, use_spa = FALSE, return_evidence = TRUE)
  std_zc <- .call_region(.region_multi, use_spa = FALSE, return_evidence = TRUE,
                         z_scale = 0.5)
  # No error (Z_standard always set), and z_scale shrinks the test-driving Z.
  expect_equal(std_zc$evidence$test_units$Z_standard,
               0.5 * std$evidence$test_units$Z_standard)
  expect_gte(std_zc$result$GLOW_Omni, std$result$GLOW_Omni)
})

test_that("use_spa = TRUE errors on a continuous trait", {
  set.seed(99)
  nm_cont <- fit_null_model(.X, rnorm(.n), trait = "continuous",
                            sample_id = .sids)
  expect_error(.call_region(.region_multi, null_model = nm_cont, use_spa = TRUE),
               "binary")
})

test_that("invalid use_spa is rejected", {
  expect_error(.call_region(.region_multi, use_spa = "yes"),         "use_spa")
  expect_error(.call_region(.region_multi, use_spa = c(TRUE, FALSE)), "use_spa")
  expect_error(.call_region(.region_multi, use_spa = NA),            "use_spa")
})


# ---- region_summary: the 7 gene-driver summary columns (Q1a) ----

.SUMMARY_COLS <- c("n_collapsed_units", "min_MAF", "max_MAF",
                   "min_MAC", "max_MAC", "min_p_standard", "min_p_SPA")

test_that("region_summary = FALSE (default) leaves the canonical row unchanged", {
  res <- .call_region(.region_multi)                       # default FALSE
  expect_false(any(.SUMMARY_COLS %in% names(res$result)))
})

test_that("region_summary = TRUE splices the 7 columns after cMAC, before GLOW_*", {
  res <- .call_region(.region_multi, region_summary = TRUE)
  expect_equal(res$status, "ok")
  expect_true(all(.SUMMARY_COLS %in% names(res$result)))

  nms <- names(res$result)
  # Spliced contiguously immediately after cMAC, immediately before the GLOW_*.
  expect_equal(which(nms == "n_collapsed_units"), which(nms == "cMAC") + 1L)
  expect_equal(which(nms == .GLOW_COLS[1L]), which(nms == "min_p_SPA") + 1L)
})

test_that("region_summary values bit-match a hand replication of the gene-driver formulas", {
  # Recompute the gene driver's exact path (02-run-glow-chr.R:170-183, :292-306):
  #   post-collapse MAF/MAC from input$G; per-unit single-variant min p-values
  #   2*pnorm(-abs(Z)) on standard (Z_standard) and SPA (Z_SPA, the binary
  #   dispatch); n_collapsed_units = sum(input$is_collapsed).
  vset  <- extract_variant_set(.gds, .region_multi, .spec,
                               sample_id = .nm$sample_id,
                               annotation_names = .feats, verbose = 0)
  input <- prepare_glow_input(vset, B_func = .b_func, PI_models = .pi_models,
                              PI_features = .feats, reference_medians = .meds,
                              ld_threshold = 0.95, mac_threshold = 11L,
                              collapse_method = "mean", verbose = 0)
  Z_SPA      <- compute_score_stats(input$G, .nm, verbose = 0)$Zscores  # binary dispatch
  Z_standard <- compute_score_stats(input$G, .nm, use_spa = FALSE,
                                    verbose = 0)$Zscores
  MAF <- colMeans(input$G) / 2
  MAC <- colSums(input$G)
  expected <- c(
    n_collapsed_units = sum(input$is_collapsed),
    min_MAF = min(MAF), max_MAF = max(MAF),
    min_MAC = min(MAC), max_MAC = max(MAC),
    min_p_standard = min(2 * pnorm(-abs(Z_standard))),
    min_p_SPA      = min(2 * pnorm(-abs(Z_SPA)))
  )

  res <- .call_region(.region_multi, region_summary = TRUE)
  got <- as.numeric(res$result[, .SUMMARY_COLS])
  expect_equal(got, as.numeric(expected))
})

test_that("min_p_SPA is NA for a continuous trait (no SPA channel)", {
  set.seed(7)
  nm_cont <- fit_null_model(.X, rnorm(.n), trait = "continuous",
                            sample_id = .sids)
  res <- .call_region(.region_multi, null_model = nm_cont, region_summary = TRUE)
  expect_equal(res$status, "ok")
  expect_true(is.na(res$result$min_p_SPA))
  expect_false(is.na(res$result$min_p_standard))
})
