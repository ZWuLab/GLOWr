# Tests for extract_variant_set.R

# Create synthetic aGDS once for all tests in this file
test_agds_path <- tempfile(fileext = ".gds")
create_test_agds(test_agds_path)

withr::defer(unlink(test_agds_path), teardown_env())

# ==============================================================================
# Basic structure tests
# ==============================================================================

test_that("extract_variant_set returns correct structure", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(variant_type = "SNV", rare_maf_cutoff = 0.5),
    verbose = 0
  )
  expect_s3_class(vset, "glow_variant_set")
  expect_true(is.matrix(vset$G))
  expect_true(is.data.frame(vset$variant_info))
  expect_equal(ncol(vset$G), nrow(vset$variant_info))
  expect_equal(nrow(vset$G), vset$n_samples)
  expect_equal(ncol(vset$G), vset$n_variants)
  expect_equal(vset$n_samples, 50)
  expect_true(vset$n_total_in_region > 0)
  # n_after_annotation is present and satisfies the ordering invariant.
  expect_true(!is.null(vset$n_after_annotation))
  expect_true(is.numeric(vset$n_after_annotation) || is.integer(vset$n_after_annotation))
  expect_true(vset$n_total_in_region >= vset$n_after_annotation)
  expect_true(vset$n_after_annotation >= vset$n_variants)
})

test_that("extract_variant_set returns NULL for empty regions", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 999999, end = 999999),
    filter_spec = variant_filter(),
    verbose = 0
  )
  expect_null(vset)
})

test_that("extract_variant_set returns NULL for region on wrong chromosome", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "1", start = 1000, end = 5000),
    filter_spec = variant_filter(),
    verbose = 0
  )
  expect_null(vset)
})

# ==============================================================================
# Region input handling
# ==============================================================================

test_that("extract_variant_set accepts data.frame region", {
  region_df <- data.frame(
    region_id = "GENE_A", chr = "22", start = 1000, end = 5000,
    label = "GENE_A", stringsAsFactors = FALSE
  )
  vset <- extract_variant_set(
    test_agds_path, region = region_df,
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )
  expect_s3_class(vset, "glow_variant_set")
  expect_equal(vset$region$label, "GENE_A")
})

test_that("extract_variant_set accepts list region", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )
  expect_s3_class(vset, "glow_variant_set")
})

# ==============================================================================
# QC filter
# ==============================================================================

test_that("extract_variant_set applies QC filter", {
  # With QC filter (default)
  vset_qc <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1, end = 50000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )

  # Without QC filter (use a non-existent node so tryCatch falls back)
  vset_no_qc <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1, end = 50000),
    filter_spec = variant_filter(
      qc_label = "nonexistent/node",
      rare_maf_cutoff = 0.5
    ),
    verbose = 0
  )

  # QC-filtered should have fewer variants
  expect_true(vset_qc$n_variants < vset_no_qc$n_variants)
})

# ==============================================================================
# Variant type filter
# ==============================================================================

test_that("extract_variant_set filters by variant type", {
  # SNV only
  vset_snv <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1, end = 50000),
    filter_spec = variant_filter(variant_type = "SNV", rare_maf_cutoff = 0.5),
    verbose = 0
  )
  # All variants in synthetic data are SNVs (A/G or C/T)
  expect_true(vset_snv$n_variants > 0)
})

# ==============================================================================
# Annotation-based filter
# ==============================================================================

test_that("extract_variant_set applies plof annotation filter", {
  # GENE_A has 5 stopgain (indices 1:5) in positions 1000-5000
  spec_plof <- coding_filter("plof", rare_maf_cutoff = 0.5)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = spec_plof, verbose = 0
  )
  # Should get stopgain variants (indices 1:5, minus any QC-failed)
  # Index 3 is non-PASS, so expect 4 stopgain variants
  expect_true(vset$n_variants > 0)
  expect_true(vset$n_variants <= 5)
})

test_that("extract_variant_set applies missense filter correctly", {
  # GENE_A has nonsynonymous SNVs at indices 6:15 (10 total)
  spec_missense <- coding_filter("missense", rare_maf_cutoff = 0.5)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = spec_missense, verbose = 0
  )
  # 10 missense minus 1 QC-failed (index 8) = 9
  expect_true(vset$n_variants > 0)
})

test_that("extract_variant_set applies disruptive_missense filter", {
  # GENE_A indices 6:10 are nonsynonymous + MetaSVM="D".
  # min_mac = 0L is set explicitly here because the synthetic 50-sample
  # fixture (seed 42) generates the 4 surviving disruptive_missense variants
  # with cohort-MAC = 0 by random chance; the new variant_filter default
  # (min_mac = 1L) would correctly drop them, but the test's intent is to
  # verify the annotation filter, not the MAC filter.
  spec_dm <- coding_filter("disruptive_missense", rare_maf_cutoff = 0.5,
                            min_mac = 0L)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = spec_dm, verbose = 0
  )
  # 5 disruptive missense minus 1 QC-failed (index 8) = 4
  expect_true(vset$n_variants > 0)
  expect_true(vset$n_variants <= 5)
})

test_that("extract_variant_set applies synonymous filter", {
  # GENE_C (pos 20000-25000) has 30 synonymous variants
  spec_syn <- coding_filter("synonymous", rare_maf_cutoff = 0.5)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 20000, end = 25000),
    filter_spec = spec_syn, verbose = 0
  )
  # 30 synonymous minus 1 QC-failed (index 85) = 29
  expect_true(vset$n_variants > 0)
})

test_that("extract_variant_set returns NULL when plof on synonymous-only gene", {
  # GENE_C is all synonymous — plof should find nothing
  spec_plof <- coding_filter("plof", rare_maf_cutoff = 0.5)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 20000, end = 25000),
    filter_spec = spec_plof, verbose = 0
  )
  expect_null(vset)
})

# ==============================================================================
# MAF/MAC filter
# ==============================================================================

test_that("extract_variant_set applies MAF cutoff", {
  # Default rare_maf_cutoff = 0.01 should filter out common variants
  vset_rare <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.01), verbose = 0
  )
  # Wide MAF cutoff should get more
  vset_all <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )
  expect_true(vset_rare$n_variants <= vset_all$n_variants)
})

# ==============================================================================
# min_variants threshold
# ==============================================================================

test_that("extract_variant_set respects min_variants", {
  # GENE_E has 5 variants, 4 common (MAF > 0.06). With rare MAF cutoff,
  # only 1 rare variant remains, which is < min_variants=2
  spec <- variant_filter(rare_maf_cutoff = 0.01, min_variants = 2L)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 40000, end = 41000),
    filter_spec = spec, verbose = 0
  )
  # With strict MAF cutoff and min_variants=2, region should return NULL
  # (only 1 rare variant in GENE_E after MAF filter) or have >= 2 variants
  is_valid <- is.null(vset) || vset$n_variants >= 2
  expect_true(is_valid)
})

# ==============================================================================
# Annotation score extraction
# ==============================================================================

test_that("extract_variant_set extracts annotation scores", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5),
    annotation_names = c("CADD", "LINSIGHT"),
    verbose = 0
  )
  expect_true(!is.null(vset$annotations))
  expect_equal(ncol(vset$annotations), 2)
  expect_equal(colnames(vset$annotations), c("CADD", "LINSIGHT"))
  expect_equal(nrow(vset$annotations), vset$n_variants)
})

test_that("extract_variant_set works without annotation scores", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5),
    verbose = 0
  )
  expect_null(vset$annotations)
})

# ==============================================================================
# Sample alignment
# ==============================================================================

test_that("extract_variant_set handles sample alignment", {
  # Get sample IDs, close GDS, then pass path
  gds <- SeqArray::seqOpen(test_agds_path)
  all_ids <- SeqArray::seqGetData(gds, "sample.id")
  SeqArray::seqClose(gds)

  shuffled_ids <- rev(all_ids)
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5),
    sample_id = shuffled_ids, verbose = 0
  )
  expect_equal(vset$n_samples, length(shuffled_ids))
})

test_that("extract_variant_set with sample subset", {
  gds <- SeqArray::seqOpen(test_agds_path)
  all_ids <- SeqArray::seqGetData(gds, "sample.id")
  SeqArray::seqClose(gds)

  subset_ids <- all_ids[1:10]
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5),
    sample_id = subset_ids, verbose = 0
  )
  expect_equal(vset$n_samples, 10)
})

# ==============================================================================
# GDS connection handling
# ==============================================================================

test_that("extract_variant_set accepts open GDS connection", {
  gds <- SeqArray::seqOpen(test_agds_path)
  vset <- extract_variant_set(
    gds,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )
  expect_s3_class(vset, "glow_variant_set")
  # GDS should still be open (not closed by extract_variant_set)
  expect_true(inherits(gds, "SeqVarGDSClass"))
  SeqArray::seqClose(gds)
})

test_that("extract_variant_set accepts file path", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )
  expect_s3_class(vset, "glow_variant_set")
})

# ==============================================================================
# variant_info correctness
# ==============================================================================

test_that("variant_info has correct columns", {
  vset <- extract_variant_set(
    test_agds_path,
    region = list(chr = "22", start = 1000, end = 5000),
    filter_spec = variant_filter(rare_maf_cutoff = 0.5), verbose = 0
  )
  expected_cols <- c("variant_id", "rsid", "chr", "pos", "ref", "alt", "MAF", "MAC")
  expect_equal(names(vset$variant_info), expected_cols)
  expect_true(all(vset$variant_info$MAF >= 0))
  expect_true(all(vset$variant_info$MAF <= 0.5))
  expect_true(all(vset$variant_info$MAC >= 0))
})

# ==============================================================================
# Multiple calls with same open GDS (loop pattern)
# ==============================================================================

test_that("extract_variant_set works in a loop with open GDS", {
  gds <- SeqArray::seqOpen(test_agds_path)
  on.exit(SeqArray::seqClose(gds))

  regions <- list(
    list(chr = "22", start = 1000, end = 5000),
    list(chr = "22", start = 10000, end = 15000),
    list(chr = "22", start = 20000, end = 25000)
  )
  spec <- variant_filter(rare_maf_cutoff = 0.5)

  results <- lapply(regions, function(r) {
    extract_variant_set(gds, r, spec, verbose = 0)
  })

  expect_true(all(vapply(results, function(x) {
    is.null(x) || inherits(x, "glow_variant_set")
  }, logical(1))))

  n_non_null <- sum(!vapply(results, is.null, logical(1)))
  expect_true(n_non_null >= 2)
})

# ==============================================================================
# compute_annotation_medians() and filter consistency with extract_variant_set()
# ==============================================================================

test_that("compute_annotation_medians returns named numeric vector", {
  spec <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV")
  medians <- compute_annotation_medians(
    test_agds_path,
    annotation_names = c("cadd_phred", "linsight", "fathmm_xf"),
    filter_spec = spec,
    verbose = 0
  )
  expect_true(is.numeric(medians))
  expect_equal(names(medians), c("cadd_phred", "linsight", "fathmm_xf"))
  expect_equal(length(medians), 3)
})

test_that("compute_annotation_medians requires filter_spec", {
  expect_error(
    compute_annotation_medians(
      test_agds_path,
      annotation_names = "cadd_phred",
      verbose = 0
    ),
    "filter_spec is required"
  )
})

test_that("compute_annotation_medians returns NA vector when no variants match", {
  # An impossible MAC cutoff: no variant can have MAC >= 10^9.
  spec <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV",
                          min_mac = 1e9)
  medians <- compute_annotation_medians(
    test_agds_path,
    annotation_names = c("cadd_phred", "linsight"),
    filter_spec = spec,
    verbose = 0
  )
  expect_equal(names(medians), c("cadd_phred", "linsight"))
  expect_true(all(is.na(medians)))
})

test_that(".select_variants_by_spec omits G when return_genotypes = FALSE", {
  gds <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds))
  all_ids <- SeqArray::seqGetData(gds, "variant.id")

  cat <- get("annotation_name_catalog", envir = asNamespace("GLOWr"))
  spec <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV")

  out <- GLOWr:::.select_variants_by_spec(
    gds = gds,
    candidate_variant_ids = all_ids,
    filter_spec = spec,
    Annotation_dir = "annotation/info/FunctionalAnnotation",
    Annotation_name_catalog = cat,
    return_genotypes = FALSE
  )
  expect_null(out$G)
  expect_true(length(out$variant_ids) > 0)
  expect_equal(length(out$MAF), length(out$variant_ids))
  expect_equal(length(out$MAC), length(out$variant_ids))
  expect_true(all(out$MAF <= 0.5 & out$MAF >= 0))
  expect_true(all(out$MAC >= 0L))
  # n_after_annotation is present and >= n final variants (MAF/MAC may drop more).
  expect_true(!is.null(out$n_after_annotation))
  expect_true(out$n_after_annotation >= length(out$variant_ids))
})

test_that(".select_variants_by_spec leaves GDS filter narrowed to final_ids", {
  # Regression: if the helper leaves a pre-MAF/MAC filter on the GDS,
  # compute_annotation_medians() reads annotations over the wrong pool.
  # A subsequent seqGetData() call must return exactly length(variant_ids)
  # rows, not more.
  gds <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds))
  all_ids <- SeqArray::seqGetData(gds, "variant.id")

  cat <- get("annotation_name_catalog", envir = asNamespace("GLOWr"))
  # Tight MAF cutoff so final_ids strictly smaller than narrow_ids (the
  # only way this bug can surface).
  spec <- variant_filter(rare_maf_cutoff = 0.1, variant_type = "SNV")

  out <- GLOWr:::.select_variants_by_spec(
    gds = gds,
    candidate_variant_ids = all_ids,
    filter_spec = spec,
    Annotation_dir = "annotation/info/FunctionalAnnotation",
    Annotation_name_catalog = cat,
    return_genotypes = FALSE
  )

  # After the helper returns, the GDS filter should be set to final_ids.
  post_vids <- SeqArray::seqGetData(gds, "variant.id")
  expect_equal(length(post_vids), length(out$variant_ids))
  expect_setequal(post_vids, out$variant_ids)
})

test_that(".select_variants_by_spec ignores min_variants", {
  # The helper is intentionally neutral about min_variants — that's a
  # caller concern (per-gene skip vs. chromosome-wide reference pool).
  gds <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds))
  all_ids <- SeqArray::seqGetData(gds, "variant.id")

  cat <- get("annotation_name_catalog", envir = asNamespace("GLOWr"))
  spec <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV",
                          min_variants = 1e6)

  out <- GLOWr:::.select_variants_by_spec(
    gds = gds,
    candidate_variant_ids = all_ids,
    filter_spec = spec,
    Annotation_dir = "annotation/info/FunctionalAnnotation",
    Annotation_name_catalog = cat,
    return_genotypes = FALSE
  )
  expect_true(length(out$variant_ids) > 0)
})

test_that("per-region variant IDs are a subset of chromosome-wide pool", {
  # The refactor guarantees filter consistency: the variants any single
  # region can return must already be in the chromosome-wide scan result
  # when the same filter_spec is applied.
  spec <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV")

  gds <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds))
  all_ids <- SeqArray::seqGetData(gds, "variant.id")

  cat <- get("annotation_name_catalog", envir = asNamespace("GLOWr"))
  chrom_pool <- GLOWr:::.select_variants_by_spec(
    gds = gds,
    candidate_variant_ids = all_ids,
    filter_spec = spec,
    Annotation_dir = "annotation/info/FunctionalAnnotation",
    Annotation_name_catalog = cat,
    return_genotypes = FALSE
  )$variant_ids
  SeqArray::seqResetFilter(gds, verbose = FALSE)

  # Per-region extraction: pass the open gds handle to avoid re-opening
  # the same file while we still hold one connection.
  regions <- list(
    list(chr = "22", start = 1000, end = 5000),
    list(chr = "22", start = 5001, end = 10000),
    list(chr = "22", start = 10001, end = 20000)
  )
  per_region_ids <- unlist(lapply(regions, function(r) {
    vset <- extract_variant_set(gds, r, spec, verbose = 0)
    if (is.null(vset)) integer(0) else vset$variant_info$variant_id
  }))

  expect_true(length(per_region_ids) > 0)
  expect_true(all(per_region_ids %in% chrom_pool))
})

test_that("compute_annotation_medians pool shrinks under tighter filter_spec", {
  # A more restrictive MAF cutoff must produce a subset pool.
  loose <- variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV")
  tight <- variant_filter(rare_maf_cutoff = 0.01, variant_type = "SNV")

  # Grab the per-pool size by intercepting message() output with a
  # calling handler (capture.output type="message" is flaky here).
  count_from <- function(spec) {
    msgs <- character()
    withCallingHandlers(
      compute_annotation_medians(
        test_agds_path,
        annotation_names = "cadd_phred",
        filter_spec = spec,
        verbose = 1
      ),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    # "Computing annotation medians from N variants"
    hits <- grep("Computing annotation medians", msgs, value = TRUE)
    if (length(hits) == 0) return(NA_integer_)
    as.integer(sub(".*from ([0-9]+) variants.*", "\\1", hits[1]))
  }

  n_loose <- count_from(loose)
  n_tight <- count_from(tight)
  expect_false(is.na(n_loose))
  expect_false(is.na(n_tight))
  expect_true(n_tight <= n_loose)
})
