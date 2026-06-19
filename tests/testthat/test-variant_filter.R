# Tests for variant_filter.R

test_that("variant_filter creates valid spec", {
  spec <- variant_filter(variant_type = "SNV", rare_maf_cutoff = 0.05)
  expect_s3_class(spec, "glow_variant_filter")
  expect_equal(spec$variant_type, "SNV")
  expect_equal(spec$rare_maf_cutoff, 0.05)
  expect_null(spec$annotation_clauses)
})

test_that("variant_filter validates inputs", {
  expect_error(variant_filter(variant_type = "invalid"))
  expect_error(variant_filter(rare_maf_cutoff = -1))
  expect_error(variant_filter(annotation_clauses = list(list())))
})

test_that("variant_filter stores annotation clauses", {
  clauses <- list(
    list("GENCODE.EXONIC.Category" = "stopgain"),
    list("GENCODE.Category" = "splicing")
  )
  spec <- variant_filter(annotation_clauses = clauses)
  expect_equal(length(spec$annotation_clauses), 2)
  expect_true("GENCODE.EXONIC.Category" %in% spec$required_annotation_fields)
  expect_true("GENCODE.Category" %in% spec$required_annotation_fields)
})

test_that("coding_filter produces correct masks for all 7 categories", {
  categories <- c("plof", "plof_ds", "missense", "disruptive_missense",
                   "synonymous", "ptv", "ptv_ds")
  for (cat in categories) {
    spec <- coding_filter(cat)
    expect_s3_class(spec, "glow_variant_filter")
    expect_equal(attr(spec, "category_name"), cat)
    expect_true(length(spec$annotation_clauses) > 0)
  }
})

test_that("coding_filter plof has correct fields", {
  spec <- coding_filter("plof")
  fields <- spec$required_annotation_fields
  expect_true("GENCODE.EXONIC.Category" %in% fields)
  expect_true("GENCODE.Category" %in% fields)
  expect_false("MetaSVM" %in% fields)
})

test_that("coding_filter plof has 6 clauses", {
  spec <- coding_filter("plof")
  expect_equal(length(spec$annotation_clauses), 6)
})

test_that("coding_filter disruptive_missense requires MetaSVM", {
  spec <- coding_filter("disruptive_missense")
  expect_true("MetaSVM" %in% spec$required_annotation_fields)
  expect_equal(length(spec$annotation_clauses), 1)
})

test_that("coding_filter plof_ds is plof + disruptive_missense", {
  plof <- coding_filter("plof")
  plof_ds <- coding_filter("plof_ds")
  expect_equal(length(plof_ds$annotation_clauses),
               length(plof$annotation_clauses) + 1)
  expect_true("MetaSVM" %in% plof_ds$required_annotation_fields)
})

test_that("coding_filter ptv has 4 clauses (no ncRNA)", {
  spec <- coding_filter("ptv")
  expect_equal(length(spec$annotation_clauses), 4)
  # Check that ncRNA categories are NOT present
  cats <- unlist(lapply(spec$annotation_clauses, function(cl) cl[["GENCODE.Category"]]))
  expect_false(any(grepl("ncRNA", cats)))
})

test_that("coding_filter ptv_ds is ptv + disruptive_missense", {
  ptv <- coding_filter("ptv")
  ptv_ds <- coding_filter("ptv_ds")
  expect_equal(length(ptv_ds$annotation_clauses),
               length(ptv$annotation_clauses) + 1)
})

test_that("coding_filter rejects unknown categories", {
  expect_error(coding_filter("nonexistent"), "Unknown coding category")
})

test_that("print.glow_variant_filter works", {
  spec <- coding_filter("plof")
  expect_output(print(spec), "plof")
  expect_output(print(spec), "GLOWr Variant Filter")
})

test_that("print.glow_variant_filter shows no-annotation case", {
  spec <- variant_filter()
  expect_output(print(spec), "none")
})

test_that("variant_filter / coding_filter default min_mac is 1L", {
  # Variants monomorphic in the analysis cohort (cohort-MAC = 0) carry no
  # information for any association test; the default excludes them. See
  # the 2026-05-14 WGS-run closeout for the motivating incident.
  expect_equal(variant_filter()$min_mac, 1L)
  expect_equal(coding_filter("plof")$min_mac, 1L)
  # Explicit override remains supported.
  expect_equal(variant_filter(min_mac = 0L)$min_mac, 0L)
  expect_equal(variant_filter(min_mac = 5L)$min_mac, 5L)
})
