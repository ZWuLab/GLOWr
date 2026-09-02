# Tests for annotation_predicate() and predicate evaluation in variant filtering
#
# These cover the 2026-08-21 filter-grammar extension: a clause condition may be
# an annotation_predicate() (nonempty/empty, in/not_in, gt/ge/lt/le) as well as
# an atomic vector of accepted values. The backward-compatibility tests are the
# important half -- every existing config and built-in coding mask uses the
# atomic %in% path and must be untouched.

test_agds_pred_path <- tempfile(fileext = ".gds")
create_test_agds(test_agds_pred_path)
withr::defer(unlink(test_agds_pred_path), teardown_env())

# The fixture's whole-chromosome region (125 variants, chr 22).
pred_region <- list(chr = "22", start = 1L, end = 1e9, label = "ALL")

# Extract with a single-clause filter, returning n_after_annotation (the count
# after QC + type + annotation clause, before MAF/MAC) so the tests isolate the
# clause from frequency filtering.
n_anno <- function(clauses) {
  vset <- extract_variant_set(
    test_agds_pred_path, region = pred_region,
    filter_spec = variant_filter(annotation_clauses = clauses,
                                 rare_maf_cutoff = 0.5, min_mac = 0L,
                                 min_variants = 1L),
    verbose = 0
  )
  if (is.null(vset)) 0L else vset$n_after_annotation
}

# Ground truth read straight from the fixture, so the expectations below are
# COMPUTED from the data rather than hard-coded. The clause is applied after the
# QC and variant-type filters, so every expected count is taken on that subset.
.fixture_truth <- local({
  g <- SeqArray::seqOpen(test_agds_pred_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(g), add = TRUE)
  keep <- (SeqArray::seqGetData(g, "annotation/filter") == "PASS") &
    GLOWr:::.is_snv_from_alleles(SeqArray::seqGetData(g, "allele"))
  list(
    keep = keep,
    n_keep = sum(keep),
    category = SeqArray::seqGetData(
      g, "annotation/info/FunctionalAnnotation/genecode_comprehensive_category"),
    cadd = SeqArray::seqGetData(
      g, "annotation/info/FunctionalAnnotation/cadd_phred"),
    genehancer = SeqArray::seqGetData(
      g, "annotation/info/FunctionalAnnotation/genehancer")
  )
})

# Count of QC-passing SNVs whose category is in `vals`.
.n_cat <- function(vals) {
  sum(.fixture_truth$keep & .fixture_truth$category %in% vals)
}

# ==============================================================================
# Constructor
# ==============================================================================

test_that("annotation_predicate builds and prints valid predicates", {
  p <- annotation_predicate("nonempty")
  expect_s3_class(p, "glow_annotation_predicate")
  expect_equal(p$op, "nonempty")
  expect_null(p$value)
  expect_output(print(p), "nonempty")

  p2 <- annotation_predicate("gt", 20)
  expect_equal(p2$op, "gt")
  expect_equal(p2$value, 20)
  expect_output(print(p2), "gt")

  expect_equal(annotation_predicate("in", c("D", "T"))$value, c("D", "T"))
})

test_that("annotation_predicate rejects malformed input", {
  expect_error(annotation_predicate("between", 1), "op must be one of")
  expect_error(annotation_predicate("nonempty", 1), "takes no value")
  expect_error(annotation_predicate("empty", "x"), "takes no value")
  expect_error(annotation_predicate("gt"), "single non-missing numeric")
  expect_error(annotation_predicate("gt", c(1, 2)), "single non-missing numeric")
  expect_error(annotation_predicate("gt", NA_real_), "single non-missing numeric")
  expect_error(annotation_predicate("gt", "20"), "single non-missing numeric")
  expect_error(annotation_predicate("in"), "non-empty atomic vector")
  expect_error(annotation_predicate("in", character(0)), "non-empty atomic vector")
})

test_that("variant_filter rejects a clause condition that is neither vector nor predicate", {
  expect_error(
    variant_filter(annotation_clauses = list(list("CADD" = NULL))),
    "atomic vector of accepted values or an annotation_predicate"
  )
  expect_error(
    variant_filter(annotation_clauses = list(list("CADD" = list(op = "gt")))),
    "atomic vector of accepted values or an annotation_predicate"
  )
  # A valid predicate passes validation and is retained in the spec.
  spec <- variant_filter(
    annotation_clauses = list(list("CADD" = annotation_predicate("gt", 20)))
  )
  expect_s3_class(spec$annotation_clauses[[1]][["CADD"]],
                  "glow_annotation_predicate")
  expect_equal(spec$required_annotation_fields, "CADD")
})

# ==============================================================================
# Backward compatibility -- the atomic %in% path is unchanged
# ==============================================================================

test_that("atomic conditions behave exactly as before", {
  expect_equal(n_anno(list(list("GENCODE.Category" = "upstream"))),
               .n_cat("upstream"))
  expect_equal(n_anno(list(list("GENCODE.Category" = "downstream"))),
               .n_cat("downstream"))
  # Multi-value membership OR-s within one field.
  expect_equal(
    n_anno(list(list("GENCODE.Category" = c("upstream", "downstream")))),
    .n_cat(c("upstream", "downstream")))
  # Clauses are OR-ed: two single-value clauses == one two-value clause.
  expect_equal(
    n_anno(list(list("GENCODE.Category" = "upstream"),
                list("GENCODE.Category" = "downstream"))),
    .n_cat(c("upstream", "downstream")))
  # Sanity: the fixture actually contains both categories, so the test has bite.
  expect_gt(.n_cat("upstream"), 0L)
  expect_gt(.n_cat("downstream"), 0L)
})

test_that("built-in coding masks still resolve and are unaffected", {
  for (cat in c("plof", "plof_ds", "missense", "disruptive_missense",
                "synonymous", "ptv", "ptv_ds")) {
    spec <- coding_filter(cat, rare_maf_cutoff = 0.5, min_mac = 0L,
                          min_variants = 1L)
    expect_s3_class(spec, "glow_variant_filter")
    vset <- extract_variant_set(test_agds_pred_path, region = pred_region,
                                filter_spec = spec, verbose = 0)
    # Every mask must run; the fixture is built so each selects >= 1 variant.
    expect_false(is.null(vset), info = cat)
  }
})

# ==============================================================================
# Predicate evaluation
# ==============================================================================

test_that("nonempty and empty partition the QC-passing variants", {
  n_nonempty <- n_anno(list(list("GENCODE.Category" = annotation_predicate("nonempty"))))
  n_empty    <- n_anno(list(list("GENCODE.Category" = annotation_predicate("empty"))))

  # Both counts match the fixture, and together they cover every QC-passing SNV.
  expect_equal(n_nonempty,
               sum(.fixture_truth$keep & nzchar(.fixture_truth$category)))
  expect_equal(n_empty,
               sum(.fixture_truth$keep & !nzchar(.fixture_truth$category)))
  expect_equal(n_nonempty + n_empty, .fixture_truth$n_keep)
  expect_gt(n_nonempty, 0L)
})

test_that("nonempty selects a sparsely-populated field correctly", {
  # GeneHancer is populated for only a handful of fixture variants -- exactly the
  # sparse-membership field the noncoding masks key on, and exactly the test the
  # old %in%-only grammar could not express.
  keep <- .fixture_truth$keep
  gh <- .fixture_truth$genehancer
  n_present <- sum(keep & nzchar(gh))

  expect_equal(n_anno(list(list("GeneHancer" = annotation_predicate("nonempty")))),
               n_present)
  expect_equal(n_anno(list(list("GeneHancer" = annotation_predicate("empty")))),
               .fixture_truth$n_keep - n_present)

  # The field is genuinely sparse in this fixture, so the two predicates are
  # doing real work rather than trivially selecting everything or nothing.
  expect_gt(n_present, 0L)
  expect_lt(n_present, .fixture_truth$n_keep)
})

test_that("numeric predicates select complementary sets and match a manual count", {
  keep <- .fixture_truth$keep
  cadd <- .fixture_truth$cadd
  expect_equal(n_anno(list(list("CADD" = annotation_predicate("gt", 20)))),
               sum(keep & cadd > 20))
  expect_equal(n_anno(list(list("CADD" = annotation_predicate("le", 20)))),
               sum(keep & cadd <= 20))
  expect_equal(n_anno(list(list("CADD" = annotation_predicate("ge", 20)))),
               sum(keep & cadd >= 20))
  expect_equal(n_anno(list(list("CADD" = annotation_predicate("lt", 20)))),
               sum(keep & cadd < 20))
  # gt and le partition the QC-passing SNVs (the fixture's CADD has no NA).
  expect_equal(n_anno(list(list("CADD" = annotation_predicate("gt", 20)))) +
                 n_anno(list(list("CADD" = annotation_predicate("le", 20)))),
               .fixture_truth$n_keep)
  # The threshold has bite in this fixture (neither side is empty).
  expect_gt(sum(keep & cadd > 20), 0L)
  expect_gt(sum(keep & cadd <= 20), 0L)
})

test_that("in and not_in are complementary and in matches the atomic path", {
  n_in  <- n_anno(list(list("GENCODE.Category" = annotation_predicate("in", "upstream"))))
  n_out <- n_anno(list(list("GENCODE.Category" = annotation_predicate("not_in", "upstream"))))
  expect_equal(n_in, .n_cat("upstream"))
  expect_equal(n_in + n_out, .fixture_truth$n_keep)
  # "in" must agree with the plain atomic condition.
  expect_equal(n_in, n_anno(list(list("GENCODE.Category" = "upstream"))))
})

test_that("predicates AND with other terms inside a clause", {
  # upstream AND CADD > 20 must be a subset of upstream.
  n_up <- n_anno(list(list("GENCODE.Category" = "upstream")))
  n_up_hi <- n_anno(list(list("GENCODE.Category" = "upstream",
                              "CADD" = annotation_predicate("gt", 20))))
  expect_lte(n_up_hi, n_up)
  expect_equal(n_up_hi, sum(.fixture_truth$keep &
                              .fixture_truth$category == "upstream" &
                              .fixture_truth$cadd > 20))

  # Mixed clause types OR correctly: (upstream) OR (CADD > 20) is at least
  # as large as either alone.
  n_or <- n_anno(list(list("GENCODE.Category" = "upstream"),
                      list("CADD" = annotation_predicate("gt", 20))))
  expect_gte(n_or, n_up)
  expect_gte(n_or, n_anno(list(list("CADD" = annotation_predicate("gt", 20)))))
})

test_that(".apply_annotation_predicate never returns NA", {
  vals_num <- c(1, NA, 30)
  vals_chr <- c("a", NA, "")
  for (p in list(annotation_predicate("gt", 20), annotation_predicate("le", 20),
                 annotation_predicate("in", 30), annotation_predicate("not_in", 30))) {
    r <- GLOWr:::.apply_annotation_predicate(vals_num, p)
    expect_false(anyNA(r))
    expect_length(r, 3L)
  }
  for (p in list(annotation_predicate("nonempty"), annotation_predicate("empty"))) {
    r <- GLOWr:::.apply_annotation_predicate(vals_chr, p)
    expect_false(anyNA(r))
    expect_length(r, 3L)
  }
  # NA fails a positive numeric test and passes its negation.
  expect_equal(GLOWr:::.apply_annotation_predicate(vals_num, annotation_predicate("gt", 20)),
               c(FALSE, FALSE, TRUE))
  expect_equal(GLOWr:::.apply_annotation_predicate(vals_chr, annotation_predicate("nonempty")),
               c(TRUE, FALSE, FALSE))
  expect_equal(GLOWr:::.apply_annotation_predicate(vals_chr, annotation_predicate("empty")),
               c(FALSE, TRUE, TRUE))
})
