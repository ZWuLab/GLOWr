# ==============================================================================
# Tests for glow_drilldown_df.R
# ==============================================================================
#
# Tests for glow_drilldown_df() (region evidence packet -> flat per-variant df).
# Builds a tiny synthetic evidence packet `e` with 4 post-filter variants:
#   variants 1,2 -> collapsed unit 1; variant 3 -> unit 2 (its own);
#   variant 4 -> LD-removed (ld_kept = FALSE, collapse_group = NA).


make_evidence <- function(with_spa = FALSE, with_mac = TRUE) {
  pf <- data.frame(
    variant_id     = 1:4,
    chr            = rep(1L, 4),
    pos            = c(100L, 200L, 300L, 400L),
    ref            = c("A", "C", "G", "T"),
    alt            = c("T", "G", "C", "A"),
    MAF            = c(0.01, 0.02, 0.10, 0.05),
    MAC            = c(5, 10, 50, 25),
    ld_kept        = c(TRUE, TRUE, TRUE, FALSE),
    collapse_group = c(1L, 1L, 2L, NA_integer_),
    rsid           = c("rs1", "rs2", "rs3", "rs4"),
    stringsAsFactors = FALSE
  )
  # 2 post-collapse units. Unit 1 = variants 1,2 (collapsed); unit 2 = variant 3.
  tu <- list(
    component_post_filter_idx = list(c(1L, 2L), 3L),
    is_collapsed    = c(TRUE, FALSE),
    MAF             = c(0.03, 0.10),
    Z_standard      = c(2.5, -1.2),
    pvalue_standard = c(0.012, 0.23),
    B               = c(0.8, 0.4),
    PI              = c(0.7, 0.3),
    weights = matrix(c(1.0, 2.0, 0.5, 1.5), nrow = 2L, byrow = TRUE,
                     dimnames = list(c("equal", "optimal"), NULL))
  )
  if (with_mac) tu$MAC <- c(15, 50)
  if (with_spa) {
    tu$Z_SPA      <- c(2.4, -1.1)
    tu$pvalue_SPA <- c(0.016, 0.27)
  }
  list(post_filter = pf, test_units = tu)
}


test_that("glow_drilldown_df returns one row per post-filter variant", {
  e <- make_evidence()
  d <- glow_drilldown_df(e)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), nrow(e$post_filter))
})

test_that("glow_drilldown_df has the expected column schema", {
  e <- make_evidence()
  d <- glow_drilldown_df(e)
  base_cols <- c("rsid", "chr", "pos", "ref", "alt", "MAF", "MAC",
                 "ld_kept", "unit_idx", "unit_n_components",
                 "unit_is_collapsed", "unit_MAF", "unit_MAC",
                 "B", "PI", "Z_standard", "pvalue_standard",
                 "Z_SPA", "pvalue_SPA")
  weight_cols <- c("weight_equal", "weight_optimal")
  expect_equal(colnames(d), c(base_cols, weight_cols))
})

test_that("glow_drilldown_df repeats unit-level fields across components", {
  e <- make_evidence()
  d <- glow_drilldown_df(e)
  # variants 1 and 2 share unit 1
  expect_equal(d$unit_idx[1:2], c(1L, 1L))
  expect_equal(d$Z_standard[1], d$Z_standard[2])
  expect_equal(d$Z_standard[1], 2.5)
  expect_equal(d$unit_n_components[1:2], c(2L, 2L))
  expect_true(d$unit_is_collapsed[1])
  # variant 3 is its own unit
  expect_equal(d$unit_idx[3], 2L)
  expect_equal(d$Z_standard[3], -1.2)
  expect_equal(d$unit_n_components[3], 1L)
  expect_false(d$unit_is_collapsed[3])
})

test_that("glow_drilldown_df makes LD-removed variants NA in unit-level cols", {
  e <- make_evidence()
  d <- glow_drilldown_df(e)
  # variant 4 (row 4) was LD-removed
  expect_false(d$ld_kept[4])
  expect_true(is.na(d$unit_idx[4]))
  expect_true(is.na(d$Z_standard[4]))
  expect_true(is.na(d$pvalue_standard[4]))
  expect_true(is.na(d$B[4]))
  expect_true(is.na(d$PI[4]))
  expect_true(is.na(d$unit_MAF[4]))
  expect_true(is.na(d$weight_equal[4]))
  expect_true(is.na(d$weight_optimal[4]))
  # per-variant identity columns are still present (not unit-level)
  expect_equal(d$rsid[4], "rs4")
  expect_equal(d$pos[4], 400L)
})

test_that("glow_drilldown_df maps weights by scheme rowname", {
  e <- make_evidence()
  d <- glow_drilldown_df(e)
  # weights matrix (rows = schemes, cols = units):
  #   equal   row = (1.0, 2.0) for units 1,2
  #   optimal row = (0.5, 1.5) for units 1,2
  expect_equal(d$weight_equal[1:2], c(1.0, 1.0))   # both in unit 1
  expect_equal(d$weight_equal[3], 2.0)             # unit 2
  expect_equal(d$weight_optimal[1], 0.5)           # unit 1
  expect_equal(d$weight_optimal[3], 1.5)           # unit 2
})

test_that("glow_drilldown_df fills Z_SPA/pvalue_SPA with NA when absent", {
  e <- make_evidence(with_spa = FALSE)
  d <- glow_drilldown_df(e)
  expect_true(all(is.na(d$Z_SPA)))
  expect_true(all(is.na(d$pvalue_SPA)))
})

test_that("glow_drilldown_df carries Z_SPA/pvalue_SPA when present", {
  e <- make_evidence(with_spa = TRUE)
  d <- glow_drilldown_df(e)
  expect_equal(d$Z_SPA[1], 2.4)
  expect_equal(d$Z_SPA[3], -1.1)
  expect_true(is.na(d$Z_SPA[4]))  # LD-removed
})

test_that("glow_drilldown_df fills unit_MAC with NA when tu$MAC absent", {
  e <- make_evidence(with_mac = FALSE)
  d <- glow_drilldown_df(e)
  expect_true(all(is.na(d$unit_MAC)))
})
