# ==============================================================================
# Tests for compute_ld_scores()
# ==============================================================================
#
# Structural / contract tests on the synthetic single-chromosome test aGDS
# (helper-create-test-agds.R: 125 variants on chr22). Numerical fidelity to the
# validated research script is checked separately by the Phase-2 reproduction of
# the cached chr22 LD scores (not a unit test; needs the cohort GDS).

skip_if_not_installed("SeqArray")

test_agds_path <- tempfile(fileext = ".gds")
create_test_agds(test_agds_path)
withr::defer(unlink(test_agds_path), teardown_env())

# Expected biallelic-SNV count, using GLOWr's own predicate (the function's filter).
.gds <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
.na  <- SeqArray::seqGetData(.gds, "$num_allele")
.alle <- SeqArray::seqGetData(.gds, "allele")
.all_vid <- SeqArray::seqGetData(.gds, "variant.id")
.expected_snv <- sum(.na == 2L & .is_snv_from_alleles(.alle))
SeqArray::seqClose(.gds)


test_that("compute_ld_scores returns the documented schema over biallelic SNVs", {
  ld <- compute_ld_scores(test_agds_path, verbose = 0)
  expect_s3_class(ld, "data.frame")
  expect_identical(names(ld), c("chr", "pos", "ref", "alt", "variant_id", "ld"))
  expect_equal(nrow(ld), .expected_snv)
  expect_true(all(ld$chr == "22"))                       # single chromosome
  expect_true(all(nchar(ld$ref) == 1L & nchar(ld$alt) == 1L))
})

test_that("LD scores are finite and bounded (sign depends on real LD)", {
  ld <- compute_ld_scores(test_agds_path, verbose = 0)
  expect_true(all(is.finite(ld$ld)))
  # With the LDSC-unbiased estimator r2tilde = r^2 - (1 - r^2)/(n - 2), UNLINKED
  # variants contribute slightly NEGATIVE r2tilde. The synthetic fixture has no
  # real LD, so scores summed over many unlinked neighbours can be negative;
  # that is correct. Positivity on real (LD-bearing) data is validated by the
  # Phase-2 chr22 reproduction. Here we only require finite, bounded values.
  expect_true(all(abs(ld$ld) < nrow(ld)))
})

test_that("output is ordered by position", {
  ld <- compute_ld_scores(test_agds_path, verbose = 0)
  expect_false(is.unsorted(ld$pos))
})

test_that("the segment size is a compute knob that does not change results", {
  a <- compute_ld_scores(test_agds_path, segment = 2000L, verbose = 0)
  b <- compute_ld_scores(test_agds_path, segment = 7L,    verbose = 0)
  expect_equal(a$ld, b$ld, tolerance = 1e-10)
})

test_that("variant_id restricts the variant set", {
  subset_vid <- .all_vid[seq(1L, length(.all_vid), by = 2L)]   # every other id
  ld <- compute_ld_scores(test_agds_path, variant_id = subset_vid, verbose = 0)
  expect_true(all(ld$variant_id %in% subset_vid))
  expect_lte(nrow(ld), .expected_snv)
  expect_gt(nrow(ld), 0L)
})
