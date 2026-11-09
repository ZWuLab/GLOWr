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

test_that("a neighbor exactly on the window boundary counts, whatever the segment size", {
  # Direct computation of ell_i = sum_j r2tilde_ij 1(|pos_i - pos_j| <= window) over the
  # fixture's biallelic SNVs, with the function's own imputation and standardization rules
  # (a zero-variance column standardizes to zeros, so its correlations are 0).
  g <- SeqArray::seqOpen(test_agds_path, readonly = TRUE)
  na <- SeqArray::seqGetData(g, "$num_allele")
  SeqArray::seqSetFilter(g, variant.sel = (na == 2L), verbose = FALSE)
  snv <- .is_snv_from_alleles(SeqArray::seqGetData(g, "allele"))
  SeqArray::seqSetFilter(g, variant.sel = snv, action = "intersect", verbose = FALSE)
  pos <- SeqArray::seqGetData(g, "position")
  D   <- SeqArray::seqGetData(g, "$dosage")
  SeqArray::seqClose(g)   # the function opens its own handle; gdsfmt refuses a second one
  o <- order(pos); pos <- pos[o]; D <- D[, o, drop = FALSE]
  n <- nrow(D)
  cm <- colMeans(D, na.rm = TRUE)
  nas <- which(is.na(D), arr.ind = TRUE); if (nrow(nas)) D[nas] <- cm[nas[, 2L]]
  csd <- sqrt((colSums(D * D) - n * cm * cm) / (n - 1)); csd[!is.finite(csd) | csd == 0] <- 1
  Ds <- sweep(sweep(D, 2L, cm, "-"), 2L, csd, "/")
  R  <- crossprod(Ds) / (n - 1)
  R2adj <- R * R - (1 - R * R) / (n - 2)
  gaps <- diff(pos); gaps <- gaps[gaps > 0]
  # Windows that hit a pairwise distance exactly, and one between two distances
  windows <- c(gaps[1L], gaps[length(gaps) %/% 2L], max(gaps), sort(unique(gaps))[2L] + 0.5)
  for (W in windows) {
    ref <- rowSums(R2adj * (abs(outer(pos, pos, "-")) <= W))
    for (seg in c(1L, 2L, 3L, 2000L)) {
      ld <- compute_ld_scores(test_agds_path, window = W, segment = seg, verbose = 0)
      expect_equal(ld$pos, pos)
      expect_equal(ld$ld, ref, tolerance = 1e-10,
                   label = sprintf("window %s, segment %d", format(W), seg))
    }
  }
})

test_that("n_cores spreads the segments over workers without changing the table", {
  skip_on_os("windows")
  # Two workers: R CMD check caps a test at two simultaneous processes (_R_CHECK_LIMIT_CORES_).
  a <- compute_ld_scores(test_agds_path, segment = 7L, verbose = 0)
  b <- compute_ld_scores(test_agds_path, segment = 7L, n_cores = 2L, verbose = 0)
  expect_identical(a, b)
  # more workers than segments: capped to the one segment, same table
  d <- compute_ld_scores(test_agds_path, segment = 2000L, n_cores = 2L, verbose = 0)
  expect_identical(d$ld, compute_ld_scores(test_agds_path, segment = 2000L, verbose = 0)$ld)
})

test_that("the file's storage order does not change the table (positions unsorted on disk)", {
  # A GDS written from a VCF whose rows are not in position order keeps that order as its
  # storage order, and a filter by variant.id returns the dosage columns in storage order,
  # not in the order asked for. The function must put them in position order itself
  set.seed(11)
  n_s <- 40L; n_v <- 30L
  pos_sorted <- sort(sample(1000:9000, n_v))
  pos_disk   <- pos_sorted[sample(n_v)]                       # the order on disk
  maf  <- runif(n_v, 0.1, 0.5)
  geno <- vapply(seq_len(n_v), function(v) rbinom(n_s, 2L, maf[v]), integer(n_s))  # n_s x n_v
  gt   <- c("0/0", "0/1", "1/1")
  vcf <- tempfile(fileext = ".vcf"); gds <- tempfile(fileext = ".gds")
  withr::defer(unlink(c(vcf, gds)))
  rows <- vapply(seq_len(n_v), function(v)
    paste(c("22", pos_disk[v], ".", "A", "G", ".", "PASS", ".", "GT", gt[geno[, v] + 1L]),
          collapse = "\t"), character(1))
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=22>",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       sprintf("s%02d", seq_len(n_s))), collapse = "\t"),
               rows), vcf)
  SeqArray::seqVCF2GDS(vcf, gds, verbose = FALSE)
  g <- SeqArray::seqOpen(gds, readonly = TRUE)
  expect_true(is.unsorted(SeqArray::seqGetData(g, "position")))   # the premise of the test
  SeqArray::seqClose(g)
  # the direct computation, in position order, from the simulated genotypes (the alternate
  # allele count; the function reads the reference dosage, and r^2 is invariant to the flip)
  o <- order(pos_disk); posp <- pos_disk[o]; D <- geno[, o, drop = FALSE]
  n <- nrow(D); cm <- colMeans(D)
  csd <- sqrt((colSums(D * D) - n * cm * cm) / (n - 1)); csd[!is.finite(csd) | csd == 0] <- 1
  Ds <- sweep(sweep(D, 2L, cm, "-"), 2L, csd, "/")
  R <- crossprod(Ds) / (n - 1); R2adj <- R * R - (1 - R * R) / (n - 2)
  W <- 2500
  ref <- rowSums(R2adj * (abs(outer(posp, posp, "-")) <= W))
  for (seg in c(1L, 7L, 2000L)) {
    ld <- compute_ld_scores(gds, window = W, segment = seg, verbose = 0)
    expect_equal(ld$pos, posp)
    expect_equal(ld$ld, ref, tolerance = 1e-10, label = sprintf("segment %d", seg))
  }
  if (.Platform$OS.type != "windows") {
    ld2 <- compute_ld_scores(gds, window = W, segment = 7L, n_cores = 2L, verbose = 0)
    expect_equal(ld2$ld, ref, tolerance = 1e-10)
  }
})
