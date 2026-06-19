# ==============================================================================
# Unit Tests for marginal_scan.R functions
# ==============================================================================
#
# Tests for:
#   - extract_pheno_covar_gds()
#   - marginal_scan()
#   - annotate_gds_marginal()

library(testthat)

# Test GDS path
test_gds_path <- file.path(
  "..", "..", "..", "..", "data", "local", "large-data", "test",
  "marginal_test_chr22.gds"
)
# Normalize path
test_gds_path <- normalizePath(test_gds_path, mustWork = FALSE)

skip_if_no_test_gds <- function() {
  if (!file.exists(test_gds_path)) {
    skip("Test GDS file not found")
  }
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    skip("SeqArray not installed")
  }
}


# ============================================================================
# Phase 3: extract_pheno_covar_gds() tests
# ============================================================================

test_that("extract_pheno_covar_gds handles PLINK encoding", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, pheno_coding = "plink",
                                    verbose = 0)
  expect_true(all(pheno$Y %in% c(0, 1)))
  expect_equal(pheno$trait, "binary")
  expect_equal(pheno$n_total, 500)
  expect_true(pheno$n_excluded >= 0)
  expect_equal(length(pheno$sample_id), length(pheno$Y))
  # Our test subset has no missing phenotype (all are 1 or 2)
  expect_equal(pheno$n_excluded, 0)
  expect_equal(sum(pheno$Y == 1), 100)  # 100 cases
  expect_equal(sum(pheno$Y == 0), 400)  # 400 controls
})

test_that("extract_pheno_covar_gds reads covariates", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  expect_true(is.matrix(pheno$X))
  expect_equal(ncol(pheno$X), 1)
  expect_true(all(pheno$X[, 1] %in% c(0, 1)))
  expect_equal(colnames(pheno$X), "sex")
})

test_that("extract_pheno_covar_gds excludes missing phenotype", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, verbose = 0)
  expect_equal(length(pheno$Y) + pheno$n_excluded, pheno$n_total)
})

test_that("extract_pheno_covar_gds errors on non-existent field", {
  skip_if_no_test_gds()
  expect_error(
    extract_pheno_covar_gds(test_gds_path, pheno_name = "nonexistent",
                             verbose = 0),
    "not found"
  )
})

test_that("extract_pheno_covar_gds errors on non-existent covariate", {
  skip_if_no_test_gds()
  expect_error(
    extract_pheno_covar_gds(test_gds_path, covar_names = "nonexistent",
                             verbose = 0),
    "not found"
  )
})

test_that("extract_pheno_covar_gds returns NULL X when no covariates", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, verbose = 0)
  expect_null(pheno$X)
})


# ============================================================================
# Phase 4: marginal_scan() standard method tests
# ============================================================================

test_that("marginal_scan produces correct output structure", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  results <- marginal_scan(test_gds_path, nm, use_SPA = FALSE, verbose = 0)

  expect_s3_class(results, "glow_marginal_scan")
  expect_true(all(c("chr", "pos", "ref", "alt", "variant_id", "MAF",
                     "MAC", "score", "var_score", "Z", "pvalue") %in%
                    names(results)))
  expect_true(all(results$MAF >= 0 & results$MAF <= 0.5))
  expect_true(all(results$pvalue >= 0 & results$pvalue <= 1, na.rm = TRUE))
})

test_that("marginal_scan Z-scores match getZ_marg_score", {
  skip_if_no_test_gds()

  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  # Method 1: marginal_scan
  scan_results <- marginal_scan(test_gds_path, nm, use_SPA = FALSE,
                                 mac_cutoff = 0, verbose = 0)

  # Method 2: getZ_marg_score (direct matrix input)
  # Load G from GDS, align samples to match null model order
  gds <- SeqArray::seqOpen(test_gds_path)
  SeqArray::seqSetFilter(gds, sample.id = pheno$sample_id, verbose = FALSE)
  G_raw <- SeqArray::seqGetData(gds, "$dosage")
  gds_ids <- SeqArray::seqGetData(gds, "sample.id")
  SeqArray::seqClose(gds)

  # Reorder G rows to match null model sample order
  reorder <- match(pheno$sample_id, gds_ids)
  G_aligned <- G_raw[reorder, , drop = FALSE]

  # Impute NAs with column means (same as marginal_scan default)
  for (j in seq_len(ncol(G_aligned))) {
    na_idx <- is.na(G_aligned[, j])
    if (any(na_idx)) {
      G_aligned[na_idx, j] <- mean(G_aligned[, j], na.rm = TRUE)
    }
  }

  direct_results <- getZ_marg_score(G_aligned, null_model = nm)

  # Compare Z-scores on non-zero-variance variants
  # (Near-monomorphic variants may differ due to zero-variance thresholding)
  non_zero <- scan_results$var_score > 1e-6
  expect_true(sum(non_zero) > 150)  # Most variants should be non-trivial
  expect_equal(scan_results$Z[non_zero], direct_results$Zscores[non_zero],
               tolerance = 1e-10)
})

test_that("marginal_scan gives same results regardless of chunk_size", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  results_10 <- marginal_scan(test_gds_path, nm, chunk_size = 10,
                               use_SPA = FALSE, mac_cutoff = 0, verbose = 0)
  results_200 <- marginal_scan(test_gds_path, nm, chunk_size = 200,
                                use_SPA = FALSE, mac_cutoff = 0, verbose = 0)

  # Compare on non-near-zero-variance variants
  non_zero <- results_200$var_score > 1e-6
  expect_equal(results_10$Z[non_zero], results_200$Z[non_zero],
               tolerance = 1e-12)
  expect_equal(results_10$pvalue[non_zero], results_200$pvalue[non_zero],
               tolerance = 1e-12)
})

test_that("marginal_scan handles MAC cutoff", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  results_no_cutoff <- marginal_scan(test_gds_path, nm, mac_cutoff = 0,
                                      use_SPA = FALSE, verbose = 0)
  results_cutoff <- marginal_scan(test_gds_path, nm, mac_cutoff = 5,
                                   use_SPA = FALSE, verbose = 0)

  expect_true(nrow(results_cutoff) <= nrow(results_no_cutoff))
  # Variants kept should all have MAC >= 5
  expect_true(all(results_cutoff$MAC >= 5))
})

test_that("marginal_scan writes CSV correctly", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  csv_path <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_path), add = TRUE)

  results <- marginal_scan(test_gds_path, nm, use_SPA = FALSE,
                            output_csv = csv_path, verbose = 0)

  expect_true(file.exists(csv_path))
  csv_data <- read.csv(csv_path)
  expect_equal(nrow(csv_data), nrow(results))
  expect_equal(csv_data$Z, results$Z, tolerance = 1e-10)
})

test_that("marginal_scan errors without sample_id", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  # Fit null model WITHOUT sample_id
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary")

  expect_error(
    marginal_scan(test_gds_path, nm, use_SPA = FALSE, verbose = 0),
    "sample_id must be provided"
  )
})

test_that("marginal_scan handles monomorphic variants correctly", {
  skip_if_no_test_gds()
  # Monomorphic variants should get near-zero variance, small Z
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  results <- marginal_scan(test_gds_path, nm, use_SPA = FALSE,
                            mac_cutoff = 0, verbose = 0)
  # Monomorphic variants (MAC=0) should exist and have near-zero Z
  mono <- results$MAC == 0
  if (any(mono)) {
    expect_true(all(abs(results$Z[mono]) < 0.01))
  }
})

test_that("marginal_scan works with chunk_size = 1", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  # Extreme chunking — one variant at a time
  results_1 <- marginal_scan(test_gds_path, nm, chunk_size = 1,
                              use_SPA = FALSE, mac_cutoff = 0, verbose = 0)
  results_all <- marginal_scan(test_gds_path, nm, chunk_size = 200,
                                use_SPA = FALSE, mac_cutoff = 0, verbose = 0)
  non_zero <- results_all$var_score > 1e-6
  expect_equal(results_1$Z[non_zero], results_all$Z[non_zero],
               tolerance = 1e-12)
})

test_that("marginal_scan missing_imputation='zero' produces valid results", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  results_zero <- marginal_scan(test_gds_path, nm, use_SPA = FALSE,
                                 missing_imputation = "zero",
                                 mac_cutoff = 0, verbose = 0)
  expect_s3_class(results_zero, "glow_marginal_scan")
  expect_true(all(results_zero$pvalue >= 0 & results_zero$pvalue <= 1))
})


# ============================================================================
# Phase 5: SPA integration tests
# ============================================================================

test_that("marginal_scan SPA produces Z_SPA and pvalue_SPA columns", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  results <- marginal_scan(test_gds_path, nm, use_SPA = TRUE, verbose = 0)
  expect_true("Z_SPA" %in% names(results))
  expect_true("pvalue_SPA" %in% names(results))
  # Standard columns should also be present
  expect_true("Z" %in% names(results))
  expect_true("pvalue" %in% names(results))
})

test_that("marginal_scan SPA Z-scores match getZ_marg_score_binary_SPA", {
  skip_if_no_test_gds()

  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  # SPA via marginal_scan
  scan_spa <- marginal_scan(test_gds_path, nm, use_SPA = TRUE,
                             mac_cutoff = 0, verbose = 0)

  # SPA via direct function
  gds <- SeqArray::seqOpen(test_gds_path)
  SeqArray::seqSetFilter(gds, sample.id = pheno$sample_id, verbose = FALSE)
  G_raw <- SeqArray::seqGetData(gds, "$dosage")
  gds_ids <- SeqArray::seqGetData(gds, "sample.id")
  SeqArray::seqClose(gds)

  reorder <- match(pheno$sample_id, gds_ids)
  G_aligned <- G_raw[reorder, , drop = FALSE]
  for (j in seq_len(ncol(G_aligned))) {
    na_idx <- is.na(G_aligned[, j])
    if (any(na_idx)) {
      G_aligned[na_idx, j] <- mean(G_aligned[, j], na.rm = TRUE)
    }
  }

  direct_spa <- getZ_marg_score_binary_SPA(G_aligned, null_model = nm)

  # Compare SPA Z-scores on non-trivial variants
  non_zero <- scan_spa$var_score > 1e-6
  expect_true(sum(non_zero) > 150)
  expect_equal(scan_spa$Z_SPA[non_zero], direct_spa$Zscores[non_zero],
               tolerance = 1e-6)
})

test_that("marginal_scan SPA errors on continuous trait", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  # Create a continuous null model with dummy Y
  nm_cont <- fit_null_model(pheno$X, rnorm(length(pheno$Y)),
                             trait = "continuous",
                             sample_id = pheno$sample_id)

  expect_error(
    marginal_scan(test_gds_path, nm_cont, use_SPA = TRUE, verbose = 0),
    "SPA.*binary"
  )
})

test_that("marginal_scan default SPA=auto uses SPA for binary", {
  skip_if_no_test_gds()
  pheno <- extract_pheno_covar_gds(test_gds_path, covar_names = "sex",
                                    verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)

  # Default use_SPA=NULL -> auto -> TRUE for binary
  results <- marginal_scan(test_gds_path, nm, verbose = 0)
  expect_true("Z_SPA" %in% names(results))
})


# ============================================================================
# Phase 6: annotate_gds_marginal() tests
# ============================================================================

test_that("annotate_gds_marginal writes and reads correctly", {
  skip_if_no_test_gds()

  # Copy test GDS to temp location
  temp_gds <- tempfile(fileext = ".gds")
  file.copy(test_gds_path, temp_gds)
  on.exit(unlink(temp_gds), add = TRUE)

  pheno <- extract_pheno_covar_gds(temp_gds, covar_names = "sex", verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  results <- marginal_scan(temp_gds, nm, use_SPA = FALSE, verbose = 0)

  # Write
  annotate_gds_marginal(temp_gds, results, verbose = 0)

  # Read back
  gds <- gdsfmt::openfn.gds(temp_gds, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)
  node <- gdsfmt::index.gdsn(gds, "annotation/info/MarginalAnalysis")
  mat <- gdsfmt::read.gdsn(node)
  attrs <- gdsfmt::get.attr.gdsn(node)

  expect_equal(ncol(mat), 6)  # MAF, MAC, score, var_score, Z, pvalue
  expect_true("column_names" %in% names(attrs))
  expect_equal(attrs$column_names,
               c("MAF", "MAC", "score", "var_score", "Z", "pvalue"))
  expect_true("analysis_date" %in% names(attrs))
})

test_that("annotate_gds_marginal respects overwrite flag", {
  skip_if_no_test_gds()

  temp_gds <- tempfile(fileext = ".gds")
  file.copy(test_gds_path, temp_gds)
  on.exit(unlink(temp_gds), add = TRUE)

  pheno <- extract_pheno_covar_gds(temp_gds, covar_names = "sex", verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  results <- marginal_scan(temp_gds, nm, use_SPA = FALSE, verbose = 0)

  # Write first time
  annotate_gds_marginal(temp_gds, results, verbose = 0)

  # Write again without overwrite -> error
  expect_error(
    annotate_gds_marginal(temp_gds, results, verbose = 0),
    "already exists"
  )

  # With overwrite=TRUE -> succeeds
  expect_no_error(
    annotate_gds_marginal(temp_gds, results, overwrite = TRUE, verbose = 0)
  )
})

test_that("annotate_gds_marginal includes SPA columns when present", {
  skip_if_no_test_gds()

  temp_gds <- tempfile(fileext = ".gds")
  file.copy(test_gds_path, temp_gds)
  on.exit(unlink(temp_gds), add = TRUE)

  pheno <- extract_pheno_covar_gds(temp_gds, covar_names = "sex", verbose = 0)
  nm <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
                       sample_id = pheno$sample_id)
  results <- marginal_scan(temp_gds, nm, use_SPA = TRUE, verbose = 0)

  annotate_gds_marginal(temp_gds, results, verbose = 0)

  gds <- gdsfmt::openfn.gds(temp_gds, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)
  node <- gdsfmt::index.gdsn(gds, "annotation/info/MarginalAnalysis")
  mat <- gdsfmt::read.gdsn(node)
  attrs <- gdsfmt::get.attr.gdsn(node)

  expect_equal(ncol(mat), 8)  # 6 standard + Z_SPA + pvalue_SPA
  expect_true("Z_SPA" %in% attrs$column_names)
  expect_true("pvalue_SPA" %in% attrs$column_names)
})
