########## Tests for prepare_B_training_data() Enhancements ##########
#
# This file tests the enhancements to prepare_B_training_data():
# 1. filter_autosomes parameter (TRUE/FALSE)
# 2. Column mapping transparency (user vs auto detection)

context("prepare_B_training_data() - filter_autosomes and column mapping transparency")

# ========== Test 1: filter_autosomes = TRUE (default, current behavior) ==========

test_that("filter_autosomes = TRUE removes non-autosomal chromosomes", {

  # Create test data with autosomes + X + Y + MT + invalid CHR
  test_data <- data.frame(
    SNP = paste0("rs", 1:10),
    CHR = c(1, 2, 3, "X", "Y", "MT", 23, 24, 25, "invalid"),
    POS = 1000 * (1:10),
    MAF = runif(10, 0.01, 0.5),
    P = runif(10, 1e-5, 0.05),
    N = rep(10000, 10),
    BETA = rnorm(10, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Apply filter_autosomes = TRUE (default)
  result <- prepare_B_training_data(
    data = test_data,
    filter_autosomes = TRUE,
    verbose = 0
  )

  # Verify only autosomes (1-22) remain
  expect_true(all(result$data$CHR >= 1 & result$data$CHR <= 22, na.rm = TRUE))

  # Verify X, Y, MT, 23, 24, 25, invalid were removed
  # Original had 10 variants, should keep only first 3 (CHR 1, 2, 3)
  expect_equal(nrow(result$data), 3)

  # Verify CHR column is numeric
  expect_true(is.numeric(result$data$CHR))

  # Verify metadata records filter_autosomes setting
  expect_equal(result$metadata$filter_autosomes, TRUE)
})


# ========== Test 2: filter_autosomes = FALSE (new behavior) ==========

test_that("filter_autosomes = FALSE retains all chromosomes", {

  # Create test data with autosomes + X + Y + MT + invalid CHR
  test_data <- data.frame(
    SNP = paste0("rs", 1:10),
    CHR = c(1, 2, 3, "X", "Y", "MT", 23, 24, 25, "invalid"),
    POS = 1000 * (1:10),
    MAF = runif(10, 0.01, 0.5),
    P = runif(10, 1e-5, 0.05),
    N = rep(10000, 10),
    BETA = rnorm(10, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Apply filter_autosomes = FALSE
  result <- prepare_B_training_data(
    data = test_data,
    filter_autosomes = FALSE,
    verbose = 0
  )

  # Verify all 10 variants are retained
  expect_equal(nrow(result$data), 10)

  # When filter_autosomes = FALSE, CHR column stays unchanged as character
  expect_true(is.character(result$data$CHR))

  # Verify all chromosome values are present (as character strings)
  chr_vals <- result$data$CHR
  expect_true("1" %in% chr_vals)
  expect_true("2" %in% chr_vals)
  expect_true("3" %in% chr_vals)
  expect_true("X" %in% chr_vals)
  expect_true("Y" %in% chr_vals)
  expect_true("MT" %in% chr_vals)

  # Verify metadata records filter_autosomes setting
  expect_equal(result$metadata$filter_autosomes, FALSE)
})


# ========== Test 3: Column mapping transparency (user vs auto) ==========

test_that("column_mappings metadata distinguishes user vs auto", {

  # Create test data with mixed column names
  test_data <- data.frame(
    my_snp_id = paste0("rs", 1:5),    # User will map this
    FREQ = runif(5, 0.01, 0.5),       # Auto-detect as MAF
    pval = runif(5, 1e-5, 0.05),      # User will map this
    sample_size = rep(10000, 5),      # Auto-detect as N
    Beta = rnorm(5, 0, 0.1),          # Auto-detect as BETA
    chr = rep(1, 5),
    pos = 1000 * (1:5),
    stringsAsFactors = FALSE
  )

  # Apply with custom column mapping for some columns
  result <- prepare_B_training_data(
    data = test_data,
    column_mapping = list(
      rsID = "my_snp_id",
      P = "pval"
    ),
    verbose = 0
  )

  # Verify column_mappings metadata exists
  expect_true(!is.null(result$metadata$column_mappings))
  expect_true(!is.null(result$metadata$mapping_sources))

  # Verify user-specified mappings are marked as "user-specified"
  expect_equal(result$metadata$mapping_sources$rsID, "user-specified")
  expect_equal(result$metadata$mapping_sources$P, "user-specified")

  # Verify auto-detected mappings are marked as "auto-detected"
  expect_equal(result$metadata$mapping_sources$MAF, "auto-detected")
  expect_equal(result$metadata$mapping_sources$N, "auto-detected")
  expect_equal(result$metadata$mapping_sources$BETA, "auto-detected")

  # Verify original column names are stored
  expect_equal(result$metadata$column_mappings$rsID, "my_snp_id")
  expect_equal(result$metadata$column_mappings$P, "pval")
  expect_equal(result$metadata$column_mappings$MAF, "FREQ")
  expect_equal(result$metadata$column_mappings$N, "sample_size")
})


# ========== Test 4: Legacy behavior replication ==========

test_that("Can replicate legacy BMD filtering with filter_autosomes = FALSE", {

  # Simulate BMD-like data with various chromosomes
  # In legacy code, all variants were kept regardless of chromosome
  test_data <- data.frame(
    SNP = paste0("rs", 1:20),
    CHR = c(rep(1:5, each=3), rep("X", 2), "Y", "MT", 23),
    POS = 1000 * (1:20),
    MAF = runif(20, 0.01, 0.5),
    P = runif(20, 1e-8, 0.05),
    N = rep(50000, 20),
    BETA = rnorm(20, 0, 0.05),
    STUDY = rep("BMD Study", 20),
    stringsAsFactors = FALSE
  )

  # Apply with filter_autosomes = FALSE (legacy behavior)
  result <- prepare_B_training_data(
    data = test_data,
    trait_type = "continuous",
    filter_autosomes = FALSE,
    verbose = 0
  )

  # Verify all 20 variants are retained
  expect_equal(nrow(result$data), 20)

  # Verify all chromosomes are preserved (CHR stays as character when filter_autosomes = FALSE)
  expect_true(is.character(result$data$CHR))

  # Verify metadata
  expect_equal(result$metadata$filter_autosomes, FALSE)
  expect_equal(result$metadata$n_variants_final, 20)
  expect_equal(result$metadata$n_variants_removed, 0)
})


# ========== Test 5: Backward compatibility ==========

test_that("Default behavior unchanged (filter_autosomes = TRUE)", {

  # Create test data
  test_data <- data.frame(
    SNP = paste0("rs", 1:8),
    CHR = c(1, 2, 3, 4, "X", "Y", "MT", 22),
    POS = 1000 * (1:8),
    MAF = runif(8, 0.01, 0.5),
    P = runif(8, 1e-5, 0.05),
    N = rep(10000, 8),
    BETA = rnorm(8, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Run with default (not specifying filter_autosomes)
  result_default <- prepare_B_training_data(
    data = test_data,
    verbose = 0
  )

  # Run with explicit filter_autosomes = TRUE
  result_explicit <- prepare_B_training_data(
    data = test_data,
    filter_autosomes = TRUE,
    verbose = 0
  )

  # Verify both produce identical results
  expect_equal(nrow(result_default$data), nrow(result_explicit$data))
  expect_equal(nrow(result_default$data), 5)  # 1, 2, 3, 4, 22

  # Verify X, Y, MT were removed in both
  expect_true(all(result_default$data$CHR %in% 1:22))
  expect_true(all(result_explicit$data$CHR %in% 1:22))
})


# ========== Test 6: Integration test with column mapping and filtering ==========

test_that("filter_autosomes works with custom column mappings", {

  # Create test data with custom column names AND various chromosomes
  test_data <- data.frame(
    variant_id = paste0("rs", 1:10),
    chromosome = c(1, 2, 3, "X", "Y", 4, 5, "MT", 22, 23),
    position = 1000 * (1:10),
    allele_freq = runif(10, 0.01, 0.5),
    pvalue = runif(10, 1e-5, 0.05),
    sample_n = rep(10000, 10),
    effect = rnorm(10, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Test with filter_autosomes = TRUE
  result_true <- prepare_B_training_data(
    data = test_data,
    column_mapping = list(
      rsID = "variant_id",
      CHR = "chromosome",
      POS = "position",
      MAF = "allele_freq",
      P = "pvalue",
      N = "sample_n",
      BETA = "effect"
    ),
    filter_autosomes = TRUE,
    verbose = 0
  )

  # Verify only autosomes retained: 1, 2, 3, 4, 5, 22 = 6 variants
  expect_equal(nrow(result_true$data), 6)
  expect_true(all(result_true$data$CHR %in% 1:22))

  # Verify column mappings were tracked
  expect_equal(result_true$metadata$column_mappings$CHR, "chromosome")
  expect_equal(result_true$metadata$mapping_sources$CHR, "user-specified")

  # Test with filter_autosomes = FALSE
  result_false <- prepare_B_training_data(
    data = test_data,
    column_mapping = list(
      rsID = "variant_id",
      CHR = "chromosome",
      POS = "position",
      MAF = "allele_freq",
      P = "pvalue",
      N = "sample_n",
      BETA = "effect"
    ),
    filter_autosomes = FALSE,
    verbose = 0
  )

  # Verify all 10 variants retained
  expect_equal(nrow(result_false$data), 10)

  # Verify metadata
  expect_equal(result_false$metadata$filter_autosomes, FALSE)
  expect_equal(result_true$metadata$filter_autosomes, TRUE)
})


# ========== Test 7: Edge case - no CHR column ==========

test_that("filter_autosomes handles missing CHR column gracefully", {

  # Create test data without CHR column
  test_data <- data.frame(
    SNP = paste0("rs", 1:5),
    MAF = runif(5, 0.01, 0.5),
    P = runif(5, 1e-5, 0.05),
    N = rep(10000, 5),
    BETA = rnorm(5, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Should work without error regardless of filter_autosomes setting
  result_true <- prepare_B_training_data(
    data = test_data,
    filter_autosomes = TRUE,
    verbose = 0
  )

  result_false <- prepare_B_training_data(
    data = test_data,
    filter_autosomes = FALSE,
    verbose = 0
  )

  # Verify both retain all variants (no filtering possible without CHR)
  expect_equal(nrow(result_true$data), 5)
  expect_equal(nrow(result_false$data), 5)
})


# ========== Test 8: Verify message differences ==========

test_that("Verbose messages differ based on filter_autosomes setting", {

  test_data <- data.frame(
    SNP = paste0("rs", 1:5),
    CHR = c(1, 2, 3, "X", "Y"),
    POS = 1000 * (1:5),
    MAF = runif(5, 0.01, 0.5),
    P = runif(5, 1e-5, 0.05),
    N = rep(10000, 5),
    BETA = rnorm(5, 0, 0.1),
    stringsAsFactors = FALSE
  )

  # Capture messages with filter_autosomes = TRUE
  expect_message(
    prepare_B_training_data(test_data, filter_autosomes = TRUE, verbose = 1),
    "filtering to autosomes"
  )

  # Capture messages with filter_autosomes = FALSE
  expect_message(
    prepare_B_training_data(test_data, filter_autosomes = FALSE, verbose = 1),
    "Skipping CHR processing"
  )
})
