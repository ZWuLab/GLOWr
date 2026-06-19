########## Test get_B data loader refactoring ##########
#
# This file tests the refactored data processing functionality:
# 1. PVALUE_MLOG is preserved as P_mlog10 (not converted to P)
# 2. CHR is left unchanged when filter_autosomes = FALSE

test_that("prepare_B_training_data preserves PVALUE_MLOG as P_mlog10", {

  # Create test data with PVALUE_MLOG column AND a P column
  # (validation requires P column to be present)
  test_data <- data.frame(
    SNP = paste0("rs", 1:10),
    CHROM = rep(1, 10),
    BP = 1:10 * 1000,
    RISK.ALLELE.FREQUENCY = runif(10, 0.01, 0.5),
    P = runif(10, 1e-8, 0.05),  # Regular P-values
    PVALUE_MLOG = runif(10, 5, 100),  # -log10(p) values
    INITIAL.SAMPLE.SIZE = rep(10000, 10),
    stringsAsFactors = FALSE
  )

  # Process the data
  result <- prepare_B_training_data(
    data = test_data,
    trait_type = "continuous",
    extract_text_values = FALSE,  # Don't extract text
    verbose = 0
  )

  # Check that P_mlog10 column exists
  expect_true("P_mlog10" %in% names(result$data))

  # Check that P_mlog10 values are preserved (not converted to P)
  # Note: Row order might change due to duplicate removal, so check set equality
  expect_setequal(result$data$P_mlog10, test_data$PVALUE_MLOG)

})


test_that("prepare_B_training_data leaves CHR unchanged when filter_autosomes=FALSE", {

  # Create test data with mixed chromosome values including X, Y, MT
  test_data <- data.frame(
    SNP = paste0("rs", 1:15),
    CHR = c("1", "2", "3", "X", "Y", "MT", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
    POS = 1:15 * 1000,
    MAF = runif(15, 0.01, 0.5),
    P = runif(15, 1e-8, 0.05),
    N = rep(10000, 15),
    stringsAsFactors = FALSE
  )

  # Process with filter_autosomes = FALSE
  result <- prepare_B_training_data(
    data = test_data,
    trait_type = "continuous",
    filter_autosomes = FALSE,  # Keep all chromosomes unchanged
    extract_text_values = FALSE,
    verbose = 0
  )

  # Check that all 15 rows are retained
  expect_equal(nrow(result$data), 15)

  # Check that CHR column retains all original values (as set, order may change)
  expect_setequal(result$data$CHR, test_data$CHR)

  # Check specifically that non-numeric CHR values are present
  expect_true("X" %in% result$data$CHR)
  expect_true("Y" %in% result$data$CHR)
  expect_true("MT" %in% result$data$CHR)

  # Check that CHR is still character type (not converted to numeric)
  expect_true(is.character(result$data$CHR))

  # Check that CHR was NOT converted to numeric
  expect_false(result$metadata$conversions_applied$CHR_to_numeric)

  # Check that filter_autosomes setting is recorded
  expect_false(result$metadata$filter_autosomes)
})


test_that("prepare_B_training_data filters to autosomes when filter_autosomes=TRUE", {

  # Create test data with mixed chromosome values
  test_data <- data.frame(
    SNP = paste0("rs", 1:15),
    CHR = c("1", "2", "3", "X", "Y", "MT", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
    POS = 1:15 * 1000,
    MAF = runif(15, 0.01, 0.5),
    P = runif(15, 1e-8, 0.05),
    N = rep(10000, 15),
    stringsAsFactors = FALSE
  )

  # Process with filter_autosomes = TRUE (default)
  result <- prepare_B_training_data(
    data = test_data,
    trait_type = "continuous",
    filter_autosomes = TRUE,  # Filter to autosomes
    extract_text_values = FALSE,
    verbose = 0
  )

  # Check that only 12 autosomal variants remain (X, Y, MT removed)
  expect_equal(nrow(result$data), 12)

  # Check that CHR is now numeric
  expect_true(is.numeric(result$data$CHR))

  # Check that CHR values are in 1-22 range (here 1-12)
  expect_true(all(result$data$CHR >= 1 & result$data$CHR <= 22))

  # Check that CHR conversion happened
  expect_true(result$metadata$conversions_applied$CHR_to_numeric)

  # Check that filter_autosomes setting is recorded
  expect_true(result$metadata$filter_autosomes)
})


test_that("prepare_B_training_data handles both PVALUE_MLOG and regular P columns", {

  # Create test data with both P and PVALUE_MLOG columns
  test_data <- data.frame(
    SNP = paste0("rs", 1:10),
    CHR = rep("1", 10),
    POS = 1:10 * 1000,
    MAF = runif(10, 0.01, 0.5),
    P = runif(10, 1e-5, 0.05),  # Regular p-values
    PVALUE_MLOG = runif(10, 50, 100),  # -log10(p) for extreme values
    N = rep(10000, 10),
    stringsAsFactors = FALSE
  )

  # Process the data
  result <- prepare_B_training_data(
    data = test_data,
    trait_type = "continuous",
    extract_text_values = FALSE,
    verbose = 0
  )

  # Check that both P and P_mlog10 columns exist
  expect_true("P" %in% names(result$data))
  expect_true("P_mlog10" %in% names(result$data))

  # Check that P values are preserved (as set, order may change due to sorting)
  expect_setequal(result$data$P, test_data$P)

  # Check that P_mlog10 values are preserved (from PVALUE_MLOG)
  expect_setequal(result$data$P_mlog10, test_data$PVALUE_MLOG)

})


test_that("standardize_column_names maps PVALUE_MLOG to P_mlog10", {

  # Create test data with various P-value column names
  test_data_1 <- data.frame(
    SNP = "rs1",
    PVALUE_MLOG = 50,
    stringsAsFactors = FALSE
  )

  test_data_2 <- data.frame(
    SNP = "rs1",
    PVALUE_MLOG10 = 50,
    stringsAsFactors = FALSE
  )

  test_data_3 <- data.frame(
    SNP = "rs1",
    MLOG10P = 50,
    stringsAsFactors = FALSE
  )

  # Test standardization
  result_1 <- standardize_column_names(test_data_1, verbose = 0)
  result_2 <- standardize_column_names(test_data_2, verbose = 0)
  result_3 <- standardize_column_names(test_data_3, verbose = 0)

  # All should have P_mlog10 column
  expect_true("P_mlog10" %in% names(result_1))
  expect_true("P_mlog10" %in% names(result_2))
  expect_true("P_mlog10" %in% names(result_3))

  # All should preserve the value
  expect_equal(result_1$P_mlog10, 50)
  expect_equal(result_2$P_mlog10, 50)
  expect_equal(result_3$P_mlog10, 50)
})


test_that("CHR handling preserves non-autosomal variants when filter_autosomes=FALSE", {

  # Create test data with specific X, Y, MT chromosomes
  test_data <- data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs4", "rs5"),
    CHR = c("X", "Y", "MT", "M", "23"),
    POS = 1:5 * 1000,
    MAF = runif(5, 0.01, 0.5),
    P = runif(5, 1e-8, 0.05),
    N = rep(10000, 5),
    stringsAsFactors = FALSE
  )

  # Process with filter_autosomes = FALSE
  result <- prepare_B_training_data(
    data = test_data,
    trait_type = "continuous",
    filter_autosomes = FALSE,
    extract_text_values = FALSE,
    verbose = 0
  )

  # Check that all 5 variants are retained
  expect_equal(nrow(result$data), 5)

  # Check that CHR values are preserved (as set, order may change)
  expect_setequal(result$data$CHR, test_data$CHR)

  # Check specifically that non-numeric values are retained
  expect_true("X" %in% result$data$CHR)
  expect_true("Y" %in% result$data$CHR)
  expect_true("MT" %in% result$data$CHR)
  expect_true("M" %in% result$data$CHR)

  # Check that CHR is still character type
  expect_true(is.character(result$data$CHR))
})
