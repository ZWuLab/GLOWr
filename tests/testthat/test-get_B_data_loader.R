########## Tests for Data Loader Functions ##########
#
# This file contains comprehensive tests for the get_B data loader module:
# - prepare_B_training_data()
# - standardize_column_names()
# - apply_qc_filters()
# - validate_B_training_data()
# - S3 methods (print, summary)

library(testthat)

# ========== Helper Functions ==========

#' Create mock GWAS data for testing
create_mock_gwas <- function(
  n = 100,
  include_rsID = TRUE,
  include_beta = TRUE,
  include_chr_pos = TRUE,
  add_duplicates = FALSE,
  add_na = FALSE
) {
  set.seed(123)

  data <- data.frame(
    MAF = runif(n, 0.001, 0.5),
    P = runif(n, 1e-8, 0.1),
    N = round(runif(n, 1000, 10000)),
    stringsAsFactors = FALSE
  )

  if (include_beta) {
    data$BETA <- rnorm(n, 0, 0.1)
  }

  if (include_rsID) {
    data$rsID <- paste0("rs", 1:n)
  }

  if (include_chr_pos) {
    data$CHR <- sample(1:22, n, replace = TRUE)
    data$POS <- round(runif(n, 1e6, 1e8))
  }

  if (add_duplicates) {
    # Add some duplicates
    data <- rbind(data, data[1:5, ])
  }

  if (add_na) {
    # Add some NAs
    data$MAF[1:3] <- NA
    data$P[4:6] <- NA
    data$N[7:9] <- NA
  }

  return(data)
}


# ========== Test standardize_column_names() ==========

test_that("standardize_column_names works with standard names", {
  data <- create_mock_gwas()

  result <- standardize_column_names(data, verbose = 0)

  # Should keep standard names unchanged
  expect_true("MAF" %in% names(result))
  expect_true("P" %in% names(result))
  expect_true("BETA" %in% names(result))
  expect_true("N" %in% names(result))
  expect_true("rsID" %in% names(result))
})


test_that("standardize_column_names handles common variations", {
  # Test various common column name formats
  test_cases <- list(
    list(
      input = c("RAF", "PVALUE", "EFFECT", "SAMPLE_SIZE"),
      expected = c("MAF", "P", "BETA", "N")
    ),
    list(
      input = c("FREQ", "P_VALUE", "EFFECT_SIZE", "INITIAL.SAMPLE.SIZE"),
      expected = c("MAF", "P", "BETA", "N")
    ),
    list(
      input = c("EAF", "P.VALUE", "EFFECT_SIZE", "N_TOTAL"),
      expected = c("MAF", "P", "BETA", "N")
    )
  )

  for (test_case in test_cases) {
    # Create data with alternative names
    data <- data.frame(
      x1 = runif(10, 0.1, 0.5),
      x2 = runif(10, 0.001, 0.1),
      x3 = rnorm(10, 0, 0.1),
      x4 = rep(1000, 10)
    )
    names(data) <- test_case$input

    result <- standardize_column_names(data, verbose = 0)

    for (expected_name in test_case$expected) {
      expect_true(expected_name %in% names(result),
                 info = paste("Expected", expected_name, "in result"))
    }
  }
})


test_that("standardize_column_names handles custom mapping", {
  data <- data.frame(
    my_freq = runif(10, 0.1, 0.5),
    my_pval = runif(10, 0.001, 0.1),
    my_effect = rnorm(10, 0, 0.1),
    my_n = rep(1000, 10)
  )

  result <- standardize_column_names(
    data,
    column_mapping = list(
      MAF = "my_freq",
      P = "my_pval",
      BETA = "my_effect",
      N = "my_n"
    ),
    verbose = 0
  )

  expect_true("MAF" %in% names(result))
  expect_true("P" %in% names(result))
  expect_true("BETA" %in% names(result))
  expect_true("N" %in% names(result))
})


test_that("prepare_B_training_data automatically swaps MAF > 0.5", {
  # Create data with some MAF > 0.5 (major allele frequencies)
  data <- data.frame(
    MAF = c(0.2, 0.7, 0.3, 0.9),  # Some > 0.5
    P = rep(0.01, 4),
    N = rep(1000, 4)
  )

  result <- prepare_B_training_data(
    data = data,
    trait_type = "binary",
    verbose = 0
  )

  # All MAF should be automatically converted to <= 0.5
  expect_true(all(result$data$MAF <= 0.5))
  expect_equal(result$data$MAF[2], 0.3)  # 1 - 0.7 = 0.3
  expect_equal(result$data$MAF[4], 0.1)  # 1 - 0.9 = 0.1
  expect_equal(result$data$MAF[1], 0.2)  # Unchanged
  expect_equal(result$data$MAF[3], 0.3)  # Unchanged
})


test_that("standardize_column_names handles case-insensitive matching", {
  data <- data.frame(
    maf = runif(10, 0.1, 0.5),
    pvalue = runif(10, 0.001, 0.1),
    beta = rnorm(10, 0, 0.1),
    sample_size = rep(1000, 10)
  )

  result <- standardize_column_names(data, verbose = 0)

  expect_true("MAF" %in% names(result))
  expect_true("P" %in% names(result))
  expect_true("BETA" %in% names(result))
  expect_true("N" %in% names(result))
})


# ========== Test apply_qc_filters() ==========

test_that("apply_qc_filters removes NA values", {
  data <- create_mock_gwas(n = 100, add_na = TRUE)

  result <- apply_qc_filters(
    data,
    filters = list(remove_na = TRUE),
    verbose = 0
  )

  # Should have removed rows with NA
  expect_true(nrow(result$data) < 100)
  expect_true(all(!is.na(result$data$MAF)))
  expect_true(all(!is.na(result$data$P)))
  expect_true(all(!is.na(result$data$N)))
})


test_that("apply_qc_filters applies sample size filter", {
  data <- create_mock_gwas(n = 100)
  data$N <- c(rep(100, 30), rep(1000, 70))  # 30 with N < 500

  result <- apply_qc_filters(
    data,
    filters = list(min_sample_size = 500),
    verbose = 0
  )

  # Should have removed 30 variants
  expect_equal(nrow(result$data), 70)
  expect_true(all(result$data$N >= 500))
})


test_that("apply_qc_filters applies p-value filter", {
  data <- create_mock_gwas(n = 100)
  data$P <- c(rep(0.001, 20), rep(0.1, 80))  # 80 with P > 0.01

  result <- apply_qc_filters(
    data,
    filters = list(max_pvalue = 0.01),
    verbose = 0
  )

  # Should have kept only 20 variants
  expect_equal(nrow(result$data), 20)
  expect_true(all(result$data$P <= 0.01))
})


test_that("apply_qc_filters applies min MAF filter", {
  data <- create_mock_gwas(n = 100)
  data$MAF <- c(
    rep(0.001, 10),  # 10 very rare
    rep(0.05, 90)    # 90 in range
  )

  result <- apply_qc_filters(
    data,
    filters = list(
      min_maf = 0.01
    ),
    verbose = 0
  )

  # Should have removed 10 variants with MAF < 0.01
  expect_equal(nrow(result$data), 90)
  expect_true(all(result$data$MAF >= 0.01))
  expect_true(all(result$data$MAF <= 0.5))
})


test_that("apply_qc_filters removes duplicates", {
  data <- create_mock_gwas(n = 100, add_duplicates = TRUE)

  n_before <- nrow(data)

  result <- apply_qc_filters(
    data,
    filters = list(remove_duplicates = TRUE),
    verbose = 0
  )

  # Should have removed duplicates
  expect_true(nrow(result$data) < n_before)
  expect_equal(length(unique(result$data$rsID)), nrow(result$data))
})


test_that("apply_qc_filters returns QC summary", {
  data <- create_mock_gwas(n = 100, add_na = TRUE, add_duplicates = TRUE)

  result <- apply_qc_filters(
    data,
    filters = list(
      min_sample_size = 5000,
      remove_na = TRUE,
      remove_duplicates = TRUE
    ),
    verbose = 0
  )

  # Check QC summary structure
  expect_true("qc_summary" %in% names(result))
  expect_s3_class(result$qc_summary, "data.frame")
  expect_true("filter" %in% names(result$qc_summary))
  expect_true("n_before" %in% names(result$qc_summary))
  expect_true("n_removed" %in% names(result$qc_summary))
  expect_true("n_after" %in% names(result$qc_summary))
})


test_that("apply_qc_filters handles missing columns gracefully", {
  # Data without rsID column
  data <- data.frame(
    MAF = runif(10, 0.1, 0.5),
    P = runif(10, 0.001, 0.1),
    N = rep(1000, 10)
  )

  # Should not fail even with remove_duplicates=TRUE
  result <- apply_qc_filters(
    data,
    filters = list(remove_duplicates = TRUE),
    verbose = 0
  )

  expect_equal(nrow(result$data), 10)
})


# ========== Test validate_B_training_data() ==========

test_that("validate_B_training_data passes for valid data", {
  data <- create_mock_gwas(n = 100)

  result <- validate_B_training_data(data, verbose = 0)

  expect_true(result$valid)
  expect_equal(length(result$errors), 0)
})


test_that("validate_B_training_data detects missing columns", {
  data <- data.frame(
    MAF = runif(10, 0.1, 0.5),
    P = runif(10, 0.001, 0.1)
    # Missing N column
  )

  result <- validate_B_training_data(data, verbose = 0)

  expect_false(result$valid)
  expect_true(any(grepl("Missing required columns", result$errors)))
})


test_that("validate_B_training_data detects invalid MAF range", {
  data <- create_mock_gwas(n = 10)
  data$MAF[1:3] <- c(-0.1, 0, 1.5)  # Invalid values

  result <- validate_B_training_data(data, verbose = 0)

  expect_false(result$valid)
  expect_true(any(grepl("MAF values must be", result$errors)))
})


test_that("validate_B_training_data detects invalid P-values", {
  data <- create_mock_gwas(n = 10)
  data$P[1:3] <- c(-0.01, 0, 1.5)  # Invalid values

  result <- validate_B_training_data(data, verbose = 0)

  expect_false(result$valid)
  expect_true(any(grepl("P-values must be", result$errors)))
})


test_that("validate_B_training_data detects invalid sample sizes", {
  data <- create_mock_gwas(n = 10)
  data$N[1:3] <- c(-100, 0, -50)  # Invalid values

  result <- validate_B_training_data(data, verbose = 0)

  expect_false(result$valid)
  expect_true(any(grepl("Sample sizes .* must be positive", result$errors)))
})


test_that("validate_B_training_data warns about small sample sizes", {
  data <- create_mock_gwas(n = 10)
  data$N <- rep(50, 10)  # Very small N

  result <- validate_B_training_data(data, verbose = 0)

  expect_true(result$valid)  # Valid but with warnings
  expect_true(any(grepl("N < 100", result$warnings)))
})


test_that("validate_B_training_data warns about duplicates", {
  data <- create_mock_gwas(n = 10, add_duplicates = TRUE)

  result <- validate_B_training_data(data, verbose = 0)

  expect_true(any(grepl("duplicate rsIDs", result$warnings)))
})


test_that("validate_B_training_data requires BETA for beta method", {
  data <- create_mock_gwas(n = 10, include_beta = FALSE)

  result <- validate_B_training_data(data, method = "beta", verbose = 0)

  expect_false(result$valid)
  expect_true(any(grepl("Missing required columns.*BETA", result$errors)))
})


test_that("validate_B_training_data returns summary statistics", {
  data <- create_mock_gwas(n = 100)

  result <- validate_B_training_data(data, verbose = 0)

  expect_true(!is.null(result$summary))
  expect_true(!is.null(result$summary$n_variants))
  expect_true(!is.null(result$summary$maf_range))
  expect_true(!is.null(result$summary$p_range))
  expect_true(!is.null(result$summary$n_range))
})


# ========== Test prepare_B_training_data() ==========

test_that("prepare_B_training_data works with data.frame input", {
  data <- create_mock_gwas(n = 100)

  result <- prepare_B_training_data(
    data = data,
    trait_type = "binary",
    verbose = 0
  )

  # Check result structure
  expect_s3_class(result, "glow_training_data")
  expect_true("data" %in% names(result))
  expect_true("metadata" %in% names(result))

  # Check data
  expect_s3_class(result$data, "data.frame")
  expect_true(nrow(result$data) > 0)

  # Check metadata
  expect_equal(result$metadata$trait_type, "binary")
  expect_true(result$metadata$n_variants_final <= result$metadata$n_variants_original)
})


test_that("prepare_B_training_data auto-detects trait type", {
  # Binary trait (small BETAs that can be negative)
  data <- create_mock_gwas(n = 100)
  data$BETA <- rnorm(100, 0, 0.2)  # Small, can be negative

  result <- prepare_B_training_data(
    data = data,
    trait_type = NULL,  # Auto-detect
    verbose = 0
  )

  expect_true(result$metadata$trait_type %in% c("binary", "continuous"))
})


test_that("prepare_B_training_data applies QC filters correctly", {
  data <- create_mock_gwas(n = 100)
  data$N <- c(rep(100, 30), rep(5000, 70))

  result <- prepare_B_training_data(
    data = data,
    qc_filters = list(min_sample_size = 1000),
    verbose = 0
  )

  # Should have filtered out N < 1000
  expect_true(all(result$data$N >= 1000))
  expect_true(result$metadata$n_variants_removed > 0)
})


test_that("prepare_B_training_data handles custom column mapping", {
  data <- data.frame(
    my_maf = runif(50, 0.01, 0.5),
    my_pval = runif(50, 1e-5, 0.1),
    my_n = rep(5000, 50),
    my_beta = rnorm(50, 0, 0.1)
  )

  result <- prepare_B_training_data(
    data = data,
    column_mapping = list(
      MAF = "my_maf",
      P = "my_pval",
      N = "my_n",
      BETA = "my_beta"
    ),
    verbose = 0
  )

  # Check that columns were mapped
  expect_true("MAF" %in% names(result$data))
  expect_true("P" %in% names(result$data))
  expect_true("N" %in% names(result$data))
  expect_true("BETA" %in% names(result$data))
})


test_that("prepare_B_training_data fails for invalid data", {
  # Data missing required columns (no MAF, no P)
  data <- data.frame(
    BETA = rnorm(10, 0, 0.1),
    N = rep(5000, 10)
  )

  expect_error(
    prepare_B_training_data(data = data, verbose = 0),
    "MAF"
  )
})


test_that("prepare_B_training_data records metadata correctly", {
  data <- create_mock_gwas(n = 100)

  result <- prepare_B_training_data(
    data = data,
    trait_type = "continuous",
    format = "dataframe",
    verbose = 0
  )

  # Check metadata fields
  expect_equal(result$metadata$trait_type, "continuous")
  expect_equal(result$metadata$format, "dataframe")
  expect_equal(result$metadata$n_variants_original, 100)
  expect_true(!is.null(result$metadata$date_prepared))
  expect_true(!is.null(result$metadata$qc_applied))
  expect_true(!is.null(result$metadata$qc_summary))
})


# ========== Test CSV File Reading ==========

test_that("prepare_B_training_data can read CSV files", {
  # Create temporary CSV file
  tmp_file <- tempfile(fileext = ".csv")
  data <- create_mock_gwas(n = 50)

  write.csv(data, tmp_file, row.names = FALSE)

  result <- prepare_B_training_data(
    data = tmp_file,
    trait_type = "binary",
    verbose = 0
  )

  # Clean up
  unlink(tmp_file)

  # Check result
  expect_s3_class(result, "glow_training_data")
  expect_true(nrow(result$data) > 0)
  expect_equal(result$metadata$format, "csv")
})


test_that("prepare_B_training_data can read TSV files", {
  # Create temporary TSV file
  tmp_file <- tempfile(fileext = ".tsv")
  data <- create_mock_gwas(n = 50)

  write.table(data, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- prepare_B_training_data(
    data = tmp_file,
    format = "auto",
    trait_type = "binary",
    verbose = 0
  )

  # Clean up
  unlink(tmp_file)

  # Check result
  expect_s3_class(result, "glow_training_data")
  expect_true(nrow(result$data) > 0)
})


test_that("prepare_B_training_data fails for non-existent file", {
  expect_error(
    prepare_B_training_data(
      data = "nonexistent_file.csv",
      verbose = 0
    ),
    "File not found"
  )
})


# ========== Test S3 Methods ==========

test_that("print.glow_training_data works", {
  data <- create_mock_gwas(n = 100)

  result <- prepare_B_training_data(
    data = data,
    trait_type = "binary",
    verbose = 0
  )

  # Should print without error
  expect_output(print(result), "GLOW Training Data")
  expect_output(print(result), "Trait type")
  expect_output(print(result), "Variants")
})


test_that("summary.glow_training_data works", {
  data <- create_mock_gwas(n = 100)

  result <- prepare_B_training_data(
    data = data,
    trait_type = "continuous",
    verbose = 0
  )

  # Should print summary without error
  expect_output(summary(result), "GLOW Training Data Summary")
  expect_output(summary(result), "QC Summary")
})


# ========== Integration Tests ==========

test_that("prepared data works with get_B", {
  # Create training data
  training_data <- create_mock_gwas(n = 50)

  prepared <- prepare_B_training_data(
    data = training_data,
    trait_type = "binary",
    verbose = 0
  )

  # Create target MAF
  target_MAF <- seq(0.01, 0.4, by = 0.05)

  # Run get_B
  B_estimates <- get_B(
    training_trait = prepared$metadata$trait_type,
    training_MAF = prepared$data$MAF,
    training_BETA = prepared$data$BETA,
    target_trait = "binary",
    target_MAF = target_MAF,
    verbose = 0
  )

  # Check output
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
})


test_that("prepared data works with get_B pvalue method", {
  # Create training data
  training_data <- create_mock_gwas(n = 50)

  prepared <- prepare_B_training_data(
    data = training_data,
    trait_type = "binary",
    verbose = 0
  )

  # Create target parameters
  target_MAF <- seq(0.01, 0.4, by = 0.05)
  target_case_prop <- 0.5

  # Run get_B with p-value method
  B_estimates <- get_B(
    training_trait = prepared$metadata$trait_type,
    training_MAF = prepared$data$MAF,
    training_P = prepared$data$P,
    training_N = prepared$data$N,
    target_trait = "binary",
    target_MAF = target_MAF,
    target_case_prop = target_case_prop,
    method = "pvalue",
    verbose = 0
  )

  # Check output
  expect_type(B_estimates, "double")
  expect_length(B_estimates, length(target_MAF))
  expect_true(all(is.finite(B_estimates)))
})


# ========== Edge Cases ==========

test_that("prepare_B_training_data handles single variant", {
  data <- create_mock_gwas(n = 1)

  result <- prepare_B_training_data(
    data = data,
    trait_type = "binary",
    verbose = 0
  )

  expect_equal(nrow(result$data), 1)
})


test_that("prepare_B_training_data handles data with minimal columns", {
  # Only required columns
  data <- data.frame(
    MAF = runif(20, 0.01, 0.5),
    P = runif(20, 1e-5, 0.1),
    N = rep(5000, 20)
  )

  result <- prepare_B_training_data(
    data = data,
    trait_type = "continuous",
    verbose = 0
  )

  expect_s3_class(result, "glow_training_data")
  expect_equal(nrow(result$data), 20)
})


test_that("prepare_B_training_data handles GWAS Catalog format columns", {
  # GWAS Catalog format
  data <- data.frame(
    RISK.ALLELE.FREQUENCY = runif(30, 0.1, 0.9),  # Some > 0.5
    P.VALUE = runif(30, 1e-8, 0.01),
    INITIAL.SAMPLE.SIZE = rep(50000, 30),
    EFFECT_SIZE = rnorm(30, 0, 0.15),
    CHR_ID = sample(1:22, 30, replace = TRUE),
    CHR_POS = round(runif(30, 1e6, 1e8)),
    SNPS = paste0("rs", 1:30)
  )

  result <- prepare_B_training_data(
    data = data,
    format = "gwas_catalog",
    verbose = 0
  )

  # Check standardization worked
  expect_true("MAF" %in% names(result$data))
  expect_true("P" %in% names(result$data))
  expect_true("N" %in% names(result$data))
  expect_true("BETA" %in% names(result$data))

  # Check MAF harmonization
  expect_true(all(result$data$MAF <= 0.5))
})


test_that("prepare_B_training_data works with very strict QC", {
  data <- create_mock_gwas(n = 1000)

  result <- prepare_B_training_data(
    data = data,
    qc_filters = list(
      min_sample_size = 9000,
      max_pvalue = 0.001,
      min_maf = 0.05
    ),
    verbose = 0
  )

  # Should have heavily filtered the data
  expect_true(result$metadata$n_variants_final < result$metadata$n_variants_original)
  # All MAF should still be <= 0.5 after automatic swapping
  expect_true(all(result$data$MAF <= 0.5))
})
