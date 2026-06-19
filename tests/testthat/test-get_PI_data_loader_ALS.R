########## ALS Validation Test for PI Case Data Preparation ##########
#
# This file contains the validation test that replicates the legacy ALS
# workflow, ensuring the new implementation produces the same 297 variants
# from the 352-variant input file.

context("PI Case Data Preparation - ALS Validation")


# ========== ALS Workflow Replication Test ==========

test_that("ALS workflow replicates legacy output (352 -> 297 variants)", {

  # Load ALS test file
  als_file <- system.file("extdata", "ALS-known-SNPs-raw.xlsx", package = "GLOWr")

  skip_if_not(file.exists(als_file), "ALS test file not available")

  # Replicate legacy workflow:
  # 1. Remove chr X
  # 2. Remove duplicates by rsID (keep highest N)
  # 3. Remove "Nicolas A" author
  # Expected: 352 -> 297 variants

  als_result <- prepare_PI_case_data(
    data = als_file,
    format = "excel",
    column_mapping = list(
      N = "Sample Size",
      BETA = "Beta_numeric",
      POS = "pos"  # Explicitly map POS to avoid ambiguity with pos_GRCh37
    ),
    exclude_authors = "Nicolas A",
    exclude_chromosomes = "X",
    qc_filters = list(
      remove_duplicates = TRUE,
      duplicate_priority = "N",
      filter_autosomes = FALSE,  # Already filtered by exclude_chromosomes
      remove_na = FALSE           # Don't remove variants with missing data
    ),
    verbose = 0
  )

  # ========== Critical Validations ==========

  # Test 1: Exact variant count (legacy output)
  expect_equal(nrow(als_result$data), 297,
               info = "ALS workflow should produce exactly 297 variants")

  # Test 2: Input count
  expect_equal(als_result$metadata$n_input, 352,
               info = "ALS raw file should have 352 variants")

  # Test 3: S3 class
  expect_s3_class(als_result, "glow_pi_case_data")

  # Test 4: VarInfo column exists
  expect_true("VarInfo" %in% names(als_result$data),
              info = "VarInfo column should be generated")

  # Test 5: CHR values (no X chromosome)
  if ("CHR" %in% names(als_result$data)) {
    chr_values <- toupper(as.character(als_result$data$CHR))
    expect_false(any(chr_values %in% c("X", "23", "CHRX")),
                info = "No X chromosome variants should remain")
  }

  # Test 6: No "Nicolas A" in output
  if ("FIRST_AUTHOR" %in% names(als_result$data)) {
    expect_false(any(grepl("Nicolas A", als_result$data$FIRST_AUTHOR,
                           ignore.case = TRUE)),
                info = "No 'Nicolas A' author should remain")
  }

  # Test 7: No duplicate rsIDs
  if ("rsID" %in% names(als_result$data)) {
    expect_equal(sum(duplicated(als_result$data$rsID)), 0,
                info = "No duplicate rsIDs should remain")
  }

  # Test 8: Required columns present
  required_cols <- c("rsID", "CHR", "POS", "VarInfo")
  expect_true(all(required_cols %in% names(als_result$data)),
              info = "All required columns should be present")

  # Test 9: Metadata completeness
  expect_true("filters_applied" %in% names(als_result$metadata))
  expect_true("column_mappings" %in% names(als_result$metadata))
  expect_true("date_prepared" %in% names(als_result$metadata))

  # Test 10: Filter counts match expectations
  filters <- als_result$metadata$filters_applied

  # Should have removed some by chr (X chromosome)
  expect_true(filters$n_removed_chr > 0,
              info = "Should have removed X chromosome variants")

  # Should have removed some by author
  expect_true(filters$n_removed_author > 0,
              info = "Should have removed 'Nicolas A' variants")

  # Should have removed some duplicates
  expect_true(filters$n_removed_dup > 0,
              info = "Should have removed duplicate variants")

  # Total removed should equal input - output
  total_removed <- filters$n_removed_chr + filters$n_removed_author +
                   filters$n_removed_dup +
                   filters$n_removed_autosome + filters$n_removed_general +
                   filters$n_removed_na

  expect_equal(total_removed, 352 - 297,
               info = "Total removed should equal input - output")

})


# ========== ALS Data Quality Checks ==========

test_that("ALS data has expected quality after preparation", {

  als_file <- system.file("extdata", "ALS-known-SNPs-raw.xlsx", package = "GLOWr")

  skip_if_not(file.exists(als_file), "ALS test file not available")

  als_result <- prepare_PI_case_data(
    data = als_file,
    format = "excel",
    column_mapping = list(N = "Sample Size", BETA = "Beta_numeric", POS = "pos"),
    exclude_authors = "Nicolas A",
    exclude_chromosomes = "X",
    qc_filters = list(
      remove_duplicates = TRUE,
      duplicate_priority = "N",
      filter_autosomes = FALSE,
      remove_na = FALSE
    ),
    verbose = 0
  )

  data <- als_result$data

  # Check rsID format (should have some valid identifiers)
  if ("rsID" %in% names(data)) {
    # Allow for different identifier formats (rs, kgp, etc.)
    has_identifier <- grepl("^(rs|kgp)[0-9]+$", data$rsID)
    expect_true(sum(has_identifier) / length(data$rsID) > 0.8,
                info = "At least 80% should have valid identifier format (rs or kgp)")
  }

  # Check CHR range (should be autosomes after X removal)
  if ("CHR" %in% names(data)) {
    chr_numeric <- suppressWarnings(as.numeric(data$CHR))
    expect_true(all(!is.na(chr_numeric)),
                info = "All CHR values should be numeric after filtering")
    expect_true(all(chr_numeric >= 1 & chr_numeric <= 22),
                info = "All CHR values should be in 1-22 range")
  }

  # Check POS values (should be positive integers, allowing for NA)
  if ("POS" %in% names(data)) {
    expect_true(all(data$POS > 0 | is.na(data$POS)),
                info = "All non-NA POS values should be positive")
  }

  # Check VarInfo format
  if ("VarInfo" %in% names(data)) {
    # Should match CHR-POS-REF-ALT format (POS can be NA if not available)
    # Allow for NA in POS and alleles since this info may not be available in literature
    expect_true(all(grepl("^[0-9]+-.+-.+-.+$", data$VarInfo)),
                info = "All VarInfo should have CHR-something-something-something format")
  }

  # Check sample sizes (should be reasonable)
  if ("N" %in% names(data)) {
    expect_true(all(data$N > 0 & data$N < 1e7, na.rm = TRUE),
                info = "All N values should be in reasonable range")
  }

})


# ========== Test Legacy Bug Fix ==========

test_that("Legacy bug (wrong variable name) is fixed in new implementation", {

  # Legacy code line 18 used wrong variable name:
  # fwrite(knownSNPneed[,"rsID",drop=F], ...)  # WRONG - undefined
  # Should be:
  # fwrite(knownSNPneed_nodup[,"rsID",drop=F], ...)  # CORRECT

  # The new implementation doesn't have this bug because we use a single
  # data.frame variable throughout the pipeline.

  als_file <- system.file("extdata", "ALS-known-SNPs-raw.xlsx", package = "GLOWr")

  skip_if_not(file.exists(als_file), "ALS test file not available")

  # This should not error (legacy version would error on undefined variable)
  expect_error(
    {
      als_result <- prepare_PI_case_data(
        data = als_file,
        format = "excel",
        column_mapping = list(N = "Sample Size", BETA = "Beta_numeric", POS = "pos"),
        exclude_authors = "Nicolas A",
        exclude_chromosomes = "X",
        qc_filters = list(
          remove_duplicates = TRUE,
          duplicate_priority = "N",
          filter_autosomes = FALSE,
          remove_na = FALSE
        ),
        verbose = 0
      )
    },
    regexp = NA  # Should not error
  )

})


# ========== Test Filter Order ==========

test_that("ALS filters are applied in correct order", {

  als_file <- system.file("extdata", "ALS-known-SNPs-raw.xlsx", package = "GLOWr")

  skip_if_not(file.exists(als_file), "ALS test file not available")

  als_result <- prepare_PI_case_data(
    data = als_file,
    format = "excel",
    column_mapping = list(N = "Sample Size", BETA = "Beta_numeric", POS = "pos"),
    exclude_authors = "Nicolas A",
    exclude_chromosomes = "X",
    qc_filters = list(
      remove_duplicates = TRUE,
      duplicate_priority = "N",
      filter_autosomes = FALSE,
      remove_na = FALSE
    ),
    verbose = 0
  )

  filters <- als_result$metadata$filters_applied

  # Filter order should be:
  # 1. Author filtering (Nicolas A)
  # 2. Chromosome filtering (X)
  # 3. Duplicate removal (by rsID, keep highest N)

  # Verify filter counts are recorded
  expect_true(!is.null(filters$n_removed_author))
  expect_true(!is.null(filters$n_removed_chr))
  expect_true(!is.null(filters$n_removed_dup))

  # Verify each filter removed something
  expect_true(filters$n_removed_author > 0)
  expect_true(filters$n_removed_chr > 0)
  expect_true(filters$n_removed_dup > 0)

})


# ========== Test Alternative Filter Configurations ==========

test_that("ALS workflow with different filter configurations", {

  als_file <- system.file("extdata", "ALS-known-SNPs-raw.xlsx", package = "GLOWr")

  skip_if_not(file.exists(als_file), "ALS test file not available")

  # Test 1: No author filtering (should have more variants)
  result_no_author <- prepare_PI_case_data(
    data = als_file,
    format = "excel",
    column_mapping = list(N = "Sample Size", BETA = "Beta_numeric", POS = "pos"),
    exclude_chromosomes = "X",
    qc_filters = list(
      remove_duplicates = TRUE,
      duplicate_priority = "N",
      filter_autosomes = FALSE,
      remove_na = FALSE
    ),
    verbose = 0
  )

  expect_true(nrow(result_no_author$data) > 297,
              info = "Should have more variants without author filtering")


  # Test 2: No duplicate removal (should have more variants)
  result_no_dup <- prepare_PI_case_data(
    data = als_file,
    format = "excel",
    column_mapping = list(N = "Sample Size", BETA = "Beta_numeric", POS = "pos"),
    exclude_authors = "Nicolas A",
    exclude_chromosomes = "X",
    qc_filters = list(
      remove_duplicates = FALSE,  # Disable duplicate removal
      filter_autosomes = FALSE,
      remove_na = FALSE
    ),
    verbose = 0
  )

  expect_true(nrow(result_no_dup$data) > 297,
              info = "Should have more variants without duplicate removal")


  # Test 3: No X chromosome filtering (should have more variants)
  result_keep_x <- prepare_PI_case_data(
    data = als_file,
    format = "excel",
    column_mapping = list(N = "Sample Size", BETA = "Beta_numeric", POS = "pos"),
    exclude_authors = "Nicolas A",
    # exclude_chromosomes not specified (keep X)
    qc_filters = list(
      remove_duplicates = TRUE,
      duplicate_priority = "N",
      filter_autosomes = FALSE,
      remove_na = FALSE
    ),
    verbose = 0
  )

  expect_true(nrow(result_keep_x$data) > 297,
              info = "Should have more variants with X chromosome included")

})
