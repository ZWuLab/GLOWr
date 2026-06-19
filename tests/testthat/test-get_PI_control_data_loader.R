# Test file for prepare_PI_control_data() function

# Helper to get test data path
get_test_plink_prefix <- function() {
  file.path(
    testthat::test_path(), "..", "..", "..", "..",
    "data", "ALS_GWAS_sample_data",
    "ALS_GWAS_chr21_hg38_rename_first800"
  )
}

# Check if SeqArray is available
skip_if_no_seqarray <- function() {
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    testthat::skip("SeqArray package not available")
  }
}

# Helper to create test GDS file
create_test_gds <- function() {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  plink_to_gds(plink_prefix, temp_gds, verbose = 0)
  return(temp_gds)
}

# Helper to create test CSV file
create_test_csv <- function() {
  temp_csv <- tempfile(fileext = ".csv")

  # Create simple test data
  test_data <- data.frame(
    CHR = c("21", "21", "21", "21", "21"),
    POS = c(100, 200, 300, 400, 500),
    REF = c("A", "G", "C", "T", "A"),
    ALT = c("G", "T", "A", "C", "G"),
    rsID = paste0("rs", 1:5),
    stringsAsFactors = FALSE
  )

  write.csv(test_data, temp_csv, row.names = FALSE)
  return(temp_csv)
}

# ========== Format Detection Tests ==========

test_that("prepare_PI_control_data auto-detects PLINK format", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Use format="auto"
  result <- prepare_PI_control_data(
    source = plink_prefix,
    format = "auto",
    verbose = 0
  )

  # Verify correct detection
  expect_equal(result$metadata$format, "plink")
  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 800)
})

test_that("prepare_PI_control_data auto-detects GDS format", {
  skip_if_no_seqarray()

  temp_gds <- create_test_gds()
  on.exit(unlink(temp_gds), add = TRUE)

  # Use format="auto" with .gds file
  result <- prepare_PI_control_data(
    source = temp_gds,
    format = "auto",
    verbose = 0
  )

  # Verify correct detection
  expect_equal(result$metadata$format, "gds")
  expect_s3_class(result, "glow_pi_control_data")
})

test_that("prepare_PI_control_data auto-detects CSV format", {
  temp_csv <- create_test_csv()
  on.exit(unlink(temp_csv), add = TRUE)

  result <- prepare_PI_control_data(
    source = temp_csv,
    format = "auto",
    verbose = 0
  )

  expect_equal(result$metadata$format, "csv")
  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 5)
})

test_that("prepare_PI_control_data auto-detects dataframe format", {
  test_df <- data.frame(
    CHR = c("1", "1", "2"),
    POS = c(100, 200, 300),
    REF = c("A", "G", "C"),
    ALT = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(
    source = test_df,
    format = "auto",
    verbose = 0
  )

  expect_equal(result$metadata$format, "dataframe")
  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 3)
})

test_that("format auto-detection fails gracefully for unknown format", {
  expect_error(
    prepare_PI_control_data(
      source = "/nonexistent/file.xyz",
      format = "auto",
      verbose = 0
    ),
    "Cannot auto-detect format"
  )
})

# ========== PLINK Loading Tests ==========

test_that("prepare_PI_control_data loads PLINK correctly", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  result <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    verbose = 0
  )

  # Verify output structure
  expect_s3_class(result, "glow_pi_control_data")
  expect_type(result, "list")
  expect_named(result, c("data", "metadata"))

  # Verify data structure
  expect_s3_class(result$data, "data.frame")
  expect_equal(nrow(result$data), 800)

  # Verify required columns exist
  expect_true(all(c("VarInfo", "CHR", "POS", "REF", "ALT") %in% names(result$data)))

  # Verify VarInfo is first column
  expect_equal(names(result$data)[1], "VarInfo")

  # Verify metadata
  expect_equal(result$metadata$n_input, 800)
  expect_equal(result$metadata$n_output, 800)
  expect_equal(result$metadata$format, "plink")
  expect_equal(result$metadata$strategy, "provided")
})

test_that("PLINK loading generates correct VarInfo format", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  result <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    verbose = 0
  )

  # Check VarInfo format (CHR-POS-REF-ALT)
  varinfo_split <- strsplit(result$data$VarInfo[1], "-")[[1]]
  expect_length(varinfo_split, 4)

  # First VarInfo should be 21-13464068-C-T
  expect_equal(result$data$VarInfo[1], "21-13464068-C-T")

  # Verify VarInfo matches component columns
  expect_equal(
    result$data$VarInfo[1],
    paste(
      result$data$CHR[1],
      result$data$POS[1],
      result$data$REF[1],
      result$data$ALT[1],
      sep = "-"
    )
  )
})

test_that("PLINK loading handles temp vs persistent GDS correctly", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Test 1: Default (temp GDS, should be deleted)
  result1 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    gds_output = NULL,  # Default: temp
    verbose = 0
  )

  # GDS should have been deleted
  if (!is.null(result1$metadata$gds_path)) {
    expect_false(file.exists(result1$metadata$gds_path))
  }
  expect_false(result1$metadata$gds_persistent)

  # Test 2: Persistent GDS
  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  result2 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    gds_output = temp_gds,
    verbose = 0
  )

  # GDS should exist
  expect_true(file.exists(temp_gds))
  expect_equal(result2$metadata$gds_path, temp_gds)
  expect_true(result2$metadata$gds_persistent)
})

# ========== GDS Loading Tests ==========

test_that("prepare_PI_control_data loads GDS correctly", {
  skip_if_no_seqarray()

  temp_gds <- create_test_gds()
  on.exit(unlink(temp_gds), add = TRUE)

  result <- prepare_PI_control_data(
    source = temp_gds,
    format = "gds",
    verbose = 0
  )

  # Verify output
  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 800)
  expect_equal(result$metadata$format, "gds")

  # Verify same variants as PLINK source
  expect_equal(result$data$VarInfo[1], "21-13464068-C-T")
})

test_that("GDS loading handles non-existent file", {
  expect_error(
    prepare_PI_control_data(
      source = "/nonexistent/file.gds",
      format = "gds",
      verbose = 0
    ),
    "GDS file not found"
  )
})

test_that("GDS loading marks GDS as persistent", {
  skip_if_no_seqarray()

  temp_gds <- create_test_gds()
  on.exit(unlink(temp_gds), add = TRUE)

  result <- prepare_PI_control_data(
    source = temp_gds,
    format = "gds",
    verbose = 0
  )

  # User-provided GDS should be marked persistent
  expect_equal(result$metadata$gds_path, temp_gds)
  expect_true(result$metadata$gds_persistent)
})

# ========== CSV Loading Tests ==========

test_that("prepare_PI_control_data loads CSV correctly", {
  temp_csv <- create_test_csv()
  on.exit(unlink(temp_csv), add = TRUE)

  result <- prepare_PI_control_data(
    source = temp_csv,
    format = "csv",
    verbose = 0
  )

  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 5)
  expect_equal(result$metadata$format, "csv")

  # Verify VarInfo generated
  expect_true("VarInfo" %in% names(result$data))
  expect_equal(result$data$VarInfo[1], "21-100-A-G")
})

test_that("CSV loading validates required columns", {
  # Create CSV with missing columns
  temp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_csv), add = TRUE)

  incomplete_data <- data.frame(
    CHR = c("1", "2"),
    POS = c(100, 200),
    # Missing REF and ALT
    stringsAsFactors = FALSE
  )

  write.csv(incomplete_data, temp_csv, row.names = FALSE)

  expect_error(
    prepare_PI_control_data(temp_csv, format = "csv", verbose = 0),
    "CSV file missing required columns"
  )
})

test_that("CSV loading adds rsID if missing", {
  temp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_csv), add = TRUE)

  # CSV without rsID column
  data_no_rsid <- data.frame(
    CHR = c("1", "2"),
    POS = c(100, 200),
    REF = c("A", "G"),
    ALT = c("G", "T"),
    stringsAsFactors = FALSE
  )

  write.csv(data_no_rsid, temp_csv, row.names = FALSE)

  result <- prepare_PI_control_data(temp_csv, format = "csv", verbose = 0)

  expect_true("rsID" %in% names(result$data))
  expect_true(all(is.na(result$data$rsID)))
})

# ========== Data.frame Loading Tests ==========

test_that("prepare_PI_control_data loads data.frame correctly", {
  test_df <- data.frame(
    CHR = c("1", "1", "2"),
    POS = c(100, 200, 300),
    REF = c("A", "G", "C"),
    ALT = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(
    source = test_df,
    format = "dataframe",
    verbose = 0
  )

  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 3)
  expect_equal(result$metadata$format, "dataframe")

  # Verify VarInfo generated
  expect_true("VarInfo" %in% names(result$data))
  expect_equal(result$data$VarInfo[1], "1-100-A-G")
})

test_that("data.frame loading validates required columns", {
  incomplete_df <- data.frame(
    CHR = c("1", "2"),
    POS = c(100, 200)
    # Missing REF and ALT
  )

  expect_error(
    prepare_PI_control_data(incomplete_df, format = "dataframe", verbose = 0),
    "data.frame missing required columns"
  )
})

# ========== Random Sampling Tests ==========

test_that("n_controls samples correct number of variants", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Sample to 100 variants
  result <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = 42,
    verbose = 0
  )

  expect_equal(nrow(result$data), 100)
  expect_equal(result$metadata$n_input, 800)
  expect_equal(result$metadata$n_output, 100)
})

test_that("random_seed produces reproducible sampling", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Sample with same seed twice
  result1 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = 42,
    verbose = 0
  )

  result2 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = 42,
    verbose = 0
  )

  # VarInfo should be identical
  expect_identical(result1$data$VarInfo, result2$data$VarInfo)

  # Verify metadata records seed
  expect_equal(result1$metadata$random_seed, 42)
  expect_equal(result2$metadata$random_seed, 42)
})

test_that("different seeds produce different samples", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Sample with different seeds
  result1 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = 42,
    verbose = 0
  )

  result2 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = 123,
    verbose = 0
  )

  # VarInfo should differ
  expect_false(identical(result1$data$VarInfo, result2$data$VarInfo))
})

test_that("n_controls > n_variants uses all variants with warning", {
  test_df <- data.frame(
    CHR = c("1", "1", "2"),
    POS = c(100, 200, 300),
    REF = c("A", "G", "C"),
    ALT = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  # Request more variants than available
  expect_warning(
    result <- prepare_PI_control_data(
      source = test_df,
      format = "dataframe",
      n_controls = 100,  # More than 3 available
      verbose = 0
    ),
    "Requested n_controls=100 but only 3 variants available"
  )

  # Should use all 3 variants
  expect_equal(nrow(result$data), 3)
})

test_that("random sampling without seed is non-deterministic", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Sample without seed twice
  result1 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = NULL,
    verbose = 0
  )

  result2 <- prepare_PI_control_data(
    source = plink_prefix,
    format = "plink",
    n_controls = 100,
    random_seed = NULL,
    verbose = 0
  )

  # Results should differ (with high probability)
  # Note: This test has small probability of false failure
  expect_false(identical(result1$data$VarInfo, result2$data$VarInfo))
})

# ========== S3 Methods Tests ==========

test_that("print.glow_pi_control_data works", {
  test_df <- data.frame(
    CHR = c("1", "1", "2"),
    POS = c(100, 200, 300),
    REF = c("A", "G", "C"),
    ALT = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  # Capture print output
  output <- capture.output(print(result))
  output_text <- paste(output, collapse = "\n")

  # Should contain key information
  expect_match(output_text, "GLOW PI Control Data")
  expect_match(output_text, "Variants:.*3")
  expect_match(output_text, "Strategy:.*provided")
  expect_match(output_text, "Data format:.*dataframe")
})

test_that("summary.glow_pi_control_data works", {
  test_df <- data.frame(
    CHR = c("1", "1", "2"),
    POS = c(100, 200, 300),
    REF = c("A", "G", "C"),
    ALT = c("G", "T", "A"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  # Capture summary output
  output <- capture.output(summary(result))
  output_text <- paste(output, collapse = "\n")

  # Should show chromosome distribution
  expect_match(output_text, "Chromosome distribution")
  expect_match(output_text, "Data structure")
  expect_match(output_text, "Source format:.*dataframe")
})

# ========== Strategy Tests ==========

test_that("prepare_PI_control_data only accepts 'provided' strategy", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    stringsAsFactors = FALSE
  )

  expect_error(
    prepare_PI_control_data(
      test_df,
      format = "dataframe",
      strategy = "synonymous",
      verbose = 0
    ),
    "Currently only strategy='provided' is supported"
  )
})

test_that("prepare_PI_control_data default strategy is 'provided'", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  expect_equal(result$metadata$strategy, "provided")
})

# ========== Verbosity Tests ==========

test_that("prepare_PI_control_data is silent when verbose=0", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    stringsAsFactors = FALSE
  )

  output <- capture.output(
    result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0),
    type = "message"
  )

  expect_length(output, 0)
})

test_that("prepare_PI_control_data reports progress when verbose=1", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    stringsAsFactors = FALSE
  )

  # Note: capture.output(type = "message") is unreliable inside testthat;
  # use expect_message() instead
  expect_message(
    prepare_PI_control_data(test_df, format = "dataframe", verbose = 1),
    "Preparing PI Control Data"
  )
  expect_message(
    prepare_PI_control_data(test_df, format = "dataframe", verbose = 1),
    "Loading variants from data\\.frame"
  )
  expect_message(
    prepare_PI_control_data(test_df, format = "dataframe", verbose = 1),
    "PI Control Data Preparation Complete"
  )
})

test_that("prepare_PI_control_data shows detailed info when verbose=2", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  output <- capture.output(
    result <- prepare_PI_control_data(
      plink_prefix,
      format = "plink",
      verbose = 2
    ),
    type = "message"
  )

  output_text <- paste(output, collapse = "\n")
  expect_match(output_text, "temporary GDS|Saving GDS")
})

# ========== Column Ordering Tests ==========

test_that("VarInfo is first column in output", {
  test_df <- data.frame(
    ALT = c("G"),
    REF = c("A"),
    POS = c(100),
    CHR = c("1"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  expect_equal(names(result$data)[1], "VarInfo")
})

test_that("columns follow priority order", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    extra_col = c("value"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  # First 6 should be priority columns
  priority_cols <- c("VarInfo", "rsID", "CHR", "POS", "REF", "ALT")
  actual_priority <- names(result$data)[1:6]
  expect_equal(actual_priority, priority_cols)

  # Extra columns should come after
  expect_true("extra_col" %in% names(result$data))
  expect_true(which(names(result$data) == "extra_col") > 6)
})

# ========== Metadata Validation Tests ==========

test_that("metadata contains all required fields", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  # Check metadata fields
  expected_fields <- c(
    "n_input", "n_output", "format", "gds_path",
    "gds_persistent", "strategy", "random_seed", "date_prepared"
  )

  expect_true(all(expected_fields %in% names(result$metadata)))
})

test_that("date_prepared is current date", {
  test_df <- data.frame(
    CHR = c("1"),
    POS = c(100),
    REF = c("A"),
    ALT = c("G"),
    stringsAsFactors = FALSE
  )

  result <- prepare_PI_control_data(test_df, format = "dataframe", verbose = 0)

  expect_equal(result$metadata$date_prepared, Sys.Date())
})

# ========== Integration Tests ==========

test_that("full workflow: PLINK -> GDS -> sample -> output", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Full workflow with all features
  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  result <- prepare_PI_control_data(
    source = plink_prefix,
    format = "auto",
    n_controls = 50,
    gds_output = temp_gds,
    random_seed = 42,
    verbose = 1
  )

  # Verify all aspects
  expect_s3_class(result, "glow_pi_control_data")
  expect_equal(nrow(result$data), 50)
  expect_equal(result$metadata$n_input, 800)
  expect_equal(result$metadata$n_output, 50)
  expect_equal(result$metadata$format, "plink")
  expect_true(file.exists(temp_gds))
  expect_true(result$metadata$gds_persistent)
  expect_equal(result$metadata$random_seed, 42)

  # Verify data quality
  expect_true(all(c("VarInfo", "CHR", "POS", "REF", "ALT") %in% names(result$data)))
  expect_true(all(result$data$CHR == "21" | result$data$CHR == 21))
})
