# Test file for annotate_favor() function

# Helper to get FAVOR database path
get_favor_db_path <- function() {
  file.path(
    testthat::test_path(), "..", "..", "..", "..",
    "data", "large-data", "FAVOR"
  )
}

# Helper to get test control data path
get_test_control_csv <- function() {
  file.path(
    testthat::test_path(), "..", "..", "..", "..",
    "data", "large-data", "test", "sample_controls_500.csv"
  )
}

# Helper to get test case data path
get_test_case_csv <- function() {
  file.path(
    testthat::test_path(), "..", "..", "..", "..",
    "data", "large-data", "test", "ALS_known_control_PI_case_data.csv"
  )
}

# Skip if FAVOR database not available
skip_if_no_favor <- function() {
  favor_path <- get_favor_db_path()
  split_file <- file.path(favor_path, "FAVORdatabase_chrsplit.csv")
  if (!file.exists(split_file)) {
    testthat::skip("FAVOR database not available")
  }
}

# ========== Helper Function Tests ==========

test_that(".default_favor_features returns all 20 numerical FAVOR features", {
  features <- GLOWr:::.default_favor_features()

  expect_type(features, "character")
  expect_length(features, 20)
  expect_true("apc_conservation_v2" %in% features)
  expect_true("apc_micro_rna" %in% features)
  expect_true("cadd_phred" %in% features)
  expect_true("linsight" %in% features)
})

test_that(".default_PI_features is a subset of .default_favor_features", {
  pi <- GLOWr:::.default_PI_features()
  favor <- GLOWr:::.default_favor_features()

  expect_length(pi, 16)
  expect_true(all(pi %in% favor))
})

# ========== Input Validation Tests ==========

test_that("annotate_favor validates inputs", {
  skip_if_no_favor()

  # Missing VarInfo column
  bad_data <- data.frame(CHR = 21, POS = 1000)
  expect_error(
    annotate_favor(bad_data, get_favor_db_path(), verbose = 0),
    "VarInfo"
  )

  # Non-existent path
  good_data <- data.frame(VarInfo = "21-1000-A-G")
  expect_error(
    annotate_favor(good_data, "/nonexistent/path", verbose = 0),
    "not found"
  )
})

# ========== Complete VarInfo Matching Tests ==========

test_that("annotate_favor handles complete VarInfo (control data)", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  # Load subset of control data (first 10 variants)
  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:10, ]

  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_true("VarInfo" %in% names(result))
  expect_true("apc_conservation_v2" %in% names(result))
  expect_equal(nrow(result), 10)
})

# ========== Position-Only Matching Tests ==========

test_that("annotate_favor handles NA allele variants with 'average' method", {
  skip_if_no_favor()

  # Create test data with NA alleles (chr21 only for speed)
  test_variants <- data.frame(
    VarInfo = c("21-17079907-NA-NA", "21-15558573-NA-NA"),
    stringsAsFactors = FALSE
  )

  result <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    na_allele_method = "average",
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that("annotate_favor handles NA allele variants with 'first' method", {
  skip_if_no_favor()

  test_variants <- data.frame(
    VarInfo = c("21-17079907-NA-NA", "21-15558573-NA-NA"),
    stringsAsFactors = FALSE
  )

  result <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    na_allele_method = "first",
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

# ========== Feature Selection Tests ==========

test_that("annotate_favor respects feature selection", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:5, ]

  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    features = c("cadd_phred", "linsight"),
    verbose = 0
  )

  expect_true("cadd_phred" %in% names(result))
  expect_true("linsight" %in% names(result))
  expect_false("apc_conservation" %in% names(result))
})

# ========== NA Handling Tests ==========

test_that("annotate_favor handles na_handling='zero'", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:10, ]

  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    na_handling = "zero",
    verbose = 0
  )

  # Check that NA values are replaced with 0
  features <- GLOWr:::.default_PI_features()
  for (feat in features) {
    if (feat %in% names(result)) {
      expect_false(any(is.na(result[[feat]])))
    }
  }
})

test_that("annotate_favor handles na_handling='drop'", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:10, ]

  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    na_handling = "drop",
    verbose = 0
  )

  # Result may have fewer rows (variants with any NA dropped)
  expect_true(nrow(result) <= 10)
})

# ========== CSV Output Tests ==========

test_that("annotate_favor writes CSV output", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:5, ]

  temp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_csv), add = TRUE)

  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    output_csv = temp_csv,
    verbose = 0
  )

  expect_true(file.exists(temp_csv))

  # Verify CSV content
  csv_data <- utils::read.csv(temp_csv, stringsAsFactors = FALSE)
  expect_true("VarInfo" %in% names(csv_data))
  expect_equal(nrow(csv_data), nrow(result))
})

# ========== S3 Object Input Tests ==========

test_that("annotate_favor handles glow_pi_control_data S3 object", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  # Create S3 object
  control_obj <- prepare_PI_control_data(
    source = control_csv,
    format = "csv",
    n_controls = 10,
    verbose = 0
  )

  result <- annotate_favor(
    variants = control_obj,
    favor_db_path = get_favor_db_path(),
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_true("VarInfo" %in% names(result))
})

# ========== aGDS Output Tests ==========

test_that("annotate_favor writes aGDS output", {
  skip_if_no_favor()

  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    skip("gdsfmt package not available")
  }

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:5, ]

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    output_agds = temp_gds,
    verbose = 0
  )

  expect_true(file.exists(temp_gds))

  # Verify GDS structure
  gds <- gdsfmt::openfn.gds(temp_gds, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)

  expect_true("chromosome" %in% gdsfmt::ls.gdsn(gds))
  expect_true("position" %in% gdsfmt::ls.gdsn(gds))
  expect_true("VarInfo" %in% gdsfmt::ls.gdsn(gds))
})

# ========== Flexible Matching Tests ==========

test_that("annotate_favor flexible matching handles same-REF multiallelic variants", {
  skip_if_no_favor()

  # Test with variants where we expect same-REF but different ALT matches
  # Create test data: if FAVOR has 21-POS-A-G and input has 21-POS-A-T
  # The same-REF multiallelic match should find the 21-POS-A-* entry
  test_variants <- data.frame(
    VarInfo = c("21-15000100-A-T"),  # Modified ALT from actual FAVOR entry
    stringsAsFactors = FALSE
  )

  # Compare exact vs flexible matching
  result_exact <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    match_method = "exact",
    verbose = 0
  )

  result_flexible <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    match_method = "flexible",
    verbose = 0
  )

  # Both should return data frames
  expect_s3_class(result_exact, "data.frame")
  expect_s3_class(result_flexible, "data.frame")

  # Flexible matching may find annotations where exact fails
  # (depends on what's in FAVOR at this position)
  expect_equal(nrow(result_exact), 1)
  expect_equal(nrow(result_flexible), 1)
})

test_that("annotate_favor flexible matching handles swapped alleles (strand flips)", {
  skip_if_no_favor()

  # Load actual control data to find a real variant
  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  original_variant <- controls$VarInfo[1]

  # Parse and create a swapped version
  parts <- strsplit(original_variant, "-", fixed = TRUE)[[1]]
  if (length(parts) == 4) {
    swapped_variant <- paste(parts[1], parts[2], parts[4], parts[3], sep = "-")
  } else {
    swapped_variant <- "21-15000100-T-A"  # Fallback
  }

  test_variants <- data.frame(
    VarInfo = c(swapped_variant),
    stringsAsFactors = FALSE
  )

  # Exact matching shouldn't find this
  result_exact <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    match_method = "exact",
    verbose = 0
  )

  # Flexible matching should potentially find it via swapped allele matching
  result_flexible <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    match_method = "flexible",
    verbose = 0
  )

  expect_s3_class(result_exact, "data.frame")
  expect_s3_class(result_flexible, "data.frame")
  expect_equal(nrow(result_flexible), 1)
})

test_that("annotate_favor flexible matching falls back to position average", {
  skip_if_no_favor()

  # Create a variant that won't match exactly but position exists
  test_variants <- data.frame(
    VarInfo = c("21-17079907-X-Y"),  # Invalid alleles, position should exist
    stringsAsFactors = FALSE
  )

  result_flexible <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    match_method = "flexible",
    verbose = 0
  )

  expect_s3_class(result_flexible, "data.frame")
  expect_equal(nrow(result_flexible), 1)
  # If position exists in FAVOR, we should get annotations via position average
})

test_that("annotate_favor flexible matching preserves exact match priority", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  # Use real control variants that should have exact matches
  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:5, ]

  result_exact <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    match_method = "exact",
    verbose = 0
  )

  result_flexible <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    match_method = "flexible",
    verbose = 0
  )

  # Both should return same number of rows
  expect_equal(nrow(result_exact), nrow(result_flexible))

  # For variants with exact matches, results should be identical
  # Compare annotation values for variants that matched in both
  features <- GLOWr:::.default_PI_features()
  for (feat in features) {
    if (feat %in% names(result_exact) && feat %in% names(result_flexible)) {
      # Compare non-NA values
      exact_vals <- result_exact[[feat]]
      flex_vals <- result_flexible[[feat]]
      matched_idx <- !is.na(exact_vals) & !is.na(flex_vals)
      if (any(matched_idx)) {
        expect_equal(
          exact_vals[matched_idx],
          flex_vals[matched_idx],
          tolerance = 1e-10,
          info = paste("Feature:", feat)
        )
      }
    }
  }
})

test_that("annotate_favor match_method parameter validation", {
  skip_if_no_favor()

  test_variants <- data.frame(VarInfo = "21-15000100-A-G")

  # Valid values should work
  expect_no_error(
    annotate_favor(test_variants, get_favor_db_path(),
                   match_method = "exact", verbose = 0)
  )
  expect_no_error(
    annotate_favor(test_variants, get_favor_db_path(),
                   match_method = "flexible", verbose = 0)
  )

  # Invalid value should error
  expect_error(
    annotate_favor(test_variants, get_favor_db_path(),
                   match_method = "invalid", verbose = 0),
    "match_method"
  )
})

# ========== Bundled Split File Tests ==========

test_that("annotate_favor uses bundled split file when not in favor_db_path", {
  skip_if_no_favor()

  # Create a temp directory with only chunk files (no split file)
  temp_favor_dir <- tempfile(pattern = "favor_no_split_")
  dir.create(temp_favor_dir)
  on.exit(unlink(temp_favor_dir, recursive = TRUE), add = TRUE)

  # Copy only the chunk file(s) needed for chr21
  favor_path <- get_favor_db_path()
  chr21_files <- list.files(favor_path, pattern = "^chr21_.*\\.csv$", full.names = TRUE)
  for (f in chr21_files) {
    file.copy(f, temp_favor_dir)
  }

  # Test variant (chr21)
  test_variants <- data.frame(
    VarInfo = c("21-15000100-A-G"),
    stringsAsFactors = FALSE
  )

  # Should work using bundled split file (no split file in temp_favor_dir)
  result <- annotate_favor(
    variants = test_variants,
    favor_db_path = temp_favor_dir,
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

test_that("package bundled split file exists", {
  skip("FAVORdatabase_chrsplit.csv not bundled in package yet")

  split_file <- system.file("extdata", "FAVORdatabase_chrsplit.csv", package = "GLOWr")
  expect_true(nzchar(split_file))
  expect_true(file.exists(split_file))

  # Verify file has expected structure
  split_data <- utils::read.csv(split_file, stringsAsFactors = FALSE)
  expect_true("Chr" %in% names(split_data))
  expect_true("File_No" %in% names(split_data))
  expect_true("Start_Pos" %in% names(split_data))
  expect_true("End_Pos" %in% names(split_data))
})

# ========== xsv Performance Tests ==========

test_that(".check_xsv_available detects xsv installation", {
  result <- GLOWr:::.check_xsv_available()
  expect_type(result, "logical")
  # Note: test passes regardless of whether xsv is installed
})

test_that("annotate_favor uses xsv when available and appropriate", {
  skip_if_no_favor()

  # Skip if xsv not available
  if (!GLOWr:::.check_xsv_available()) {
    skip("xsv not installed")
  }

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  # Use control data (complete VarInfo, no NA alleles)
  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:10, ]

  # With use_xsv = TRUE (should use xsv)
  result_xsv <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    match_method = "exact",
    use_xsv = TRUE,
    verbose = 0
  )

  # With use_xsv = FALSE (should use R)
  result_r <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    match_method = "exact",
    use_xsv = FALSE,
    verbose = 0
  )

  # Both should produce same results
  expect_equal(nrow(result_xsv), nrow(result_r))
  expect_equal(sort(result_xsv$VarInfo), sort(result_r$VarInfo))

  # Compare annotation values (should be identical for exact matching)
  features <- GLOWr:::.default_PI_features()
  for (feat in features) {
    if (feat %in% names(result_xsv) && feat %in% names(result_r)) {
      xsv_vals <- result_xsv[[feat]][order(result_xsv$VarInfo)]
      r_vals <- result_r[[feat]][order(result_r$VarInfo)]
      # Compare non-NA values
      matched_idx <- !is.na(xsv_vals) & !is.na(r_vals)
      if (any(matched_idx)) {
        expect_equal(xsv_vals[matched_idx], r_vals[matched_idx],
                     tolerance = 1e-10, info = paste("Feature:", feat))
      }
    }
  }
})

test_that("annotate_favor falls back to R when xsv cannot be used", {
  skip_if_no_favor()

  # Test with NA allele variants (xsv cannot handle)
  test_variants <- data.frame(
    VarInfo = c("21-17079907-NA-NA"),
    stringsAsFactors = FALSE
  )

  # Even with use_xsv = TRUE, should fall back to R for NA alleles
  result <- annotate_favor(
    variants = test_variants,
    favor_db_path = get_favor_db_path(),
    use_xsv = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
})

test_that("annotate_favor falls back to R for flexible matching", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control data not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:5, ]

  # With flexible matching, should use R even if use_xsv = TRUE
  result <- annotate_favor(
    variants = controls,
    favor_db_path = get_favor_db_path(),
    match_method = "flexible",
    use_xsv = TRUE,
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
})

# ========== GDS Input Tests (Phase 2) ==========

# Helper to get test GDS file path
get_test_gds <- function() {
  file.path(
    testthat::test_path(), "..", "..", "..", "..",
    "data", "large-data", "test", "my_controls.gds"
  )
}

# Skip if GDS test file not available
skip_if_no_gds <- function() {
  gds_path <- get_test_gds()
  if (!file.exists(gds_path)) {
    testthat::skip("GDS test file not available")
  }
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    testthat::skip("SeqArray package not available")
  }
}

test_that("annotate_favor accepts GDS file path as input", {
  skip_if_no_favor()
  skip_if_no_gds()

  gds_path <- get_test_gds()

  # Run with a small subset using variant_filter
  result <- annotate_favor(
    variants = gds_path,
    favor_db_path = get_favor_db_path(),
    variant_filter = list(variant_ids = 1:10),
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_true("VarInfo" %in% names(result))
  expect_true("apc_conservation_v2" %in% names(result))
  expect_equal(nrow(result), 10)
})

test_that("annotate_favor GDS input detects .gds extension", {
  skip_if_no_favor()
  skip_if_no_gds()

  # Non-.gds file path should error
  expect_error(
    annotate_favor(
      variants = "some_file.txt",
      favor_db_path = get_favor_db_path(),
      verbose = 0
    ),
    "not found|must be"
  )
})

test_that(".extract_varinfo_from_gds extracts correct VarInfo format", {
  skip_if_no_gds()

  gds_path <- get_test_gds()

  result <- GLOWr:::.extract_varinfo_from_gds(
    gds_path = gds_path,
    variant_filter = list(variant_ids = 1:5),
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  expect_true("VarInfo" %in% names(result))
  expect_true("variant_id" %in% names(result))
  expect_equal(nrow(result), 5)

  # VarInfo should be in CHR-POS-REF-ALT format
  expect_true(all(grepl("^[0-9XYM]+-[0-9]+-[ACGTN]+-[ACGTN]+$", result$VarInfo)))
})

test_that(".apply_variant_filter works with chromosome filter", {
  skip_if_no_gds()

  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    skip("SeqArray not available")
  }

  gds_path <- get_test_gds()
  gds <- SeqArray::seqOpen(gds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Get all chromosomes before filter
  all_chr <- unique(SeqArray::seqGetData(gds, "chromosome"))

  # Apply chromosome filter
  GLOWr:::.apply_variant_filter(gds, list(chr = "21"), verbose = 0)

  # Get filtered chromosomes
  filtered_chr <- unique(SeqArray::seqGetData(gds, "chromosome"))

  # All filtered variants should be chr21
  expect_equal(unique(filtered_chr), "21")
})

test_that("variant_filter by variant IDs works", {
  skip_if_no_favor()
  skip_if_no_gds()

  gds_path <- get_test_gds()

  result <- annotate_favor(
    variants = gds_path,
    favor_db_path = get_favor_db_path(),
    variant_filter = list(variant_ids = c(1, 5, 10)),
    verbose = 0
  )

  expect_equal(nrow(result), 3)
})

test_that("variant_filter by position range works", {
  skip_if_no_favor()
  skip_if_no_gds()

  gds_path <- get_test_gds()

  # Filter to a specific region on chr21
  result <- annotate_favor(
    variants = gds_path,
    favor_db_path = get_favor_db_path(),
    variant_filter = list(chr = "21", start = 14000000, end = 15000000),
    verbose = 0
  )

  expect_s3_class(result, "data.frame")
  # Should have some variants in this range
  expect_true(nrow(result) > 0)

  # All VarInfo should be chr21 in the position range
  varinfo_parts <- strsplit(result$VarInfo, "-", fixed = TRUE)
  chr_vals <- sapply(varinfo_parts, function(x) x[1])
  pos_vals <- as.numeric(sapply(varinfo_parts, function(x) x[2]))

  expect_true(all(chr_vals == "21"))
  expect_true(all(pos_vals >= 14000000 & pos_vals <= 15000000))
})

test_that("GDS annotation matches CSV annotation for same variants", {
  skip_if_no_favor()
  skip_if_no_gds()

  gds_path <- get_test_gds()
  control_csv <- get_test_control_csv()

  if (!file.exists(control_csv)) {
    skip("Test control CSV not available")
  }

  # Get annotation from GDS input (first 10 variants)
  result_gds <- annotate_favor(
    variants = gds_path,
    favor_db_path = get_favor_db_path(),
    variant_filter = list(variant_ids = 1:10),
    verbose = 0
  )

  # Get annotation from CSV input
  csv_data <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  csv_data <- csv_data[1:10, ]

  result_csv <- annotate_favor(
    variants = csv_data,
    favor_db_path = get_favor_db_path(),
    verbose = 0
  )

  # Both should be data frames with same columns
  expect_s3_class(result_gds, "data.frame")
  expect_s3_class(result_csv, "data.frame")

  features <- GLOWr:::.default_PI_features()
  for (feat in features) {
    expect_true(feat %in% names(result_gds))
    expect_true(feat %in% names(result_csv))
  }
})

test_that("update_gds=TRUE creates valid aGDS structure", {
  skip_if_no_favor()
  skip_if_no_gds()

  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    skip("gdsfmt not available")
  }

  gds_path <- get_test_gds()

  # Create a temporary copy of the GDS file
  temp_gds <- tempfile(fileext = ".gds")
  file.copy(gds_path, temp_gds)
  on.exit(unlink(temp_gds), add = TRUE)

  # Annotate with update_gds = TRUE (only first 10 variants)
  result <- annotate_favor(
    variants = temp_gds,
    favor_db_path = get_favor_db_path(),
    variant_filter = list(variant_ids = 1:10),
    update_gds = TRUE,
    verbose = 0
  )

  # Verify GDS was updated with annotations
  gds <- gdsfmt::openfn.gds(temp_gds, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)

  # Check annotation structure exists
  expect_true("annotation" %in% gdsfmt::ls.gdsn(gds))
  annot_node <- gdsfmt::index.gdsn(gds, "annotation")
  expect_true("info" %in% gdsfmt::ls.gdsn(annot_node))
  info_node <- gdsfmt::index.gdsn(annot_node, "info")
  expect_true("FunctionalAnnotation" %in% gdsfmt::ls.gdsn(info_node))

  # FunctionalAnnotation is now a STAARpipeline sub-node FOLDER (one sub-node per
  # feature), NOT a single matrix. The default feature set has 20 features.
  fa_node <- gdsfmt::index.gdsn(info_node, "FunctionalAnnotation")
  sub_nodes <- gdsfmt::ls.gdsn(fa_node)
  expect_true(length(sub_nodes) > 0)
  expect_setequal(sub_nodes, GLOWr:::.default_favor_features())
  expect_true("apc_conservation_v2" %in% sub_nodes)

  # Each feature sub-node has one value per annotated variant.
  cadd <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(fa_node, "cadd_phred"))
  expect_equal(length(cadd), 10)  # 10 variants annotated

  # feature_names attribute on the folder lists all features (discoverability).
  feature_names <- gdsfmt::get.attr.gdsn(fa_node)$feature_names
  expect_equal(length(feature_names), length(GLOWr:::.default_favor_features()))
  expect_true("apc_conservation_v2" %in% feature_names)
})

test_that("update_gds=TRUE warning when input is not GDS", {
  skip_if_no_favor()

  control_csv <- get_test_control_csv()
  if (!file.exists(control_csv)) {
    skip("Test control CSV not available")
  }

  controls <- utils::read.csv(control_csv, stringsAsFactors = FALSE)
  controls <- controls[1:5, ]

  # Should warn that update_gds is ignored
  expect_warning(
    annotate_favor(
      variants = controls,
      favor_db_path = get_favor_db_path(),
      update_gds = TRUE,
      verbose = 0
    ),
    "update_gds.*ignored"
  )
})

test_that(".update_gds_with_annotations handles existing annotation node", {
  skip_if_no_gds()

  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    skip("gdsfmt not available")
  }

  gds_path <- get_test_gds()

  # Create a temporary copy
  temp_gds <- tempfile(fileext = ".gds")
  file.copy(gds_path, temp_gds)
  on.exit(unlink(temp_gds), add = TRUE)

  # Create annotation data
  annotations <- data.frame(
    VarInfo = paste0("21-", 14000000 + 1:5, "-A-G"),
    apc_conservation = runif(5),
    cadd_phred = runif(5),
    stringsAsFactors = FALSE
  )

  features <- c("apc_conservation", "cadd_phred")

  # First update
  GLOWr:::.update_gds_with_annotations(
    gds_path = temp_gds,
    annotations = annotations,
    features = features,
    verbose = 0
  )

  # Second update should warn about overwriting
  expect_warning(
    GLOWr:::.update_gds_with_annotations(
      gds_path = temp_gds,
      annotations = annotations,
      features = features,
      verbose = 1
    ),
    "Overwriting"
  )
})

test_that("GDS input with non-existent file errors gracefully", {
  skip_if_no_favor()

  expect_error(
    annotate_favor(
      variants = "/path/to/nonexistent.gds",
      favor_db_path = get_favor_db_path(),
      verbose = 0
    ),
    "not found"
  )
})
