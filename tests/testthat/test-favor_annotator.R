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

test_that(".default_favor_features covers the full FAVOR Essential DB content", {
  features <- GLOWr:::.default_favor_features()

  expect_type(features, "character")
  expect_length(features, 30)
  expect_false(anyDuplicated(features) > 0)

  # 20 numeric scores (17 aPCs + CADD/LINSIGHT/FATHMM-XF)
  expect_true("apc_conservation_v2" %in% features)
  expect_true("apc_micro_rna" %in% features)
  expect_true("cadd_phred" %in% features)
  expect_true("linsight" %in% features)

  # 10 categorical fields. These are the regression guard for the 2026-08-21
  # fix: without them an aGDS built with the default cannot run GLOW's own
  # built-in coding masks, which read GENCODE category/exonic category and
  # MetaSVM.
  for (f in c("genecode_comprehensive_category",
              "genecode_comprehensive_exonic_category",
              "genecode_comprehensive_info",
              "genecode_comprehensive_exonic_info",
              "metasvm_pred", "genehancer", "cage_tc", "cage_promoter",
              "rdhs", "rsid")) {
    expect_true(f %in% features, info = f)
  }
})

test_that("the default feature set covers every field the built-in coding masks read", {
  features <- GLOWr:::.default_favor_features()
  catalog <- GLOWr:::annotation_name_catalog
  masks <- GLOWr:::.coding_annotation_masks()

  mask_fields <- unique(unlist(lapply(masks, function(mask) {
    unlist(lapply(mask, names))
  })))
  # Map catalog alias (e.g. "GENCODE.Category") to the aGDS node name.
  nodes <- vapply(mask_fields, function(f) {
    idx <- which(catalog$name == f)
    if (length(idx)) catalog$dir[idx[1]] else f
  }, character(1))

  expect_true(all(nodes %in% features),
              info = paste("missing:", paste(setdiff(nodes, features), collapse = ", ")))
})

test_that("every catalog name resolves to a node the annotator writes by default", {
  # The three aPC aliases inherited from STAARpipeline-Tutorial pointed at nodes
  # the FAVOR databases carry only under versioned names (plan 11, corrected
  # 2026-09-13); every catalog dir must be a default feature.
  catalog <- GLOWr:::annotation_name_catalog
  features <- GLOWr:::.default_favor_features()
  expect_true(all(catalog$dir %in% features),
              info = paste("not written by default:",
                           paste(setdiff(catalog$dir, features), collapse = ", ")))
  alias <- function(nm) catalog$dir[match(nm, catalog$name)]
  expect_equal(alias("aPC.Protein"), "apc_protein_function_v3")
  expect_equal(alias("aPC.Conservation"), "apc_conservation_v2")
  expect_equal(alias("aPC.LocalDiversity"), "apc_local_nucleotide_diversity_v3")
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

# ==============================================================================
# Self-contained tests of the corrected matching and both database backends
# ==============================================================================
# Matching, corrected design"; gates G1 and G2 of
# synthetic (decision D2: no FAVOR content ships with the package). One small
# chr22 database is written once as FAVOR v1 CSV chunks and once as a FAVOR 2.0
# style Parquet file with nested columns, so the two backends read the same
# content.

# Positions of the synthetic database (all on chr22):
#   1000: REF G, SNV rows G-A, G-C, G-T              (complete position)
#   2000: REF A, SNV rows A-C, A-G, A-T and indel A-AT
#   3000: REF C, one SNV row C-T                     (sparse: unmatched_alt)
#   4000: REF T, rows T-A (all features missing), T-C, T-G
#   5000: REF C, SNV rows C-A, C-G, C-T and deletion CAT-C
#   6000: no rows                                     (uncovered)
#   7000: REF A, the indel row A-AT only              (no SNV row: a position-only
#                                                      key is uncovered, decision D10)
# The rsid column imitates FAVOR: dbSNP's rsID sits on the specific REF-ALT row
# dbSNP knows and the other rows at the position carry none (decision D9).
.syn_rows <- function() {
  d <- data.frame(
    position = c(1000L, 1000L, 1000L, 2000L, 2000L, 2000L, 2000L, 3000L,
                 4000L, 4000L, 4000L, 5000L, 5000L, 5000L, 5000L, 7000L),
    ref_vcf  = c("G", "G", "G", "A", "A", "A", "A", "C", "T", "T", "T",
                 "C", "C", "C", "CAT", "A"),
    alt_vcf  = c("A", "C", "T", "C", "G", "T", "AT", "T", "A", "C", "G",
                 "A", "G", "T", "C", "AT"),
    stringsAsFactors = FALSE)
  # Dyadic values, exactly representable in float32 and in decimal text, so the
  # CSV and Parquet copies hold bitwise-identical numbers.
  d$cadd_phred <- c(1.25, 2.5, 3.75, 4.125, 4.25, 4.375, 9.5, 5.5,
                    NA, 6.125, 6.25, 7.125, 7.25, 7.375, 8.5, 10.5)
  d$linsight   <- c(0.5, 0.25, 0.125, 0.5, 0.5, 0.5, NA, 0.75,
                    NA, 0.25, 0.25, 0.5, 0.5, 0.5, NA, NA)
  d$genecode_comprehensive_category <- c(rep("exonic", 3), rep("intronic", 4), "UTR3",
                                         NA, "intergenic", "intergenic",
                                         rep("exonic", 4), "intronic")
  d$metasvm_pred <- c("D", "T", "", rep("", 13))
  d$rsid <- c("rs1000", "", "rs1003", "rs2001", "", "", "rs2007", "rs3000",
              "", "rs4002", "", "rs5001", "rs5002", "rs5003", "rs5004", "rs7007")
  d$variant_vcf <- paste("22", d$position, d$ref_vcf, d$alt_vcf, sep = "-")
  d
}

# Write the rows as a one-chunk FAVOR v1 CSV database (plus a second, empty-range
# chunk so that chunk selection is exercised) and its split table.
.syn_csv_db <- function(dir) {
  d <- .syn_rows()
  cols <- c("variant_vcf", "chromosome", "position", "ref_vcf", "alt_vcf",
            "cadd_phred", "linsight", "genecode_comprehensive_category", "metasvm_pred",
            "rsid")
  d$chromosome <- 22L
  data.table::fwrite(d[d$position < 4000, cols], file.path(dir, "chr22_1.csv"))
  data.table::fwrite(d[d$position >= 4000, cols], file.path(dir, "chr22_2.csv"))
  data.table::fwrite(
    data.frame(Chr = c(22L, 22L), File_No = 1:2, Start_Pos = c(1L, 3501L),
               End_Pos = c(3500L, 1e7)),
    file.path(dir, "FAVORdatabase_chrsplit.csv"))
  dir
}

# Write the same rows as a FAVOR 2.0 style Parquet file with nested columns,
# the scores stored as float32 (as FAVOR 2.0 stores them) and one row per row
# group, so that a position's rows straddle row groups, some groups hold no
# queried row, and the reader's reopening after eight groups is exercised.
# `drop` names features whose leaf columns are left out of the file.
.syn_parquet_db <- function(dir, drop = character(0), rows_per_group = 1L) {
  d <- .syn_rows()
  # arrow::write_parquet() rejects nested data.frame columns, so the struct
  # columns are built as StructArrays (main.cadd.phred is a struct in a struct).
  cadd <- data.frame(row = seq_len(nrow(d))); cadd$cadd <- data.frame(phred = d$cadd_phred)
  cadd$row <- NULL
  cols <- list(
    position = d$position, ref_vcf = d$ref_vcf, alt_vcf = d$alt_vcf,
    variant_vcf = d$variant_vcf, linsight = d$linsight,
    main    = arrow::StructArray$create(cadd),
    gencode = arrow::StructArray$create(
      data.frame(region_type = d$genecode_comprehensive_category, stringsAsFactors = FALSE)),
    dbnsfp  = arrow::StructArray$create(
      data.frame(metasvm_pred = d$metasvm_pred, stringsAsFactors = FALSE)),
    dbsnp   = arrow::StructArray$create(
      data.frame(rsid = d$rsid, stringsAsFactors = FALSE)))
  types <- list(
    position = arrow::int64(), ref_vcf = arrow::string(), alt_vcf = arrow::string(),
    variant_vcf = arrow::string(), linsight = arrow::float32(),
    main = arrow::struct(cadd = arrow::struct(phred = arrow::float32())),
    gencode = arrow::struct(region_type = arrow::string()),
    dbnsfp = arrow::struct(metasvm_pred = arrow::string()),
    dbsnp = arrow::struct(rsid = arrow::string()))
  omit <- c(linsight = "linsight", cadd_phred = "main", genecode_comprehensive_category = "gencode",
            metasvm_pred = "dbnsfp")[drop]
  cols <- cols[setdiff(names(cols), omit)]; types <- types[names(cols)]
  tb <- do.call(arrow::arrow_table, cols)$cast(do.call(arrow::schema, types))
  arrow::write_parquet(tb, file.path(dir, "chromosome_22.parquet"), chunk_size = rows_per_group)
  dir
}

.syn_feats <- function() {
  c("cadd_phred", "linsight", "genecode_comprehensive_category", "metasvm_pred")
}

# One input per outcome. Expected outcome under flexible and exact matching,
# the database row that must supply the annotation, the input's rsID and the
# expected rsID check (decision D9), and the outcome under
# rsid_policy = "require" (a transformed match checked differs or favor_none is
# withheld as rsid_conflict; exact matches and chip_none checks are kept).
.syn_inputs <- function() {
  data.frame(
    VarInfo = c("22-1000-G-A",      # exact
                "22-1000-A-G",      # swapped
                "22-1000-C-T",      # flipped (reverse strand of G-A); rsID differs
                "22-1000-T-C",      # flipped_swapped; no rsID on the input
                "22-1000-A-C",      # flipped_swapped: read as G-T (design's example)
                "22-1000-C-G",      # palindromic containing r: swapped to G-C, no rsID in FAVOR
                "22-1000-A-T",      # palindromic lacking r: unmatched_ref
                "22-2000-A-AT",     # indel, exact
                "22-5000-CAT-C",    # deletion, exact; no rsID on the input
                "22-5000-C-CAT",    # reciprocal of a deletion: never swapped
                "22-5000-GAT-G",    # indel whose anchor differs from r
                "22-3000-C-A",      # sparse position: unmatched_alt
                "22-4000-T-A",      # matched row whose features (and rsid) are all missing
                "22-6000-A-G",      # uncovered
                "22-2000-NA-NA",    # position-only, averaged over the 3 SNV rows
                "22-6000-NA-NA",    # position-only, uncovered
                "22-7000-NA-NA",    # position-only at an indel-only position: uncovered (D10)
                "22-7000-A-AT",     # the indel row itself: exact
                "22-1000-G-<DEL>",  # symbolic allele
                "22-1000-G-A,C",    # multiallelic record
                "22-1000-G-G",      # identical alleles
                "22-1000-g-a",      # lower case: the lookup key is upper-cased; rsID differs
                "chr22-1000-G-A",   # chr prefix; rsID in upper case
                "22-1000-G-A",      # duplicate of the first row
                "23-100-A-G"),      # chromosome without a database
    rsID = c("rs1000", "rs1000", "rs9999", NA, "rs1003", "rs1000", "rs1000", "rs2007",
             "", NA, NA, "rs3000", "rs4000", "rs6000", "rs2000", NA, NA, "rs7007",
             NA, NA, NA, "rs7777", "RS1000", "rs1000", "rs23"),
    flexible = c("exact", "swapped", "flipped", "flipped_swapped", "flipped_swapped",
                 "swapped", "unmatched_ref", "exact", "exact", "unmatched_alt",
                 "unmatched_ref", "unmatched_alt", "exact", "uncovered", "position_only",
                 "uncovered", "uncovered", "exact", "unsupported_allele", "unsupported_allele",
                 "unsupported_allele", "exact", "exact", "exact", "no_database"),
    exact = c("exact", "unmatched_ref", "unmatched_ref", "unmatched_ref", "unmatched_ref",
              "unmatched_ref", "unmatched_ref", "exact", "exact", "unmatched_alt",
              "unmatched_ref", "unmatched_alt", "exact", "uncovered", "position_only",
              "uncovered", "uncovered", "exact", "unsupported_allele", "unsupported_allele",
              "unsupported_allele", "exact", "exact", "exact", "no_database"),
    require = c("exact", "swapped", "rsid_conflict", "flipped_swapped", "flipped_swapped",
                "rsid_conflict", "unmatched_ref", "exact", "exact", "unmatched_alt",
                "unmatched_ref", "unmatched_alt", "exact", "uncovered", "position_only",
                "uncovered", "uncovered", "exact", "unsupported_allele", "unsupported_allele",
                "unsupported_allele", "exact", "exact", "exact", "no_database"),
    key_flexible = c("22-1000-G-A", "22-1000-G-A", "22-1000-G-A", "22-1000-G-A",
                     "22-1000-G-T", "22-1000-G-C", NA, "22-2000-A-AT", "22-5000-CAT-C",
                     NA, NA, NA, "22-4000-T-A", NA, NA, NA, NA, "22-7000-A-AT", NA, NA, NA,
                     "22-1000-G-A", "22-1000-G-A", "22-1000-G-A", NA),
    rsid_check = c("same", "same", "differs", "chip_none", "same", "favor_none", NA,
                   "same", "chip_none", NA, NA, NA, "favor_none", NA, NA, NA, NA, "same",
                   NA, NA, NA, "differs", "same", "same", NA),
    stringsAsFactors = FALSE)
}

.syn_annotate <- function(db, method, format, rsid_policy = "record", ...) {
  # "record" by default here, so that the outcome tables of the fixture (which
  # describe the lookup alone) apply; the package default is "require".
  inp <- .syn_inputs()
  suppressWarnings(annotate_favor(
    variants = inp[, c("VarInfo", "rsID")], favor_db_path = db,
    features = .syn_feats(), match_method = method, favor_db_format = format,
    verbose = 0, rsid_policy = rsid_policy, ...))
}

test_that("the 48 SNV combinations resolve as the design's table states", {
  bases <- c("A", "C", "G", "T")
  fr <- do.call(rbind, lapply(seq_along(bases), function(i)
    data.frame(position = i, ref_vcf = bases[i], alt_vcf = setdiff(bases, bases[i]))))
  cases <- expand.grid(r = seq_along(bases), a1 = bases, a2 = bases,
                       stringsAsFactors = FALSE)
  cases <- cases[cases$a1 != cases$a2, ]
  pal <- GLOWr:::.favor_complement(cases$a1) == cases$a2
  # The design's rule, written out independently of the classifier: the tier
  # each case must get, and the REF-ALT of the row it must point at.
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  rule <- function(r, a1, a2, method) {
    if (a1 == r) return(c("exact", a1, a2))
    if (method == "exact") return(c("unmatched_ref", NA, NA))
    if (a2 == r) return(c("swapped", a2, a1))
    if (comp[[a1]] == r) return(c("flipped", comp[[a1]], comp[[a2]]))
    if (comp[[a2]] == r) return(c("flipped_swapped", comp[[a2]], comp[[a1]]))
    c("unmatched_ref", NA, NA)
  }
  for (method in c("flexible", "exact")) {
    cl <- GLOWr:::.favor_classify(cases$r, cases$a1, cases$a2, rep(TRUE, nrow(cases)),
                                  fr$position, fr$ref_vcf, fr$alt_vcf, method)
    want <- t(mapply(rule, bases[cases$r], cases$a1, cases$a2, MoreArgs = list(method = method)))
    # Case by case: the tier, and the row's REF and ALT.
    expect_equal(cl$tier, unname(want[, 1]), info = method)
    got_ref <- ifelse(is.na(cl$row), NA, fr$ref_vcf[cl$row])
    got_alt <- ifelse(is.na(cl$row), NA, fr$alt_vcf[cl$row])
    expect_equal(got_ref, unname(want[, 2]), info = method)
    expect_equal(got_alt, unname(want[, 3]), info = method)
    # The totals of the design's table.
    tab <- table(pal, factor(cl$tier, GLOWr:::.favor_tier_levels()))
    if (method == "flexible") {
      expect_equal(as.integer(tab["FALSE", c("exact", "swapped", "flipped", "flipped_swapped")]),
                   c(8L, 8L, 8L, 8L))
      expect_equal(as.integer(tab["TRUE", c("exact", "swapped", "unmatched_ref")]),
                   c(4L, 4L, 8L))
      expect_equal(sum(tab[, "unmatched_alt"]), 0L)   # unreachable at a complete position
    } else {
      expect_equal(sum(tab[, "exact"]), 12L)
      expect_equal(sum(tab[, "unmatched_ref"]), 36L)
    }
  }
})

test_that("at a sparse position a transformation that finds the reference but no row is unmatched_alt", {
  # Only A-G exists. Input T-G: REF as given is not the reference base A, the
  # complement A-C is the reference orientation but no row carries A-C. The
  # outcome names the alternate allele as the reason (unmatched_alt); under
  # exact matching no transformation is enabled, so it is unmatched_ref.
  cl <- GLOWr:::.favor_classify(1L, "T", "G", TRUE, 1L, "A", "G", "flexible")
  expect_equal(cl$tier, "unmatched_alt")
  cl <- GLOWr:::.favor_classify(1L, "T", "G", TRUE, 1L, "A", "G", "exact")
  expect_equal(cl$tier, "unmatched_ref")
  # C-G (palindromic, containing neither A nor its complement T): unmatched_ref.
  cl <- GLOWr:::.favor_classify(1L, "C", "G", TRUE, 1L, "A", "G", "flexible")
  expect_equal(cl$tier, "unmatched_ref")
})

test_that("a malformed position is refused rather than truncated", {
  db <- .syn_csv_db(withr::local_tempdir())
  for (bad in c("22-101.9-A-C", "22-0-A-C", "22-1e3-A-C", "22-99999999999-A-C")) {
    expect_error(annotate_favor(data.frame(VarInfo = c("22-1000-G-A", bad)), favor_db_path = db,
                                features = .syn_feats(), use_xsv = FALSE, verbose = 0),
                 "malformed position", info = bad)
  }
  # A key with the wrong number of parts is still an unsupported input, not an error.
  res <- annotate_favor(data.frame(VarInfo = c("22-1000-G-A", "22-1000-G")), favor_db_path = db,
                        features = .syn_feats(), use_xsv = FALSE, verbose = 0)
  expect_equal(res$match_tier, c("exact", "unsupported_allele"))
})

test_that("every outcome is recorded with its database key (CSV backend)", {
  db <- .syn_csv_db(withr::local_tempdir())
  inp <- .syn_inputs()
  for (method in c("flexible", "exact")) {
    res <- .syn_annotate(db, method, "csv", use_xsv = FALSE)
    expect_equal(res$VarInfo, inp$VarInfo)                      # identity and order kept
    expect_equal(res$match_tier, inp[[method]], info = method)
  }
  res <- .syn_annotate(db, "flexible", "csv", use_xsv = FALSE)
  expect_equal(res$favor_key, inp$key_flexible)
  expect_equal(res$rsid_check, inp$rsid_check)
  rows <- .syn_rows()
  # A flipped variant stores its own row's value, not the mean of the position.
  i <- which(inp$VarInfo == "22-1000-C-T")
  expect_equal(res$cadd_phred[i], 1.25)
  expect_equal(res$genecode_comprehensive_category[i], "exonic")
  expect_equal(res$metasvm_pred[i], "D")
  # A matched row whose features are all missing stays matched, not re-averaged.
  j <- which(inp$VarInfo == "22-4000-T-A")
  expect_equal(res$match_tier[j], "exact")
  expect_true(is.na(res$cadd_phred[j]) && is.na(res$linsight[j]))
  # Position-only average: numeric mean over the three SNV rows only (the indel
  # row A-AT, cadd 9.5, is excluded; decision D10), string = first non-empty.
  k <- which(inp$VarInfo == "22-2000-NA-NA")
  expect_equal(res$cadd_phred[k], mean(c(4.125, 4.25, 4.375)))
  expect_equal(res$linsight[k], 0.5)
  expect_equal(res$genecode_comprehensive_category[k], "intronic")
  # A position covered only by an indel row is uncovered for a position-only key,
  # while the indel itself matches by its exact key.
  expect_equal(res$match_tier[inp$VarInfo == "22-7000-NA-NA"], "uncovered")
  expect_equal(res$match_tier[inp$VarInfo == "22-7000-A-AT"], "exact")
  expect_equal(res$cadd_phred[inp$VarInfo == "22-7000-A-AT"], 10.5)
  # Duplicates each receive the annotation; non-matches carry NA.
  expect_equal(res$cadd_phred[inp$VarInfo == "22-1000-G-A"], c(1.25, 1.25))
  expect_true(all(is.na(res$cadd_phred[res$match_tier %in%
    c("unmatched_alt", "unmatched_ref", "uncovered", "no_database", "unsupported_allele")])))
  # Features keep their type even when every value is missing.
  expect_type(res$metasvm_pred, "character")
  expect_type(res$cadd_phred, "double")
})

test_that("the rsID check records the evidence, and rsid_policy = 'require' withholds unconfirmed transformed matches", {
  db <- .syn_csv_db(withr::local_tempdir())
  inp <- .syn_inputs()
  rec <- .syn_annotate(db, "flexible", "csv", use_xsv = FALSE)
  req <- .syn_annotate(db, "flexible", "csv", use_xsv = FALSE, rsid_policy = "require")
  # "record" keeps every match and records the check; "require" withholds the
  # transformed matches checked differs or favor_none, and nothing else.
  expect_equal(rec$match_tier, inp$flexible)
  expect_equal(req$match_tier, inp$require)
  expect_equal(req$rsid_check, inp$rsid_check)      # the check itself does not change
  w <- which(req$match_tier == "rsid_conflict")
  expect_equal(inp$VarInfo[w], c("22-1000-C-T", "22-1000-C-G"))
  expect_equal(req$rsid_check[w], c("differs", "favor_none"))
  # A withheld row carries no annotation, and favor_key names the withheld row.
  for (f in .syn_feats()) expect_true(all(is.na(req[[f]][w])), info = f)
  expect_equal(req$favor_key[w], rec$favor_key[w])
  # Exact matches are never withheld, whatever the check; chip_none passes through.
  e <- which(inp$VarInfo == "22-1000-g-a")
  expect_equal(req$match_tier[e], "exact"); expect_equal(req$rsid_check[e], "differs")
  e <- which(inp$VarInfo == "22-1000-T-C")
  expect_equal(req$match_tier[e], "flipped_swapped"); expect_equal(req$rsid_check[e], "chip_none")
  # Every kept value equals the "record" run's.
  kept <- req$match_tier != "rsid_conflict"
  for (f in .syn_feats()) expect_equal(req[[f]][kept], rec[[f]][kept], info = f)
  # Non-matches and position-only keys carry no check.
  expect_true(all(is.na(rec$rsid_check[rec$match_tier %in% c("position_only", "uncovered",
    "unmatched_alt", "unmatched_ref", "no_database", "unsupported_allele")])))
  # Without an rsID column every check is chip_none and "require" withholds nothing.
  none <- suppressWarnings(annotate_favor(inp[, "VarInfo", drop = FALSE], db, features = .syn_feats(),
                                          match_method = "flexible", rsid_policy = "require",
                                          use_xsv = FALSE, verbose = 0))
  expect_equal(none$match_tier, inp$flexible)
  expect_true(all(none$rsid_check[!is.na(none$rsid_check)] == "chip_none"))
  # Exact mode records the check too, and the rsID is never a lookup key: the
  # input rs9999 at 22-1000-C-T is not found by its rsID.
  ex <- .syn_annotate(db, "exact", "csv", use_xsv = FALSE)
  expect_equal(ex$rsid_check[inp$VarInfo == "22-1000-G-A"], c("same", "same"))
  expect_equal(ex$match_tier[inp$VarInfo == "22-1000-C-T"], "unmatched_ref")
  # The check is made whether or not rsid is a requested feature.
  expect_false("rsid" %in% names(rec))
  # The summary message reports the policy and the counts.
  expect_message(suppressWarnings(annotate_favor(inp[1:3, c("VarInfo", "rsID")], db,
    features = .syn_feats(), match_method = "flexible", rsid_policy = "require",
    use_xsv = FALSE, verbose = 1)), "rsID check \\(policy require, source rsID column\\): same=2, differs=1")
})

test_that("the xsv position join returns the same rows as the fread path", {
  skip_if_not(GLOWr:::.check_xsv_available(), "xsv not installed")
  db <- .syn_csv_db(withr::local_tempdir())
  a <- .syn_annotate(db, "flexible", "csv", use_xsv = FALSE)
  b <- .syn_annotate(db, "flexible", "csv", use_xsv = TRUE)
  expect_identical(b, a)
  # Both chunks contribute: rows at 1000 (chunk 1) and 5000 (chunk 2) resolve.
  expect_false(is.na(b$cadd_phred[b$VarInfo == "22-5000-CAT-C"]))
})

test_that("the Parquet backend returns what the CSV backend returns (G1, synthetic)", {
  skip_if_not_installed("arrow")
  csv_db <- .syn_csv_db(withr::local_tempdir())
  pq_db  <- .syn_parquet_db(withr::local_tempdir())
  for (method in c("flexible", "exact")) {
    a <- .syn_annotate(csv_db, method, "csv", use_xsv = FALSE)
    b <- .syn_annotate(pq_db, method, "parquet")
    expect_identical(b$match_tier, a$match_tier, info = method)
    expect_identical(b$favor_key, a$favor_key, info = method)
    expect_identical(b$rsid_check, a$rsid_check, info = method)
    for (f in .syn_feats()) expect_identical(b[[f]], a[[f]], info = paste(method, f))
  }
  # The rsID policy acts the same on both backends (the Parquet rsID is dbsnp.rsid).
  a <- .syn_annotate(csv_db, "flexible", "csv", use_xsv = FALSE, rsid_policy = "require")
  b <- .syn_annotate(pq_db, "flexible", "parquet", rsid_policy = "require")
  expect_identical(b$match_tier, a$match_tier)
  expect_identical(b$match_tier, .syn_inputs()$require)
  # Auto-detection picks each backend, and a mixed directory is refused.
  expect_identical(.syn_annotate(pq_db, "flexible", "auto")$match_tier,
                   .syn_inputs()$flexible)
  file.copy(file.path(csv_db, "chr22_1.csv"), pq_db)
  expect_error(.syn_annotate(pq_db, "flexible", "auto"), "both Parquet and CSV")
})

test_that("a chromosome without a database is 'no_database', with a warning", {
  db <- .syn_csv_db(withr::local_tempdir())
  expect_warning(
    res <- annotate_favor(data.frame(VarInfo = c("23-100-A-G", "22-1000-G-A")),
                          favor_db_path = db, features = .syn_feats(), verbose = 0),
    "no data for chromosome 23")
  expect_equal(res$match_tier, c("no_database", "exact"))
})

test_that("a missing chunk file, or one without the key columns, is fatal", {
  db <- .syn_csv_db(withr::local_tempdir())
  unlink(file.path(db, "chr22_2.csv"))
  expect_error(annotate_favor(data.frame(VarInfo = "22-5000-C-A"), favor_db_path = db,
                              features = .syn_feats(), verbose = 0),
               "chunk file not found")
  writeLines(c("variant_vcf,cadd_phred", "22-5000-C-A,1"), file.path(db, "chr22_2.csv"))
  expect_error(annotate_favor(data.frame(VarInfo = "22-5000-C-A"), favor_db_path = db,
                              features = .syn_feats(), verbose = 0),
               "lacks position")
})

test_that("na_allele_method = 'first' takes every feature from the first row", {
  db <- .syn_csv_db(withr::local_tempdir())
  res <- annotate_favor(data.frame(VarInfo = "22-1000-NA-NA"), favor_db_path = db,
                        features = .syn_feats(), na_allele_method = "first",
                        use_xsv = FALSE, verbose = 0)
  expect_equal(res$match_tier, "position_only")
  expect_equal(res$cadd_phred, 1.25)
  expect_equal(res$metasvm_pred, "D")
})

test_that("na_handling acts on features only and never changes the outcome", {
  db <- .syn_csv_db(withr::local_tempdir())
  inp <- data.frame(VarInfo = c("22-1000-G-A", "22-6000-A-G"))
  z <- annotate_favor(inp, favor_db_path = db, features = .syn_feats(),
                      na_handling = "zero", use_xsv = FALSE, verbose = 0)
  expect_equal(z$match_tier, c("exact", "uncovered"))
  expect_equal(z$cadd_phred, c(1.25, 0))
  expect_true(is.na(z$genecode_comprehensive_category[2]))  # strings keep NA
  d <- annotate_favor(inp, favor_db_path = db, features = .syn_feats(),
                      na_handling = "drop", use_xsv = FALSE, verbose = 0)
  expect_equal(d$VarInfo, "22-1000-G-A")
})

test_that("the nested FAVOR 2.0 fields serialize into the v1 layouts", {
  fr <- data.frame(position = 1:3)
  fr$.gh_id <- c("GH22I017115", NA, "GH22I000001")
  fr$.gh_score <- c(0.72, NA, 1.2)
  fr$.gh_targets <- list(data.frame(gene = c("TMEM121B", "IL17RA"), score = c(12.18, 3)),
                         NULL, data.frame(gene = character(0), score = numeric(0)))
  fr$.gei_tx <- list(
    data.frame(gene = c("CECR2", "CECR2"), transcript_id = c("ENST1.6", "ENST2.1"),
               location = c("exon8", "exon8"), hgvsc = c("c.G389A", "c.G812A"),
               hgvsp = c("p.R130H", "p.R271H")),
    data.frame(gene = "UNKNOWN", transcript_id = "", location = "", hgvsc = NA, hgvsp = ""),
    NULL)
  fr$.gi_genes <- list(c("ENST00000651146.1:c.*975C>T", "ENST00000465611.1:c.*2092C>T)", "GAB4"),
                       c("TCN2", "PES1"), character(0))
  s <- GLOWr:::.favor2_serialize(fr, 1:3, c("genehancer", "genecode_comprehensive_exonic_info",
                                            "genecode_comprehensive_info"))
  expect_equal(s$genehancer, c(
    "Name=0.72;genehancer_id=GH22I017115;connected_gene=TMEM121B;score=12.18;connected_gene=IL17RA;score=3.00",
    NA, "Name=1.20;genehancer_id=GH22I000001"))
  expect_equal(s$genecode_comprehensive_exonic_info, c(
    "CECR2:ENST1.6:exon8:c.G389A:p.R130H,CECR2:ENST2.1:exon8:c.G812A:p.R271H,",
    "UNKNOWN", NA))
  expect_equal(s$genecode_comprehensive_info, c("GAB4", "TCN2,PES1", NA))
  # STAARpipeline's enhancer masks read the first connected gene this way.
  expect_equal(strsplit(strsplit(s$genehancer[1], "=")[[1]][4], ";")[[1]][1], "TMEM121B")
})

test_that("a FAVOR 2.0 file lacking a requested feature's leaves skips that feature", {
  skip_if_not_installed("arrow")
  db <- .syn_parquet_db(withr::local_tempdir(), drop = "linsight")
  inp <- .syn_inputs()
  expect_message(
    res <- suppressWarnings(annotate_favor(inp[, "VarInfo", drop = FALSE], favor_db_path = db,
                          features = .syn_feats(), favor_db_format = "parquet",
                          match_method = "flexible", verbose = 1)),
    "lacks the columns of linsight")
  expect_false("linsight" %in% names(res))
  expect_equal(res$match_tier, inp$flexible)            # matching is unaffected
  expect_equal(res$cadd_phred[inp$VarInfo == "22-1000-C-T"], 1.25)
  # The key columns stay mandatory.
  db2 <- withr::local_tempdir()
  arrow::write_parquet(arrow::arrow_table(position = 1:2, ref_vcf = c("A", "C")),
                       file.path(db2, "chromosome_22.parquet"))
  expect_error(annotate_favor(inp[1, "VarInfo", drop = FALSE], favor_db_path = db2,
                              features = .syn_feats(), favor_db_format = "parquet", verbose = 0),
               "lacks the key columns")
})

test_that("the xsv path accepts a database directory whose path contains spaces", {
  skip_if_not(GLOWr:::.check_xsv_available(), "xsv not available")
  db <- file.path(withr::local_tempdir(), "favor db with spaces")
  dir.create(db)
  .syn_csv_db(db)
  inp <- .syn_inputs()
  a <- .syn_annotate(db, "flexible", "csv", use_xsv = TRUE)
  b <- .syn_annotate(db, "flexible", "csv", use_xsv = FALSE)
  expect_equal(a$match_tier, inp$flexible)
  expect_equal(a, b)
})

test_that("the backend is resolved from the directory or refused", {
  d <- withr::local_tempdir()
  expect_error(GLOWr:::.favor_resolve_backend(d, "auto"), "No FAVOR database files")
  file.create(file.path(d, "chr1_1.csv"))
  expect_equal(GLOWr:::.favor_resolve_backend(d, "auto"), "csv")
  expect_equal(GLOWr:::.favor_resolve_backend(d, "parquet"), "parquet")
})
