# Test file for plink_to_gds() function

# Helper to get test data path
get_test_plink_prefix <- function() {
  file.path(
    testthat::test_path(), "..", "..", "..", "..",
    "data", "ALS_GWAS_sample_data",
    "ALS_GWAS_chr21_hg38_rename_first800"
  )
}

# Check if SeqArray is available (required for these tests)
skip_if_no_seqarray <- function() {
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    testthat::skip("SeqArray package not available")
  }
}

# ========== Basic Functionality Tests ==========

test_that("plink_to_gds converts PLINK files to GDS", {
  skip_if_no_seqarray()

  # Get test data prefix
  plink_prefix <- get_test_plink_prefix()

  # Skip if test data not available
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Create temporary GDS output
  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Convert PLINK to GDS
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_gds,
    verbose = 0
  )

  # Verify file was created
  expect_true(file.exists(temp_gds))

  # Verify function returns output path
  expect_equal(result, temp_gds)

  # Open GDS and verify contents
  gds <- SeqArray::seqOpen(temp_gds, readonly = TRUE)

  tryCatch({
    # Verify 800 variants (as expected from test data)
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    expect_equal(n_variants, 800)

    # Verify all variants are from chr21
    chr <- SeqArray::seqGetData(gds, "chromosome")
    expect_true(all(chr == "21" | chr == 21))

    # Verify basic data structures
    expect_type(SeqArray::seqGetData(gds, "position"), "integer")
    expect_type(SeqArray::seqGetData(gds, "$ref"), "character")

  }, finally = {
    SeqArray::seqClose(gds)
  })
})

test_that("plink_to_gds creates GDS with correct first variant", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  plink_to_gds(plink_prefix, temp_gds, verbose = 0)

  # Open GDS and check first variant
  gds <- SeqArray::seqOpen(temp_gds, readonly = TRUE)

  tryCatch({
    chr <- SeqArray::seqGetData(gds, "chromosome")[1]
    pos <- SeqArray::seqGetData(gds, "position")[1]
    ref <- SeqArray::seqGetData(gds, "$ref")[1]
    alt_raw <- SeqArray::seqGetData(gds, "$alt")

    # Handle alt (may be list or character)
    if (is.list(alt_raw)) {
      alt <- alt_raw[[1]][1]
    } else {
      alt <- strsplit(as.character(alt_raw[1]), ",")[[1]][1]
    }

    # Expected first variant: 21-13464068-C-T
    expect_true(chr == "21" | chr == 21)
    expect_equal(pos, 13464068)
    expect_equal(ref, "C")
    expect_equal(alt, "T")

  }, finally = {
    SeqArray::seqClose(gds)
  })
})

# ========== Input Validation Tests ==========

test_that("plink_to_gds handles missing files with informative error", {
  skip_if_no_seqarray()

  # Test with non-existent prefix
  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  expect_error(
    plink_to_gds(
      plink_prefix = "/nonexistent/path/to/data",
      output_gds = temp_gds,
      verbose = 0
    ),
    "Missing required PLINK files"
  )
})

test_that("plink_to_gds handles missing single PLINK file", {
  skip_if_no_seqarray()

  # Create temp directory with incomplete PLINK files
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  temp_prefix <- file.path(temp_dir, "incomplete")
  temp_gds <- tempfile(fileext = ".gds")

  # Create only .bed file (missing .bim and .fam)
  file.create(paste0(temp_prefix, ".bed"))

  expect_error(
    plink_to_gds(temp_prefix, temp_gds, verbose = 0),
    "Missing required PLINK files"
  )
})

test_that("plink_to_gds validates input types", {
  skip_if_no_seqarray()

  temp_gds <- tempfile(fileext = ".gds")

  # Invalid plink_prefix type
  expect_error(
    plink_to_gds(plink_prefix = 123, output_gds = temp_gds),
    "plink_prefix must be a single character string"
  )

  # Invalid output_gds type
  expect_error(
    plink_to_gds(plink_prefix = "test", output_gds = c("a", "b")),
    "output_gds must be a single character string"
  )
})

test_that("plink_to_gds checks output directory exists", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  expect_error(
    plink_to_gds(
      plink_prefix = plink_prefix,
      output_gds = "/nonexistent/dir/output.gds",
      verbose = 0
    ),
    "Output directory does not exist"
  )
})

# ========== Extension Stripping Tests ==========

test_that("plink_to_gds strips .bed extension from prefix", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Provide prefix with .bed extension
  result <- plink_to_gds(
    plink_prefix = paste0(plink_prefix, ".bed"),
    output_gds = temp_gds,
    verbose = 0
  )

  expect_true(file.exists(temp_gds))
  expect_equal(result, temp_gds)
})

test_that("plink_to_gds strips .bim extension from prefix", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Provide prefix with .bim extension
  result <- plink_to_gds(
    plink_prefix = paste0(plink_prefix, ".bim"),
    output_gds = temp_gds,
    verbose = 0
  )

  expect_true(file.exists(temp_gds))
})

test_that("plink_to_gds strips .fam extension from prefix", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Provide prefix with .fam extension
  result <- plink_to_gds(
    plink_prefix = paste0(plink_prefix, ".fam"),
    output_gds = temp_gds,
    verbose = 0
  )

  expect_true(file.exists(temp_gds))
})

# ========== Verbosity Tests ==========

test_that("plink_to_gds is silent when verbose=0", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Capture messages
  output <- capture.output(
    result <- plink_to_gds(plink_prefix, temp_gds, verbose = 0),
    type = "message"
  )

  # Should produce no messages
  expect_length(output, 0)
})

test_that("plink_to_gds reports counts when verbose=1", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Capture messages
  output <- capture.output(
    result <- plink_to_gds(plink_prefix, temp_gds, verbose = 1),
    type = "message"
  )

  # Should report conversion and counts
  output_text <- paste(output, collapse = "\n")
  expect_match(output_text, "Converting PLINK to GDS")
  expect_match(output_text, "800 variants")
  expect_match(output_text, "samples")
})

test_that("plink_to_gds warns when overwriting existing GDS", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Create initial GDS
  plink_to_gds(plink_prefix, temp_gds, verbose = 0)

  # Overwrite should produce warning
  expect_warning(
    plink_to_gds(plink_prefix, temp_gds, verbose = 1),
    "already exists and will be overwritten"
  )
})

# ========== Edge Cases ==========

test_that("plink_to_gds handles paths with spaces", {
  skip_if_no_seqarray()
  skip_on_ci()  # Skip on CI due to file path complications

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Create temp directory with space in name
  temp_dir <- file.path(tempdir(), "test dir with spaces")
  dir.create(temp_dir, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  temp_gds <- file.path(temp_dir, "output.gds")

  result <- plink_to_gds(plink_prefix, temp_gds, verbose = 0)

  expect_true(file.exists(temp_gds))
})

test_that("plink_to_gds handles relative paths", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Use current directory for output
  temp_gds <- "test_output.gds"
  on.exit(unlink(temp_gds), add = TRUE)

  result <- plink_to_gds(plink_prefix, temp_gds, verbose = 0)

  expect_true(file.exists(temp_gds))
})

# ========== Chromosome Split Tests ==========

test_that("plink_to_gds with split_by_chr creates per-chromosome files", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Create temp output directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Test data is chr21 only, so should create one file
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,
    split_by_chr = TRUE,
    verbose = 0
  )

  # Should return character vector
  expect_type(result, "character")
  expect_length(result, 1)

  # File should be named with chr21
  expect_true(grepl("_chr21\\.gds$", result))

  # File should exist and contain correct variants
  expect_true(file.exists(result))

  gds <- SeqArray::seqOpen(result, readonly = TRUE)
  tryCatch({
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    expect_equal(n_variants, 800)
  }, finally = {
    SeqArray::seqClose(gds)
  })
})

test_that("plink_to_gds with split_by_chr and specific chromosomes", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Request chr21 (exists) and chr22 (doesn't exist in test data)
  expect_warning(
    result <- plink_to_gds(
      plink_prefix = plink_prefix,
      output_gds = temp_dir,
      split_by_chr = TRUE,
      chromosomes = c(21, 22),
      verbose = 0
    ),
    "not found in data"
  )

  # Should only create chr21 file
  expect_length(result, 1)
  expect_true(grepl("_chr21\\.gds$", result))
})

test_that("plink_to_gds with custom output_prefix", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,
    split_by_chr = TRUE,
    output_prefix = "custom_prefix",
    verbose = 0
  )

  # File should use custom prefix
  expect_true(grepl("^custom_prefix_chr", basename(result)))
})

# ========== Region Filter Tests ==========

test_that("plink_to_gds with single region (shorthand format)", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Filter to a region that contains variants
  # Test data is chr21, positions around 13464068 to 14454530
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,
    regions = list(chr = "21", start = 13400000, end = 14000000),
    verbose = 0
  )

  expect_type(result, "character")
  expect_length(result, 1)

  # File should be named with region coordinates
  expect_true(grepl("_chr21_13400000_14000000\\.gds$", result))
  expect_true(file.exists(result))

  # Verify variants are within region
  gds <- SeqArray::seqOpen(result, readonly = TRUE)
  tryCatch({
    positions <- SeqArray::seqGetData(gds, "position")
    expect_true(all(positions >= 13400000 & positions <= 14000000))
    expect_true(length(positions) > 0)
  }, finally = {
    SeqArray::seqClose(gds)
  })
})

test_that("plink_to_gds with multiple regions (list-of-lists format)", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Multiple regions on chr21
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,
    regions = list(
      list(chr = "21", start = 13400000, end = 13600000),
      list(chr = "21", start = 14000000, end = 14200000)
    ),
    verbose = 0
  )

  # Should create two files (or fewer if some regions have no variants)
  expect_type(result, "character")
  expect_true(length(result) >= 1)

  # All files should exist
  expect_true(all(file.exists(result)))
})

# ========== Mutual Exclusivity Tests ==========

test_that("plink_to_gds errors when both split_by_chr and regions specified", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  expect_error(
    plink_to_gds(
      plink_prefix = plink_prefix,
      output_gds = temp_dir,
      split_by_chr = TRUE,
      regions = list(chr = "21", start = 1, end = 1000000),
      verbose = 0
    ),
    "mutually exclusive"
  )
})

# ========== Region Validation Tests ==========

test_that("plink_to_gds validates region structure", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Missing 'end'
  expect_error(
    plink_to_gds(
      plink_prefix = plink_prefix,
      output_gds = temp_dir,
      regions = list(chr = "21", start = 1000000),
      verbose = 0
    ),
    "chr, start, end"
  )

  # start > end
  expect_error(
    plink_to_gds(
      plink_prefix = plink_prefix,
      output_gds = temp_dir,
      regions = list(chr = "21", start = 5000000, end = 1000000),
      verbose = 0
    ),
    "cannot be greater than end"
  )
})

test_that("plink_to_gds validates split_by_chr type", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  expect_error(
    plink_to_gds(
      plink_prefix = plink_prefix,
      output_gds = temp_dir,
      split_by_chr = "yes",
      verbose = 0
    ),
    "split_by_chr must be a single logical value"
  )
})

# ========== Directory Creation Tests ==========

test_that("plink_to_gds creates output directory in split mode if parent exists", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Create parent directory but not the target
  parent_dir <- tempfile()
  dir.create(parent_dir)
  on.exit(unlink(parent_dir, recursive = TRUE), add = TRUE)

  new_dir <- file.path(parent_dir, "new_output")
  expect_false(dir.exists(new_dir))

  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = new_dir,
    split_by_chr = TRUE,
    verbose = 0
  )

  # Directory should have been created
  expect_true(dir.exists(new_dir))
  expect_true(length(result) >= 1)
})

# ========== Directory Auto-Filename Tests ==========

test_that("plink_to_gds auto-generates filename when output_gds is directory", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Create temp output directory
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Provide directory instead of file path
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,  # Directory, not file path
    verbose = 0
  )

  # Should return auto-generated filename
  expect_type(result, "character")
  expect_length(result, 1)

  # Filename should be based on PLINK prefix
  expect_true(grepl("ALS_GWAS_chr21_hg38_rename_first800\\.gds$", result))
  expect_true(file.exists(result))

  # File should be in the provided directory
  expect_equal(dirname(result), temp_dir)

  # Verify contents
  gds <- SeqArray::seqOpen(result, readonly = TRUE)
  tryCatch({
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    expect_equal(n_variants, 800)
  }, finally = {
    SeqArray::seqClose(gds)
  })
})

# ========== Backward Compatibility Tests ==========

test_that("plink_to_gds default mode unchanged with new parameters", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Call with new parameters at default values
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_gds,
    split_by_chr = FALSE,
    chromosomes = NULL,
    regions = NULL,
    output_prefix = NULL,
    verbose = 0
  )

  # Should return single path (not vector)
  expect_type(result, "character")
  expect_length(result, 1)
  expect_equal(result, temp_gds)
  expect_true(file.exists(temp_gds))

  # Verify contents match expected
  gds <- SeqArray::seqOpen(temp_gds, readonly = TRUE)
  tryCatch({
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    expect_equal(n_variants, 800)
  }, finally = {
    SeqArray::seqClose(gds)
  })
})

# ========== Version 2: chr.conv, parallel, ... Parameter Tests ==========

test_that(".convert_chromosome_names applies correct conversion", {
  # Test with chr.conv = TRUE (should convert)
  expect_equal(
    GLOWr:::.convert_chromosome_names(c("1", "22", "23", "24", "25", "26"), TRUE),
    c("1", "22", "X", "Y", "XY", "MT")
  )

  # Test with chr.conv = FALSE (should not convert)
  expect_equal(
    GLOWr:::.convert_chromosome_names(c("1", "22", "23", "24"), FALSE),
    c("1", "22", "23", "24")
  )

  # Test mixed - only numeric sex/mito chromosomes converted
  expect_equal(
    GLOWr:::.convert_chromosome_names(c("21", "X", "23"), TRUE),
    c("21", "X", "X")  # "23" -> "X", "X" stays "X"
  )
})

test_that("plink_to_gds passes chr.conv to seqBED2GDS", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  # Test with chr.conv = TRUE (default)
  temp_gds_true <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds_true), add = TRUE)

  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_gds_true,
    chr.conv = TRUE,
    verbose = 0
  )

  expect_true(file.exists(temp_gds_true))

  # Test with chr.conv = FALSE
  temp_gds_false <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds_false), add = TRUE)

  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_gds_false,
    chr.conv = FALSE,
    verbose = 0
  )

  expect_true(file.exists(temp_gds_false))

  # Both should have 800 variants (test data is chr21, no 23/24/25/26)
  gds1 <- SeqArray::seqOpen(temp_gds_true, readonly = TRUE)
  gds2 <- SeqArray::seqOpen(temp_gds_false, readonly = TRUE)
  tryCatch({
    expect_equal(
      length(SeqArray::seqGetData(gds1, "variant.id")),
      length(SeqArray::seqGetData(gds2, "variant.id"))
    )
  }, finally = {
    SeqArray::seqClose(gds1)
    SeqArray::seqClose(gds2)
  })
})

test_that("plink_to_gds passes parallel parameter to seqBED2GDS", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Test with parallel = FALSE (default)
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_gds,
    parallel = FALSE,
    verbose = 0
  )

  expect_true(file.exists(temp_gds))

  # Verify contents correct
  gds <- SeqArray::seqOpen(temp_gds, readonly = TRUE)
  tryCatch({
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    expect_equal(n_variants, 800)
  }, finally = {
    SeqArray::seqClose(gds)
  })
})

test_that("plink_to_gds passes ... to seqBED2GDS", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_gds), add = TRUE)

  # Test passing additional parameters via ...
  # Using optimize = FALSE as a test parameter
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_gds,
    verbose = 0,
    optimize = FALSE  # This should be passed to seqBED2GDS via ...
  )

  expect_true(file.exists(temp_gds))

  # Verify contents correct
  gds <- SeqArray::seqOpen(temp_gds, readonly = TRUE)
  tryCatch({
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    expect_equal(n_variants, 800)
  }, finally = {
    SeqArray::seqClose(gds)
  })
})

test_that("plink_to_gds split mode with chr.conv applies conversion to chromosome names", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Test with chr.conv = TRUE (default) in split mode
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,
    split_by_chr = TRUE,
    chr.conv = TRUE,
    verbose = 0
  )

  # Test data is chr21 only, should create one file
  expect_length(result, 1)
  expect_true(file.exists(result))
  expect_true(grepl("_chr21\\.gds$", result))
})

test_that("plink_to_gds regions mode with chr.conv applies conversion", {
  skip_if_no_seqarray()

  plink_prefix <- get_test_plink_prefix()
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    skip("Test data not available")
  }

  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Test with chr.conv = TRUE in regions mode
  result <- plink_to_gds(
    plink_prefix = plink_prefix,
    output_gds = temp_dir,
    regions = list(chr = "21", start = 13400000, end = 14000000),
    chr.conv = TRUE,
    verbose = 0
  )

  expect_length(result, 1)
  expect_true(file.exists(result))
})
