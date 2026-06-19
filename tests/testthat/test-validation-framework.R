#' Unit Tests for GLOWr Validation Framework
#'
#' This file tests the validation framework functions to ensure they work
#' correctly before using them to validate actual GLOWr implementations.
#'
#' File Log (reverse chronological order):
#' - 2025-10-18: Created by r-developer - Initial validation framework tests

# Test helper-validation.R functions
test_that("load_legacy_glow() loads legacy functions successfully", {
  # Skip if legacy directory not available
  legacy_dir <- file.path(Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"), "legacy-materials/code/GLOW_R_pacakge/GLOW/R")
  skip_if_not(dir.exists(legacy_dir), "Legacy directory not found")

  # Load legacy code
  expect_silent(legacy_env <- load_legacy_glow(verbose = FALSE))
  expect_true(is.environment(legacy_env))

  # Check that expected functions are loaded
  expected_functions <- c(
    "select_best_model",  # from helpers_optimalWeights.R
    "LASSOmodel",         # from helpers_optimalWeights.R
    "GLMmodel"            # from helpers_optimalWeights.R
  )

  for (fn in expected_functions) {
    expect_true(fn %in% ls(legacy_env),
                info = paste("Function", fn, "should be loaded"))
  }
})


test_that("generate_test_data() creates valid test data", {
  # Simple continuous data
  data1 <- generate_test_data(n = 100, p = 10, binary = FALSE, seed = 123)

  expect_type(data1, "list")
  expect_named(data1, c("G", "X", "Y", "true_B", "true_PI", "causal_idx", "maf", "description"))

  # Check dimensions
  expect_equal(nrow(data1$G), 100)
  expect_equal(ncol(data1$G), 10)
  expect_equal(nrow(data1$X), 100)
  expect_equal(length(data1$Y), 100)
  expect_equal(length(data1$true_B), 10)
  expect_equal(length(data1$true_PI), 10)

  # Check genotype values are valid (0, 1, 2)
  expect_true(all(data1$G %in% c(0, 1, 2)))

  # Check MAF in valid range
  expect_true(all(data1$maf > 0 & data1$maf < 1))

  # Check causal variants
  expect_true(all(data1$true_PI %in% c(0, 1)))
  expect_equal(sum(data1$true_PI), length(data1$causal_idx))

  # Binary data
  data2 <- generate_test_data(n = 200, p = 5, binary = TRUE, seed = 456)

  expect_true(all(data2$Y %in% c(0, 1)))

  # Rare variants
  data3 <- generate_test_data(n = 100, p = 10, rare = TRUE, seed = 789)

  expect_true(all(data3$maf < 0.01))

  # Correlated variants
  data4 <- generate_test_data(n = 150, p = 8, ld_structure = TRUE,
                               ld_strength = 0.5, seed = 101112)

  # Check that correlation structure exists (not a strict requirement, but likely)
  # Calculate correlation between adjacent variants
  if (ncol(data4$G) > 1) {
    cors <- sapply(1:(ncol(data4$G)-1), function(i) {
      cor(data4$G[, i], data4$G[, i+1])
    })
    # With LD structure, we expect some correlation (not necessarily all positive)
    expect_true(any(abs(cors) > 0.1))
  }
})


test_that("compare_with_legacy() correctly compares numeric vectors", {
  # Identical vectors
  v1 <- c(1.0, 2.0, 3.0)
  v2 <- c(1.0, 2.0, 3.0)

  result <- compare_with_legacy(v1, v2, tolerance = 1e-10, test_name = "identical")

  expect_true(result$passed)
  expect_equal(result$max_diff, 0)
  expect_equal(result$mean_diff, 0)

  # Nearly identical vectors (within tolerance)
  v3 <- c(1.0, 2.0, 3.0)
  v4 <- c(1.0 + 1e-12, 2.0, 3.0)

  result2 <- compare_with_legacy(v3, v4, tolerance = 1e-10, test_name = "near_identical")

  expect_true(result2$passed)
  expect_true(result2$max_diff < 1e-10)

  # Vectors exceeding tolerance
  v5 <- c(1.0, 2.0, 3.0)
  v6 <- c(1.0, 2.0, 3.1)

  result3 <- compare_with_legacy(v5, v6, tolerance = 1e-10, test_name = "different")

  expect_false(result3$passed)
  expect_true(result3$max_diff > 1e-10)
  expect_equal(result3$max_diff, 0.1)
})


test_that("compare_with_legacy() correctly compares matrices", {
  # Identical matrices
  m1 <- matrix(c(1, 2, 3, 4), 2, 2)
  m2 <- matrix(c(1, 2, 3, 4), 2, 2)

  result <- compare_with_legacy(m1, m2, tolerance = 1e-10, test_name = "matrix_identical")

  expect_true(result$passed)
  expect_equal(result$max_diff, 0)

  # Different matrices
  m3 <- matrix(c(1, 2, 3, 4), 2, 2)
  m4 <- matrix(c(1, 2, 3, 4.001), 2, 2)

  result2 <- compare_with_legacy(m3, m4, tolerance = 1e-10, test_name = "matrix_different")

  expect_false(result2$passed)
  expect_equal(result2$max_diff, 0.001)

  # Dimension mismatch
  m5 <- matrix(c(1, 2, 3, 4), 2, 2)
  m6 <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)

  result3 <- compare_with_legacy(m5, m6, tolerance = 1e-10, test_name = "matrix_dim_mismatch")

  expect_false(result3$passed)
  expect_match(result3$details, "Dimension mismatch")
})


test_that("compare_with_legacy() correctly compares lists", {
  # Identical lists
  list1 <- list(a = c(1, 2, 3), b = matrix(1:4, 2, 2))
  list2 <- list(a = c(1, 2, 3), b = matrix(1:4, 2, 2))

  result <- compare_with_legacy(list1, list2, tolerance = 1e-10, test_name = "list_identical")

  expect_true(result$passed)

  # Different list values
  list3 <- list(a = c(1, 2, 3), b = matrix(1:4, 2, 2))
  list4 <- list(a = c(1, 2, 3.001), b = matrix(1:4, 2, 2))

  result2 <- compare_with_legacy(list3, list4, tolerance = 1e-10, test_name = "list_different")

  expect_false(result2$passed)
  expect_true(length(result2$mismatches) > 0)

  # Different list names
  list5 <- list(a = c(1, 2, 3))
  list6 <- list(b = c(1, 2, 3))

  result3 <- compare_with_legacy(list5, list6, tolerance = 1e-10, test_name = "list_names")

  expect_false(result3$passed)
  expect_match(result3$details, "List names differ")
})


test_that("compare_with_legacy() handles NULL values", {
  # Both NULL
  result1 <- compare_with_legacy(NULL, NULL, test_name = "both_null")
  expect_true(result1$passed)
  expect_match(result1$details, "NULL")

  # One NULL
  result2 <- compare_with_legacy(NULL, c(1, 2, 3), test_name = "one_null")
  expect_false(result2$passed)
  expect_match(result2$details, "NULL")
})


test_that("compare_with_legacy() detects type mismatches", {
  v1 <- c(1, 2, 3)
  m1 <- matrix(1:3, 3, 1)

  result <- compare_with_legacy(v1, m1, test_name = "type_mismatch")

  expect_false(result$passed)
  expect_match(result$details, "Type mismatch")
})


test_that("validation_report() generates report file", {
  # Create some test results
  test_results <- list(
    test1 = list(
      passed = TRUE,
      max_diff = 0,
      mean_diff = 0,
      details = "Passed",
      mismatches = list()
    ),
    test2 = list(
      passed = FALSE,
      max_diff = 0.001,
      mean_diff = 0.0005,
      details = "Failed: max diff exceeded",
      mismatches = list(position = 5)
    )
  )

  # Create temporary report directory
  temp_dir <- tempdir()
  report_dir <- file.path(temp_dir, "validation_reports")
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

  # Generate report
  report_path <- validation_report(
    test_results,
    report_name = "test_report",
    report_dir = report_dir,
    format = "txt"
  )

  # Check file was created
  expect_true(file.exists(report_path))

  # Read and check content
  report_content <- readLines(report_path)

  expect_true(any(grepl("GLOWr Validation Report", report_content)))
  expect_true(any(grepl("Total tests: 2", report_content)))
  expect_true(any(grepl("Passed: 1", report_content)))
  expect_true(any(grepl("Failed: 1", report_content)))

  # Clean up
  unlink(report_path)
})


test_that("Test fixtures load successfully", {
  fixture_dir <- "../../inst/validation/fixtures"

  # Skip if fixtures not available
  skip_if_not(dir.exists(fixture_dir), "Fixture directory not found")

  # List all fixtures
  fixtures <- list.files(fixture_dir, pattern = "\\.rds$", full.names = TRUE)

  # Skip if no fixtures
  skip_if(length(fixtures) == 0, "No fixtures found")

  # Expected fixtures
  expected_fixtures <- c(
    "test_data_simple.rds",
    "test_data_correlated.rds",
    "test_data_binary.rds",
    "test_data_continuous.rds",
    "test_data_rare.rds",
    "test_data_sparse.rds"
  )

  # Test each fixture
  for (fixture_name in expected_fixtures) {
    fixture_path <- file.path(fixture_dir, fixture_name)

    if (file.exists(fixture_path)) {
      # Load fixture
      data <- readRDS(fixture_path)

      # Check structure (fixture_name context for debugging)
      expect_type(data, "list")
      expect_named(data, c("G", "X", "Y", "true_B", "true_PI", "causal_idx", "maf", "description"))

      # Check dimensions are consistent
      n <- nrow(data$G)
      p <- ncol(data$G)

      expect_equal(nrow(data$X), n)
      expect_equal(length(data$Y), n)
      expect_equal(length(data$true_B), p)
      expect_equal(length(data$true_PI), p)
      expect_equal(length(data$maf), p)

      # Check genotypes are valid
      expect_true(all(data$G %in% c(0, 1, 2)))

      # Check description exists
      expect_type(data$description, "character")
      expect_true(nchar(data$description) > 0)

      message(sprintf("✓ Fixture %s: %s", fixture_name, data$description))
    } else {
      message(sprintf("⚠ Fixture %s not found", fixture_name))
    }
  }
})


test_that("Test fixtures have expected properties", {
  fixture_dir <- "../../inst/validation/fixtures"
  skip_if_not(dir.exists(fixture_dir), "Fixture directory not found")

  # Test simple fixture
  simple_path <- file.path(fixture_dir, "test_data_simple.rds")
  if (file.exists(simple_path)) {
    simple <- readRDS(simple_path)
    expect_equal(nrow(simple$G), 100)
    expect_equal(ncol(simple$G), 10)
  }

  # Test binary fixture
  binary_path <- file.path(fixture_dir, "test_data_binary.rds")
  if (file.exists(binary_path)) {
    binary <- readRDS(binary_path)
    expect_true(all(binary$Y %in% c(0, 1)))
    expect_true(sum(binary$Y) > 0)  # Has some cases
    expect_true(sum(binary$Y) < length(binary$Y))  # Has some controls
  }

  # Test rare variant fixture
  rare_path <- file.path(fixture_dir, "test_data_rare.rds")
  if (file.exists(rare_path)) {
    rare <- readRDS(rare_path)
    expect_true(all(rare$maf < 0.01))
    expect_true(all(rare$Y %in% c(0, 1)))  # Should be binary
  }

  # Test sparse fixture
  sparse_path <- file.path(fixture_dir, "test_data_sparse.rds")
  if (file.exists(sparse_path)) {
    sparse <- readRDS(sparse_path)
    expect_equal(ncol(sparse$G), 50)  # Should have 50 variants
    expect_true(sum(sparse$true_PI) <= 2)  # Sparse: 1-2 causal
  }
})


test_that("Validation framework integration test", {
  # This test demonstrates a complete validation workflow

  # Skip if legacy code not available
  legacy_dir <- file.path(Sys.getenv("GLOW_LEGACY_ROOT", unset = "/nonexistent"), "legacy-materials/code/GLOW_R_pacakge/GLOW/R")
  skip_if_not(dir.exists(legacy_dir), "Legacy directory not found")

  # 1. Load legacy code
  legacy_env <- load_legacy_glow(verbose = FALSE)
  expect_true(is.environment(legacy_env))

  # 2. Generate test data
  test_data <- generate_test_data(n = 50, p = 5, seed = 999)

  # 3. Create some mock results to compare
  # (In real validation, these would be actual function outputs)
  glowr_result <- list(
    stat = 1.234,
    pval = 0.456,
    scores = c(0.1, 0.2, 0.3, 0.4, 0.5)
  )

  legacy_result <- list(
    stat = 1.234,
    pval = 0.456,
    scores = c(0.1, 0.2, 0.3, 0.4, 0.5)
  )

  # 4. Compare results
  comparison <- compare_with_legacy(glowr_result, legacy_result, tolerance = 1e-10)

  expect_true(comparison$passed)

  # 5. Generate validation report
  test_results <- list(
    integration_test = comparison
  )

  temp_dir <- tempdir()
  report_path <- validation_report(
    test_results,
    report_name = "integration_test",
    report_dir = temp_dir
  )

  expect_true(file.exists(report_path))

  # Clean up
  unlink(report_path)

  message("✓ Integration test completed successfully")
})
