########## Tests for extract_sample_size_from_text() ##########
#
# This file tests the enhanced sample size extraction function with intelligent
# group detection and summing capabilities.

test_that("extract_sample_size_from_text handles single numeric values", {
  # Direct numeric input
  expect_equal(
    GLOWr:::extract_sample_size_from_text(10000, verbose = FALSE),
    10000
  )

  # Single number as text
  expect_equal(
    GLOWr:::extract_sample_size_from_text("10000", verbose = FALSE),
    10000
  )

  # Single number with commas
  expect_equal(
    GLOWr:::extract_sample_size_from_text("10,000", verbose = FALSE),
    10000
  )

  # Single number with description
  expect_equal(
    GLOWr:::extract_sample_size_from_text("394,929 European ancestry individuals", verbose = FALSE),
    394929
  )
})


test_that("extract_sample_size_from_text handles gender-stratified samples (SUMMING)", {
  # Main test case from user report
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text(
      "76,067 European ancestry women, 66,420 European ancestry men",
      verbose = TRUE  # Should warn
    )
  )
  expect_equal(result, 142487)

  # Additional gender patterns
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("10,000 women, 15,000 men", verbose = FALSE)),
    25000
  )

  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("5000 females and 4000 males", verbose = FALSE)),
    9000
  )

  # Should still work with "female" and "male" singular
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("100 female, 200 male participants", verbose = FALSE)),
    300
  )
})


test_that("extract_sample_size_from_text handles case-control studies (SUMMING)", {
  # Explicit plus operator
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("5,000 cases + 5,000 controls", verbose = FALSE)),
    10000
  )

  # Cases and controls with "and"
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("10,000 cases and 15,000 controls", verbose = FALSE)),
    25000
  )

  # Just "cases" and "controls" keywords
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("2,500 cases, 7,500 controls", verbose = FALSE)),
    10000
  )

  # Singular forms
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("100 case and 300 control samples", verbose = FALSE)),
    400
  )
})


test_that("extract_sample_size_from_text handles 'and' conjunction between large numbers", {
  # "and" between two large numbers should trigger summing
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("10,000 and 15,000 participants", verbose = FALSE)
  )
  expect_equal(result, 25000)

  # "and" with smaller numbers should NOT trigger summing (likely not groups)
  # Should take max instead
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("50 and 100", verbose = FALSE)
  )
  expect_equal(result, 100)
})


test_that("extract_sample_size_from_text takes maximum when no summing indicators", {
  # Multiple numbers without clear indicators -> take max
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("Study with 10,000 in 2019", verbose = FALSE)
  )
  expect_equal(result, 10000)

  # Year should not be summed
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("5000 participants, data from 2020", verbose = FALSE)
  )
  expect_equal(result, 5000)
})


test_that("extract_sample_size_from_text handles NA and empty strings", {
  expect_true(is.na(GLOWr:::extract_sample_size_from_text(NA, verbose = FALSE)))
  expect_true(is.na(GLOWr:::extract_sample_size_from_text("", verbose = FALSE)))
  expect_true(is.na(GLOWr:::extract_sample_size_from_text("   ", verbose = FALSE)))
  expect_true(is.na(GLOWr:::extract_sample_size_from_text("No numbers here", verbose = FALSE)))
})


test_that("extract_sample_size_from_text handles vectorized input", {
  input <- c(
    "10,000",
    "5,000 cases + 5,000 controls",
    "76,067 women, 66,420 men",
    "394,929 individuals"
  )

  result <- suppressWarnings(GLOWr:::extract_sample_size_from_text(input, verbose = FALSE))

  expect_equal(result[1], 10000)
  expect_equal(result[2], 10000)
  expect_equal(result[3], 142487)
  expect_equal(result[4], 394929)
})


test_that("extract_sample_size_from_text produces warnings when verbose=TRUE", {
  # Should warn about summing
  expect_warning(
    GLOWr:::extract_sample_size_from_text("5,000 women, 4,000 men", verbose = TRUE),
    "Sample size extraction: Summed"
  )

  # "and" between large numbers triggers summing (not max warning)
  expect_warning(
    GLOWr:::extract_sample_size_from_text("Study with 10,000 and 20,000", verbose = TRUE),
    "'and' conjunction between large numbers"
  )

  # Multiple large numbers WITHOUT indicators should warn about taking max
  expect_warning(
    GLOWr:::extract_sample_size_from_text("Study 10,000 participants 20,000", verbose = TRUE),
    "Multiple large numbers found"
  )
})


test_that("extract_sample_size_from_text suppresses warnings when verbose=FALSE", {
  # Should NOT warn when verbose = FALSE
  expect_no_warning(
    GLOWr:::extract_sample_size_from_text("5,000 women, 4,000 men", verbose = FALSE)
  )

  expect_no_warning(
    GLOWr:::extract_sample_size_from_text("Study with 10,000 in 2019", verbose = FALSE)
  )
})


test_that("extract_sample_size_from_text handles edge cases", {
  # Very large numbers
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("1,000,000 participants", verbose = FALSE)),
    1000000
  )

  # Multiple commas
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("100,000 women, 200,000 men", verbose = FALSE)),
    300000
  )

  # Mixed formats
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("5000 cases + 10,000 controls", verbose = FALSE)),
    15000
  )
})


test_that("extract_sample_size_from_text handles complex real-world examples", {
  # Example 1: Multi-ancestry studies
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text(
      "Meta-analysis of 76,067 European women and 66,420 European men",
      verbose = FALSE
    )
  )
  expect_equal(result, 142487)

  # Example 2: Case-control with detailed description
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text(
      "Discovery cohort: 2,678 cases and 8,534 controls",
      verbose = FALSE
    )
  )
  expect_equal(result, 11212)

  # Example 3: Study with year but no summing indicators
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text(
      "Cohort study of 50,000 participants (2015-2020)",
      verbose = FALSE
    )
  )
  expect_equal(result, 50000)  # Should take 50,000, not sum with years
})


test_that("extract_sample_size_from_text prioritizes correct summing indicators", {
  # Explicit + should take precedence
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("5,000 + 10,000", verbose = FALSE)
  )
  expect_equal(result, 15000)

  # Gender indicators should trigger summing
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("Data from 2020 with 10,000 women, 15,000 men", verbose = FALSE)
  )
  expect_equal(result, 25000)  # Should sum 10,000 and 15,000, ignore 2020
})


test_that("extract_sample_size_from_text handles small numbers correctly", {
  # Small numbers with "and" should NOT be summed (likely not sample size groups)
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("Between 50 and 100", verbose = FALSE)
  )
  expect_equal(result, 100)

  # But small numbers with explicit indicators SHOULD be summed
  result <- suppressWarnings(
    GLOWr:::extract_sample_size_from_text("50 cases + 100 controls", verbose = FALSE)
  )
  expect_equal(result, 150)
})


test_that("extract_sample_size_from_text backward compatibility", {
  # Ensure old behavior is preserved for simple cases
  expect_equal(
    GLOWr:::extract_sample_size_from_text("10000", verbose = FALSE),
    10000
  )

  # Old "+" behavior should still work
  expect_equal(
    suppressWarnings(GLOWr:::extract_sample_size_from_text("5000 + 5000", verbose = FALSE)),
    10000
  )
})


test_that("extract_sample_size_from_text warning messages contain useful information", {
  # Check that warnings include the original text and reason
  expect_warning(
    GLOWr:::extract_sample_size_from_text("5,000 women, 4,000 men", verbose = TRUE),
    "Original text.*women.*men"
  )

  expect_warning(
    GLOWr:::extract_sample_size_from_text("5,000 women, 4,000 men", verbose = TRUE),
    "Reason:.*gender"
  )

  expect_warning(
    GLOWr:::extract_sample_size_from_text("5,000 women, 4,000 men", verbose = TRUE),
    "Please verify"
  )
})
