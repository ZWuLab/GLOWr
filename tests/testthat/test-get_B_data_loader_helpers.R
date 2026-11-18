########## Unit Tests for the Helpers of prepare_B_training_data() ##########
#
# Tests of the internal helpers in R/get_B_data_loader_helpers.R that the two data loaders
# share (prepare_B_training_data() and prepare_PI_case_data()).

filter_test_data <- function() {
  data.frame(
    rsID = paste0("rs", 1:6),
    N = c(1000, 100000, 250000, 1000, 5000, NA),
    MAF = c(0.5, 0.25, 0.1, 0.05, 0.01, 0.2),
    TRAIT = c("ALS", "ALS", "height", "007", "7", "ALS"),
    stringsAsFactors = FALSE
  )
}

test_that(".filter_values_for_column converts text for a numeric column only", {
  f <- GLOWr:::.filter_values_for_column
  expect_identical(f(c("1e3", "100000", "0.50"), c(1, 2), "N"), c(1000, 100000, 0.5))
  expect_identical(f(c("1", NA), c(1, 2), "N"), c(1, NA))
  # A character column keeps its text, so a label such as "007" is not turned into 7.
  expect_identical(f(c("007", "7"), c("a", "b"), "TRAIT"), c("007", "7"))
  # Numeric values are returned as given.
  expect_identical(f(c(1000, 5000), c(1, 2), "N"), c(1000, 5000))
  expect_identical(f(1:5, c("1", "X"), "CHR"), 1:5)
  expect_error(f(c("1000", "many"), c(1, 2), "N"),
               "column 'N' has value\\(s\\) that are not numbers \\('many'\\)")
})

test_that("an inclusion filter given as text matches a numeric column", {
  d <- filter_test_data()
  # 100000 as text: R alone compares it with "1e+05" and matches nothing.
  expect_false(100000 %in% "100000")
  kept <- GLOWr:::.apply_general_filters(d, include = list(N = c("100000", "1e3")), verbose = 0)
  expect_identical(kept$rsID, c("rs1", "rs2", "rs4"))
  # The same filter given as numbers keeps the same rows.
  kept_num <- GLOWr:::.apply_general_filters(d, include = list(N = c(100000, 1000)), verbose = 0)
  expect_identical(kept_num$rsID, kept$rsID)
  kept_maf <- GLOWr:::.apply_general_filters(d, include = list(MAF = "0.50"), verbose = 0)
  expect_identical(kept_maf$rsID, "rs1")
})

test_that("an exclusion filter given as text matches a numeric column", {
  d <- filter_test_data()
  left <- GLOWr:::.apply_general_filters(d, exclude = list(N = "250000"), verbose = 0)
  # The row with a missing N is kept: an exclusion removes explicit matches only.
  expect_identical(left$rsID, c("rs1", "rs2", "rs4", "rs5", "rs6"))
})

test_that("text labels of a character column are still compared as text", {
  d <- filter_test_data()
  kept <- GLOWr:::.apply_general_filters(d, include = list(TRAIT = "007"),
                                         pattern_matching = FALSE, verbose = 0)
  expect_identical(kept$rsID, "rs4")
})

test_that("a filter value that is not a number stops for a numeric column", {
  d <- filter_test_data()
  expect_error(GLOWr:::.apply_general_filters(d, include = list(N = "large"), verbose = 0),
               "column 'N' has value\\(s\\) that are not numbers")
  expect_error(GLOWr:::.apply_general_filters(d, exclude = list(MAF = c("0.1", "rare")), verbose = 0),
               "column 'MAF' has value\\(s\\) that are not numbers \\('rare'\\)")
})

test_that("text values for a 64-bit integer column are compared digit for digit", {
  skip_if_not_installed("bit64")
  f <- GLOWr:::.filter_values_for_column
  # two neighbours beyond 2^53, which a double cannot tell apart
  ids <- bit64::as.integer64(c("9007199254740992", "9007199254740993"))
  expect_true(as.numeric("9007199254740992") == as.numeric("9007199254740993"))
  expect_identical(ids %in% f("9007199254740993", ids, "tag"), c(FALSE, TRUE))
  expect_identical(ids %in% f(" 9007199254740992 ", ids, "tag"), c(TRUE, FALSE))
  expect_error(f(c("9007199254740993", "1e3"), ids, "tag"),
               "column 'tag' has value\\(s\\) that are not whole numbers in digits \\('1e3'\\)")
  # values that are not text are returned as given
  expect_identical(f(ids[2L], ids, "tag"), ids[2L])

  d <- data.frame(rsID = c("rs1", "rs2"), N = c(1000, 2000), stringsAsFactors = FALSE)
  d$tag <- ids
  kept <- GLOWr:::.apply_general_filters(d, include = list(tag = "9007199254740993"), verbose = 0)
  expect_identical(kept$rsID, "rs2")
  left <- GLOWr:::.apply_general_filters(d, exclude = list(tag = "9007199254740993"), verbose = 0)
  expect_identical(left$rsID, "rs1")
})

test_that("the public B loader keeps the right row of a 64-bit integer column", {
  skip_if_not_installed("bit64")
  d <- data.frame(rsID = c("rs1", "rs2"), CHR = c(1, 2), POS = c(101, 202), MAF = c(.1, .2),
                  BETA = c(.1, .2), P = c(.001, .002), N = c(1e4, 1e4))
  d$tag <- bit64::as.integer64(c("9007199254740992", "9007199254740993"))
  kept <- prepare_B_training_data(data = d, trait_type = "binary",
                                  filter_include = list(tag = "9007199254740993"), verbose = 0)
  expect_identical(kept$data$rsID, "rs2")
  left <- prepare_B_training_data(data = d, trait_type = "binary",
                                  filter_exclude = list(tag = "9007199254740993"), verbose = 0)
  expect_identical(left$data$rsID, "rs1")
})

test_that("text values for a logical column are compared as logical values", {
  f <- GLOWr:::.filter_values_for_column
  expect_identical(f(c("true", "FALSE", "T", "F", "True", NA), c(TRUE, FALSE), "KEEP"),
                   c(TRUE, FALSE, TRUE, FALSE, TRUE, NA))
  expect_error(f(c("true", "yes"), c(TRUE, FALSE), "KEEP"),
               "column 'KEEP' has value\\(s\\) that are not true or false \\('yes'\\)")
  expect_identical(f(TRUE, c(TRUE, FALSE), "KEEP"), TRUE)

  d <- filter_test_data()
  d$KEEP <- c(TRUE, FALSE, TRUE, FALSE, NA, TRUE)
  expect_false(TRUE %in% "true")      # what R alone does
  kept <- GLOWr:::.apply_general_filters(d, include = list(KEEP = "true"), verbose = 0)
  expect_identical(kept$rsID, c("rs1", "rs3", "rs6"))
  # the same filter given as a logical value keeps the same rows
  expect_identical(GLOWr:::.apply_general_filters(d, include = list(KEEP = TRUE), verbose = 0)$rsID, kept$rsID)
  left <- GLOWr:::.apply_general_filters(d, exclude = list(KEEP = "false"), verbose = 0)
  expect_identical(left$rsID, c("rs1", "rs3", "rs5", "rs6"))
})

test_that("a value beyond the range of a 64-bit integer stops, and matches no other row", {
  skip_if_not_installed("bit64")
  f <- GLOWr:::.filter_values_for_column
  ids <- bit64::as.integer64(c("9223372036854775807", NA, "42"))     # the largest value, a missing one, an ordinary one
  # (review 2026-11-16, second round) Converted without a check, the first of these became the
  # largest value and the second a missing value, so each selected a row it does not name.
  expect_error(f("9223372036854775808", ids, "tag"), "not within the range of a 64-bit integer \\('9223372036854775808'\\)")
  expect_error(f("-9223372036854775809", ids, "tag"), "not within the range of a 64-bit integer")
  # (third review, 2026-11-16) zero written with a sign or with leading zeros is zero, not a
  # value beyond the range
  zeros <- bit64::as.integer64(c("0", "1"))
  for (z in c("0", "+0", "-0", "000", "-0000"))
    expect_identical(as.character(f(z, zeros, "id")), "0", info = z)
  expect_identical(as.character(f("-007", bit64::as.integer64(c("-7", "7")), "id")), "-7")
  d <- data.frame(rsID = c("rsMAX", "rsNA", "rsOK"), N = c(1, 2, 3), stringsAsFactors = FALSE); d$tag <- ids
  expect_error(GLOWr:::.apply_general_filters(d, include = list(tag = "9223372036854775808"), verbose = 0), "64-bit integer")
  expect_error(GLOWr:::.apply_general_filters(d, exclude = list(tag = "-9223372036854775809"), verbose = 0), "64-bit integer")
  # the two ends of the range, leading zeros and a plus sign are values like any other
  expect_identical(GLOWr:::.apply_general_filters(d, include = list(tag = "9223372036854775807"), verbose = 0)$rsID, "rsMAX")
  expect_identical(GLOWr:::.apply_general_filters(d, include = list(tag = c("0042", "+42")), verbose = 0)$rsID, "rsOK")
  expect_identical(GLOWr:::.apply_general_filters(d, exclude = list(tag = "42"), verbose = 0)$rsID, c("rsMAX", "rsNA"))
  expect_identical(as.character(f("-9223372036854775807", ids, "tag")), "-9223372036854775807")
})
