# ==============================================================================
# Tests for glow_test(), print.glow_test_result, as.data.frame.glow_test_result
# ==============================================================================
#
# Tests the GLOW variant-set test runner that wraps GLOW_Omni and individual
# tests, returning structured S3 result objects.

# --- Helper: create mock marg_score_stats for unit testing ---
# Uses identity M_Z and M_s to avoid needing real genotype data.
.make_mock_marg_stats <- function(p, seed = 42) {
  set.seed(seed)
  list(
    Zscores = rnorm(p),
    M_Z = diag(p),
    M_s = diag(p),
    s0 = 1
  )
}


# ==================== S3 class and field structure ====================

test_that("glow_test returns correct S3 class with all fields", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 42)
  set.seed(42)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI,
                      region_info = list(chr = "22", start = 1000,
                                         end = 2000, label = "GENE1"),
                      verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_true("pvalues" %in% names(result))
  expect_true("statistics" %in% names(result))
  expect_true("region_info" %in% names(result))
  expect_true("variant_summary" %in% names(result))
  expect_true("settings" %in% names(result))
  expect_true("raw" %in% names(result))

  # All p-values in [0, 1]
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
})


# ==================== Omni mode ====================

test_that("glow_test omni mode returns 16 p-values with correct names", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 43)
  set.seed(43)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "omni", verbose = 0)

  expect_equal(length(result$pvalues), 16)
  expect_equal(length(result$statistics), 16)

  # Check key standardized names exist
  expect_true("GLOW_Omni" %in% names(result$pvalues))
  expect_true("GLOW_Burden_BE" %in% names(result$pvalues))
  expect_true("GLOW_SKAT_BE_N" %in% names(result$pvalues))
  expect_true("GLOW_Fisher_BE_N" %in% names(result$pvalues))
  expect_true("GLOW_BSF_Omni" %in% names(result$pvalues))
  expect_true("GLOW_SNV_CCT" %in% names(result$pvalues))

  # Statistics vector has same names as pvalues
  expect_identical(names(result$statistics), names(result$pvalues))
})

test_that("glow_test omni stores raw GLOW_Omni output", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 44)
  set.seed(44)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "omni", verbose = 0)

  expect_true("omni" %in% names(result$raw))
  expect_true(is.matrix(result$raw$omni$STAT))
  expect_true(is.matrix(result$raw$omni$PVAL))
  expect_equal(nrow(result$raw$omni$PVAL), 16)
})


# ==================== as.data.frame ====================

test_that("glow_test as.data.frame produces single row with all columns", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 45)
  set.seed(45)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI,
                      region_info = list(chr = "22", start = 1, end = 100,
                                         label = "TEST"),
                      variant_summary = list(n_original = 20,
                                             n_after_collapse = 10,
                                             cMAC = 50),
                      verbose = 0)
  df <- as.data.frame(result)

  expect_equal(nrow(df), 1)

  # Region info columns
  expect_true("label" %in% names(df))
  expect_true("chr" %in% names(df))
  expect_true("start" %in% names(df))
  expect_true("end" %in% names(df))
  expect_equal(df$label, "TEST")
  expect_equal(df$chr, "22")

  # Variant summary columns
  expect_true("n_variants" %in% names(df))
  expect_true("cMAC" %in% names(df))
  expect_equal(df$n_variants, 20)
  expect_equal(df$cMAC, 50)

  # P-value columns
  expect_true("GLOW_Omni" %in% names(df))
  expect_true("GLOW_Burden_BE" %in% names(df))
})

test_that("glow_test as.data.frame handles NULL region/summary gracefully", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 46)
  set.seed(46)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, verbose = 0)
  df <- as.data.frame(result)

  expect_equal(nrow(df), 1)
  expect_true(is.na(df$label))
  expect_true(is.na(df$chr))
  expect_true(is.na(df$n_variants))
  expect_true(is.na(df$cMAC))

  # P-values should still be present
  expect_true("GLOW_Omni" %in% names(df))
})


# ==================== Individual test selection ====================

test_that("glow_test handles individual test selection (burden)", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 47)
  set.seed(47)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "burden", verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_true(length(result$pvalues) > 0)
  # All p-values should have Burden-related names (df_Inf_ or GLOW_Burden)
  expect_true(all(grepl("Burden|Inf", names(result$pvalues))))
  expect_true("burden" %in% names(result$raw))
  expect_null(result$raw$omni)

  # P-values in [0, 1]
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
})

test_that("glow_test handles individual test selection (skat)", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 48)
  set.seed(48)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "skat", verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_true(length(result$pvalues) > 0)
  # All names should have SKAT-related names (df_1_ or GLOW_SKAT)
  expect_true(all(grepl("SKAT|df_1", names(result$pvalues))))
  expect_true("skat" %in% names(result$raw))
})

test_that("glow_test handles multiple individual tests", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 49)
  set.seed(49)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI,
                      tests = c("burden", "skat"),
                      verbose = 0)

  expect_s3_class(result, "glow_test_result")
  # Should have both burden and skat p-values
  expect_true("burden" %in% names(result$raw))
  expect_true("skat" %in% names(result$raw))
  expect_true(any(grepl("Burden|Inf", names(result$pvalues))))
  expect_true(any(grepl("SKAT|df_1", names(result$pvalues))))
})


# ==================== Input validation ====================

test_that("glow_test validates B length mismatch", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 50)

  expect_error(
    glow_test(marg_stats, B = rep(0.5, p - 1), PI = runif(p), verbose = 0),
    "B must have length"
  )
})

test_that("glow_test validates PI length mismatch", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 51)

  expect_error(
    glow_test(marg_stats, B = rep(0.5, p), PI = runif(p + 2), verbose = 0),
    "PI must have length"
  )
})

test_that("glow_test validates PI range", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 52)

  expect_error(
    glow_test(marg_stats, B = rep(0.5, p), PI = rep(-0.1, p), verbose = 0),
    "PI values must be in"
  )
  expect_error(
    glow_test(marg_stats, B = rep(0.5, p), PI = rep(1.5, p), verbose = 0),
    "PI values must be in"
  )
})

test_that("glow_test validates NA in B", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 53)
  B_na <- rep(0.5, p)
  B_na[3] <- NA

  expect_error(
    glow_test(marg_stats, B = B_na, PI = runif(p, 0.1, 0.9), verbose = 0),
    "B must not contain NA"
  )
})

test_that("glow_test validates missing marg_score_stats fields", {
  p <- 10
  bad_stats <- list(Zscores = rnorm(p))  # missing M_Z, M_s, s0

  expect_error(
    glow_test(bad_stats, B = rep(0.5, p), PI = runif(p), verbose = 0),
    "marg_score_stats is missing required fields"
  )
})

test_that("glow_test validates invalid test names", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 54)

  expect_error(
    glow_test(marg_stats, B = rep(0.5, p), PI = runif(p, 0.1, 0.9),
              tests = "bogus", verbose = 0),
    "Invalid test"
  )
})


# ==================== Print method ====================

test_that("glow_test print method works with region label", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 55)
  set.seed(55)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI,
                      region_info = list(label = "TEST_GENE"),
                      verbose = 0)

  expect_output(print(result), "GLOW Test Result")
  expect_output(print(result), "TEST_GENE")
  expect_output(print(result), "Top p-values")
})

test_that("glow_test print method works without region info", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 56)
  set.seed(56)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, verbose = 0)

  # Should print without error
  expect_output(print(result), "GLOW Test Result")
  expect_output(print(result), "omni")
})

test_that("glow_test print shows variant summary when provided", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 57)
  set.seed(57)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI,
                      variant_summary = list(n_original = 20,
                                             n_after_collapse = 10,
                                             cMAC = 50),
                      verbose = 0)

  expect_output(print(result), "Variants: 10")
  expect_output(print(result), "Cumulative MAC: 50")
})


# ==================== Settings stored correctly ====================

test_that("glow_test stores settings correctly", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 58)
  set.seed(58)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "omni", verbose = 0)
  expect_equal(result$settings$tests_run, "omni")

  result2 <- glow_test(marg_stats, B, PI,
                       tests = c("burden", "fisher"), verbose = 0)
  expect_equal(result2$settings$tests_run, c("burden", "fisher"))
})


# ==================== Verbose messaging ====================

test_that("glow_test produces message when verbose >= 1", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 59)
  set.seed(59)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  expect_message(
    glow_test(marg_stats, B, PI, verbose = 1),
    "GLOW test: p=10"
  )
})

test_that("glow_test is silent when verbose = 0", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 60)
  set.seed(60)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  expect_silent(glow_test(marg_stats, B, PI, verbose = 0))
})


# ==================== Single-variant edge case ====================

test_that("glow_test handles single variant (p=1) correctly", {
  marg_stats <- .make_mock_marg_stats(1, seed = 70)
  B <- 0.5
  PI <- 0.3

  result <- glow_test(marg_stats, B, PI, verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_equal(length(result$pvalues), 16)

  # All p-values should equal the marginal p-value
  expected_pval <- 2 * pnorm(-abs(marg_stats$Zscores))
  expect_true(all(abs(result$pvalues - expected_pval) < 1e-12))

  # Settings should indicate single variant
  expect_true(result$settings$single_variant)
})

test_that("glow_test single variant as.data.frame works", {
  marg_stats <- .make_mock_marg_stats(1, seed = 71)

  result <- glow_test(marg_stats, B = 0.5, PI = 0.3,
                      region_info = list(label = "SINGLE"),
                      verbose = 0)
  df <- as.data.frame(result)
  expect_equal(nrow(df), 1)
  expect_true("GLOW_Omni" %in% names(df))
  expect_equal(df$label, "SINGLE")
})


# ==================== Custom test_specs (Phase 2) ====================

test_that("glow_test with custom 2-family test_specs returns dynamic names", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 80)
  set.seed(80)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  # Custom spec: only SKAT + Burden (no Fisher)
  my_specs <- list(
    list(family = "SKAT", g = function(x) x^2, df = 1),
    list(family = "Burden", g = function(x) x, df = Inf)
  )

  result <- glow_test(marg_stats, B, PI,
                      test_specs = my_specs,
                      verbose = 0)

  expect_s3_class(result, "glow_test_result")

  # SKAT(5) + Burden(3) + BSF_Omni + SNV_CCT + Omni = 11
  expect_equal(length(result$pvalues), 11)
  expect_equal(length(result$statistics), 11)

  # Should have SKAT and Burden names, but NOT Fisher
  expect_true(any(grepl("^GLOW_SKAT_", names(result$pvalues))))
  expect_true(any(grepl("^GLOW_Burden_", names(result$pvalues))))
  expect_false(any(grepl("^GLOW_Fisher_", names(result$pvalues))))

  # Should have omnibus names
  expect_true("GLOW_BSF_Omni" %in% names(result$pvalues))
  expect_true("GLOW_SNV_CCT" %in% names(result$pvalues))
  expect_true("GLOW_Omni" %in% names(result$pvalues))

  # All p-values valid
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))

  # Statistics names match p-value names
  expect_identical(names(result$statistics), names(result$pvalues))
})

test_that("glow_test custom specs as.data.frame has dynamic columns", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 81)
  set.seed(81)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  my_specs <- list(
    list(family = "SKAT", g = function(x) x^2, df = 1),
    list(family = "Burden", g = function(x) x, df = Inf)
  )

  result <- glow_test(marg_stats, B, PI, test_specs = my_specs,
                      region_info = list(label = "CUSTOM"),
                      verbose = 0)
  df <- as.data.frame(result)

  expect_equal(nrow(df), 1)
  expect_equal(df$label, "CUSTOM")

  # Should have GLOW_Omni but not 16 p-value columns
  expect_true("GLOW_Omni" %in% names(df))
  n_pval_cols <- sum(grepl("^GLOW_", names(df)))
  expect_equal(n_pval_cols, 11)
})


# ==================== include_snv_cct = FALSE (Phase 2) ====================

test_that("glow_test with include_snv_cct=FALSE omits SNV and Omni", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 82)
  set.seed(82)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI,
                      include_snv_cct = FALSE,
                      verbose = 0)

  expect_s3_class(result, "glow_test_result")

  # 13 individual + BSF_Omni = 14 (no SNV_CCT, no Omni)
  expect_equal(length(result$pvalues), 14)

  expect_true("GLOW_BSF_Omni" %in% names(result$pvalues))
  expect_false("GLOW_SNV_CCT" %in% names(result$pvalues))
  expect_false("GLOW_Omni" %in% names(result$pvalues))

  # All p-values valid
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
})


# ==================== Single-variant with custom specs (Phase 2) ====================

test_that("glow_test single variant with custom specs returns correct names", {
  marg_stats <- .make_mock_marg_stats(1, seed = 83)
  B <- 0.5
  PI <- 0.3

  # Custom 2-family spec
  my_specs <- list(
    list(family = "SKAT", g = function(x) x^2, df = 1),
    list(family = "Burden", g = function(x) x, df = Inf)
  )

  result <- glow_test(marg_stats, B, PI,
                      test_specs = my_specs,
                      verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_true(result$settings$single_variant)

  # SKAT(5) + Burden(3) + BSF_Omni + SNV_CCT + Omni = 11
  expect_equal(length(result$pvalues), 11)

  # All p-values should equal the marginal p-value
  expected_pval <- 2 * pnorm(-abs(marg_stats$Zscores))
  expect_true(all(abs(result$pvalues - expected_pval) < 1e-12))

  # Should have SKAT and Burden names, but NOT Fisher
  expect_true(any(grepl("^GLOW_SKAT_", names(result$pvalues))))
  expect_true(any(grepl("^GLOW_Burden_", names(result$pvalues))))
  expect_false(any(grepl("^GLOW_Fisher_", names(result$pvalues))))

  # Omni names present
  expect_true("GLOW_Omni" %in% names(result$pvalues))
})

test_that("glow_test single variant with include_snv_cct=FALSE", {
  marg_stats <- .make_mock_marg_stats(1, seed = 84)
  B <- 0.5
  PI <- 0.3

  result <- glow_test(marg_stats, B, PI,
                      include_snv_cct = FALSE,
                      verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_true(result$settings$single_variant)

  # 13 individual + BSF_Omni = 14 (no SNV_CCT, no Omni)
  expect_equal(length(result$pvalues), 14)

  expect_true("GLOW_BSF_Omni" %in% names(result$pvalues))
  expect_false("GLOW_SNV_CCT" %in% names(result$pvalues))
  expect_false("GLOW_Omni" %in% names(result$pvalues))

  # All p-values should equal the marginal p-value
  expected_pval <- 2 * pnorm(-abs(marg_stats$Zscores))
  expect_true(all(abs(result$pvalues - expected_pval) < 1e-12))
})


# ==================== return_weights ====================

# Expected 13 scheme rownames produced by run_bsf_tests with default_test_specs().
.expected_default_scheme_names <- c(
  "SKAT_BE_N", "SKAT_BE_sparse", "SKAT_APR_N", "SKAT_APR_sparse", "SKAT_equ",
  "Burden_BE", "Burden_APR", "Burden_equ",
  "Fisher_BE_N", "Fisher_BE_sparse", "Fisher_APR_N", "Fisher_APR_sparse",
  "Fisher_equ"
)

test_that("glow_test default does not attach weights field", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 101)
  set.seed(101)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "omni", verbose = 0)
  expect_false("weights" %in% names(result))
})

test_that("glow_test(return_weights = TRUE) attaches weight matrix with expected shape", {
  p <- 10
  marg_stats <- .make_mock_marg_stats(p, seed = 102)
  set.seed(102)
  B <- rep(0.5, p)
  PI <- runif(p, 0.1, 0.9)

  result <- glow_test(marg_stats, B, PI, tests = "omni",
                      return_weights = TRUE, verbose = 0)

  expect_true("weights" %in% names(result))
  expect_true(is.matrix(result$weights))
  expect_equal(nrow(result$weights), 13L)
  expect_equal(ncol(result$weights), p)
  expect_setequal(rownames(result$weights), .expected_default_scheme_names)
})

test_that("glow_test single-variant + return_weights attaches uniform 13x1 matrix of 1s", {
  marg_stats <- .make_mock_marg_stats(1L, seed = 103)
  result <- glow_test(marg_stats, B = 0.5, PI = 0.5,
                      return_weights = TRUE, verbose = 0)

  expect_true("weights" %in% names(result))
  expect_true(is.matrix(result$weights))
  expect_equal(nrow(result$weights), 13L)
  expect_equal(ncol(result$weights), 1L)
  expect_true(all(result$weights == 1))
  expect_setequal(rownames(result$weights), .expected_default_scheme_names)
})

test_that("return_weights default-off works in single-variant edge case", {
  marg_stats <- .make_mock_marg_stats(1L, seed = 104)
  result <- glow_test(marg_stats, B = 0.5, PI = 0.5, verbose = 0)
  expect_false("weights" %in% names(result))
})
