# ==============================================================================
# Phase 3: Post-Retirement Correctness Tests
# ==============================================================================
#
# Replaces Phase 2.5 equivalence gate tests. The old BSF_test, BSF_cctP_test,
# BSF_test_byP, BSF_cctP_test_byP functions have been removed (Phase 3).
# These tests verify that run_bsf_tests and the public API (GLOW_Omni,
# GLOW_Omni_byP, glow_test) produce correct, consistent results.
#
# Structure:
#   1. run_bsf_tests correctness (structure, dimensions, valid p-values)
#   2. Name correctness (16 hierarchical names for default specs)
#   3. GLOW_Omni / GLOW_Omni_byP consistency
#   4. Edge cases (all-zero B, extreme PI, m=2, single variant via glow_test)
#   5. Multiple random seeds
#   6. Integration (glow_test result structure and values)


# --- Helper: generate test data ---
.make_test_data <- function(m, seed) {
  set.seed(seed)
  n <- 100

  G <- matrix(rbinom(n * m, 2, 0.2), n, m)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.5)
  B <- rnorm(m, 0, 0.1)
  PI <- runif(m, 0.2, 0.8)

  marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")
  Bstar <- sqrt(diag(marg_stats$M_s)) * B / marg_stats$s0

  list(
    Zscores = marg_stats$Zscores,
    M = marg_stats$M_Z,
    Bstar = Bstar,
    PI = PI,
    B = B,
    marg_stats = marg_stats
  )
}


# ==============================================================================
# 1. run_bsf_tests CORRECTNESS (default specs, 16 rows)
# ==============================================================================

test_that("run_bsf_tests with default specs produces 16-row output", {
  td <- .make_test_data(m = 10, seed = 1001)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  # 16 rows: 5 SKAT + 3 Burden + 5 Fisher + BSF_Omni + SNV_CCT + Omni
  expect_equal(nrow(result$PVAL), 16)
  expect_equal(nrow(result$STAT), 16)
  expect_equal(ncol(result$PVAL), 1)
  expect_equal(ncol(result$STAT), 1)

  # All p-values in [0, 1]
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))

  # All statistics finite
  expect_true(all(is.finite(result$STAT)))

  # test_names field present and correct length
  expect_equal(length(result$test_names), 16)
})


test_that("run_bsf_tests with include_snv_cct=FALSE produces 14-row output", {
  td <- .make_test_data(m = 10, seed = 1002)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = FALSE
  )

  # 14 rows: 13 individual + BSF_Omni (no SNV_CCT, no Omni)
  expect_equal(nrow(result$PVAL), 14)
  expect_equal(nrow(result$STAT), 14)

  # Last row is BSF_Omni
  expect_equal(rownames(result$PVAL)[14], "BSF_Omni")

  # No SNV_CCT or Omni rows
  expect_false("SNV_CCT" %in% rownames(result$PVAL))
  expect_false("Omni" %in% rownames(result$PVAL))

  # All p-values valid
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("GLOW_Omni_byP produces same results as GLOW_Omni on same data", {
  td <- .make_test_data(m = 12, seed = 1003)

  # Reconstruct p-values from Z-scores
  Pvalues <- 2 * pnorm(-abs(td$Zscores))
  Zsigns <- sign(td$Zscores)
  Zsigns[Zsigns == 0] <- 1

  # GLOW_Omni path
  result_omni <- GLOW_Omni(td$marg_stats, td$B, td$PI)

  # GLOW_Omni_byP path
  result_byP <- GLOW_Omni_byP(
    Pvalues = Pvalues,
    Zsigns = Zsigns,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI
  )

  expect_equal(nrow(result_omni$PVAL), 16)
  expect_equal(nrow(result_byP$PVAL), 16)

  # Numerically equal p-values and statistics (not bit-identical because
  # GLOW_Omni_byP reconstructs Z-scores from p-values, introducing minor
  # floating-point differences from the direct Z-score path in GLOW_Omni)
  expect_equal(
    as.vector(unname(result_byP$PVAL)),
    as.vector(unname(result_omni$PVAL)),
    tolerance = 1e-12
  )
  expect_equal(
    as.vector(unname(result_byP$STAT)),
    as.vector(unname(result_omni$STAT)),
    tolerance = 1e-12
  )
})


# ==============================================================================
# 2. NAME CORRECTNESS: 16 hierarchical names for default specs
# ==============================================================================

test_that("default spec produces the expected 16 hierarchical names", {
  td <- .make_test_data(m = 8, seed = 2001)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores, M = td$M,
    Bstar = td$Bstar, PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expected_names <- c(
    "df_1_SKAT_BE_N", "df_1_SKAT_APR_N",
    "df_1_SKAT_BE_sparse", "df_1_SKAT_APR_sparse",
    "df_1_SKAT_equ",
    "df_Inf_Burden_BE", "df_Inf_Burden_APR", "df_Inf_Burden_equ",
    "df_2_Fisher_BE_N", "df_2_Fisher_APR_N",
    "df_2_Fisher_BE_sparse", "df_2_Fisher_APR_sparse",
    "df_2_Fisher_equ",
    "BSF_Omni", "SNV_CCT", "Omni"
  )

  expect_equal(rownames(result$PVAL), expected_names)
  expect_equal(length(rownames(result$PVAL)), 16)
})


test_that("glow_test output names have GLOW_ prefix", {
  td <- .make_test_data(m = 8, seed = 2002)

  result <- glow_test(td$marg_stats, td$B, td$PI, verbose = 0)

  expected_glow_names <- c(
    "GLOW_SKAT_BE_N", "GLOW_SKAT_APR_N",
    "GLOW_SKAT_BE_sparse", "GLOW_SKAT_APR_sparse",
    "GLOW_SKAT_equ",
    "GLOW_Burden_BE", "GLOW_Burden_APR", "GLOW_Burden_equ",
    "GLOW_Fisher_BE_N", "GLOW_Fisher_APR_N",
    "GLOW_Fisher_BE_sparse", "GLOW_Fisher_APR_sparse",
    "GLOW_Fisher_equ",
    "GLOW_BSF_Omni", "GLOW_SNV_CCT", "GLOW_Omni"
  )

  expect_equal(names(result$pvalues), expected_glow_names)
  expect_equal(names(result$statistics), expected_glow_names)
})


test_that("all 16 names are unique (no collisions)", {
  td <- .make_test_data(m = 8, seed = 2003)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores, M = td$M,
    Bstar = td$Bstar, PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(length(result$test_names), length(unique(result$test_names)))
})


# ==============================================================================
# 3. EDGE CASES
# ==============================================================================

test_that("all-zero B vector produces valid results", {
  set.seed(3001)
  n <- 100; m <- 8
  G <- matrix(rbinom(n * m, 2, 0.2), n, m)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.5)
  PI <- runif(m, 0.2, 0.8)

  marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")
  B_zero <- rep(0, m)
  Bstar_zero <- sqrt(diag(marg_stats$M_s)) * B_zero / marg_stats$s0

  result <- GLOWr:::run_bsf_tests(
    Zscores = marg_stats$Zscores, M = marg_stats$M_Z,
    Bstar = Bstar_zero, PI = PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("extreme PI (near 0) produces valid results", {
  set.seed(3002)
  n <- 100; m <- 8
  G <- matrix(rbinom(n * m, 2, 0.2), n, m)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.5)
  B <- rnorm(m, 0, 0.1)

  marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")
  Bstar <- sqrt(diag(marg_stats$M_s)) * B / marg_stats$s0
  PI_low <- rep(0.01, m)

  result <- GLOWr:::run_bsf_tests(
    Zscores = marg_stats$Zscores, M = marg_stats$M_Z,
    Bstar = Bstar, PI = PI_low,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("extreme PI (near 1) produces valid results", {
  set.seed(3003)
  n <- 100; m <- 8
  G <- matrix(rbinom(n * m, 2, 0.2), n, m)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.5)
  B <- rnorm(m, 0, 0.1)

  marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")
  Bstar <- sqrt(diag(marg_stats$M_s)) * B / marg_stats$s0
  PI_high <- rep(0.99, m)

  result <- GLOWr:::run_bsf_tests(
    Zscores = marg_stats$Zscores, M = marg_stats$M_Z,
    Bstar = Bstar, PI = PI_high,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


test_that("single variant (p=1) via glow_test produces valid results", {
  # Single-variant case is handled by glow_test, which assigns the marginal
  # p-value to all output slots.
  set.seed(3004)
  marg_stats <- list(
    Zscores = rnorm(1),
    M_Z = matrix(1, 1, 1),
    M_s = matrix(1, 1, 1),
    s0 = 1
  )
  B <- 0.5
  PI <- 0.3

  result <- glow_test(marg_stats, B, PI, verbose = 0)

  expect_s3_class(result, "glow_test_result")
  expect_equal(length(result$pvalues), 16)
  expect_true(result$settings$single_variant)

  # All p-values should equal the marginal p-value
  expected_pval <- 2 * pnorm(-abs(marg_stats$Zscores))
  expect_true(all(abs(result$pvalues - expected_pval) < 1e-12))

  # Names follow the naming convention
  expect_true("GLOW_Omni" %in% names(result$pvalues))
  expect_true("GLOW_BSF_Omni" %in% names(result$pvalues))
  expect_true("GLOW_SNV_CCT" %in% names(result$pvalues))
  expect_true("GLOW_SKAT_BE_N" %in% names(result$pvalues))
  expect_true("GLOW_Burden_BE" %in% names(result$pvalues))
  expect_true("GLOW_Fisher_BE_N" %in% names(result$pvalues))
})


test_that("m=2 (minimal multi-variant) produces valid 16-row results", {
  set.seed(3005)
  n <- 100; m <- 2
  G <- matrix(rbinom(n * m, 2, 0.2), n, m)
  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.5)
  B <- rnorm(m, 0, 0.1)
  PI <- runif(m, 0.2, 0.8)

  marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")
  Bstar <- sqrt(diag(marg_stats$M_s)) * B / marg_stats$s0

  result <- GLOWr:::run_bsf_tests(
    Zscores = marg_stats$Zscores, M = marg_stats$M_Z,
    Bstar = Bstar, PI = PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_equal(nrow(result$STAT), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


# ==============================================================================
# 4. INTEGRATION: glow_test result structure, names, and values
# ==============================================================================

test_that("glow_test omni returns valid complete results", {
  td <- .make_test_data(m = 10, seed = 4001)

  result <- glow_test(td$marg_stats, td$B, td$PI, verbose = 0)

  # 16 p-values
  expect_equal(length(result$pvalues), 16)

  # All valid
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
  expect_true(all(is.finite(result$statistics)))

  # GLOW_Omni is the last entry
  expect_equal(names(result$pvalues)[16], "GLOW_Omni")
})


test_that("glow_test result structure is complete", {
  td <- .make_test_data(m = 10, seed = 4002)

  result <- glow_test(
    td$marg_stats, td$B, td$PI,
    region_info = list(chr = "22", start = 1000, end = 2000, label = "TEST_GENE"),
    variant_summary = list(n_original = 15, n_after_collapse = 10, cMAC = 42),
    verbose = 0
  )

  # S3 class
  expect_s3_class(result, "glow_test_result")

  # All required fields present
  expect_true(all(c("pvalues", "statistics", "region_info",
                    "variant_summary", "settings", "raw") %in% names(result)))

  # Region info preserved
  expect_equal(result$region_info$label, "TEST_GENE")
  expect_equal(result$region_info$chr, "22")

  # Variant summary preserved
  expect_equal(result$variant_summary$n_original, 15)
  expect_equal(result$variant_summary$cMAC, 42)

  # Settings correct
  expect_equal(result$settings$tests_run, "omni")

  # Raw output present and has correct dimensions
  expect_true("omni" %in% names(result$raw))
  expect_equal(nrow(result$raw$omni$PVAL), 16)

  # as.data.frame works
  df <- as.data.frame(result)
  expect_equal(nrow(df), 1)
  expect_true("GLOW_Omni" %in% names(df))
  expect_equal(df$label, "TEST_GENE")
  n_pval_cols <- sum(grepl("^GLOW_", names(df)))
  expect_equal(n_pval_cols, 16)
})


test_that("glow_test pvalues and statistics have identical names", {
  td <- .make_test_data(m = 10, seed = 4003)
  result <- glow_test(td$marg_stats, td$B, td$PI, verbose = 0)
  expect_identical(names(result$pvalues), names(result$statistics))
})


# ==============================================================================
# 5. MULTIPLE RANDOM SEEDS (consistency across different data)
# ==============================================================================

test_that("run_bsf_tests produces valid results across seed=5001 (m=8)", {
  td <- .make_test_data(m = 8, seed = 5001)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores, M = td$M,
    Bstar = td$Bstar, PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


test_that("run_bsf_tests produces valid results across seed=5002 (m=15)", {
  td <- .make_test_data(m = 15, seed = 5002)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores, M = td$M,
    Bstar = td$Bstar, PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


test_that("run_bsf_tests produces valid results across seed=5003 (m=25)", {
  td <- .make_test_data(m = 25, seed = 5003)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores, M = td$M,
    Bstar = td$Bstar, PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


test_that("run_bsf_tests produces valid results across seed=5004 (m=5, small)", {
  td <- .make_test_data(m = 5, seed = 5004)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores, M = td$M,
    Bstar = td$Bstar, PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})


test_that("run_bsf_tests with correlated genotypes (seed=5005)", {
  set.seed(5005)
  n <- 100; m <- 10

  # Generate correlated genotypes via multivariate normal
  library(MASS)
  cor_mat <- matrix(0.3, m, m)
  diag(cor_mat) <- 1
  Z_corr <- mvrnorm(n, mu = rep(0, m), Sigma = cor_mat)
  G <- matrix(as.integer(pnorm(Z_corr) < 0.2), n, m) +
    matrix(as.integer(pnorm(Z_corr) < 0.1), n, m)

  X <- matrix(rnorm(n * 2), n, 2)
  Y <- rbinom(n, 1, 0.5)
  B <- rnorm(m, 0, 0.1)
  PI <- runif(m, 0.2, 0.8)

  marg_stats <- getZ_marg_score(G, X, Y, trait = "binary")
  Bstar <- sqrt(diag(marg_stats$M_s)) * B / marg_stats$s0

  result <- GLOWr:::run_bsf_tests(
    Zscores = marg_stats$Zscores, M = marg_stats$M_Z,
    Bstar = Bstar, PI = PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(nrow(result$PVAL), 16)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
  expect_true(all(is.finite(result$STAT)))
})
