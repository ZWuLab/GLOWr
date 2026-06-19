# ==============================================================================
# Tests for glow_composition.R
# ==============================================================================
#
# Tests for call_weight_fn(), run_bsf_tests(), and generate_test_names().
# The CRITICAL test is numerical equivalence between run_bsf_tests() and
# the existing GLOW_Omni() pipeline.


# --- Helper: create test data for composition tests ---
.make_composition_test_data <- function(m = 10, seed = 100) {
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


# ==================== call_weight_fn() ====================

test_that("call_weight_fn dispatches to function with specific formals", {
  # Function that takes only 'p'
  simple_fn <- function(p) list(equ = rep(1, p))
  all_args <- list(g = function(x) x, Bstar = c(0.5), PI = c(0.5),
                   M = matrix(1), p = 3, is.posi.wts = TRUE)

  result <- GLOWr:::call_weight_fn(simple_fn, all_args)
  expect_named(result, "equ")
  expect_equal(result$equ, rep(1, 3))
})

test_that("call_weight_fn passes all args to function with ...", {
  # Function with ... receives everything
  dotfun <- function(...) {
    args <- list(...)
    list(got = names(args))
  }
  all_args <- list(a = 1, b = 2, c = 3)

  result <- GLOWr:::call_weight_fn(dotfun, all_args)
  expect_true(all(c("a", "b", "c") %in% result$got))
})

test_that("call_weight_fn works with Optimal_Weights_M signature", {
  # Optimal_Weights_M(g, Bstar, PI, M, is.posi.wts = TRUE)
  # Use p=2 to avoid diag() scalar issue with p=1
  all_args <- list(
    g = function(x) x,
    Bstar = c(0.5, 0.5),
    PI = c(0.5, 0.5),
    M = diag(2),
    p = 2,
    is.posi.wts = TRUE
  )

  result <- GLOWr:::call_weight_fn(Optimal_Weights_M, all_args)
  # Burden identity function -> returns wts_BE, wts_APE
  expect_true("wts_BE" %in% names(result))
  expect_true("wts_APE" %in% names(result))
})


# ==================== generate_test_names() ====================

test_that("generate_test_names produces 16 names for default specs", {
  specs <- default_test_specs()
  names <- GLOWr:::generate_test_names(specs, include_snv_cct = TRUE)

  expect_equal(length(names), 16)

  # Check structure: SKAT(5) + Burden(3) + Fisher(5) + BSF_Omni + SNV_CCT + Omni
  # Names include df_ prefix from multi_SgZ_test convention
  # SKAT (df=1): BE_N, APR_N, BE_sparse, APR_sparse, equ
  expect_true("df_1_SKAT_BE_N" %in% names)
  expect_true("df_1_SKAT_APR_N" %in% names)
  expect_true("df_1_SKAT_BE_sparse" %in% names)
  expect_true("df_1_SKAT_APR_sparse" %in% names)
  expect_true("df_1_SKAT_equ" %in% names)

  # Burden (df=Inf): BE, APR, equ
  expect_true("df_Inf_Burden_BE" %in% names)
  expect_true("df_Inf_Burden_APR" %in% names)
  expect_true("df_Inf_Burden_equ" %in% names)

  # Fisher (df=2): BE_N, APR_N, BE_sparse, APR_sparse, equ
  expect_true("df_2_Fisher_BE_N" %in% names)
  expect_true("df_2_Fisher_APR_N" %in% names)
  expect_true("df_2_Fisher_BE_sparse" %in% names)
  expect_true("df_2_Fisher_APR_sparse" %in% names)
  expect_true("df_2_Fisher_equ" %in% names)

  # CCT rows
  expect_true("BSF_Omni" %in% names)
  expect_true("SNV_CCT" %in% names)
  expect_true("Omni" %in% names)
})

test_that("generate_test_names respects include_snv_cct = FALSE", {
  specs <- default_test_specs()
  names_without <- GLOWr:::generate_test_names(specs, include_snv_cct = FALSE)

  # Without SNV_CCT: 13 individual + BSF_Omni = 14
  expect_equal(length(names_without), 14)
  expect_true("BSF_Omni" %in% names_without)
  expect_false("SNV_CCT" %in% names_without)
  expect_false("Omni" %in% names_without)
})

test_that("generate_test_names rejects duplicate family names", {
  specs <- list(
    list(family = "SKAT", g = function(x) x^2, df = 1),
    list(family = "SKAT", g = function(x) x, df = Inf)
  )

  expect_error(
    GLOWr:::generate_test_names(specs),
    "Duplicate family names"
  )
})

test_that("generate_test_names works with custom 2-family spec", {
  specs <- list(
    GLOWr:::make_test_spec("TestA", function(x) x, Inf),
    GLOWr:::make_test_spec("TestB", function(x) x^2, 1)
  )

  names <- GLOWr:::generate_test_names(specs, include_snv_cct = TRUE)

  # TestA (Burden identity, df=Inf) -> 2 optimal + 1 equal = 3
  # TestB (SKAT-like, df=1) -> 4 optimal + 1 equal = 5
  # + BSF_Omni + SNV_CCT + Omni = 11
  expect_equal(length(names), 11)
  expect_true(any(grepl("TestA_", names)))
  expect_true(any(grepl("TestB_", names)))
  expect_true(any(grepl("^df_Inf_TestA_", names)))
  expect_true(any(grepl("^df_1_TestB_", names)))
})

test_that("generate_test_names all names are unique", {
  specs <- default_test_specs()
  names <- GLOWr:::generate_test_names(specs, include_snv_cct = TRUE)
  expect_equal(length(names), length(unique(names)))
})


# ==================== run_bsf_tests() ====================

test_that("run_bsf_tests returns correct structure", {
  td <- .make_composition_test_data(m = 8, seed = 200)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_type(result, "list")
  expect_true("STAT" %in% names(result))
  expect_true("PVAL" %in% names(result))
  expect_true("test_names" %in% names(result))

  # 16 rows: 13 individual + BSF_Omni + SNV_CCT + Omni
  expect_equal(nrow(result$STAT), 16)
  expect_equal(nrow(result$PVAL), 16)
  expect_equal(ncol(result$STAT), 1)
  expect_equal(ncol(result$PVAL), 1)

  # All p-values in [0, 1]
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))

  # test_names matches rownames
  expect_equal(result$test_names, rownames(result$PVAL))
})

test_that("run_bsf_tests with include_snv_cct=FALSE omits SNV rows", {
  td <- .make_composition_test_data(m = 8, seed = 201)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = FALSE
  )

  # 13 individual + BSF_Omni = 14
  expect_equal(nrow(result$PVAL), 14)
  expect_true("BSF_Omni" %in% rownames(result$PVAL))
  expect_false("SNV_CCT" %in% rownames(result$PVAL))
  expect_false("Omni" %in% rownames(result$PVAL))
})

test_that("run_bsf_tests rejects duplicate family names", {
  td <- .make_composition_test_data(m = 5, seed = 202)

  bad_specs <- list(
    list(family = "SKAT", g = function(x) x^2, df = 1),
    list(family = "SKAT", g = function(x) x, df = Inf)
  )

  expect_error(
    GLOWr:::run_bsf_tests(
      Zscores = td$Zscores,
      M = td$M,
      Bstar = td$Bstar,
      PI = td$PI,
      test_specs = bad_specs
    ),
    "Duplicate family names"
  )
})

test_that("run_bsf_tests works with custom 2-family spec", {
  td <- .make_composition_test_data(m = 6, seed = 203)

  # Only Burden + SKAT (no Fisher)
  specs <- list(
    GLOWr:::make_test_spec("Burden", function(x) x, Inf),
    GLOWr:::make_test_spec("SKAT", function(x) x^2, 1)
  )

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = specs,
    include_snv_cct = TRUE
  )

  # Burden(3) + SKAT(5) + BSF_Omni + SNV_CCT + Omni = 11
  expect_equal(nrow(result$PVAL), 11)
  expect_true(all(result$PVAL >= 0 & result$PVAL <= 1))
})


# ==================== CRITICAL: Numerical equivalence ====================

test_that("run_bsf_tests produces identical p-values to GLOW_Omni", {
  # This is the most important test in this file.
  # run_bsf_tests with default_test_specs() must produce numerically
  # identical results to the existing GLOW_Omni pipeline.

  td <- .make_composition_test_data(m = 10, seed = 300)

  # --- Run existing GLOW_Omni ---
  omni_result <- GLOW_Omni(
    marg_score_stats = td$marg_stats,
    B = td$B,
    PI = td$PI
  )

  # --- Run new run_bsf_tests ---
  comp_result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  # Both should have 16 rows
  expect_equal(nrow(omni_result$PVAL), 16)
  expect_equal(nrow(comp_result$PVAL), 16)

  # P-values must be numerically identical (same code paths, same
  # floating-point operations). Compare unname'd values since the row
  # names intentionally differ (new: hierarchical, old: legacy).
  for (i in 1:16) {
    expect_identical(
      unname(comp_result$PVAL[i, 1]),
      unname(omni_result$PVAL[i, 1]),
      info = paste0("Row ", i, " p-value mismatch: ",
                    "new=", comp_result$PVAL[i, 1],
                    " old=", omni_result$PVAL[i, 1])
    )
  }

  # Statistics must be numerically identical
  for (i in 1:16) {
    expect_identical(
      unname(comp_result$STAT[i, 1]),
      unname(omni_result$STAT[i, 1]),
      info = paste0("Row ", i, " stat mismatch: ",
                    "new=", comp_result$STAT[i, 1],
                    " old=", omni_result$STAT[i, 1])
    )
  }
})

test_that("run_bsf_tests equivalence holds for different data sizes", {
  # Test with m = 5
  td5 <- .make_composition_test_data(m = 5, seed = 301)

  omni5 <- GLOW_Omni(td5$marg_stats, td5$B, td5$PI)
  comp5 <- GLOWr:::run_bsf_tests(
    Zscores = td5$Zscores, M = td5$M,
    Bstar = td5$Bstar, PI = td5$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_identical(as.vector(comp5$PVAL), as.vector(omni5$PVAL))
  expect_identical(as.vector(comp5$STAT), as.vector(omni5$STAT))

  # Test with m = 20
  td20 <- .make_composition_test_data(m = 20, seed = 302)

  omni20 <- GLOW_Omni(td20$marg_stats, td20$B, td20$PI)
  comp20 <- GLOWr:::run_bsf_tests(
    Zscores = td20$Zscores, M = td20$M,
    Bstar = td20$Bstar, PI = td20$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_identical(as.vector(comp20$PVAL), as.vector(omni20$PVAL))
  expect_identical(as.vector(comp20$STAT), as.vector(omni20$STAT))
})

test_that("run_bsf_tests equivalence holds with correlated genotypes", {
  set.seed(303)
  n <- 100
  m <- 8

  # Generate correlated genotypes
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

  omni <- GLOW_Omni(marg_stats, B, PI)
  comp <- GLOWr:::run_bsf_tests(
    Zscores = marg_stats$Zscores, M = marg_stats$M_Z,
    Bstar = Bstar, PI = PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_identical(as.vector(comp$PVAL), as.vector(omni$PVAL))
  expect_identical(as.vector(comp$STAT), as.vector(omni$STAT))
})


# ==================== Name correctness ====================

test_that("run_bsf_tests row names match generate_test_names", {
  td <- .make_composition_test_data(m = 8, seed = 400)

  result <- GLOWr:::run_bsf_tests(
    Zscores = td$Zscores,
    M = td$M,
    Bstar = td$Bstar,
    PI = td$PI,
    test_specs = default_test_specs(),
    include_snv_cct = TRUE
  )

  expected_names <- GLOWr:::generate_test_names(
    default_test_specs(),
    include_snv_cct = TRUE
  )

  expect_equal(result$test_names, expected_names)
})
